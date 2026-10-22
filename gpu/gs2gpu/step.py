"""One GS2 time step: ``fields::advance`` -> ``advance_implicit``/``advance_local``
-> ``timeadv`` (spec 01 sections 1.6-1.10 and 3.2-3.3, spec 10 section C).

The step composes the kernel modules:

=====================  ===================================================
stage                  module
=====================  ===================================================
flow shear             ``exb.exb_shear``
explicit / NL terms    ``explicit.add_explicit_terms`` (+ ``nonlinear``)
source                 ``source.get_source_term_on_local_domain``
parallel solve         ``invert_rhs.invert_rhs`` (+ linked BCs, parity)
hyperviscosity         ``hyper.hyper_diff``
collisions             ``collisions.vspace_derivatives`` / ``collisions_advance``
field solve            ``fields.getfield`` / ``fields_local.getfield_local``
zonal-flow removal     ``zonal.remove_zonal_flows``
=====================  ===================================================

Everything that depends on ``code_dt`` lives in a :class:`DtSet`, rebuilt by
:meth:`Model.build_dtset` when the time step changes. Host-side control
(resets, counters, time block) is plain Python; the arrays stay on the device.
"""

from __future__ import annotations

import dataclasses
from dataclasses import dataclass
from functools import partial
from typing import Any, NamedTuple, Optional

import jax
import jax.numpy as jnp
import numpy as np

from . import collisions as coll
from . import collisions_init as coll_init
from . import coeffs as coeffs_mod
from . import exb as exb_mod
from . import explicit as explicit_mod
from . import fields as fields_mod
from . import fields_local as fields_local_mod
from . import hyper as hyper_mod
from . import invert_rhs as invert_mod
from . import nonlinear as nl_mod
from . import parity as parity_mod
from . import response as response_mod
from . import source as source_mod
from . import zonal as zonal_mod
from .gadjust import FROM_G, TO_G, GAdjustData, g_adjust_direction
from .gates import check_gates
from .io.export import Export
from .reset import ReinitParams, ResetCounters
from .types import Coeffs, Dims, Fields, TimeBlock, fields_from_group, to_jax


def _flag(group, name: str, default: bool = False) -> bool:
    if group is None or name not in group:
        return default
    return bool(np.asarray(group[name]).item())


@dataclass(frozen=True)
class RunParams:
    """Run-control switches read by the step and the driver."""

    fphi: float
    fapar: float
    fbpar: float
    immediate_reset: bool
    nstep: int
    max_sim_time: float
    ncheck_stop: int
    nonlin: bool
    local_fields: bool
    remove_zonal_flows: bool
    force_maxwell_reinit: bool
    def_parity: bool
    colls: bool
    split_collisions: bool
    exb_active: bool

    @classmethod
    def from_export(cls, ex: Export, local_fields: bool) -> "RunParams":
        rp = ex["/run_parameters"]
        fl = ex["/fields"]
        cg = ex["/collisions"] if "collisions" in ex.root.groups else None
        st = ex["/state_t0"]
        colls = _flag(cg, "colls")
        return cls(
            fphi=float(rp["fphi"]), fapar=float(rp["fapar"]), fbpar=float(rp["fbpar"]),
            immediate_reset=_flag(rp, "immediate_reset", True),
            nstep=int(rp["nstep"]), max_sim_time=float(rp["max_sim_time"]),
            ncheck_stop=int(rp["ncheck_stop"]),
            nonlin=_flag(ex["/nonlinear_terms"], "nonlin"),
            local_fields=local_fields,
            remove_zonal_flows=_flag(fl, "remove_zonal_flows_switch"),
            force_maxwell_reinit=_flag(fl, "force_maxwell_reinit", True),
            def_parity=_flag(ex["/dist_fn"], "def_parity"),
            colls=colls,
            split_collisions=_flag(cg, "split_collisions"),
            exb_active=("kx_shift" in st) or ("theta0_shift" in st),
        )


class DtSet(NamedTuple):
    """Every array that depends on ``code_dt`` (spec 10 section E). A pytree, so it can
    be passed to jitted functions."""

    code_dt: float
    coeffs: Coeffs                         # dist_fn coefficients incl. g_h (device)
    coll: Optional[Any]                    # collisions_init.CollCoeffs (device) or None
    aminv: Any                             # response.Aminv (implicit or converted local)


@dataclass
class State:
    """The evolving state of the run (the Fortran module arrays and SAVE counters)."""

    g: jnp.ndarray
    gnew: jnp.ndarray
    old: Fields                            # phi, apar, bpar
    new: Fields                            # phinew, aparnew, bparnew
    time: TimeBlock
    reset: bool
    xstate: Any                            # explicit.ExplicitState (gexp_1..3, nl istep_last)
    exb: Optional[Any]                     # exb.ExbState or None
    vnmult: np.ndarray
    timeadv_istep_last: int = -1
    counters: ResetCounters = dataclasses.field(default_factory=ResetCounters)
    hypervisc_filter: Optional[jnp.ndarray] = None

    def replace(self, **kw) -> "State":
        return dataclasses.replace(self, **kw)


class Model:
    """Static data and parameters of every kernel, built once from the export."""

    def __init__(self, ex: Export, *, check: bool = True):
        if check:
            check_gates(ex)
        self.ex = ex
        self.dims = Dims.from_export(ex)
        self.reinit = ReinitParams.from_export(ex)
        # dt-dependent coefficient builders
        self.coeff_inputs = coeffs_mod.CoeffInputs.from_export(ex)
        self.coeff_params = coeffs_mod.CoeffParams.from_export(ex)
        # source, parallel solve (+ linked BCs, parity)
        self.src_params = source_mod.SourceParams.from_export(ex)
        self.src_data = to_jax(source_mod.SourceData.from_export(ex))
        self.inv_params = invert_mod.InvertParams.from_export(ex)
        self.inv_data = to_jax(invert_mod.InvertData.from_export(ex))
        # fields
        self.fld_params = fields_mod.FieldParams.from_export(ex)
        self.fld_data = to_jax(fields_mod.FieldData.from_export(ex, self.fld_params))
        local = self.fld_params.field_option == "local"
        self.run = RunParams.from_export(ex, local_fields=local)
        # explicit / nonlinear terms
        self.nl_params = nl_mod.NLParams.from_export(ex)
        self.nl_data = to_jax(nl_mod.NLData.from_export(ex)) if self.run.nonlin else None
        self.x_params = explicit_mod.ExplicitParams.from_export(ex)
        # hyperviscosity
        self.hyp_params = hyper_mod.HyperParams.from_export(ex)
        self.hyp_data = to_jax(hyper_mod.HyperData.from_export(ex, self.hyp_params))
        # collisions
        if self.run.colls:
            self.coll_params = coll_init.CollParams.from_export(ex)
            self.coll_host = coll_init.CollData.from_export(ex, self.coll_params)
            self.coll_data = to_jax(self.coll_host.device())
        else:
            self.coll_params = self.coll_host = self.coll_data = None
        # split_collisions with inactive collisions still re-solves the fields from h
        # (see _refield_from_h); it needs g_adjust's data.
        self.gadj = None
        if self.run.split_collisions and not self.run.colls:
            self.gadj = to_jax(GAdjustData.from_export(ex))
        # flow shear, zonal flows
        if self.run.exb_active:
            self.exb_params = exb_mod.ExbParams.from_export(ex)
            self.exb_data = exb_mod.ExbData.from_export(ex)
        else:
            self.exb_params = self.exb_data = None
        self.zon_params = zonal_mod.ZonalParams.from_export(ex)
        self.zon_data = to_jax(zonal_mod.ZonalData.from_export(ex))

    # ------------------------------------------------------------------ dt set
    def dtset_from_export(self) -> DtSet:
        """The dt-dependent arrays exactly as GS2 built them at the initial dt."""
        ex = self.ex
        c = Coeffs.from_export(ex)
        if c.g_h is None:
            c = c._replace(g_h=invert_mod.init_homogeneous_g(to_jax(c), self.inv_data,
                                                              self.inv_params))
        cc = coll_init.CollCoeffs.from_export(ex) if self.run.colls else None
        if self.run.local_fields:
            aminv = fields_local_mod.local_aminv_from_export(ex)
        else:
            aminv = response_mod.Aminv.from_export(ex)
        return DtSet(code_dt=c.code_dt, coeffs=to_jax(c), coll=to_jax(cc) if cc else None,
                     aminv=to_jax(aminv))

    def build_dtset(self, code_dt: float, vnmult) -> DtSet:
        """Rebuild every dt-dependent array for ``code_dt`` (spec 10 section E.3-E.5)."""
        c = coeffs_mod.compute_coeffs(code_dt, self.coeff_inputs, self.coeff_params)
        c = to_jax(c)
        c = c._replace(g_h=invert_mod.init_homogeneous_g(c, self.inv_data, self.inv_params))
        cc = None
        if self.run.colls:
            cc = to_jax(coll_init.compute_coll_coeffs(code_dt, vnmult, self.coll_host,
                                                      self.coll_params))
        partial_set = DtSet(code_dt=code_dt, coeffs=c, coll=cc, aminv=None)
        st, dims = _statics(self), self.dims
        apply_hyper = hyper_mod.hyper_applies(0, self.hyp_params)
        aminv = response_mod.build_aminv(
            lambda fnew, args: _probe(fnew, args[0], args[1], st, dims, apply_hyper),
            self.fld_data, self.fld_params, code_dt=code_dt,
            timeadv_args=(partial_set, _data(self)))
        return DtSet(code_dt=code_dt, coeffs=c, coll=cc, aminv=aminv)

    # ------------------------------------------------------------------- state
    def initial_state(self) -> State:
        """The JAX initial state: ``/state_t0`` (after GS2's initialize_diagnostics)."""
        ex = self.ex
        st = ex["/state_t0"]
        dims = self.dims
        old = to_jax(fields_from_group(st, dims, new=False))
        new = to_jax(fields_from_group(st, dims, new=True))
        exb = exb_mod.ExbState.from_export(ex, self.exb_params) if self.run.exb_active else None
        return State(
            g=jnp.asarray(st["g"]), gnew=jnp.asarray(st["gnew"]), old=old, new=new,
            time=TimeBlock.from_group(st), reset=bool(np.asarray(st["reset"]).item()),
            xstate=explicit_mod.ExplicitState.from_export(ex), exb=exb,
            vnmult=np.asarray(st["vnmult"], dtype=np.float64),
            timeadv_istep_last=0,   # the response matrix was built with istep = 0 (09 section 3.13)
        )


# ---------------------------------------------------------------------------- timeadv core
class _Statics(NamedTuple):
    """Hashable parameters of the jitted step kernels (static under ``jax.jit``)."""

    src: Any
    inv: Any
    hyp: Any
    coll: Any
    fld: Any
    run: RunParams


class _Data(NamedTuple):
    """Device arrays of the jitted step kernels (traced arguments)."""

    src: Any
    inv: Any
    hyp: Any
    coll: Any
    fld: Any


def _statics(model: Model) -> _Statics:
    return _Statics(model.src_params, model.inv_params, model.hyp_params, model.coll_params,
                    model.fld_params, model.run)


def _data(model: Model) -> _Data:
    return _Data(model.src_data, model.inv_data, model.hyp_data, model.coll_data,
                 model.fld_data)


def _pad_ab(ab) -> jnp.ndarray:
    """AB coefficients as a fixed-length (3,) array; zero padding is exact (source.py)."""
    out = np.zeros(3, np.float64)
    out[:len(ab)] = ab
    return jnp.asarray(out)


@partial(jax.jit, static_argnames=("st", "explicit_on", "apply_hyper", "mode"))
def _invert_and_vspace(g, old: Fields, new: Fields, gnew_in, gexp, ab, code_time, code_dt,
                       coeffs, coll_coeffs, data: _Data, st: _Statics, explicit_on: bool,
                       apply_hyper: bool, mode):
    """``invert_rhs`` + hyper + collisions + parity of ``timeadv`` (dist_fn.fpp:4068-4091).

    ``explicit_on`` is ``istep > 0`` (the only property of ``istep`` the source uses);
    ``apply_hyper`` is ``hyper.hyper_applies(istep, ...)``."""
    istep = 1 if explicit_on else 0
    src = source_mod.get_source_term_on_local_domain(
        g, old, new, gexp, ab, istep, coeffs, data.src, st.src, code_time=code_time)
    gnew = invert_mod.invert_rhs(src, new.phi, new.bpar, coeffs, data.inv, st.inv,
                                 gnew_in=gnew_in)
    hfilter = None
    if apply_hyper:
        gnew, hfilter = hyper_mod.hyper_diff(gnew, new.phi, code_dt, data.hyp, st.hyp)
    if st.run.colls and not st.run.split_collisions:
        gnew = coll.vspace_derivatives(gnew, g, old, new, new.apar, mode, coll_coeffs,
                                       data.coll, st.coll)
    if st.run.def_parity:
        gnew = parity_mod.enforce_parity(gnew, data.inv, st.inv)
    return gnew, hfilter


@partial(jax.jit, static_argnames=("fp", "local"))
def _getfield_jit(new: Fields, old: Fields, gnew, aminv, fld_data, fp, local: bool) -> Fields:
    """The field solve and update of advance_implicit (fields_implicit.f90:352-356) or
    advance_local (fields_local.f90:4313-4320)."""
    if local:
        u = fields_local_mod.getfield_local(new, gnew, aminv, fld_data, fp)
        return fields_local_mod.field_update_local(u, old, fp)
    u = fields_mod.getfield(new, gnew, aminv, fld_data, fp)
    return fields_mod.field_update(u, old, fp)


def _probe(fields_new: Fields, dts: DtSet, data: _Data, st: _Statics, dims: Dims,
           apply_hyper: bool):
    """``timeadv(phi=0, apar=0, bpar=0, phinew, aparnew, bparnew, istep=0)`` from ``g = 0``:
    one response-matrix probe (fields_implicit.f90:650-719; spec 10 section C, G-27/G-28).

    Pure JAX with every array an argument, so ``response.build_aminv`` can vmap it over
    probes without baking the model's arrays into the compiled program as constants."""
    zero_g = jnp.zeros(dims.g_shape, jnp.complex128)
    zero_f = jnp.zeros(dims.field_shape, jnp.complex128)
    old = Fields(zero_f, zero_f, zero_f)
    gnew, _ = _invert_and_vspace(
        zero_g, old, fields_new, zero_g, None, _pad_ab(()), 0.0, dts.code_dt, dts.coeffs,
        dts.coll, data, st, explicit_on=False, apply_hyper=apply_hyper, mode=None)
    return gnew


def timeadv_probe(model: Model, dts: DtSet, fields_new: Fields):
    """One response-matrix probe for ``model`` (see :func:`_probe`)."""
    return _probe(fields_new, dts, _data(model), _statics(model), model.dims,
                  hyper_mod.hyper_applies(0, model.hyp_params))


def timeadv(model: Model, dts: DtSet, S: State, istep: int, mode) -> tuple[State, bool]:
    """``timeadv`` (dist_fn.fpp:4014-4116) with ``phi.. = S.old``, ``phinew.. = S.new``.

    Returns ``(state, aborted)``; ``aborted`` is ``reset .and. immediate_reset`` at
    dist_fn.fpp:4066 (``timeadv::istep_last`` is then not updated)."""
    xo = explicit_mod.add_explicit_terms(S.xstate, S.g, S.old, istep, S.time, S.reset,
                                         model.nl_data, model.nl_params, model.x_params)
    S = S.replace(xstate=xo.state, time=xo.time, reset=xo.reset)
    if xo.abort_timeadv:
        return S, True
    gexp = None
    if model.run.nonlin:
        gexp = (S.xstate.gexp_1, S.xstate.gexp_2, S.xstate.gexp_3)
    ab = source_mod.ab_coefficients(S.time, istep) if istep > 0 else ()
    gnew, hfilter = _invert_and_vspace(
        S.g, S.old, S.new, S.gnew, gexp, _pad_ab(ab), S.time.code_time, S.time.code_dt,
        dts.coeffs, dts.coll, _data(model), _statics(model), explicit_on=istep > 0,
        apply_hyper=hyper_mod.hyper_applies(istep, model.hyp_params), mode=mode)
    S = S.replace(gnew=gnew, timeadv_istep_last=istep,
                  hypervisc_filter=hfilter if hfilter is not None else S.hypervisc_filter)
    return S, False


def _getfield(model: Model, dts: DtSet, S: State) -> Fields:
    """The field solve and update of advance_implicit (fields_implicit.f90:352-356) or
    advance_local (fields_local.f90:4313-4320)."""
    return _getfield_jit(S.new, S.old, S.gnew, dts.aminv, model.fld_data, model.fld_params,
                         model.run.local_fields)


@partial(jax.jit, static_argnames=("fphi", "fbpar", "fp"))
def _refield_from_h(gnew, new: Fields, gadj, fld_data, fphi: float, fbpar: float, fp):
    """The field re-solve of ``collisions_advance`` without collisions
    (dist_fn.fpp:4161-4185): ``h = g_adjust(gnew, from_g)``; fields from
    ``calculate_potentials_from_nonadiabatic_dfn(h)``; ``gnew = g_adjust(h, to_g)``
    with the new fields."""
    h = g_adjust_direction(gnew, new.phi, new.bpar, FROM_G, fphi, fbpar, gadj)
    fields = fields_mod.calculate_potentials_from_nonadiabatic_dfn(h, fld_data, fp)
    if not isinstance(fields, Fields):
        fields = Fields(*fields)
    return g_adjust_direction(h, fields.phi, fields.bpar, TO_G, fphi, fbpar, gadj), fields


def potentials_from_h(model: Model):
    """``calculate_potentials_from_nonadiabatic_dfn`` as a callable for collisions_advance."""
    return lambda h: fields_mod.calculate_potentials_from_nonadiabatic_dfn(
        h, model.fld_data, model.fld_params)


def advance(model: Model, dts: DtSet, S: State, istep: int) -> tuple[State, bool]:
    """``fields::advance(istep)`` for implicit and local fields (spec 01 section 3.2).

    Returns ``(state, reset)``: ``reset`` means an immediate CFL reset was requested in the
    predictor; the state then has ``gnew == g`` and the caller resets dt and redoes the step."""
    run = model.run
    if run.exb_active:
        gnew, new, exb = exb_mod.exb_shear(S.gnew, S.new, S.exb, istep, S.time,
                                           model.exb_data, model.exb_params)
        S = S.replace(gnew=gnew, new=new, exb=exb)
    S = S.replace(g=S.gnew, old=S.new)
    S, aborted = timeadv(model, dts, S, istep, None)
    if S.reset:
        return S, True
    new = _getfield(model, dts, S)
    if run.remove_zonal_flows:
        new = new._replace(phi=zonal_mod.remove_zonal_flows(new.phi, model.zon_data,
                                                            model.zon_params))
    S = S.replace(new=new)
    S, _ = timeadv(model, dts, S, istep, 1)
    if run.split_collisions:                     # collisions_advance (fields_implicit.f90:369)
        if run.colls:
            gnew, new = coll.collisions_advance(S.gnew, S.g, S.old, S.new,
                                                potentials_from_h(model), dts.coll,
                                                model.coll_data, model.coll_params,
                                                diagnostics=1)
        else:
            # GS2 quirk: collisions_advance only tests split_collisions, so with
            # collisions inactive (colls = F) vspace_derivatives is a no-op but the
            # fields are still re-solved from h, without the fl_avg term
            # (dist_fn.fpp:4161-4185, le_derivatives.fpp:70).
            gnew, new = _refield_from_h(S.gnew, S.new, model.gadj, model.fld_data, run.fphi,
                                        run.fbpar, model.fld_params)
        S = S.replace(gnew=gnew, new=new)
    return S, False


def restore_after_reset(model: Model, S: State, gnew_saved) -> State:
    """State restoration after a dt change (gs2_init.fpp:292-375; spec 01 section 3.9):
    ``g = gnew = saved gnew``; fields recomputed from g (force_maxwell_reinit, default) or
    kept; ``phi = phinew``. gexp_*, vnmult, kx_shift and the counters are unchanged."""
    g = gnew_saved
    if model.run.force_maxwell_reinit:
        new = fields_mod.get_init_field(g, model.fld_data, model.fld_params)
    else:
        new = S.new
    return S.replace(g=g, gnew=g, new=new, old=new)
