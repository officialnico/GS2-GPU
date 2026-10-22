"""Source term (right-hand side ``g_work``) of the implicit GK equation.

Port of GS2 8.2.1 (owner spec ``02_source_term.md``, sections 3.2-3.5 and 3.7; critique
``10_critique.md`` G-14/G-30).  Line numbers refer to commit 4d8c94bcf:

* :func:`get_source_term_on_local_domain` mirrors ``get_source_term_on_local_domain``
  (src/dist_fn.fpp:5363-5420) with ``get_source_term`` (4856-5130) and its contained
  ``set_source`` (5017-5128) inlined: time-decentred potentials, ``set_source``
  (incl. the parallel flow-shear drive and the explicit Adams-Bashforth terms with
  ``gexp_1..3``), ``phiext_full``, the ``b*g + a*g`` "matrix multiplication" and the
  totally-trapped-particle (ttp) overwrite;
* :func:`ab_coefficients` / :func:`adams_bashforth_coefficients` mirror
  ``get_adams_bashforth_coefficients`` (src/gs2_time.f90:250-312) and
  ``time_history_available`` (234-244) on the host;
* :func:`sourcefac` mirrors the ``sourcefac`` block of ``invert_rhs``
  (src/dist_fn.fpp:6166-6176).

``opt_source = .true.`` (``get_source_term_opt``) is GATED, as is
``source_option = 'phiext_full'`` with ``phi_ext /= 0`` (G-14); ``phiext_full`` with
``phi_ext == 0`` is supported (numerically identical to ``full``).

Round-off: the Fortran association order is kept everywhere (fully parenthesised);
mixed real/complex products use ``numerics.rmul/imul/cmul``; real products feeding a sum
are wrapped in ``numerics.ob`` (strict mode = no FMA contraction); the one division
(``.../spec(is)%stm`` in ``set_source``) goes through ``coeffs.exact_div`` because XLA
turns division by a broadcast vector into a multiplication by its reciprocal.  The
``b*g + a*g`` products mirror the FCMLA instructions gfortran emits for that vectorised
loop (:func:`fcmla_mul`).  With these, strict mode reproduces a gfortran build of the
verbatim Fortran expressions with the reference flags bit for bit (checked on synthetic
inputs, see the B2 report).
"""

from __future__ import annotations

import cmath
import math
from typing import NamedTuple, Optional, Sequence

import jax.numpy as jnp
import numpy as np
from jax import lax

from .coeffs import _f64, _pad_theta, compress_axes, exact_div
from .io.export import Export
from .numerics import cmul, imul, ob, rmul
from .types import Coeffs, Fields, TimeBlock

# dist_fn.fpp:130-131
SOURCE_OPTION_FULL = 1
SOURCE_OPTION_PHIEXT_FULL = 5
SOURCE_OPTION_HOMOGENEOUS = 6
_SOURCE_OPTION_NAMES = {"default": SOURCE_OPTION_FULL, "full": SOURCE_OPTION_FULL,
                        "phiext_full": SOURCE_OPTION_PHIEXT_FULL,
                        "homogeneous": SOURCE_OPTION_HOMOGENEOUS}

# gs2_time.f90:19 (dt_not_set) and Fortran epsilon(0.0) under -fdefault-real-8.
DT_NOT_SET = -1.0
EPSILON = float(np.finfo(np.float64).eps)


# ============================================================================ host side
def time_history_available(code_dt_prev1: float, code_dt_prev2: float) -> int:
    """Number of time levels of history; mirrors ``time_history_available``
    (src/gs2_time.f90:234-244): 3 if ``code_dt_prev2`` is set, else 2 if
    ``code_dt_prev1`` is set, else 1 (sentinel ``dt_not_set = -1.0``)."""
    if code_dt_prev2 != DT_NOT_SET:
        return 3
    if code_dt_prev1 != DT_NOT_SET:
        return 2
    return 1


def adams_bashforth_coefficients(code_dt: float, code_dt_prev1: float,
                                 code_dt_prev2: float, maximum_order: Optional[int] = None,
                                 target_dt: Optional[float] = None) -> tuple:
    """Variable-step Adams-Bashforth coefficients (newest first), host float64.

    Literal mirror of ``get_adams_bashforth_coefficients`` (src/gs2_time.f90:250-312;
    spec 01 section 3.10, 02 section 3.3), including the optional ``maximum_order`` and
    ``target_dt`` arguments.  ``x**2`` is ``x*x``; ``/3.`` and ``/2.`` are true
    divisions.  With equal steps the result is *not* exactly (23/12, -4/3, 5/12).
    """
    order = time_history_available(code_dt_prev1, code_dt_prev2)
    if maximum_order is not None:
        order = min(order, int(maximum_order))
    this_dt = float(code_dt if target_dt is None else target_dt)
    p1 = float(code_dt_prev1)
    p2 = float(code_dt_prev2)
    if order == 3:
        # (((this_dt + p1)**2/p1)*((this_dt+p1)/3. + p2/2.) - p1*(p1/3. + p2/2.))
        #   / (this_dt*(p1 + p2))                                       gs2_time.f90:287-291
        q = this_dt + p1
        c1 = ((((q * q) / p1) * ((q / 3.0) + (p2 / 2.0))) - (p1 * ((p1 / 3.0) + (p2 / 2.0)))) \
            / (this_dt * (p1 + p2))
        # - this_dt*(this_dt/3. + (p1+p2)/2.)/(p1*p2)                   :293-295
        c2 = -((this_dt * ((this_dt / 3.0) + ((p1 + p2) / 2.0))) / (p1 * p2))
        # this_dt*(this_dt/3. + p1/2.)/(p2*(p1+p2))                    :297-299
        c3 = (this_dt * ((this_dt / 3.0) + (p1 / 2.0))) / (p2 * (p1 + p2))
        return (c1, c2, c3)
    if order == 2:
        c1 = ((0.5 * this_dt) + p1) / p1                               # :303
        c2 = -((0.5 * this_dt) / p1)                                   # :305
        return (c1, c2)
    return (1.0,)                                                      # :309


def ab_coefficients(time: TimeBlock, istep: Optional[int] = None, *,
                    maximum_order: Optional[int] = None,
                    target_dt: Optional[float] = None) -> tuple:
    """AB coefficients for the source of step ``istep`` (host; spec 01 section 3.10).

    Mirrors ``get_adams_bashforth_coefficients()`` as called by
    ``get_source_term_on_local_domain`` (src/dist_fn.fpp:5388-5389), from the
    ``code_dt, code_dt_prev1, code_dt_prev2`` of ``time``.  GS2 evaluates them only
    when ``istep > 0`` (``include_explicit``, :5382); for ``istep <= 0`` this returns
    ``()`` (no explicit terms, e.g. in the response-matrix probes).  Both ``timeadv``
    calls of a step use the same values.
    """
    if istep is not None and int(istep) <= 0:
        return ()
    return adams_bashforth_coefficients(time.code_dt, time.code_dt_prev1, time.code_dt_prev2,
                                        maximum_order=maximum_order, target_dt=target_dt)


def sourcefac(code_time, params: "SourceParams"):
    """Time factor of the external source; mirrors ``invert_rhs``
    (src/dist_fn.fpp:6166-6176)::

        time > t0:  source0*exp(-zi*omega0*time + gamma0*time)
        else:       (0.5 - 0.5*cos(pi*time/t0))*exp(-zi*omega0*time + gamma0*time)

    ``-zi*omega0*time + gamma0*time`` lowers to ``(gamma0*time, -(omega0*time))``.
    Host Python complex for a Python/NumPy ``code_time`` (``cmath.exp`` may differ from
    gfortran's ``cexp`` in the last bit); a ``jnp`` value for a traced one.  Only used
    multiplied by ``phi_ext``, which must be 0 (G-14).
    """
    t0, omega0, gamma0, source0 = params.t0, params.omega0, params.gamma0, params.source0
    if isinstance(code_time, (int, float, np.floating, np.integer)):
        time = float(code_time)
        e = cmath.exp(complex(gamma0 * time, -(omega0 * time)))
        fac = source0 if time > t0 else (0.5 - (0.5 * math.cos((math.pi * time) / t0)))
        return complex(fac * e.real, fac * e.imag)
    time = jnp.asarray(code_time, jnp.float64)
    mag = jnp.exp(gamma0 * time)
    ph = -(omega0 * time)
    e_re, e_im = mag * jnp.cos(ph), mag * jnp.sin(ph)
    ramp = 0.5 - (0.5 * jnp.cos(exact_div(math.pi * time, _f64(t0))))
    fac = jnp.where(time > t0, source0, ramp)
    return lax.complex(fac * e_re, fac * e_im)


# ============================================================================ containers
def _scalar(group, name: str, default=None):
    if name in group:
        return float(np.asarray(group[name]).reshape(()))
    if default is not None:
        return float(default)
    raise KeyError(f"{group.path}/{name} missing from the export")


def _flag(group, name: str, default: bool = False) -> bool:
    if name in group:
        return bool(np.asarray(group[name]).reshape(()))
    return default


def _source_option(df) -> int:
    if "source_option_switch" in df:
        return int(np.asarray(df["source_option_switch"]).reshape(()))
    name = str(df.attrs.get("source_option", "default")).strip().lower()
    if name not in _SOURCE_OPTION_NAMES:
        raise NotImplementedError(f"source_option = {name!r}")
    return _SOURCE_OPTION_NAMES[name]


class SourceParams(NamedTuple):
    """Static switches of the source (hashable; ``jit`` static argument).

    * ``source_option``: ``source_option_switch`` (1 full, 5 phiext_full, 6 homogeneous)
    * ``nonlin``, ``split_nonlinear``: ``include_explicit = istep > 0 .and. nonlin .and.
      .not. split_nonlinear`` (dist_fn.fpp:5382-5386)
    * ``trapped_grid``: ``nlambda > ng2``; ``has_ttp``: trapped grid with at least one
      ``can_be_ttp(il) .and. is_ttp(ig,il)`` point at ``ig < ntgrid``
    * ``apar_terms`` / ``bpar_terms``: ``has_apar`` / ``has_bpar``.  When false, GS2
      still evaluates the corresponding terms but they are exactly +-0 (the absent
      fields are identically zero), so skipping them changes no value (only the sign
      of some zeros)
    * ``phiext_modes``: ``phiext_full`` and some ``aky < epsilon(0.0)``
    * ``t0, omega0, gamma0, source0``: ``source_knobs`` (only for :func:`sourcefac`)
    """

    source_option: int
    nonlin: bool
    split_nonlinear: bool
    trapped_grid: bool
    has_ttp: bool
    apar_terms: bool
    bpar_terms: bool
    phiext_modes: bool
    t0: float
    omega0: float
    gamma0: float
    source0: float

    @classmethod
    def from_export(cls, ex: Export) -> "SourceParams":
        df = ex["/dist_fn"]
        if _flag(df, "opt_source"):
            raise NotImplementedError("opt_source = .true. (get_source_term_opt) is GATED")
        option = _source_option(df)
        if option not in (SOURCE_OPTION_FULL, SOURCE_OPTION_PHIEXT_FULL,
                          SOURCE_OPTION_HOMOGENEOUS):
            raise NotImplementedError(f"source_option_switch = {option}")
        phi_ext = _scalar(df, "phi_ext", 0.0)
        if option == SOURCE_OPTION_PHIEXT_FULL and phi_ext != 0.0:
            raise NotImplementedError("source_option = 'phiext_full' with phi_ext /= 0 is "
                                      "GATED (10_critique G-14)")
        le = ex["/le_grids"]
        trapped = ex.nlambda > int(le["ng2"])
        ttp = np.asarray(le["can_be_ttp"], bool)[:, None] & np.asarray(le["is_ttp"], bool)
        rp = ex["/run_parameters"]
        nl = ex["/nonlinear_terms"]
        aky = np.asarray(ex["/kt_grids/aky"], np.float64)
        return cls(
            source_option=option,
            nonlin=_flag(nl, "nonlin"),
            split_nonlinear=_flag(nl, "split_nonlinear"),
            trapped_grid=bool(trapped),
            has_ttp=bool(trapped and ttp[:, :-1].any()),
            apar_terms=_flag(rp, "has_apar"),
            bpar_terms=_flag(rp, "has_bpar"),
            phiext_modes=bool(option == SOURCE_OPTION_PHIEXT_FULL and np.any(aky < EPSILON)),
            t0=_scalar(df, "t0", 100.0), omega0=_scalar(df, "omega0", 0.0),
            gamma0=_scalar(df, "gamma0", 0.0), source0=_scalar(df, "source0", 1.0),
        )


class SourceData(NamedTuple):
    """dt-independent arrays read by the source (NumPy; ``to_jax`` before use).

    * ``aj0, aj1``: G6; ``vperp2``: G6 with ``Y, X`` compressed when constant;
      ``vpac``: G7 with ``Y, X`` compressed when constant (see ``coeffs.compress_axes``)
    * ``nonmaxw_corr, ufac``: ``(S, E)``; ``D_res``: K2 ``(Y, X)`` (frozen at the
      initial dt, 02 pitfall 7); ``itor_over_B``: ``(T,)``; ``wunits``: ``(Y,)``
    * species ``(S,)``: ``zstm, stm, tz, fexp, bkdiff``
    * 0-d: ``fphi, fapar, fbpar, omprimfac, g_exb, phi_ext``
    * masks: ``lmask (L,)`` = ``il <= lmax``; ``ttp_mask (L, T)`` = ``can_be_ttp(il)
      .and. is_ttp(ig, il)``; ``kwork_filter`` K2 (value during the time loop);
      ``aky_small (Y,)`` = ``aky < epsilon(0.0)``
    """

    aj0: np.ndarray
    aj1: np.ndarray
    vperp2: np.ndarray
    vpac: np.ndarray
    nonmaxw_corr: np.ndarray
    ufac: np.ndarray
    D_res: np.ndarray
    itor_over_B: np.ndarray
    wunits: np.ndarray
    zstm: np.ndarray
    stm: np.ndarray
    tz: np.ndarray
    fexp: np.ndarray
    bkdiff: np.ndarray
    fphi: np.ndarray
    fapar: np.ndarray
    fbpar: np.ndarray
    omprimfac: np.ndarray
    g_exb: np.ndarray
    phi_ext: np.ndarray
    lmask: np.ndarray
    ttp_mask: np.ndarray
    kwork_filter: np.ndarray
    aky_small: np.ndarray

    @classmethod
    def from_export(cls, ex: Export) -> "SourceData":
        arrays = ex["/dist_fn_arrays"]
        df = ex["/dist_fn"]
        sp = ex["/species"]
        le = ex["/le_grids"]
        rp = ex["/run_parameters"]
        kt = ex["/kt_grids"]
        hyper = ex["/hyper"] if "hyper" in ex.root.groups else None
        # GS2 quirk: D_res contains code_dt but init_hyper is guarded by `initialized`
        # (hyper.f90:370) and level 19 is never re-entered, so D_res stays at its
        # initial-dt value for the whole run, also after dt resets (02 pitfall 7).
        if hyper is not None and "D_res" in hyper:
            d_res = np.asarray(hyper["D_res"], np.float64)
        else:                                   # hyper.f90:409, 615: D_res = 0 when unused
            d_res = np.zeros((ex.naky, ex.ntheta0))
        lmax = int(le["lmax"])
        il = np.arange(1, ex.nlambda + 1)
        aky = np.asarray(kt["aky"], np.float64)
        return cls(
            aj0=np.asarray(arrays["aj0"], np.float64),
            aj1=np.asarray(arrays["aj1"], np.float64),
            vperp2=compress_axes(arrays["vperp2"], (3, 4)),
            vpac=compress_axes(arrays["vpac"], (3, 4)),
            nonmaxw_corr=np.asarray(sp["nonmaxw_corr"], np.float64),
            ufac=np.asarray(df["ufac"], np.float64),
            D_res=d_res,
            itor_over_B=np.asarray(ex["/theta_grid/itor_over_B"], np.float64),
            wunits=np.asarray(ex["/gs2_time/wunits"], np.float64),
            zstm=np.asarray(sp["spec_zstm"], np.float64),
            stm=np.asarray(sp["spec_stm"], np.float64),
            tz=np.asarray(sp["spec_tz"], np.float64),
            fexp=np.atleast_1d(np.asarray(df["fexp"], np.float64)),
            bkdiff=np.atleast_1d(np.asarray(df["bkdiff"], np.float64)),
            fphi=_f64(_scalar(rp, "fphi")), fapar=_f64(_scalar(rp, "fapar")),
            fbpar=_f64(_scalar(rp, "fbpar")),
            omprimfac=_f64(_scalar(df, "omprimfac")), g_exb=_f64(_scalar(df, "g_exb")),
            phi_ext=_f64(_scalar(df, "phi_ext", 0.0)),
            lmask=il <= lmax,
            ttp_mask=np.asarray(le["can_be_ttp"], bool)[:, None] & np.asarray(le["is_ttp"], bool),
            kwork_filter=np.asarray(kt["kwork_filter"], bool),
            aky_small=aky < EPSILON,
        )


# ============================================================================== kernel
def include_explicit(istep, ab, params: SourceParams):
    """``include_explicit`` of dist_fn.fpp:5382-5386: ``istep > 0 .and. nonlin .and.
    .not. split_nonlinear``.  A Python bool for a Python ``istep`` (then a nonlinear
    ``istep > 0`` without AB coefficients is an error); for a traced ``istep`` a traced
    bool (``False`` if no coefficients were passed)."""
    if not params.nonlin or params.split_nonlinear:
        return False
    no_ab = ab is None or len(ab) == 0
    if isinstance(istep, (int, np.integer)):
        if istep <= 0:
            return False
        if no_ab:
            raise ValueError(f"nonlinear source at istep={istep} needs the AB coefficients "
                             "(ab_coefficients(time, istep))")
        return True
    if no_ab:
        return False
    return jnp.asarray(istep) > 0


def _add_if(value, term, flag):
    """``value + term`` where the explicit terms are included (static or traced flag)."""
    if flag is True:
        return value + term
    return jnp.where(flag, value + term, value)


def explicit_sum(gexp, ab):
    """``ab(1)*gexp_1 [+ ab(2)*gexp_2 [+ ab(3)*gexp_3]]`` on cells ``ig < ntgrid``, G7
    ``(..., 2, N)``; the ``select case (size(ab_coefficients))`` of set_source
    (src/dist_fn.fpp:5106-5127) and of the ttp block (4986-5010), left to right."""
    k = len(ab)
    if k not in (1, 2, 3):
        raise ValueError(f"expected 1-3 Adams-Bashforth coefficients, got {k}")
    total = rmul(jnp.asarray(gexp[0])[..., :-1], ab[0])
    if k >= 2:
        total = total + rmul(jnp.asarray(gexp[1])[..., :-1], ab[1])
    if k >= 3:
        total = total + rmul(jnp.asarray(gexp[2])[..., :-1], ab[2])
    return total


def fcmla_mul(c, x):
    """Fortran ``c*x`` (complex*complex) as gfortran vectorises it for ``-mcpu=apple-m1``.

    In vectorised loops gfortran -O3 emits the ARMv8.3 ``FCMLA #0, #90`` pair for a
    complex product even with ``-ffp-contract=off``: the ``c.re`` products are rounded,
    the ``c.im`` products are fused into them::

        re = fma(-c.im, x.im, c.re*x.re) ;  im = fma(c.im, x.re, c.re*x.im)

    with ``c`` the left operand as written in the Fortran.  Verified bit for bit against
    the verbatim matrix-term loop of get_source_term (dist_fn.fpp:4932-4944) compiled
    with the reference flags (``objs/dist_fn.o`` of the reference build contains 48
    ``fcmla``), and by the strict-mode trace tests (all ``source`` probes of all cases
    bitwise); scalar (non-vectorised) loops such as the ttp block use the naive
    product (``numerics.cmul``).  In strict mode the ``c.re`` products sit behind a barrier and
    XLA:CPU fuses the remaining ``t - c.im*x.im`` into an FMA, which reproduces FCMLA
    exactly (tested); in default mode XLA contracts as it likes (round-off level).
    """
    t_re = ob(c.real * x.real)
    t_im = ob(c.real * x.imag)
    return lax.complex(t_re - c.imag * x.imag, t_im + c.imag * x.real)


def _f6(field):
    """F3 ``(Y, X, T)`` -> G6 broadcast ``(1, 1, 1, Y, X, T)``."""
    return jnp.asarray(field)[None, None, None]


def _lo(x):
    """G6 grid values -> cell values at ``ig`` with a length-1 sign axis."""
    return x[..., None, :-1]


def _hi(x):
    """G6 grid values -> cell values at ``ig + 1`` with a length-1 sign axis."""
    return x[..., None, 1:]


def _s7(v):
    return v[:, None, None, None, None, None, None]


def gyro_potentials(fields_old: Fields, fields_new: Fields, data: SourceData,
                    params: SourceParams):
    """``(phigavg, apargavg)`` on grid points, G6 complex (``apargavg`` None when
    ``not apar_terms``); mirrors src/dist_fn.fpp:4891-4896::

        phigavg  = (fexp*phi + (1-fexp)*phinew)*aj0*fphi
                 + (fexp*bpar + (1-fexp)*bparnew)*aj1*fbpar*2.0*vperp2*tz
        apargavg = (fexp*apar + (1-fexp)*aparnew)*aj0*fapar
    """
    fe = data.fexp[:, None, None, None, None, None]
    ofe = 1.0 - fe
    phi_c = rmul(_f6(fields_old.phi), fe) + rmul(_f6(fields_new.phi), ofe)
    phigavg = rmul(rmul(phi_c, data.aj0), data.fphi)
    if params.bpar_terms:
        bpar_c = rmul(_f6(fields_old.bpar), fe) + rmul(_f6(fields_new.bpar), ofe)
        bpart = rmul(rmul(rmul(rmul(rmul(bpar_c, data.aj1), data.fbpar), 2.0), data.vperp2),
                     data.tz[:, None, None, None, None, None])
        phigavg = phigavg + bpart
    apargavg = None
    if params.apar_terms:
        apar_c = rmul(_f6(fields_old.apar), fe) + rmul(_f6(fields_new.apar), ofe)
        apargavg = rmul(rmul(apar_c, data.aj0), data.fapar)
    return phigavg, apargavg


def set_source(phigavg, apargavg, fields_old: Fields, fields_new: Fields, coeffs: Coeffs,
               data: SourceData, params: SourceParams):
    """``set_source`` without the explicit terms (src/dist_fn.fpp:5017-5103), on cells
    ``ig < ntgrid``, G7 complex ``(S, E, L, Y, X, 2, N)``::

        source = (-2.0*vpar*nm*phi_m
                  - zstm*vpac*nm*((aj0(ig+1) + aj0(ig))*0.5*apar_m + D_res*apar_p)
                  - zi*wdrift*nm*phi_p)
                 + zi*(wstar + vpac*code_dt*wunits*ufac
                       - 2.0*omprimfac*vpac*nm*code_dt*wunits*g_exb*itor_over_B/stm)
                   *(phi_p - apar_p*stm*vpac)

    with ``phi_p = bdfac_p*phigavg(ig+1) + bdfac_m*phigavg(ig)``, ``phi_m =
    phigavg(ig+1) - phigavg(ig)``, ``apar_p`` likewise and ``apar_m = bdfac_p*(aparnew -
    apar)(ig+1)*fapar + bdfac_m*(aparnew - apar)(ig)*fapar``; ``bdfac_p = 1+bd*(3-2*isgn)``.
    The ``g_exb`` (parallel flow shear) drive uses the raw ``g_exb`` (02 section 3.8).
    """
    dt = jnp.asarray(coeffs.code_dt, jnp.float64)
    nm7 = data.nonmaxw_corr[:, :, None, None, None, None, None]
    sgn = 3.0 - 2.0 * jnp.asarray([1.0, 2.0])                      # 3.-2.*real(isgn)
    bdfac_p = (1.0 + data.bkdiff[:, None] * sgn)[:, None, None, None, None, :, None]
    bdfac_m = (1.0 - data.bkdiff[:, None] * sgn)[:, None, None, None, None, :, None]
    phi_p = rmul(_hi(phigavg), bdfac_p) + rmul(_lo(phigavg), bdfac_m)
    phi_m = _hi(phigavg) - _lo(phigavg)
    vpar = coeffs.vpar[..., :-1]
    vpac = data.vpac[..., :-1]
    wdrift = coeffs.wdrift[..., :-1]
    term1 = -rmul(phi_m, (2.0 * vpar) * nm7)
    term3 = imul(wdrift * nm7, phi_p)
    if apargavg is not None:
        apar_p = rmul(_hi(apargavg), bdfac_p) + rmul(_lo(apargavg), bdfac_m)
        dapar = jnp.asarray(fields_new.apar) - jnp.asarray(fields_old.apar)
        dapar = dapar[None, None, None, :, :, None, :]
        apar_m = rmul(rmul(dapar[..., 1:], bdfac_p), data.fapar) \
            + rmul(rmul(dapar[..., :-1], bdfac_m), data.fapar)
        aj0_c = (_hi(data.aj0) + _lo(data.aj0)) * 0.5
        d_res = data.D_res[None, None, None, :, :, None, None]
        xa = rmul(apar_m, aj0_c) + rmul(apar_p, d_res)
        term2 = rmul(xa, (_s7(data.zstm) * vpac) * nm7)
        src = (term1 - term2) - term3
        chi = phi_p - rmul(rmul(apar_p, _s7(data.stm)), vpac)
    else:
        src = term1 - term3
        chi = phi_p
    wunits = data.wunits[None, None, None, :, None, None, None]
    wstar = coeffs.wstar[:, :, None, :, None, None, None]
    ufac = data.ufac[:, :, None, None, None, None, None]
    # GS2 quirk: the parallel flow-shear drive uses the raw g_exb, neither multiplied by
    # g_exbfac nor gated by g_exb_start_time/timestep (dist_fn.fpp:5101; 02 section 3.8).
    shear = (((((((2.0 * data.omprimfac) * vpac) * nm7) * dt) * wunits) * data.g_exb)
             * data.itor_over_B[:-1])
    w = (wstar + ob(((vpac * dt) * wunits) * ufac)) - exact_div(shear, _s7(data.stm))
    return src + imul(w, chi)


def get_source_term_on_local_domain(g, fields_old: Fields, fields_new: Fields,
                                    gexp: Optional[Sequence], ab: Sequence, istep,
                                    coeffs: Coeffs, data: SourceData, params: SourceParams,
                                    *, code_time=None, aparnew_module=None):
    """The source ``g_work`` (G7 complex ``(S, E, L, Y, X, 2, T)``) of one ``invert_rhs``.

    Mirrors ``get_source_term_on_local_domain`` (src/dist_fn.fpp:5363-5420), i.e.
    ``get_source_term`` (4856-5130) for both signs with ``opt_source = .false.``.

    * ``g``: G7, the module ``g`` (g^n = the step's starting ``gnew``; 0 in probes);
    * ``fields_old``: ``phi, apar, bpar`` (t^n); ``fields_new``: ``phinew, aparnew,
      bparnew`` (pass 1: equal to ``fields_old``; pass 2: after the field solve);
    * ``gexp``: ``(gexp_1, gexp_2, gexp_3)`` G7 (``nl_center``-ed) or None (linear);
    * ``ab``: AB coefficients from :func:`ab_coefficients` (length 1-3, or ``()``);
      padding with zeros is exact up to signed zeros;
    * ``istep``: GS2 step index (0 in response probes); Python int or traced;
    * ``coeffs``: dt-dependent ``vpar, wdrift, wdriftttp, wstar, a, b`` and ``code_dt``;
    * ``code_time``: ``gs2_time::code_time`` (only used by :func:`sourcefac`);
    * ``aparnew_module``: ignored (interface symmetry): inside ``set_source`` the name
      ``aparnew`` is the *dummy* argument of ``get_source_term``, i.e.
      ``fields_new.apar``.

    Per element on cells ``ig < ntgrid`` (02 sections 3.4-3.5)::

        source = ((set_source [+ explicit] [- phiext]) + coef1*g(ig)) + coef2*g(ig+1)

    with ``(coef1, coef2) = (b, a)`` for ``isgn = 1`` and ``(a, b)`` for ``isgn = 2``;
    ``set_source`` only where ``il <= lmax`` (else 0); ``homogeneous`` keeps only the
    matrix terms; the matrix terms only where ``il <= lmax``.  At ttp points
    (``isgn = 2``, full/phiext_full) the value is overwritten with
    ``(g*a - zi*wdriftttp*nm*phigavg) + zi*wstar*phigavg [- phiext] [+ explicit]``.
    Entries ``ig = ntgrid`` and every ``kwork_filter`` mode are 0.
    """
    del aparnew_module
    option = params.source_option
    if option not in (SOURCE_OPTION_FULL, SOURCE_OPTION_PHIEXT_FULL, SOURCE_OPTION_HOMOGENEOUS):
        raise NotImplementedError(f"source_option_switch = {option}")
    explicit = include_explicit(istep, ab, params)
    if explicit is not False and gexp is None:
        raise ValueError("include_explicit needs gexp = (gexp_1, gexp_2, gexp_3)")
    g = jnp.asarray(g)
    a = coeffs.a
    b = coeffs.b
    half_dt = 0.5 * jnp.asarray(coeffs.code_dt, jnp.float64)
    lmask = data.lmask[None, None, :, None, None, None, None]
    driven = option in (SOURCE_OPTION_FULL, SOURCE_OPTION_PHIEXT_FULL)
    phiext = option == SOURCE_OPTION_PHIEXT_FULL and params.phiext_modes
    if phiext:
        if code_time is None:
            raise ValueError("phiext_full needs code_time (sourcefac)")
        sfac = sourcefac(code_time, params)

    # ---- select case (source_option_switch), dist_fn.fpp:4899-4927 ------------------
    if driven:
        phigavg, apargavg = gyro_potentials(fields_old, fields_new, data, params)
        src = set_source(phigavg, apargavg, fields_old, fields_new, coeffs, data, params)
        if explicit is not False:                                        # :5106-5127
            expl = explicit_sum(gexp, ab)
            src = _add_if(src, rmul(expl, half_dt), explicit)
        if phiext:                                                       # :4917-4920
            wdrift = coeffs.wdrift[..., :-1]
            nm7 = data.nonmaxw_corr[:, :, None, None, None, None, None]
            ext = rmul(imul(((wdrift * nm7) * 2.0) * data.phi_ext, sfac), _lo(data.aj0))
            small = data.aky_small[None, None, None, :, None, None, None]
            src = jnp.where(small, src - ext, src)
        body = jnp.where(lmask, src, 0.0 + 0.0j)
    else:                                                                # homogeneous
        shape = g.shape[:-1] + (g.shape[-1] - 1,)
        body = jnp.zeros(shape, jnp.complex128)

    # ---- matrix multiplication with g^n, il <= lmax (dist_fn.fpp:4932-4944) ----------
    # (vectorised by gfortran with FCMLA complex products, see fcmla_mul)
    coef1 = jnp.concatenate([b[..., 0:1, :-1], a[..., 1:2, :-1]], axis=-2)
    coef2 = jnp.concatenate([a[..., 0:1, :-1], b[..., 1:2, :-1]], axis=-2)
    out = (body + fcmla_mul(coef1, g[..., :-1])) + fcmla_mul(coef2, g[..., 1:])
    out = jnp.where(lmask, out, body)

    # ---- totally trapped particles, isgn = 2 (dist_fn.fpp:4957-5012) -----------------
    if params.has_ttp and driven:
        nm6 = data.nonmaxw_corr[:, :, None, None, None, None]
        pg = phigavg[..., :-1]
        wdttp = coeffs.wdriftttp[..., 1, :-1]
        wstar6 = coeffs.wstar[:, :, None, :, None, None]
        ttp_val = (cmul(g[..., 1, :-1], a[..., 1, :-1]) - imul(wdttp * nm6, pg)) \
            + imul(wstar6, pg)
        if phiext:
            ext = rmul(imul(((wdttp * nm6) * 2.0) * data.phi_ext, sfac), data.aj0[..., :-1])
            small = data.aky_small[None, None, None, :, None, None]
            ttp_val = jnp.where(small, ttp_val - ext, ttp_val)
        if explicit is not False:
            # GS2 quirk: nl_center gives ~2x the cell value but the unmodified grid value
            # at ttp points; both get the same 0.5*code_dt (02 pitfall 10).  The set_source
            # addition was overwritten above, so the explicit terms enter exactly once.
            ttp_val = _add_if(ttp_val, rmul(expl[..., 1, :], half_dt), explicit)
        ttp = data.ttp_mask[None, None, :, None, None, :-1]
        out = jnp.stack([out[..., 0, :], jnp.where(ttp, ttp_val, out[..., 1, :])], axis=-2)

    # ---- kwork_filter modes are skipped: g_work stays 0 (dist_fn.fpp:5401) ------------
    kwork = data.kwork_filter[None, None, None, :, :, None, None]
    out = jnp.where(kwork, 0.0 + 0.0j, out)
    return _pad_theta(out)
