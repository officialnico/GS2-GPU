"""Explicit (nonlinear) history ``gexp_1..3`` and the CFL / time-step control signals.

Mirrors (GS2 8.2.1, commit 4d8c94bcf):

* ``nonlinear_terms::add_explicit_terms`` / ``add_explicit``   src/nonlinear_terms.fpp:407-612
* ``nonlinear_terms::estimate_error``                          src/nonlinear_terms.fpp:620-672
* ``get_max_vel``, ``set_cfl_time_step_and_check_if_too_large``,
  ``nb_check_time_step_too_large``                             src/nonlinear_terms.fpp:1191-1252
* ``gs2_time::{time_history_available, get_adams_bashforth_coefficients,
  check_time_step_too_large, save_dt_cfl}``                    src/gs2_time.f90:154-183, 234-312
* ``gs2_reinit::check_time_step`` (end-of-step test)           src/gs2_reinit.f90:228-254

Specs: 05 sections 4.1, 4.9-4.11; 02 section 3.7; 01 sections 1.14, 3.1, 3.5.

The array work (the NL term, the order-based error reductions) runs on the device; the
bookkeeping (istep_last, the max-velocity components, ``code_dt_cfl``, ``reset``) is host
Python float64 mirroring the Fortran literally (00_decisions D9). Reading the CFL velocities
back to the host once per step is inherent: GS2 makes a host-level branch on them.

Control flow of one step with ``immediate_reset = .true.`` (default)::

    out = add_explicit_terms(state, g, fields_old, istep, time, reset, ...)   # pass 1
    if out.abort_timeadv: reset_time_step(...), then redo the step: add_explicit_terms
        is called again with out.state and returns it unchanged (NL not recomputed,
        history not reshuffled: istep == istep_last).
    out2 = add_explicit_terms(out.state, ..., istep, ...)                       # pass 2: no-op

With ``immediate_reset = .false.`` the CFL check is deferred to ``check_time_step`` at the end of
the step (nonblocking reduction in GS2), which also performs the dt-increase test.
"""

from __future__ import annotations

from typing import NamedTuple, Optional

import jax
import jax.numpy as jnp
import numpy as np

from . import nonlinear as nlm
from .io.export import Export
from .numerics import rmul
from .types import Fields, TimeBlock

#: ``epsilon(0.0)`` (nonlinear_terms.fpp:518).
EPS = nlm.EPS
#: ``dt_cfl_default_large`` (nonlinear_terms.fpp:90).
DT_CFL_DEFAULT_LARGE = 1.0e8
#: ``gs2_time::dt_not_set`` (gs2_time.f90:19).
DT_NOT_SET = -1.0


def _item(value):
    return np.asarray(value).item()


def _flag(value) -> bool:
    return bool(np.asarray(value).item())


def _f64(x) -> np.float64:
    return np.float64(x)


# --------------------------------------------------------------------------------------------
# containers
# --------------------------------------------------------------------------------------------


class ExplicitParams(NamedTuple):
    """Knobs of the explicit-term bookkeeping and the CFL control (host scalars)."""

    nonlin: bool
    split_nonlinear: bool
    zip: bool
    use_cfl_limit: bool
    use_2d_cfl: bool
    use_order_based_error: bool
    error_target: float
    istep_error_start: int
    cfl: float
    immediate_reset: bool
    delt_adj: float
    delt_cushion: float

    @classmethod
    def from_export(cls, ex: Export) -> "ExplicitParams":
        """From /nonlinear_terms (09 section 3.12), /run_parameters and /gs2_reinit (3.5).
        Every value is required: a missing knob is an exporter bug, not a default."""
        nt = ex["/nonlinear_terms"]
        reinit = ex["/gs2_reinit"]
        return cls(
            nonlin=_flag(nt["nonlin"]),
            split_nonlinear=_flag(nt["split_nonlinear"]),
            zip=_flag(nt["zip"]),
            use_cfl_limit=_flag(nt["use_cfl_limit"]),
            use_2d_cfl=_flag(nt["use_2d_cfl"]),
            use_order_based_error=_flag(nt["use_order_based_error"]),
            error_target=float(_item(nt["error_target"])),
            istep_error_start=int(_item(nt["istep_error_start"])),
            cfl=float(_item(nt["cfl"])),
            immediate_reset=_flag(ex["/run_parameters/immediate_reset"]),
            delt_adj=float(_item(reinit["delt_adj"])),
            delt_cushion=float(_item(reinit["delt_cushion"])),
        )


class ExplicitState(NamedTuple):
    """State carried between calls (survives dt resets, spec 01 section 1.12).

    ``gexp_1..3``: G7 complex device arrays (None for linear runs: GS2 does not allocate them,
    dist_fn.fpp:3970-3980); ``istep_last``: ``nonlinear_terms::istep_last`` (init 0);
    ``max_vel_components``: the module array (init 0, :323), read again by the deferred
    check; ``cfl_pending``: a nonblocking reduction was started (``cfl_req_hand /= null``).
    """

    gexp_1: Optional[jnp.ndarray]
    gexp_2: Optional[jnp.ndarray]
    gexp_3: Optional[jnp.ndarray]
    istep_last: int = 0
    max_vel_components: tuple = (0.0, 0.0, 0.0)
    cfl_pending: bool = False

    @classmethod
    def from_group(cls, grp, max_vel_components=(0.0, 0.0, 0.0)) -> "ExplicitState":
        """From a group holding ``gexp_1..3`` and the time block's ``nl_istep_last``
        (``/state_t0`` or a trace record such as ``advance_begin``, 09a section 5.2). The
        ``gexp_*`` are None when GS2 did not allocate them (linear runs: listed as absent)."""
        gexp = []
        for k in (1, 2, 3):
            value = grp.get(f"gexp_{k}")
            gexp.append(None if value is None else jnp.asarray(value))
        return cls(gexp_1=gexp[0], gexp_2=gexp[1], gexp_3=gexp[2],
                   istep_last=int(_item(grp["nl_istep_last"])),
                   max_vel_components=tuple(float(v) for v in max_vel_components),
                   cfl_pending=False)

    @classmethod
    def from_export(cls, ex: Export) -> "ExplicitState":
        """The t = 0 state: ``/state_t0`` (10_critique G-05) with the module
        ``/nonlinear_terms/max_vel_components`` (0 after init, nonlinear_terms.fpp:323)."""
        st = ex["/state_t0"]
        if _flag(st["nl_initializing"]):
            # add_explicit's `initializing` branch (nonlinear_terms.fpp:466-476) is not ported:
            # it is dead once fields::init_fields has called nl_finish_init (spec 05 section 2.1).
            raise ValueError("export has nonlinear_terms::initializing = .true. at t = 0")
        mvc = np.asarray(ex["/nonlinear_terms/max_vel_components"], np.float64)
        return cls.from_group(st, max_vel_components=mvc)


class ExplicitOutcome(NamedTuple):
    """Result of one ``add_explicit_terms`` call.

    ``time``: the time block with ``code_dt_cfl``/``user_dt_cfl`` updated when
    ``save_dt_cfl`` ran; ``reset``: ``run_parameters::reset`` after the call (unchanged unless
    ``check_time_step_too_large`` ran); ``computed``: the NL term was evaluated (first call of
    a step); ``returned_early``: ``add_explicit`` returned at :585-587; ``abort_timeadv``:
    ``reset .and. immediate_reset`` -- timeadv returns before invert_rhs (dist_fn.fpp:4066).
    """

    state: ExplicitState
    time: TimeBlock
    reset: bool
    computed: bool
    returned_early: bool
    abort_timeadv: bool


# --------------------------------------------------------------------------------------------
# gs2_time helpers (host)
# --------------------------------------------------------------------------------------------


def time_history_available(time: TimeBlock) -> int:
    """``gs2_time::time_history_available`` (gs2_time.f90:234-244)."""
    if time.code_dt_prev2 != DT_NOT_SET:
        return 3
    if time.code_dt_prev1 != DT_NOT_SET:
        return 2
    return 1


def ab_coefficients(time: TimeBlock, maximum_order: Optional[int] = None,
                    target_dt: Optional[float] = None) -> tuple:
    """``gs2_time::get_adams_bashforth_coefficients`` (gs2_time.f90:250-312), literal float64.

    Variable-step AB coefficients of order ``min(time_history_available, maximum_order)``;
    ``this_dt = target_dt`` or ``code_dt``. ``(x)**2`` is ``x*x`` in gfortran; ``/3.`` and ``/2.``
    are true divisions. (B2 consumes the same numbers in the source term, spec 02 section 3.3;
    this copy serves ``estimate_error``.)
    """
    order = time_history_available(time)
    if maximum_order is not None:
        order = min(order, maximum_order)
    this_dt = _f64(time.code_dt if target_dt is None else target_dt)
    p1 = _f64(time.code_dt_prev1)
    p2 = _f64(time.code_dt_prev2)
    if order == 3:
        q = this_dt + p1
        c1 = (((q * q) / p1) * ((this_dt + p1) / 3.0 + p2 / 2.0)
              - p1 * (p1 / 3.0 + p2 / 2.0)) / (this_dt * (p1 + p2))
        c2 = -(this_dt * (this_dt / 3.0 + (p1 + p2) / 2.0)) / (p1 * p2)
        c3 = this_dt * (this_dt / 3.0 + p1 / 2.0) / (p2 * (p1 + p2))
        return (float(c1), float(c2), float(c3))
    if order == 2:
        c1 = (0.5 * this_dt + p1) / p1
        c2 = -0.5 * this_dt / p1
        return (float(c1), float(c2))
    return (1.0,)


def get_max_vel(components, use_2d_cfl: bool) -> float:
    """``get_max_vel`` (nonlinear_terms.fpp:1214-1224)."""
    c = [float(v) for v in components]
    max_vel_cfl = c[0] + c[1] if use_2d_cfl else max(c[0], c[1])
    return max(max_vel_cfl, c[2])


def set_cfl_time_step_and_check_if_too_large(components, time: TimeBlock, xparams: ExplicitParams):
    """``set_cfl_time_step_and_check_if_too_large`` (nonlinear_terms.fpp:1229-1252).

    ``dt_cfl = 1e8`` if ``max_vel <= 0`` else ``1/max_vel``; ``save_dt_cfl`` sets
    ``code_dt_cfl = user_dt_cfl = dt_cfl`` (gs2_time.f90:177-183);
    ``reset = code_dt > code_dt_cfl`` (check_time_step_too_large, gs2_time.f90:154-158).
    Returns ``(time, reset)``.
    """
    max_vel = get_max_vel(components, xparams.use_2d_cfl)
    dt_cfl = DT_CFL_DEFAULT_LARGE if max_vel <= 0.0 else 1.0 / max_vel
    time = time._replace(code_dt_cfl=dt_cfl, user_dt_cfl=dt_cfl)
    return time, bool(time.code_dt > time.code_dt_cfl)


# --------------------------------------------------------------------------------------------
# order-based error estimate (device reductions + host scalars)
# --------------------------------------------------------------------------------------------


@jax.jit
def _error_reductions(g1, g2, g3, c3, c2):
    """``max|s3 - s2|`` and ``max|s2|`` of estimate_error (nonlinear_terms.fpp:649-661), with
    ``s3 = (c3(1)*g1 + c3(2)*g2) + c3(3)*g3`` and ``s2 = c2(1)*g1 + c2(2)*g2``."""
    s3 = (rmul(g1, c3[0]) + rmul(g2, c3[1])) + rmul(g3, c3[2])
    s2 = rmul(g1, c2[0]) + rmul(g2, c2[1])
    return jnp.max(jnp.abs(s3 - s2)), jnp.max(jnp.abs(s2))


def estimate_error(g1, g2, g3, time: TimeBlock, xparams: ExplicitParams) -> float:
    """``estimate_error`` (nonlinear_terms.fpp:620-672): AB3 vs AB2 relative difference.

    Without three history levels it returns ``(error_target*cfl)*cfl`` (:634-637).
    """
    if time_history_available(time) != 3:
        return float((_f64(xparams.error_target) * _f64(xparams.cfl)) * _f64(xparams.cfl))
    c3 = jnp.asarray(ab_coefficients(time, maximum_order=3), jnp.float64)
    c2 = jnp.asarray(ab_coefficients(time, maximum_order=2), jnp.float64)
    e1, e2 = _error_reductions(g1, g2, g3, c3, c2)
    with np.errstate(divide="ignore", invalid="ignore"):
        return float(_f64(np.asarray(e1)) / _f64(np.asarray(e2)))


# --------------------------------------------------------------------------------------------
# add_explicit_terms
# --------------------------------------------------------------------------------------------


def _default_kernel(g, fields, data, params):
    return nlm.nl_source_jit(g, fields, data, params)


def add_explicit_terms(state: ExplicitState, g, fields: Fields, istep: int, time: TimeBlock,
                       reset: bool, data, params, xparams: ExplicitParams,
                       *, kernel=None) -> ExplicitOutcome:
    """``add_explicit_terms`` + ``add_explicit`` (nonlinear_terms.fpp:407-612).

    Called once per ``timeadv`` (both passes, and again in the redone step after an
    immediate reset). ``g`` and ``fields`` are the start-of-step ``g`` and ``phi, apar, bpar``.

    * linear: ``save_dt_cfl(1e8)`` on every call (:428-430);
    * nonlinear, ``istep == 0``: nothing (response-matrix probes, :425);
    * nonlinear, ``istep /= istep_last`` (:495): ``istep_last = istep``; history shuffle
      ``gexp_3 <- gexp_2``, ``gexp_2 <- gexp_1`` (:505-506); ``max_vel_components = EPS``
      (:518); ``gexp_1 = nl_center(add_nl(g, fields), bd)`` (:520-523); CFL limit
      (:532-534); order-based error (:540-556); reduction + ``set_cfl...`` now
      (``immediate_reset``) or at the end of the step (:561-583); early return if reset
      (:585-587);
    * otherwise (second pass, redone step): nothing.

    ``zip`` and ``split_nonlinear`` are GATED (``NotImplementedError``). The ``initializing``
    branch (:466-476) is dead while stepping (spec 05 section 2.1) and not reproduced.
    ``data``/``params`` are ``nonlinear.NLData`` (device arrays) / ``nonlinear.NLParams``;
    ``kernel(g, fields, data, params) -> (gexp_1, max_vel_components)`` defaults to
    ``nonlinear.nl_source_jit``.
    """
    reset = bool(reset)
    if not xparams.nonlin:
        time = time._replace(code_dt_cfl=DT_CFL_DEFAULT_LARGE, user_dt_cfl=DT_CFL_DEFAULT_LARGE)
        return ExplicitOutcome(state, time, reset, False, False, reset and xparams.immediate_reset)
    if xparams.zip:
        raise NotImplementedError("nonlinear_terms_knobs: zip = .true. is GATED (10_critique G-08)")
    if xparams.split_nonlinear:
        raise NotImplementedError("nonlinear_terms_knobs: split_nonlinear = .true. is GATED")
    if istep == 0:
        return ExplicitOutcome(state, time, reset, False, False, reset and xparams.immediate_reset)

    computed = False
    if istep != state.istep_last:
        computed = True
        gexp_3 = state.gexp_2                                  # copy(g2, g3)  (:505)
        gexp_2 = state.gexp_1                                  # copy(g1, g2)  (:506)
        kern = _default_kernel if kernel is None else kernel
        gexp_1, mv = kern(g, fields, data, params)             # add_nl + nl_center (:520-523)
        mv = np.asarray(mv, dtype=np.float64)
        mvc = [EPS, EPS, EPS]                                  # :518
        mvc[0], mvc[1] = float(mv[0]), float(mv[1])            # add_nl :913-914
        if not xparams.use_cfl_limit:                          # :532-534
            mvc[0] = mvc[1] = EPS
        if xparams.use_order_based_error and istep > xparams.istep_error_start:   # :540-556
            error = estimate_error(gexp_1, gexp_2, gexp_3, time, xparams)
            with np.errstate(divide="ignore", invalid="ignore"):
                target_dt = _f64(time.code_dt) * np.sqrt(_f64(xparams.error_target) / _f64(error))
                mvc[2] = float(_f64(1.0) / target_dt)
        pending = state.cfl_pending
        if xparams.use_cfl_limit or xparams.use_order_based_error:   # :561
            if xparams.immediate_reset:
                # max_allreduce is the identity on one rank (and exact in general).
                time, reset = set_cfl_time_step_and_check_if_too_large(mvc, time, xparams)
            else:
                pending = True                                 # nb_max_allreduce (:573)
        else:
            time, reset = set_cfl_time_step_and_check_if_too_large(mvc, time, xparams)   # :582
        state = ExplicitState(gexp_1=gexp_1, gexp_2=gexp_2, gexp_3=gexp_3, istep_last=istep,
                              max_vel_components=tuple(mvc), cfl_pending=pending)
        if reset:                                              # :585-587
            return ExplicitOutcome(state, time, reset, True, True,
                                   reset and xparams.immediate_reset)
    # zip block (:599-609) is GATED.
    state = state._replace(istep_last=istep)                   # :611
    return ExplicitOutcome(state, time, reset, computed, False, reset and xparams.immediate_reset)


# --------------------------------------------------------------------------------------------
# deferred CFL check and end-of-step test
# --------------------------------------------------------------------------------------------


def nb_check_time_step_too_large(state: ExplicitState, time: TimeBlock, reset: bool,
                                 xparams: ExplicitParams):
    """``nb_check_time_step_too_large`` (nonlinear_terms.fpp:1191-1197).

    No-op with ``immediate_reset``; otherwise finish the (pending) reduction and call
    ``set_cfl_time_step_and_check_if_too_large`` with the module ``max_vel_components`` --
    also when no reduction was pending (``finish_nonblocking_max_vel_reduction`` returns early
    but the CFL routine still runs). Returns ``(state, time, reset)``.
    """
    if xparams.immediate_reset:
        return state, time, bool(reset)
    state = state._replace(cfl_pending=False)
    time, reset = set_cfl_time_step_and_check_if_too_large(state.max_vel_components, time, xparams)
    return state, time, reset


def check_time_step(state: ExplicitState, time: TimeBlock, exit: bool, xparams: ExplicitParams):
    """``gs2_reinit::check_time_step`` (gs2_reinit.f90:228-254), end of every step.

    ``reset = .false.``; nothing more if exiting; finish a deferred CFL check; then
    ``reset`` if ``code_dt*1.0 > code_dt_cfl`` (too large) or
    ``code_dt*1.0 < min(code_dt_max, (code_dt_cfl/delt_adj)/delt_cushion)`` (too small: the
    only place dt is increased). The dummy ``reset`` aliases ``run_parameters::reset``, which
    ``nb_check`` assigns; the two tests recompute it, so the alias is harmless. Returns
    ``(state, time, reset)``; ``reset_time_step`` itself belongs to the integration layer.
    """
    reset = False
    if exit:
        return state, time, reset
    state, time, reset = nb_check_time_step_too_large(state, time, reset, xparams)
    fac = 1.0
    if time.code_dt * fac > time.code_dt_cfl:
        reset = True
    if time.code_dt * fac < min(time.code_dt_max,
                                time.code_dt_cfl / xparams.delt_adj / xparams.delt_cushion):
        reset = True
    return state, time, reset
