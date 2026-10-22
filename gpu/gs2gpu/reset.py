"""Host-side time-step control: ``update_time``, ``check_time_step`` and the
scalar part of ``reset_time_step``.

Literal float64 ports of src/gs2_time.f90:161-174 (update_time),
src/gs2_reinit.f90:228-254 (check_time_step) and src/gs2_reinit.f90:104-225
(reset_time_step; the coefficient rebuild and state restoration it triggers
are done by ``gs2gpu.step``). Spec 01 sections 1.11, 3.4-3.6.
"""

from __future__ import annotations

from dataclasses import dataclass

from .types import TimeBlock

TIMESTEP_CHANGE_ATTEMPTS_LIMIT = 10  # gs2_reinit.f90:126


@dataclass(frozen=True)
class ReinitParams:
    """``gs2_reinit`` knobs (export group ``/gs2_reinit``)."""

    delt_adj: float
    delt_cushion: float
    delt_minimum: float
    dt0: float
    abort_rapid_time_step_change: bool
    in_memory: bool
    immediate_reset: bool

    @classmethod
    def from_export(cls, ex) -> "ReinitParams":
        g = ex["/gs2_reinit"]
        return cls(
            delt_adj=float(g["delt_adj"]), delt_cushion=float(g["delt_cushion"]),
            delt_minimum=float(g["delt_minimum"]), dt0=float(g["dt0"]),
            abort_rapid_time_step_change=bool(g["abort_rapid_time_step_change"]),
            in_memory=bool(g["in_memory"]),
            immediate_reset=bool(ex["/run_parameters/immediate_reset"]),
        )


@dataclass
class ResetCounters:
    """SAVE variables of ``reset_time_step`` (gs2_reinit.f90:120-121)."""

    istep_last: int = -1
    nconsec: int = 0


def update_time(t: TimeBlock) -> TimeBlock:
    """``update_time`` (gs2_time.f90:161-174)."""
    return t._replace(
        code_dt_prev2=t.code_dt_prev1,
        code_dt_prev1=t.code_dt,
        code_dt_old=t.code_dt,
        code_time=t.code_time + t.code_dt,
        user_time=t.user_time + t.user_dt,
    )


def save_dt_cfl(t: TimeBlock, delt_cfl: float) -> TimeBlock:
    """``save_dt_cfl`` (gs2_time.f90:177-183)."""
    return t._replace(code_dt_cfl=float(delt_cfl), user_dt_cfl=float(delt_cfl))


def check_time_step(t: TimeBlock, p: ReinitParams, exit_: bool) -> bool:
    """``check_time_step`` (gs2_reinit.f90:228-254), after any deferred CFL
    result has been folded into ``t.code_dt_cfl`` by the caller."""
    if exit_:
        return False
    fac = 1.0
    reset = False
    if t.code_dt * fac > t.code_dt_cfl:
        reset = True
    if t.code_dt * fac < min(t.code_dt_max, t.code_dt_cfl / p.delt_adj / p.delt_cushion):
        reset = True
    return reset


@dataclass(frozen=True)
class ResetOutcome:
    """Result of the scalar part of ``reset_time_step``."""

    time: TimeBlock
    exit: bool
    rebuild: bool          # True: coefficients must be rebuilt and the state restored
    reason: str = ""
    # True for the two exits after GS2 has zeroed gnew and torn the fields down
    # (gs2_reinit.f90:163-166; fields.f90:618-620): the final diagnostics then see
    # gnew = 0 and zero fields (g is kept).
    torn_down: bool = False


def reset_time_step(t: TimeBlock, istep: int, p: ReinitParams,
                    counters: ResetCounters) -> ResetOutcome:
    """Scalar part of ``reset_time_step`` (gs2_reinit.f90:104-225).

    Mutates ``counters`` like the Fortran SAVE variables. On ``rebuild`` the
    caller rebuilds every dt-dependent array for ``outcome.time.code_dt`` and
    restores the state (gs2_init.fpp:292-375), then the step continues.
    """
    if counters.istep_last + 1 == istep:
        counters.nconsec += 1
    else:
        counters.nconsec = 0

    if counters.nconsec > 4 and p.abort_rapid_time_step_change:
        return ResetOutcome(t, exit=True, rebuild=False,
                            reason="Time step changing rapidly.  Abort run.")

    if t.code_dt / p.delt_adj <= t.code_dt_min:
        return ResetOutcome(t._replace(code_dt=t.code_dt_min), exit=True, rebuild=False,
                            reason="Time step wants to fall below delt_min.  Abort run.")

    fac = 1.0
    code_dt = t.code_dt
    attempts = 0
    while attempts < TIMESTEP_CHANGE_ATTEMPTS_LIMIT:
        if code_dt * fac > t.code_dt_cfl:
            code_dt = code_dt / p.delt_adj                                   # reduce_time_step
        elif code_dt * fac < min(t.code_dt_max, t.code_dt_cfl / p.delt_adj / p.delt_cushion):
            code_dt = min(code_dt * p.delt_adj, t.code_dt_max)               # increase_time_step
        else:
            break
        attempts += 1

    if attempts == TIMESTEP_CHANGE_ATTEMPTS_LIMIT:
        # GS2 leaves code_dt changed but user_dt untouched and aborts.
        return ResetOutcome(t._replace(code_dt=code_dt), exit=True, rebuild=False,
                            reason="Attempting to change the timestep too much in one go",
                            torn_down=True)

    t = t._replace(code_dt=code_dt, user_dt=code_dt)                         # save_dt
    if t.code_dt <= t.code_dt_min:
        return ResetOutcome(t._replace(code_dt=t.code_dt_min), exit=True, rebuild=False,
                            reason="Time step wants to fall below delt_min.  Abort run.",
                            torn_down=True)

    counters.istep_last = istep
    return ResetOutcome(t, exit=False, rebuild=True)


def max_sim_time_reached(t: TimeBlock, max_sim_time: float) -> bool:
    """Stopping test of gs2_main.fpp:833 (uses the post-reset user_dt)."""
    return t.user_time + t.user_dt > max_sim_time


__all__ = ["ReinitParams", "ResetCounters", "ResetOutcome", "update_time", "save_dt_cfl",
           "check_time_step", "reset_time_step", "max_sim_time_reached"]
