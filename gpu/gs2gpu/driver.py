"""The GS2 main loop on JAX (``evolve_equations``, src/gs2_main.fpp:649-894;
spec 01 sections 1.5 and 3.1).

``run`` loads a ``<run>.gs2export.nc`` produced by ``bin/gs2_export``, copies
the exporter's ``<run>.out.nc`` (GS2's static part and the istep-0 record) to
the output path, and advances the run for ``nstep`` steps (or until an exit
condition), appending one record per ``nwrite`` steps exactly as GS2 does.
"""

from __future__ import annotations

import logging
import shutil
import time as walltime
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Optional

import jax.numpy as jnp
import numpy as np

from . import explicit as explicit_mod
from .diagnostics import Diagnostics
from .io.export import Export, load_export
from .reset import max_sim_time_reached, reset_time_step, update_time
from .step import DtSet, Model, State, advance, restore_after_reset
from .types import Fields

log = logging.getLogger("gs2gpu")


@dataclass
class RunResult:
    """Outcome of :func:`run`."""

    out_path: Optional[Path]
    steps: int
    exit_reason: str
    final_state: State
    resets: list = field(default_factory=list)      # (istep, kind, old_dt, new_dt)
    wall_seconds: float = 0.0


class DtCache:
    """Dt-dependent sets keyed by the exact (code_dt, vnmult) bits.

    GS2 rebuilds the coefficients on every reset; the rebuild is deterministic, so
    reusing a set built for the same dt gives the same numbers. The export's own set
    (GS2's matrices) is used whenever the run returns to the initial dt. At most
    ``max_sets`` rebuilt sets are kept (least recently used evicted; with a
    non-power-of-2 ``delt_adj`` the reachable dt values are path dependent)."""

    def __init__(self, model: Model, max_sets: int = 6):
        self.model = model
        self.max_sets = max_sets
        initial = model.dtset_from_export()
        vn = tuple(float(v) for v in model.initial_state().vnmult)
        self._initial_key = self._key(initial.code_dt, vn)
        self._sets = {self._initial_key: initial}
        self.initial = initial

    @staticmethod
    def _key(code_dt, vnmult):
        return (np.float64(code_dt).tobytes(), tuple(np.float64(v).tobytes() for v in vnmult))

    def get(self, code_dt: float, vnmult) -> DtSet:
        key = self._key(code_dt, tuple(float(v) for v in vnmult))
        if key in self._sets:
            self._sets[key] = self._sets.pop(key)          # most recently used last
            return self._sets[key]
        t0 = walltime.perf_counter()
        self._sets[key] = self.model.build_dtset(code_dt, vnmult)
        log.info("rebuilt dt-dependent arrays for code_dt=%r in %.1f s", code_dt,
                 walltime.perf_counter() - t0)
        rebuilt = [k for k in self._sets if k != self._initial_key]
        while len(rebuilt) > self.max_sets:
            del self._sets[rebuilt.pop(0)]
        return self._sets[key]


def _do_reset(model: Model, cache: DtCache, S: State, istep: int):
    """``reset_time_step`` + rebuild + state restoration (gs2_reinit.f90:104-225,
    gs2_init.fpp:292-375). Returns ``(state, dtset, exit, reason)``."""
    outcome = reset_time_step(S.time, istep, model.reinit, S.counters)
    S = S.replace(time=outcome.time)
    if outcome.exit:
        if outcome.torn_down:
            # GS2 has already zeroed gnew (gs2_reinit.f90:163) and torn the fields down
            # (fields.f90:618-620) before these two exits; the final diagnostics see that.
            zg = jnp.zeros_like(S.gnew)
            zf = jnp.zeros(model.dims.field_shape, jnp.complex128)
            zero = Fields(zf, zf, zf)
            S = S.replace(gnew=zg, old=zero, new=zero)
        return S, None, True, outcome.reason
    dts = cache.get(outcome.time.code_dt, S.vnmult)
    S = restore_after_reset(model, S, S.gnew)
    return S, dts, False, ""


def run(export: str | Path | Export, out_path: str | Path | None = None, *,
        nstep: Optional[int] = None, check_gates: bool = True,
        progress: Optional[Callable[[int, State], None]] = None,
        checkpoint_path: str | Path | None = None, checkpoint_every: int = 0,
        restart: str | Path | None = None) -> RunResult:
    """Run the case of ``export`` with GS2's numerics.

    ``out_path``: where to write ``<run>.out.nc`` (default: next to the export file,
    named ``<run>.jax.out.nc``). ``nstep`` overrides the input's ``nstep`` (the total
    step count of the run, also when restarting).

    ``checkpoint_every = N`` writes a checkpoint to ``checkpoint_path`` (default
    ``<out>.ckpt``) after every N-th step and at the end; ``restart = path`` resumes a
    run from such a checkpoint, appending to the same ``out_path``, as if it had never
    stopped (same step counter, same dt history, same diagnostics records).
    """
    from . import checkpoint as ckpt

    t_start = walltime.perf_counter()
    ex = export if isinstance(export, Export) else load_export(export,
                                                              allow_gates=not check_gates)
    model = Model(ex, check=check_gates)
    run_name = ex.path.name.replace(".gs2export.nc", "")
    if out_path is None:
        out_path = ex.path.with_name(f"{run_name}.jax.out.nc")
    out_path = Path(out_path)
    if checkpoint_path is None:
        checkpoint_path = out_path.with_name(out_path.name + ".ckpt")
    checkpoint_path = Path(checkpoint_path)
    # GS2 looks for <run_name>.stop in the run directory (job_manage.fpp:151); accept it
    # next to the export (the run directory of gs2_export) or next to the output.
    stop_files = {ex.path.with_name(f"{run_name}.stop"), out_path.with_name(f"{run_name}.stop")}
    rp = ex["/run_parameters"]
    avail_cpu_time = float(rp["avail_cpu_time"])
    margin_cpu_time = float(rp["margin_cpu_time"])

    cache = DtCache(model)
    if restart is not None:
        saved = ckpt.load(restart)
        if Path(saved["export_path"]).resolve() != ex.path.resolve():
            log.warning("checkpoint %s was written for %s, restarting with %s", restart,
                        saved["export_path"], ex.path)
        if not out_path.exists():
            raise FileNotFoundError(f"restart needs the run's existing output {out_path}")
        S = saved["state"]
        first = saved["istep"] + 1
        dts = cache.get(S.time.code_dt, S.vnmult)
        diag = Diagnostics.from_export(ex, out_path)
        diag.state = saved["diag_state"]
        log.info("restarting %s after step %d", run_name, saved["istep"])
    else:
        src_out = ex.path.with_name(f"{run_name}.out.nc")
        if src_out.exists() and src_out.resolve() != out_path.resolve():
            shutil.copyfile(src_out, out_path)
        S = model.initial_state()
        first = 1
        dts = cache.initial
        diag = Diagnostics.from_export(ex, out_path)

    def save_checkpoint(at_step: int):
        ckpt.save(checkpoint_path, S, at_step, export_path=ex.path, diag_state=diag.state)

    nstep_run = model.run.nstep if nstep is None else int(nstep)
    exit_ = False
    reason = "nstep"
    clean_exit = True        # False after a reset_time_step abort: no final checkpoint
    resets = []
    istep = first - 1
    completed = True
    for istep in range(first, nstep_run + 1):
        S, reset = advance(model, dts, S, istep)
        if reset:                                            # gs2_main.fpp:739-752
            old_dt = S.time.code_dt
            S, new_dts, exit_, why = _do_reset(model, cache, S, istep)
            S = S.replace(reset=False)
            if not exit_:
                dts = new_dts
                resets.append((istep, "immediate", old_dt, S.time.code_dt))
                S, again = advance(model, dts, S, istep)
                if again:
                    raise RuntimeError(f"istep {istep}: reset requested again after a reset")
            else:
                reason, clean_exit = why, False
        if exit_:
            completed = False
            break
        S = S.replace(time=update_time(S.time))              # gs2_main.fpp:771
        gnew, g, exit_ = diag.step(istep, S.gnew, S.g, S.old, S.new, S.time.code_dt,
                                   S.time.user_time, exit_)
        S = S.replace(gnew=gnew, g=g)
        if exit_:
            reason = "diagnostics"
        if not exit_:                                        # gs2_main.fpp:792-823
            xstate, time_block, reset = explicit_mod.check_time_step(
                S.xstate, S.time, exit_, model.x_params)
            S = S.replace(xstate=xstate, time=time_block)
            if reset:
                old_dt = S.time.code_dt
                S, new_dts, exit_, why = _do_reset(model, cache, S, istep)
                S = S.replace(reset=False)
                if exit_:
                    reason, clean_exit = why, False
                else:
                    dts = new_dts
                    resets.append((istep, "end_of_step", old_dt, S.time.code_dt))
            if (not exit_ and model.run.ncheck_stop > 0 and istep % model.run.ncheck_stop == 0
                    and any(f.exists() for f in stop_files)):
                exit_, reason = True, "stop file"
            if (not exit_ and walltime.perf_counter() - t_start
                    >= avail_cpu_time - margin_cpu_time):          # job_manage.fpp checktime
                exit_, reason = True, "avail_cpu_time"
        if not exit_ and max_sim_time_reached(S.time, model.run.max_sim_time):
            exit_, reason = True, "max_sim_time"
        if progress is not None:
            progress(istep, S)
        if checkpoint_every and istep % checkpoint_every == 0 and not exit_:
            save_checkpoint(istep)
        if exit_:
            completed = False
            break
    last = istep if not completed else nstep_run              # gs2_main.fpp:877-885
    if checkpoint_every and clean_exit:
        save_checkpoint(last)                                 # before the final forced write
    elif checkpoint_every:
        log.warning("run aborted in reset_time_step at step %d (%s); the last periodic "
                    "checkpoint is the restart point", istep, reason)
    gnew, g = diag.finalize(last, S.gnew, S.g, S.old, S.new, S.time.code_dt,
                            S.time.user_time, exit_)
    S = S.replace(gnew=gnew, g=g)
    return RunResult(out_path=out_path, steps=last, exit_reason=reason, final_state=S,
                     resets=resets, wall_seconds=walltime.perf_counter() - t_start)
