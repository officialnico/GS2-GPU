"""End-to-end replication: run every case with the JAX driver and compare its
``<run>.out.nc`` with GS2's own output of the same full run (00_decisions D1).

Every variable of GS2's file that the JAX run writes is compared over the
common time records with a max-abs-relative tolerance (1e-8). Variables that
are zero up to round-off in GS2 (the particle flux with adiabatic electrons,
~1e-24; the energy exchange; the zonal density) cannot meet a relative test;
for those the reference is GS2's own round-off spread, measured by rerunning
GS2 with a different data layout (gpu/scripts/make_noise_refdata.sh), and the
JAX difference must stay within NOISE_FACTOR times that spread.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from conftest import all_cases, case_dir, have
from gs2gpu.compare import compare_files, structure_differences
from gs2gpu.driver import run

CASES = [c for c in all_cases() if have(c, "out")] or ["<no-refdata>"]
RTOL = 1e-8
ATOL = 1e-20
NOISE_FACTOR = 20.0


def gs2_noise(case: str) -> dict:
    """Per-variable max |GS2(layout B) - GS2(layout A)| (empty if not generated)."""
    other = case_dir(case) / "gs2_layout" / f"{case}.out.nc"
    if not other.exists():
        return {}
    return {d.name: d.absdiff for d in compare_files(case_dir(case) / f"{case}.out.nc", other)}


RESTART_CASES = [c for c in ("lin_single_ae", "nl_box_ae_dt_decrease_deferred") if c in CASES]


@pytest.mark.parametrize("case", RESTART_CASES or ["<no-refdata>"])
def test_restart_is_seamless(case, tmp_path):
    """Stopping after step K (a write step) and resuming from the checkpoint gives an
    out.nc bitwise identical to the uninterrupted run."""
    import netCDF4

    export = case_dir(case) / "export" / f"{case}.gs2export.nc"
    if not export.exists():
        pytest.skip(f"no export for {case}")
    full = tmp_path / "full.out.nc"
    whole = run(export, full)
    k = (whole.steps // 2) // 10 * 10 or 10
    part = tmp_path / "part.out.nc"
    run(export, part, nstep=k, checkpoint_every=k)
    resumed = run(export, part, restart=part.with_name(part.name + ".ckpt"))
    assert resumed.steps == whole.steps
    with netCDF4.Dataset(full) as a, netCDF4.Dataset(part) as b:
        for name, var in a.variables.items():
            x, y = np.asarray(var[...]), np.asarray(b.variables[name][...])
            assert x.shape == y.shape, name
            if x.dtype.kind in "fc":
                assert np.array_equal(x, y, equal_nan=True), name
            else:
                assert np.array_equal(x, y), name


def test_no_checkpoint_after_reset_abort(tmp_path):
    """A run that aborts inside reset_time_step (here at step 1, before any step
    completed) must not leave a checkpoint labelled as a completed step."""
    case = "nl_box_ae_reset_abort"
    export = case_dir(case) / "export" / f"{case}.gs2export.nc"
    if not export.exists():
        pytest.skip(f"no export for {case}")
    out = tmp_path / f"{case}.out.nc"
    result = run(export, out, checkpoint_every=1)
    assert "too much" in result.exit_reason
    assert not out.with_name(out.name + ".ckpt").exists()


@pytest.mark.parametrize("case", CASES)
def test_run_matches_gs2(case, tmp_path):
    export = case_dir(case) / "export" / f"{case}.gs2export.nc"
    if not export.exists():
        pytest.skip(f"no export for {case}")
    out = tmp_path / f"{case}.out.nc"
    result = run(export, out)
    ref = case_dir(case) / f"{case}.out.nc"
    diffs = compare_files(ref, out)
    assert diffs, "nothing compared"
    structure = structure_differences(ref, out)
    assert not structure, f"{case}: " + "; ".join(structure)
    noise = gs2_noise(case)

    def ok(d):
        if d.ok(RTOL, ATOL):
            return True
        n = noise.get(d.name)
        return n is not None and math.isfinite(d.absdiff) and d.absdiff <= NOISE_FACTOR * n

    bad = [d for d in diffs if not ok(d)]
    report = "\n".join(f"  {d.name}: rel={d.err:.3e} abs={d.absdiff:.3e} scale={d.scale:.3e} "
                       f"gs2-noise={noise.get(d.name, float('nan')):.3e}"
                       for d in sorted(bad, key=lambda d: -d.err if math.isfinite(d.err)
                                       else -math.inf)[:15])
    signal = [d for d in diffs if math.isfinite(d.err) and d.err <= RTOL]
    worst = max((d.err for d in signal), default=0.0)
    ratios = [d.absdiff / noise[d.name] for d in diffs
              if d not in signal and noise.get(d.name) and math.isfinite(d.absdiff)]
    print(f"{case}: {len(diffs)} variables, worst relative {worst:.3e} on {len(signal)}, "
          f"{len(ratios)} round-off-level variables at <= {max(ratios, default=0):.2f}x GS2's "
          f"own noise, {result.steps} steps, {len(result.resets)} resets, "
          f"{result.wall_seconds:.1f} s")
    assert not bad, f"{case}: {len(bad)} variables differ from GS2:\n{report}"
