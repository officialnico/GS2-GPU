"""Whole-step tests: one JAX ``advance`` from GS2's state at ``advance_begin`` must
reproduce GS2's state at ``advance_end`` (00_decisions D1: 1e-12), for every
traced step of every case. Steps where GS2 performed an immediate reset are
checked for the reset decision here and end to end in test_end_to_end.py.
"""

from __future__ import annotations

from functools import lru_cache

import numpy as np
import pytest

from conftest import RTOL_STEP, all_cases, assert_close, export_of, have, trace_of
from gs2gpu import exb as exb_mod
from gs2gpu import explicit as explicit_mod
from gs2gpu.step import Model, State, advance
from gs2gpu.types import TimeBlock, fields_from_group, to_jax

CASES = [c for c in all_cases() if have(c, "trace")] or ["<no-refdata>"]


@lru_cache(maxsize=2)
def model_of(case: str) -> Model:
    return Model(export_of(case))


@lru_cache(maxsize=2)
def dtset_of(case: str):
    return model_of(case).dtset_from_export()


def state_from_record(model: Model, grp, istep: int) -> State:
    """The state at an ``advance_begin`` probe (spec 09 section 5.3 table)."""
    new = to_jax(fields_from_group(grp, model.dims, new=True))
    exb = None
    if model.run.exb_active:
        exb = exb_mod.ExbState.from_group(grp, model.exb_params, istep_last=istep - 1)
    import jax.numpy as jnp
    return State(
        g=jnp.asarray(grp["gnew"]), gnew=jnp.asarray(grp["gnew"]), old=new, new=new,
        time=TimeBlock.from_group(grp), reset=bool(np.asarray(grp["reset"]).item()),
        xstate=explicit_mod.ExplicitState.from_group(grp), exb=exb,
        vnmult=np.asarray(grp["vnmult"], dtype=np.float64), timeadv_istep_last=istep - 1,
    )


@pytest.mark.parametrize("case", CASES)
def test_advance_matches_trace(case):
    model = model_of(case)
    tr = trace_of(case)
    base = dtset_of(case)
    dtsets = {base.code_dt: base}
    checked = 0
    for istep in tr.steps:
        attempts = tr.attempts(istep)
        begin = tr.load(tr.record("advance_begin", istep, attempt=1))
        S = state_from_record(model, begin, istep)
        dt = S.time.code_dt
        if dt not in dtsets:
            dtsets[dt] = model.build_dtset(dt, S.vnmult)
        S2, reset = advance(model, dtsets[dt], S, istep)
        if len(attempts) > 1:
            assert reset, f"{case} istep={istep}: GS2 reset immediately, JAX did not"
            checked += 1
            continue
        end = tr.load(tr.record("advance_end", istep, attempt=attempts[-1]))
        if bool(np.asarray(end["reset"]).item()):
            # GS2 requested a reset and never redid the step (reset_time_step aborted)
            assert reset, f"{case} istep={istep}: GS2 requested a reset, JAX did not"
            checked += 1
            continue
        assert not reset, f"{case} istep={istep}: JAX requested a reset, GS2 did not"
        what = f"{case} istep={istep}"
        assert_close(end["gnew"], np.asarray(S2.gnew), RTOL_STEP, f"{what} gnew")
        for name, value in zip(("phinew", "aparnew", "bparnew"), S2.new):
            if name in end:
                assert_close(end[name], np.asarray(value), RTOL_STEP, f"{what} {name}")
        checked += 1
    assert checked or len(tr.steps) == 0
