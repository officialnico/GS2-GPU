"""Loader tests: every reference case loads, has canonical shapes, passes the
gates, and its trace is indexable (spec 09 / 09a, 11_api.md section 2)."""

from __future__ import annotations

import numpy as np
import pytest

from conftest import all_cases, export_of, have, trace_of
from gs2gpu.gates import check_gates, violations
from gs2gpu.types import (Coeffs, Dims, Geometry, KGrids, PointArrays, Species, VGrids,
                          fields_from_group)

CASES = [c for c in all_cases() if have(c)] or ["<no-refdata>"]


@pytest.mark.parametrize("case", CASES)
def test_export_loads_with_canonical_shapes(case):
    ex = export_of(case)
    dims = Dims.from_export(ex)
    assert dims.nth == 2 * dims.ntgrid + 1
    g7 = dims.g_shape
    f3 = dims.field_shape
    pts = PointArrays.from_export(ex)
    assert pts.vpa.shape == g7 and pts.vpa.dtype == np.float64
    assert pts.aj0.shape == g7[:5] + (dims.nth,)
    coeffs = Coeffs.from_export(ex)
    for name in ("a", "b", "r", "ainv"):
        arr = getattr(coeffs, name)
        assert arr.shape == g7 and arr.dtype == np.complex128, name
    assert coeffs.wstar.shape == (dims.nspec, dims.negrid, dims.naky)
    geo = Geometry.from_export(ex)
    assert geo.bmag.shape == (dims.nth,)
    kg = KGrids.from_export(ex)
    assert kg.kperp2.shape == f3 and kg.kwork_filter.dtype == bool
    vg = VGrids.from_export(ex)
    assert vg.wl.shape == (dims.nlambda, dims.nth) and vg.forbid.dtype == bool
    sp = Species.from_export(ex)
    assert sp.z.shape == (dims.nspec,)
    state = ex["/state_t0"]
    assert state["gnew"].shape == g7 and state["gnew"].dtype == np.complex128
    fields = fields_from_group(state, dims)
    assert fields.phi.shape == f3


@pytest.mark.parametrize("case", CASES)
def test_species_group_alias(case):
    ex = export_of(case)
    # The file stores the group as /species_ (netCDF-4 name clash with the
    # root dimension 'species'); both spellings must work.
    assert np.array_equal(ex["/species/spec_z"], ex["/species_/spec_z"])
    assert ex["/species"].path in ("/species", "/species_")


@pytest.mark.parametrize("case", CASES)
def test_gates_pass(case):
    ex = export_of(case)
    assert violations(ex) == []
    check_gates(ex)


@pytest.mark.parametrize("case", CASES)
def test_trace_index(case):
    tr = trace_of(case)
    assert tr.steps, "trace has no step records"
    first = tr.steps[0]
    rec = tr.record("advance_begin", first)
    nxt = tr.next(rec, "advance_end")
    assert nxt is not None and nxt.seq > rec.seq
    grp = tr.load(rec)
    assert grp["gnew"].shape == export_of(case).g_shape
    # Invariant of spec 09 section 6.2: step_end of step k == advance_begin of step k+1.
    for k, k1 in zip(tr.steps, tr.steps[1:]):
        if k1 != k + 1:
            continue
        end = tr.load(tr.record("step_end", k))["gnew"]
        begin = tr.load(tr.record("advance_begin", k1, attempt=1))["gnew"]
        assert np.array_equal(end, begin), (k, k1)
