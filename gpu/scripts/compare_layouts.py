#!/usr/bin/env python3
"""Compare two GS2 export files of the same case run with different `layout` strings.

The canonical index expansion (gpu/spec/09_export_infra.md §2) must make every exported array
independent of the layout. Arrays computed point by point (grids, geometry, Bessel functions,
vpa/vpar/wdrift, the invert_rhs coefficients a/b/r/ainv/g_h, collision LE matrices, ...) must be
bitwise identical. Arrays that GS2 obtains through velocity-space sums (gamtot*, the field-solve
matrices, the conservation vectors, t = 0 fields and everything downstream of them) are summed in
iglo order, which depends on the layout, so they may differ at round-off level; they are reported
with their relative difference (must stay below --rtol).

Usage: compare_layouts.py A.gs2export.nc B.gs2export.nc [--rtol 1e-13]
Exit status 0 if the static set is bitwise identical and the rest is within rtol.
"""
import argparse
import sys

import numpy as np
import netCDF4

# Expected to depend on the layout by construction
LAYOUT_DEFINED = {"/gs2_layouts/iglo_of", "/gs2_layouts/ile_of", "/gs2_layouts/g_lo_ik_ord",
                  "/gs2_layouts/g_lo_it_ord", "/gs2_layouts/g_lo_il_ord", "/gs2_layouts/g_lo_ie_ord",
                  "/gs2_layouts/g_lo_is_ord", "/gs2_layouts/le_lo_ig_ord", "/gs2_layouts/le_lo_it_ord",
                  "/gs2_layouts/le_lo_ik_ord", "/gs2_layouts/le_lo_is_ord", "/nonlinear_terms/accelerated"}
# Velocity-space sums (iglo / ile order) and everything computed from them
SUM_DERIVED_PREFIXES = ("/dist_fn/gamtot", "/dist_fn/inv_phi_denominator_g", "/dist_fn/inv_bpar_denominator_g",
                        "/dist_fn/awgt", "/fields_implicit/aminv_c", "/fields_local/resp_",
                        "/collisions/z0le", "/collisions/s0le", "/collisions/w0le", "/collisions/bz0le",
                        "/collisions/bs0le", "/collisions/bw0le", "/collisions/z0", "/collisions/s0",
                        "/collisions/w0", "/collisions/bz0", "/collisions/bs0", "/collisions/bw0",
                        "/state_t0", "/state_t0_prediag", "/dt_ladder")


def walk(grp, prefix=""):
    yield prefix or "/", grp
    for name, sub in grp.groups.items():
        yield from walk(sub, f"{prefix}/{name}")


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("a")
    p.add_argument("b")
    p.add_argument("--rtol", type=float, default=1e-13)
    args = p.parse_args(argv)
    da, db = netCDF4.Dataset(args.a), netCDF4.Dataset(args.b)
    for d in (da, db):
        for _, g in walk(d):
            g.set_auto_mask(False)
    print(f"layouts: {da.getncattr('layout')} vs {db.getncattr('layout')}")
    ga, gb = dict(walk(da)), dict(walk(db))
    n_static = n_layout = n_sum_same = 0
    bad_static, sum_diffs, bad_sum, missing = [], [], [], []
    for gpath, grp in ga.items():
        if gpath not in gb:
            missing.append(gpath)
            continue
        for name, va in grp.variables.items():
            full = f"{gpath.rstrip('/')}/{name}"
            if name not in gb[gpath].variables:
                missing.append(full)
                continue
            a = np.asarray(va[...])
            b = np.asarray(gb[gpath].variables[name][...])
            same = a.shape == b.shape and a.dtype == b.dtype and a.tobytes() == b.tobytes()
            if full in LAYOUT_DEFINED or gpath == "/" and name in ("glo",):
                n_layout += 1
                continue
            if any(full.startswith(pfx) for pfx in SUM_DERIVED_PREFIXES):
                if same:
                    n_sum_same += 1
                else:
                    scale = float(np.max(np.abs(b))) if b.size and b.dtype.kind == "f" else 1.0
                    rel = float(np.max(np.abs(a - b))) / (scale or 1.0) if b.dtype.kind == "f" else np.inf
                    sum_diffs.append((full, rel))
                    if not rel <= args.rtol:
                        bad_sum.append((full, rel))
                continue
            n_static += 1
            if not same:
                bad_static.append(full)
    print(f"  static per-point arrays/scalars bitwise identical: {n_static - len(bad_static)}/{n_static}")
    for f in bad_static:
        print(f"    DIFFERENT: {f}")
    print(f"  layout-defined variables (iglo_of, orders, ...): {n_layout} (not compared)")
    print(f"  velocity-sum derived: {n_sum_same} bitwise identical, {len(sum_diffs)} differ at round-off:")
    for f, r in sorted(sum_diffs, key=lambda x: -x[1])[:12]:
        print(f"    {f:55s} max|diff|/max|ref| = {r:.2e}")
    if len(sum_diffs) > 12:
        print(f"    ... {len(sum_diffs) - 12} more")
    for m in missing:
        print(f"    MISSING in B: {m}")
    ok = not bad_static and not bad_sum and not missing
    print(f"compare_layouts: {'PASS' if ok else 'FAIL'}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
