#!/usr/bin/env python3
"""Bitwise comparison of two netCDF files (all groups, all variables).

Used to prove that the GPU-port edits of the Fortran sources do not change GS2's numerics
(gpu/spec/09_export_infra.md §4.4) and, by make_refdata.sh, that a traced gs2_export run
follows the same trajectory as bin/gs2.

Usage:
    compare_netcdf.py A.nc B.nc [--attributes] [--ignore-var NAME ...] [--ignore-att NAME ...]
                                [--records-of-b] [--quiet]

* Every variable present in either file is compared: same dimensions, same dtype, same shape
  and identical bytes (so -0.0 vs +0.0 and NaN payloads count as differences).
* --records-of-b: compare only the time records (dimension `t`) of B whose time value occurs
  (bitwise) in A, against those records of A. Used to compare a short traced run (B) with a full
  run (A); B's final forced record (step N not a multiple of nwrite) is skipped.
* --time-only: compare only variables whose first dimension is the time dimension `t`.
* Global/variable attributes are compared only with --attributes; date_created, id and the
  git/build attributes are always ignored.
Exit status 0 if identical, 1 otherwise.
"""
import argparse
import sys

import numpy as np
import netCDF4

ALWAYS_IGNORED_ATTS = {"date_created", "id", "creation_date", "git_hash", "git_state",
                       "build_config", "software_version"}


def raw(arr):
    a = np.ascontiguousarray(np.asarray(arr))
    if a.dtype.kind in "OU":
        return a.astype(str).tobytes()
    return a.tobytes()


def walk(grp, prefix=""):
    yield prefix or "/", grp
    for name, sub in grp.groups.items():
        yield from walk(sub, f"{prefix}/{name}")


def time_index_map(da, db):
    """Records of B whose t value occurs exactly once (bitwise) in A, and the matching A indices."""
    ta = np.asarray(da.variables["t"][:], dtype=np.float64)
    tb = np.asarray(db.variables["t"][:], dtype=np.float64)
    ia, ib = [], []
    for j, v in enumerate(tb):
        hits = np.nonzero(ta.view(np.uint64) == np.float64(v).view(np.uint64))[0]
        if len(hits) == 1:
            ia.append(int(hits[0]))
            ib.append(j)
    return np.array(ia, dtype=int), np.array(ib, dtype=int), tb


def compare(path_a, path_b, attributes=False, ignore_var=(), ignore_att=(), records_of_b=False,
            quiet=False, time_only=False):
    da = netCDF4.Dataset(path_a)
    db = netCDF4.Dataset(path_b)
    for d in (da, db):
        d.set_auto_mask(False)
        d.set_auto_scale(False)
    ignored_atts = ALWAYS_IGNORED_ATTS | set(ignore_att)
    ndiff = 0
    nvars = 0
    msgs = []

    tidx = bidx = None
    if records_of_b:
        tidx, bidx, tb = time_index_map(da, db)
        print(f"compare_netcdf: {len(bidx)} of {len(tb)} time records of B matched in A "
              f"(t = {', '.join(f'{x:g}' for x in tb[bidx])})")
        if len(bidx) == 0:
            msgs.append(f"records: none of B's t values {tb} found in A")
            ndiff += 1

    groups_a = dict(walk(da))
    groups_b = dict(walk(db))
    for gpath in sorted(set(groups_a) | set(groups_b)):
        ga, gb = groups_a.get(gpath), groups_b.get(gpath)
        if ga is None or gb is None:
            msgs.append(f"group {gpath}: only in {'B' if ga is None else 'A'}")
            ndiff += 1
            continue
        if attributes:
            for att in sorted(set(ga.ncattrs()) | set(gb.ncattrs())):
                if att in ignored_atts:
                    continue
                va = ga.getncattr(att) if att in ga.ncattrs() else None
                vb = gb.getncattr(att) if att in gb.ncattrs() else None
                if va is None or vb is None or raw(va) != raw(vb):
                    msgs.append(f"attribute {gpath}:{att} differs")
                    ndiff += 1
        for vname in sorted(set(ga.variables) | set(gb.variables)):
            full = f"{gpath.rstrip('/')}/{vname}"
            if vname in ignore_var or full in ignore_var:
                continue
            va, vb = ga.variables.get(vname), gb.variables.get(vname)
            if time_only and ((va is None or "t" not in va.dimensions[:1]) and
                              (vb is None or "t" not in vb.dimensions[:1])):
                continue
            if va is None or vb is None:
                msgs.append(f"variable {full}: only in {'B' if va is None else 'A'}")
                ndiff += 1
                continue
            nvars += 1
            if va.dimensions != vb.dimensions or va.dtype != vb.dtype:
                msgs.append(f"variable {full}: dims/dtype {va.dimensions}/{va.dtype} vs "
                            f"{vb.dimensions}/{vb.dtype}")
                ndiff += 1
                continue
            a = va[...]
            b = vb[...]
            if records_of_b and tidx is not None and len(va.dimensions) > 0 \
                    and va.dimensions[0] == "t":
                a = np.asarray(a)[tidx]
                b = np.asarray(b)[bidx]
            if np.shape(a) != np.shape(b):
                msgs.append(f"variable {full}: shape {np.shape(a)} vs {np.shape(b)}")
                ndiff += 1
                continue
            if raw(a) != raw(b):
                an = np.asarray(a)
                bn = np.asarray(b)
                extra = ""
                if an.dtype.kind in "fc":
                    with np.errstate(invalid="ignore"):
                        d = np.abs(an.astype(np.complex128) - bn.astype(np.complex128))
                        extra = f" (max|diff|={np.nanmax(d):.3e}, max|A|={np.nanmax(np.abs(an)):.3e})"
                msgs.append(f"variable {full}: values differ{extra}")
                ndiff += 1
            if attributes:
                for att in sorted(set(va.ncattrs()) | set(vb.ncattrs())):
                    if att in ignored_atts:
                        continue
                    xa = va.getncattr(att) if att in va.ncattrs() else None
                    xb = vb.getncattr(att) if att in vb.ncattrs() else None
                    if xa is None or xb is None or raw(xa) != raw(xb):
                        msgs.append(f"attribute {full}:{att} differs")
                        ndiff += 1
    da.close()
    db.close()
    if not quiet or ndiff:
        for m in msgs:
            print("  " + m)
    print(f"compare_netcdf: {nvars} variables compared, {ndiff} difference(s): "
          f"{'IDENTICAL' if ndiff == 0 else 'DIFFERENT'}  [{path_a} vs {path_b}]")
    return ndiff == 0


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("a")
    p.add_argument("b")
    p.add_argument("--attributes", action="store_true")
    p.add_argument("--ignore-var", nargs="*", default=[])
    p.add_argument("--ignore-att", nargs="*", default=[])
    p.add_argument("--records-of-b", action="store_true")
    p.add_argument("--time-only", action="store_true")
    p.add_argument("--quiet", action="store_true")
    args = p.parse_args(argv)
    ok = compare(args.a, args.b, args.attributes, args.ignore_var, args.ignore_att,
                 args.records_of_b, args.quiet, args.time_only)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
