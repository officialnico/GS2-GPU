"""Compare two GS2 ``<run>.out.nc`` files variable by variable.

Typical use: GS2's own output of a case against the output the JAX stepper
wrote for the same case::

    python -m gs2gpu.compare gpu/refdata/lin_single_ae/lin_single_ae.out.nc jax.out.nc

For every variable present in both files the error is
``max|got - ref| / max|ref|`` over the records both files contain (time
variables are cut to the common number of records; the ``t`` arrays
themselves are compared too). Complex variables (trailing ``ri``) are
compared as complex numbers.
"""

from __future__ import annotations

import argparse
import math
import sys
from dataclasses import dataclass
from pathlib import Path

import netCDF4
import numpy as np


@dataclass
class VarDiff:
    name: str
    shape: tuple
    nrec: int | None
    err: float
    scale: float

    absdiff: float = 0.0

    def ok(self, rtol: float, atol: float = 0.0) -> bool:
        return self.err <= rtol or self.absdiff <= atol


def _complexify(var, data):
    if var.dimensions and var.dimensions[-1] == "ri" and data.shape[-1] == 2:
        return data[..., 0] + 1j * data[..., 1]
    return data


def _read(var):
    data = np.asarray(var[...])
    if data.dtype.kind in "SUO":
        return None
    return _complexify(var, data.astype(np.complex128 if data.dtype.kind == "c" else data.dtype))


def compare_files(ref_path: str | Path, got_path: str | Path, names=None,
                  max_records: int | None = None) -> list[VarDiff]:
    """Per-variable differences between two out.nc files."""
    out = []
    with netCDF4.Dataset(ref_path) as ref, netCDF4.Dataset(got_path) as got:
        ref.set_auto_mask(False)
        got.set_auto_mask(False)
        nref = len(ref.dimensions["t"]) if "t" in ref.dimensions else 0
        ngot = len(got.dimensions["t"]) if "t" in got.dimensions else 0
        nrec = min(nref, ngot)
        if max_records is not None:
            nrec = min(nrec, max_records)
        for name, rvar in ref.variables.items():
            if names and name not in names:
                continue
            if name not in got.variables:
                continue
            gvar = got.variables[name]
            r = _read(rvar)
            g = _read(gvar)
            if r is None or g is None:
                continue
            timed = bool(rvar.dimensions) and rvar.dimensions[0] == "t"
            if timed:
                r = r[:nrec]
                g = g[:nrec]
            if r.shape != g.shape:
                out.append(VarDiff(name, r.shape, nrec if timed else None, math.inf, math.nan,
                                   math.inf))
                continue
            if r.size == 0:
                continue
            scale = float(np.max(np.abs(r)))
            diff = float(np.max(np.abs(g.astype(np.complex128) - r.astype(np.complex128))))
            err = diff / scale if scale > 0 else diff
            out.append(VarDiff(name, r.shape, nrec if timed else None, err, scale, diff))
    return out


def structure_differences(ref_path: str | Path, got_path: str | Path) -> list[str]:
    """Structural mismatches :func:`compare_files` cannot see because it compares the
    common records only: different numbers of time records, and variables of the
    reference that the candidate lacks."""
    problems = []
    with netCDF4.Dataset(ref_path) as ref, netCDF4.Dataset(got_path) as got:
        nref = len(ref.dimensions["t"]) if "t" in ref.dimensions else 0
        ngot = len(got.dimensions["t"]) if "t" in got.dimensions else 0
        if nref != ngot:
            problems.append(f"time records: reference {nref}, candidate {ngot}")
        missing = sorted(set(ref.variables) - set(got.variables))
        if missing:
            problems.append(f"variables missing from the candidate: {', '.join(missing)}")
    return problems


def growth_rates(path: str | Path):
    """``omega_average`` at the last record, shape (naky, ntheta0), complex:
    real part = frequency, imaginary part = growth rate."""
    with netCDF4.Dataset(path) as d:
        d.set_auto_mask(False)
        w = np.asarray(d["omega_average"][-1])
    return w[..., 0] + 1j * w[..., 1]


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("reference")
    parser.add_argument("candidate")
    parser.add_argument("--rtol", type=float, default=1e-8)
    parser.add_argument("--atol", type=float, default=1e-20,
                        help="absolute floor: variables that are round-off noise around zero "
                             "(e.g. the particle flux with adiabatic electrons) pass if "
                             "max|got - ref| <= atol")
    parser.add_argument("--vars", nargs="*", default=None)
    parser.add_argument("--max-records", type=int, default=None)
    parser.add_argument("--all", action="store_true", help="also print passing variables")
    args = parser.parse_args(argv)
    diffs = compare_files(args.reference, args.candidate, args.vars, args.max_records)
    structure = [] if args.vars or args.max_records else structure_differences(
        args.reference, args.candidate)
    for problem in structure:
        print(f"FAIL {problem}")
    worst = 0.0
    nfail = 0
    for d in sorted(diffs, key=lambda d: -d.err if math.isfinite(d.err) else -math.inf):
        bad = not d.ok(args.rtol, args.atol)
        nfail += bad
        if math.isfinite(d.err):
            worst = max(worst, d.err)
        if bad or args.all:
            rec = f" records={d.nrec}" if d.nrec is not None else ""
            print(f"{'FAIL' if bad else 'ok  '} {d.name:40s} err={d.err:.3e} "
                  f"abs={d.absdiff:.3e} scale={d.scale:.3e} shape={d.shape}{rec}")
    print(f"{len(diffs)} variables compared, {nfail} above rtol={args.rtol:g} "
          f"(atol={args.atol:g}), worst relative={worst:.3e}")
    return 1 if (nfail or structure) else 0


if __name__ == "__main__":
    sys.exit(main())
