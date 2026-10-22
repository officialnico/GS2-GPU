"""Command line: ``gs2gpu run <run>.gs2export.nc`` and ``gs2gpu compare ref new``.

Typical workflow::

    mpirun -np 1 bin/gs2_export case.in          # Fortran front end (CPU, seconds)
    gs2gpu run case.gs2export.nc -o case.jax.out.nc --checkpoint-every 500
    gs2gpu run case.gs2export.nc -o case.jax.out.nc --restart case.jax.out.nc.ckpt   # resume
    gs2gpu compare case_gs2.out.nc case.jax.out.nc
"""

from __future__ import annotations

import argparse
import logging
import sys


def _run(args) -> int:
    from .driver import run

    def progress(istep, state):
        if args.progress and istep % args.progress == 0:
            print(f"step {istep:6d}  t = {state.time.user_time:.6g}  dt = {state.time.code_dt:.6g}",
                  flush=True)

    result = run(args.export, args.out, nstep=args.nstep, check_gates=not args.no_gates,
                 progress=progress, checkpoint_path=args.checkpoint,
                 checkpoint_every=args.checkpoint_every, restart=args.restart)
    print(f"{result.steps} steps, exit: {result.exit_reason}, "
          f"{len(result.resets)} time-step changes, {result.wall_seconds:.1f} s")
    for istep, kind, old, new in result.resets:
        print(f"  step {istep}: {kind} reset, dt {old:.6g} -> {new:.6g}")
    print(f"output: {result.out_path}")
    return 0


def _bench(args) -> int:
    from .bench import bench

    r = bench(args.export, steps=args.steps, warmup=args.warmup)
    print(f"device            {r.device}")
    print(f"g shape           {r.g_shape}  ({r.g_bytes / 2**20:.1f} MiB per distribution function)")
    print(f"warm-up           {r.warmup_s:.1f} s ({args.warmup} steps, includes compilation)")
    print(f"time per step     {r.ms_per_step:.2f} ms over {r.steps} steps ({r.resets} dt resets)")
    if r.peak_bytes is not None:
        print(f"peak device mem   {r.peak_bytes / 2**30:.2f} GiB")
    return 0


def _compare(args) -> int:
    from .compare import main as compare_main

    argv = [args.reference, args.candidate, "--rtol", str(args.rtol)]
    if args.all:
        argv.append("--all")
    return compare_main(argv)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(prog="gs2gpu", description="GS2 time stepping on JAX")
    parser.add_argument("-v", "--verbose", action="store_true")
    sub = parser.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("run", help="run a case from its gs2_export file")
    p.add_argument("export", help="<run>.gs2export.nc written by bin/gs2_export")
    p.add_argument("-o", "--out", default=None, help="output <run>.out.nc path")
    p.add_argument("--nstep", type=int, default=None, help="override nstep")
    p.add_argument("--progress", type=int, default=0, help="print every N steps")
    p.add_argument("--checkpoint-every", type=int, default=0, metavar="N",
                   help="write a restart checkpoint every N steps and at the end")
    p.add_argument("--checkpoint", default=None, help="checkpoint path (default <out>.ckpt)")
    p.add_argument("--restart", default=None, metavar="CKPT",
                   help="resume from a checkpoint, appending to the same --out file")
    p.add_argument("--no-gates", action="store_true", help="debugging only: skip the gate check")
    p.set_defaults(func=_run)
    b = sub.add_parser("bench", help="time the step on the current device")
    b.add_argument("export")
    b.add_argument("--steps", type=int, default=50)
    b.add_argument("--warmup", type=int, default=3)
    b.set_defaults(func=_bench)
    c = sub.add_parser("compare", help="compare two out.nc files")
    c.add_argument("reference")
    c.add_argument("candidate")
    c.add_argument("--rtol", type=float, default=1e-8)
    c.add_argument("--all", action="store_true")
    c.set_defaults(func=_compare)
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO if args.verbose else logging.WARNING)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
