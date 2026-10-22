"""Time the JAX step on the current device: ``gs2gpu bench <run>.gs2export.nc``.

Runs a few warm-up steps (compilation), then times ``steps`` steps of the main
loop without diagnostics and reports milliseconds per step, the size of one
distribution function and, on GPUs, the device memory in use.
"""

from __future__ import annotations

import time
from dataclasses import dataclass

import jax
import numpy as np

from .driver import DtCache, _do_reset
from .io.export import load_export
from .reset import update_time
from .step import Model, advance


@dataclass
class BenchResult:
    device: str
    g_shape: tuple
    g_bytes: int
    warmup_s: float
    ms_per_step: float
    steps: int
    resets: int
    peak_bytes: int | None


def bench(export_path, steps: int = 50, warmup: int = 3) -> BenchResult:
    ex = load_export(export_path)
    model = Model(ex)
    cache = DtCache(model)
    dts = cache.initial
    S = model.initial_state()
    resets = 0

    def one(S, dts, istep):
        nonlocal resets
        S, reset = advance(model, dts, S, istep)
        if reset:
            S, new_dts, exit_, _ = _do_reset(model, cache, S, istep)
            S = S.replace(reset=False)
            if exit_:
                raise RuntimeError(f"run exits at step {istep} (time step limits)")
            dts = new_dts
            resets += 1
            S, _ = advance(model, dts, S, istep)
        return S.replace(time=update_time(S.time)), dts

    t0 = time.perf_counter()
    istep = 0
    for _ in range(warmup):
        istep += 1
        S, dts = one(S, dts, istep)
    jax.block_until_ready(S.gnew)
    warm = time.perf_counter() - t0
    t0 = time.perf_counter()
    for _ in range(steps):
        istep += 1
        S, dts = one(S, dts, istep)
    jax.block_until_ready(S.gnew)
    elapsed = time.perf_counter() - t0
    dev = jax.devices()[0]
    peak = None
    try:
        stats = dev.memory_stats() or {}
        peak = stats.get("peak_bytes_in_use")
    except Exception:  # CPU devices have no memory stats
        pass
    shape = model.dims.g_shape
    return BenchResult(device=f"{dev.platform}:{dev.device_kind}", g_shape=shape,
                       g_bytes=int(np.prod(shape)) * 16, warmup_s=warm,
                       ms_per_step=1e3 * elapsed / steps, steps=steps, resets=resets,
                       peak_bytes=peak)
