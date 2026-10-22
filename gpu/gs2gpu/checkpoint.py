"""Checkpoint and restart of a JAX run.

A checkpoint holds the complete :class:`gs2gpu.step.State` (distribution
functions, fields, explicit-term history, the time block, the flow-shear
bookkeeping and every SAVE counter of 00_decisions D7/Q-10), the step it was
taken after, and the diagnostics state, so a restarted run continues exactly
as if it had never stopped. Arrays are stored as NumPy (host) copies; the file
is a pickle and must only be loaded from trusted locations (your own runs).
"""

from __future__ import annotations

import dataclasses
import os
import pickle
from pathlib import Path
from typing import Any, Optional

import jax
import jax.numpy as jnp
import numpy as np

from .step import State

FORMAT = "gs2gpu-checkpoint-1"


def _to_host(tree):
    return jax.tree_util.tree_map(
        lambda x: np.asarray(x) if isinstance(x, jax.Array) else x, tree)


def _to_device(tree):
    return jax.tree_util.tree_map(
        lambda x: jnp.asarray(x) if isinstance(x, np.ndarray) else x, tree)


def save(path: str | Path, state: State, istep: int, *, export_path: str | Path,
         diag_state: Any = None, extra: Optional[dict] = None) -> Path:
    """Write a checkpoint atomically (``path`` is replaced only once complete)."""
    path = Path(path)
    payload = {
        "format": FORMAT,
        "istep": int(istep),
        "export_path": str(export_path),
        "state": {f.name: _to_host(getattr(state, f.name)) for f in dataclasses.fields(state)},
        "diag_state": _to_host(diag_state),
        "extra": extra or {},
    }
    tmp = path.with_suffix(path.suffix + ".tmp")
    with open(tmp, "wb") as fh:
        pickle.dump(payload, fh, protocol=pickle.HIGHEST_PROTOCOL)
    os.replace(tmp, path)
    return path


def load(path: str | Path) -> dict:
    """Read a checkpoint written by :func:`save`; returns a dict with ``state`` (a
    :class:`State` on device), ``istep``, ``export_path``, ``diag_state``, ``extra``."""
    with open(path, "rb") as fh:
        payload = pickle.load(fh)
    if payload.get("format") != FORMAT:
        raise ValueError(f"{path}: not a {FORMAT} checkpoint")
    fields = _to_device(payload["state"])
    payload["state"] = State(**fields)
    payload["diag_state"] = _to_device(payload["diag_state"])
    return payload
