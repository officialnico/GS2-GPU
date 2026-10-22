"""Shared containers passed between the kernels of the GS2 JAX port.

Every container is built from the export file (``gs2gpu.io.export``) with a
``from_export`` constructor and holds NumPy arrays in the canonical layout of
spec 09 section 0.2 (NumPy C order; Fortran order reversed):

=====  ======================================  =============================
tag    NumPy shape                             Fortran
=====  ======================================  =============================
G7     ``(S, E, L, Y, X, 2, T)``               ``g(ig, isgn, it, ik, il, ie, is)``
G6     ``(S, E, L, Y, X, T)``                  ``x(ig, it, ik, il, ie, is)``
F3     ``(Y, X, T)``                           ``phi(ig, it, ik)``
K2     ``(Y, X)``                              ``x(it, ik)``
LT     ``(L, T)``                              ``wl(ig, il)``, ``forbid(ig, il)``
=====  ======================================  =============================

with S = nspec, E = negrid, L = nlambda, Y = naky, X = ntheta0, T = 2*ntgrid+1
(theta index 0 is Fortran ``ig = -ntgrid``), sign index 0/1 = Fortran isgn 1/2.

Containers are ``NamedTuple``s, which JAX treats as pytrees; ``to_jax``
moves every array leaf to the default device. Static sizes and switches live
in ``Dims`` (a frozen dataclass, hashable, safe as a ``jit`` static argument).
Arrays that GS2 did not allocate are ``None``.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import NamedTuple, Optional

import jax
import jax.numpy as jnp
import numpy as np

from .io.export import Export


def to_jax(tree):
    """Copy every NumPy leaf of ``tree`` to the default JAX device."""
    return jax.tree_util.tree_map(
        lambda x: jnp.asarray(x) if isinstance(x, np.ndarray) else x, tree)


def _opt(ex: Export, path: str):
    """``ex[path]`` or None if GS2 did not allocate it (listed as absent)."""
    group, name = path.rsplit("/", 1)
    grp = ex[group] if group.strip("/") else ex.root
    if name in grp:
        return grp[name]
    if grp.is_absent(name):
        return None
    raise KeyError(f"{path} missing from export and not listed as absent")


def _bool(value) -> bool:
    return bool(np.asarray(value).item()) if value is not None else False


@dataclass(frozen=True)
class Dims:
    """Static sizes and a few structural switches (hashable)."""

    ntgrid: int
    nth: int
    ntheta0: int
    naky: int
    nlambda: int
    negrid: int
    nspec: int
    ng2: int
    nxi: int
    box: bool
    nonlin: bool
    has_phi: bool
    has_apar: bool
    has_bpar: bool

    @property
    def g_shape(self) -> tuple[int, ...]:
        return (self.nspec, self.negrid, self.nlambda, self.naky, self.ntheta0, 2, self.nth)

    @property
    def field_shape(self) -> tuple[int, int, int]:
        return (self.naky, self.ntheta0, self.nth)

    @classmethod
    def from_export(cls, ex: Export) -> "Dims":
        return cls(
            ntgrid=ex.ntgrid, nth=ex.nth, ntheta0=ex.ntheta0, naky=ex.naky,
            nlambda=ex.nlambda, negrid=ex.negrid, nspec=ex.nspec,
            ng2=int(ex["/le_grids/ng2"]), nxi=int(ex["/le_grids/nxi"]),
            box=_bool(ex["/kt_grids/box"]),
            nonlin=_bool(ex["/nonlinear_terms/nonlin"]),
            has_phi=_bool(ex["/run_parameters/has_phi"]),
            has_apar=_bool(ex["/run_parameters/has_apar"]),
            has_bpar=_bool(ex["/run_parameters/has_bpar"]),
        )


class Geometry(NamedTuple):
    """``/theta_grid`` arrays, all shape ``(T,)`` unless noted."""

    theta: np.ndarray
    delthet: np.ndarray
    bmag: np.ndarray
    gradpar: np.ndarray
    itor_over_B: np.ndarray
    gbdrift: np.ndarray
    gbdrift0: np.ndarray
    cvdrift: np.ndarray
    cvdrift0: np.ndarray
    cdrift: np.ndarray
    cdrift0: np.ndarray
    gds2: np.ndarray
    gds21: np.ndarray
    gds22: np.ndarray
    grho: np.ndarray
    jacob: np.ndarray
    shat: float
    kxfac: float
    qval: float
    drhodpsi: float
    field_line_average_weight: float

    @classmethod
    def from_export(cls, ex: Export) -> "Geometry":
        g = ex["/theta_grid"]
        return cls(**{name: g[name] for name in cls._fields})


class KGrids(NamedTuple):
    """``/kt_grids`` arrays: ``aky (Y,)``, ``akx (X,)``, ``theta0`` K2,
    ``kperp2``/``inv_kperp2`` F3, ``kwork_filter`` K2 bool (value after init_fields)."""

    aky: np.ndarray
    akx: np.ndarray
    ikx: np.ndarray
    theta0: np.ndarray
    kperp2: np.ndarray
    inv_kperp2: np.ndarray
    kwork_filter: np.ndarray

    @classmethod
    def from_export(cls, ex: Export) -> "KGrids":
        g = ex["/kt_grids"]
        return cls(**{name: g[name] for name in cls._fields})


class VGrids(NamedTuple):
    """``/le_grids`` velocity grids and masks.

    ``energy, w, speed``: ``(S, E)``; ``al``: ``(L,)``; ``wl``, ``forbid``,
    ``is_ttp``, ``is_bounce_point``, ``is_lower_bounce_point``,
    ``is_upper_bounce_point``: ``(L, T)``; ``jend``: ``(T,)`` int (Fortran
    values); ``can_be_ttp``, ``il_is_wfb``, ``il_is_passing``,
    ``il_is_trapped``: ``(L,)`` bool; ``is_passing_hybrid_electron``:
    ``(S, Y, L)`` bool.
    """

    energy: np.ndarray
    w: np.ndarray
    speed: np.ndarray
    al: np.ndarray
    wl: np.ndarray
    jend: np.ndarray
    forbid: np.ndarray
    is_ttp: np.ndarray
    is_bounce_point: np.ndarray
    is_lower_bounce_point: np.ndarray
    is_upper_bounce_point: np.ndarray
    can_be_ttp: np.ndarray
    il_is_wfb: np.ndarray
    il_is_passing: np.ndarray
    il_is_trapped: np.ndarray
    is_passing_hybrid_electron: np.ndarray
    passing_wfb: bool
    trapped_wfb: bool
    mixed_wfb: bool

    @classmethod
    def from_export(cls, ex: Export) -> "VGrids":
        g = ex["/le_grids"]
        kw = {name: g[name] for name in cls._fields}
        for flag in ("passing_wfb", "trapped_wfb", "mixed_wfb"):
            kw[flag] = _bool(kw[flag])
        return cls(**kw)


class Species(NamedTuple):
    """``/species``: per-species vectors ``(S,)`` (``spec_<component>``) and
    ``nonmaxw_corr``, ``f0_values``, ``dlnf0drho`` of shape ``(S, E)``."""

    z: np.ndarray
    mass: np.ndarray
    dens: np.ndarray
    temp: np.ndarray
    tprim: np.ndarray
    fprim: np.ndarray
    uprim: np.ndarray
    vnewk: np.ndarray
    stm: np.ndarray
    zstm: np.ndarray
    tz: np.ndarray
    smz: np.ndarray
    zt: np.ndarray
    nu_h: np.ndarray
    bess_fac: np.ndarray
    type: np.ndarray
    nonmaxw_corr: np.ndarray
    f0_values: np.ndarray
    dlnf0drho: np.ndarray

    @classmethod
    def from_export(cls, ex: Export) -> "Species":
        g = ex["/species"]
        kw = {}
        for name in cls._fields:
            kw[name] = g[name] if name in ("nonmaxw_corr", "f0_values", "dlnf0drho") \
                else g[f"spec_{name}"]
        return cls(**kw)


class PointArrays(NamedTuple):
    """dt-independent per-point arrays: ``vpa, vpac`` G7 real;
    ``vperp2, aj0, aj1`` G6 real (``/dist_fn_arrays``)."""

    vpa: np.ndarray
    vpac: np.ndarray
    vperp2: np.ndarray
    aj0: np.ndarray
    aj1: np.ndarray

    @classmethod
    def from_export(cls, ex: Export) -> "PointArrays":
        g = ex["/dist_fn_arrays"]
        return cls(**{name: g[name] for name in cls._fields})


class Coeffs(NamedTuple):
    """dt-dependent ``dist_fn`` coefficients (rebuilt on every dt change).

    ``vpar, wdrift, wdriftttp``: G7 real; ``a, b, r, ainv, g_h``: G7 complex;
    ``wstar``: ``(S, E, Y)`` real (Fortran ``wstar(ik, ie, is)``). ``code_dt``
    is the time step they were built for. ``g_h`` is None unless allocated.
    """

    code_dt: float
    vpar: np.ndarray
    wdrift: np.ndarray
    wdriftttp: np.ndarray
    wstar: np.ndarray
    a: np.ndarray
    b: np.ndarray
    r: np.ndarray
    ainv: np.ndarray
    g_h: Optional[np.ndarray]

    @classmethod
    def from_export(cls, ex: Export, rung: Optional[int] = None) -> "Coeffs":
        """Coefficients at the initial dt, or at dt-ladder rung ``rung`` (1-based)."""
        if rung is None:
            dist_fn = ex["/dist_fn"]
            arrays = ex["/dist_fn_arrays"]
            code_dt = float(ex["/state_t0/code_dt"])
        else:
            base = ex[f"/dt_ladder/k{rung}"]
            dist_fn = base["dist_fn"]
            arrays = base["dist_fn_arrays"]
            code_dt = float(base.attrs["code_dt"])
        return cls(
            code_dt=code_dt,
            vpar=dist_fn["vpar"], wdrift=dist_fn["wdrift"], wdriftttp=dist_fn["wdriftttp"],
            wstar=arrays["wstar"], a=dist_fn["a"], b=dist_fn["b"], r=dist_fn["r"],
            ainv=dist_fn["ainv"], g_h=dist_fn.get("g_h"),
        )


class Fields(NamedTuple):
    """Potentials, each F3 complex ``(Y, X, T)``; absent fields are zeros."""

    phi: jnp.ndarray
    apar: jnp.ndarray
    bpar: jnp.ndarray


class TimeBlock(NamedTuple):
    """The ``gs2_time`` scalars (spec 09 section 7.6), host float64."""

    code_dt: float
    user_dt: float
    code_dt_prev1: float
    code_dt_prev2: float
    code_dt_old: float
    code_time: float
    user_time: float
    code_dt_cfl: float
    user_dt_cfl: float
    code_dt_min: float
    code_dt_max: float

    @classmethod
    def from_group(cls, grp) -> "TimeBlock":
        return cls(**{name: float(grp[name]) for name in cls._fields})


def fields_from_group(grp, dims: Dims, new: bool = True) -> Fields:
    """``Fields`` from a group holding ``phinew/aparnew/bparnew`` (or
    ``phi/apar/bpar`` with ``new=False``); missing fields become zeros."""
    suffix = "new" if new else ""
    out = []
    for name in ("phi", "apar", "bpar"):
        key = f"{name}{suffix}"
        value = grp.get(key)
        out.append(np.zeros(dims.field_shape, np.complex128) if value is None else value)
    return Fields(*out)
