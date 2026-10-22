"""dt-dependent ``dist_fn`` coefficients (rebuilt on every change of ``code_dt``).

Port of the level-31 coefficient initialisation of GS2 8.2.1 (owner specs
``02_source_term.md`` section 3.1 and ``03_invert_rhs.md`` section 3.3; the critique
``10_critique.md`` A.8 / G-30 confirms both). Line numbers refer to commit 4d8c94bcf:

==========  ==============================================  ===============================
array       Fortran                                         here
==========  ==============================================  ===============================
vpar        ``init_vpar`` (src/dist_fn.fpp:1681-1729)       :func:`init_vpar`
wdrift      ``init_wdrift`` (src/dist_fn.fpp:1405-1495),    :func:`init_wdrift`
wdriftttp   ``wdrift_func`` (1500-1530),                    :func:`wdrift_func`,
            ``wcoriolis_func`` (1534-1557)                  :func:`wcoriolis_func`
wstar       ``init_wstar`` (src/dist_fn.fpp:1732-1750)      :func:`init_wstar`
a, b, r,    ``init_invert_rhs`` (src/dist_fn.fpp:1826-1929) :func:`init_invert_rhs`
ainv
==========  ==============================================  ===============================

``g_h`` (``init_homogeneous_g``) belongs to the parallel-solve module and is returned
as ``None``.  ``ufac`` and ``hyper::D_res`` are *not* rebuilt on a dt change in GS2
(02 section 1.1); they live in :mod:`gs2gpu.source`.

Round-off notes (measured on XLA:CPU, JAX 0.12, see the module tests):

* Every expression keeps the Fortran association order (fully parenthesised in the
  code below).  Products that feed a sum are wrapped in :func:`gs2gpu.numerics.ob`, so
  strict mode (``numerics.set_strict(True)``) forbids FMA contraction exactly where a
  ``-ffp-contract=off`` gfortran build cannot fuse.
* XLA rewrites ``x / y`` into ``x * (1/y)`` when ``y`` is a constant or a broadcast of
  a lower-rank array (e.g. a per-species or per-theta vector divided into a G7 array),
  which differs from the Fortran division in ~30 % of the cases.  All such divisions go
  through :func:`exact_div`, which hides the broadcast divisor behind an optimisation
  barrier (always, not only in strict mode; the divisors here are small).  This also
  holds when the inputs are closed over as constants under ``jax.jit`` (tested); the
  complex denominators handed to ``numerics.recip`` get the same treatment
  (:func:`_divisor`).
* XLA:CPU flushes subnormal numbers to zero (gfortran does not); irrelevant for these
  coefficients.
"""

from __future__ import annotations

from typing import NamedTuple

import jax.numpy as jnp
import numpy as np
from jax import lax

from .io.export import Export
from .numerics import cmul, ob, radd, recip
from .types import Coeffs

# Fortran sentinel of dist_fn_knobs::tpdriftknob (dist_fn.fpp:1286).
TPDRIFTKNOB_SENTINEL = -9.9e9


# ----------------------------------------------------------------------------- helpers
def exact_div(num, den):
    """IEEE ``num / den`` that XLA cannot turn into ``num * (1/den)``.

    The divisor is broadcast to the result shape and passed through
    ``lax.optimization_barrier`` so the algebraic simplifier no longer sees a
    broadcast (or constant) divisor.  Use it for every division whose divisor has a
    lower rank than the result.  (Division by an exact power of two needs no guard:
    ``x/2 == x*0.5`` bit for bit.)
    """
    num = jnp.asarray(num)
    den = jnp.asarray(den)
    shape = jnp.broadcast_shapes(num.shape, den.shape)
    den = lax.optimization_barrier(jnp.broadcast_to(den, shape))
    return jnp.broadcast_to(num, shape) / den


def _complex(re, im):
    """``lax.complex`` with NumPy broadcasting of the two parts."""
    re, im = jnp.broadcast_arrays(jnp.asarray(re), jnp.asarray(im))
    return lax.complex(re, im)


def _divisor(re, im):
    """Complex denominator for ``numerics.recip/cdiv`` with both parts materialised.

    ``cdiv`` divides by the real and imaginary parts of its denominator
    (``bi/br``, ``br/bi``); if one of them is a broadcast of a lower-rank array (e.g.
    the sign-independent real part of ``ainv``'s denominator), XLA rewrites that
    division into a reciprocal multiply (1 ulp off).  The barrier hides the broadcast.
    """
    re, im = jnp.broadcast_arrays(jnp.asarray(re), jnp.asarray(im))
    re, im = lax.optimization_barrier((re, im))
    return lax.complex(re, im)


def _pad_theta(x):
    """Append the ``ig = ntgrid`` entry (0, from ``zero_array``) on the last axis."""
    pad = [(0, 0)] * (x.ndim - 1) + [(0, 1)]
    return jnp.pad(x, pad)


def compress_axes(arr: np.ndarray, axes: tuple[int, ...]) -> np.ndarray:
    """Drop redundant axes (kept with length 1) if ``arr`` is constant along them.

    GS2 stores e.g. ``vpa(ig, isgn, iglo)`` for every ``(it, ik)`` although the values
    do not depend on them (02 section 5).  The compressed array broadcasts back to the
    original bit for bit; if any value differs the array is returned unchanged.
    """
    arr = np.asarray(arr)
    index = tuple(slice(0, 1) if ax in axes else slice(None) for ax in range(arr.ndim))
    small = arr[index]
    if np.array_equal(np.broadcast_to(small, arr.shape), arr):
        return np.ascontiguousarray(small)
    return arr


def _f64(value) -> np.ndarray:
    """0-d float64 array (a traced leaf under jit, never a baked-in constant)."""
    return np.asarray(value, dtype=np.float64)


def _scalar(ex: Export, path: str) -> float:
    value = ex[path]
    return float(np.asarray(value).reshape(()))


def _bool(ex: Export, path: str) -> bool:
    return bool(np.asarray(ex[path]).reshape(()))


def pitch_angle_classes(ex: Export) -> tuple[np.ndarray, np.ndarray]:
    """``(il_is_wfb, il_is_trapped)`` as ``(L,)`` bool: the exported arrays, else their
    definitions ``il == ng2+1`` and ``il > ng2+1`` (le_grids.f90:453-491)."""
    le = ex["/le_grids"]
    if "il_is_wfb" in le and "il_is_trapped" in le:
        return np.asarray(le["il_is_wfb"], bool), np.asarray(le["il_is_trapped"], bool)
    ng2 = int(le["ng2"])
    il = np.arange(1, ex.nlambda + 1)
    return il == ng2 + 1, il > ng2 + 1


def lower_bounce_points(ex: Export) -> np.ndarray:
    """``is_lower_bounce_point(ig, il)`` as ``(L, T)`` bool: the exported array, else the
    elemental function (le_grids.f90:2775-2819; 03 section 3.1)."""
    le = ex["/le_grids"]
    if "is_lower_bounce_point" in le:
        return np.asarray(le["is_lower_bounce_point"], bool)
    bounce = np.asarray(le["is_bounce_point"], bool)
    forbid = np.asarray(le["forbid"], bool)
    il_is_wfb, _ = pitch_angle_classes(ex)
    trapped_wfb = bool(np.asarray(le.get("trapped_wfb", False)).reshape(()))
    lower = np.zeros_like(bounce)
    lower[:, 0] = bounce[:, 0]
    lower[:, 1:-1] = bounce[:, 1:-1] & forbid[:, :-2]
    lower &= ~(il_is_wfb[:, None] & (not trapped_wfb))
    return lower


# ------------------------------------------------------------------------- containers
class CoeffParams(NamedTuple):
    """Static switches of the coefficient build (hashable; ``jit`` static argument).

    ``trapped_grid``: ``le_grids::grid_has_trapped_particles()`` = ``nlambda > ng2``
    (le_grids.f90:447-450); gates the trapped-particle part of ``init_invert_rhs``.
    ``g_shape``: canonical G7 shape ``(S, E, L, Y, X, 2, T)``.
    """

    trapped_grid: bool
    g_shape: tuple

    @classmethod
    def from_export(cls, ex: Export) -> "CoeffParams":
        ng2 = int(ex["/le_grids/ng2"])
        trapped = ex.nlambda > ng2
        if "grid_has_trapped_particles" in ex["/le_grids"]:
            flag = _bool(ex, "/le_grids/grid_has_trapped_particles")
            if flag != trapped:
                raise ValueError("grid_has_trapped_particles disagrees with nlambda > ng2")
        # Gate Q-07 (00_decisions D7): wdrift_func/wcoriolis_func divide by shat for ky = 0.
        shat = _scalar(ex, "/theta_grid/shat")
        if shat == 0.0 and np.any(np.asarray(ex["/kt_grids/aky"]) == 0.0):
            raise NotImplementedError("shat == 0 with an aky == 0 mode is GATED (Q-07)")
        if np.iscomplexobj(ex["/dist_fn"].get("wdrift", np.zeros(0))):
            raise NotImplementedError("LOWFLOW (complex wdrift) is GATED")
        return cls(trapped_grid=bool(trapped), g_shape=tuple(ex.g_shape))


class CoeffInputs(NamedTuple):
    """dt-independent inputs of the coefficient build (NumPy; ``to_jax`` before use).

    Shapes (canonical, spec 09 section 0.2); ``vpa``, ``vpac`` (G7) and ``vperp2``
    (G6) are stored with their redundant ``Y, X`` axes compressed to length 1 when
    the export proves them constant (:func:`compress_axes`).

    * ``vpa, vpac``: ``(S, E, L, 1|Y, 1|X, 2, T)``; ``vperp2``: ``(S, E, L, 1|Y, 1|X, T)``
    * theta vectors ``(T,)``: ``delthet, gradpar, gbdrift, gbdrift0, cvdrift, cvdrift0,
      cdrift, cdrift0``
    * ``aky (Y,)``, ``akx (X,)``, ``theta0 (Y, X)``, ``wunits, tunits (Y,)``
    * species ``(S,)``: ``zstm, stm, tz, fexp, bkdiff``; ``dlnf0drho (S, E)``
    * 0-d: ``shat, mach, tpdriftknob`` (effective value); ``drift_knob (L,)`` is the
      per-pitch-angle ``scaling_knob`` of ``init_wdrift`` (tpdriftknob for wfb and
      trapped ``il``, driftknob otherwise, dist_fn.fpp:1431-1436)
    * masks: ``forbid, is_ttp, is_lower_bounce_point (L, T)``, ``can_be_ttp (L,)``
    """

    vpa: np.ndarray
    vpac: np.ndarray
    vperp2: np.ndarray
    delthet: np.ndarray
    gradpar: np.ndarray
    gbdrift: np.ndarray
    gbdrift0: np.ndarray
    cvdrift: np.ndarray
    cvdrift0: np.ndarray
    cdrift: np.ndarray
    cdrift0: np.ndarray
    shat: np.ndarray
    aky: np.ndarray
    akx: np.ndarray
    theta0: np.ndarray
    wunits: np.ndarray
    tunits: np.ndarray
    zstm: np.ndarray
    stm: np.ndarray
    tz: np.ndarray
    fexp: np.ndarray
    bkdiff: np.ndarray
    dlnf0drho: np.ndarray
    mach: np.ndarray
    tpdriftknob: np.ndarray
    drift_knob: np.ndarray
    forbid: np.ndarray
    is_ttp: np.ndarray
    can_be_ttp: np.ndarray
    is_lower_bounce_point: np.ndarray

    @classmethod
    def from_export(cls, ex: Export) -> "CoeffInputs":
        arrays = ex["/dist_fn_arrays"]
        tg = ex["/theta_grid"]
        kt = ex["/kt_grids"]
        le = ex["/le_grids"]
        sp = ex["/species"]
        df = ex["/dist_fn"]
        gt = ex["/gs2_time"]

        driftknob = float(np.asarray(df["driftknob"]).reshape(()))
        tpdriftknob = float(np.asarray(df["tpdriftknob"]).reshape(()))
        if tpdriftknob == TPDRIFTKNOB_SENTINEL:          # dist_fn.fpp:1286
            tpdriftknob = driftknob
        il_is_wfb, il_is_trapped = pitch_angle_classes(ex)
        # init_wdrift: "Note here we always treat the wfb as trapped" (dist_fn.fpp:1431-1436)
        drift_knob = np.where(il_is_wfb | il_is_trapped, tpdriftknob, driftknob)

        def theta(name):
            return np.asarray(tg[name], np.float64)

        return cls(
            vpa=compress_axes(arrays["vpa"], (3, 4)),
            vpac=compress_axes(arrays["vpac"], (3, 4)),
            vperp2=compress_axes(arrays["vperp2"], (3, 4)),
            delthet=theta("delthet"), gradpar=theta("gradpar"),
            gbdrift=theta("gbdrift"), gbdrift0=theta("gbdrift0"),
            cvdrift=theta("cvdrift"), cvdrift0=theta("cvdrift0"),
            cdrift=theta("cdrift"), cdrift0=theta("cdrift0"),
            shat=_f64(tg["shat"]),
            aky=np.asarray(kt["aky"], np.float64), akx=np.asarray(kt["akx"], np.float64),
            theta0=np.asarray(kt["theta0"], np.float64),
            wunits=np.asarray(gt["wunits"], np.float64),
            tunits=np.asarray(gt["tunits"], np.float64),
            zstm=np.asarray(sp["spec_zstm"], np.float64),
            stm=np.asarray(sp["spec_stm"], np.float64),
            tz=np.asarray(sp["spec_tz"], np.float64),
            fexp=np.atleast_1d(np.asarray(df["fexp"], np.float64)),
            bkdiff=np.atleast_1d(np.asarray(df["bkdiff"], np.float64)),
            dlnf0drho=np.asarray(sp["dlnf0drho"], np.float64),
            mach=_f64(df["mach"]),
            tpdriftknob=_f64(tpdriftknob),
            drift_knob=np.asarray(drift_knob, np.float64),
            forbid=np.asarray(le["forbid"], bool),
            is_ttp=np.asarray(le["is_ttp"], bool),
            can_be_ttp=np.asarray(le["can_be_ttp"], bool),
            is_lower_bounce_point=lower_bounce_points(ex),
        )


# ---------------------------------------------------------------------------- kernels
# Broadcast helpers: canonical G6 axes (S, E, L, Y, X, T), G7 axes (S, E, L, Y, X, 2, T).
def _s6(v):
    return v[:, None, None, None, None, None]


def _s7(v):
    return v[:, None, None, None, None, None, None]


def _y6(v):
    return v[None, None, None, :, None, None]


def _x6(v):
    return v[None, None, None, None, :, None]


def _l6(v):
    return v[None, None, :, None, None, None]


def _lt7(m):
    """``(L, T')`` le_grids mask -> G7 broadcast ``(1, 1, L, 1, 1, 1, T')``."""
    return m[None, None, :, None, None, None, :]


def init_vpar(code_dt, inputs: CoeffInputs, params: CoeffParams):
    """``vpar`` (G7 real); mirrors ``init_vpar`` (src/dist_fn.fpp:1681-1729).

    ``vpar(ig,isgn) = zstm*tunits*code_dt*0.5/delthet(ig)*(|gradpar(ig)|+|gradpar(ig+1)|)
    *vpac(ig,isgn)`` for ``ig < ntgrid`` (left to right, 02 section 3.1); ``vpar(ntgrid) = 0``.
    """
    pref = ((inputs.zstm[:, None] * inputs.tunits[None, :]) * code_dt) * 0.5     # (S, Y)
    fac = exact_div(pref[:, :, None], inputs.delthet[:-1])                       # (S, Y, N)
    fac = fac * (jnp.abs(inputs.gradpar[:-1]) + jnp.abs(inputs.gradpar[1:]))
    vpar = fac[:, None, None, :, None, None, :] * inputs.vpac[..., :-1]
    shape = params.g_shape[:-1] + (params.g_shape[-1] - 1,)
    return _pad_theta(jnp.broadcast_to(vpar, shape))


def wdrift_func(code_dt, inputs: CoeffInputs):
    """Magnetic drift at grid points, G6 ``(S, E, L, Y, X, T)``; mirrors ``wdrift_func``
    (src/dist_fn.fpp:1500-1530), evaluated for every point (the caller masks)::

        aky == 0: akx/shat*(cvdrift0*vpa2 + gbdrift0*0.5*vperp2)*code_dt/2.0
        else:     ((cvdrift + theta0*cvdrift0)*vpa2
                   + (gbdrift + theta0*gbdrift0)*0.5*vperp2)*code_dt*wunits

    with ``vpa2 = vpa(ig,1)**2`` (= ``vpa*vpa``) and the exact ``aky == 0.0`` test.
    """
    vperp2 = inputs.vperp2
    vpa1 = inputs.vpa[..., 0, :]
    vpa2 = vpa1 * vpa1
    th0 = inputs.theta0[None, None, None, :, :, None]
    # zonal branch
    kx_over_shat = _x6(exact_div(inputs.akx, inputs.shat))
    zonal = ob(inputs.cvdrift0 * vpa2) + ob((inputs.gbdrift0 * 0.5) * vperp2)
    zonal = ((kx_over_shat * zonal) * code_dt) / 2.0
    # finite-ky branch
    cv = inputs.cvdrift + ob(th0 * inputs.cvdrift0)
    gb = inputs.gbdrift + ob(th0 * inputs.gbdrift0)
    finite = ob(cv * vpa2) + ob((gb * 0.5) * vperp2)
    finite = (finite * code_dt) * _y6(inputs.wunits)
    return jnp.where(_y6(inputs.aky) == 0.0, zonal, finite)


def wcoriolis_func(code_dt, inputs: CoeffInputs):
    """Coriolis drift at grid points, G6; mirrors ``wcoriolis_func``
    (src/dist_fn.fpp:1534-1557)::

        aky == 0: mach*vpa(ig,1)*cdrift0*code_dt*akx/(2.*shat*stm)
        else:     mach*vpa(ig,1)*(cdrift + theta0*cdrift0)*code_dt*wunits/stm
    """
    vpa1 = inputs.vpa[..., 0, :]
    mv = inputs.mach * vpa1
    th0 = inputs.theta0[None, None, None, :, :, None]
    num0 = ((mv * inputs.cdrift0) * code_dt) * _x6(inputs.akx)
    zonal = exact_div(num0, _s6((2.0 * inputs.shat) * inputs.stm))
    cd = inputs.cdrift + ob(th0 * inputs.cdrift0)
    num1 = ((mv * cd) * code_dt) * _y6(inputs.wunits)
    finite = exact_div(num1, _s6(inputs.stm))
    return jnp.where(_y6(inputs.aky) == 0.0, zonal, finite)


def init_wdrift(code_dt, inputs: CoeffInputs, params: CoeffParams):
    """``(wdrift, wdriftttp)``, both G7 real; mirrors ``init_wdrift``
    (src/dist_fn.fpp:1405-1495).

    Pass 1 (grid values, skipped at forbidden points, which stay 0)::

        wdrift(ig,1) = wd*scaling_knob + wc ;  wdrift(ig,2) = wd*scaling_knob - wc

    pass 2 (``can_be_ttp(il)`` and ``is_ttp(ig,il)`` only, all ``ig``)::

        wdriftttp(ig,1) = wd*tpdriftknob + wc ;  wdriftttp(ig,2) = wd*tpdriftknob - wc

    pass 3: ``wdrift(ig) = 0.5*(wdrift(ig) + wdrift(ig+1))`` for ``ig < ntgrid`` in
    ascending order, i.e. from the *uncentred* ``ig+1`` (02 section 3.1, pitfall 6);
    ``wdrift(ntgrid)`` keeps its grid value; ``wdriftttp`` is not centred.
    """
    wd = wdrift_func(code_dt, inputs)
    wc = wcoriolis_func(code_dt, inputs)
    shape = params.g_shape
    # GS2 quirk: the Coriolis drift is not scaled by driftknob/tpdriftknob ("Is there a
    # reason we don't scale this by scaling_knob?", dist_fn.fpp:1442-1443, 1473-1474).
    wdk = ob(wd * _l6(inputs.drift_knob))
    grid = jnp.stack([wdk + wc, wdk - wc], axis=-2)
    forbid = inputs.forbid[None, None, :, None, None, None, :]
    grid = jnp.where(forbid, 0.0, jnp.broadcast_to(grid, shape))
    wdt = ob(wd * inputs.tpdriftknob)
    ttp_vals = jnp.stack([wdt + wc, wdt - wc], axis=-2)
    ttp = (inputs.can_be_ttp[:, None] & inputs.is_ttp)[None, None, :, None, None, None, :]
    wdriftttp = jnp.where(ttp, jnp.broadcast_to(ttp_vals, shape), 0.0)
    # GS2 quirk: plain 0.5 average, not bakdif-weighted ("This should be weighted by
    # bakdif", dist_fn.fpp:1481); forbidden grid values enter as 0; wdriftttp uncentred.
    centred = 0.5 * (grid[..., :-1] + grid[..., 1:])
    wdrift = jnp.concatenate([centred, grid[..., -1:]], axis=-1)
    return wdrift, wdriftttp


def init_wstar(code_dt, inputs: CoeffInputs):
    """``wstar`` with NumPy shape ``(S, E, Y)`` (Fortran ``wstar(ik, ie, is)``); mirrors
    ``init_wstar`` (src/dist_fn.fpp:1732-1750): ``- code_dt*wunits(ik)*dlnf0drho(ie,is)``.
    """
    return -((code_dt * inputs.wunits[None, None, :]) * inputs.dlnf0drho[:, :, None])


def init_invert_rhs(vpar, wdrift, wdriftttp, inputs: CoeffInputs, params: CoeffParams):
    """``(a, b, r, ainv)``, G7 complex; mirrors ``init_invert_rhs``
    (src/dist_fn.fpp:1826-1929, without the ``init_homogeneous_g`` call).

    For ``ig < ntgrid``, both ``isgn``, with ``wd = wdrift(ig,isgn)``,
    ``vp = vpar(ig,1)`` (sigma = 1 for both signs) and the GCC component lowering of
    03 section 3.3 (``zi*wd = (0, wd)``)::

        ainv = 1.0/(1.0 + bd + (1.0-fexp)*tz*(zi*wd*(1.0+bd) + 2.0*vp))
        r    = (1.0 - bd + (1.0-fexp)*tz*(zi*wd*(1.0-bd) - 2.0*vp))*ainv
        a    = 1.0 + bd + fexp*tz*(-zi*wd*(1.0+bd) - 2.0*vp)
        b    = 1.0 - bd + fexp*tz*(-zi*wd*(1.0-bd) + 2.0*vp)

    then, only on a trapped grid: ``r = ainv = 0`` where ``forbid(ig) .or. forbid(ig+1)``;
    ``ainv = 1.0 + ainv`` for ``isgn == 1`` where ``is_lower_bounce_point(ig+1)``; at
    ``is_ttp(ig)``: ``ainv = 1/(1 + zi*(1-fexp)*tz*wdttp)``, ``a = 1 - zi*fexp*tz*wdttp``,
    ``r = 0``.  Entries ``ig = ntgrid`` are 0.  The complex division is GCC's Smith
    method (``numerics.recip``).
    """
    wd = wdrift[..., :-1]
    vp = vpar[..., 0:1, :-1]
    bd = _s7(inputs.bkdiff)
    fe = _s7(inputs.fexp)
    tz = _s7(inputs.tz)
    opb = 1.0 + bd
    omb = 1.0 - bd
    c1 = (1.0 - fe) * tz
    c2 = fe * tz
    v2 = 2.0 * vp
    wp = wd * opb
    wm = wd * omb
    c1v2 = ob(c1 * v2)
    c2v2 = ob(c2 * v2)
    ainv = recip(_divisor(opb + c1v2, c1 * wp))
    r = cmul(_complex(omb - c1v2, c1 * wm), ainv)
    a = _complex(opb - c2v2, -(c2 * wp))
    b = _complex(omb + c2v2, -(c2 * wm))
    if params.trapped_grid:
        forbid = inputs.forbid
        fz = _lt7(forbid[:, :-1] | forbid[:, 1:])
        r = jnp.where(fz, 0.0 + 0.0j, r)                                  # :1901-1904
        ainv = jnp.where(fz, 0.0 + 0.0j, ainv)
        sigma1 = jnp.asarray([True, False])[:, None]
        lb = _lt7(inputs.is_lower_bounce_point[:, 1:]) & sigma1           # :1910-1914
        ainv = jnp.where(lb, radd(1.0, ainv), ainv)
        ttp = _lt7(inputs.is_ttp[:, :-1])                                 # :1917-1921
        wdttp = wdriftttp[..., :-1]
        ones = jnp.ones(jnp.broadcast_shapes(c1.shape, wdttp.shape))
        ainv_ttp = recip(_divisor(ones, c1 * wdttp))
        a_ttp = _complex(ones, -(c2 * wdttp))
        ainv = jnp.where(ttp, ainv_ttp, ainv)
        a = jnp.where(ttp, a_ttp, a)
        r = jnp.where(ttp, 0.0 + 0.0j, r)
    shape = params.g_shape[:-1] + (params.g_shape[-1] - 1,)
    out = []
    for x in (a, b, r, ainv):
        out.append(_pad_theta(jnp.broadcast_to(x, shape)))
    return tuple(out)


def compute_coeffs(code_dt, inputs: CoeffInputs, params: CoeffParams) -> Coeffs:
    """All dt-dependent ``dist_fn`` coefficients at ``code_dt`` (``g_h = None``).

    Mirrors the dt-dependent part of ``init_dist_fn_level_3`` (src/dist_fn.fpp:960-1013):
    ``init_vpar`` (1681), ``init_wdrift`` (1405), ``init_wstar`` (1732) and
    ``init_invert_rhs`` (1826) in that order.  ``jax.jit``-compatible with ``params``
    static; ``code_dt`` may be a Python float or a traced scalar.
    """
    dt = jnp.asarray(code_dt, dtype=jnp.float64)
    vpar = init_vpar(dt, inputs, params)
    wdrift, wdriftttp = init_wdrift(dt, inputs, params)
    wstar = init_wstar(dt, inputs)
    a, b, r, ainv = init_invert_rhs(vpar, wdrift, wdriftttp, inputs, params)
    return Coeffs(code_dt=code_dt, vpar=vpar, wdrift=wdrift, wdriftttp=wdriftttp,
                  wstar=wstar, a=a, b=b, r=r, ainv=ainv, g_h=None)
