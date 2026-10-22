"""Tests of gs2gpu.coeffs (dt-dependent dist_fn coefficients).

Reference-data tests (skip without ``gpu/refdata``): ``compute_coeffs`` against the
exported coefficients at the initial dt (``/dist_fn/*``, ``/dist_fn_arrays/wstar``) and at
every dt-ladder rung (``types.Coeffs.from_export(ex, rung=k)``), for every case, at
``RTOL_KERNEL``; the same in strict mode (no FMA, should be bitwise with the
``-ffp-contract=off`` reference build).

Synthetic tests (always run): a small synthetic export with physically consistent
le_grids masks (two magnetic wells -> ttp points, bounce points, forbidden regions) is
built in memory; ``compute_coeffs`` (strict mode, jitted) must equal *bit for bit* a
literal scalar transcription of the Fortran loops (``ref_coeffs`` below, plain Python
floats, Fortran association order; complex division by GCC's Smith method).  The
synthetic builder and the scalar reference are reused by ``test_source.py``.
"""

from __future__ import annotations

import contextlib

from pathlib import Path

import jax
import jax.numpy as jnp
import numpy as np
import pytest

from conftest import RTOL_KERNEL, assert_close, cases_where, export_of
from gs2gpu import coeffs as co
from gs2gpu import numerics
from gs2gpu.io.export import Export, Group
from gs2gpu.types import Coeffs, to_jax

COEFF_NAMES = ("vpar", "wdrift", "wdriftttp", "wstar", "a", "b", "r", "ainv")
EPS = float(np.finfo(np.float64).eps)


# =============================================================================== utils
@contextlib.contextmanager
def strict_mode(flag: bool = True):
    """Switch numerics strict mode (no FMA contraction) for code traced inside."""
    old = numerics.strict()
    numerics.set_strict(flag)
    try:
        yield
    finally:
        numerics.set_strict(old)


def jit_coeffs(params):
    """A fresh jitted ``compute_coeffs`` (new trace, so the current strict flag applies)."""
    return jax.jit(lambda dt, inputs: co.compute_coeffs(dt, inputs, params))


def bitwise_equal(ref, got) -> bool:
    """Exact equality of values (``-0.0 == 0.0``), real and imaginary parts."""
    ref = np.asarray(ref)
    got = np.asarray(got)
    return ref.shape == got.shape and bool(np.all(ref == got))


def strict_rtol(ex) -> float:
    """Tolerance of the strict-mode tests: bitwise (0) against the ``-ffp-contract=off``
    reference build (00_decisions D2; measured exact), else 1e-15."""
    return 0.0 if "-ffp-contract=off" in str(ex.attrs.get("build_config", "")) else 1e-15


def ladder_rungs(ex) -> list[int]:
    if "dt_ladder" not in ex.root.groups:
        return []
    out = []
    for name in ex["/dt_ladder"].groups:
        if name.startswith("k") and name[1:].isdigit():
            out.append(int(name[1:]))
    return sorted(out)


# ===================================================================== synthetic export
def _group(path: str, arrays: dict, groups: dict | None = None, attrs: dict | None = None):
    name = path.rsplit("/", 1)[-1] or "/"
    return Group(name=name, path=path, arrays=dict(arrays), groups=dict(groups or {}),
                 attrs=dict(attrs or {}))


def le_masks(al, bmag, ng2, ntgrid, trapped_wfb=False):
    """forbid, is_bounce_point, is_lower/upper_bounce_point, can_be_ttp, is_ttp (L, T)
    from their definitions (le_grids.f90:2715-2946; 03 section 3.1), Fortran loops."""
    nl, nt = len(al), len(bmag)
    fuzz = 10 * EPS
    forbid = np.zeros((nl, nt), bool)
    bounce = np.zeros((nl, nt), bool)
    for il in range(nl):
        for ig in range(nt):
            forbid[il, ig] = 1.0 - al[il] * bmag[ig] < -fuzz
            if il + 1 >= ng2 + 2:
                bounce[il, ig] = abs(1.0 - al[il] * bmag[ig]) <= fuzz
    lower = np.zeros_like(forbid)
    upper = np.zeros_like(forbid)
    for il in range(nl):
        is_wfb = il + 1 == ng2 + 1
        for ig in range(nt):
            if not bounce[il, ig] or (is_wfb and not trapped_wfb):
                continue
            if ig == 0:
                lower[il, ig], upper[il, ig] = True, False
            elif ig == nt - 1:
                lower[il, ig], upper[il, ig] = False, True
            else:
                lower[il, ig] = forbid[il, ig - 1]
                upper[il, ig] = forbid[il, ig + 1]
    can_be_ttp = np.zeros(nl, bool)
    is_ttp = np.zeros((nl, nt), bool)
    if nl > ng2:
        ittp = np.full(nt, nl + 1)                       # Fortran il values
        for ig in range(1, nt - 1):
            for il in range(ng2 + 1, nl + 1):
                if forbid[il - 1, ig - 1] and forbid[il - 1, ig + 1] and not forbid[il - 1, ig]:
                    ittp[ig] = il
                    break
        for il in range(ng2 + 1, nl + 1):
            can_be_ttp[il - 1] = bool(np.any(il >= ittp))
            for ig in range(nt):
                if il >= ittp[ig] and not forbid[il - 1, ig]:
                    is_ttp[il - 1, ig] = True
    return forbid, bounce, lower, upper, can_be_ttp, is_ttp


def make_synthetic_export(seed: int = 0, *, nspec: int = 2, negrid: int = 3, ntgrid: int = 4,
                          naky: int = 3, ntheta0: int = 4, ng2: int = 3, trapped: bool = True,
                          zonal: bool = True, box: bool = True, nonlin: bool = False,
                          has_apar: bool = True, has_bpar: bool = True, mach: float = 0.07,
                          g_exb: float = 0.05, driftknob: float = 1.1,
                          tpdriftknob: float = 0.9, source_option: int = 1,
                          nonmaxwell: bool = True, code_dt: float = 0.05,
                          ladder: tuple = (0.025,), bkdiff: tuple = (0.05, 0.05),
                          fexp: tuple = (0.45, 0.4)) -> Export:
    """In-memory Export with every array coeffs/source read (spec 09 names).

    Geometry: a double-well ``bmag`` (local maximum at theta = 0, minima at +-pi/2), so a
    trapped grid has ttp points, bounce points and separate allowed wells.  Velocity
    arrays (``vpa, vpac, vperp2``) follow init_vpa_vpac / init_vperp2 and are stored
    with their redundant ky/kx axes, as GS2 does.  The dt-dependent coefficient groups
    (``/dist_fn``, ``/dist_fn_arrays/wstar``, ``/dt_ladder/k<k>``) are filled with the
    scalar reference :func:`ref_coeffs`.
    """
    rng = np.random.default_rng(seed)
    nth = 2 * ntgrid + 1
    S, E, Y, X, T = nspec, negrid, naky, ntheta0, nth
    tu = np.linspace(-np.pi, np.pi, nth)
    theta = tu + 0.03 * np.sin(3.0 * tu)
    bmag = 1.0 - 0.2 * np.cos(tu) + 0.15 * np.cos(2.0 * tu)
    bmag = 0.5 * (bmag + bmag[::-1])                          # exactly symmetric
    bmax = bmag.max()
    al = list(np.linspace(0.08, 0.9, ng2) / bmax)
    if trapped:
        bset = np.unique(np.round(bmag[ntgrid:], 14))[::-1]   # bmax first (the wfb)
        al += [1.0 / b for b in bset]
    al = np.asarray(al)
    L = len(al)
    lmax = L - 1 if trapped else L
    forbid, bounce, lower, upper, can_be_ttp, is_ttp = le_masks(al, bmag, ng2, ntgrid)
    il1 = np.arange(1, L + 1)
    energy = np.outer(1.0 + 0.1 * np.arange(S), np.linspace(0.4, 3.1, E))
    # init_vpa_vpac / init_vperp2 (dist_fn.fpp:1560-1678), (S, E, L, T)
    vpa1 = np.sqrt(energy[:, :, None, None] * np.maximum(0.0, 1.0 - al[None, None, :, None]
                                                          * bmag[None, None, None, :]))
    vpa1 = np.where(forbid[None, None], 0.0, vpa1)
    vpa = np.stack([vpa1, -vpa1], axis=-2)                    # (S, E, L, 2, T)
    vpac = np.zeros_like(vpa)
    fz = forbid[:, :-1] | forbid[:, 1:]
    vpac[..., :-1] = np.where(fz[None, None, :, None, :], 0.0,
                              0.5 * (vpa[..., :-1] + vpa[..., 1:]))
    vperp2 = np.where(~forbid[None, None], (bmag[None, None, None, :] * al[None, None, :, None])
                      * energy[:, :, None, None], 0.0)

    def g7(x):          # (S, E, L, 2, T) -> G7 with redundant ky, kx axes
        return np.ascontiguousarray(np.broadcast_to(x[:, :, :, None, None], (S, E, L, Y, X, 2, T)))

    def g6(x):          # (S, E, L, T) -> G6
        return np.ascontiguousarray(np.broadcast_to(x[:, :, :, None, None], (S, E, L, Y, X, T)))

    aky = np.array([0.0, 0.35, 0.7][:Y]) if zonal else np.linspace(0.2, 0.7, Y)
    akx = np.concatenate([[0.0], np.sort(rng.uniform(-1.0, 1.0, X - 1))]) if box \
        else rng.uniform(-1.0, 1.0, X)
    shat = 0.8
    theta0 = rng.uniform(-1.5, 1.5, (Y, X))
    delthet = np.concatenate([np.diff(theta), [0.0]])
    kwork = np.zeros((Y, X), bool)
    if box:
        kwork[0, 0] = True

    def smooth(scale=1.0, offset=0.0):
        c = rng.normal(size=3)
        return offset + scale * (c[0] + c[1] * np.cos(tu) + c[2] * np.sin(2.0 * tu))

    tg = _group("/theta_grid", {
        "theta": theta, "delthet": delthet, "bmag": bmag,
        "gradpar": 0.4 + 0.05 * np.cos(tu) + 0.01 * rng.normal(size=T),
        "gbdrift": smooth(), "gbdrift0": smooth(0.5), "cvdrift": smooth(),
        "cvdrift0": smooth(0.5), "cdrift": smooth(0.3), "cdrift0": smooth(0.2),
        "itor_over_B": smooth(0.2, 1.0), "shat": shat,
    })
    kt = _group("/kt_grids", {"aky": aky, "akx": akx, "theta0": theta0, "kwork_filter": kwork,
                              "box": box})
    le = _group("/le_grids", {
        "al": al, "ng2": ng2, "nlambda": L, "lmax": lmax, "energy": energy,
        "forbid": forbid, "is_bounce_point": bounce, "is_lower_bounce_point": lower,
        "is_upper_bounce_point": upper, "can_be_ttp": can_be_ttp, "is_ttp": is_ttp,
        "il_is_wfb": il1 == ng2 + 1, "il_is_trapped": il1 > ng2 + 1,
        "il_is_passing": il1 <= ng2, "grid_has_trapped_particles": L > ng2,
        "nxi": max(2 * L - 1, 2 * ng2),
    })
    zstm = rng.uniform(0.5, 1.5, S)
    stm = rng.uniform(0.5, 1.5, S)
    tz = rng.uniform(0.5, 1.5, S)
    nmc = rng.uniform(0.8, 1.2, (S, E)) if nonmaxwell else np.ones((S, E))
    sp = _group("/species", {"spec_zstm": zstm, "spec_stm": stm, "spec_tz": tz,
                             "spec_zt": 1.0 / tz, "dlnf0drho": rng.uniform(-3.0, -1.0, (S, E)),
                             "nonmaxw_corr": nmc})
    df = _group("/dist_fn", {
        "fexp": np.asarray(fexp[:S], float), "bkdiff": np.asarray(bkdiff[:S], float),
        "driftknob": driftknob, "tpdriftknob": tpdriftknob, "mach": mach,
        "g_exb": g_exb, "omprimfac": 1.0, "phi_ext": 0.0, "t0": 100.0, "omega0": 0.0,
        "gamma0": 0.0, "source0": 1.0, "source_option_switch": source_option,
        "opt_source": False, "ufac": rng.uniform(-0.3, 0.3, (S, E)),
    })
    da = _group("/dist_fn_arrays", {
        "vpa": g7(vpa), "vpac": g7(vpac), "vperp2": g6(vperp2),
        "aj0": rng.uniform(0.2, 1.0, (S, E, L, Y, X, T)),
        "aj1": rng.uniform(0.1, 0.5, (S, E, L, Y, X, T)),
    })
    # GS2 has tunits = 1 (wstar_units is gated); non-unit values make the vpar chain
    # order observable in the bitwise tests.
    gt = _group("/gs2_time", {"wunits": aky / 2.0, "tunits": 1.0 + 0.1 * rng.uniform(size=Y)})
    rp = _group("/run_parameters", {        # non-unit factors make association order visible
        "fphi": 0.95, "fapar": 0.8 if has_apar else 0.0, "fbpar": 0.9 if has_bpar else -1.0,
        "has_phi": True, "has_apar": has_apar, "has_bpar": has_bpar})
    hy = _group("/hyper", {"D_res": rng.uniform(0.0, 0.02, (Y, X))})
    nl = _group("/nonlinear_terms", {"nonlin": nonlin, "split_nonlinear": False})
    st = _group("/state_t0", {"code_dt": code_dt})
    root = _group("/", {}, {g.name: g for g in (tg, kt, le, sp, df, da, gt, rp, hy, nl, st)})
    root.path = "/"
    dims = {"theta": T, "kx": X, "ky": Y, "lambda": L, "egrid": E, "species": S}
    ex = Export(path=Path("<synthetic>"), root=root, dims=dims, attrs={"gate_status": "ok"})
    # dt-dependent groups from the scalar reference
    ref = ref_coeffs(ex, code_dt)
    for name in ("vpar", "wdrift", "wdriftttp", "a", "b", "r", "ainv"):
        df.arrays[name] = ref[name]
    da.arrays["wstar"] = ref["wstar"]
    rungs = {}
    for k, dt in enumerate(ladder, start=1):
        rk = ref_coeffs(ex, dt)
        rungs[f"k{k}"] = _group(f"/dt_ladder/k{k}", {}, {
            "dist_fn": _group(f"/dt_ladder/k{k}/dist_fn",
                              {n: rk[n] for n in ("vpar", "wdrift", "wdriftttp", "a", "b", "r",
                                                   "ainv")}),
            "dist_fn_arrays": _group(f"/dt_ladder/k{k}/dist_fn_arrays", {"wstar": rk["wstar"]}),
        }, attrs={"code_dt": dt})
    if rungs:
        root.groups["dt_ladder"] = _group("/dt_ladder", {}, rungs)
    return ex


# ========================================================== scalar Fortran transcription
def cdiv_gcc(a: complex, b: complex) -> complex:
    """GCC's Smith division (ties take the second branch), 03 section 3.0."""
    ar, ai, br, bi = a.real, a.imag, b.real, b.imag
    if abs(br) < abs(bi):
        ratio = br / bi
        div = br * ratio + bi
        return complex((ar * ratio + ai) / div, (ai * ratio - ar) / div)
    ratio = bi / br
    div = bi * ratio + br
    return complex((ai * ratio + ar) / div, (ai - ar * ratio) / div)


def ref_coeffs(ex: Export, code_dt: float) -> dict:
    """Literal scalar transcription of init_vpar, init_wdrift (+ wdrift_func,
    wcoriolis_func), init_wstar and init_invert_rhs (src/dist_fn.fpp:1405-1929).

    Every operand is a Python float/complex (``tolist``): Python evaluates ``*``/``/``
    left to right like Fortran, never contracts to FMA, and (3.14) mixes real and
    complex operands componentwise, as GCC lowers them.  NumPy complex scalars are
    avoided (their multiply contracts to FMA), and the complex division uses
    :func:`cdiv_gcc` (CPython's algorithm differs).  Fortran's unary minus binds after
    ``*``; the value is the same (IEEE rounding is sign-symmetric).
    """
    tg, kt, le, sp, df = (ex[p] for p in ("/theta_grid", "/kt_grids", "/le_grids", "/species",
                                          "/dist_fn"))
    da, gt = ex["/dist_fn_arrays"], ex["/gs2_time"]
    S, E, L, Y, X, _, T = da["vpa"].shape
    vpa, vpac, vperp2 = (np.asarray(da[n]).tolist() for n in ("vpa", "vpac", "vperp2"))
    lst = lambda grp, n: np.asarray(grp[n]).tolist()          # noqa: E731
    delthet, gradpar = lst(tg, "delthet"), lst(tg, "gradpar")
    gbdrift, gbdrift0 = lst(tg, "gbdrift"), lst(tg, "gbdrift0")
    cvdrift, cvdrift0 = lst(tg, "cvdrift"), lst(tg, "cvdrift0")
    cdrift, cdrift0, shat = lst(tg, "cdrift"), lst(tg, "cdrift0"), float(tg["shat"])
    aky, akx, theta0 = lst(kt, "aky"), lst(kt, "akx"), lst(kt, "theta0")
    wunits, tunits = lst(gt, "wunits"), lst(gt, "tunits")
    zstm, stm, tz = lst(sp, "spec_zstm"), lst(sp, "spec_stm"), lst(sp, "spec_tz")
    dlnf0drho = lst(sp, "dlnf0drho")
    fexp, bkdiff = np.atleast_1d(df["fexp"]).tolist(), np.atleast_1d(df["bkdiff"]).tolist()
    driftknob, tpdriftknob, mach = float(df["driftknob"]), float(df["tpdriftknob"]), float(df["mach"])
    if tpdriftknob == -9.9e9:
        tpdriftknob = driftknob
    forbid, is_ttp, can_be_ttp = lst(le, "forbid"), lst(le, "is_ttp"), lst(le, "can_be_ttp")
    lower = lst(le, "is_lower_bounce_point")
    il_is_wfb, il_is_trapped = lst(le, "il_is_wfb"), lst(le, "il_is_trapped")
    trapped_grid = L > int(le["ng2"])
    zi = 1j

    vpar = np.zeros((S, E, L, Y, X, 2, T))
    wdrift = np.zeros_like(vpar)
    wdriftttp = np.zeros_like(vpar)
    a = np.zeros(vpar.shape, complex)
    b = np.zeros_like(a)
    r = np.zeros_like(a)
    ainv = np.zeros_like(a)
    wstar = np.zeros((S, E, Y))

    def wdrift_func(s, e, l, y, x, ig):                    # dist_fn.fpp:1500-1530
        vperp2_local = vperp2[s][e][l][y][x][ig]
        vpa2_local = vpa[s][e][l][y][x][0][ig] * vpa[s][e][l][y][x][0][ig]
        if aky[y] == 0.0:
            return akx[x]/shat \
                * (cvdrift0[ig]*vpa2_local
                   + gbdrift0[ig]*0.5*vperp2_local) \
                * code_dt/2.0
        return ((cvdrift[ig] + theta0[y][x]*cvdrift0[ig])
                * vpa2_local
                + (gbdrift[ig] + theta0[y][x]*gbdrift0[ig])
                * 0.5*vperp2_local) \
            * code_dt*wunits[y]

    def wcoriolis_func(s, e, l, y, x, ig):                 # dist_fn.fpp:1534-1557
        if aky[y] == 0.0:
            return mach*vpa[s][e][l][y][x][0][ig] \
                * cdrift0[ig] * code_dt * akx[x]/(2.*shat*stm[s])
        return mach*vpa[s][e][l][y][x][0][ig] \
            * (cdrift[ig] + theta0[y][x]*cdrift0[ig])*code_dt*wunits[y]/stm[s]

    for s in range(S):
        for e in range(E):
            for l in range(L):
                for y in range(Y):
                    for x in range(X):
                        # ---- init_vpar (dist_fn.fpp:1715-1725)
                        for isg in range(2):
                            for ig in range(T - 1):
                                vpar[s, e, l, y, x, isg, ig] = \
                                    zstm[s]*tunits[y]*code_dt \
                                    * 0.5/delthet[ig] \
                                    * (abs(gradpar[ig]) + abs(gradpar[ig + 1])) \
                                    * vpac[s][e][l][y][x][isg][ig]
                        # ---- init_wdrift (dist_fn.fpp:1423-1491)
                        knob = tpdriftknob if (il_is_wfb[l] or il_is_trapped[l]) else driftknob
                        w = [[0.0] * T, [0.0] * T]
                        for ig in range(T):
                            if forbid[l][ig]:
                                continue
                            wdv = wdrift_func(s, e, l, y, x, ig)
                            wcv = wcoriolis_func(s, e, l, y, x, ig)
                            w[0][ig] = wdv * knob + wcv
                            w[1][ig] = wdv * knob - wcv
                        if can_be_ttp[l]:
                            for ig in range(T):
                                if not is_ttp[l][ig]:
                                    continue
                                wdv = wdrift_func(s, e, l, y, x, ig)
                                wcv = wcoriolis_func(s, e, l, y, x, ig)
                                wdriftttp[s, e, l, y, x, 0, ig] = wdv * tpdriftknob + wcv
                                wdriftttp[s, e, l, y, x, 1, ig] = wdv * tpdriftknob - wcv
                        for ig in range(T - 1):          # in place, ascending
                            for isg in range(2):
                                w[isg][ig] = 0.5*(w[isg][ig] + w[isg][ig + 1])
                        wdrift[s, e, l, y, x] = w
    for s in range(S):                                      # init_wstar (:1746)
        for e in range(E):
            for y in range(Y):
                wstar[s, e, y] = - code_dt*wunits[y] * dlnf0drho[s][e]
    wd_l, wdttp_l, vpar_l = wdrift.tolist(), wdriftttp.tolist(), vpar.tolist()
    # ---- init_invert_rhs (dist_fn.fpp:1862-1925)
    for s in range(S):
        for e in range(E):
            for l in range(L):
                for y in range(Y):
                    for x in range(X):
                        for isg in range(2):
                            for ig in range(T - 1):
                                wd = wd_l[s][e][l][y][x][isg][ig]
                                wdttp = wdttp_l[s][e][l][y][x][isg][ig]
                                vp = vpar_l[s][e][l][y][x][0][ig]
                                bd = bkdiff[s]
                                ai = cdiv_gcc(1.0 + 0j, 1.0 + bd
                                              + (1.0-fexp[s])*tz[s]*(zi*wd*(1.0+bd) + 2.0*vp))
                                rr = (1.0 - bd
                                      + (1.0-fexp[s])*tz[s]*(zi*wd*(1.0-bd) - 2.0*vp)) \
                                    * ai
                                aa = 1.0 + bd \
                                    + fexp[s]*tz[s]*(-(zi*wd*(1.0+bd)) - 2.0*vp)
                                bb = 1.0 - bd \
                                    + fexp[s]*tz[s]*(-(zi*wd*(1.0-bd)) + 2.0*vp)
                                if trapped_grid:
                                    if forbid[l][ig] or forbid[l][ig + 1]:
                                        rr = 0.0
                                        ai = 0.0
                                    if isg == 0 and lower[l][ig + 1]:
                                        ai = 1.0 + ai
                                    if is_ttp[l][ig]:
                                        ai = cdiv_gcc(1.0 + 0j, 1.0 + zi*(1.0-fexp[s])*tz[s]*wdttp)
                                        aa = 1.0 - zi*fexp[s]*tz[s]*wdttp
                                        rr = 0.0
                                p = (s, e, l, y, x, isg, ig)
                                ainv[p], r[p], a[p], b[p] = ai, rr, aa, bb
    return {"vpar": vpar, "wdrift": wdrift, "wdriftttp": wdriftttp, "wstar": wstar,
            "a": a, "b": b, "r": r, "ainv": ainv}


# ==================================================================== synthetic tests
SYNTH_CONFIGS = {
    "trapped_zonal_box": dict(),
    "untrapped_finite_ky": dict(trapped=False, zonal=False, box=False, mach=0.0, seed=1),
    "maxwell_mach0_knobs": dict(mach=0.0, nonmaxwell=False, tpdriftknob=-9.9e9, seed=2,
                                bkdiff=(0.0, 0.0), fexp=(0.5, 0.5)),
    "one_species": dict(nspec=1, negrid=2, ntgrid=3, naky=2, ntheta0=3, seed=3),
}


@pytest.fixture(scope="module", params=sorted(SYNTH_CONFIGS))
def synth(request):
    cfg = dict(SYNTH_CONFIGS[request.param])
    seed = cfg.pop("seed", 0)
    return request.param, make_synthetic_export(seed, **cfg)


def test_compress_axes():
    x = np.broadcast_to(np.arange(6.0).reshape(2, 1, 3), (2, 4, 3))
    small = co.compress_axes(x, (1,))
    assert small.shape == (2, 1, 3) and np.array_equal(np.broadcast_to(small, x.shape), x)
    y = np.array(x)
    y[1, 2, 0] += 1e-16 + 1.0
    assert co.compress_axes(y, (1,)).shape == y.shape


def test_exact_div_defeats_reciprocal_rewrite():
    rng = np.random.default_rng(5)
    num = rng.standard_normal((7, 300))
    den = rng.uniform(1.0, 2.0, 7)
    got = jax.jit(lambda n, d: co.exact_div(n, d[:, None]))(num, den)
    assert bitwise_equal(num / den[:, None], got)


def test_coeffs_synthetic_strict_bitwise(synth):
    name, ex = synth
    inputs, params = co.CoeffInputs.from_export(ex), co.CoeffParams.from_export(ex)
    if name != "untrapped_finite_ky":   # redundant ky/kx axes must have been compressed
        assert inputs.vpa.shape[3:5] == (1, 1) and inputs.vperp2.shape[3:5] == (1, 1)
    ref = Coeffs.from_export(ex)
    with strict_mode():
        got = jit_coeffs(params)(ref.code_dt, to_jax(inputs))
    for field in COEFF_NAMES:
        assert bitwise_equal(getattr(ref, field), getattr(got, field)), \
            f"{name}: {field} differs: {numerics.max_abs_rel(getattr(ref, field), getattr(got, field)):.3e}"
    assert got.g_h is None


def test_coeffs_synthetic_default_mode(synth):
    """Non-strict (FMA allowed) jit and eager evaluation agree with the reference to ~1 ulp."""
    name, ex = synth
    inputs, params = co.CoeffInputs.from_export(ex), co.CoeffParams.from_export(ex)
    ref = Coeffs.from_export(ex)
    got_jit = jit_coeffs(params)(ref.code_dt, to_jax(inputs))
    got_eager = co.compute_coeffs(ref.code_dt, to_jax(inputs), params)
    for field in COEFF_NAMES:
        assert_close(getattr(ref, field), getattr(got_jit, field), 1e-14, f"{name} jit {field}")
        assert_close(getattr(ref, field), getattr(got_eager, field), 1e-14, f"{name} eager {field}")


def test_coeffs_synthetic_ladder(synth):
    name, ex = synth
    inputs, params = to_jax(co.CoeffInputs.from_export(ex)), co.CoeffParams.from_export(ex)
    fn = jit_coeffs(params)
    for k in ladder_rungs(ex):
        ref = Coeffs.from_export(ex, rung=k)
        got = fn(ref.code_dt, inputs)                      # traced dt: no retrace per rung
        for field in COEFF_NAMES:
            assert_close(getattr(ref, field), getattr(got, field), 1e-14, f"{name} k{k} {field}")


def test_coeffs_masks_semantics():
    """Spot checks of the structure the source and invert steps rely on (02/03)."""
    ex = make_synthetic_export(0)
    inputs, params = co.CoeffInputs.from_export(ex), co.CoeffParams.from_export(ex)
    c = co.compute_coeffs(0.05, to_jax(inputs), params)
    vpar, a, b, r, ainv = (np.asarray(getattr(c, n)) for n in ("vpar", "a", "b", "r", "ainv"))
    assert np.array_equal(vpar[..., 1, :], -vpar[..., 0, :])            # 02 section 3.1
    for x in (vpar, a, b, r, ainv):
        assert np.all(x[..., -1] == 0)                                   # ig = ntgrid
    le = ex["/le_grids"]
    fz = le["forbid"][:, :-1] | le["forbid"][:, 1:]
    fz7 = np.broadcast_to(fz[None, None, :, None, None, None, :], r[..., :-1].shape)
    ttp7 = np.broadcast_to(le["is_ttp"][None, None, :, None, None, None, :-1], fz7.shape)
    lb7 = np.zeros(fz7.shape, bool)
    lb7[..., 0, :] = le["is_lower_bounce_point"][None, None, :, None, None, 1:]
    assert np.all(r[..., :-1][fz7 & ~ttp7] == 0)
    assert np.all(ainv[..., :-1][fz7 & ~ttp7 & ~lb7] == 0)
    assert np.all(ainv[..., :-1][fz7 & ~ttp7 & lb7] == 1)              # 03 section 3.3 facts
    assert ttp7.any() and np.all(a[..., :-1][ttp7].real == 1.0)
    wdt = np.asarray(c.wdriftttp)
    assert np.all(wdt[~np.broadcast_to(le["is_ttp"][None, None, :, None, None, None, :], wdt.shape)] == 0)


def test_derived_mask_fallbacks():
    """If the exporter omitted the derived masks, their definitions are recomputed."""
    ex = make_synthetic_export(0, ladder=())
    le = ex["/le_grids"]
    want_lb = le["is_lower_bounce_point"].copy()
    want_wfb, want_tr = le["il_is_wfb"].copy(), le["il_is_trapped"].copy()
    assert want_lb.any()
    for name in ("is_lower_bounce_point", "il_is_wfb", "il_is_trapped"):
        del le.arrays[name]
    assert np.array_equal(co.lower_bounce_points(ex), want_lb)
    wfb, tr = co.pitch_angle_classes(ex)
    assert np.array_equal(wfb, want_wfb) and np.array_equal(tr, want_tr)


def test_q07_gate():
    ex = make_synthetic_export(0, ladder=())
    ex["/theta_grid"].arrays["shat"] = 0.0
    with pytest.raises(NotImplementedError):
        co.CoeffParams.from_export(ex)


# ===================================================================== refdata tests
CASES = cases_where(lambda ex: True)


def _compare(ex, ref: Coeffs, got: Coeffs, what: str, rtol=RTOL_KERNEL) -> dict:
    errs = {}
    for field in COEFF_NAMES:
        errs[field] = assert_close(getattr(ref, field), getattr(got, field), rtol,
                                   f"{what} {field}")
    return errs


@pytest.mark.parametrize("case", CASES)
def test_coeffs_initial_dt(case):
    ex = export_of(case)
    inputs, params = to_jax(co.CoeffInputs.from_export(ex)), co.CoeffParams.from_export(ex)
    ref = Coeffs.from_export(ex)
    got = jit_coeffs(params)(ref.code_dt, inputs)
    errs = _compare(ex, ref, got, f"{case} dt={ref.code_dt!r}")
    print(f"{case} initial dt: " + ", ".join(f"{k}={v:.2e}" for k, v in errs.items()))


@pytest.mark.parametrize("case", CASES)
def test_coeffs_dt_ladder(case):
    ex = export_of(case)
    rungs = ladder_rungs(ex)
    if not rungs:
        pytest.skip(f"{case}: no dt ladder in the export")
    inputs, params = to_jax(co.CoeffInputs.from_export(ex)), co.CoeffParams.from_export(ex)
    fn = jit_coeffs(params)
    for k in rungs:
        ref = Coeffs.from_export(ex, rung=k)
        got = fn(ref.code_dt, inputs)
        errs = _compare(ex, ref, got, f"{case} rung k{k} dt={ref.code_dt!r}")
        print(f"{case} k{k}: " + ", ".join(f"{n}={v:.2e}" for n, v in errs.items()))


@pytest.mark.parametrize("case", CASES)
def test_coeffs_strict_exactness(case):
    """Strict mode (no FMA) against the -ffp-contract=off reference build: bitwise."""
    ex = export_of(case)
    inputs, params = to_jax(co.CoeffInputs.from_export(ex)), co.CoeffParams.from_export(ex)
    refs = [Coeffs.from_export(ex)] + [Coeffs.from_export(ex, rung=k) for k in ladder_rungs(ex)]
    with strict_mode():
        fn = jit_coeffs(params)
        for ref in refs:
            got = fn(ref.code_dt, inputs)
            errs = _compare(ex, ref, got, f"{case} strict dt={ref.code_dt!r}",
                            rtol=strict_rtol(ex))
            print(f"{case} strict dt={ref.code_dt!r}: "
                  + ", ".join(f"{n}={v:.2e}" for n, v in errs.items()))
