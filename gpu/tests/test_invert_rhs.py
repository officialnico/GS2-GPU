"""Tests for the implicit parallel solve (spec 03; spec 07 sections 3.6, 3.8).

Part 1 (always runs): synthetic grids built with the le_grids formulas
(s-alpha-like ``bmag``, GS2-style lambda grid with passing, wfb, trapped and
totally trapped pitch angles, box twist-and-shift topology), compared with a
**literal scalar transcription of the Fortran** (one IEEE operation per Python
operator, so no FMA): ``init_homogeneous_g``, ``invert_rhs_1`` (with
``apply_parallel_boundary_conditions``, ``self_periodic``, ``set_h_zero``,
``enforce_parity``) and ``apply_linked_boundary_conditions`` with an emulation
of the ``links_p/links_h/wfb_p/wfb_h`` fills into ``g_adj`` (slot layout of
``setup_connected_bc_redistribute*``). In strict mode (``numerics.set_strict``,
no FMA in XLA) the JAX kernels must agree exactly; in the default (jit, FMA)
mode to round-off.

Part 2 (needs ``gpu/refdata``): the probe comparisons of spec 11 section 6 for
every case, every traced step and both passes, and ``init_homogeneous_g``
against the exported ``g_h`` at the initial dt and every dt-ladder rung (D1
tolerance). GS2 is not reproduced bit for bit: the ``-ffp-contract=off``
reference build still fuses the complex products that GCC vectorises into
ARMv8.3 ``fcmla`` instructions. ``test_g_h_rounding_model`` shows that this
alone explains every last-bit difference of ``g_h``.
"""

from __future__ import annotations

import itertools
from functools import lru_cache

import jax
import numpy as np
import pytest

from conftest import RTOL_KERNEL, all_cases, assert_close, export_of, have, trace_of

import gs2gpu  # noqa: F401  (x64)
from gs2gpu import invert_rhs as ir
from gs2gpu import linked as lk
from gs2gpu import numerics
from gs2gpu import parity as par
from gs2gpu.types import Coeffs, Dims, fields_from_group, to_jax

LINKED, ZERO_BC, SELF_PERIODIC, ALT_ZERO = (ir.BOUNDARY_LINKED, ir.BOUNDARY_ZERO,
                                            ir.BOUNDARY_SELF_PERIODIC,
                                            ir.BOUNDARY_ALTERNATE_ZERO)

# =============================================================================
# Part 1a: scalar complex arithmetic on (re, im) float pairs (no FMA possible)
# =============================================================================
ZERO = (0.0, 0.0)
ONE = (1.0, 0.0)


def c_(z):
    return (float(z.real), float(z.imag))


def z_(c):
    return complex(c[0], c[1])


def c_add(a, b):
    return (a[0] + b[0], a[1] + b[1])


def c_sub(a, b):
    return (a[0] - b[0], a[1] - b[1])


def c_neg(a):
    return (-a[0], -a[1])


def c_mul(a, b):
    return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])


def c_rmul(a, r):
    return (a[0] * r, a[1] * r)


def c_div(a, b):
    """GCC's Smith division (spec 03 section 3.0)."""
    ar, ai = a
    br, bi = b
    if abs(br) < abs(bi):
        r = br / bi
        d = br * r + bi
        return ((ar * r + ai) / d, (ai * r - ar) / d)
    r = bi / br
    d = bi * r + br
    return ((ai * r + ar) / d, (ai - ar * r) / d)


def r_sub(r, a):
    """Fortran ``r - z`` with real ``r``."""
    return (r - a[0], -a[1])


def np_cdiv(a, b):
    """Vectorised Smith division (NumPy, only used to build synthetic coefficients)."""
    ar, ai, br, bi = a.real, a.imag, b.real, b.imag
    with np.errstate(divide="ignore", invalid="ignore"):
        t = np.abs(br) < np.abs(bi)
        r1 = br / bi
        d1 = br * r1 + bi
        r2 = bi / br
        d2 = bi * r2 + br
        re = np.where(t, (ar * r1 + ai) / d1, (ai * r2 + ar) / d2)
        im = np.where(t, (ai * r1 - ar) / d1, (ai - ar * r2) / d2)
    return re + 1j * im


# =============================================================================
# Part 1b: synthetic grids (le_grids.f90 formulas, literal)
# =============================================================================
class Grid:
    """Synthetic run: grids, masks, topology and knobs (all NumPy, canonical layouts)."""


def _le_flags(al, bmag, ng2, trapped_wfb):
    """forbid, is_bounce_point, lower/upper bounce points, is_ttp: literal le_grids loops."""
    nlambda, nth = len(al), len(bmag)
    fuzz = 10 * np.finfo(np.float64).eps
    forbid = np.zeros((nlambda, nth), bool)
    bounce = np.zeros((nlambda, nth), bool)
    for il in range(1, nlambda + 1):
        for t in range(nth):
            forbid[il - 1, t] = 1.0 - al[il - 1] * bmag[t] < -fuzz
            if il >= ng2 + 2:
                bounce[il - 1, t] = abs(1.0 - al[il - 1] * bmag[t]) <= fuzz
    lower = np.zeros_like(bounce)
    upper = np.zeros_like(bounce)
    for il in range(1, nlambda + 1):
        for t in range(nth):
            if not bounce[il - 1, t] or (il == ng2 + 1 and not trapped_wfb):
                continue
            if t == 0:
                lower[il - 1, t], upper[il - 1, t] = True, False
            elif t == nth - 1:
                lower[il - 1, t], upper[il - 1, t] = False, True
            else:
                fm, fp = forbid[il - 1, t - 1], forbid[il - 1, t + 1]
                lower[il - 1, t] = fm or (fm and fp)
                upper[il - 1, t] = fp or (fm and fp)
    ittp = np.full(nth, nlambda + 1)
    if nlambda > ng2:
        for t in range(1, nth - 1):
            for il in range(ng2 + 1, nlambda + 1):
                if forbid[il - 1, t - 1] and forbid[il - 1, t + 1] and not forbid[il - 1, t]:
                    ittp[t] = il
                    break
    ttp = np.zeros((nlambda, nth), bool)
    if nlambda > ng2:
        for il in range(ng2 + 1, nlambda + 1):
            for t in range(nth):
                if il >= ittp[t] and not forbid[il - 1, t]:
                    ttp[il - 1, t] = True
    return forbid, bounce, lower, upper, ttp


def _itleft_itright(jshift0, aky):
    """compute_itleft_and_itright (dist_fn.fpp:2441-2513), literal; 1-based values."""
    naky = len(aky)

    def build(ntheta0):
        itleft = np.full((naky, ntheta0), -1, np.int32)
        itright = np.full((naky, ntheta0), -1, np.int32)
        for ik in range(1, naky + 1):
            for it in range(1, ntheta0 + 1):
                it0 = it - ntheta0 - 1 if it > (ntheta0 + 1) // 2 else it - 1
                if ik == 1:
                    if aky[0] != 0.0 and naky == 1:
                        itl, itr = it0 + jshift0, it0 - jshift0
                    else:
                        itl, itr = it0, it0
                else:
                    itl, itr = it0 + (ik - 1) * jshift0, it0 - (ik - 1) * jshift0
                for v, arr in ((itl, itleft), (itr, itright)):
                    if 0 <= v < (ntheta0 + 1) // 2:
                        arr[ik - 1, it - 1] = v + 1
                    elif (ntheta0 + 1) // 2 < v + ntheta0 + 1 <= ntheta0:
                        arr[ik - 1, it - 1] = v + ntheta0 + 1
                    else:
                        arr[ik - 1, it - 1] = -1
        return itleft, itright

    return build


def _count_links(itc):
    """count_links_one_side (dist_fn.fpp:2608-2640), literal."""
    naky, ntheta0 = itc.shape
    links = np.zeros((naky, ntheta0), np.int32)
    for ik in range(naky):
        for it in range(ntheta0):
            cur = it
            while True:
                nxt = itc[ik, cur] - 1
                if cur == nxt or nxt < 0:
                    break
                links[ik, it] += 1
                cur = nxt
    return links


def make_grid(seed=0, *, bos=LINKED, wfbbc="mixed", nonad_zero=True, zero_forbid=False,
              def_parity=False, even=True, hybrid=False, nperiod=1, ntheta=8, ng2=2,
              nspec=2, negrid=2, jshift0=2, wfb=1.0, fphi=1.0, fbpar=1.0):
    rng = np.random.default_rng(seed)
    g = Grid()
    g.bos = bos
    g.linked = bos == LINKED
    g.passing_wfb, g.trapped_wfb = wfbbc == "passing", wfbbc == "trapped"
    g.mixed_wfb = not (g.passing_wfb or g.trapped_wfb)
    g.nonad_zero, g.zero_forbid, g.def_parity, g.even = nonad_zero, zero_forbid, def_parity, even
    g.wfb, g.fphi, g.fbpar = wfb, fphi, fbpar
    g.ng2 = ng2
    g.ntgrid = ntheta // 2 + (nperiod - 1) * ntheta
    nth = 2 * g.ntgrid + 1
    t = np.arange(nth)
    k = (t - g.ntgrid) % ntheta
    kk = np.minimum(k, ntheta - k)                          # exactly symmetric and periodic
    eps = 0.18
    g.bmag = 1.0 / (1.0 + eps * np.cos(kk * (2 * np.pi / ntheta)))
    half = ntheta // 2
    bvals = 1.0 / (1.0 + eps * np.cos(np.arange(half, -1, -1) * (2 * np.pi / ntheta)))
    al_pass = np.sort(rng.uniform(0.05, 0.95, ng2)) / bvals[0]
    g.al = np.concatenate([al_pass, 1.0 / bvals])            # wfb = 1/bmax ... ttp = 1/bmin
    nlambda = len(g.al)
    g.forbid, g.bounce, g.lbp, g.ubp, g.ttp = _le_flags(g.al, g.bmag, ng2, g.trapped_wfb)
    assert g.ttp.any() and g.forbid.any()

    if g.linked:
        naky, ntheta0 = 4, 7
        g.aky = np.arange(naky) * 0.25
        g.box = True
        g.itleft, g.itright = _itleft_itright(jshift0, g.aky)(ntheta0)
        g.l_links, g.r_links = _count_links(g.itleft), _count_links(g.itright)
    else:
        naky, ntheta0 = 3, 2
        g.aky = np.array([0.0, 0.3, 0.5])                   # include a ky = 0 mode
        g.box = False
        g.itleft = g.itright = None
        g.l_links = g.r_links = np.zeros((naky, ntheta0), np.int32)
    g.kwork = np.zeros((naky, ntheta0), bool)
    if g.box:
        g.kwork[0, 0] = True
    g.ieqzip = np.zeros((naky, ntheta0), bool)
    S, E, L, Y, X, T = nspec, negrid, nlambda, naky, ntheta0, nth
    g.shape7 = (S, E, L, Y, X, 2, T)
    il = np.arange(1, L + 1)
    g.stype = np.array([1, 4 if hybrid else 2][:S])
    g.hybrid = ((g.stype == 4)[:, None, None] & (g.aky != 0.0)[None, :, None]
                & (il <= ng2 + 1)[None, None, :])
    g.zt = rng.uniform(0.5, 2.0, S)
    g.nonmaxw = rng.uniform(0.8, 1.2, (S, E))
    g.vperp2 = rng.uniform(0.1, 2.0, (S, E, L, Y, X, T))
    g.aj0 = rng.uniform(0.3, 1.0, (S, E, L, Y, X, T))
    g.aj1 = rng.uniform(0.1, 0.5, (S, E, L, Y, X, T))
    # mixed-wfb end factors: functions of (is, ik, itl, ie) / (is, ik, itr, ie), pre-gathered (G-04)
    ends_l = rng.uniform(0.3, 1.0, (4, S, Y, X, E))
    g.wfb_end = {}
    if g.linked:
        lm = np.zeros((Y, X), int)
        rm = np.zeros((Y, X), int)
        for iy, ix in itertools.product(range(Y), range(X)):
            cur = ix
            while g.itleft[iy, cur] >= 0 and g.itleft[iy, cur] - 1 != cur:
                cur = g.itleft[iy, cur] - 1
            lm[iy, ix] = cur
            cur = ix
            while g.itright[iy, cur] >= 0 and g.itright[iy, cur] - 1 != cur:
                cur = g.itright[iy, cur] - 1
            rm[iy, ix] = cur
        g.leftmost, g.rightmost = lm, rm
        yy = np.arange(Y)[:, None]
        g.end_vals = ends_l                                  # indexed [k, is, ik, it_end, ie]
        g.wfb_end = {
            "wfb_aj0l": ends_l[0][:, yy, lm, :], "wfb_aj1l": ends_l[1][:, yy, lm, :],
            "wfb_aj0r": ends_l[2][:, yy, rm, :], "wfb_aj1r": ends_l[3][:, yy, rm, :],
            "wfb_vperp2left": rng.uniform(0.1, 1.0, (S, E)),
            "wfb_vperp2right": rng.uniform(0.1, 1.0, (S, E)),
        }

    # ---- invert_rhs coefficients (03 section 3.3 component form) at some dt
    vp = rng.uniform(0.05, 0.6, (S, E, L, Y, X, T))
    wd = rng.normal(0.0, 0.05, (S, E, L, Y, X, 2, T))
    wdttp = rng.normal(0.0, 0.05, (S, E, L, Y, X, 2, T))
    fe = rng.uniform(0.3, 0.5, S)[:, None, None, None, None, None, None]
    bd = rng.uniform(0.0, 0.1, S)[:, None, None, None, None, None, None]
    tz = rng.uniform(0.5, 2.0, S)[:, None, None, None, None, None, None]
    v2 = (2.0 * vp)[..., None, :]
    c1 = (1.0 - fe) * tz
    opb, omb = 1.0 + bd, 1.0 - bd
    D = (opb + c1 * v2) + 1j * (c1 * (wd * opb))
    ainv = np_cdiv(np.ones_like(D), D)
    N = (omb - c1 * v2) + 1j * (c1 * (wd * omb))
    r = (N * ainv)
    fb = g.forbid[None, None, :, None, None, None, :]
    fb_next = np.zeros_like(g.forbid)
    fb_next[:, :-1] = g.forbid[:, 1:]
    fbn = fb_next[None, None, :, None, None, None, :]
    zero_mask = (fb | fbn)
    r = np.where(zero_mask, 0.0, r)
    ainv = np.where(zero_mask, 0.0, ainv)
    lb_next = np.zeros_like(g.lbp)
    lb_next[:, :-1] = g.lbp[:, 1:]
    sig1 = np.zeros((2,), bool)
    sig1[0] = True
    m = lb_next[None, None, :, None, None, None, :] & sig1[None, None, None, None, None, :, None]
    ainv = np.where(m, 1.0 + ainv, ainv)
    ttp = g.ttp[None, None, :, None, None, None, :]
    ainv = np.where(ttp, np_cdiv(np.ones_like(D), 1.0 + 1j * (c1 * wdttp)), ainv)
    r = np.where(ttp, 0.0, r)
    r[..., -1] = 0.0
    ainv[..., -1] = 0.0
    g.r = r.astype(np.complex128)
    g.ainv = ainv.astype(np.complex128)

    def cplx(*shape):
        return rng.normal(size=shape) + 1j * rng.normal(size=shape)

    g.source = cplx(*g.shape7)
    g.source[..., -1] = 0.0
    g.phinew = cplx(Y, X, T)
    g.bparnew = cplx(Y, X, T)
    g.gnew_in = cplx(*g.shape7)
    return g


def build_from_grid(g):
    S, E, L, Y, X, _, T = g.shape7
    return ir.build(
        ntgrid=g.ntgrid, ng2=g.ng2, aky=g.aky, kwork_filter=g.kwork, ieqzip=g.ieqzip,
        forbid=g.forbid, is_ttp=g.ttp, is_lower_bounce_point=g.lbp,
        is_upper_bounce_point=g.ubp, is_passing_hybrid_electron=g.hybrid,
        passing_wfb=g.passing_wfb, trapped_wfb=g.trapped_wfb, mixed_wfb=g.mixed_wfb,
        boundary_option_switch=g.bos, nonad_zero=g.nonad_zero, zero_forbid=g.zero_forbid,
        def_parity=g.def_parity, even=g.even, wfb=g.wfb, fphi=g.fphi, fbpar=g.fbpar,
        vperp2=g.vperp2, aj0=g.aj0, aj1=g.aj1, zt=g.zt, nonmaxw_corr=g.nonmaxw,
        l_links=g.l_links, r_links=g.r_links, itleft=g.itleft, itright=g.itright,
        wfb_factors=g.wfb_end if g.linked else None, box=g.box)


# =============================================================================
# Part 1c: literal scalar transcription of the Fortran
# =============================================================================
def _lines(g):
    S, E, L, Y, X, _, T = g.shape7
    return itertools.product(range(S), range(E), range(L), range(Y), range(X))


def ref_init_homogeneous_g(g, r):
    """init_homogeneous_g (dist_fn.fpp:1937-2125), one line at a time."""
    S, E, L, Y, X, _, T = g.shape7
    out = np.zeros(g.shape7, np.complex128)
    uph_b = g.bos in (SELF_PERIODIC, LINKED)
    for s, e, l0, y, x in _lines(g):
        il = l0 + 1
        uph = uph_b or g.aky[y] == 0.0
        ilmin = 1 if uph else g.ng2 + 1
        if il < ilmin:
            continue
        is_wfb, is_passing, is_trapped = il == g.ng2 + 1, il <= g.ng2, il > g.ng2 + 1
        gh = [[ZERO, ZERO] for _ in range(T)]
        g2 = [[ZERO, ZERO] for _ in range(T)]
        if uph and is_passing:
            gh[0][0] = ONE
            gh[T - 1][1] = ONE
        if is_wfb and g.mixed_wfb:
            gh[0][0] = (g.wfb, 0.0)
            gh[T - 1][1] = (g.wfb, 0.0)
        elif is_wfb and g.trapped_wfb:
            gh[T - 1][1] = ONE
        elif is_wfb and g.passing_wfb:
            gh[0][0] = ONE
            gh[T - 1][1] = ONE
        if is_trapped:
            for t in range(T):
                if g.ttp[l0, t]:
                    continue
                if g.ubp[l0, t]:
                    g2[t][1] = ONE
        rr = r[s, e, l0, y, x]
        for t in range(T - 2, -1, -1):
            gh[t][1] = c_add(c_neg(c_mul(gh[t + 1][1], c_(rr[1, t]))), g2[t][1])
        if is_trapped:
            for t in range(T - 1):
                if g.ttp[l0, t]:
                    continue
                if g.lbp[l0, t + 1]:
                    g2[t + 1][0] = gh[t + 1][1]
        if is_wfb and g.trapped_wfb:
            gh[0][0] = gh[0][1]
        for t in range(T - 1):
            if g.ttp[l0, t]:
                continue
            gh[t + 1][0] = c_add(c_neg(c_mul(gh[t][0], c_(rr[0, t]))), g2[t + 1][0])
        for t in range(T):
            out[s, e, l0, y, x, 0, t] = z_(gh[t][0])
            out[s, e, l0, y, x, 1, t] = z_(gh[t][1])
    return out


def _neighbor(g, l0, y, x):
    il = l0 + 1
    if il > g.ng2 + 1 or (il == g.ng2 + 1 and g.trapped_wfb):
        return False, False, False
    left = g.itleft[y, x] >= 0
    right = g.itright[y, x] >= 0
    return left or right, left, right


def _save_h(g, l0, y, x):
    """set_save_h (dist_fn.fpp:2570-2604) on top of compute_connections (2518-2565)."""
    nb, left, right = _neighbor(g, l0, y, x)
    s1 = left and g.aky[y] != 0.0
    s2 = right and g.aky[y] != 0.0
    if l0 + 1 == g.ng2 + 1 and nb and (not g.passing_wfb) and g.aky[y] != 0.0:
        s1 = s2 = True
    return s1, s2


def _get_adjust(g, s, e, l0, y, x, t, facphi, facbpar):
    """get_adjust (dist_fn_arrays.fpp:178-188): 2*vperp2*aj1*bpar*facbpar + zt*phi*aj0*facphi*nonmaxw."""
    a = c_rmul(c_rmul(c_(g.bparnew[y, x, t]), (2.0 * g.vperp2[s, e, l0, y, x, t])
                      * g.aj1[s, e, l0, y, x, t]), facbpar)
    b = c_rmul(c_rmul(c_rmul(c_rmul(c_(g.phinew[y, x, t]), g.zt[s]), g.aj0[s, e, l0, y, x, t]),
                      facphi), g.nonmaxw[s, e])
    return c_add(a, b)


def ref_invert_rhs_1(g, g_h):
    """invert_rhs_1 (dist_fn.fpp:5440-5707), one line at a time, then enforce_parity."""
    S, E, L, Y, X, _, T = g.shape7
    nt = g.ntgrid
    gnew = g.gnew_in.copy()
    for s, e, l0, y, x in _lines(g):
        if g.kwork[y, x] or g.ieqzip[y, x]:
            continue
        il = l0 + 1
        if g.hybrid[s, y, l0]:
            facphi, facbpar = -g.fphi, -g.fbpar
            phi_factor = g.zt[s] * g.nonmaxw[s, e] * facphi if facphi != 0 else 0.0
            bpar_factor = 2.0 * facbpar if facbpar != 0 else 0.0
            for t in range(T):
                adj = c_add(c_rmul(c_(g.bparnew[y, x, t]), (bpar_factor * g.vperp2[s, e, l0, y, x, t])
                                   * g.aj1[s, e, l0, y, x, t]),
                            c_rmul(c_rmul(c_(g.phinew[y, x, t]), phi_factor),
                                   g.aj0[s, e, l0, y, x, t]))
                gnew[s, e, l0, y, x, 0, t] = z_(adj)
                gnew[s, e, l0, y, x, 1, t] = z_(adj)
            continue
        speriod = False
        if g.bos == SELF_PERIODIC:
            uph = True
        elif g.bos == LINKED:
            uph = True
            speriod = g.aky[y] == 0.0
        else:
            uph = False
        uph = uph or g.aky[y] == 0.0
        is_wfb, is_passing, is_trapped = il == g.ng2 + 1, il <= g.ng2, il > g.ng2 + 1
        is_pl = is_passing or (is_wfb and g.passing_wfb)
        if g.bos == LINKED:
            nb = _neighbor(g, l0, y, x)[0]
            is_sp = ((speriod and is_passing) or (speriod and is_wfb and not g.trapped_wfb)
                     or (is_wfb and not nb and g.mixed_wfb))
        else:
            is_sp = (uph and is_pl) or (is_wfb and g.mixed_wfb)
        G = [[ZERO, ZERO] for _ in range(T)]
        src = [[c_(g.source[s, e, l0, y, x, 0, t]), c_(g.source[s, e, l0, y, x, 1, t])]
               for t in range(T)]
        A = g.ainv[s, e, l0, y, x]
        R = g.r[s, e, l0, y, x]
        H = g_h[s, e, l0, y, x]
        if is_pl and not is_sp:
            is_left = g.l_links[y, x] == 0
            is_right = g.r_links[y, x] == 0
            if (is_left or is_right) and g.nonad_zero:
                if is_left:
                    G[0][0] = c_neg(_get_adjust(g, s, e, l0, y, x, 0, g.fphi, g.fbpar))
                if is_right:
                    G[T - 1][1] = c_neg(_get_adjust(g, s, e, l0, y, x, T - 1, g.fphi, g.fbpar))
        for t in range(T - 2, -1, -1):
            G[t][1] = c_add(c_neg(c_mul(G[t + 1][1], c_(R[1, t]))), c_mul(c_(A[1, t]), src[t][1]))
        if is_trapped:
            for t in range(T - 1):
                if g.ttp[l0, t]:
                    continue
                if g.lbp[l0, t + 1]:
                    src[t][0] = G[t + 1][1]
        if is_wfb and g.trapped_wfb:
            G[0][0] = G[0][1]
        for t in range(T - 1):
            if not g.ttp[l0, t]:
                G[t + 1][0] = c_mul(c_(A[0, t]), src[t][0])
                G[t + 1][0] = c_sub(G[t + 1][0], c_mul(G[t][0], c_(R[0, t])))
            else:
                if not g.forbid[l0, t]:
                    G[t][0] = G[t][1]
        if is_sp:
            if is_wfb and g.nonad_zero and not g.passing_wfb:
                def adj_sp(t):
                    a = c_rmul(c_rmul(c_(g.bparnew[y, x, t]), (2.0 * g.vperp2[s, e, l0, y, x, t])
                                      * g.aj1[s, e, l0, y, x, t]), g.fbpar)
                    b = c_rmul(c_rmul(c_rmul(c_rmul(c_(g.phinew[y, x, t]), g.zt[s]),
                                             g.aj0[s, e, l0, y, x, t]), g.nonmaxw[s, e]), g.fphi)
                    return c_add(a, b)
                dadj = c_sub(adj_sp(0), adj_sp(T - 1))
            else:
                dadj = ZERO
            if H[0, T - 1] != 1.0:
                b1 = c_div(c_sub(c_sub(G[T - 1][0], G[0][0]), dadj), r_sub(1.0, c_(H[0, T - 1])))
                for t in range(T):
                    G[t][0] = c_add(G[t][0], c_mul(b1, c_(H[0, t])))
            if H[1, 0] != 1.0:
                b1 = c_div(c_add(c_sub(G[0][1], G[T - 1][1]), dadj), r_sub(1.0, c_(H[1, 0])))
                for t in range(T):
                    G[t][1] = c_add(G[t][1], c_mul(b1, c_(H[1, t])))
        if is_trapped:
            beta1 = ZERO
            for t in range(T - 2, -1, -1):
                if g.forbid[l0, t] or g.ttp[l0, t]:
                    beta1 = ZERO
                    continue
                elif g.forbid[l0, t + 1]:
                    beta1 = c_div(c_sub(G[t][0], G[t][1]), r_sub(1.0, c_(H[0, t])))
                G[t][0] = c_add(G[t][0], c_mul(beta1, c_(H[0, t])))
                G[t][1] = c_add(G[t][1], c_mul(beta1, c_(H[1, t])))
        if is_wfb and g.trapped_wfb:
            beta1 = c_div(c_sub(G[T - 1][0], G[T - 1][1]), r_sub(1.0, c_(H[0, T - 1])))
            for t in range(T):
                G[t][0] = c_add(G[t][0], c_mul(beta1, c_(H[0, t])))
                G[t][1] = c_add(G[t][1], c_mul(beta1, c_(H[1, t])))
        if g.zero_forbid:
            for t in range(T):
                if g.forbid[l0, t]:
                    G[t][0] = ZERO
                    G[t][1] = ZERO
        for t in range(T):
            gnew[s, e, l0, y, x, 0, t] = z_(G[t][0])
            gnew[s, e, l0, y, x, 1, t] = z_(G[t][1])
    return ref_enforce_parity(g, gnew)


def _parity_conn(g, y, x):
    """get_parity_conn (dist_fn.fpp:4742-4793), literal walk."""
    nl, nr = g.l_links[y, x], g.r_links[y, x]
    tmp = x
    if nl == nr:
        return x
    if nl > nr:
        for _ in range(nl):
            tmp = g.itleft[y, tmp] - 1
            if g.r_links[y, tmp] == nl:
                break
    else:
        for _ in range(nr):
            tmp = g.itright[y, tmp] - 1
            if g.l_links[y, tmp] == nr:
                break
    return tmp


def ref_enforce_parity(g, gnew):
    """enforce_parity (dist_fn.fpp:4796-4853)."""
    if not g.def_parity:
        return gnew
    mult = 1 if g.even else -1
    S, E, L, Y, X, _, T = g.shape7
    nt = g.ntgrid
    out = gnew.copy()
    for s, e, l0, y, x in _lines(g):
        if g.linked:
            conn = _parity_conn(g, y, x)
            for ig in range(-nt, nt + 1):
                out[s, e, l0, y, x, 0, -ig + nt] = gnew[s, e, l0, y, conn, 1, ig + nt]
        else:
            for j in range(1, nt + 1):
                out[s, e, l0, y, x, 0, nt - j] = z_(c_rmul(c_(gnew[s, e, l0, y, x, 1, nt + j]),
                                                           float(mult)))
                out[s, e, l0, y, x, 0, nt + j] = z_(c_rmul(c_(gnew[s, e, l0, y, x, 1, nt - j]),
                                                           float(mult)))
    if g.linked and mult != 1:
        for s, e, l0, y, x in _lines(g):
            for t in range(T):
                out[s, e, l0, y, x, 0, t] = z_(c_rmul(c_(out[s, e, l0, y, x, 0, t]), float(mult)))
    return out


def ref_linked(g, gnew, g_h):
    """apply_linked_boundary_conditions (dist_fn.fpp:5919-6149) with the g_adj fills
    emulated from setup_connected_bc_redistribute(_mixed_wfb) (2759-3224)."""
    S, E, L, Y, X, _, T = g.shape7
    ncell = g.l_links + g.r_links + 1
    if (ncell == 1).all():
        return gnew                                           # no_connections
    nslot = 2 * int(ncell.max()) + 2
    g_adj = {}

    def slot(s, e, l0, y, x):
        key = (s, e, l0, y, x)
        if key not in g_adj:
            g_adj[key] = [[ZERO, ZERO] for _ in range(nslot + 1)]   # 1-based slots
        return g_adj[key]

    # init_homogeneous_g: fill(links_h, g_h, g_adj); fill(wfb_h, ...);  invert_rhs: links_p, wfb_p
    for s, e, l0, y, x in _lines(g):
        il = l0 + 1
        if il <= g.ng2 or (il == g.ng2 + 1 and g.passing_wfb):
            if ncell[y, x] == 1:
                continue
            cur = x
            for j in range(1, g.r_links[y, x] + 1):
                nxt = g.itright[y, cur] - 1
                slot(s, e, l0, y, nxt)[j][0] = c_(gnew[s, e, l0, y, x, 0, T - 1])
                if g.l_links[y, x] > 0:
                    slot(s, e, l0, y, nxt)[2 * g.l_links[y, x] + j][0] = c_(g_h[s, e, l0, y, x, 0, T - 1])
                cur = nxt
            cur = x
            for j in range(1, g.l_links[y, x] + 1):
                nxt = g.itleft[y, cur] - 1
                slot(s, e, l0, y, nxt)[j][1] = c_(gnew[s, e, l0, y, x, 1, 0])
                if g.r_links[y, x] > 0:
                    slot(s, e, l0, y, nxt)[2 * g.r_links[y, x] + j][1] = c_(g_h[s, e, l0, y, x, 1, 0])
                cur = nxt
        if g.mixed_wfb and il == g.ng2 + 1 and ncell[y, x] > 1:
            n = ncell[y, x]
            star = g.leftmost[y, x]
            for _ in range(n):
                slot(s, e, l0, y, star)[g.r_links[y, x] + 1][0] = c_(gnew[s, e, l0, y, x, 0, T - 1])
                slot(s, e, l0, y, star)[2 * n - g.r_links[y, x]][0] = c_(g_h[s, e, l0, y, x, 0, T - 1])
                star = g.itright[y, star] - 1
            star = g.rightmost[y, x]
            for _ in range(n):
                slot(s, e, l0, y, star)[g.l_links[y, x] + 1][1] = c_(gnew[s, e, l0, y, x, 1, 0])
                slot(s, e, l0, y, star)[2 * n - g.l_links[y, x]][1] = c_(g_h[s, e, l0, y, x, 1, 0])
                star = g.itleft[y, star] - 1

    out = gnew.copy()
    for s, e, l0, y, x in _lines(g):
        s1, s2 = _save_h(g, l0, y, x)
        if not (s1 or s2) or g.kwork[y, x]:
            continue
        n = ncell[y, x]
        if n == 1 or g.hybrid[s, y, l0]:
            continue
        ga = slot(s, e, l0, y, x)
        H = g_h[s, e, l0, y, x]
        if l0 + 1 == g.ng2 + 1 and g.mixed_wfb:
            if g.nonad_zero:
                itl, itr = g.leftmost[y, x], g.rightmost[y, x]
                vpl, vpr = g.wfb_end["wfb_vperp2left"][s, e], g.wfb_end["wfb_vperp2right"][s, e]
                aj0l, aj1l = g.end_vals[0][s, y, itl, e], g.end_vals[1][s, y, itl, e]
                aj0r, aj1r = g.end_vals[2][s, y, itr, e], g.end_vals[3][s, y, itr, e]
                adjl = c_add(c_rmul(c_rmul(c_(g.bparnew[y, itl, 0]), (2.0 * vpl) * aj1l), g.fbpar),
                             c_rmul(c_rmul(c_rmul(c_rmul(c_(g.phinew[y, itl, 0]), g.zt[s]), aj0l),
                                           g.nonmaxw[s, e]), g.fphi))
                adjr = c_add(c_rmul(c_rmul(c_(g.bparnew[y, itr, T - 1]), (2.0 * vpr) * aj1r), g.fbpar),
                             c_rmul(c_rmul(c_rmul(c_rmul(c_(g.phinew[y, itr, T - 1]), g.zt[s]), aj0r),
                                           g.nonmaxw[s, e]), g.fphi))
                dadj = c_sub(adjl, adjr)
            else:
                dadj = ZERO
            for sg, save, links, sign in ((0, s1, g.l_links[y, x], -1), (1, s2, g.r_links[y, x], 1)):
                if not save:
                    continue
                prod = ONE
                for j in range(n + 1, 2 * n + 1):
                    prod = c_mul(prod, ga[j][sg])
                facd = c_div(ONE, r_sub(1.0, prod))
                b0 = ZERO
                fac = ONE
                for i in range(n - 1, 0, -1):
                    fac = c_mul(fac, ga[n + i + 1][sg])
                    b0 = c_add(b0, c_mul(fac, ga[n + 1 - i][sg]))
                b0 = c_add(b0, ga[1][sg])
                b0 = c_sub(b0, dadj) if sign < 0 else c_add(b0, dadj)
                b0 = c_mul(b0, facd)
                for i in range(1, links + 1):
                    b0 = c_add(c_mul(b0, ga[n + i][sg]), ga[n + 1 - i][sg])
                for t in range(T):
                    out[s, e, l0, y, x, sg, t] = z_(c_add(c_(out[s, e, l0, y, x, sg, t]),
                                                          c_mul(b0, c_(H[sg, t]))))
        else:
            for sg, save, links in ((0, s1, g.l_links[y, x]), (1, s2, g.r_links[y, x])):
                if not (save and links > 0):
                    continue
                nl = 2 * links - 1                                   # n_links
                b0 = ga[1][sg]
                fac = ONE
                for i in range(2, links + 1):
                    fac = c_mul(fac, ga[nl + 2 - i][sg])
                    b0 = c_add(b0, c_mul(ga[i][sg], fac))
                for t in range(T):
                    out[s, e, l0, y, x, sg, t] = z_(c_add(c_(out[s, e, l0, y, x, sg, t]),
                                                          c_mul(b0, c_(H[sg, t]))))
    return out


# =============================================================================
# Part 1d: synthetic tests
# =============================================================================
SYNTH = {
    "linked_mixed": dict(),
    "linked_mixed_np2": dict(nperiod=2, ntheta=8),
    "linked_mixed_nonadF": dict(nonad_zero=False),
    "linked_passing": dict(wfbbc="passing"),
    "linked_trapped": dict(wfbbc="trapped"),
    "linked_zero_forbid": dict(zero_forbid=True),
    "linked_parity_even": dict(def_parity=True, even=True),
    "linked_parity_odd_negshear": dict(def_parity=True, even=False, jshift0=-2),
    "linked_hybrid": dict(hybrid=True, wfbbc="passing", fbpar=-1.0),
    "linked_jshift1": dict(jshift0=1, wfb=0.9),
    "zero_mixed_np2": dict(bos=ZERO_BC, nperiod=2),
    "altzero_passing": dict(bos=ALT_ZERO, wfbbc="passing"),
    "zero_trapped_parity_odd": dict(bos=ZERO_BC, wfbbc="trapped", def_parity=True, even=False),
    "zero_hybrid": dict(bos=ZERO_BC, hybrid=True),
    "selfperiodic_np2": dict(bos=SELF_PERIODIC, nperiod=2, nonad_zero=True),
    "selfperiodic_nonadF_trapped": dict(bos=SELF_PERIODIC, wfbbc="trapped", nonad_zero=False),
}


@lru_cache(maxsize=None)
def _synth(name):
    g = make_grid(seed=sum(map(ord, name)), **SYNTH[name])
    data, params = build_from_grid(g)
    gh_ref = ref_init_homogeneous_g(g, g.r)
    out1_ref = ref_invert_rhs_1(g, gh_ref)
    out_ref = ref_linked(g, out1_ref, gh_ref) if g.linked else out1_ref
    return g, data, params, gh_ref, out1_ref, out_ref


def _coeffs(g, g_h=None):
    zeros = np.zeros(g.shape7)
    return Coeffs(code_dt=0.05, vpar=zeros, wdrift=zeros, wdriftttp=zeros,
                  wstar=np.zeros((g.shape7[0], g.shape7[1], g.shape7[3])), a=g.ainv, b=g.ainv,
                  r=g.r, ainv=g.ainv, g_h=g_h)


def _exact(ref, got, what):
    got = np.asarray(got)
    bad = ref != got
    assert not bad.any(), (f"{what}: {bad.sum()} of {bad.size} values differ; max |diff| "
                           f"{np.max(np.abs(got - ref)):.3e} (max |ref| {np.max(np.abs(ref)):.3e})")


@pytest.fixture
def strict():
    numerics.set_strict(True)
    yield
    numerics.set_strict(False)


@pytest.mark.parametrize("name", sorted(SYNTH))
def test_synthetic_strict_exact(name, strict):
    """Strict mode (no FMA, eager): JAX == literal Fortran transcription, bit for bit."""
    g, data, params, gh_ref, out1_ref, out_ref = _synth(name)
    d = to_jax(data)
    gh = ir.init_homogeneous_g(_coeffs(g), d, params)
    _exact(gh_ref, gh, f"{name} init_homogeneous_g")
    co = _coeffs(g, gh_ref)
    out1 = ir.invert_rhs_1(g.source, g.phinew, g.bparnew, co, d, params, gnew_in=g.gnew_in)
    _exact(out1_ref, out1, f"{name} invert_rhs_1")
    out = lk.apply_linked_boundary_conditions(out1_ref, g.phinew, g.bparnew, co, d, params)
    _exact(out_ref, out, f"{name} linked")
    both = ir.invert_rhs(g.source, g.phinew, g.bparnew, co, d, params, gnew_in=g.gnew_in)
    _exact(out_ref, both, f"{name} invert_rhs")


@pytest.mark.parametrize("name", sorted(SYNTH))
def test_synthetic_jit(name):
    """Default mode under jit (XLA may contract FMAs): round-off agreement."""
    g, data, params, gh_ref, out1_ref, out_ref = _synth(name)
    d = to_jax(data)
    f_gh = jax.jit(ir.init_homogeneous_g, static_argnames=("params",))
    f_inv = jax.jit(ir.invert_rhs_1, static_argnames=("params",))
    f_lnk = jax.jit(lk.apply_linked_boundary_conditions, static_argnames=("params",))
    gh = f_gh(_coeffs(g), d, params=params)
    assert_close(gh_ref, gh, 1e-13, f"{name} g_h")
    co = to_jax(_coeffs(g, gh_ref))
    out1 = f_inv(g.source, g.phinew, g.bparnew, co, d, params=params, gnew_in=g.gnew_in)
    assert_close(out1_ref, out1, 1e-12, f"{name} invert_rhs_1")
    out = f_lnk(out1_ref, g.phinew, g.bparnew, co, d, params=params)
    assert_close(out_ref, out, 1e-12, f"{name} linked")


def test_synthetic_vmap_over_impulses():
    """The kernels vmap over a leading batch of (source, fields, gnew_in), as a response
    build does with its impulses (spec 04 section 5.4)."""
    g, data, params, gh_ref, out1_ref, out_ref = _synth("linked_mixed")
    d = to_jax(data)
    co = to_jax(_coeffs(g, gh_ref))
    rng = np.random.default_rng(5)
    nb = 3
    src = np.stack([g.source] + [g.source * rng.normal() for _ in range(nb - 1)])
    phi = np.stack([g.phinew] + [g.phinew * rng.normal() for _ in range(nb - 1)])
    bpar = np.stack([g.bparnew] * nb)
    gin = np.stack([g.gnew_in] * nb)

    def one(s, p, b, gi):
        return ir.invert_rhs(s, p, b, co, d, params, gnew_in=gi)

    got = jax.jit(jax.vmap(one))(src, phi, bpar, gin)
    assert_close(out_ref, got[0], 1e-12, "vmap element 0")
    for i in range(1, nb):
        assert_close(np.asarray(one(src[i], phi[i], bpar[i], gin[i])), got[i], 1e-12,
                     f"vmap element {i}")


def test_synthetic_parity_alone(strict):
    """enforce_parity on its own (the timeadv call site, dist_fn.fpp:4091)."""
    for name in ("linked_parity_even", "linked_parity_odd_negshear", "zero_trapped_parity_odd"):
        g, data, params, *_ = _synth(name)
        got = par.enforce_parity(g.gnew_in, to_jax(data), params)
        _exact(ref_enforce_parity(g, g.gnew_in), got, f"{name} parity")


def test_synthetic_g_h_none_and_passthrough_checks():
    """g_h=None is computed on the fly; gnew_in is required with kwork_filter cells."""
    g, data, params, gh_ref, out1_ref, out_ref = _synth("linked_mixed")
    d = to_jax(data)
    out = ir.invert_rhs(g.source, g.phinew, g.bparnew, _coeffs(g), d, params, gnew_in=g.gnew_in)
    assert_close(out_ref, out, 1e-12, "g_h computed on the fly")
    with pytest.raises(ValueError, match="gnew_in"):
        ir.invert_rhs_1(g.source, g.phinew, g.bparnew, _coeffs(g, gh_ref), d, params)
    # kwork_filter override: all cells solved (initial response-matrix build)
    g2 = make_grid(seed=sum(map(ord, "linked_mixed")), **SYNTH["linked_mixed"])
    g2.kwork[:] = False
    ref = ref_invert_rhs_1(g2, gh_ref)
    got = ir.invert_rhs_1(g.source, g.phinew, g.bparnew, _coeffs(g, gh_ref), d, params,
                          gnew_in=g.gnew_in, kwork_filter=np.zeros_like(g.kwork))
    assert_close(ref, got, 1e-12, "kwork_filter override")


def test_synthetic_topology():
    """Chains, parity partners and bounce-point helpers against literal walks."""
    for name in ("linked_mixed", "linked_parity_odd_negshear", "linked_jshift1"):
        g, data, params, *_ = _synth(name)
        ld = data.linked
        for c in range(ld.chain_n.shape[0]):
            n = ld.chain_n[c]
            cells = ld.chain_x[c, :n]
            y = ld.chain_y[c]
            assert g.leftmost[y, cells[0]] == cells[0] and g.rightmost[y, cells[-1]] == cells[-1]
            assert list(ld.chain_xr[c, :n]) == list(cells[::-1])
            for k, x in enumerate(cells):
                assert g.l_links[y, x] == k and g.r_links[y, x] == n - 1 - k
        pt = par.parity_partner(g.l_links, g.r_links, g.itleft, g.itright, True)
        for y, x in itertools.product(range(g.aky.size), range(g.kwork.shape[1])):
            assert pt[y, x] == _parity_conn(g, y, x)
        bp = g.bounce
        assert (ir.lower_bounce_points(bp, g.forbid, g.trapped_wfb, g.ng2) == g.lbp).all()
        assert (ir.upper_bounce_points(bp, g.forbid, g.trapped_wfb, g.ng2) == g.ubp).all()


def test_gates():
    g = make_grid(seed=1)
    for knob in ("start_from_previous_solution", "boundary_off_grid", "exponential_boundary"):
        with pytest.raises(NotImplementedError, match=knob):
            ir.build(**_build_kwargs(g), **{knob: True})
    with pytest.raises(NotImplementedError, match="linked"):
        ir.build(**(_build_kwargs(g) | {"box": False}))


def _fake_export(g, *, with_checks=True):
    """An in-memory ``io.export.Export`` with the spec 09 names (canonical NumPy layouts)."""
    from pathlib import Path

    from gs2gpu.io.export import Export, Group

    S, E, L, Y, X, _, T = g.shape7

    def grp(path, arrays):
        return Group(name=path.strip("/"), path=path, arrays=dict(arrays))

    dfn = {"boundary_option_switch": g.bos, "nonad_zero": g.nonad_zero,
           "zero_forbid": g.zero_forbid, "def_parity": g.def_parity, "even": g.even,
           "wfb": g.wfb, "start_from_previous_solution": False, "boundary_off_grid": False,
           "exponential_boundary": False, "l_links": g.l_links, "r_links": g.r_links,
           "r": g.r, "ainv": g.ainv, "a": g.ainv, "b": g.ainv, "vpar": np.zeros(g.shape7),
           "wdrift": np.zeros(g.shape7), "wdriftttp": np.zeros(g.shape7)}
    if g.linked:
        dfn |= {"itleft": g.itleft, "itright": g.itright, "no_connections": False}
        dfn |= {k: np.asarray(v) for k, v in g.wfb_end.items()}
    le = {"ng2": g.ng2, "passing_wfb": g.passing_wfb, "trapped_wfb": g.trapped_wfb,
          "mixed_wfb": g.mixed_wfb, "forbid": g.forbid, "is_ttp": g.ttp,
          "is_bounce_point": g.bounce}
    if with_checks:
        le |= {"is_lower_bounce_point": g.lbp, "is_upper_bounce_point": g.ubp,
               "is_passing_hybrid_electron": g.hybrid}
    root = Group(name="/", path="/", groups={
        "dist_fn": grp("/dist_fn", dfn),
        "le_grids": grp("/le_grids", le),
        "kt_grids": grp("/kt_grids", {"aky": g.aky, "kwork_filter": g.kwork, "box": g.box}),
        "run_parameters": grp("/run_parameters", {"fphi": g.fphi, "fbpar": g.fbpar,
                                                  "ieqzip": g.ieqzip}),
        "species": grp("/species", {"spec_zt": g.zt, "nonmaxw_corr": g.nonmaxw,
                                    "spec_type": g.stype}),
        "dist_fn_arrays": grp("/dist_fn_arrays", {"vperp2": g.vperp2, "aj0": g.aj0,
                                                  "aj1": g.aj1,
                                                  "wstar": np.zeros((S, E, Y))}),
        "state_t0": grp("/state_t0", {"code_dt": 0.05}),
    })
    dims = {"theta": T, "kx": X, "ky": Y, "lambda": L, "egrid": E, "species": S}
    return Export(path=Path("<fake>"), root=root, dims=dims, attrs={})


@pytest.mark.parametrize("name", ["linked_mixed", "linked_hybrid", "zero_trapped_parity_odd",
                                  "selfperiodic_np2"])
@pytest.mark.parametrize("with_checks", [True, False])
def test_from_export_fake(name, with_checks):
    """from_export reads the spec 09 names and reproduces ``build`` (also when the
    optional check arrays are absent and the masks are recomputed)."""
    g, data, params, gh_ref, out1_ref, out_ref = _synth(name)
    ex = _fake_export(g, with_checks=with_checks)
    d2, p2 = ir.from_export(ex)
    assert p2 == params
    assert ir.InvertParams.from_export(ex) == params
    flat1, tree1 = jax.tree_util.tree_flatten(data)
    flat2, tree2 = jax.tree_util.tree_flatten(d2)
    assert tree1 == tree2
    for a, b in zip(flat1, flat2):
        assert np.array_equal(np.asarray(a), np.asarray(b))
    if g.linked:
        assert lk.LinkedData.from_export(ex) is not None
    if g.def_parity:
        assert np.array_equal(par.ParityData.from_export(ex).par_it, data.parity.par_it)
    co = Coeffs.from_export(ex)._replace(g_h=gh_ref)
    got = ir.invert_rhs(g.source, g.phinew, g.bparnew, to_jax(co), to_jax(d2), p2,
                        gnew_in=g.gnew_in)
    assert_close(out_ref, got, 1e-12, f"{name} via from_export")


def _build_kwargs(g):
    return dict(
        ntgrid=g.ntgrid, ng2=g.ng2, aky=g.aky, kwork_filter=g.kwork, ieqzip=g.ieqzip,
        forbid=g.forbid, is_ttp=g.ttp, is_lower_bounce_point=g.lbp,
        is_upper_bounce_point=g.ubp, is_passing_hybrid_electron=g.hybrid,
        passing_wfb=g.passing_wfb, trapped_wfb=g.trapped_wfb, mixed_wfb=g.mixed_wfb,
        boundary_option_switch=g.bos, nonad_zero=g.nonad_zero, zero_forbid=g.zero_forbid,
        def_parity=g.def_parity, even=g.even, wfb=g.wfb, fphi=g.fphi, fbpar=g.fbpar,
        vperp2=g.vperp2, aj0=g.aj0, aj1=g.aj1, zt=g.zt, nonmaxw_corr=g.nonmaxw,
        l_links=g.l_links, r_links=g.r_links, itleft=g.itleft, itright=g.itright,
        wfb_factors=g.wfb_end if g.linked else None, box=g.box)


# =============================================================================
# Part 2: reference data (gpu/refdata), spec 11 section 6
# =============================================================================
# Every case exercises invert_rhs_1 (01 §C); the linked/parity tests reduce to identities or
# skip where the option is off. Parametrize without loading the exports at collection time.
ALL_CASES = [c for c in all_cases() if have(c)] or ["<no-refdata>"]

# probes that dump gnew / the new fields (spec 09 section 5.3); after_getfield holds the increment
GNEW_TAGS = ("after_parity", "after_vspace_derivatives", "after_hyper", "after_invert_rhs",
             "after_timeadv", "after_exb_shear", "advance_begin")
FIELD_TAGS = ("after_remove_zonal_flows", "after_field_update", "after_exb_shear", "advance_begin")


@lru_cache(maxsize=4)
def _setup(case):
    ex = export_of(case)
    data, params = ir.from_export(ex)
    return ex, to_jax(data), params, Dims.from_export(ex)


@lru_cache(maxsize=4)
def _jitted():
    return (jax.jit(ir.invert_rhs_1, static_argnames=("params",)),
            jax.jit(lk.apply_linked_boundary_conditions, static_argnames=("params",)),
            jax.jit(par.enforce_parity, static_argnames=("params",)),
            jax.jit(ir.init_homogeneous_g, static_argnames=("params",)))


def _before(tr, rec, tags):
    """Closest earlier record of the same (istep, attempt) whose tag is in ``tags``."""
    best = None
    for r in tr.records:
        if r.istep == rec.istep and r.attempt == rec.attempt and r.seq < rec.seq and r.tag in tags:
            if best is None or r.seq > best.seq:
                best = r
    if best is None:
        raise AssertionError(f"no record among {tags} before {rec}")
    return best


def _ladder_rungs(ex):
    if "dt_ladder" not in ex.root.groups:
        return []
    return sorted(int(k[1:]) for k in ex["/dt_ladder"].groups if k.startswith("k"))


_COEFF_CACHE: dict = {}


def _coeffs_for(case, ex, code_dt, data, params):
    """Coefficients at ``code_dt``: the export, a ladder rung with the same bit pattern, or
    the B1 recompute (``gs2gpu.coeffs``) plus ``init_homogeneous_g``; None if unavailable."""
    key = (case, float(code_dt))
    if key in _COEFF_CACHE:
        return _COEFF_CACHE[key]
    found = None
    base = Coeffs.from_export(ex)
    if base.code_dt == code_dt:
        found = base
    for k in _ladder_rungs(ex):
        if found is not None:
            break
        c = Coeffs.from_export(ex, rung=k)
        if c.code_dt == code_dt:
            found = c
    if found is None:
        try:
            from gs2gpu import coeffs as b1
            c = b1.compute_coeffs(code_dt, b1.CoeffInputs.from_export(ex),
                                  b1.CoeffParams.from_export(ex))
            found = c._replace(g_h=np.asarray(ir.init_homogeneous_g(c, data, params)))
        except Exception:  # noqa: BLE001  (B1 missing or failing: skip those records)
            found = None
    if found is not None:
        found = to_jax(found)
        if found.g_h is None:
            found = found._replace(g_h=ir.init_homogeneous_g(found, data, params))
    _COEFF_CACHE[key] = found
    return found


@pytest.mark.parametrize("case", ALL_CASES)
def test_masks_vs_export(case):
    """Recomputed masks/topology agree with the exporter's check arrays (spec 09 section 7.7)."""
    ex, data, params, dims = _setup(case)
    le = ex["/le_grids"]
    ng2 = int(le["ng2"])
    tw = bool(np.asarray(le["trapped_wfb"]).item())
    if "is_bounce_point" in le:
        bp = np.asarray(le["is_bounce_point"], bool)
        fb = np.asarray(le["forbid"], bool)
        for name, fn in (("is_lower_bounce_point", ir.lower_bounce_points),
                         ("is_upper_bounce_point", ir.upper_bounce_points)):
            if name in le:
                assert (fn(bp, fb, tw, ng2) == np.asarray(le[name], bool)).all(), name
    if "is_passing_hybrid_electron" in le and "spec_type" in ex["/species"]:
        stype = np.asarray(ex["/species/spec_type"])
        il = np.arange(1, ex.nlambda + 1)
        aky = np.asarray(ex["/kt_grids/aky"])
        hyb = ((stype == 4)[:, None, None] & (aky != 0.0)[None, :, None]
               & (il <= ng2 + 1)[None, None, :])
        assert (hyb == np.asarray(le["is_passing_hybrid_electron"], bool)).all()
    if not params.linked:
        return
    dfn = ex["/dist_fn"]
    il = np.arange(1, ex.nlambda + 1)
    ld = data.linked
    if "save_h" in dfn and ld is not None:
        save = np.asarray(dfn["save_h"], bool)                      # (S, E, L, Y, X, 2)
        # reconstruct save_h from the kernel masks: open masks need l_links/r_links > 0
        itl = np.asarray(dfn["itleft"])
        itr = np.asarray(dfn["itright"])
        aky = np.asarray(ex["/kt_grids/aky"])
        excl = (il > ng2 + 1) | ((il == ng2 + 1) & tw)
        left = (~excl)[:, None, None] & (itl >= 0)[None]
        right = (~excl)[:, None, None] & (itr >= 0)[None]
        nb = left | right
        s1 = left & (aky != 0)[None, :, None]
        s2 = right & (aky != 0)[None, :, None]
        pw = bool(np.asarray(le["passing_wfb"]).item())
        w = (il == ng2 + 1)[:, None, None] & nb & (not pw) & (aky != 0)[None, :, None]
        s1, s2 = s1 | w, s2 | w
        assert (save[..., 0] == s1[None, None]).all() and (save[..., 1] == s2[None, None]).all()
        if "connections_neighbor" in dfn:
            assert (np.asarray(dfn["connections_neighbor"], bool) == nb[None, None]).all()
    if "n_links" in dfn:
        nl = np.asarray(dfn["n_links"])                              # (Y, X, 2)
        for i, name in ((0, "l_links"), (1, "r_links")):
            links = np.asarray(dfn[name])
            assert (nl[..., i] == np.where(links == 0, 0, 2 * links - 1)).all()
    if ld is not None and "itleftmost" in dfn:
        lm = np.asarray(dfn["itleftmost"]) - 1
        rm = np.asarray(dfn["itrightmost"]) - 1
        for c in range(np.asarray(ld.chain_n).shape[0]):
            n = int(ld.chain_n[c])
            cells = np.asarray(ld.chain_x)[c, :n]
            y = int(ld.chain_y[c])
            assert (lm[y, cells] == cells[0]).all() and (rm[y, cells] == cells[-1]).all()
            if "n_cells" in dfn:
                assert (np.asarray(dfn["n_cells"])[y, cells] == n).all()


@pytest.mark.parametrize("case", ALL_CASES)
def test_init_homogeneous_g(case):
    """g_h recomputed from r vs the exported g_h, at the initial dt and every ladder rung."""
    ex, data, params, dims = _setup(case)
    f_gh = _jitted()[3]
    rungs = [None] + _ladder_rungs(ex)
    for rung in rungs:
        co = Coeffs.from_export(ex, rung=rung)
        if co.g_h is None:
            pytest.skip(f"{case}: g_h not exported")
        got = f_gh(to_jax(co._replace(g_h=None)), data, params=params)
        assert_close(co.g_h, got, RTOL_KERNEL, f"{case} g_h rung={rung} dt={co.code_dt}")


@pytest.mark.parametrize("case", ALL_CASES)
def test_invert_rhs_1(case):
    """invert_rhs_1 on the probe ``source`` vs ``after_invert_rhs_1``, every step and pass."""
    ex, data, params, dims = _setup(case)
    tr = trace_of(case)
    f_inv = _jitted()[0]
    recs = tr.find("source")
    if not recs:
        pytest.skip(f"{case}: no source probes in the trace")
    checked = 0
    for rec in recs:
        out_rec = tr.next(rec, "after_invert_rhs_1")
        assert out_rec is not None, rec
        src = tr.load(rec)["g_work"]
        ref = tr.load(out_rec)["gnew"]
        fields = fields_from_group(tr.load(_before(tr, rec, FIELD_TAGS)), dims)
        gnew_in = tr.load(_before(tr, rec, GNEW_TAGS))["gnew"]
        code_dt = float(tr.load(rec)["code_dt"])
        co = _coeffs_for(case, ex, code_dt, data, params)
        if co is None:
            continue
        got = f_inv(src, fields.phi, fields.bpar, co, data, params=params, gnew_in=gnew_in)
        assert_close(ref, got, RTOL_KERNEL,
                     f"{case} istep={rec.istep} attempt={rec.attempt} pass={rec.pass_}")
        checked += 1
    assert checked > 0, f"{case}: no record had coefficients available"


@pytest.mark.parametrize("case", ALL_CASES)
def test_linked_boundary_conditions(case):
    """apply_linked_boundary_conditions on ``after_invert_rhs_1`` vs ``after_invert_rhs``
    (the identity for non-linked cases)."""
    ex, data, params, dims = _setup(case)
    tr = trace_of(case)
    f_lnk = _jitted()[1]
    recs = tr.find("after_invert_rhs_1")
    if not recs:
        pytest.skip(f"{case}: no after_invert_rhs_1 probes")
    checked = 0
    for rec in recs:
        out_rec = tr.next(rec, "after_invert_rhs")
        assert out_rec is not None, rec
        gin = tr.load(rec)["gnew"]
        ref = tr.load(out_rec)["gnew"]
        fields = fields_from_group(tr.load(_before(tr, rec, FIELD_TAGS)), dims)
        co = _coeffs_for(case, ex, float(tr.load(rec)["code_dt"]), data, params)
        if co is None:
            continue
        got = f_lnk(gin, fields.phi, fields.bpar, co, data, params=params)
        assert_close(ref, got, RTOL_KERNEL,
                     f"{case} istep={rec.istep} attempt={rec.attempt} pass={rec.pass_}")
        checked += 1
    assert checked > 0


@pytest.mark.parametrize("case", ALL_CASES)
def test_strict_first_step(case, strict):
    """Strict mode (no FMA contraction by XLA, eager) on the first traced step: D1 tolerance
    for g_h, invert_rhs_1 and the linked BCs.

    Bitwise agreement is *not* expected. Despite ``-ffp-contract=off`` the reference build
    (GCC 15, ``-mcpu=apple-m1``) vectorises some complex products into ARMv8.3 ``fcmla``
    pairs, which fuse the multiply-adds per lane (25 objects, e.g. 6 in ``invert_rhs``, 4 in
    ``init_homogeneous_g``, 4 in ``apply_linked_boundary_conditions``); see
    :func:`test_g_h_rounding_model`. ``-mcpu=apple-m1+nofcma`` removes them."""
    ex, data, params, dims = _setup(case)
    tr = trace_of(case)
    base = Coeffs.from_export(ex)
    if base.g_h is None:
        pytest.skip(f"{case}: g_h not exported")
    gh = ir.init_homogeneous_g(to_jax(base._replace(g_h=None)), data, params)
    assert_close(base.g_h, gh, RTOL_KERNEL, f"{case} g_h strict")
    co = to_jax(base)
    recs = [r for r in tr.find("source") if float(tr.load(r)["code_dt"]) == base.code_dt]
    if not recs:
        pytest.skip(f"{case}: no source probe at the initial dt")
    for rec in recs[:2]:                                      # both passes of the first step
        fields = fields_from_group(tr.load(_before(tr, rec, FIELD_TAGS)), dims)
        gnew_in = tr.load(_before(tr, rec, GNEW_TAGS))["gnew"]
        out1 = ir.invert_rhs_1(tr.load(rec)["g_work"], fields.phi, fields.bpar, co, data,
                               params, gnew_in=gnew_in)
        r1 = tr.next(rec, "after_invert_rhs_1")
        assert_close(tr.load(r1)["gnew"], out1, RTOL_KERNEL,
                     f"{case} invert_rhs_1 strict pass={rec.pass_}")
        out = lk.apply_linked_boundary_conditions(tr.load(r1)["gnew"], fields.phi, fields.bpar,
                                                  co, data, params)
        assert_close(tr.load(tr.next(rec, "after_invert_rhs"))["gnew"], out, RTOL_KERNEL,
                     f"{case} linked strict pass={rec.pass_}")


_FMA = np.frompyfunc(__import__("math").fma, 3, 1)


def _fma(a, b, c):
    return np.asarray(_FMA(a, b, c), np.float64)


def _fcmla_product(n, m):
    """``n*m`` as GCC emits it with ``fcmla #0`` then ``fcmla #90`` into a zero accumulator
    (ARM ARM FCMLA: each lane is one fused multiply-add)."""
    re1 = _fma(n.real, m.real, np.zeros_like(n.real))
    im1 = _fma(n.real, m.imag, np.zeros_like(n.real))
    out = np.empty(np.broadcast(n, m).shape, np.complex128)
    out.real = _fma(n.imag, -m.imag, re1)
    out.imag = _fma(n.imag, m.real, im1)
    return out


@pytest.mark.parametrize("case", ALL_CASES)
def test_g_h_rounding_model(case, strict):
    """Pins down every last-bit difference of g_h against GS2 (dist_fn.fpp:2071-2114).

    Sign 2 (a scalar loop in the reference objects) equals strict JAX bit for bit. Sign 1
    (vectorised with ``fcmla``) equals an fcmla emulation of the same recurrence bit for bit,
    started from GS2's own ``g_h(-ntgrid,1)`` and lower-bounce values; every value equals
    either that emulation or the plain product (strict JAX), so the association order of
    the port is the Fortran one."""
    ex, data, params, dims = _setup(case)
    base = Coeffs.from_export(ex)
    if base.g_h is None:
        pytest.skip(f"{case}: g_h not exported")
    ref = np.asarray(base.g_h)
    got = np.asarray(ir.init_homogeneous_g(to_jax(base._replace(g_h=None)), data, params))
    _exact(ref[..., 1, :], got[..., 1, :], f"{case} g_h sign 2")
    r1 = np.moveaxis(np.asarray(base.r)[..., 0, :], -1, 0)          # (T, S, E, L, Y, X)
    ref1 = np.moveaxis(ref[..., 0, :], -1, 0)
    ref2 = np.moveaxis(ref[..., 1, :], -1, 0)
    ttp = np.asarray(data.ttp).T[:, None, None, :, None, None]
    lbsrc = np.asarray(data.lbsrc).T[:, None, None, :, None, None]
    emu = np.empty_like(ref1)
    emu[0] = ref1[0]
    for t in range(ref1.shape[0] - 1):
        g2 = np.where(lbsrc[t], ref2[t + 1], 0.0)
        nxt = g2 - _fcmla_product(r1[t], emu[t])                    # n = r (Vn), m = g_h (Vm)
        emu[t + 1] = np.where(ttp[t], 0.0, nxt)
    emu = np.moveaxis(emu, 0, -1)
    skip = np.asarray(data.gh_skip)[None, None, :, :, None, None]
    emu = np.where(skip, 0.0, emu)
    plain = got[..., 0, :]
    ok = (ref[..., 0, :] == emu) | (ref[..., 0, :] == plain)
    assert ok.all(), (f"{case}: {int((~ok).sum())} sign-1 g_h values match neither the fcmla "
                      f"emulation nor the plain product")


@pytest.mark.parametrize("case", ALL_CASES)
def test_parity(case):
    """enforce_parity at the timeadv call site (dist_fn.fpp:4091) vs ``after_parity``."""
    ex, data, params, dims = _setup(case)
    if not params.def_parity:
        pytest.skip(f"{case}: def_parity is off")
    tr = trace_of(case)
    f_par = _jitted()[2]
    recs = tr.find("after_parity")
    if not recs:
        pytest.skip(f"{case}: no after_parity probes")
    for rec in recs:
        src_rec = _before(tr, rec, ("after_vspace_derivatives", "after_hyper", "after_invert_rhs"))
        got = f_par(tr.load(src_rec)["gnew"], data, params=params)
        assert_close(tr.load(rec)["gnew"], got, RTOL_KERNEL,
                     f"{case} istep={rec.istep} attempt={rec.attempt} pass={rec.pass_}")
