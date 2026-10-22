"""Tests of the nonlinear term (gs2gpu.nonlinear) and the explicit-term / CFL bookkeeping
(gs2gpu.explicit). Specs 05, 02 section 3.7, 01 sections 1.14 / 3.1 / 3.5.

Two families:

* synthetic tests (always run): analytic Poisson brackets of Fourier modes against an
  independent direct-sum implementation, the GS2 index / normalisation conventions, literal
  loop versions of the loads and of ``nl_center``, the host bookkeeping, and the FFT
  round-off against FFTW itself (the library GS2 links, called through ctypes with GS2's plan
  flags; skipped if libfftw3 cannot be loaded);
* reference-data tests (skip without ``gpu/refdata``): ``add_nl`` vs probe
  ``nl_before_center``, ``nl_center`` vs ``nl_after_center``, the CFL signals vs ``nl_cfl``,
  the ``gexp`` history vs ``after_add_explicit_terms`` (both passes, every attempt), and the
  end-of-step / deferred CFL decisions of the dt-change cases.

Measured errors are printed (``pytest -s``) as ``[nl-err] ...`` lines.
"""

from __future__ import annotations

import ctypes
import ctypes.util
from functools import lru_cache

import jax
import jax.numpy as jnp
import numpy as np
import pytest

from conftest import RTOL_KERNEL, assert_close, cases_where, export_of, trace_of

import gs2gpu  # noqa: F401  (x64)
from gs2gpu import explicit as xp
from gs2gpu import nonlinear as nlm
from gs2gpu.gadjust import GAdjustData
from gs2gpu.numerics import max_abs_rel
from gs2gpu.types import Dims, Fields, TimeBlock, fields_from_group, to_jax

EPS = np.finfo(np.float64).eps


def _report(what, err):
    print(f"[nl-err] {what}: {err:.3e}")


# ============================================================================================
# FFTW oracle: GS2's transforms with the FFTW library and plan flags GS2 uses
# (FFTW_ESTIMATE + FFTW_UNALIGNED, externals/utils/fft_work.fpp:239-245; in-place 1D c2c along
# x, out-of-place 1D c2r/r2c along y, gs2_transforms.fpp:286-314), and add_nl in NumPy with
# one rounding per operation (the reference is built with -ffp-contract=off, 00_decisions D2).
# ============================================================================================

_FFTW_FORWARD, _FFTW_BACKWARD = -1, 1
_FFTW_FLAGS = (1 << 6) | (1 << 1)          # FFTW_ESTIMATE | FFTW_UNALIGNED


@lru_cache(maxsize=1)
def _fftw():
    for cand in ("/opt/homebrew/opt/fftw/lib/libfftw3.3.dylib", "/opt/homebrew/lib/libfftw3.dylib",
                 ctypes.util.find_library("fftw3")):
        if not cand:
            continue
        try:
            lib = ctypes.CDLL(cand)
        except OSError:
            continue
        P, I = ctypes.c_void_p, ctypes.c_int
        IP = ctypes.POINTER(I)
        lib.fftw_plan_many_dft.restype = P
        lib.fftw_plan_many_dft.argtypes = [I, IP, I, P, IP, I, I, P, IP, I, I, I, ctypes.c_uint]
        for name in ("fftw_plan_many_dft_c2r", "fftw_plan_many_dft_r2c"):
            getattr(lib, name).restype = P
            getattr(lib, name).argtypes = [I, IP, I, P, IP, I, I, P, IP, I, I, ctypes.c_uint]
        for name in ("fftw_execute_dft", "fftw_execute_dft_c2r", "fftw_execute_dft_r2c"):
            getattr(lib, name).restype = None
            getattr(lib, name).argtypes = [P, P, P]
        lib.fftw_destroy_plan.restype = None
        lib.fftw_destroy_plan.argtypes = [P]
        return lib
    return None


def _need_fftw():
    lib = _fftw()
    if lib is None:
        pytest.skip("libfftw3 not loadable (FFTW comparison only)")
    return lib


def _ints(*v):
    return (ctypes.c_int * len(v))(*v)


def _fftw_c2c_inplace(x, sign):
    lib = _need_fftw()
    x = np.ascontiguousarray(x, np.complex128).copy()
    n = x.shape[-1]
    plan = lib.fftw_plan_many_dft(1, _ints(n), x.size // n, x.ctypes.data, _ints(n), 1, n,
                                  x.ctypes.data, _ints(n), 1, n, sign, _FFTW_FLAGS)
    data = np.ascontiguousarray(x).copy()        # planning may scribble on the array
    lib.fftw_execute_dft(plan, data.ctypes.data, data.ctypes.data)
    lib.fftw_destroy_plan(plan)
    return data


def _fftw_c2r(x, n):
    lib = _need_fftw()
    x = np.ascontiguousarray(x, np.complex128).copy()
    nc = n // 2 + 1
    how = x.size // nc
    din = np.zeros_like(x)
    dout = np.zeros(x.shape[:-1] + (n,))
    plan = lib.fftw_plan_many_dft_c2r(1, _ints(n), how, din.ctypes.data, _ints(nc), 1, nc,
                                      dout.ctypes.data, _ints(n), 1, n, _FFTW_FLAGS)
    out = np.zeros(x.shape[:-1] + (n,))
    lib.fftw_execute_dft_c2r(plan, x.ctypes.data, out.ctypes.data)
    lib.fftw_destroy_plan(plan)
    return out


def _fftw_r2c(x):
    lib = _need_fftw()
    x = np.ascontiguousarray(x, np.float64).copy()
    n = x.shape[-1]
    nc = n // 2 + 1
    how = x.size // n
    din = np.zeros_like(x)
    dout = np.zeros(x.shape[:-1] + (nc,), np.complex128)
    plan = lib.fftw_plan_many_dft_r2c(1, _ints(n), how, din.ctypes.data, _ints(n), 1, n,
                                      dout.ctypes.data, _ints(nc), 1, nc, _FFTW_FLAGS)
    out = np.zeros(x.shape[:-1] + (nc,), np.complex128)
    lib.fftw_execute_dft_r2c(plan, x.ctypes.data, out.ctypes.data)
    lib.fftw_destroy_plan(plan)
    return out


def _xla_c2c_inplace(x, sign):
    x = jnp.asarray(x)
    if sign == _FFTW_BACKWARD:
        return np.asarray(jnp.conj(jnp.fft.fft(jnp.conj(x), axis=-1)))
    return np.asarray(jnp.fft.fft(x, axis=-1))


def _xla_c2r(x, n):
    return np.asarray(jnp.fft.irfft(jnp.asarray(x), n=n, axis=-1, norm="forward"))


def _xla_r2c(x):
    return np.asarray(jnp.fft.rfft(jnp.asarray(x), axis=-1))


#: FFT primitives used by the replica: FFTW (GS2's library) or XLA's (to isolate the FFTs).
_ORACLE_FFT = {"c2c": _fftw_c2c_inplace, "c2r": _fftw_c2r, "r2c": _fftw_r2c}


class _oracle_ffts:
    """Context manager: run the replica with XLA's FFTs instead of FFTW's."""

    def __enter__(self):
        self.saved = dict(_ORACLE_FFT)
        _ORACLE_FFT.update(c2c=_xla_c2c_inplace, c2r=_xla_c2r, r2c=_xla_r2c)

    def __exit__(self, *exc):
        _ORACLE_FFT.update(self.saved)


def _cplx(re, im):
    out = np.empty(np.broadcast(re, im).shape, np.complex128)
    out.real = re
    out.imag = im
    return out


def _np_rmul(z, r):
    return _cplx(z.real * r, z.imag * r)


def oracle_transform2(G, p: nlm.NLParams, forbid=None):
    """transform2_5d with FFTW (GS2's data flow; forbidden lines are not gathered)."""
    S, E, L, Y, X, _, T = G.shape
    nx, ny = p.nx, p.ny
    half = np.where(np.arange(Y) == 0, 1.0, 0.5)[None, None, None, :, None, None, None]
    Gh = _np_rmul(G, half)
    xxf = np.zeros((S, E, L, 2, T, Y, nx), np.complex128)          # lines (ik, ig, isgn, il, ...)
    J = nlm.kx_fft_index(X, nx)
    xxf[..., J] = np.transpose(Gh, (0, 1, 2, 5, 6, 3, 4))
    if forbid is not None:
        xxf = np.where(forbid[None, None, :, None, :, None, None], 0.0, xxf)
    xxf = _ORACLE_FFT["c2c"](xxf.reshape(-1, nx), _FFTW_BACKWARD).reshape(xxf.shape)
    fft = np.zeros((S, E, L, 2, T, nx, ny // 2 + 1), np.complex128)
    fft[..., :Y] = np.swapaxes(xxf, -1, -2)
    fft[..., 0] = fft[..., 0].real              # FFTW's c2r ignores Im(DC) (explicit for XLA)
    return _ORACLE_FFT["c2r"](fft.reshape(-1, ny // 2 + 1), ny).reshape(S, E, L, 2, T, nx, ny)


def oracle_inverse2(f, p: nlm.NLParams):
    """inverse2_5d with FFTW, including the scale of gs2_transforms.fpp:968-983."""
    S, E, L, _, T, nx, ny = f.shape
    Y, X = p.naky, p.ntheta0
    fft = _ORACLE_FFT["r2c"](f.reshape(-1, ny)).reshape(S, E, L, 2, T, nx, ny // 2 + 1)
    xxf = np.ascontiguousarray(np.swapaxes(fft[..., :Y], -1, -2))  # (..., Y, nx)
    xxf = _ORACLE_FFT["c2c"](xxf.reshape(-1, nx), _FFTW_FORWARD).reshape(xxf.shape)
    G = xxf[..., nlm.kx_fft_index(X, nx)]                          # (S,E,L,2,T,Y,X)
    G = np.transpose(G, (0, 1, 2, 5, 6, 3, 4))
    scale = 1.0 / float(nx * ny) if p.accelerated else (1.0 / float(nx)) * (1.0 / float(ny))
    fac = np.where(np.arange(Y) == 0, scale, scale * 2.0)[None, None, None, :, None, None, None]
    return _np_rmul(G, fac)


def _np_load(fields, k6, D, p):
    """load_k*_phi/bpar/apar, one rounding per operation."""
    aj0, aj1, vperp2 = D.gadj.aj0, D.gadj.aj1, D.gadj.vperp2
    if p.phi_on:
        c = k6 * aj0
        ph = fields.phi[None, None, None]
        g6 = _cplx((-(c * ph.imag)) * p.fphi, (c * ph.real) * p.fphi)
    else:
        g6 = np.zeros(aj0.shape, np.complex128)
    if p.bpar_on:
        c = (((k6 * aj1) * 2.0) * vperp2) * D.tz[:, None, None, None, None, None]
        b = fields.bpar[None, None, None]
        fr, fi = (-(c * b.imag)) * p.fbpar, (c * b.real) * p.fbpar
        g6 = _cplx(g6.real + fr, g6.imag + fi)
    g7 = np.broadcast_to(g6[..., None, :], aj0.shape[:-1] + (2, aj0.shape[-1]))
    if p.apar_on:
        c = (k6 * aj0) * D.stm[:, None, None, None, None, None]
        a = fields.apar[None, None, None]
        fr = ((-(c * a.imag)) * p.fapar)[..., None, :]
        fi = ((c * a.real) * p.fapar)[..., None, :]
        g7 = _cplx(g7.real - fr * D.vpa, g7.imag - fi * D.vpa)
    return np.array(g7)


def _np_g_adjust(g, phi, bpar, facphi, facbpar, gd):
    """g_adjust_floats (dist_fn_arrays.fpp:117-175), one rounding per operation."""
    pf = (gd.zt[:, None] * gd.nonmaxw_corr) * facphi if facphi != 0 else np.zeros_like(gd.nonmaxw_corr)
    bf = 2.0 * facbpar if facbpar != 0 else 0.0
    r1 = (bf * gd.vperp2) * gd.aj1
    b = bpar[None, None, None]
    ph = phi[None, None, None]
    pf6 = pf[:, :, None, None, None, None]
    adj = _cplx(r1 * b.real + ((pf6 * ph.real) * gd.aj0), r1 * b.imag + ((pf6 * ph.imag) * gd.aj0))
    out = _cplx(g.real + adj.real[..., None, :], g.imag + adj.imag[..., None, :])
    mask = gd.adjust_mask[None, None, :, None, None, None, :]
    return np.where(mask, out, g)


def oracle_add_nl(g, fields, D, p, forbid=None):
    """add_nl (nonlinear_terms.fpp:788-929) with FFTW transforms and NumPy arithmetic."""
    akx6 = D.akx[None, None, None, None, :, None]
    aky6 = D.aky[None, None, None, :, None, None]
    ba = oracle_transform2(_np_load(fields, akx6, D, p), p, forbid)
    mvy = np.max(np.abs(ba)) * p.cfly
    h = _np_g_adjust(g, fields.phi, fields.bpar, p.fphi, p.fbpar, D.gadj)
    k = D.aky[None, None, None, :, None, None, None]
    gb = oracle_transform2(_cplx(-h.imag * k, h.real * k), p, forbid)
    bracket = ba * gb
    ba = oracle_transform2(_np_load(fields, aky6, D, p), p, forbid)
    mvx = np.max(np.abs(ba)) * p.cflx
    k = D.akx[None, None, None, None, :, None, None]
    gb = oracle_transform2(_cplx(-h.imag * k, h.real * k), p, forbid)
    bracket = bracket - ba * gb
    nl = _np_rmul(oracle_inverse2(bracket, p), p.kxfac)
    return nl, np.array([mvx, mvy, EPS])


def oracle_nl_center(g, bd, D, p):
    """nl_center as the literal in-place ascending loop of nonlinear_terms.fpp:705-758."""
    g = np.array(g, np.complex128, copy=True)
    S, E, L, Y, X, _, T = g.shape
    for il in range(L):            # cell_zero = trapped .and. (il > jend(ig) .or. il > jend(ig+1))
        if p.nl_forbid_force_zero:
            fb = np.asarray(D.forbid)[il]
            g[:, :, il, :, :, :, fb] = 0.0
        for ig in range(T - 1):
            if D.is_ttp[il, ig]:
                continue
            if D.cell_zero[il, ig]:
                g[:, :, il, :, :, :, ig] = 0.0
            else:
                a1, b1 = g[:, :, il, :, :, 0, ig + 1], g[:, :, il, :, :, 0, ig]
                a2, b2 = g[:, :, il, :, :, 1, ig + 1], g[:, :, il, :, :, 1, ig]
                g[:, :, il, :, :, 0, ig] = _cplx((1. + bd) * a1.real + (1. - bd) * b1.real,
                                                 (1. + bd) * a1.imag + (1. - bd) * b1.imag)
                g[:, :, il, :, :, 1, ig] = _cplx((1. - bd) * a2.real + (1. + bd) * b2.real,
                                                 (1. - bd) * a2.imag + (1. + bd) * b2.imag)
    return g


# ============================================================================================
# synthetic setups
# ============================================================================================


def _signed_kx(X, nx):
    J = nlm.kx_fft_index(X, nx)
    return np.where(J <= nx // 2, J, J - nx)


def _hermitian_ky0(F, kxi):
    """Make the ky = 0 row (axis -2 = Y, axis -1 = X) of ``F`` Hermitian in kx."""
    F = np.array(F, copy=True)
    pos = {int(k): i for i, k in enumerate(kxi)}
    for it, k in enumerate(kxi):
        jt = pos[-int(k)]
        if jt == it:
            F[..., 0, it] = F[..., 0, it].real
        elif jt > it:
            F[..., 0, jt] = np.conj(F[..., 0, it])
    return F


def synthetic(seed=0, S=2, E=2, L=5, Y=4, X=7, ntgrid=3, nx=12, ny=12, em=True, trapped=True,
              bd=0.05, kxfac=1.37, hermitian=True, accelerated=False, include=(True, True, True)):
    """A small but complete nonlinear setup with GS2-like structure (forbidden regions with
    exact zeros, ttp points, jend), random state; returns (g, Fields, NLData (NumPy), NLParams)."""
    rng = np.random.default_rng(seed)
    T = 2 * ntgrid + 1
    dkx, dky = 0.31, 0.27
    kxi = _signed_kx(X, nx)
    akx, aky = kxi * dkx, np.arange(Y) * dky
    if trapped:
        ng2 = L - 2
        jend = np.array([max(L - abs(t - ntgrid), ng2) for t in range(T)])   # all il allowed at centre
        jend[0] = jend[-1] = ng2
    else:
        jend = np.zeros(T, int)
    il = np.arange(1, L + 1)[:, None]
    forbid = (il > jend[None, :]) if trapped else np.zeros((L, T), bool)
    is_ttp = np.zeros((L, T), bool)
    if trapped:
        is_ttp[L - 1, ntgrid] = True            # totally trapped at the outboard midplane
        is_ttp[L - 2, ntgrid - 1] = True
    cell_zero = trapped & ((il > jend[None, :-1]) | (il > jend[None, 1:]))
    shape6 = (S, E, L, Y, X, T)
    ok6 = ~forbid[None, None, :, None, None, :]
    aj0 = np.where(ok6, rng.uniform(0.3, 1.0, shape6), 0.0)
    aj1 = np.where(ok6, rng.uniform(0.1, 0.5, shape6), 0.0)
    vperp2 = np.where(ok6, rng.uniform(0.2, 2.0, shape6), 0.0)
    vpa1 = np.where(ok6, rng.uniform(0.1, 1.5, shape6), 0.0)
    vpa = np.stack([vpa1, -vpa1], axis=-2)
    gadj = GAdjustData(vperp2=vperp2, aj0=aj0, aj1=aj1, zt=np.array([1.0, -2.0])[:S],
                       nonmaxw_corr=rng.uniform(0.9, 1.1, (S, E)),
                       adjust_mask=~(trapped & (il > jend[None, :])))
    data = nlm.NLData(akx=akx, aky=aky, vpa=vpa, stm=np.array([1.0, 7.5])[:S],
                      tz=np.array([1.0, -0.5])[:S], forbid=forbid, is_ttp=is_ttp,
                      cell_zero=cell_zero, gadj=gadj)
    cfl = 0.5
    params = nlm.NLParams(
        nonlin=True, nx=nx, ny=ny, naky=Y, ntheta0=X, accelerated=accelerated, kxfac=kxfac,
        cfl=cfl, cflx=(float(np.max(akx)) * 0.5) / cfl, cfly=(float(np.max(aky)) * 0.5) / cfl,
        fphi=1.0, fapar=1.0 if em else 0.0, fbpar=1.0 if em else 0.0,
        phi_on=include[0], apar_on=em and include[1], bpar_on=em and include[2],
        nl_forbid_force_zero=True, trapped=trapped, bd=bd)

    def cfield():
        f = rng.standard_normal((Y, X, T)) + 1j * rng.standard_normal((Y, X, T))
        if hermitian:
            f = np.moveaxis(_hermitian_ky0(np.moveaxis(f, -1, 0), kxi), 0, -1)
        return f
    phi = cfield()
    apar = 0.3 * cfield() if em else np.zeros((Y, X, T), np.complex128)
    bpar = 0.2 * cfield() if em else np.zeros((Y, X, T), np.complex128)
    g = rng.standard_normal((S, E, L, Y, X, 2, T)) + 1j * rng.standard_normal((S, E, L, Y, X, 2, T))
    if hermitian:
        g = np.moveaxis(_hermitian_ky0(np.moveaxis(g, (-2, -1), (0, 1)), kxi), (0, 1), (-2, -1))
    return g, Fields(phi, apar, bpar), data, params


def _jfields(f: Fields):
    return Fields(*(jnp.asarray(x) for x in f))


def _eager(fn, *args, **kw):
    """Run op-by-op (no XLA fusion, hence no FMA contraction)."""
    with jax.disable_jit():
        return fn(*args, **kw)


# ============================================================================================
# synthetic tests: conventions, analytic brackets, literal loops
# ============================================================================================


@pytest.mark.parametrize("ntheta0,nx", [(11, 16), (7, 12), (5, 8), (21, 32), (1, 4)])
def test_kx_fft_index(ntheta0, nx):
    """J(it) of gs2_transforms.fpp:485-489 (1-based it, it0 = it + nx - ntheta0 past the middle)."""
    J = nlm.kx_fft_index(ntheta0, nx)
    for it in range(1, ntheta0 + 1):
        it0 = it + nx - ntheta0 if it > (ntheta0 + 1) // 2 else it    # Fortran integer division
        assert J[it - 1] == it0 - 1
    assert len(set(J.tolist())) == ntheta0


@pytest.mark.parametrize("nx,ny,Y,X", [(16, 16, 6, 11), (12, 12, 4, 7), (12, 15, 5, 7), (8, 8, 3, 5)])
def test_transform_roundtrip(nx, ny, Y, X):
    """inverse2(transform2(G)) = G for ky > 0 and the Hermitian projection of the ky = 0 row
    (FFTW's c2r drops Im of the x-transformed DC bin, spec 05 P6)."""
    g, _, _, p = synthetic(seed=1, Y=Y, X=X, nx=nx, ny=ny, hermitian=False)
    G = jnp.asarray(g)
    back = np.asarray(nlm.inverse2_jit(nlm.transform2_jit(G, p), p))
    kxi = _signed_kx(X, nx)
    pos = {int(k): i for i, k in enumerate(kxi)}
    proj = np.array(g, copy=True)
    for it, k in enumerate(kxi):
        proj[:, :, :, 0, it] = 0.5 * (g[:, :, :, 0, it] + np.conj(g[:, :, :, 0, pos[-int(k)]]))
    err = max_abs_rel(proj, back)
    _report(f"round trip nx={nx} ny={ny}", err)
    assert err < 1e-14


def _direct_bracket(g, fields: Fields, D, p):
    """Independent pseudo-spectral reference by direct summation (no FFT): GS2 coefficients
    -> real space on the nx*ny grid -> bracket -> GS2 coefficients."""
    Y, X, nx, ny = p.naky, p.ntheta0, p.nx, p.ny
    kxi = _signed_kx(X, nx)
    ex = np.exp(2j * np.pi * np.outer(kxi, np.arange(nx)) / nx)       # (X, nx)
    ey = np.exp(2j * np.pi * np.outer(np.arange(Y), np.arange(ny)) / ny)

    def real(F7):
        """G7 coefficients -> (S,E,L,2,T,nx,ny) real values with GS2's convention:
        Re sum_kx F[0,kx] e^{i kx x} + sum_{kx, ky>0} 2 Re((F/2) e^{i(kx x + ky y)})."""
        F = np.transpose(F7, (0, 1, 2, 5, 6, 3, 4))                         # (..., Y, X)
        f0 = np.einsum("...x,xi->...i", F[..., 0, :], ex).real[..., None]
        fk = np.einsum("...yx,xi,yj->...ij", F[..., 1:, :], ex, ey[1:]).real
        return f0 + fk

    k7x = D.akx[None, None, None, None, :, None, None]
    k7y = D.aky[None, None, None, :, None, None, None]
    # chi itself (no i k): chi = J0 phi fphi + 2 vperp2 J1 tz bpar fbpar - J0 stm vpa apar fapar
    gd = D.gadj
    chi6_phi = (gd.aj0 * fields.phi[None, None, None]) * p.fphi if p.phi_on else 0.0
    chi6_b = ((((gd.aj1 * 2.0) * gd.vperp2) * D.tz[:, None, None, None, None, None])
              * fields.bpar[None, None, None]) * p.fbpar if p.bpar_on else 0.0
    chi7 = np.broadcast_to(np.asarray(chi6_phi + chi6_b)[..., None, :], g.shape).astype(np.complex128)
    if p.apar_on:
        chi7 = chi7 - ((gd.aj0 * D.stm[:, None, None, None, None, None])
                       * fields.apar[None, None, None] * p.fapar)[..., None, :] * D.vpa
    h = _np_g_adjust(g, fields.phi, fields.bpar, p.fphi, p.fbpar, gd)
    cx, cy = real(1j * k7x * chi7), real(1j * k7y * chi7)
    hx, hy = real(1j * k7x * h), real(1j * k7y * h)
    B = cx * hy - cy * hx
    Gout = np.einsum("...ij,xi,yj->...yx", B, np.conj(ex), np.conj(ey)) / (nx * ny)
    Gout[..., 1:, :] *= 2.0
    nl = p.kxfac * np.transpose(Gout, (0, 1, 2, 5, 6, 3, 4))
    return nl, np.array([np.max(np.abs(cy)) * p.cflx, np.max(np.abs(cx)) * p.cfly, EPS])


@pytest.mark.parametrize("em", [False, True])
@pytest.mark.parametrize("nx,ny,Y,X", [(16, 16, 6, 11), (12, 15, 5, 7)])
def test_add_nl_vs_direct_sum(em, nx, ny, Y, X):
    """The full add_nl (loads, g_adjust, i*k, padding, 0.5/2 convention, bracket sign, scale,
    kxfac, CFL maxima) against an FFT-free direct-summation pseudo-spectral bracket."""
    g, f, D, p = synthetic(seed=2, Y=Y, X=X, nx=nx, ny=ny, em=em)
    ref, ref_mv = _direct_bracket(g, f, D, p)
    nl, mv = nlm.add_nl_jit(jnp.asarray(g), _jfields(f), to_jax(D), p)
    err = max_abs_rel(ref, np.asarray(nl))
    _report(f"add_nl vs direct sum em={em} nx={nx} ny={ny}", err)
    assert err < 1e-13
    np.testing.assert_allclose(np.asarray(mv), ref_mv, rtol=1e-13)


def test_single_mode_bracket_analytic():
    """{cos(k1.r), cos(k2.r)} = (kx1 ky2 - ky1 kx2) sin(k1.r) sin(k2.r) in GS2 coefficients."""
    S = E = L = 1
    Y, X, T, nx, ny = 6, 11, 1, 16, 16
    g, f, D, p = synthetic(seed=3, S=1, E=1, L=1, Y=Y, X=X, ntgrid=0, nx=nx, ny=ny, em=False,
                           trapped=False, kxfac=1.0)
    gd = D.gadj._replace(aj0=np.ones_like(D.gadj.aj0), zt=np.array([0.0]))  # chi = phi, h = g
    D = D._replace(gadj=gd)
    kxi = _signed_kx(X, nx)
    it = {int(k): i for i, k in enumerate(kxi)}
    # chi = cos(kx1 x + ky1 y) with (kx1, ky1) = (2, 1): GS2 coefficient 1 at (ik=1, kx=2)
    phi = np.zeros((Y, X, T), np.complex128)
    phi[1, it[2]] = 1.0
    # h = cos(kx2 x + ky2 y) with (kx2, ky2) = (-1, 2)
    gg = np.zeros((S, E, L, Y, X, 2, T), np.complex128)
    gg[:, :, :, 2, it[-1]] = 1.0
    nl, mv = nlm.add_nl_jit(jnp.asarray(gg), _jfields(Fields(phi, 0 * phi, 0 * phi)), to_jax(D), p)
    nl = np.asarray(nl)[0, 0, 0, :, :, 0, 0]
    kx1, ky1, kx2, ky2 = 2 * 0.31, 1 * 0.27, -1 * 0.31, 2 * 0.27
    amp = kx1 * ky2 - ky1 * kx2
    # sin a sin b = (cos(a-b) - cos(a+b))/2 ; a-b: (3, -1) -> cos((-3, 1)): coefficient +amp/2 at
    # (ik=1, kx=-3); a+b: (1, 3): coefficient -amp/2 at (ik=3, kx=1)
    expected = np.zeros((Y, X), np.complex128)
    expected[1, it[-3]] = 0.5 * amp
    expected[3, it[1]] = -0.5 * amp
    assert np.max(np.abs(nl - expected)) < 1e-14 * abs(amp) * 10
    # CFL maxima: |d_x chi| max = |kx1|, |d_y chi| max = |ky1| (cos has unit amplitude)
    np.testing.assert_allclose(np.asarray(mv)[:2], [ky1 * p.cflx, kx1 * p.cfly], rtol=1e-14)


def test_loads_literal():
    """load_k{x,y}_{phi,bpar,apar}: bitwise against the one-rounding-per-op NumPy version when
    run op by op; within 2 ulp-scale when XLA may contract FMAs."""
    g, f, D, p = synthetic(seed=4, em=True)
    for name, k6 in (("x", D.akx[None, None, None, None, :, None]),
                     ("y", D.aky[None, None, None, :, None, None])):
        ref = _np_load(f, k6, D, p)
        fn = nlm.load_kx_chi if name == "x" else nlm.load_ky_chi
        got = np.asarray(_eager(fn, _jfields(f), to_jax(D), p))
        assert np.array_equal(ref, got), f"load_k{name} differs (op by op)"
        got_jit = np.asarray(jax.jit(fn, static_argnums=2)(_jfields(f), to_jax(D), p))
        assert max_abs_rel(ref, got_jit) < 1e-15


def test_load_switches():
    """phi off -> zero_array before the bpar/apar loads; include_* switches respected."""
    g, f, D, p = synthetic(seed=5, em=True, include=(False, True, False))
    got = np.asarray(nlm.load_kx_chi(_jfields(f), to_jax(D), p))
    ref = _np_load(f, D.akx[None, None, None, None, :, None], D, p)
    assert np.array_equal(got, ref)
    p0 = p._replace(apar_on=False)
    assert not np.any(np.asarray(nlm.load_kx_chi(_jfields(f), to_jax(D), p0)))


@pytest.mark.parametrize("bd", [0.0, 0.05])
@pytest.mark.parametrize("force_zero", [True, False])
def test_nl_center_literal(bd, force_zero):
    """nl_center vs the literal in-place ascending loop (bitwise op by op)."""
    g, f, D, p = synthetic(seed=6, bd=bd)
    p = p._replace(nl_forbid_force_zero=force_zero)
    ref = oracle_nl_center(g, bd, D, p)
    got = np.asarray(_eager(nlm.nl_center, jnp.asarray(g), bd, to_jax(D), p))
    assert np.array_equal(ref, got)
    got_jit = np.asarray(nlm.nl_center_jit(jnp.asarray(g), bd, to_jax(D), p))
    assert max_abs_rel(ref, got_jit) < 1e-15
    # structure: ttp points keep the grid value, ig = ntgrid keeps the grid value
    ttp = np.argwhere(D.is_ttp)
    for il, ig in ttp:
        np.testing.assert_array_equal(got[:, :, il, :, :, :, ig], g[:, :, il, :, :, :, ig])
    ok = ~D.forbid[:, -1]
    np.testing.assert_array_equal(got[:, :, ok, :, :, :, -1], g[:, :, ok, :, :, :, -1])


def test_jit_traces_and_shapes():
    g, f, D, p = synthetic(seed=7)
    out, mv = nlm.nl_source_jit(jnp.asarray(g), _jfields(f), to_jax(D), p)
    assert out.shape == g.shape and out.dtype == jnp.complex128
    assert mv.shape == (3,) and float(mv[2]) == EPS
    ba = nlm.transform2_jit(jnp.asarray(g), p)
    assert ba.shape == (2, 2, 5, 2, 7, 12, 12) and ba.dtype == jnp.float64


def test_forbidden_lines_are_zero():
    """Forbidden lines give exact zeros whether or not they are transformed (spec 05 4.5.3)."""
    g, f, D, p = synthetic(seed=8)
    nl, _ = nlm.add_nl_jit(jnp.asarray(g), _jfields(f), to_jax(D), p)
    nl = np.asarray(nl)
    fb = D.forbid
    assert fb.any()
    for il in range(fb.shape[0]):
        assert np.all(nl[:, :, il, :, :, :, fb[il]] == 0)


# ============================================================================================
# FFT round-off against FFTW (the library GS2 links)
# ============================================================================================


@pytest.mark.parametrize("nx,ny,Y,X", [(16, 16, 6, 11), (8, 8, 3, 5), (12, 15, 5, 7), (48, 48, 16, 31)])
def test_transforms_vs_fftw(nx, ny, Y, X):
    _need_fftw()
    g, f, D, p = synthetic(seed=9, Y=Y, X=X, nx=nx, ny=ny, L=3, E=2, S=1, ntgrid=2, hermitian=False)
    ref = oracle_transform2(g, p)
    got = np.asarray(nlm.transform2_jit(jnp.asarray(g), p))
    err = max_abs_rel(ref, got)
    _report(f"transform2 JAX vs FFTW nx={nx} ny={ny}", err)
    assert err < 1e-14
    refi = oracle_inverse2(ref, p)
    goti = np.asarray(nlm.inverse2_jit(jnp.asarray(ref), p))
    erri = max_abs_rel(refi, goti)
    _report(f"inverse2 JAX vs FFTW nx={nx} ny={ny}", erri)
    assert erri < 1e-14
    for accel in (True,):
        pa = p._replace(accelerated=accel)
        np.testing.assert_allclose(np.asarray(nlm.inverse2_jit(jnp.asarray(ref), pa)),
                                   oracle_inverse2(ref, pa), rtol=0, atol=1e-14 * np.max(np.abs(refi)))


@pytest.mark.parametrize("em", [False, True])
@pytest.mark.parametrize("nx,ny,Y,X", [(16, 16, 6, 11), (8, 8, 3, 5)])
def test_add_nl_vs_fftw_oracle(em, nx, ny, Y, X):
    """Full add_nl + nl_center against the FFTW/NumPy replica of GS2's arithmetic: the residual
    is FFT-library (DUCC vs FFTW) plus FMA round-off."""
    _need_fftw()
    g, f, D, p = synthetic(seed=10, Y=Y, X=X, nx=nx, ny=ny, em=em)
    ref, ref_mv = oracle_add_nl(g, f, D, p, forbid=D.forbid)
    nl, mv = nlm.add_nl_jit(jnp.asarray(g), _jfields(f), to_jax(D), p)
    err = max_abs_rel(ref, np.asarray(nl))
    _report(f"add_nl JAX(jit) vs FFTW replica em={em} nx={nx}", err)
    assert err < 1e-13
    err_mv = float(np.max(np.abs(np.asarray(mv) - ref_mv) / ref_mv))
    _report(f"max_vel_components JAX vs FFTW replica em={em} nx={nx}", err_mv)
    assert err_mv < 1e-14
    refc = oracle_nl_center(ref, p.bd, D, p)
    gotc = np.asarray(nlm.nl_center_jit(nl, p.bd, to_jax(D), p))
    errc = max_abs_rel(refc, gotc)
    _report(f"nl_center(add_nl) JAX vs FFTW replica em={em} nx={nx}", errc)
    assert errc < 1e-13


@pytest.mark.parametrize("em", [False, True])
def test_add_nl_arithmetic_bitwise(em):
    """With the replica's FFTs swapped for XLA's, the replica (a one-rounding-per-operation
    transcription of the Fortran) and the JAX kernel run op by op agree BITWISE, and strict-mode
    jit (numerics.set_strict) equals op by op: every non-FFT operation is in GS2's order and
    every FMA-contractible product is protected. Hence all JAX-vs-GS2 differences in add_nl
    come from the FFT library (and from XLA's FMA contraction outside strict mode)."""
    from gs2gpu import numerics
    g, f, D, p = synthetic(seed=12, em=em, nx=16, ny=16, Y=6, X=11)
    Dj, G, F = to_jax(D), jnp.asarray(g), _jfields(f)
    with _oracle_ffts():
        ref, ref_mv = oracle_add_nl(g, f, D, p, forbid=D.forbid)
    eager, mv = _eager(nlm.add_nl, G, F, Dj, p)
    assert np.array_equal(ref, np.asarray(eager))
    assert np.array_equal(ref_mv, np.asarray(mv))
    numerics.set_strict(True)
    try:
        strict = jax.jit(lambda a, b, c: nlm.nl_source(a, b, c, p))    # fresh trace
        s, _ = strict(G, F, Dj)
    finally:
        numerics.set_strict(False)
    assert np.array_equal(oracle_nl_center(ref, p.bd, D, p), np.asarray(s))


# ============================================================================================
# explicit-term bookkeeping and CFL control (synthetic)
# ============================================================================================


def _time(code_dt=0.05, prev1=-1.0, prev2=-1.0, cfl=1e8, dt_max=0.05):
    return TimeBlock(code_dt=code_dt, user_dt=code_dt, code_dt_prev1=prev1, code_dt_prev2=prev2,
                     code_dt_old=0.0, code_time=0.0, user_time=0.0, code_dt_cfl=cfl,
                     user_dt_cfl=cfl, code_dt_min=1e-6, code_dt_max=dt_max)


def _xparams(**kw):
    base = dict(nonlin=True, split_nonlinear=False, zip=False, use_cfl_limit=True, use_2d_cfl=True,
                use_order_based_error=False, error_target=0.1, istep_error_start=30, cfl=0.5,
                immediate_reset=True, delt_adj=2.0, delt_cushion=1.5)
    base.update(kw)
    return xp.ExplicitParams(**base)


def _stub_kernel(mvx, mvy):
    calls = []

    def kernel(g, fields, data, params):
        calls.append(1)
        return jnp.full(g.shape, complex(len(calls), 0.0)), jnp.asarray([mvx, mvy, EPS])
    kernel.calls = calls
    return kernel


def _state0(shape=(1, 1, 1, 2, 3, 2, 3)):
    z = jnp.zeros(shape, jnp.complex128)
    return xp.ExplicitState(z, z, z, istep_last=0)


def test_history_shuffle_and_second_pass():
    """gexp_3 <- gexp_2 <- gexp_1 <- NL once per istep; the second pass is a no-op."""
    kern = _stub_kernel(1.0, 2.0)
    st = _state0()
    g = jnp.zeros(st.gexp_1.shape, jnp.complex128)
    f = Fields(*(jnp.zeros((2, 3, 3), jnp.complex128),) * 3)
    t = _time()
    for istep in (1, 2, 3):
        out = xp.add_explicit_terms(st, g, f, istep, t, False, None, None, _xparams(), kernel=kern)
        assert out.computed and not out.reset and not out.abort_timeadv
        assert out.state.istep_last == istep
        assert float(out.state.gexp_1[0, 0, 0, 0, 0, 0, 0].real) == istep
        assert float(out.state.gexp_2[0, 0, 0, 0, 0, 0, 0].real) == istep - 1
        assert float(out.state.gexp_3[0, 0, 0, 0, 0, 0, 0].real) == max(istep - 2, 0)
        assert out.time.code_dt_cfl == 1.0 / (1.0 + 2.0) == out.time.user_dt_cfl
        out2 = xp.add_explicit_terms(out.state, g, f, istep, out.time, False, None, None,
                                     _xparams(), kernel=kern)
        assert not out2.computed and out2.state is not None
        for a, b in zip(out2.state[:3], out.state[:3]):
            assert a is b
        assert out2.time == out.time
        st, t = out2.state, out2.time
    assert len(kern.calls) == 3


def test_immediate_reset_and_redo():
    """code_dt > dt_cfl: reset, early return, abort timeadv; the redone step reuses gexp."""
    kern = _stub_kernel(30.0, 20.0)                 # dt_cfl = 1/50 = 0.02 < 0.05
    st = _state0()
    g = jnp.zeros(st.gexp_1.shape, jnp.complex128)
    f = Fields(*(jnp.zeros((2, 3, 3), jnp.complex128),) * 3)
    out = xp.add_explicit_terms(st, g, f, 1, _time(), False, None, None, _xparams(), kernel=kern)
    assert out.reset and out.returned_early and out.abort_timeadv and out.computed
    assert out.time.code_dt_cfl == 1.0 / 50.0
    assert out.state.istep_last == 1
    # reset_time_step halves dt (integration layer); gs2_main sets reset = .false.; redo
    t2 = out.time._replace(code_dt=0.0125, user_dt=0.0125)
    for _ in range(2):                              # both passes of the redone step
        out2 = xp.add_explicit_terms(out.state, g, f, 1, t2, False, None, None, _xparams(),
                                     kernel=kern)
        assert not out2.computed and not out2.reset and not out2.abort_timeadv
        assert out2.state.gexp_1 is out.state.gexp_1 and out2.time == t2
    assert len(kern.calls) == 1


def test_deferred_reset_path():
    """immediate_reset = F: no CFL update inside the step; check_time_step folds it in."""
    kern = _stub_kernel(30.0, 20.0)
    st = _state0()
    g = jnp.zeros(st.gexp_1.shape, jnp.complex128)
    f = Fields(*(jnp.zeros((2, 3, 3), jnp.complex128),) * 3)
    xpar = _xparams(immediate_reset=False)
    t = _time(cfl=0.123)
    out = xp.add_explicit_terms(st, g, f, 1, t, False, None, None, xpar, kernel=kern)
    assert out.computed and not out.reset and not out.abort_timeadv and out.state.cfl_pending
    assert out.time == t                             # code_dt_cfl untouched
    assert out.state.max_vel_components == (30.0, 20.0, EPS)
    st2, t2, reset = xp.check_time_step(out.state, out.time, False, xpar)
    assert t2.code_dt_cfl == 1.0 / 50.0 and reset and not st2.cfl_pending
    # with immediate_reset the end-of-step check does not touch code_dt_cfl
    st3, t3, reset3 = xp.check_time_step(out.state, t, False, _xparams())
    assert t3 == t and not reset3                    # 0.05 < 0.123, and 0.05 >= min(0.05, ...)


def test_check_time_step_increase():
    """dt is increased when code_dt < min(code_dt_max, (dt_cfl/delt_adj)/delt_cushion)."""
    st = _state0()
    xpar = _xparams()
    t = _time(code_dt=0.0125, cfl=1.0, dt_max=0.1)
    assert xp.check_time_step(st, t, False, xpar)[2]                 # 0.0125 < min(0.1, 1/3)
    t = _time(code_dt=0.1, cfl=0.31, dt_max=0.1)
    assert not xp.check_time_step(st, t, False, xpar)[2]             # 0.1 = dt_max
    t = _time(code_dt=0.05, cfl=0.12, dt_max=0.1)
    assert xp.check_time_step(st, t, False, xpar)[2] is False        # 0.05 >= 0.12/2/1.5 = 0.04
    assert not xp.check_time_step(st, _time(code_dt=1.0, cfl=0.1), True, xpar)[2]   # exiting


def test_cfl_limit_off_and_1d():
    kern = _stub_kernel(30.0, 20.0)
    st = _state0()
    g = jnp.zeros(st.gexp_1.shape, jnp.complex128)
    f = Fields(*(jnp.zeros((2, 3, 3), jnp.complex128),) * 3)
    out = xp.add_explicit_terms(st, g, f, 1, _time(), False, None, None,
                                _xparams(use_cfl_limit=False), kernel=kern)
    assert out.time.code_dt_cfl == 1.0 / (EPS + EPS) and not out.reset          # P12
    out = xp.add_explicit_terms(st, g, f, 1, _time(), False, None, None,
                                _xparams(use_2d_cfl=False), kernel=kern)
    assert out.time.code_dt_cfl == 1.0 / 30.0
    tz = xp.set_cfl_time_step_and_check_if_too_large((0.0, 0.0, 0.0), _time(), _xparams())[0]
    assert tz.code_dt_cfl == 1e8


def test_linear_and_istep0():
    st = xp.ExplicitState(None, None, None, 0)
    out = xp.add_explicit_terms(st, None, None, 5, _time(cfl=3.0), False, None, None,
                                _xparams(nonlin=False))
    assert out.time.code_dt_cfl == 1e8 == out.time.user_dt_cfl and not out.computed
    kern = _stub_kernel(1.0, 1.0)
    out = xp.add_explicit_terms(_state0(), None, None, 0, _time(cfl=3.0), False, None, None,
                                _xparams(), kernel=kern)
    assert out.time.code_dt_cfl == 3.0 and not kern.calls
    with pytest.raises(NotImplementedError):
        xp.add_explicit_terms(_state0(), None, None, 1, _time(), False, None, None, _xparams(zip=True))


def test_ab_coefficients():
    t = _time(code_dt=0.1)
    assert xp.time_history_available(t) == 1 and xp.ab_coefficients(t) == (1.0,)
    t = _time(code_dt=0.1, prev1=0.1)
    assert xp.time_history_available(t) == 2
    np.testing.assert_allclose(xp.ab_coefficients(t), (1.5, -0.5), rtol=0, atol=2 * EPS)
    t = _time(code_dt=0.1, prev1=0.1, prev2=0.1)
    assert xp.time_history_available(t) == 3
    np.testing.assert_allclose(xp.ab_coefficients(t), (23 / 12, -4 / 3, 5 / 12), rtol=4 * EPS)
    assert len(xp.ab_coefficients(t, maximum_order=2)) == 2
    # literal formula with unequal steps (gs2_time.f90:282-291)
    dt, p1, p2 = 0.05, 0.1, 0.2
    c = xp.ab_coefficients(_time(code_dt=dt, prev1=p1, prev2=p2))
    c1 = (((dt + p1) ** 2 / p1) * ((dt + p1) / 3. + p2 / 2.) - p1 * (p1 / 3. + p2 / 2.)) / (dt * (p1 + p2))
    assert c[0] == c1
    assert abs(sum(c) - 1.0) < 1e-14                    # consistency


def test_order_based_error():
    """estimate_error / target_dt (nonlinear_terms.fpp:540-556, 620-672)."""
    rng = np.random.default_rng(11)
    shape = (1, 1, 2, 2, 3, 2, 3)
    g1, g2, g3 = (jnp.asarray(rng.standard_normal(shape) + 1j * rng.standard_normal(shape))
                  for _ in range(3))
    xpar = _xparams(use_order_based_error=True, istep_error_start=2)
    t1 = _time(code_dt=0.1, prev1=0.1)
    assert xp.estimate_error(g1, g2, g3, t1, xpar) == (0.1 * 0.5) * 0.5
    t = _time(code_dt=0.1, prev1=0.1, prev2=0.08)
    c3, c2 = xp.ab_coefficients(t, 3), xp.ab_coefficients(t, 2)
    a1, a2, a3 = (np.asarray(x) for x in (g1, g2, g3))
    s3 = c3[0] * a1 + c3[1] * a2 + c3[2] * a3
    s2 = c2[0] * a1 + c2[1] * a2
    ref = np.max(np.abs(s3 - s2)) / np.max(np.abs(s2))
    err = xp.estimate_error(g1, g2, g3, t, xpar)
    assert abs(err - ref) <= 1e-14 * ref
    # inside add_explicit: component 3 = 1/target_dt, only for istep > istep_error_start
    st = xp.ExplicitState(g1, g2, g3, istep_last=2)
    kern = lambda g, f, d, p: (g1, jnp.asarray([1.0, 1.0, EPS]))  # noqa: E731
    out = xp.add_explicit_terms(st, None, None, 3, t, False, None, None, xpar, kernel=kern)
    e = xp.estimate_error(g1, g1, g2, t, xpar)
    target = 0.1 * np.sqrt(0.1 / e)
    assert out.state.max_vel_components[2] == 1.0 / target
    assert out.time.code_dt_cfl == 1.0 / max(2.0, 1.0 / target)
    out = xp.add_explicit_terms(st, None, None, 3, t, False, None, None,
                                xpar._replace(istep_error_start=5), kernel=kern)
    assert out.computed and out.state.max_vel_components[2] == EPS


# ============================================================================================
# reference data
# ============================================================================================


def _is_nl(ex) -> bool:
    return bool(np.asarray(ex["/nonlinear_terms/nonlin"]).item())


NL_CASES = cases_where(_is_nl)


def _g7_of(grp):
    """The array of the nl_before_center / nl_after_center probes: ``g1`` (09a section 3.6)."""
    return grp["g1"]


def _vec3_of(grp):
    """The array of the nl_cfl probe: ``max_vel_components`` (09a section 3.6)."""
    return np.asarray(grp["max_vel_components"], np.float64)


def _bool(v) -> bool:
    return bool(np.asarray(v).item())


@lru_cache(maxsize=4)
def _setup(case):
    ex = export_of(case)
    data = to_jax(nlm.NLData.from_export(ex))
    return ex, Dims.from_export(ex), data, nlm.NLParams.from_export(ex), xp.ExplicitParams.from_export(ex)


def _start_of_step(tr, rec):
    """The g / old fields add_nl reads, in rec's attempt: ``after_exb_shear`` when flow shear
    is active, else ``advance_begin`` (after_exb_shear is only recorded with kx_shift /
    theta0_shift allocated, 09a section 3.6)."""
    start = tr.previous(rec, "after_exb_shear") or tr.previous(rec, "advance_begin")
    assert start is not None, f"no advance_begin before {rec.path}"
    return tr.load(start)


def _nl_inputs(case, tr, rec, dims):
    grp = _start_of_step(tr, rec)
    g = jnp.asarray(grp["gnew"])
    return g, _jfields(fields_from_group(grp, dims))


def _need_records(tr, tag):
    """Records every nonlinear trace must contain (an empty list is an exporter bug)."""
    recs = tr.find(tag)
    if not recs and not tr.find("step_end"):
        pytest.skip(f"{tr.path.name}: GS2 aborted before completing a step (no '{tag}'; the "
                    "exit path is covered by test_end_to_end)")
    assert recs, f"{tr.path}: no '{tag}' records"
    return recs


@pytest.mark.parametrize("case", NL_CASES)
def test_cflx_cfly_vs_export(case):
    """cflx/cfly recomputed (nonlinear_terms.fpp:357-358) == exported, bitwise (10 A.5 BOTH)."""
    ex, dims, data, p, xpar = _setup(case)
    nt = ex["/nonlinear_terms"]
    for name in ("cflx", "cfly"):
        assert float(np.asarray(nt[name]).item()) == getattr(p, name), name
    assert xpar.cfl == p.cfl


@pytest.mark.parametrize("case", NL_CASES)
def test_add_nl_vs_trace(case):
    """add_nl(g, phi, apar, bpar) vs probe nl_before_center (attempt with the NL evaluation)."""
    ex, dims, data, p, xpar = _setup(case)
    tr = trace_of(case)
    worst = 0.0
    for rec in _need_records(tr, "nl_before_center"):
        g, f = _nl_inputs(case, tr, rec, dims)
        ref = _g7_of(tr.load(rec))
        nl, _ = nlm.add_nl_jit(g, f, data, p)
        err = assert_close(ref, np.asarray(nl), RTOL_KERNEL, f"{case} istep={rec.istep} a={rec.attempt}")
        worst = max(worst, err)
    _report(f"{case} add_nl vs nl_before_center (max over steps)", worst)


@pytest.mark.parametrize("case", NL_CASES)
def test_add_nl_fftw_replica_vs_trace(case):
    """Diagnostic: the FFTW/NumPy replica of add_nl against GS2. Expected ~0 (same FFT
    library, plans and association order); separates FFT-library from algorithm errors."""
    _need_fftw()
    ex, dims, data, p, xpar = _setup(case)
    tr = trace_of(case)
    D = nlm.NLData.from_export(ex)
    worst = 0.0
    for rec in _need_records(tr, "nl_before_center")[:3]:
        grp = _start_of_step(tr, rec)
        f = fields_from_group(grp, dims)
        ref = _g7_of(tr.load(rec))
        nl, _ = oracle_add_nl(np.asarray(grp["gnew"]), f, D, p, forbid=D.forbid)
        worst = max(worst, max_abs_rel(ref, nl))
    _report(f"{case} FFTW replica vs nl_before_center", worst)
    assert worst <= RTOL_KERNEL


@pytest.mark.parametrize("case", NL_CASES)
def test_nl_center_vs_trace(case):
    """nl_center(nl_before_center) vs nl_after_center."""
    ex, dims, data, p, xpar = _setup(case)
    tr = trace_of(case)
    worst = 0.0
    for rec in _need_records(tr, "nl_after_center"):
        before = tr.previous(rec, "nl_before_center")
        assert before is not None
        got = nlm.nl_center_jit(jnp.asarray(_g7_of(tr.load(before))), p.bd, data, p)
        err = assert_close(_g7_of(tr.load(rec)), np.asarray(got), RTOL_KERNEL,
                           f"{case} istep={rec.istep} a={rec.attempt}")
        worst = max(worst, err)
    _report(f"{case} nl_center vs nl_after_center", worst)


@pytest.mark.parametrize("case", NL_CASES)
def test_cfl_vs_trace(case):
    """max_vel_components, code_dt_cfl and reset vs probe nl_cfl."""
    ex, dims, data, p, xpar = _setup(case)
    tr = trace_of(case)
    recs = tr.find("nl_cfl")
    if not recs:
        if not xpar.immediate_reset and (xpar.use_cfl_limit or xpar.use_order_based_error):
            pytest.skip("deferred CFL check: no nl_cfl probe (see test_end_of_step_cfl_vs_trace)")
        pytest.skip("trace has no 'nl_cfl' records")
    worst = 0.0
    for rec in recs:
        grp = tr.load(rec)
        g, f = _nl_inputs(case, tr, rec, dims)
        start = tr.load(tr.previous(rec, "advance_begin"))
        st = xp.ExplicitState.from_group(start)
        t0 = TimeBlock.from_group(start)
        out = xp.add_explicit_terms(st, g, f, rec.istep, t0, _bool(start["reset"]), data, p, xpar)
        ref = _vec3_of(grp)
        got = np.asarray(out.state.max_vel_components)
        err = float(np.max(np.abs(got - ref) / np.abs(ref)))
        worst = max(worst, err)
        assert err <= RTOL_KERNEL, f"{case} istep={rec.istep}: max_vel_components {got} vs {ref}"
        ref_cfl = float(grp["code_dt_cfl"])
        assert abs(out.time.code_dt_cfl - ref_cfl) <= RTOL_KERNEL * ref_cfl
        knife = abs(out.time.code_dt - ref_cfl) / out.time.code_dt
        if knife < 1e-10:
            print(f"[nl-err] WARNING {case} istep={rec.istep}: knife-edge CFL decision "
                  f"|code_dt-dt_cfl|/code_dt = {knife:.2e} (spec 05 P9)")
        assert out.reset == _bool(grp["reset"]), f"{case} istep={rec.istep}: reset"
    _report(f"{case} max_vel_components vs nl_cfl", worst)


def _attempt_records(tr, istep, attempt, tag):
    return [r for r in tr.find(tag, istep=istep, attempt=attempt)]


@pytest.mark.parametrize("case", NL_CASES)
def test_explicit_bookkeeping_vs_trace(case):
    """gexp_1..3 / istep_last / code_dt_cfl / reset vs after_add_explicit_terms, both passes,
    every attempt; the second pass and a redone step must not reshuffle."""
    ex, dims, data, p, xpar = _setup(case)
    tr = trace_of(case)
    _need_records(tr, "after_add_explicit_terms")
    worst = 0.0
    carried = None
    for istep in tr.steps:
        for attempt in tr.attempts(istep):
            begin = _attempt_records(tr, istep, attempt, "advance_begin")
            if not begin:
                continue
            bgrp = tr.load(begin[0])
            st = xp.ExplicitState.from_group(bgrp)
            if carried is not None and carried.istep_last == istep:
                # Redone step (immediate reset). GS2 invariant: it starts from the history of
                # the aborted attempt (no reshuffle), nl_istep_last == istep.
                assert st.istep_last == istep, f"{case} istep={istep} a={attempt}: nl_istep_last"
                prev = _attempt_records(tr, istep, attempt - 1, "after_add_explicit_terms")
                pgrp = tr.load(prev[-1])
                for k in (1, 2, 3):
                    assert np.array_equal(pgrp[f"gexp_{k}"], bgrp[f"gexp_{k}"]), \
                        f"{case} istep={istep}: GS2 history changed between attempts (gexp_{k})"
                # ours: the carried state (gexp_1 is our own NL term: round-off level)
                assert_close(bgrp["gexp_1"], np.asarray(carried.gexp_1), RTOL_KERNEL,
                             f"{case} istep={istep} carried gexp_1")
                for k in (2, 3):
                    assert np.array_equal(bgrp[f"gexp_{k}"], np.asarray(carried[k - 1])), f"gexp_{k}"
                st = carried
            time = TimeBlock.from_group(bgrp)
            reset = _bool(bgrp["reset"])
            recs = _attempt_records(tr, istep, attempt, "after_add_explicit_terms")
            for rec in recs:
                start = tr.previous(rec, "after_exb_shear") or begin[0]
                sgrp = tr.load(start)
                g = jnp.asarray(sgrp["gnew"])
                f = _jfields(fields_from_group(sgrp, dims))
                out = xp.add_explicit_terms(st, g, f, istep, time, reset, data, p, xpar)
                grp = tr.load(rec)
                what = f"{case} istep={istep} a={attempt} pass={rec.pass_}"
                worst = max(worst, assert_close(grp["gexp_1"], np.asarray(out.state.gexp_1),
                                                RTOL_KERNEL, what + " gexp_1"))
                for k in (2, 3):
                    ref = grp[f"gexp_{k}"]
                    got = np.asarray(getattr(out.state, f"gexp_{k}"))
                    if out.computed:
                        # pure copies of the previous history: bitwise
                        assert np.array_equal(ref, got), what + f" gexp_{k} (copy)"
                    else:
                        assert_close(ref, got, RTOL_KERNEL, what + f" gexp_{k}")
                assert out.time.code_dt_cfl == pytest.approx(float(grp["code_dt_cfl"]), rel=RTOL_KERNEL)
                assert out.reset == _bool(grp["reset"]), what + " reset"
                assert out.state.istep_last == int(np.asarray(grp["nl_istep_last"]).item()), \
                    what + " istep_last"
                if rec.pass_ == 2 or attempt > 1:
                    assert not out.computed, what + ": NL recomputed outside the first call"
                st, time, reset = out.state, out.time, out.reset
                if out.abort_timeadv:
                    later = [r for r in recs if r.seq > rec.seq]
                    assert not later, what + ": timeadv continued after an immediate reset"
                    break
            carried = st
            # attempt 2 of an immediately reset step has no nl_* probes
            if attempt > 1:
                for tag in ("nl_before_center", "nl_after_center", "nl_cfl"):
                    assert not tr.find(tag, istep=istep, attempt=attempt), f"{case}: {tag} in retry"
    _report(f"{case} gexp_1 vs after_add_explicit_terms", worst)


def _end_of_step_records(tr, istep):
    """(step_end record, after_reset_time_step record after it or None) of the last attempt."""
    attempt = tr.last_attempt(istep)
    ends = tr.find("step_end", istep=istep, attempt=attempt)
    if not ends:
        return None, None
    end = ends[-1]
    resets = [r for r in tr.find("after_reset_time_step", istep=istep) if r.key > end.key]
    return end, (resets[0] if resets else None)


@pytest.mark.parametrize("case", NL_CASES)
def test_end_of_step_cfl_vs_trace(case):
    """check_time_step at the end of every step, driven by OUR CFL value of that step:
    immediate_reset = T -> the code_dt_cfl saved by our add_explicit; = F -> the deferred
    fold-in of our max_vel_components (nb_check). Compared with GS2's code_dt_cfl after the
    check (in after_reset_time_step, or first in advance_begin(k+1), 09a section 3.5) and with
    GS2's dt-change decision (an after_reset_time_step record after step_end)."""
    ex, dims, data, p, xpar = _setup(case)
    tr = trace_of(case)
    _need_records(tr, "step_end")
    nchecked, worst = 0, 0.0
    for istep in tr.steps:
        end, reset_rec = _end_of_step_records(tr, istep)
        if end is None:
            continue
        # this step's NL evaluation: its first attempt
        bgrp = tr.load(tr.find("advance_begin", istep=istep)[0])
        st = xp.ExplicitState.from_group(bgrp)
        t0 = TimeBlock.from_group(bgrp)
        first_call = tr.find("after_add_explicit_terms", istep=istep)[0]
        sgrp = _start_of_step(tr, first_call)
        out = xp.add_explicit_terms(st, jnp.asarray(sgrp["gnew"]),
                                    _jfields(fields_from_group(sgrp, dims)), istep, t0,
                                    _bool(bgrp["reset"]), data, p, xpar)
        assert out.computed
        t_end = TimeBlock.from_group(tr.load(end))
        if xpar.immediate_reset:
            # nothing between add_explicit and check_time_step touches code_dt_cfl
            assert float(t_end.code_dt_cfl) == pytest.approx(out.time.code_dt_cfl, rel=RTOL_KERNEL)
            t_end = t_end._replace(code_dt_cfl=out.time.code_dt_cfl, user_dt_cfl=out.time.user_dt_cfl)
        st2, t2, reset = xp.check_time_step(out.state, t_end, False, xpar)
        if reset_rec is not None:
            ref_cfl = float(tr.load(reset_rec)["code_dt_cfl"])
        else:
            nxt = tr.find("advance_begin", istep=istep + 1)
            if not nxt:
                continue                  # last traced step: GS2's value is not recorded
            ref_cfl = float(tr.load(nxt[0])["code_dt_cfl"])
        err = abs(t2.code_dt_cfl - ref_cfl) / ref_cfl
        worst = max(worst, err)
        assert err <= RTOL_KERNEL, f"{case} istep={istep}: code_dt_cfl {t2.code_dt_cfl} vs {ref_cfl}"
        assert reset == (reset_rec is not None), f"{case} istep={istep}: end-of-step reset"
        nchecked += 1
    assert nchecked, "no complete steps in the trace"
    _report(f"{case} end-of-step code_dt_cfl ({nchecked} steps, decisions equal)", worst)


@pytest.mark.parametrize("case", NL_CASES)
def test_dt_change_cases_exercise_resets(case):
    """The dt-change cases must change dt within the traced steps (gpu/cases/README.md):
    immediate -> a retried step (attempt 2); deferred -> end-of-step resets and never a retry
    (impossible with immediate_reset = F); increase -> an end-of-step reset that raises dt."""
    if not case.startswith("nl_box_ae_dt_"):
        pytest.skip("not a dt-change case")
    tr = trace_of(case)
    retries = [s for s in tr.steps if len(tr.attempts(s)) > 1]
    end_resets = [s for s in tr.steps if _end_of_step_records(tr, s)[1] is not None]
    if case.endswith("immediate"):
        assert retries, "no immediate reset in the traced steps"
    elif case.endswith("deferred"):
        assert end_resets and not retries
    else:
        assert end_resets, "no end-of-step reset in the traced steps"
        s = end_resets[0]
        dt_before = float(tr.load(_end_of_step_records(tr, s)[0])["code_dt"])
        dt_after = float(tr.load(_end_of_step_records(tr, s)[1])["code_dt"])
        assert dt_after > dt_before


@pytest.mark.parametrize("case", NL_CASES)
def test_ab_order_buildup_vs_trace(case):
    """Fresh start: AB order 1, 2, 3 on steps 1, 2, 3 (spec 02 section 3.7) -- the history
    length seen by every attempt of step k is min(k, 3) and a reset does not change it."""
    tr = trace_of(case)
    recs = _need_records(tr, "advance_begin")
    t0 = TimeBlock.from_group(tr.load(recs[0]))
    if recs[0].istep != 1 or t0.code_dt_prev1 != xp.DT_NOT_SET:
        pytest.skip("trace does not start at a fresh step 1")
    for rec in recs:
        t = TimeBlock.from_group(tr.load(rec))
        assert xp.time_history_available(t) == min(rec.istep, 3), f"istep={rec.istep} a={rec.attempt}"
        if rec.istep <= 3:
            assert len(xp.ab_coefficients(t)) == min(rec.istep, 3)
