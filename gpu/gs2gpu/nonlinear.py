"""The nonlinear E x B term of GS2 (spec 05): loads, GS2 FFTs, Poisson bracket, theta centring.

Mirrors (line numbers of GS2 8.2.1, commit 4d8c94bcf):

* ``nonlinear_terms::add_nl``               src/nonlinear_terms.fpp:788-929
* ``load_k{x,y}_{phi,apar,bpar}``          src/nonlinear_terms.fpp:932-1134
* ``calculate_bracket``                     src/nonlinear_terms.fpp:1141-1187
* ``nl_center``                             src/nonlinear_terms.fpp:677-760
* ``gs2_transforms::transform2_5d``         src/gs2_transforms.fpp:919-955 (x c2c BACKWARD, then y c2r)
* ``gs2_transforms::inverse2_5d``           src/gs2_transforms.fpp:958-985 (y r2c, then x c2c FORWARD)
* ``transform2_5d_accel/inverse2_5d_accel`` src/gs2_transforms.fpp:990-1352 (2D FFTW plans;
  identical semantics, scale ``1/real(nx*ny)``, externals/utils/fft_work.fpp:220,481)

Layouts. Spectral arrays are canonical G7 ``(S, E, L, Y, X, 2, T)`` (Fortran
``g(ig, isgn, it, ik, il, ie, is)``); fields are F3 ``(Y, X, T)``. Real-space arrays use the
"NL layout" ``(S, E, L, 2, T, nx, ny)`` (spec 05 section 3.3), i.e. NumPy index
``[is, ie, il, isgn, ig, ix, iy]`` of GS2's ``ba(iy, iyxf)``.

GS2 Fourier convention (gs2_transforms.fpp:929-942): ``F_k = f_k`` for ky = 0 and ``2 f_k`` for
ky > 0, i.e. ``f(x, y) = sum_kx F[0, kx] e^{i kx x} + sum_{kx, ky>0} 2 Re((F/2) e^{i(kx x + ky y)})``.
kx index ``it`` sits at FFT position ``J(it) = it-1`` (kx >= 0) or ``it-1+nx-ntheta0`` (kx < 0)
(gs2_transforms.fpp:485-489); ky index ``ik`` sits at ``ik-1``; everything else is zero padding
(the 2/3 rule). FFTW's c2r ignores the imaginary part of the ky = 0 (and Nyquist) bin; this is
reproduced explicitly (spec 05 section 4.5.1, P6).

Round-off. The arithmetic keeps the Fortran association order; products that feed an addition
go through ``numerics.ob`` / ``numerics.rmul`` so that ``numerics.set_strict(True)`` removes
XLA's FMA contraction (the reference build uses ``-ffp-contract=off``, 00_decisions D2). The
FFTs are XLA's (DUCC on CPU) instead of FFTW's: same transforms, different butterflies, about
1 ulp of the largest element per transform (measured in tests/test_nonlinear.py).
"""

from __future__ import annotations

from typing import NamedTuple

import jax
import jax.numpy as jnp
import numpy as np
from jax import lax

from .gadjust import GAdjustData, g_adjust
from .io.export import Export
from .numerics import ob, rmul
from .types import Fields

#: ``epsilon(0.0)`` under ``-fdefault-real-8`` (spec 05 build facts).
EPS = float(np.finfo(np.float64).eps)


def _item(value):
    return np.asarray(value).item()


def _flag(value) -> bool:
    return bool(np.asarray(value).item())


# --------------------------------------------------------------------------------------------
# static parameters and arrays
# --------------------------------------------------------------------------------------------


class NLParams(NamedTuple):
    """Static switches and scalars of the nonlinear term (hashable; a ``jit`` static argument).

    ``phi_on`` = ``has_phi .and. include_phi`` (load_k*_phi, nonlinear_terms.fpp:946/979),
    ``apar_on`` = ``include_apar .and. has_apar`` (:1015/:1048),
    ``bpar_on`` = ``include_bpar .and. has_bpar`` (:1081/:1116).
    ``cflx = maxval(akx)*0.5/cfl``, ``cfly = maxval(aky)*0.5/cfl`` (:357-358).
    ``bd = bkdiff(1)`` is the argument ``nl_center`` receives (dist_fn.fpp:4061).
    """

    nonlin: bool
    nx: int
    ny: int
    naky: int
    ntheta0: int
    accelerated: bool
    kxfac: float
    cfl: float
    cflx: float
    cfly: float
    fphi: float
    fapar: float
    fbpar: float
    phi_on: bool
    apar_on: bool
    bpar_on: bool
    nl_forbid_force_zero: bool
    trapped: bool
    bd: float

    @classmethod
    def from_export(cls, ex: Export) -> "NLParams":
        """Static NL parameters (export groups /nonlinear_terms, /kt_grids, /run_parameters,
        /theta_grid, /le_grids, /dist_fn). Raises ``NotImplementedError`` for the GATED
        options ``zip`` and ``split_nonlinear`` (10_critique G-08, B.4)."""
        nt = ex["/nonlinear_terms"]
        kt = ex["/kt_grids"]
        rp = ex["/run_parameters"]
        nonlin = _flag(nt["nonlin"])
        if nonlin and _flag(nt["zip"]):
            raise NotImplementedError("nonlinear_terms_knobs: zip = .true. is GATED (10_critique G-08)")
        if nonlin and _flag(nt["split_nonlinear"]):
            raise NotImplementedError("nonlinear_terms_knobs: split_nonlinear = .true. is GATED "
                                      "(spec 05 section 4.13)")
        akx = np.asarray(kt["akx"], dtype=np.float64)
        aky = np.asarray(kt["aky"], dtype=np.float64)
        nx = int(_item(kt["nx"]))
        ny = int(_item(kt["ny"]))
        naky, ntheta0 = ex.naky, ex.ntheta0
        cfl = float(_item(nt["cfl"]))
        if nonlin:
            if nx == 0 or ny == 0:
                # gs2_transforms.fpp:106-108 aborts here.
                raise ValueError("nonlinear run with nx or ny = 0 (grid_option must be 'box')")
            if ntheta0 > nx or naky > ny // 2 + 1:
                raise ValueError(f"FFT grid too small: ntheta0={ntheta0} nx={nx} naky={naky} ny={ny}")
            # nonlinear_terms.fpp:357-358, evaluated left to right.
            cfly = (float(np.max(aky)) * 0.5) / cfl
            cflx = (float(np.max(akx)) * 0.5) / cfl
        else:
            cflx = cfly = 0.0
        has_phi = _flag(rp["has_phi"])
        has_apar = _flag(rp["has_apar"])
        has_bpar = _flag(rp["has_bpar"])
        trapped = _flag(ex["/le_grids/grid_has_trapped_particles"])     # le_grids.f90:443-446
        bkdiff = np.atleast_1d(np.asarray(ex["/dist_fn/bkdiff"], dtype=np.float64))
        return cls(
            nonlin=nonlin, nx=nx, ny=ny, naky=naky, ntheta0=ntheta0,
            accelerated=_flag(nt["accelerated"]),
            kxfac=float(_item(ex["/theta_grid/kxfac"])),
            cfl=cfl, cflx=cflx, cfly=cfly,
            fphi=float(_item(rp["fphi"])), fapar=float(_item(rp["fapar"])),
            fbpar=float(_item(rp["fbpar"])),
            phi_on=has_phi and _flag(nt["include_phi"]),
            apar_on=has_apar and _flag(nt["include_apar"]),
            bpar_on=has_bpar and _flag(nt["include_bpar"]),
            nl_forbid_force_zero=_flag(nt["nl_forbid_force_zero"]),
            trapped=trapped,
            bd=float(bkdiff[0]),
        )


class NLData(NamedTuple):
    """Arrays read by the nonlinear term (dt independent).

    ``akx (X,)``, ``aky (Y,)``; ``vpa`` G7 real; ``stm, tz (S,)``; ``forbid, is_ttp (L, T)`` bool;
    ``cell_zero (L, T-1)`` bool = ``trapped .and. (il > jend(ig) .or. il > jend(ig+1))``
    (nl_center, nonlinear_terms.fpp:741); ``gadj`` holds ``vperp2, aj0, aj1`` (G6) and the
    g_adjust inputs (shared, not duplicated).
    """

    akx: np.ndarray
    aky: np.ndarray
    vpa: np.ndarray
    stm: np.ndarray
    tz: np.ndarray
    forbid: np.ndarray
    is_ttp: np.ndarray
    cell_zero: np.ndarray
    gadj: GAdjustData

    @classmethod
    def from_export(cls, ex: Export) -> "NLData":
        kt = ex["/kt_grids"]
        le = ex["/le_grids"]
        sp = ex["/species"]                     # stored as /species_ (09a section 3.1); loader alias
        trapped = _flag(le["grid_has_trapped_particles"])
        jend = np.asarray(le["jend"]).astype(np.int64)               # (T,), Fortran values
        il = np.arange(1, ex.nlambda + 1)[:, None]                    # Fortran il, (L, 1)
        cell_zero = trapped & ((il > jend[None, :-1]) | (il > jend[None, 1:]))
        return cls(
            akx=np.asarray(kt["akx"], dtype=np.float64),
            aky=np.asarray(kt["aky"], dtype=np.float64),
            vpa=np.asarray(ex["/dist_fn_arrays/vpa"], dtype=np.float64),
            stm=np.asarray(sp["spec_stm"], dtype=np.float64),
            tz=np.asarray(sp["spec_tz"], dtype=np.float64),
            forbid=np.asarray(le["forbid"], dtype=bool),
            is_ttp=np.asarray(le["is_ttp"], dtype=bool),
            cell_zero=np.asarray(cell_zero, dtype=bool),
            gadj=GAdjustData.from_export(ex),
        )


# --------------------------------------------------------------------------------------------
# GS2 transforms (spec 05 section 5)
# --------------------------------------------------------------------------------------------


def kx_fft_index(ntheta0: int, nx: int) -> np.ndarray:
    """0-based FFT position ``J`` of each 0-based kx index (gs2_transforms.fpp:485-489:
    ``it0 = it + nx - ntheta0`` if ``it > (ntheta0+1)/2`` else ``it``)."""
    it = np.arange(1, ntheta0 + 1)
    return np.where(it > (ntheta0 + 1) // 2, it + nx - ntheta0, it) - 1


def _to_xy(a):
    """G7 ``(S,E,L,Y,X,2,T)`` -> ``(S,E,L,2,T,X,Y)`` (FFT axes last)."""
    return jnp.transpose(a, (0, 1, 2, 5, 6, 4, 3))


def _from_xy(a):
    """``(S,E,L,2,T,X,Y)`` -> G7 ``(S,E,L,Y,X,2,T)``."""
    return jnp.transpose(a, (0, 1, 2, 6, 5, 3, 4))


def _c2c_backward(a, axis):
    """Unnormalised c2c with ``exp(+2 pi i j k / n)`` (FFTW_BACKWARD = FFT_TO_REAL_SPACE).

    Written as ``conj(fft(conj(a)))``: conjugation is exact, and XLA's inverse FFT would apply
    a ``1/n`` that ``norm="forward"`` multiplies back (two roundings unless n is a power of 2).
    """
    return jnp.conj(jnp.fft.fft(jnp.conj(a), axis=axis))


def transform2_xy(gxy, params: NLParams):
    """``transform2_5d`` on an array with the spectral axes last: ``(..., X, Y)`` -> ``(..., nx, ny)``.

    Mirrors gs2_transforms.fpp:919-955 + transform_x5d (830-855) + transform_y5d (875-899):
    halve ky > 0 (exact; GS2 does it in place, P7 -- the input here is not modified), gather
    into the zero-padded ``xxf`` (J mapping), c2c BACKWARD along x (length nx), gather into
    ``fft`` (rows ky >= naky zero), c2r along y (length ny). Forbidden lines are transformed
    too (their loads are exact zeros, spec 05 section 4.5.3).
    """
    nx, ny = params.nx, params.ny
    ntheta0, naky = gxy.shape[-2], gxy.shape[-1]
    half = jnp.asarray(np.where(np.arange(naky) == 0, 1.0, 0.5))
    gh = rmul(gxy, half)                                           # :947-951 (exact)
    npos = (ntheta0 + 1) // 2                                     # kx >= 0 count
    lead = gh.shape[:-2]
    pad_x = jnp.zeros(lead + (nx - ntheta0, naky), gh.dtype)
    xxf = jnp.concatenate([gh[..., :npos, :], pad_x, gh[..., npos:, :]], axis=-2)
    xk = _c2c_backward(xxf, axis=-2)                               # xf_fft (:853)
    nyc = ny // 2 + 1
    fft = jnp.concatenate([xk, jnp.zeros(lead + (nx, nyc - naky), xk.dtype)], axis=-1)
    # FFTW c2r ignores Im of the DC bin (and of the Nyquist bin for even ny) -- P6.
    fft = fft.at[..., 0].set(lax.complex(fft[..., 0].real, jnp.zeros_like(fft[..., 0].real)))
    if ny % 2 == 0 and naky > ny // 2:
        k = ny // 2
        fft = fft.at[..., k].set(lax.complex(fft[..., k].real, jnp.zeros_like(fft[..., k].real)))
    return jnp.fft.irfft(fft, n=ny, axis=-1, norm="forward")      # yf_fft (:897), unnormalised


def inverse2_xy(f, naky: int, ntheta0: int, params: NLParams):
    """``inverse2_5d`` with the spectral axes last: ``(..., nx, ny)`` real -> ``(..., X, Y)``.

    Mirrors gs2_transforms.fpp:958-985 + inverse_y5d (902-916) + inverse_x5d (858-872):
    r2c along y (unnormalised, ``exp(-i)``), keep ky < naky, c2c FORWARD along x, keep the kx of
    the J mapping, then ``g*(scale*2.0)`` for ky > 0 and ``g*scale`` for ky = 0 with
    ``scale = (1/nx)*(1/ny)`` (:968) or, accelerated, ``1/real(nx*ny)`` (fft_work.fpp:220,481;
    gs2_transforms.fpp:1326-1335).
    """
    nx, ny = params.nx, params.ny
    yk = jnp.fft.rfft(f, axis=-1)[..., :naky]                     # yb_fft (:912)
    zk = jnp.fft.fft(yk, axis=-2)                                  # xb_fft (:868)
    npos = (ntheta0 + 1) // 2
    zk = jnp.concatenate([zk[..., :npos, :], zk[..., nx - ntheta0 + npos:, :]], axis=-2)
    if params.accelerated:
        scale = 1.0 / float(nx * ny)
    else:
        scale = (1.0 / float(nx)) * (1.0 / float(ny))
    fac = jnp.asarray(np.where(np.arange(naky) == 0, scale, scale * 2.0))
    return rmul(zk, fac)                                           # :977-983


def transform2(g, params: NLParams):
    """G7 spectral -> real space ``(S, E, L, 2, T, nx, ny)`` (transform2_5d, gs2_transforms.fpp:919-955)."""
    return transform2_xy(_to_xy(g), params)


def inverse2(f, params: NLParams):
    """Real space ``(S, E, L, 2, T, nx, ny)`` -> G7 (inverse2_5d, gs2_transforms.fpp:958-985)."""
    return _from_xy(inverse2_xy(f, params.naky, params.ntheta0, params))


# --------------------------------------------------------------------------------------------
# the loads of i*k*chi (nonlinear_terms.fpp:932-1134)
# --------------------------------------------------------------------------------------------


def _load_chi(fields: Fields, k6, data: NLData, params: NLParams):
    """``g1 = i*k*chi`` on G7, ``k6`` = akx or aky broadcast to G6.

    ``zi*k*x*...*z`` evaluates as ``(0, c)*(z_r, z_i) = (0*z_r - c*z_i, 0*z_i + c*z_r)``, whose
    ``0*`` terms are exact (signed) zeros, so it is written ``(-(c*z_i), c*z_r)``.
    Order phi -> bpar -> apar (nonlinear_terms.fpp:816-818, 863-865):

    * phi  (:946-957): ``fac = zi*k*aj0*phi*fphi`` = ``(-(c*phi_i), c*phi_r)*fphi``, ``c = k*aj0``;
      ``g1(:,1) = g1(:,2) = fac``; ``zero_array(g1)`` if ``.not.(has_phi .and. include_phi)``.
    * bpar (:1081-1097): ``fac = zi*k*aj1*2.0*vperp2*tz*bpar*fbpar``, ``c = (((k*aj1)*2.0)*vperp2)*tz``;
      ``g1 = g1 + fac`` on both signs.
    * apar (:1015-1029): ``fac = zi*k*aj0*stm*apar*fapar``, ``c = (k*aj0)*stm``;
      ``g1(:,isgn) = g1(:,isgn) - fac*vpa(:,isgn)``.
    """
    gadj = data.gadj
    shape6 = gadj.aj0.shape
    if params.phi_on:
        c = k6 * gadj.aj0
        p = fields.phi[None, None, None]
        g6 = lax.complex(ob(-(c * p.imag) * params.fphi), ob((c * p.real) * params.fphi))
    else:
        g6 = jnp.zeros(shape6, jnp.complex128)
    if params.bpar_on:
        tz = data.tz[:, None, None, None, None, None]
        c = (((k6 * gadj.aj1) * 2.0) * gadj.vperp2) * tz
        b = fields.bpar[None, None, None]
        fac = lax.complex(ob(-(c * b.imag) * params.fbpar), ob((c * b.real) * params.fbpar))
        g6 = g6 + fac
    g7 = jnp.broadcast_to(g6[..., None, :], shape6[:-1] + (2, shape6[-1]))
    if params.apar_on:
        stm = data.stm[:, None, None, None, None, None]
        c = (k6 * gadj.aj0) * stm
        a = fields.apar[None, None, None]
        fac = lax.complex(-(c * a.imag) * params.fapar, (c * a.real) * params.fapar)
        g7 = g7 - rmul(fac[..., None, :], data.vpa)
    return g7


def load_kx_chi(fields: Fields, data: NLData, params: NLParams):
    """``i*kx*chi`` (load_kx_phi/bpar/apar, nonlinear_terms.fpp:932-962, 1067-1099, 1001-1031)."""
    return _load_chi(fields, data.akx[None, None, None, None, :, None], data, params)


def load_ky_chi(fields: Fields, data: NLData, params: NLParams):
    """``i*ky*chi`` (load_ky_phi/bpar/apar, nonlinear_terms.fpp:965-996, 1102-1134, 1034-1064)."""
    return _load_chi(fields, data.aky[None, None, None, :, None, None], data, params)


def _times_i(h, k7):
    """Fortran ``h*zi*k`` = ``(h*zi)*k`` = ``(-h_i*k, h_r*k)`` (nonlinear_terms.fpp:845, 888)."""
    return lax.complex(-h.imag * k7, h.real * k7)


# --------------------------------------------------------------------------------------------
# add_nl and nl_center
# --------------------------------------------------------------------------------------------


def add_nl(g, fields: Fields, data: NLData, params: NLParams, adjust: bool = True):
    """The nonlinear term at theta grid points; mirrors ``add_nl`` (nonlinear_terms.fpp:788-929).

    ``g`` is the start-of-step ``g`` (G7; after exb_shear, fields_implicit.f90:337), ``fields``
    the start-of-step ``phi, apar, bpar`` (F3). Returns ``(nl, max_vel_components)``:

    * ``nl`` (G7) = ``kxfac * inverse2( d_x chi * d_y h - d_y chi * d_x h )`` where
      ``h = g_adjust(g, phi, bpar, from_g_gs2)`` (:835-837, dist_fn_arrays.fpp:101-175, skipped
      if ``adjust`` is False) and the derivatives are ``i k`` products (:816-818, 845, 863-865,
      888);
    * ``max_vel_components = [max|d_y chi|*cflx, max|d_x chi|*cfly, EPS]`` (:825-830, 872-877,
      913-914; the third entry is the ``epsilon(0.0)`` stored at :518 before the call).

    The halving GS2 applies in place to its ``transform2`` inputs (P7) is not visible here.
    """
    akx7 = data.akx[None, None, None, None, :, None, None]
    aky7 = data.aky[None, None, None, :, None, None, None]
    # (A)-(B): d_x chi in real space (:816-830)
    ba = transform2(load_kx_chi(fields, data, params), params)
    max_vel_y = jnp.max(jnp.abs(ba)) * params.cfly
    # (C)-(D): h (g_wesson) and d_y h (:835-854)
    if adjust:
        h = g_adjust(g, fields.phi, fields.bpar, params.fphi, params.fbpar, data.gadj)
    else:
        h = g
    gb = transform2(_times_i(h, aky7), params)
    # (E): bracket = ba*gb (:858, 1141-1165); stored before (J), so never part of an FMA
    bracket = ob(ba * gb)
    # (F)-(G): d_y chi (:863-877)
    ba = transform2(load_ky_chi(fields, data, params), params)
    max_vel_x = jnp.max(jnp.abs(ba)) * params.cflx
    # (H)-(I): d_x h (:885-896)
    gb = transform2(_times_i(h, akx7), params)
    # (J): bracket = bracket - ba*gb (:900, 1166-1185)
    bracket = bracket - ob(ba * gb)
    # (K)-(L): back to spectral space and scale by kxfac (:906-923)
    nl = rmul(inverse2(bracket, params), params.kxfac)
    mvc = jnp.stack([max_vel_x, max_vel_y, jnp.asarray(EPS, jnp.float64)])
    return nl, mvc


def nl_center(g1, bd, data: NLData, params: NLParams):
    """Grid values -> (2x) upwinded cell centres; mirrors ``nl_center`` (nonlinear_terms.fpp:677-760).

    * if ``nl_forbid_force_zero``: ``g1 = 0`` where ``forbid(ig, il)`` (all ig, :718-724);
    * for ``ig = -ntgrid .. ntgrid-1`` (the ascending in-place loop only ever reads the
      unmodified ``ig+1``): keep the grid value at ttp points (:739); ``0`` if
      ``trapped .and. (il > jend(ig) .or. il > jend(ig+1))`` (:741-747); else
      ``out(ig,1) = (1+bd)*g(ig+1,1) + (1-bd)*g(ig,1)``,
      ``out(ig,2) = (1-bd)*g(ig+1,2) + (1+bd)*g(ig,2)`` (:753-754);
    * ``ig = ntgrid`` keeps its (possibly zeroed) grid value.
    """
    if params.nl_forbid_force_zero:
        forbid = jnp.asarray(data.forbid)[None, None, :, None, None, None, :]
        g1 = jnp.where(forbid, jnp.zeros((), g1.dtype), g1)
    p = 1.0 + bd
    m = 1.0 - bd
    lo = g1[..., :-1]
    hi = g1[..., 1:]
    c1 = rmul(hi[..., 0:1, :], p) + rmul(lo[..., 0:1, :], m)
    c2 = rmul(hi[..., 1:2, :], m) + rmul(lo[..., 1:2, :], p)
    cells = jnp.concatenate([c1, c2], axis=-2)
    ttp = jnp.asarray(data.is_ttp)[None, None, :, None, None, None, :-1]
    zero = jnp.asarray(data.cell_zero)[None, None, :, None, None, None, :]
    cells = jnp.where(zero, jnp.zeros((), g1.dtype), cells)
    cells = jnp.where(ttp, lo, cells)
    return jnp.concatenate([cells, g1[..., -1:]], axis=-1)


def nl_source(g, fields: Fields, data: NLData, params: NLParams, bd=None):
    """``add_nl`` followed by ``nl_center`` (add_explicit, nonlinear_terms.fpp:520-523).

    Returns ``(gexp_1, max_vel_components)``; ``bd`` defaults to ``params.bd`` (``bkdiff(1)``).
    """
    nl, mvc = add_nl(g, fields, data, params)
    return nl_center(nl, params.bd if bd is None else bd, data, params), mvc


add_nl_jit = jax.jit(add_nl, static_argnames=("params", "adjust"))
nl_center_jit = jax.jit(nl_center, static_argnames=("params",))
nl_source_jit = jax.jit(nl_source, static_argnames=("params",))
transform2_jit = jax.jit(transform2, static_argnames=("params",))
inverse2_jit = jax.jit(inverse2, static_argnames=("params",))
