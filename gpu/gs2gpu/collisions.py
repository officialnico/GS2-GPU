"""Collision operator time step: ``vspace_derivatives`` and ``collisions_advance`` (spec 06 section 4).

Default ``use_le_layout = .true.`` path (le_derivatives.fpp:67-112; spec 10 G-24 gates the other).
One application, for every spatial point ``(is, ik, it, ig)`` independently (collisions.fpp:3035-3090)::

    h    = g_adjust(gnew, phinew, bparnew, from_g)          (if adjust)
    gle  = gather(g2le, h)                                  LE6 (S, Y, X, T, E+1, nxi+1), zero-padded
    gle  = ediffuse solve -> its conservation               (has_diffuse)
    gle  = drag (reads aparnew) -> Lorentz solve -> its conservation   (has_lorentz)
    h    = scatter(g2le, gle, h)                            (allowed points only)
    gnew = g_adjust(h, phinew, bparnew, to_g)               (if adjust)

The tridiagonal solves are GS2's own forward elimination / back substitution, batched over all
independent lines with ``lax.scan`` along the solved axis (energy for the energy diffusion, xi for the
pitch-angle scattering). The skip masks (``vnew`` = 0 species, ``kwork_filter``, ``ieqzip``; spec 10
G-01) are explicit ``jnp.where``. Every real x complex product is componentwise (``numerics.rmul``)
and complex / real is componentwise (gfortran's lowering, spec 06 section 8.6). The conservation
moments are summed in the Fortran order (``CollParams.exact_sums``). Heating (``c_rate``) is not
ported; it does not change the evolution (spec 10 G-26).
"""

from __future__ import annotations

from functools import partial
from typing import Callable, Optional

import jax
import jax.numpy as jnp
from jax import lax

from .collisions_init import (  # noqa: F401  (re-exported: the module's public interface)
    CollCoeffs, CollData, CollInitData, CollParams, compute_coll_coeffs, raw_from_export)
from .gadjust import FROM_G, TO_G, g_adjust_direction
from .numerics import ob, rmul
from .types import Fields

STAGE_TAGS = ("coll_h_in", "coll_gle_in", "coll_after_ediffuse", "coll_after_ediffuse_conserve",
              "coll_after_drag", "coll_after_lorentz", "coll_after_lorentz_conserve",
              "coll_gle_out", "coll_h_out")


# --------------------------------------------------------------------------- helpers
def _pad_xi(a, nx: int):
    """Pad the xi axis of a coefficient array from ``nxi_lim`` to ``nxi + 1`` rows with zeros
    (they differ only on grids without trapped particles)."""
    a = jnp.asarray(a)
    if a.shape[-1] == nx:
        return a
    return jnp.pad(a, [(0, 0)] * (a.ndim - 1) + [(0, nx - a.shape[-1])])


def _cdiv_r(z, r):
    """Fortran ``complex / real``: componentwise (spec 06 section 8.6)."""
    return lax.complex(z.real / r, z.imag / r)


def _bt(mask_t_nx):
    """(T, NX) mask -> broadcastable against LE6 (S, Y, X, T, E, NX)."""
    return jnp.asarray(mask_t_nx)[None, None, None, :, None, :]


def _take_xi(g, idx_t):
    """``g[..., idx[t]]`` for LE6 ``g`` and a per-theta row index ``idx (T,)``: (S, Y, X, T, E, 1)."""
    idx = jnp.broadcast_to(jnp.asarray(idx_t)[None, None, None, :, None, None],
                           g.shape[:-1] + (1,))
    return jnp.take_along_axis(g, idx, axis=-1, mode="clip")


def _sey(a):
    """(S, E, Y) velocity factor -> broadcastable against LE6 (S, Y, X, T, E, XI)."""
    return jnp.transpose(jnp.asarray(a), (0, 2, 1))[:, :, None, None, :, None]


def _tridiag(gs, qs, bs, cs):
    """GS2's tridiagonal sweeps along axis 0 (collisions.fpp:3958-3967 and 4063-4072)::

        delta(1) = g(1);  delta(i+1) = g(i+1) - q(i+1)*delta(i)
        x(n+1) = 0;       x(i) = (delta(i) - c(i)*x(i+1)) * beta(i),   i = n, ..., 1

    ``gs`` complex ``(n, ...)``; ``qs, bs, cs`` real, same shape. One ``lax.scan`` per sweep, batched
    over every other axis.
    """
    d0 = gs[0]

    def fwd(dprev, xs):
        gi, qi = xs
        di = gi - rmul(dprev, qi)
        return di, di

    _, drest = lax.scan(fwd, d0, (gs[1:], qs[1:]))
    delta = jnp.concatenate([d0[None], drest], axis=0)

    def bwd(gnext, xs):
        di, ci, bi = xs
        gi = rmul(di - rmul(gnext, ci), bi)
        return gi, gi

    _, out = lax.scan(bwd, jnp.zeros_like(d0), (delta, cs, bs), reverse=True)
    return out


def _moment(term, params: CollParams):
    """Velocity sum ``v = v + term(ixi, ie)`` of the conservation passes, (S, Y, X, T, E, XI) -> (S, Y, X, T).

    Fortran order (collisions.fpp:3288, 3317, 3370, 3570, 3596, 3622):
    ``v = 0``, ``do ie; do ixi``: one rounded complex add per term. With ``exact_sums = False`` a tree
    reduction (round-off level differences, spec 06 section 7).
    """
    if not params.exact_sums:
        return jnp.sum(term, axis=(-2, -1))
    flat = term.reshape(term.shape[:-2] + (-1,))          # ie outer, ixi inner (C order)
    flat = jnp.moveaxis(flat, -1, 0)

    def body(acc, x):
        return acc + x, None

    acc, _ = lax.scan(body, jnp.zeros(flat.shape[1:], flat.dtype), flat, unroll=8)
    return acc


def _set_energy_rows(gle, rows):
    """``gle`` with its first ``negrid`` energy rows replaced (row negrid+1 kept)."""
    E = rows.shape[4]
    return jnp.concatenate([rows, gle[..., E:, :]], axis=4)


# --------------------------------------------------------------------------- gather / scatter
def g2le(h, data: CollData, params: CollParams):
    """``gle = 0; gather(g2le, h, gle)`` (le_derivatives.fpp:79-83; le_grids.f90:3379-3599).

    G7 ``h (S, E, L, Y, X, 2, T)`` -> LE6 ``(S, Y, X, T, E+1, nxi+1)``; only allowed points are copied,
    everything else (forbidden rows, energy row negrid+1) is 0.
    """
    S, E, L, Y, X, _, T = h.shape
    NX = params.nx
    ht = jnp.transpose(h, (0, 3, 4, 6, 1, 2, 5)).reshape(S, Y, X, T, E, 2 * L)
    idx = jnp.broadcast_to(jnp.asarray(data.le_cidx)[None, None, None, :, None, :],
                           (S, Y, X, T, E, NX))
    core = jnp.take_along_axis(ht, idx, axis=-1, mode="clip")
    core = jnp.where(_bt(data.le_valid), core, jnp.zeros((), core.dtype))
    return jnp.pad(core, ((0, 0),) * 4 + ((0, 1), (0, 0)))


def le2g(gle, h, data: CollData, params: CollParams):
    """``scatter(g2le, gle, h)`` (le_derivatives.fpp:94): allowed points of ``h`` take the LE values,
    forbidden points keep theirs."""
    S, E, L, Y, X, _, T = h.shape
    lt = gle[..., :E, :]
    idx = jnp.transpose(jnp.asarray(data.g_ixi), (2, 0, 1)).reshape(T, 2 * L)
    idx = jnp.broadcast_to(idx[None, None, None, :, None, :], (S, Y, X, T, E, 2 * L))
    hn = jnp.take_along_axis(lt, idx, axis=-1, mode="clip")
    hn = hn.reshape(S, Y, X, T, E, L, 2).transpose(0, 4, 5, 1, 2, 6, 3)
    allowed = jnp.asarray(data.g_allowed)[None, None, :, None, None, None, :]
    return jnp.where(allowed, hn, h)


# --------------------------------------------------------------------------- the stages
def solfp_ediffuse(gle, coeffs: CollCoeffs, data: CollData, params: CollParams):
    """Energy-diffusion implicit solve, solfp_ediffuse_le_layout (src/collisions.fpp:4037-4076).

    Per line ``(is, ik, it, ig, ixi)`` with ``ixi <= max(2*jend(ig), 2*ng2)``, skipped where
    ``vnew(ik,1,is) ~ 0`` (unless force_collisions), ``kwork_filter`` or ``ieqzip``; sets
    ``gle(ixi, negrid+1) = 0`` on solved lines.
    """
    E, NX = params.negrid, params.nx
    eq = jnp.moveaxis(_pad_xi(coeffs.eqle, NX), 4, 0)
    eb = jnp.moveaxis(_pad_xi(coeffs.ebetaale, NX), 4, 0)
    ec = jnp.moveaxis(_pad_xi(coeffs.ec1le, NX), 4, 0)
    gs = jnp.moveaxis(gle[..., :E, :], 4, 0)                       # (E, S, Y, X, T, NX)
    new = jnp.moveaxis(_tridiag(gs, eq, eb, ec), 0, 4)
    new = _set_energy_rows(jnp.zeros_like(gle), new)               # gle(ixi, negrid+1) = 0
    mask = jnp.asarray(data.apply)[:, :, :, None, None, None] & _bt(data.ed_rows)
    return jnp.where(mask, new, gle)


def conserve_diffuse(gle, coeffs: CollCoeffs, data: CollData, params: CollParams):
    """Energy-diffusion momentum/energy restoring terms, conserve_diffuse_le_layout
    (src/collisions.fpp:3541-3628); modes with ``kwork_filter`` or ``ieqzip`` are skipped (no vnew test)::

        v0y0 = sum ((vnE*aj0le)*gle)*wxi,         vnE = (vnmult(2)*vnew_E)*w;    gle -= bz0le*v0y0
        v1y1 = sum ((vpa_aj0_le*delnu)*gle)*wxi,  delnu = (vnmult(2)*delvnew)*w; gle -= bs0le*v1y1
        v2y2 = sum ((delnu*vperp_aj1le)*gle)*wxi;                                gle -= bw0le*v2y2
    """
    p = params
    E, XI, NX = p.negrid, p.nxi_lim, p.nx
    vnm2 = jnp.asarray(coeffs.vnmult)[1]
    W = jnp.asarray(data.w)[:, :, None]
    wxi = jnp.asarray(data.wxi)[:, :XI][None, None, None, :, None, :]
    act = jnp.asarray(data.active)[None, :, :, None, None, None]

    def moment(g, fac):
        return _moment(rmul(rmul(g[..., :E, :XI], fac), wxi), p)

    def restore(g, v, vec):
        return jnp.where(act, g - rmul(v[..., None, None], _pad_xi(vec, NX)), g)

    g = gle
    vnE = _sey((vnm2 * jnp.asarray(data.vnew_E)) * W)
    g = restore(g, moment(g, vnE * jnp.asarray(data.aj0le)[..., :E, :XI]), coeffs.bz0le)
    delnu = _sey((vnm2 * jnp.asarray(data.delvnew)) * W)
    g = restore(g, moment(g, jnp.asarray(data.vpa_aj0_le)[..., :E, :XI] * delnu), coeffs.bs0le)
    g = restore(g, moment(g, delnu * jnp.asarray(data.vperp_aj1le)[..., :E, :XI]), coeffs.bw0le)
    return g


def add_drag(gle, aparnew, coeffs: CollCoeffs, data: CollData, params: CollParams):
    """Electron-ion drag (resistivity) term of solfp1_le_layout (src/collisions.fpp:3056-3085)::

        tmp = ((((vnmult(1)*vnewk)*code_dt)*kperp2)*aparnew) / ((((-(z*dens))*beta)*stm)*energy**1.5)
        gle(ixi, ie) += tmp * vpa_aj0_le(ixi, ie)          ixi = 1..nxi_lim, ie = 1..negrid

    Electron-like species only; ``kwork_filter``/``ieqzip`` modes skipped; no vnew test. ``aparnew`` is
    the module ``aparnew`` GS2 holds at this point (spec 10 G-27, G-33).
    """
    p = params
    E, NX = p.negrid, p.nx
    vnm1 = jnp.asarray(coeffs.vnmult)[0]
    # GS2 quirk: no tunits factor here, while the z0 conservation vector has one
    # (collisions.fpp:3073-3075 vs 740)
    num = (vnm1 * jnp.asarray(data.vnewk)) * coeffs.code_dt               # (S,)
    num = num[:, None, None, None] * jnp.asarray(data.kperp2)[None]       # (S, Y, X, T)
    za = rmul(jnp.asarray(aparnew)[None], num)
    # XLA's algebraic simplifier rewrites ``A / broadcast(B)`` into ``A * broadcast(1/B)`` (a second
    # rounding, <= 1 ulp, that gfortran does not do); in strict mode the barrier hides the broadcast.
    den = ob(jnp.broadcast_to(jnp.asarray(data.drag_den)[:, None, None, None, :], za.shape + (E,)))
    tmp = _cdiv_r(za[..., None], den)                                    # (S, Y, X, T, E)
    vaj = _pad_xi(data.vpa_aj0_le, NX)[..., :E, :]
    new = gle[..., :E, :] + rmul(tmp[..., None], vaj)
    m = (jnp.asarray(data.is_electron)[:, None, None]
         & jnp.asarray(data.active)[None])[:, :, :, None, None, None]
    return _set_energy_rows(gle, jnp.where(m, new, gle[..., :E, :]))


def solfp_lorentz(gle, coeffs: CollCoeffs, data: CollData, params: CollParams):
    """Pitch-angle scattering implicit solve, solfp_lorentz_le_layout (src/collisions.fpp:3897-3980).

    Per line ``(is, ik, it, ig, ie)`` (same skips as the energy solve), with ``JE = max(2*jend, 2*ng2+1)``:
    average the two wfb vpar=0 copies (``vpar_zero_mean``), zero ``gle(JE)``, remove the wfb from the
    system (``special_wfb_lorentz``), sweep, reinsert, copy the bounce point ``gle(2*jend) = gle(jend)``.
    The sweep runs over all ``nxi+1`` rows: ``qle, betaale, c1le`` are 0 beyond the system size ``te2``
    and the rows beyond hold exact zeros, so the extra rows come out as (signed) zeros and the system
    rows are unchanged (spec 06 section 7). On grids without trapped particles this is the intended
    ``2*ng2`` system (spec 06 section 8.1).
    """
    p = params
    E, NX, ng2 = p.negrid, p.nx, p.ng2
    g0 = gle[..., :E, :]
    g = g0
    k = jnp.arange(NX)
    cj, cj2, je = (jnp.asarray(a) for a in (data.lz_cj, data.lz_cj2, data.lz_je))

    def onehot(idx):
        return k[None, :] == idx[:, None]

    if p.any_wfb_average:
        s = _take_xi(g, cj) + _take_xi(g, cj2)
        avg = _cdiv_r(s, 2.0)                                     # (g(je) + g(2je)) / 2.0
        g = jnp.where(_bt(onehot(cj) & jnp.asarray(data.lz_avg)[:, None]), avg, g)
    g = jnp.where(_bt(onehot(je)), jnp.zeros((), g.dtype), g)   # gle(je, ie) = 0
    if p.any_special_shift:
        gwfb = g[..., ng2:ng2 + 1]
        g = jnp.where(_bt(data.lz_shift), jnp.roll(g, -1, axis=-1), g)
        g = jnp.where(_bt(data.lz_shift_zero), jnp.zeros((), g.dtype), g)
    q = jnp.moveaxis(_pad_xi(coeffs.qle, NX), -1, 0)
    b = jnp.moveaxis(_pad_xi(coeffs.betaale, NX), -1, 0)
    c = jnp.moveaxis(_pad_xi(coeffs.c1le, NX), -1, 0)
    g = jnp.moveaxis(_tridiag(jnp.moveaxis(g, -1, 0), q, b, c), 0, -1)
    if p.any_special_shift:
        g = jnp.where(_bt(data.lz_unshift), jnp.roll(g, 1, axis=-1), g)
        g = jnp.where(_bt((k[None, :] == ng2) & jnp.asarray(data.lz_special)[:, None]), gwfb, g)
    if p.trapped:
        # GS2 quirk: the isgn=1 copy of every bounce point is overwritten by the scattered isgn=2
        # value; its own input is discarded except through the wfb average (collisions.fpp:3976)
        g = jnp.where(_bt(onehot(cj2) & jnp.asarray(data.lz_bounce)[:, None]), _take_xi(g, cj), g)
    g = jnp.where(jnp.asarray(data.apply)[:, :, :, None, None, None], g, g0)
    return _set_energy_rows(gle, g)


def conserve_lorentz(gle, coeffs: CollCoeffs, data: CollData, params: CollParams):
    """Lorentz momentum/energy restoring terms, conserve_lorentz_le_layout
    (src/collisions.fpp:3256-3376), ``conservative = .true.``; ``kwork_filter``/``ieqzip`` skipped::

        [drag] v0y0 = sum (((vpa_aj0_le*gle)*w)*wxi);                          gle -= z0le*v0y0
        v1y1 = sum ((((vpdiffle*nud)*aj0le)*gle)*wxi), nud = ((speed*vnmult(1))*vnew_D)*w;  gle -= s0le*v1y1
        v2y2 = sum (((nud*vperp_aj1le)*gle)*wxi),     nud = (vnmult(1)*vnew_D)*w;          gle -= w0le*v2y2
    """
    p = params
    E, XI, NX = p.negrid, p.nxi_lim, p.nx
    vnm1 = jnp.asarray(coeffs.vnmult)[0]
    W = jnp.asarray(data.w)
    wxi = jnp.asarray(data.wxi)[:, :XI][None, None, None, :, None, :]
    act = jnp.asarray(data.active)[None, :, :, None, None, None]

    def restore(g, v, vec):
        return jnp.where(act, g - rmul(v[..., None, None], _pad_xi(vec, NX)), g)

    g = gle
    if p.drag:
        term = rmul(rmul(rmul(g[..., :E, :XI], jnp.asarray(data.vpa_aj0_le)[..., :E, :XI]),
                         W[:, None, None, None, :, None]), wxi)
        g = restore(g, _moment(term, p), coeffs.z0le)
    vnew_D = jnp.asarray(data.vnew_D)
    nud = _sey(((jnp.asarray(data.speed)[None, :, None] * vnm1) * vnew_D) * W[:, :, None])
    fac = (jnp.asarray(data.vpdiffle)[None, None, None, :, None, :] * nud) \
        * jnp.asarray(data.aj0le)[..., :E, :XI]
    g = restore(g, _moment(rmul(rmul(g[..., :E, :XI], fac), wxi), p), coeffs.s0le)
    nud = _sey((vnm1 * vnew_D) * W[:, :, None])
    fac = nud * jnp.asarray(data.vperp_aj1le)[..., :E, :XI]
    g = restore(g, _moment(rmul(rmul(g[..., :E, :XI], fac), wxi), p), coeffs.w0le)
    return g


def solfp1(gle, aparnew, coeffs: CollCoeffs, data: CollData, params: CollParams, *,
           return_stages: bool = False):
    """solfp1_le_layout (src/collisions.fpp:3035-3090): energy diffusion, its conservation, drag,
    Lorentz, its conservation. With ``return_stages`` also returns the state after each stage, keyed
    by the probe tags of spec 09 section 5.2 (a tag is recorded even when its optional stage is off)."""
    p = params
    stages = {}
    if p.has_diffuse:
        gle = solfp_ediffuse(gle, coeffs, data, p)
        stages["coll_after_ediffuse"] = gle
        if p.conserve_moments:
            gle = conserve_diffuse(gle, coeffs, data, p)
        stages["coll_after_ediffuse_conserve"] = gle
    if p.has_lorentz:
        if p.drag:
            gle = add_drag(gle, aparnew, coeffs, data, p)
        stages["coll_after_drag"] = gle
        gle = solfp_lorentz(gle, coeffs, data, p)
        stages["coll_after_lorentz"] = gle
        if p.conserve_moments:
            gle = conserve_lorentz(gle, coeffs, data, p)
        stages["coll_after_lorentz_conserve"] = gle
    return (gle, stages) if return_stages else gle


# --------------------------------------------------------------------------- public kernels
def vspace_derivatives(gnew, g, fields_old: Fields, fields_new: Fields, aparnew,
                       diagnostics: Optional[int], coll_coeffs: CollCoeffs, data: CollData,
                       params: CollParams, *, return_stages: bool = False):
    """Collisions applied to ``gnew``: vspace_derivatives (src/le_derivatives.fpp:15-174), le layout
    branch (67-112), with solfp1_le_layout (src/collisions.fpp:3035-3090).

    ``gnew`` G7 complex; ``fields_new`` holds ``phinew, bparnew`` for the g <-> h transforms;
    ``aparnew`` (F3) is the field the drag term reads (the module ``aparnew``: ``apar`` in pass 1, the
    solved one in pass 2, the probe during a response build; spec 10 G-27). ``g``, ``fields_old`` and
    ``diagnostics`` only feed the heating diagnostic, which is inert and not ported (spec 10 G-26).
    Returns the new ``gnew`` (and, with ``return_stages``, a dict of the probe-point states of spec 09
    section 5.2: ``coll_h_in`` ... ``coll_h_out``).
    """
    del g, fields_old, diagnostics          # heating only (c_rate), not ported
    p = params
    if not p.colls:                          # le path: nothing at all without collisions (LD:68)
        return (gnew, {}) if return_stages else gnew
    h = gnew
    if p.adjust:
        h = g_adjust_direction(h, fields_new.phi, fields_new.bpar, FROM_G, p.fphi, p.fbpar,
                               data.gadj)
    gle = g2le(h, data, p)
    gle_out, stages = solfp1(gle, aparnew, coll_coeffs, data, p, return_stages=True)
    h_out = le2g(gle_out, h, data, p)
    out = h_out
    if p.adjust:
        out = g_adjust_direction(h_out, fields_new.phi, fields_new.bpar, TO_G, p.fphi, p.fbpar,
                                 data.gadj)
    if not return_stages:
        return out
    stages = {"coll_h_in": h, "coll_gle_in": gle, **stages, "coll_gle_out": gle_out,
              "coll_h_out": h_out}
    return out, stages


def collisions_advance(gnew, g, fields_old: Fields, fields_new: Fields,
                       potentials: Callable[..., Fields], coll_coeffs: CollCoeffs, data: CollData,
                       params: CollParams, diagnostics: Optional[int] = 1):
    """Split collisions after the second ``timeadv``: collisions_advance (src/dist_fn.fpp:4130-4167)::

        vspace_derivatives(gnew, ...)                      (drag reads the solved aparnew)
        h = g_adjust(gnew, phinew, bparnew, from_g)
        phinew, aparnew, bparnew = calculate_potentials_from_nonadiabatic_dfn(h)
        gnew = g_adjust(h, phinew, bparnew, to_g)          (with the new fields)

    ``potentials(h) -> Fields`` is ``calculate_potentials_from_nonadiabatic_dfn`` (spec 04 section 3.9,
    owner B4). Returns ``(gnew, fields)``; a no-op returning its inputs unless ``split_collisions``.
    The g_adjust calls here do not depend on ``adjust`` (dist_fn.fpp:4164-4166).
    """
    p = params
    if not p.split_collisions:
        return gnew, fields_new
    gnew = vspace_derivatives(gnew, g, fields_old, fields_new, fields_new.apar, diagnostics,
                              coll_coeffs, data, p)
    h = g_adjust_direction(gnew, fields_new.phi, fields_new.bpar, FROM_G, p.fphi, p.fbpar,
                           data.gadj)
    new = potentials(h)
    if not isinstance(new, Fields):
        new = Fields(*new)
    gnew = g_adjust_direction(h, new.phi, new.bpar, TO_G, p.fphi, p.fbpar, data.gadj)
    return gnew, new


@partial(jax.jit, static_argnames=("diagnostics", "params", "return_stages"))
def _vspace_jit(gnew, g, fields_old, fields_new, aparnew, diagnostics, coll_coeffs, data, params,
                return_stages=False):
    return vspace_derivatives(gnew, g, fields_old, fields_new, aparnew, diagnostics, coll_coeffs,
                              data, params, return_stages=return_stages)


def vspace_derivatives_jit(gnew, g, fields_old: Fields, fields_new: Fields, aparnew,
                           diagnostics: Optional[int], coll_coeffs: CollCoeffs, data: CollData,
                           params: CollParams, *, return_stages: bool = False):
    """``jax.jit`` of ``vspace_derivatives`` (``params``/``diagnostics`` static; ``data.init`` dropped
    so the host-only rebuild arrays never go to the device)."""
    return _vspace_jit(gnew, g, fields_old, fields_new, aparnew, diagnostics, coll_coeffs,
                       data._replace(init=None), params, return_stages=return_stages)
