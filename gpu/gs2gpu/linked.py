"""Linked (twist-and-shift) parallel boundary conditions (spec 03 section 3.11).

Mirrors ``apply_linked_boundary_conditions`` (src/dist_fn.fpp:5919-6149) with
the communication patterns ``links_p``/``links_h``
(``setup_connected_bc_redistribute``, 2759-3022) and ``wfb_p``/``wfb_h``
(``setup_connected_bc_redistribute_mixed_wfb``, 3028-3224) replaced by
gathers along the supercell chains. Line numbers refer to GS2 8.2.1, commit
4d8c94bcf (the specs' numbering).

Chains (supercells) are built on the host from the exported ``itright``,
``l_links``, ``r_links`` (``compute_itleft_and_itright`` 2441-2513,
``count_links`` 2608-2672) and padded to a rectangular ``(C, Nmax)`` table.
Only chains with ``N >= 2`` cells take part (``ncell == 1`` cells are skipped
at dist_fn.fpp:6002, and ky = 0 box cells are self-linked with ``N = 1``).

Per cell at chain position ``k`` (0-based from the left) with endpoint values
``p1[m] = gnew(ntgrid, 1, cell m)``, ``h1[m] = g_h(ntgrid, 1, cell m)``,
``p2[m] = gnew(-ntgrid, 2, cell m)``, ``h2[m] = g_h(-ntgrid, 2, cell m)``,
all taken from the **pre-correction** ``gnew`` (the fills precede the loop):

* passing (and passing wfb), open chain (6115-6144): sign 1 if
  ``save_h(1) .and. l_links > 0``::

      b0 = p1[k-1]; fac = 1
      for i = 2..k: fac = fac*h1[k-i+1]; b0 = b0 + p1[k-i]*fac
      gnew(:, 1) += b0*g_h(:, 1)

  nearest neighbour first, each cell's own sum (not a Horner recurrence).
  Sign 2 is the mirror image in chain order reversed (right to left).
* mixed wfb, cyclic chain (6013-6114): ``facd = 1/(1 - prod(h))``,
  ``b0 = sum``, ``A0 = ((b0 + p_last) -/+ dadj)*facd`` and the cascade
  ``A <- A*h[m] + p[m]`` over the upstream cells (a forward recurrence whose
  prefixes are shared by the chain).

Both signs are written once, in "upstream order": sign 1 uses the chain left
to right, sign 2 the chain right to left; with that relabelling the Fortran
loops of the two signs are identical (spec 03 section 3.11.2-3.11.3; the
slot bookkeeping of ``g_adj`` was checked against the redistribute set-up).
"""

from __future__ import annotations

from typing import NamedTuple, Optional

import jax.numpy as jnp
import numpy as np
from jax import lax

from .numerics import cmul, recip, rmul, rsub


class LinkedData(NamedTuple):
    """Static chain tables and masks for the linked boundary conditions.

    ``chain_y`` ``(C,)``: 0-based ik of chain ``c``; ``chain_x`` /
    ``chain_xr`` ``(C, Nmax)``: 0-based it of the cells left-to-right /
    right-to-left (padding repeats a valid cell; padded positions are never
    used); ``chain_n`` ``(C,)``: number of cells ``N >= 2``.
    ``cell_chain, cell_pos, cell_posr`` ``(Y, X)``: chain id and position from
    the left / right of every cell (0 for cells in no chain; masked).
    ``m_open1/2`` ``(L, Y, X)``: the passing open-chain correction applies to
    sign 1/2; ``m_cyc1/2`` ``(L, Y, X)``: the mixed-wfb cyclic correction
    applies (only at ``il = ng2+1``). The masks already contain
    ``any(save_h)``, ``.not. kwork_filter`` and ``ncell > 1`` (6000-6003);
    ``hybrid`` ``(S, Y, L)`` is the passing-hybrid skip (6005).
    Mixed-wfb ``dadj`` inputs (dist_fn.fpp:6013-6047, per G-04 pre-gathered at
    the chain ends): ``aj0l, aj1l, aj0r, aj1r`` ``(S, E, C)``, ``vperp2left``,
    ``vperp2right`` ``(S, E)``, ``zt`` ``(S,)``, ``nonmaxw_corr`` ``(S, E)``;
    ``None`` unless ``mixed_wfb .and. nonad_zero``.
    """

    chain_y: np.ndarray
    chain_x: np.ndarray
    chain_xr: np.ndarray
    chain_n: np.ndarray
    cell_chain: np.ndarray
    cell_pos: np.ndarray
    cell_posr: np.ndarray
    m_open1: np.ndarray
    m_open2: np.ndarray
    m_cyc1: np.ndarray
    m_cyc2: np.ndarray
    hybrid: np.ndarray
    aj0l: Optional[np.ndarray]
    aj1l: Optional[np.ndarray]
    aj0r: Optional[np.ndarray]
    aj1r: Optional[np.ndarray]
    vperp2left: Optional[np.ndarray]
    vperp2right: Optional[np.ndarray]
    zt: Optional[np.ndarray]
    nonmaxw_corr: Optional[np.ndarray]

    @classmethod
    def from_export(cls, ex) -> Optional["LinkedData"]:
        """Linked tables of an export (None unless linked with connections)."""
        from .invert_rhs import InvertData  # local import: invert_rhs imports this module

        return InvertData.from_export(ex).linked


# --------------------------------------------------------------------------- host side
def build_chains(itleft, itright, l_links, r_links):
    """Supercell chains ``[(ik, [it_0, ..., it_{N-1}]), ...]`` (0-based), every cell once.

    Heads are the cells with ``l_links == 0``; a chain follows ``itright``
    ``r_links(head)`` times (count_links_one_side, dist_fn.fpp:2608-2640;
    get_leftmost_it/get_rightmost_it, 3889-3920). The links are checked for
    consistency (``l_links`` increases by one per hop, ``itleft`` points back).
    """
    itleft = np.asarray(itleft)
    itright = np.asarray(itright)
    l_links = np.asarray(l_links)
    r_links = np.asarray(r_links)
    naky, ntheta0 = itright.shape
    chains = []
    seen = np.zeros((naky, ntheta0), bool)
    for ik in range(naky):
        for it in range(ntheta0):
            if int(l_links[ik, it]) != 0:
                continue
            n = int(r_links[ik, it]) + 1
            cells = [it]
            cur = it
            for k in range(1, n):
                nxt = int(itright[ik, cur]) - 1          # Fortran value is 1-based
                if nxt < 0 or nxt == cur:
                    raise ValueError(f"linked chain broken at ik={ik + 1} it={cur + 1}")
                if int(l_links[ik, nxt]) != k or int(itleft[ik, nxt]) - 1 != cur:
                    raise ValueError(f"inconsistent links at ik={ik + 1} it={nxt + 1}")
                cells.append(nxt)
                cur = nxt
            if int(r_links[ik, cur]) != 0:
                raise ValueError(f"chain from ik={ik + 1} it={it + 1} does not end at r_links=0")
            for c in cells:
                if seen[ik, c]:
                    raise ValueError(f"cell ik={ik + 1} it={c + 1} in two chains")
                seen[ik, c] = True
            chains.append((ik, cells))
    if not seen.all():
        raise ValueError("some (ik, it) cells belong to no chain")
    return chains


def build_linked(*, itleft, itright, l_links, r_links, aky, kwork_filter, ng2: int,
                 nlambda: int, passing_wfb: bool, trapped_wfb: bool, mixed_wfb: bool,
                 nonad_zero: bool, hybrid, zt=None, nonmaxw_corr=None, wfb_factors=None
                 ) -> tuple[Optional[LinkedData], int]:
    """Host construction of :class:`LinkedData`; returns ``(data or None, Nmax)``.

    ``hybrid`` is ``is_passing_hybrid_electron`` ``(S, Y, L)``.
    ``wfb_factors`` is a dict with the exported ``wfb_aj0l, wfb_aj1l, wfb_aj0r,
    wfb_aj1r`` ``(S, Y, X, E)`` and ``wfb_vperp2left, wfb_vperp2right``
    ``(S, E)`` (needed only for ``mixed_wfb .and. nonad_zero``).
    Masks mirror compute_connections (2518-2565) and set_save_h (2570-2604).
    """
    itleft = np.asarray(itleft)
    itright = np.asarray(itright)
    l_links = np.asarray(l_links)
    r_links = np.asarray(r_links)
    aky = np.asarray(aky)
    kwork_filter = np.asarray(kwork_filter, bool)
    naky, ntheta0 = itright.shape

    chains = [ch for ch in build_chains(itleft, itright, l_links, r_links) if len(ch[1]) >= 2]
    if not chains:
        return None, 0
    nmax = max(len(cells) for _, cells in chains)
    nch = len(chains)
    chain_y = np.zeros(nch, np.int32)
    chain_x = np.zeros((nch, nmax), np.int32)
    chain_xr = np.zeros((nch, nmax), np.int32)
    chain_n = np.zeros(nch, np.int32)
    cell_chain = np.zeros((naky, ntheta0), np.int32)
    cell_pos = np.zeros((naky, ntheta0), np.int32)
    cell_posr = np.zeros((naky, ntheta0), np.int32)
    for c, (ik, cells) in enumerate(chains):
        n = len(cells)
        chain_y[c] = ik
        chain_n[c] = n
        padded = cells + [cells[-1]] * (nmax - n)
        chain_x[c] = padded
        rev = cells[::-1]
        chain_xr[c] = rev + [rev[-1]] * (nmax - n)
        for k, it in enumerate(cells):
            cell_chain[ik, it] = c
            cell_pos[ik, it] = k
            cell_posr[ik, it] = n - 1 - k

    il = np.arange(1, nlambda + 1)
    is_passing = il <= ng2
    is_wfb = il == ng2 + 1
    is_trapped = il > ng2 + 1
    # compute_connections: trapped (and trapped wfb) have no connections
    excluded = is_trapped | (is_wfb & trapped_wfb)
    left = (~excluded)[:, None, None] & (itleft >= 0)[None]       # iglo_left >= 0
    right = (~excluded)[:, None, None] & (itright >= 0)[None]
    neighbor = left | right
    aky_nz = (aky != 0.0)[None, :, None]
    save1 = left & aky_nz
    save2 = right & aky_nz
    wfb_all = is_wfb[:, None, None] & neighbor & (not passing_wfb) & aky_nz
    save1 = save1 | wfb_all
    save2 = save2 | wfb_all
    ncell = l_links + r_links + 1
    base = (save1 | save2) & (~kwork_filter)[None] & (ncell > 1)[None]
    cyc_l = (is_wfb & mixed_wfb)[:, None, None]
    m_open1 = base & ~cyc_l & save1 & (l_links > 0)[None]
    m_open2 = base & ~cyc_l & save2 & (r_links > 0)[None]
    m_cyc1 = base & cyc_l & save1
    m_cyc2 = base & cyc_l & save2
    # the open-chain correction only exists for passing-like pitch angles
    passing_like = is_passing | (is_wfb & passing_wfb)
    assert not (m_open1 | m_open2)[~passing_like].any()

    aj0l = aj1l = aj0r = aj1r = vp2l = vp2r = None
    zt_ = nm_ = None
    if mixed_wfb and nonad_zero and bool((m_cyc1 | m_cyc2).any()):
        if wfb_factors is None:
            raise ValueError("mixed wfb with nonad_zero needs the exported wfb end factors "
                             "(/dist_fn/wfb_aj0l ... wfb_vperp2right)")
        x_left = chain_x[:, 0]
        x_right = chain_x[np.arange(nch), chain_n - 1]

        def at(name, xs):
            arr = np.asarray(wfb_factors[name])                     # (S, Y, X, E)
            return np.ascontiguousarray(arr[:, chain_y, xs, :].transpose(0, 2, 1))  # (S, E, C)

        # G-04: the factors are already evaluated at itl/itr; reading them at
        # the chain's own end cells is exact either way.
        aj0l, aj1l = at("wfb_aj0l", x_left), at("wfb_aj1l", x_left)
        aj0r, aj1r = at("wfb_aj0r", x_right), at("wfb_aj1r", x_right)
        vp2l = np.asarray(wfb_factors["wfb_vperp2left"], np.float64)
        vp2r = np.asarray(wfb_factors["wfb_vperp2right"], np.float64)
        zt_ = np.asarray(zt, np.float64)
        nm_ = np.asarray(nonmaxw_corr, np.float64)

    data = LinkedData(
        chain_y=chain_y, chain_x=chain_x, chain_xr=chain_xr, chain_n=chain_n,
        cell_chain=cell_chain, cell_pos=cell_pos, cell_posr=cell_posr,
        m_open1=m_open1, m_open2=m_open2, m_cyc1=m_cyc1, m_cyc2=m_cyc2,
        hybrid=np.asarray(hybrid, bool),
        aj0l=aj0l, aj1l=aj1l, aj0r=aj0r, aj1r=aj1r, vperp2left=vp2l, vperp2right=vp2r,
        zt=zt_, nonmaxw_corr=nm_)
    return data, nmax


# --------------------------------------------------------------------------- device side
def _open_chain(p, h):
    """Passing open-chain amplitude ``b0`` for every chain position (6121-6143).

    ``p``, ``h``: ``(..., C, Nmax)`` endpoint values in upstream order. For
    position ``k``: ``b0 = p[k-1]``, then for ``i = 2..k``
    ``fac = fac*h[k-i+1]; b0 = b0 + p[k-i]*fac``. Position 0 and padded
    positions return garbage (masked by the caller).
    """
    nmax = p.shape[-1]
    pos = jnp.arange(nmax)
    b0 = jnp.take(p, jnp.clip(pos - 1, 0, nmax - 1), axis=-1)
    fac = jnp.ones_like(b0)                                  # fac = 1.0

    def body(i, carry):
        fac, b0 = carry
        active = i <= pos
        hv = jnp.take(h, jnp.clip(pos - i + 1, 0, nmax - 1), axis=-1)
        pv = jnp.take(p, jnp.clip(pos - i, 0, nmax - 1), axis=-1)
        fac_n = cmul(fac, hv)                                # fac = fac * g_adj(n+2-i)
        b0_n = b0 + cmul(pv, fac_n)                          # b0 = b0 + g_adj(i)*fac
        return jnp.where(active, fac_n, fac), jnp.where(active, b0_n, b0)

    if nmax > 2:
        fac, b0 = lax.fori_loop(2, nmax, body, (fac, b0))
    return b0


def _cyclic_chain(p, h, n, dadj, plus: bool):
    """Mixed-wfb cyclic amplitudes for every chain position (6049-6114).

    ``p``, ``h``: ``(S, E, C, Nmax)`` in upstream order; ``n``: ``(C,)``;
    ``dadj``: ``(S, E, C)``; ``plus`` selects ``+dadj`` (sign 2) or ``-dadj``
    (sign 1). Returns ``(S, E, C, Nmax)``: the amplitude after ``k`` cascade
    steps for position ``k``.
    """
    nmax = p.shape[-1]
    pos = jnp.arange(nmax)
    valid = pos[None, :] < n[:, None]                        # (C, Nmax)
    one = jnp.ones(p.shape[:-1], p.dtype)

    # facd = 1./(1. - product(h)) ; PRODUCT is sequential, starting from 1
    def prod_step(prod, x):
        hm, vm = x
        return jnp.where(vm, cmul(prod, hm), prod), None

    prod, _ = lax.scan(prod_step, one, (jnp.moveaxis(h, -1, 0), valid.T))
    facd = recip(rsub(1.0, prod))

    # b0 = 0. ; fac = 1 ; do i = ncell-1, 1, -1: fac = fac*h[i]; b0 = b0 + fac*p[i-1]
    def b_step(carry, i):
        fac, b0 = carry
        active = i <= n - 1                                  # (C,)
        fac_n = cmul(fac, jnp.take(h, i, axis=-1))
        b0_n = b0 + cmul(fac_n, jnp.take(p, i - 1, axis=-1))
        return (jnp.where(active, fac_n, fac), jnp.where(active, b0_n, b0)), None

    (_, b0), _ = lax.scan(b_step, (one, jnp.zeros_like(one)), jnp.arange(nmax - 1, 0, -1))
    idx = jnp.broadcast_to((n - 1)[None, None, :, None], p.shape[:-1] + (1,))
    p_last = jnp.take_along_axis(p, idx, axis=-1)[..., 0]
    s = b0 + p_last
    s = (s + dadj) if plus else (s - dadj)
    a0 = cmul(s, facd)                                       # (b0 + g_adj(1) -/+ dadj)*facd

    # cascade: b0 = b0*h[m] + p[m] for the upstream cells m = 0..k-1
    def c_step(a, x):
        hm, pm = x
        return cmul(a, hm) + pm, a

    a_last, amps = lax.scan(c_step, a0, (jnp.moveaxis(h[..., :-1], -1, 0),
                                         jnp.moveaxis(p[..., :-1], -1, 0)))
    amps = jnp.concatenate([amps, a_last[None]], axis=0)     # (Nmax, S, E, C)
    return jnp.moveaxis(amps, 0, -1)


def _wfb_dadj(phinew, bparnew, ld: LinkedData, params):
    """``dadj = adjl - adjr`` of the mixed wfb in a chain (dist_fn.fpp:6013-6047), ``(S, E, C)``.

    ``adjl = 2.0*vperp2left*aj1l*bparnew(-ntgrid,itl,ik)*fbpar
    + spec(is)%zt*phinew(-ntgrid,itl,ik)*aj0l*nonmaxw_corr(ie,is)*fphi``
    (left to right), ``adjr`` likewise at ``(ntgrid, itr)``.
    """
    nchain = ld.chain_n.shape[0]
    cy = jnp.asarray(ld.chain_y)
    xl = jnp.asarray(ld.chain_x)[:, 0]
    xr = jnp.asarray(ld.chain_x)[jnp.arange(nchain), jnp.asarray(ld.chain_n) - 1]
    zt = jnp.asarray(ld.zt)[:, None, None]
    nm = jnp.asarray(ld.nonmaxw_corr)[:, :, None]

    def adj(field_b, field_p, x, t, vperp2_end, aj0, aj1):
        bp = field_b[cy, x, t][None, None, :]
        ph = field_p[cy, x, t][None, None, :]
        coef = (2.0 * jnp.asarray(vperp2_end)[:, :, None]) * aj1
        term_b = rmul(rmul(bp, coef), params.fbpar)
        term_p = rmul(rmul(rmul(rmul(ph, zt), aj0), nm), params.fphi)
        return term_b + term_p

    nth = phinew.shape[-1]
    adjl = adj(bparnew, phinew, xl, 0, ld.vperp2left, jnp.asarray(ld.aj0l), jnp.asarray(ld.aj1l))
    adjr = adj(bparnew, phinew, xr, nth - 1, ld.vperp2right, jnp.asarray(ld.aj0r),
               jnp.asarray(ld.aj1r))
    return adjl - adjr


def apply_linked_boundary_conditions(gnew, phinew, bparnew, coeffs, data, params):
    """Mirrors ``apply_linked_boundary_conditions`` (src/dist_fn.fpp:5919-6149).

    ``gnew``: G7 output of ``invert_rhs_1``; ``phinew``, ``bparnew``: F3
    fields of this ``timeadv`` call; ``coeffs``: ``types.Coeffs`` (uses
    ``g_h``; computed with ``invert_rhs.init_homogeneous_g`` if None);
    ``data``: ``InvertData`` (or a :class:`LinkedData`), ``params``:
    ``InvertParams``. Returns the corrected ``gnew``. A no-op unless
    ``boundary_option = 'linked'`` with connections (``no_connections``
    return, 5948). ``start_from_previous_solution`` is gated (refused by
    ``InvertParams.from_export``). ``ieqzip`` is not consulted here, as in GS2.
    """
    if not params.linked or params.no_connections:
        return gnew
    ld = data.linked if hasattr(data, "linked") else data
    if ld is None or ld.chain_n.shape[0] == 0:
        return gnew
    g_h = coeffs.g_h
    if g_h is None:
        if not hasattr(data, "linked"):
            raise ValueError("coeffs.g_h is None; pass InvertData so it can be computed")
        from .invert_rhs import init_homogeneous_g
        g_h = init_homogeneous_g(coeffs, data, params)
    gnew = jnp.asarray(gnew)
    g_h = jnp.asarray(g_h)
    phinew = jnp.asarray(phinew)
    bparnew = jnp.asarray(bparnew)
    nth = gnew.shape[-1]
    ng2 = params.ng2
    nlambda = gnew.shape[2]

    cy = jnp.asarray(ld.chain_y)[:, None]
    cx = jnp.asarray(ld.chain_x)
    cxr = jnp.asarray(ld.chain_xr)
    n = jnp.asarray(ld.chain_n)
    # endpoint gathers (S, E, L, C, Nmax): the links_p/wfb_p and links_h/wfb_h fills
    p1 = gnew[:, :, :, cy, cx, 0, nth - 1]
    h1 = g_h[:, :, :, cy, cx, 0, nth - 1]
    p2 = gnew[:, :, :, cy, cxr, 1, 0]
    h2 = g_h[:, :, :, cy, cxr, 1, 0]

    cc = jnp.asarray(ld.cell_chain)
    cp = jnp.asarray(ld.cell_pos)
    cpr = jnp.asarray(ld.cell_posr)
    s_, e_ = gnew.shape[0], gnew.shape[1]
    ny, nx = gnew.shape[3], gnew.shape[4]

    pieces1, pieces2 = [], []
    n_open = params.n_open_l
    if n_open > 0:
        b1 = _open_chain(p1[:, :, :n_open], h1[:, :, :n_open])
        b2 = _open_chain(p2[:, :, :n_open], h2[:, :, :n_open])
        pieces1.append(b1[:, :, :, cc, cp])                   # (S, E, n_open, Y, X)
        pieces2.append(b2[:, :, :, cc, cpr])
    done = n_open
    if params.mixed_wfb and params.has_wfb:
        assert done == ng2
        if params.nonad_zero and ld.aj0l is not None:
            dadj = _wfb_dadj(phinew, bparnew, ld, params)
        else:
            dadj = jnp.zeros((s_, e_, n.shape[0]), gnew.dtype)   # dadj = 0.0
        a1 = _cyclic_chain(p1[:, :, ng2], h1[:, :, ng2], n, dadj, plus=False)
        a2 = _cyclic_chain(p2[:, :, ng2], h2[:, :, ng2], n, dadj, plus=True)
        pieces1.append(a1[:, :, cc, cp][:, :, None])
        pieces2.append(a2[:, :, cc, cpr][:, :, None])
        done += 1
    if done < nlambda:
        pad = jnp.zeros((s_, e_, nlambda - done, ny, nx), gnew.dtype)
        pieces1.append(pad)
        pieces2.append(pad)
    amp1 = jnp.concatenate(pieces1, axis=2)                   # (S, E, L, Y, X)
    amp2 = jnp.concatenate(pieces2, axis=2)

    # GS2 quirk: ieqzip cells are not skipped here although invert_rhs_1 skips them
    # (dist_fn.fpp:5996-6010; 03 §6.7; eqzip_option is gated, G-07)
    not_hyb = ~jnp.asarray(ld.hybrid).transpose(0, 2, 1)[:, None, :, :, None]   # (S,1,L,Y,1)
    mask1 = (jnp.asarray(ld.m_open1) | jnp.asarray(ld.m_cyc1))[None, None] & not_hyb
    mask2 = (jnp.asarray(ld.m_open2) | jnp.asarray(ld.m_cyc2))[None, None] & not_hyb
    g1 = gnew[..., 0, :]
    g2 = gnew[..., 1, :]
    g1 = jnp.where(mask1[..., None], g1 + cmul(amp1[..., None], g_h[..., 0, :]), g1)
    g2 = jnp.where(mask2[..., None], g2 + cmul(amp2[..., None], g_h[..., 1, :]), g2)
    return jnp.stack([g1, g2], axis=-2)
