"""Implicit response matrix ``aminv``: storage, construction by unit probes, inversion.

Mirrors ``init_response_matrix`` / ``init_response_row`` / ``init_inverse_matrix``
(src/fields_implicit.f90:459-1071; spec 04 sections 3.12-3.15 and 5.4-5.5; spec 10 G-17, G-18,
G-27, G-28).

Index orders per class ``c`` (``N = N_c`` cells, ``nf`` present fields, ``nth`` theta points,
``R = N*nf*nth``), 0-based:

* **K** (``am`` / ``f_lo``): ``K(ig, f, n) = f + nf*(igp + nth*n)`` (field fastest);
* **B** (``aminv`` / ``fl``): ``B(ig, f, n) = igp + nth*(f + nf*n)`` (theta fastest);
* skip set ``S``: ``{igp = nth-1, n < N-1}`` (all ``f``), the duplicated right-boundary points.

``AMT[m, j, k] = am(k+1, idx(j+1, m))`` is the response (equation ``k``, K order) of member ``m``
to the unit perturbation ``j`` (K order); :func:`gauss_jordan` is the literal pivot-free
Gauss-Jordan of ``init_inverse_matrix`` on it, and the runtime matrix is
``aminv_c[m, B(o), B(i)] = AIT[m, K(o), K(i)]`` (:class:`Aminv`).
"""

from __future__ import annotations

import math
from typing import Callable, NamedTuple, Optional

import jax
import jax.numpy as jnp
import numpy as np
from jax import lax

from .fields import (FIELDOPT_LOCAL, FieldData, FieldParams, Topology, _first, _opt, _scalar,
                     getan, getfieldeq, stack_present)
from .io.export import Export
from .numerics import cmul, ob, strict
from .types import Fields


class Aminv(NamedTuple):
    """Inverse response matrices, one per class: ``mats[c]`` is ``(M_c, R_c, R_c)`` complex,
    ``[member, output row (B order), input column (B order)]`` -- the rows
    ``aminv(c)%dcell(:)%supercell(:)`` of fields_implicit.f90:24-47 (spec 09 section 7.4).
    The class/member order is that of ``FieldData.topo``. ``code_dt`` is the time step the
    matrices were built for (NaN if unknown)."""

    mats: tuple
    code_dt: float = math.nan

    @classmethod
    def from_export(cls, ex: Export, rung: Optional[int] = None) -> "Aminv":
        """``/fields_implicit/aminv_c<c>`` (or ``/dt_ladder/k<rung>/fields_implicit/...``).

        For ``field_option='local'`` runs this returns the exported local matrices
        (``/fields_local/resp_ik<ik>_is<is>``) converted to the same class/B layout
        (:func:`gs2gpu.fields_local.local_aminv_from_export`)."""
        fsw = _scalar(_first(ex, "/fields/fieldopt_switch"), int)
        base_path = "/fields_implicit" if rung is None else f"/dt_ladder/k{rung}/fields_implicit"
        if fsw == FIELDOPT_LOCAL or (_opt(ex, base_path + "/aminv_c1") is None
                                     and _opt(ex, "/fields_local") is not None):
            from .fields_local import local_aminv_from_export
            return local_aminv_from_export(ex, rung=rung)
        nclass = int(np.asarray(_first(ex, "/dist_fn/N_class")).size)
        mats = tuple(np.asarray(ex[f"{base_path}/aminv_c{c + 1}"], np.complex128)
                     for c in range(nclass))
        if rung is None:
            code_dt = _scalar(_first(ex, "/state_t0/code_dt", "/gs2_time/code_dt"), float)
        else:
            code_dt = _scalar(ex[f"/dt_ladder/k{rung}"].attrs.get("code_dt"), float)
        return cls(mats=mats, code_dt=math.nan if code_dt is None else float(code_dt))


# --------------------------------------------------------------------------- index maps (host)
def k_index(N: int, nf: int, nth: int):
    """``(n, igp, f)`` of every K index ``j`` (arrays of length ``R``)."""
    j = np.arange(N * nf * nth)
    return j // (nf * nth), (j // nf) % nth, j % nf


def skip_set(N: int, nf: int, nth: int) -> np.ndarray:
    """``iskip``/``jskip`` of init_inverse_matrix (fields_implicit.f90:770-774, 816-821) in K
    order: ``N > 1 .and. ig == ntgrid .and. n < N``."""
    n, igp, _ = k_index(N, nf, nth)
    return (N > 1) & (igp == nth - 1) & (n < N - 1)


def k_of_b(N: int, nf: int, nth: int) -> np.ndarray:
    """``kofb[b] = K`` index of the point with B index ``b`` (fields_implicit.f90:1036-1063)."""
    b = np.arange(N * nf * nth)
    igp = b % nth
    f = (b // nth) % nf
    n = b // (nf * nth)
    return (f + nf * (igp + nth * n)).astype(np.int32)


def reorder_k_to_b(AIT, nf: int, nth: int, N: int):
    """``aminv_c[m, B(o), B(i)] = AIT[m, K(o), K(i)]`` (fields_implicit.f90:1018-1063)."""
    kofb = k_of_b(N, nf, nth)
    return AIT[:, kofb][:, :, kofb]


def reorder_b_to_k(mat, nf: int, nth: int, N: int):
    """Inverse of :func:`reorder_k_to_b`."""
    kofb = k_of_b(N, nf, nth)
    bofk = np.argsort(kofb)
    return mat[:, bofk][:, :, bofk]


# --------------------------------------------------------------------------- inversion
def _smith_factors(b):
    """Divisor part of GCC's Smith division (numerics.cdiv; spec 04 P2), hoisted as GCC does
    for a fixed divisor: branch ``t = |br| < |bi|``, ratio ``r`` and denominator ``den``."""
    br, bi = b.real, b.imag
    t = jnp.abs(br) < jnp.abs(bi)
    r = jnp.where(t, br / bi, bi / br)
    den = jnp.where(t, ob(br * r) + bi, ob(bi * r) + br)
    return t, r, den


def _bdiv(x, d):
    """``x / d``. XLA's algebraic simplifier rewrites ``A / broadcast(B)`` into
    ``A * broadcast(1/B)`` (not IEEE-literal, <= 1.5 ulp); in numerics strict mode the
    broadcast divisor is materialised behind an optimization barrier to block that."""
    if strict() and jnp.shape(d) != jnp.shape(x):
        d = ob(jnp.broadcast_to(d, jnp.shape(x)))
    return x / d


def _smith_apply(a, t, r, den):
    """``a / b`` with the factors of :func:`_smith_factors`; bitwise equal to
    ``numerics.cdiv(a, b)`` (``q*r - p == -(p - q*r)`` exactly in IEEE arithmetic)."""
    ar, ai = a.real, a.imag
    p = jnp.where(t, ar, ai)
    q = jnp.where(t, ai, ar)
    re = _bdiv(ob(p * r) + q, den)
    d = p - ob(q * r)
    im = _bdiv(jnp.where(t, -d, d), den)
    return lax.complex(re, im)


def gauss_jordan(AMT, zonal_col, nf: int, nth: int, N: int):
    """``init_inverse_matrix`` (src/fields_implicit.f90:722-914), exact transcription (spec 04
    section 3.13): pivot-free Gauss-Jordan on ``AMT[m, j, k] = am(k+1, idx(j+1, m))``.

    For each pivot ``i`` not in the skip set (ascending), with ``L = AMT[:, i, :]`` and
    ``Rr = AIT[:, i, :]`` copied first (778-807), every row ``j`` not in the skip set of a
    non-zonal member (834) gets ``fac = AMT[j,i]/L[i]``; ``AMT[j,i] = fac``;
    ``AMT[j,k] -= L[k]*fac``; ``AIT[j,:] = AIT[i,:]/L[i]`` if ``j == i`` else
    ``AIT[j,:] - (Rr*L[j])/L[i]`` (836-853). Zonal members (``aky == 0 .and. akx == 0``)
    get ``AIT[j,:] = 0`` at every pivot (856). Then ``AIT[:, :, S] = 0`` (868-880) and
    ``AIT[:, j, :] = AIT[:, j+nf, :]`` for ``j`` in ``S`` (883-908).

    ``zonal_col``: ``(M, R)`` bool, the zonal flag of the cell of column ``j``.
    Returns ``AIT`` ``(M, R, R)`` in K order.
    """
    AMT = jnp.asarray(AMT, jnp.complex128)
    M, n, _ = AMT.shape
    skip = skip_set(N, nf, nth)
    piv = np.nonzero(~skip)[0].astype(np.int32)
    AIT = jnp.broadcast_to(jnp.eye(n, dtype=jnp.complex128), (M, n, n))
    nonskip = jnp.asarray(~skip)
    # GS2 quirk: members with aky == 0 .and. akx == 0 (exact) get zero inverse rows at every
    # pivot, so their fields never change -- also for such modes in non-box grids (Q-08)
    # (fields_implicit.f90:834-857)
    zc = jnp.asarray(zonal_col, bool)
    act = (nonskip[None, :] & ~zc)[:, :, None]
    zrow = (nonskip[None, :] & zc)[:, :, None]
    pivots = jnp.asarray(piv)
    czero = jnp.zeros((), jnp.complex128)

    def body(step, carry):
        A, B = carry
        i = pivots[step]
        L = lax.dynamic_index_in_dim(A, i, axis=1, keepdims=False)      # lhscol(:, m)
        Rr = lax.dynamic_index_in_dim(B, i, axis=1, keepdims=False)     # rhsrow(:, m)
        Lii = lax.dynamic_index_in_dim(L, i, axis=1, keepdims=False)    # lhscol(i, m)
        col = lax.dynamic_index_in_dim(A, i, axis=2, keepdims=False)    # am(i, jlo)
        t, r, den = _smith_factors(Lii)
        fac = _smith_apply(col, t[:, None], r[:, None], den[:, None])          # (M, n)
        An = A - cmul(L[:, None, :], fac[:, :, None])
        An = lax.dynamic_update_index_in_dim(An, fac, i, axis=2)
        Bo = B - _smith_apply(cmul(Rr[:, None, :], L[:, :, None]),
                              t[:, None, None], r[:, None, None], den[:, None, None])
        Bp = _smith_apply(Rr, t[:, None], r[:, None], den[:, None])
        Bn = lax.dynamic_update_index_in_dim(Bo, Bp, i, axis=1)
        A2 = jnp.where(act, An, A)
        B2 = jnp.where(act, Bn, jnp.where(zrow, czero, B))
        return A2, B2

    if piv.size:
        _, AIT = lax.fori_loop(0, int(piv.size), body, (AMT, AIT))
    AIT = jnp.where(jnp.asarray(skip)[None, None, :], czero, AIT)
    sk = np.nonzero(skip)[0]
    if sk.size:
        AIT = AIT.at[:, sk, :].set(AIT[:, sk + nf, :])
    return AIT


def zonal_columns(topo: Topology, c: int, nf: int, nth: int):
    """``(M, R)`` zonal flag of the cell of every K index of class ``c``."""
    M, N = topo.sizes(c)
    n_of_j = k_index(N, nf, nth)[0]
    return jnp.asarray(topo.zonal[c])[:, n_of_j]


def invert_response(AMT, topo: Topology, c: int, nf: int, nth: int):
    """``init_inverse_matrix`` for class ``c``: ``AMT`` (K order) -> ``aminv_c`` (B order)."""
    M, N = topo.sizes(c)
    AIT = gauss_jordan(AMT, zonal_columns(topo, c, nf, nth), nf, nth, N)
    return reorder_k_to_b(AIT, nf, nth, N)


def amt_from_export(ex: Export, c: int, *, layout: str = "K"):
    """``AMT`` ``(M_c, R_c, R_c)`` of class ``c`` (0-based) from the optional probe export
    ``/fields_implicit/am_c<c+1>`` (P-R1, spec 09 section 5.2), or None.

    The Fortran array is ``am(k, ilo)`` ``(R_c, R_c*M_c)``, NumPy ``[ilo, k]``; with
    ``ilo = j + R_c*m`` this reshapes to ``AMT[m, j, k]``. ``layout='B'`` if the exporter
    wrote the response index ``k`` in B order instead of GS2's K order.
    """
    am = _opt(ex, f"/fields_implicit/am_c{c + 1}")
    if am is None:
        return None
    am = np.asarray(am, np.complex128)
    R = am.shape[-1]
    AMT = am.reshape(-1, R, R)
    if layout == "B":
        nf = int(_scalar(_first(ex, "/fields_implicit/nfield"), int))
        nth = ex.nth
        N = R // (nf * nth)
        bofk = np.argsort(k_of_b(N, nf, nth))
        AMT = AMT[:, :, bofk]
    return AMT


# --------------------------------------------------------------------------- probes
def probe_table(N: int, nf: int, nth: int):
    """Probe columns ``j`` (K order) GS2 uses (fields_implicit.f90:558-633): all ``(n, ig, f)``
    except the skip set, whose responses are never read by init_inverse_matrix."""
    skip = skip_set(N, nf, nth)
    js = np.nonzero(~skip)[0].astype(np.int32)
    n, igp, f = k_index(N, nf, nth)
    return js, n[js].astype(np.int32), igp[js].astype(np.int32), f[js].astype(np.int32)


def probe_fields(n, igp, f, ikc, itc, shape, slots) -> Fields:
    """Unit field perturbation of probe ``(n, ig, f)`` for all members of a class at once
    (fields_implicit.f90:565-633): ``X_f(ig, it(m,n), ik(m)) = 1`` and, for ``n > 1`` (1-based)
    and ``ig == -ntgrid``, also ``X_f(ntgrid, it(m,n-1), ik(m)) = 1``; other fields 0."""
    Y, X, T = shape
    ikc = jnp.asarray(ikc)
    itc = jnp.asarray(itc)
    cell = jnp.zeros((Y, X), bool).at[ikc, itc[:, n]].set(True)
    prev = jnp.zeros((Y, X), bool).at[ikc, itc[:, jnp.maximum(n - 1, 0)]].set(True)
    t = jnp.arange(T)
    here = cell[..., None] & (t == igp)[None, None, :]
    left = (n > 0) & (igp == 0) & prev[..., None] & (t == T - 1)[None, None, :]
    unit = jnp.where(here | left, jnp.ones((), jnp.complex128), jnp.zeros((), jnp.complex128))
    zero = jnp.zeros((Y, X, T), jnp.complex128)
    out = [zero, zero, zero]
    for s, i in enumerate(slots):
        out[i] = jnp.where(f == s, unit, zero)
    return Fields(*out)


def _class_rows(fl, ikc, itc):
    """Residuals of one probe as ``(M, R)`` rows in K order (``init_response_row``,
    fields_implicit.f90:676-712: ``am(nidx*(nn-1) + f' + nf*igp', iflo) = fieldeq_f'``)."""
    flc = fl[:, jnp.asarray(ikc)[:, None], jnp.asarray(itc), :]     # (nf, M, N, T)
    M, N = flc.shape[1], flc.shape[2]
    return jnp.transpose(flc, (1, 2, 3, 0)).reshape(M, -1)


def _residual_fn(params: FieldParams, local: bool):
    if local:
        from .fields_local import getfieldeq_local
        return getfieldeq_local
    return getfieldeq


_NO_ARGS = object()


def response_rows(timeadv0: Callable[..., jnp.ndarray], data: FieldData,
                  params: FieldParams, c: int, *, batch_size: int = 1, local: bool = False,
                  timeadv_args=_NO_ARGS):
    """``AMT`` ``(M, R, R)`` of class ``c`` from unit probes (``init_response_row``,
    fields_implicit.f90:650-719): ``gnew = timeadv0(probe[, timeadv_args])`` (timeadv with
    ``istep = 0``, ``g = 0`` and old fields 0, fields_implicit.f90:553-557, 669), then
    ``getfieldeq(probe, gnew)`` (670). Skip-set rows are left 0 (never read)."""
    topo = data.topo
    M, N = topo.sizes(c)
    nf, nth = params.nfield, params.nth
    R = N * nf * nth
    js, ns, igps, fs = probe_table(N, nf, nth)
    ikc, itc = topo.ik[c], topo.it[c]
    shape = (params.naky, params.ntheta0, nth)
    residual = _residual_fn(params, local)

    def one(args):
        n, igp, f = args
        probe = probe_fields(n, igp, f, ikc, itc, shape, params.slots)
        g = timeadv0(probe) if timeadv_args is _NO_ARGS else timeadv0(probe, timeadv_args)
        ant = getan(g, data, params)
        fe = residual(probe, *ant, data, params)
        return _class_rows(stack_present(fe, params), ikc, itc)

    xs = (jnp.asarray(ns), jnp.asarray(igps), jnp.asarray(fs))
    rows = lax.map(one, xs, batch_size=max(1, int(batch_size)))          # (P, M, R)
    AMT = jnp.zeros((M, R, R), jnp.complex128)
    return AMT.at[:, jnp.asarray(js), :].set(jnp.transpose(rows, (1, 0, 2)))


# --------------------------------------------------------------------------- local grouping
def local_probe_index(f, p, N: int, nf: int, nth: int):
    """Call ``p`` (1-based) of field slot ``f`` in ``fm_populate`` for a supercell of ``N`` cells
    (src/fields_local.f90:2120-2167): the impulse goes to the first not yet initialised point
    of the extended domain, i.e. ``iex = p``, while ``p <= nextend = 1 + 2*ntgrid*N``.
    ``sc_iex_to_ig``/``sc_iex_to_it`` (838-855, 822-835): ``ig = ntgrid`` at ``iex = nextend``,
    else ``ig = -ntgrid + mod(iex-1, 2*ntgrid)`` in cell ``ic = 1 + (iex-1)/(2*ntgrid)``.

    Returns traced ``(active, n, igp, j)``: 0-based cell ``n``, theta index ``igp`` and the K
    index ``j = f + nf*(igp + nth*n)`` of that perturbation (``j`` = 0 when inactive)."""
    two_nt = nth - 1
    nextend = 1 + two_nt * N
    active = p <= nextend
    last = p == nextend
    n = jnp.where(last, N - 1, jnp.clip((p - 1) // two_nt, 0, N - 1))
    igp = jnp.where(last, nth - 1, (p - 1) % two_nt)
    j = f + nf * (igp + nth * n)
    return active, n, igp, jnp.where(active, j, 0)


def local_probe_fields(f, p, topo: Topology, params: FieldParams) -> Fields:
    """The fields of one ``fm_init_next_field_points`` call (src/fields_local.f90:2081-2167):
    ``field(ig, it, ik) = 1`` at extended point ``iex = p`` of **every** supercell of every ky
    that is not finished yet (empty ky = kx = 0 one included), plus the duplicate
    ``field(ntgrid, it_left, ik) = 1`` when ``ig = -ntgrid`` is not in the leftmost cell
    (2163-2167); finished supercells get nothing (they are only kwork_filter-ed, which does not
    change the other supercells' response). ``f`` is the field slot (phi, apar, bpar order)."""
    Y, X, T = params.naky, params.ntheta0, params.nth
    t = jnp.arange(T)
    unit = jnp.zeros((Y, X, T), bool)
    for c in range(topo.nclass):
        M, N = topo.sizes(c)
        active, n, igp, _ = local_probe_index(f, p, N, params.nfield, T)
        ikc = jnp.asarray(topo.ik[c])
        itc = jnp.asarray(topo.it[c])
        cell = jnp.zeros((Y, X), bool).at[ikc, itc[:, n]].set(True)
        prev = jnp.zeros((Y, X), bool).at[ikc, itc[:, jnp.maximum(n - 1, 0)]].set(True)
        here = cell[..., None] & (t == igp)[None, None, :]
        left = (n > 0) & (igp == 0) & prev[..., None] & (t == T - 1)[None, None, :]
        unit = unit | (active & (here | left))
    one = jnp.where(unit, jnp.ones((), jnp.complex128), jnp.zeros((), jnp.complex128))
    zero = jnp.zeros((Y, X, T), jnp.complex128)
    out = [zero, zero, zero]
    for s_, i in enumerate(params.slots):
        out[i] = jnp.where(f == s_, one, zero)
    return Fields(*out)


def response_rows_local(timeadv0: Callable[..., jnp.ndarray], data: FieldData,
                        params: FieldParams, *, batch_size: int = 1, timeadv_args=_NO_ARGS):
    """``AMT`` of every class (tuple of ``(M_c, R_c, R_c)``) with GS2's field_option='local'
    grouping: ``fm_populate`` (src/fields_local.f90:1973-2077) loops over the present fields
    (phi, apar, bpar) and calls ``fm_init_next_field_points`` until every point is done, i.e.
    ``p = 1 .. max(nextend)`` calls per field, each one ``timeadv(..., 0)`` (2184) with
    :func:`local_probe_fields` followed by ``getfieldeq_nogath`` (2189-2195, local
    association) and ``sc_store_fq`` of the residual as row ``iex = p`` of every active
    supercell (2199-2268, 1090-1165). The calls are independent (g and the old fields stay 0,
    the impulses are cleared after each call), so they are evaluated ``batch_size`` at a time.
    With amplitude-dependent hyperviscosity the simultaneous impulses of all supercells enter
    the shear rate exactly as in GS2 (spec 04 P9)."""
    from .fields_local import getfieldeq_local
    topo = data.topo
    nf, nth = params.nfield, params.nth
    if nth < 3:
        raise ValueError("ntgrid >= 1 required")
    C = topo.nclass
    sizes = [topo.sizes(c) for c in range(C)]
    pmax = max(1 + (nth - 1) * N for _, N in sizes)
    calls_f = np.repeat(np.arange(nf, dtype=np.int32), pmax)
    calls_p = np.tile(np.arange(1, pmax + 1, dtype=np.int32), nf)
    B = max(1, int(batch_size))
    pad = (-calls_f.size) % B
    calls_f = np.concatenate([calls_f, np.zeros(pad, np.int32)])
    calls_p = np.concatenate([calls_p, np.full(pad, pmax + 1, np.int32)])   # inactive
    nb = calls_f.size // B

    def one(f, p):
        probe = local_probe_fields(f, p, topo, params)
        g = timeadv0(probe) if timeadv_args is _NO_ARGS else timeadv0(probe, timeadv_args)
        ant = getan(g, data, params)
        fl = stack_present(getfieldeq_local(probe, *ant, data, params), params)
        return tuple(_class_rows(fl, topo.ik[c], topo.it[c]) for c in range(C))

    def body(AMTs, batch):
        fb, pb = batch
        rows = jax.vmap(one)(fb, pb)                          # tuple of (B, M_c, R_c)
        out = []
        for c in range(C):
            A = AMTs[c]
            M, N = sizes[c]
            for b in range(B):
                active, _, _, j = local_probe_index(fb[b], pb[b], N, nf, nth)
                cur = lax.dynamic_index_in_dim(A, j, axis=1, keepdims=False)
                A = lax.dynamic_update_index_in_dim(A, jnp.where(active, rows[c][b], cur),
                                                    j, axis=1)
            out.append(A)
        return tuple(out), None

    AMT0 = tuple(jnp.zeros((M, N * nf * nth, N * nf * nth), jnp.complex128) for M, N in sizes)
    xs = (jnp.asarray(calls_f.reshape(nb, B)), jnp.asarray(calls_p.reshape(nb, B)))
    AMTs, _ = lax.scan(body, AMT0, xs)
    return AMTs


def build_aminv(timeadv0: Callable[..., jnp.ndarray], data: FieldData, params: FieldParams,
                *, batch_size: int = 1, code_dt: float = math.nan, jit: bool = True,
                local: Optional[bool] = None, timeadv_args=_NO_ARGS) -> Aminv:
    """Rebuild the response matrices as GS2 does after a dt change.

    ``timeadv0(fields_new) -> gnew`` (or ``timeadv0(fields_new, timeadv_args)`` when
    ``timeadv_args`` is given; pass the dt-dependent coefficients that way so they are traced
    arguments, not closure constants) must run GS2's ``timeadv(phi=0, apar=0, bpar=0, phinew,
    aparnew, bparnew, istep=0)`` from ``g = 0`` with the current coefficients (no explicit/NL
    terms, hyper iff hyper_in_initialisation, collisions iff not split_collisions with the
    probe ``aparnew`` in the drag (G-27), parity, linked BCs; spec 10 section C).
    ``batch_size`` probes are vmapped per ``timeadv0`` evaluation (memory ~ batch_size*G7).

    * ``local=False`` (field_option 'implicit'): ``init_response_matrix`` /
      ``init_response_row`` (src/fields_implicit.f90:459-719): every probe perturbs all members
      of one class at one ``(n, ig, f)``, exactly GS2's pattern (G-18); residuals from
      getfieldeq1. Classes are processed independently (GS2 runs ``c = i_class..1``).
    * ``local=True`` (default when ``params.field_option == 'local'``): ``fm_populate``
      (src/fields_local.f90:1973-2077, called at 4133 on every (re)initialisation): one impulse
      in every unfinished supercell of every ky per ``timeadv`` call
      (:func:`response_rows_local`), residuals with the fields_local association
      (``fm_getfieldeq1_nogath``). The inversion is the implicit Gauss-Jordan instead of
      ``invert_block`` (round-off level, 00_decisions D5).

    The inversion is ``init_inverse_matrix`` (:func:`gauss_jordan`) in both cases. The current
    ``data.kwork_filter`` is used (``(1,1) = .true.`` in box runs; it only affects the zonal
    member whose inverse rows are zeroed either way, spec 10 G-17; fm_populate's filtering of
    finished supercells does not change the other supercells' rows). ``gnew`` of the evolving
    state is not touched (G-28).
    """
    if local is None:
        local = params.field_option == "local"
    nf, nth = params.nfield, params.nth
    if local:
        def run_local(d, targs):
            AMTs = response_rows_local(timeadv0, d, params, batch_size=batch_size,
                                       timeadv_args=targs)
            return tuple(invert_response(A, d.topo, c, nf, nth) for c, A in enumerate(AMTs))
        if timeadv_args is _NO_ARGS:
            run = (lambda d: run_local(d, _NO_ARGS))
            mats = (jax.jit(run) if jit else run)(data)
        else:
            mats = (jax.jit(run_local) if jit else run_local)(data, timeadv_args)
        return Aminv(mats=tuple(mats), code_dt=float(code_dt))
    mats = []
    for c in range(data.topo.nclass):
        def per_class(d, targs, c=c):
            AMT = response_rows(timeadv0, d, params, c, batch_size=batch_size, local=False,
                                timeadv_args=targs)
            return invert_response(AMT, d.topo, c, nf, nth)
        if timeadv_args is _NO_ARGS:
            run = (lambda d, f=per_class: f(d, _NO_ARGS))
            mats.append((jax.jit(run) if jit else run)(data))
        else:
            mats.append((jax.jit(per_class) if jit else per_class)(data, timeadv_args))
    return Aminv(mats=tuple(mats), code_dt=float(code_dt))


def condition_numbers(mat, topo: Topology, c: int, nf: int, nth: int) -> np.ndarray:
    """2-norm condition number of every member's reduced system (inputs/outputs outside the
    skip set; zonal members give NaN), from ``aminv_c`` (B order). Diagnostic only."""
    M, N = topo.sizes(c)
    skip = skip_set(N, nf, nth)
    kofb = k_of_b(N, nf, nth)
    keep_b = ~skip[kofb]
    out = np.full(M, np.nan)
    mat = np.asarray(mat)
    zon = np.asarray(topo.zonal[c]).any(axis=1)
    for m in range(M):
        if zon[m]:
            continue
        sub = mat[m][np.ix_(keep_b, keep_b)]
        out[m] = np.linalg.cond(sub)
    return out
