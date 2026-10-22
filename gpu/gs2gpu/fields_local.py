"""``field_option = 'local'`` field solve (src/fields_local.f90; spec 04 section 7, spec 10 G-23,
00_decisions D5).

On one MPI rank (decomp_type 3, ``pc_decomp_owncells_simplempi``, fields_local.f90:3771-3927)
every cell of a supercell owns all rows of the supercell's matrix, so the update is

* residuals from ``fm_getfieldeq1_nogath`` (fields_local.f90:3374-3443; ``has_to_gather =
  .false.``, 3788), whose phi equation is ``((antot - gamtot*gridfac1*phi) + fl_avg) +
  bpar*gamtot1`` (a different association from getfieldeq1);
* per cell ``c`` and present field ``f``: ``S_{c,f} = matmul(fq_f(:ncol), rb(f)%data)``
  (rb_mv_mult_fun, 370-384), ``T_c = ((0 - S_phi) - S_apar) - S_bpar`` (c_get_field_update,
  512-546), and ``u = ((0 + T_1) + T_2) + ...`` over non-empty cells (sc_reduce_tmpsum,
  723-774). The matmul sum is sequential in ``k`` from 0 (spec 04 Q3, measured): GS2's
  reference build reproduces bitwise with ``S = ((0 + p_1) + p_2) + ...`` and FMA-fused
  products ``p.re = fma(-a.im, x.im, a.re*x.re)``, ``p.im = fma(a.im, x.re, a.re*x.im)``
  (``a`` = matrix, ``x`` = vect; the precompiled libgfortran ``matmul_c8`` is not affected by
  ``-ffp-contract=off``). JAX has no FMA primitive, so plain products are used here (<= 1 ulp
  per product; measured <= 7e-17 relative on ``u`` in strict mode);
* duplicate boundary outputs filled by copy and empty (``ky = kx = 0``) supercells left at 0
  (fm_gather_fields, 2593-2803); ``phinew = phinew + phi`` for present fields only (4316-4320).

The exported local inverse (``/fields_local/resp_ik<ik>_is<is>``, the rowblock data after
``sc_invert_local`` written by ``sc_dump_to_file``, fields_local.f90:987-1017) is converted to the
class/B layout of :class:`gs2gpu.response.Aminv`: column ``B(ig, f, n)`` of a non-last cell's
``ig = ntgrid`` gets a zero coefficient (that input is not part of the local matvec; appending
``0*x`` at the end of a sequential sum leaves it unchanged), and the duplicate output rows are
copies of their partners. After a dt change, ``response.build_aminv`` (implicit algorithm,
local residual association) provides the matrix (D5: round-off-level deviation from GS2's
``invert_block`` inverse).
"""

from __future__ import annotations

import math
import re
from typing import Optional

import jax.numpy as jnp
import numpy as np

from .fields import (FieldData, FieldParams, Topology, _first, _fields_from_stack, _fl_avg,
                     _opt, _scalar, _seq_prod_sum, _zeros_field, gather_cells, getan,
                     scatter_cells, stack_present, topology_from_export)
from .io.export import Export
from .numerics import rmul
from .response import Aminv
from .types import Fields

_RESP = re.compile(r"^resp_ik(\d+)_is(\d+)$")


def local_supercell_order(itleftmost, naky: int, ntheta0: int):
    """``is`` numbering of ``fm_init`` (fields_local.f90:1880-1969): per ``ik``, supercells in
    order of first appearance scanning ``it = 1..ntheta0``; returns ``{(ik, is): it_leftmost}``
    (all 1-based). ``itleftmost`` is ``get_leftmost_it(it, ik)`` NumPy ``(Y, X)``, 1-based."""
    out = {}
    itleftmost = np.asarray(itleftmost)
    for ik in range(naky):
        seen = []
        for it in range(ntheta0):
            lm = int(itleftmost[ik, it])
            if lm not in seen:
                seen.append(lm)
        for is_, lm in enumerate(seen):
            out[(ik + 1, is_ + 1)] = lm
    return out


def _leftmost_from_export(ex: Export, topo: Topology):
    lm = _opt(ex, "/dist_fn/itleftmost")
    if lm is not None:
        return np.asarray(lm)
    Y, X = np.shape(topo.perm)
    out = np.zeros((Y, X), np.int32)
    for ikc, itc in zip(topo.ik, topo.it):
        for m in range(itc.shape[0]):
            out[ikc[m], itc[m]] = itc[m, 0] + 1
    return out


def ext_index_of_b(N: int, nf: int, ntgrid: int) -> tuple[np.ndarray, np.ndarray]:
    """For every B index ``b = igp + nth*(f + nf*n)`` of a supercell with ``N`` cells: the
    extended-domain index ``iex - 1 = igp + 2*ntgrid*n + nextend*f`` (sc_store_fq,
    fields_local.f90:1112; ``nextend = 1 + 2*ntgrid*N``) and whether ``b`` is a duplicate
    point (``ig = ntgrid`` of a non-last cell)."""
    nth = 2 * ntgrid + 1
    nextend = 1 + 2 * ntgrid * N
    b = np.arange(N * nf * nth)
    igp = b % nth
    f = (b // nth) % nf
    n = b // (nf * nth)
    e = igp + 2 * ntgrid * n + nextend * f
    dup = (igp == nth - 1) & (n < N - 1)
    return e.astype(np.int64), dup


def local_to_b(resp, N: int, nf: int, ntgrid: int) -> np.ndarray:
    """One supercell's local matrix ``resp[row, col]`` (extended order) -> ``(R, R)`` B layout."""
    e, dup = ext_index_of_b(N, nf, ntgrid)
    resp = np.asarray(resp, np.complex128)
    nextend = 1 + 2 * ntgrid * N
    if resp.shape != (nf * nextend, nf * nextend):
        raise ValueError(f"local response shape {resp.shape}, expected {(nf * nextend,) * 2}")
    mat = resp[np.ix_(e, e)].copy()
    mat[:, dup] = 0.0
    return mat


def b_to_local(mat, N: int, nf: int, ntgrid: int) -> np.ndarray:
    """Inverse of :func:`local_to_b` (drops the duplicate rows/columns)."""
    e, dup = ext_index_of_b(N, nf, ntgrid)
    nextend = 1 + 2 * ntgrid * N
    out = np.zeros((nf * nextend, nf * nextend), np.complex128)
    keep = ~dup
    out[np.ix_(e[keep], e[keep])] = np.asarray(mat)[np.ix_(keep, keep)]
    return out


def local_aminv_from_export(ex: Export, rung: Optional[int] = None,
                            topo: Optional[Topology] = None, transpose: bool = False) -> Aminv:
    """``/fields_local/resp_ik<ik>_is<is>`` (spec 09 section 7.5) -> :class:`Aminv` in the
    class/member order of ``topo`` and B layout.

    The netCDF ``response(ri, ax1, ax2)`` holds the Fortran array ``tmp_arr(col, row)``
    pulled from the rowblocks (``arr(cl:cu, rl:ru) = rb%data``, 1330-1335), so the NumPy array
    is ``[row (output), col (input)]``. Missing supercells (the ky = kx = 0 one is never dumped,
    1047-1057) stay 0. ``transpose=True`` swaps the interpretation (diagnostics).
    """
    if topo is None:
        topo = topology_from_export(ex)
    base = "/fields_local" if rung is None else f"/dt_ladder/k{rung}/fields_local"
    grp = ex[base]
    nf = int(_scalar(_first(ex, "/fields_implicit/nfield"), int) or
             sum(bool(np.asarray(ex[f"/run_parameters/has_{n}"]).reshape(-1)[0])
                 for n in ("phi", "apar", "bpar")))
    ntgrid = ex.ntgrid
    nth = ex.nth
    leftmost = _leftmost_from_export(ex, topo)
    order = local_supercell_order(leftmost, ex.naky, ex.ntheta0)
    where = {}
    for c, (ikc, itc) in enumerate(zip(topo.ik, topo.it)):
        for m in range(itc.shape[0]):
            where[(int(ikc[m]) + 1, int(itc[m, 0]) + 1)] = (c, m)
    mats = [np.zeros((itc.shape[0], itc.shape[1] * nf * nth, itc.shape[1] * nf * nth),
                     np.complex128) for itc in topo.it]
    for name in getattr(grp, "arrays", {}):
        mt = _RESP.match(name)
        if not mt:
            continue
        ik, is_ = int(mt.group(1)), int(mt.group(2))
        lm = _first(ex, f"{base}/it_leftmost_ik{ik}_is{is_}")
        if lm is None:
            lm = grp.meta.get(name, {}).get("it_leftmost")
        lm = int(_scalar(lm, int)) if lm is not None else order[(ik, is_)]
        c, m = where[(ik, lm)]
        N = topo.it[c].shape[1]
        resp = np.asarray(grp[name], np.complex128)
        if transpose:
            resp = resp.T
        mats[c][m] = local_to_b(resp, N, nf, ntgrid)
    if rung is None:
        code_dt = _scalar(_first(ex, "/state_t0/code_dt", "/gs2_time/code_dt"), float)
    else:
        code_dt = _scalar(ex[f"/dt_ladder/k{rung}"].attrs.get("code_dt"), float)
    return Aminv(mats=tuple(mats), code_dt=math.nan if code_dt is None else float(code_dt))


def getfieldeq_local(fields_new: Fields, antot, antota, antotp, data: FieldData,
                     params: FieldParams, kwork_filter=None):
    """``fm_getfieldeq1_nogath`` (src/fields_local.f90:3374-3443).

    ``fieldeq  = ((antot - (gamtot*gridfac1)*phi) + fl_avg) [+ bpar*gamtot1 if has_bpar]``
    ``fieldeqa = antota - (kperp2*gridfac1)*apar``
    ``fieldeqp = bpar*gridfac1 + ((antotp + bpar*gamtot2) + (0.5*phi)*gamtot1)*(beta*apfac/bmag**2)``

    GS2 skips empty and kwork_filter cells (their values are stale and never used); here they
    are computed like the others.
    """
    phi, apar, bpar = fields_new
    fl = _fl_avg(antot, data, params, kwork_filter)
    fe = fea = fep = _zeros_field(antot.shape)
    if params.has_phi:
        fe = (antot - rmul(phi, data.gamtot_gf)) + fl[..., None]
        if params.has_bpar:
            fe = fe + rmul(bpar, data.gamtot1)
    if params.has_apar:
        fea = antota - rmul(apar, data.kperp2_gf)
    if params.has_bpar:
        fep = rmul(bpar, data.gridfac1) + rmul(
            (antotp + rmul(bpar, data.gamtot2)) + rmul(rmul(phi, 0.5), data.gamtot1), data.bfac)
    return fe, fea, fep


def apply_local(fl, mats, topo: Topology, params: FieldParams):
    """Local matvec (fields_local.f90:370-384, 512-546, 723-774) on B-layout matrices.

    ``fl``: ``(nf, Y, X, T)`` residuals; returns ``(nf, Y, X, T)``.
    """
    nf, nth = params.nfield, params.nth
    outs = []
    for c, mat in enumerate(mats):
        M, N = topo.sizes(c)
        R = N * nf * nth
        flc = gather_cells(fl, topo, c).reshape(M, N, nf, nth)
        mat5 = jnp.reshape(mat, (M, R, N, nf, nth))
        emp = jnp.asarray(topo.empty[c])                    # (M, N)
        czero = jnp.zeros((), jnp.complex128)
        if params.exact:
            S = _seq_prod_sum(mat5, flc[:, None])            # (M, R, N, nf): matmul per rowblock
            Tc = jnp.zeros((M, R, N), jnp.complex128)
            for f in range(nf):                              # tmp_sum = 0; tmp_sum -= rb(f)%tmp_sum
                Tc = Tc - S[..., f]
            u = jnp.zeros((M, R), jnp.complex128)
            for n in range(N):                               # local_tmp_sum += cell%tmp_sum
                u = u + jnp.where(emp[:, None, n], czero, Tc[:, :, n])
        else:
            Tc = -jnp.einsum("mrnfk,mnfk->mrn", mat5, flc)
            u = jnp.sum(jnp.where(emp[:, None, :], czero, Tc), axis=-1)
        outs.append(u)
    return scatter_cells(outs, topo, nf, nth)


def getfield_local(fields_new: Fields, gnew, aminv, data: FieldData,
                   params: FieldParams) -> Fields:
    """``getfield_local(phinew, aparnew, bparnew)`` (src/fields_local.f90:4193-4230) with
    ``do_gather = .true.``, ``do_update = .false.`` (single rank, 4270-4272).

    Returns the increment for the present fields (``fm_gather_fields`` copies only those,
    2771-2784); absent fields are returned unchanged from ``fields_new``.
    """
    ant = getan(gnew, data, params)                 # getan_nogath == getan on one rank
    fe = getfieldeq_local(fields_new, *ant, data, params)
    fl = stack_present(fe, params)
    mats = aminv.mats if hasattr(aminv, "mats") else aminv
    u = apply_local(fl, mats, data.topo, params)
    return _fields_from_stack(u, fields_new, params)


def field_update_local(u: Fields, fields_old: Fields, params: FieldParams) -> Fields:
    """``if (has_phi) phinew = phinew + phi`` etc. (src/fields_local.f90:4316-4320)."""
    out = list(u)
    for i in params.slots:
        out[i] = u[i] + fields_old[i]
    return Fields(*out)
