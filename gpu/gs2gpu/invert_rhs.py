"""The implicit parallel solve: ``invert_rhs_1`` and ``init_homogeneous_g`` (spec 03).

Mirrors, in src/dist_fn.fpp (GS2 8.2.1, commit 4d8c94bcf line numbers, as in
the specs):

* ``invert_rhs_1`` (5440-5707) with ``apply_parallel_boundary_conditions``
  (5710-5849, the ``nonad_zero`` branch) and ``self_periodic`` (5864-5906);
* ``init_homogeneous_g`` (1937-2125);
* ``invert_rhs`` (6152-6193) minus the source (owned by ``gs2gpu.source``);
* ``set_h_zero`` and ``get_adjust`` (src/dist_fn_arrays.fpp:246-294, 178-188).

The linked boundary conditions live in ``gs2gpu.linked`` and the parity
projection in ``gs2gpu.parity``; they share :class:`InvertData` and
:class:`InvertParams`.

Layout: arrays enter and leave in the canonical G7 layout ``(S, E, L, Y, X,
2, T)``. Inside, each sign is moved to theta-first ``(T, S, E, L, Y, X)`` so
that the theta recurrences are ``lax.scan`` loops batched over every line.

Round-off: every expression keeps the Fortran association order; mixed
real/complex products are component-wise (``numerics.rmul``), complex
division is Smith's method (``numerics.cdiv``). Fortran ``-a*b + c`` is
written ``c - a*b`` (bitwise identical in IEEE arithmetic).

Supported ``boundary_option``: default (resolved by GS2), zero, unconnected,
alternate-zero (exact alias of zero, G-13), self-periodic, periodic,
kperiod=1, linked (box only). Gated (``NotImplementedError``):
``start_from_previous_solution``, ``boundary_off_grid``,
``exponential_boundary``, linked boundaries without a box grid.
"""

from __future__ import annotations

from typing import NamedTuple, Optional

import jax.numpy as jnp
import numpy as np
from jax import lax

from . import linked as _linked
from . import parity as _parity
from .numerics import cdiv, cmul, rmul, rsub

# boundary_option_switch values (src/dist_fn.fpp:121-125)
BOUNDARY_DEFAULT = 0
BOUNDARY_ZERO = 1
BOUNDARY_SELF_PERIODIC = 2
BOUNDARY_ALTERNATE_ZERO = 3
BOUNDARY_LINKED = 4

HYBRID_ELECTRON_SPECIES = 4   # src/species.f90:40


# ============================================================================ helpers (host)
def _opt_array(group, name, default=None):
    """``group[name]`` as NumPy if present (allocated), else ``default``."""
    if name in group:
        return np.asarray(group[name])
    return default


def _flag(value) -> bool:
    return bool(np.asarray(value).item()) if value is not None else False


def _scalar(group, name, default=None):
    if name in group:
        return np.asarray(group[name]).item()
    if default is None:
        raise KeyError(f"{getattr(group, 'path', '?')}/{name} missing from the export")
    return default


def lower_bounce_points(is_bounce_point, forbid, trapped_wfb: bool, ng2: int):
    """``is_lower_bounce_point(ig, il)`` for every point (src/le_grids.f90:2775-2819)."""
    bp = np.asarray(is_bounce_point, bool)
    fb = np.asarray(forbid, bool)
    nlambda, nth = bp.shape
    out = np.zeros_like(bp)
    wfb_block = (np.arange(1, nlambda + 1) == ng2 + 1) & (not trapped_wfb)
    out[:, 0] = bp[:, 0]
    out[:, 1:-1] = bp[:, 1:-1] & (fb[:, :-2] | (fb[:, :-2] & fb[:, 2:]))
    out[:, -1] = False
    out[wfb_block] = False
    return out


def upper_bounce_points(is_bounce_point, forbid, trapped_wfb: bool, ng2: int):
    """``is_upper_bounce_point(ig, il)`` for every point (src/le_grids.f90:2823-2867)."""
    bp = np.asarray(is_bounce_point, bool)
    fb = np.asarray(forbid, bool)
    nlambda, nth = bp.shape
    out = np.zeros_like(bp)
    wfb_block = (np.arange(1, nlambda + 1) == ng2 + 1) & (not trapped_wfb)
    out[:, 0] = False
    out[:, 1:-1] = bp[:, 1:-1] & (fb[:, 2:] | (fb[:, :-2] & fb[:, 2:]))
    out[:, -1] = bp[:, -1]
    out[wfb_block] = False
    return out


def trapped_beta_tables(forbid, is_ttp, is_trapped):
    """Static form of the trapped homogeneous correction loop (dist_fn.fpp:5658-5671).

    Simulates on the masks only::

        beta1 = 0
        do ig = ntgrid-1, -ntgrid, -1
           if (forbid(ig) .or. is_ttp(ig)) then: beta1 = 0; cycle
           else if (forbid(ig+1)) then: beta1 = (gnew(ig,1) - gnew(ig,2))/(1 - g_h(ig,1))
           gnew(ig,:) += beta1*g_h(ig,:)

    Returns ``top`` (theta index whose beta applies, -1 for beta = 0) and
    ``upd`` (point updated), both ``(L, T)``. Every point is updated once, by
    its own well, and a top is updated only after its beta is formed, so
    evaluating every beta from the pre-update ``gnew`` is bitwise identical.
    """
    forbid = np.asarray(forbid, bool)
    is_ttp = np.asarray(is_ttp, bool)
    nlambda, nth = forbid.shape
    top = np.full((nlambda, nth), -1, np.int32)
    upd = np.zeros((nlambda, nth), bool)
    for il in range(nlambda):
        if not is_trapped[il]:
            continue
        cur = -1
        # GS2 quirk: the loop starts at ig = ntgrid-1, so ig = ntgrid is never corrected
        # and an allowed region touching ig = ntgrid gets beta = 0 (dist_fn.fpp:5660; 03 §6.9)
        for t in range(nth - 2, -1, -1):
            if forbid[il, t] or is_ttp[il, t]:
                cur = -1
                continue
            if forbid[il, t + 1]:
                cur = t
            top[il, t] = cur
            upd[il, t] = True
    return top, upd


# ============================================================================ containers
class InvertParams(NamedTuple):
    """Static switches and sizes of the parallel solve (hashable; jit-static)."""

    ntgrid: int
    nspec: int
    negrid: int
    nlambda: int
    naky: int
    ntheta0: int
    ng2: int
    boundary_option_switch: int
    linked: bool                 # boundary_option_switch == linked
    nonad_zero: bool
    zero_forbid: bool
    def_parity: bool
    even: bool
    passing_wfb: bool
    trapped_wfb: bool
    mixed_wfb: bool
    has_wfb: bool                # nlambda > ng2 (grid_has_trapped_particles)
    has_trapped: bool            # some il > ng2+1
    has_ttp: bool                # any(is_ttp)
    any_hybrid: bool             # any(is_passing_hybrid_electron)
    any_passthrough: bool        # any(kwork_filter .or. ieqzip)
    any_self_periodic: bool
    no_connections: bool
    n_open_l: int                # pitch angles with the passing open-chain linked fix
    fphi: float
    fbpar: float
    wfb: float

    @classmethod
    def from_export(cls, ex) -> "InvertParams":
        return from_export(ex)[1]


class InvertData(NamedTuple):
    """Arrays of the parallel solve (NumPy from ``from_export``; ``to_jax`` before use).

    ``vperp2, aj0, aj1``: G6; ``zt`` ``(S,)``; ``nonmaxw_corr`` ``(S, E)``;
    ``kwork_filter, ieqzip`` ``(Y, X)`` bool; ``hybrid`` ``(S, Y, L)``
    (``is_passing_hybrid_electron``); ``forbid, ttp, ttp_copy, lbsrc`` ``(L,
    T)`` bool (``lbsrc[l, ig]``: the source at ``ig`` takes the lower-bounce
    value, ``.not. is_ttp(ig) .and. is_lower_bounce_point(ig+1)``, trapped
    only); ``beta_top`` (clipped to >= 0), ``beta_has_top``, ``beta_upd``
    ``(L, T)`` (:func:`trapped_beta_tables`); ``self_periodic, bc_left,
    bc_right`` ``(L, Y, X)`` bool; ``sp_dadj, wfb_trapped`` ``(L,)`` bool;
    for ``init_homogeneous_g``: ``gh_skip`` ``(L, Y)`` bool (``il < ilmin``),
    ``gh_bc1, gh_bc2`` ``(L, Y)`` real incoming values, ``gh_g2up`` ``(L, T)``
    real (1 at trapped non-ttp upper bounce points); ``linked``:
    ``linked.LinkedData`` or None; ``parity``: ``parity.ParityData`` or None.
    """

    vperp2: np.ndarray
    aj0: np.ndarray
    aj1: np.ndarray
    zt: np.ndarray
    nonmaxw_corr: np.ndarray
    kwork_filter: np.ndarray
    ieqzip: np.ndarray
    hybrid: np.ndarray
    forbid: np.ndarray
    ttp: np.ndarray
    ttp_copy: np.ndarray
    lbsrc: np.ndarray
    beta_top: np.ndarray
    beta_has_top: np.ndarray
    beta_upd: np.ndarray
    self_periodic: np.ndarray
    bc_left: np.ndarray
    bc_right: np.ndarray
    sp_dadj: np.ndarray
    wfb_trapped: np.ndarray
    gh_skip: np.ndarray
    gh_bc1: np.ndarray
    gh_bc2: np.ndarray
    gh_g2up: np.ndarray
    linked: Optional[_linked.LinkedData]
    parity: Optional[_parity.ParityData]

    @classmethod
    def from_export(cls, ex) -> "InvertData":
        return from_export(ex)[0]


# ============================================================================ host construction
def build(*, ntgrid: int, ng2: int, aky, kwork_filter, ieqzip, forbid, is_ttp,
          is_lower_bounce_point, is_upper_bounce_point, is_passing_hybrid_electron,
          passing_wfb: bool, trapped_wfb: bool, mixed_wfb: bool,
          boundary_option_switch: int, nonad_zero: bool, zero_forbid: bool,
          def_parity: bool, even: bool, wfb: float, fphi: float, fbpar: float,
          vperp2, aj0, aj1, zt, nonmaxw_corr, l_links=None, r_links=None,
          itleft=None, itright=None, no_connections=None, wfb_factors=None,
          box: bool = True, start_from_previous_solution: bool = False,
          boundary_off_grid: bool = False, exponential_boundary: bool = False
          ) -> tuple[InvertData, InvertParams]:
    """Host-side construction of the static masks (flags of dist_fn.fpp:5466-5558).

    All inputs are NumPy in the canonical layouts (spec 11 section 3); index
    values in ``itleft/itright`` are Fortran (1-based, -1 = none).
    """
    if start_from_previous_solution:
        raise NotImplementedError("dist_fn_knobs: start_from_previous_solution is gated (03 §3.12.1)")
    if boundary_off_grid:
        raise NotImplementedError("dist_fn_knobs: boundary_off_grid is gated (03 §3.8)")
    if exponential_boundary:
        raise NotImplementedError("dist_fn_knobs: exponential_boundary is gated (03 §3.8)")
    bos = int(boundary_option_switch)
    if bos not in (BOUNDARY_ZERO, BOUNDARY_SELF_PERIODIC, BOUNDARY_ALTERNATE_ZERO,
                   BOUNDARY_LINKED):
        raise NotImplementedError(f"boundary_option_switch = {bos} (unresolved default?)")
    linked = bos == BOUNDARY_LINKED
    if linked and not box:
        raise NotImplementedError("boundary_option = 'linked' with grid_option /= 'box' is gated (03 §4)")

    aky = np.asarray(aky, np.float64)
    kwork_filter = np.asarray(kwork_filter, bool)
    ieqzip = np.asarray(ieqzip, bool)
    forbid = np.asarray(forbid, bool)
    is_ttp = np.asarray(is_ttp, bool)
    lbp = np.asarray(is_lower_bounce_point, bool)
    ubp = np.asarray(is_upper_bounce_point, bool)
    hybrid = np.asarray(is_passing_hybrid_electron, bool)
    nlambda, nth = forbid.shape
    assert nth == 2 * ntgrid + 1
    naky, ntheta0 = kwork_filter.shape
    nspec, negrid = np.asarray(nonmaxw_corr).shape
    if l_links is None:
        l_links = np.zeros((naky, ntheta0), np.int32)
    if r_links is None:
        r_links = np.zeros((naky, ntheta0), np.int32)
    l_links = np.asarray(l_links)
    r_links = np.asarray(r_links)

    il = np.arange(1, nlambda + 1)
    is_passing = il <= ng2                                   # il_is_passing
    is_wfb = il == ng2 + 1                                   # il_is_wfb
    is_trapped = il > ng2 + 1                                # il_is_trapped
    passing_like = is_passing | (is_wfb & passing_wfb)
    aky0 = aky == 0.0

    # ---- flags (dist_fn.fpp:5503-5558)
    use_pass_homog = (bos in (BOUNDARY_SELF_PERIODIC, BOUNDARY_LINKED)) | aky0     # (Y,)
    shape_lyx = (nlambda, naky, ntheta0)
    if linked:
        itleft = np.asarray(itleft)
        itright = np.asarray(itright)
        excluded = is_trapped | (is_wfb & trapped_wfb)
        neighbor = (~excluded)[:, None, None] & ((itleft >= 0) | (itright >= 0))[None]
        speriod = aky0[None, :, None]
        self_periodic = ((speriod & is_passing[:, None, None])
                         | (speriod & (is_wfb & (not trapped_wfb))[:, None, None])
                         | ((is_wfb & mixed_wfb)[:, None, None] & ~neighbor))
    else:
        self_periodic = ((use_pass_homog[None, :, None] & passing_like[:, None, None])
                         | (is_wfb & mixed_wfb)[:, None, None])
    self_periodic = np.broadcast_to(self_periodic, shape_lyx).copy()
    apply_bc = passing_like[:, None, None] & ~self_periodic
    bc_left = apply_bc & (l_links == 0)[None]
    bc_right = apply_bc & (r_links == 0)[None]

    # ---- trapped particles (dist_fn.fpp:5604-5622, 5634-5645, 5658-5671)
    lbsrc = np.zeros((nlambda, nth), bool)
    lbsrc[:, :-1] = is_trapped[:, None] & ~is_ttp[:, :-1] & lbp[:, 1:]
    top, upd = trapped_beta_tables(forbid, is_ttp, is_trapped)
    ttp_copy = is_ttp & ~forbid

    # ---- init_homogeneous_g (dist_fn.fpp:1959-2060)
    ilmin_skip = (~use_pass_homog)[None, :] & (il < ng2 + 1)[:, None]               # (L, Y)
    gh_bc1 = np.zeros((nlambda, naky))
    gh_bc2 = np.zeros((nlambda, naky))
    pass_homog = is_passing[:, None] & use_pass_homog[None, :]
    gh_bc1[pass_homog] = 1.0
    gh_bc2[pass_homog] = 1.0
    if mixed_wfb:
        gh_bc1[is_wfb, :] = wfb
        gh_bc2[is_wfb, :] = wfb
    elif trapped_wfb:
        gh_bc2[is_wfb, :] = 1.0
    elif passing_wfb:
        gh_bc1[is_wfb, :] = 1.0
        gh_bc2[is_wfb, :] = 1.0
    gh_g2up = (is_trapped[:, None] & ~is_ttp & ubp).astype(np.float64)

    # ---- linked topology and parity map
    linked_data = None
    nmax = 0
    if linked:
        linked_data, nmax = _linked.build_linked(
            itleft=itleft, itright=itright, l_links=l_links, r_links=r_links, aky=aky,
            kwork_filter=kwork_filter, ng2=ng2, nlambda=nlambda, passing_wfb=passing_wfb,
            trapped_wfb=trapped_wfb, mixed_wfb=mixed_wfb, nonad_zero=nonad_zero,
            hybrid=hybrid, zt=zt, nonmaxw_corr=nonmaxw_corr, wfb_factors=wfb_factors)
        # GS2 sets no_connections from the passing-particle links only
        # (setup_connected_bc_redistribute, dist_fn.fpp:2870-2876); the exported value wins.
        if no_connections is None:
            no_connections = linked_data is None
        elif not bool(no_connections) and linked_data is None:
            raise ValueError("exported no_connections = F but itright/l_links/r_links "
                             "contain no chain with two or more cells")
    no_connections = bool(no_connections) if linked else True
    parity_data = None
    if def_parity:
        parity_data = _parity.build_parity(l_links, r_links, itleft, itright, linked)

    has_wfb = nlambda > ng2
    data = InvertData(
        vperp2=np.asarray(vperp2, np.float64), aj0=np.asarray(aj0, np.float64),
        aj1=np.asarray(aj1, np.float64), zt=np.asarray(zt, np.float64),
        nonmaxw_corr=np.asarray(nonmaxw_corr, np.float64),
        kwork_filter=kwork_filter, ieqzip=ieqzip, hybrid=hybrid,
        forbid=forbid, ttp=is_ttp, ttp_copy=ttp_copy, lbsrc=lbsrc,
        beta_top=np.maximum(top, 0).astype(np.int32), beta_has_top=top >= 0, beta_upd=upd,
        self_periodic=self_periodic, bc_left=bc_left, bc_right=bc_right,
        sp_dadj=is_wfb & bool(nonad_zero) & (not passing_wfb),
        wfb_trapped=is_wfb & bool(trapped_wfb),
        gh_skip=ilmin_skip, gh_bc1=gh_bc1, gh_bc2=gh_bc2, gh_g2up=gh_g2up,
        linked=linked_data, parity=parity_data)
    params = InvertParams(
        ntgrid=int(ntgrid), nspec=int(nspec), negrid=int(negrid), nlambda=int(nlambda),
        naky=int(naky), ntheta0=int(ntheta0), ng2=int(ng2), boundary_option_switch=bos,
        linked=linked, nonad_zero=bool(nonad_zero), zero_forbid=bool(zero_forbid),
        def_parity=bool(def_parity), even=bool(even), passing_wfb=bool(passing_wfb),
        trapped_wfb=bool(trapped_wfb), mixed_wfb=bool(mixed_wfb), has_wfb=bool(has_wfb),
        has_trapped=bool(is_trapped.any()), has_ttp=bool(is_ttp.any()),
        any_hybrid=bool(hybrid.any()), any_passthrough=bool((kwork_filter | ieqzip).any()),
        any_self_periodic=bool(self_periodic.any()), no_connections=no_connections,
        n_open_l=int(ng2 + (1 if (passing_wfb and has_wfb) else 0)),
        fphi=float(fphi), fbpar=float(fbpar), wfb=float(wfb))
    return data, params


def from_export(ex) -> tuple[InvertData, InvertParams]:
    """Build ``(InvertData, InvertParams)`` from a loaded export (``io.export.Export``).

    Reads the spec 09 names (``/dist_fn``, ``/le_grids``, ``/kt_grids``,
    ``/run_parameters``, ``/species``, ``/dist_fn_arrays``). The bounce-point
    masks and ``is_passing_hybrid_electron`` are recomputed (spec 03 section
    3.1, le_grids.f90:4490-4499) when the export does not carry them.
    """
    dfn = ex["/dist_fn"]
    le = ex["/le_grids"]
    kt = ex["/kt_grids"]
    rp = ex["/run_parameters"]
    sp = ex["/species"]
    dfa = ex["/dist_fn_arrays"]
    ntgrid = ex.ntgrid
    ng2 = int(_scalar(le, "ng2"))
    nlambda = ex.nlambda
    naky, ntheta0 = ex.naky, ex.ntheta0
    passing_wfb = _flag(le["passing_wfb"])
    trapped_wfb = _flag(le["trapped_wfb"])
    mixed_wfb = _flag(le["mixed_wfb"])
    forbid = np.asarray(le["forbid"], bool)
    is_ttp = np.asarray(le["is_ttp"], bool)
    lbp = _opt_array(le, "is_lower_bounce_point")
    ubp = _opt_array(le, "is_upper_bounce_point")
    if lbp is None or ubp is None:
        bp = np.asarray(le["is_bounce_point"], bool)
        lbp = lower_bounce_points(bp, forbid, trapped_wfb, ng2)
        ubp = upper_bounce_points(bp, forbid, trapped_wfb, ng2)
    aky = np.asarray(kt["aky"], np.float64)
    hyb = _opt_array(le, "is_passing_hybrid_electron")
    if hyb is None:
        stype = np.asarray(sp["spec_type"])
        il = np.arange(1, nlambda + 1)
        hyb = ((stype == HYBRID_ELECTRON_SPECIES)[:, None, None] & (aky != 0.0)[None, :, None]
               & (il <= ng2 + 1)[None, None, :])
    kwork = np.asarray(kt["kwork_filter"], bool)
    ieqzip = _opt_array(rp, "ieqzip", np.zeros((naky, ntheta0), bool)).astype(bool)
    bos = int(_scalar(dfn, "boundary_option_switch"))
    nonad_zero = _flag(_scalar(dfn, "nonad_zero"))
    box = _flag(kt["box"]) if "box" in kt else True
    wfb_factors = None
    names = ("wfb_aj0l", "wfb_aj1l", "wfb_aj0r", "wfb_aj1r", "wfb_vperp2left", "wfb_vperp2right")
    if all(name in dfn for name in names):
        wfb_factors = {name: np.asarray(dfn[name]) for name in names}
    no_conn = _opt_array(dfn, "no_connections")
    return build(
        ntgrid=ntgrid, ng2=ng2, aky=aky, kwork_filter=kwork, ieqzip=ieqzip, forbid=forbid,
        is_ttp=is_ttp, is_lower_bounce_point=lbp, is_upper_bounce_point=ubp,
        is_passing_hybrid_electron=hyb, passing_wfb=passing_wfb, trapped_wfb=trapped_wfb,
        mixed_wfb=mixed_wfb, boundary_option_switch=bos, nonad_zero=nonad_zero,
        zero_forbid=_flag(_scalar(dfn, "zero_forbid")),
        def_parity=_flag(_scalar(dfn, "def_parity")), even=_flag(_scalar(dfn, "even", True)),
        wfb=float(_scalar(dfn, "wfb")), fphi=float(_scalar(rp, "fphi")),
        fbpar=float(_scalar(rp, "fbpar")),
        vperp2=dfa["vperp2"], aj0=dfa["aj0"], aj1=dfa["aj1"], zt=sp["spec_zt"],
        nonmaxw_corr=sp["nonmaxw_corr"],
        l_links=_opt_array(dfn, "l_links"), r_links=_opt_array(dfn, "r_links"),
        itleft=_opt_array(dfn, "itleft"), itright=_opt_array(dfn, "itright"),
        no_connections=None if no_conn is None else bool(no_conn.item()),
        wfb_factors=wfb_factors, box=box,
        start_from_previous_solution=_flag(_scalar(dfn, "start_from_previous_solution", False)),
        boundary_off_grid=_flag(_scalar(dfn, "boundary_off_grid", False)),
        exponential_boundary=_flag(_scalar(dfn, "exponential_boundary", False)))


# ============================================================================ device helpers
def _tf(x):
    """``(..., T)`` -> theta-first ``(T, ...)``."""
    return jnp.moveaxis(x, -1, 0)


def _tl(x):
    """theta-first ``(T, ...)`` -> ``(..., T)``."""
    return jnp.moveaxis(x, 0, -1)


def _lt(m):
    """``(L, T)`` -> ``(T, 1, 1, L, 1, 1)`` against theta-first ``(T, S, E, L, Y, X)``."""
    return jnp.asarray(m).T[:, None, None, :, None, None]


def _l(m):
    """``(L,)`` -> ``(1, 1, L, 1, 1)`` against the batch ``(S, E, L, Y, X)``."""
    return jnp.asarray(m)[None, None, :, None, None]


def _ly(m):
    """``(L, Y)`` -> ``(1, 1, L, Y, 1)``."""
    return jnp.asarray(m)[None, None, :, :, None]


def _lyx(m):
    """``(L, Y, X)`` -> ``(1, 1, L, Y, X)``."""
    return jnp.asarray(m)[None, None]


def _end(arr6, t):
    """Value of a G6 array at theta index ``t``: ``(S, E, L, Y, X)``."""
    return jnp.asarray(arr6)[..., t]


def _field_end(field, t):
    """F3 field at theta index ``t``, broadcast against the batch: ``(1, 1, 1, Y, X)``."""
    return jnp.asarray(field)[..., t][None, None, None]


def get_adjust_end(t, phinew, bparnew, data: InvertData, params: InvertParams):
    """``get_adjust(ig, ..., fphi, fbpar, phinew, bparnew)`` at theta index ``t``
    for every line (src/dist_fn_arrays.fpp:178-188), ``(S, E, L, Y, X)``::

        2*vperp2*aj1*bpar*facbpar + zt*phi*aj0*facphi*nonmaxw_corr   (left to right)
    """
    coef = (2.0 * _end(data.vperp2, t)) * _end(data.aj1, t)
    term_b = rmul(rmul(_field_end(bparnew, t), coef), params.fbpar)
    zt = jnp.asarray(data.zt)[:, None, None, None, None]
    nm = jnp.asarray(data.nonmaxw_corr)[:, :, None, None, None]
    term_p = rmul(rmul(rmul(rmul(_field_end(phinew, t), zt), _end(data.aj0, t)),
                       params.fphi), nm)
    return term_b + term_p


def _self_periodic_adj(t, phinew, bparnew, data: InvertData, params: InvertParams):
    """``adjl``/``adjr`` of ``self_periodic`` (dist_fn.fpp:5879-5888)::

        2.0*vperp2*aj1*bparnew*fbpar + zt*phinew*aj0*nonmaxw_corr*fphi   (left to right)

    (note ``... *nonmaxw_corr*fphi`` here versus ``... *facphi*nonmaxw_corr``
    in ``get_adjust``; 03 section 6.3).
    """
    coef = (2.0 * _end(data.vperp2, t)) * _end(data.aj1, t)
    term_b = rmul(rmul(_field_end(bparnew, t), coef), params.fbpar)
    zt = jnp.asarray(data.zt)[:, None, None, None, None]
    nm = jnp.asarray(data.nonmaxw_corr)[:, :, None, None, None]
    term_p = rmul(rmul(rmul(rmul(_field_end(phinew, t), zt), _end(data.aj0, t)), nm),
                  params.fphi)
    return term_b + term_p


def set_h_zero(phinew, bparnew, data: InvertData, params: InvertParams):
    """``set_h_zero`` values (src/dist_fn_arrays.fpp:246-294), G6 ``(S, E, L, Y, X, T)``::

        facphi = -fphi ; facbpar = -fbpar
        phi_factor  = (zt*nonmaxw_corr)*facphi   (0 if facphi == 0)
        bpar_factor = 2.0*facbpar                (0 if facbpar == 0)
        adj = ((bpar_factor*vperp2)*aj1)*bpar + (phi_factor*phi)*aj0

    The caller applies it to the passing hybrid-electron lines (both signs).
    """
    facphi = -params.fphi
    facbpar = -params.fbpar
    zt = jnp.asarray(data.zt)[:, None]
    nm = jnp.asarray(data.nonmaxw_corr)
    if facphi != 0:
        phi_factor = (zt * nm) * facphi
    else:
        phi_factor = jnp.zeros_like(nm)
    bpar_factor = 2.0 * facbpar if facbpar != 0 else 0.0
    coef_b = (bpar_factor * jnp.asarray(data.vperp2)) * jnp.asarray(data.aj1)
    term_b = rmul(jnp.asarray(bparnew)[None, None, None], coef_b)
    term_p = rmul(rmul(jnp.asarray(phinew)[None, None, None],
                       phi_factor[:, :, None, None, None, None]), jnp.asarray(data.aj0))
    return term_b + term_p


def _sweep_down(c0, t, r):
    """Sign-2 recurrence (dist_fn.fpp:5596-5598 / 2071-2073), right to left::

        do ig = ntgrid-1, -ntgrid, -1:  g(ig) = -g(ig+1)*r(ig) + t(ig)

    ``c0``: ``g(ntgrid)`` (batch); ``t``, ``r``: ``(T-1, ...)`` for
    ``ig = -ntgrid..ntgrid-1``. Returns ``g`` theta-first ``(T, ...)``.
    """
    def step(c, x):
        ti, ri = x
        c = ti - cmul(c, ri)
        return c, c

    _, ys = lax.scan(step, c0, (t, r), reverse=True)
    return jnp.concatenate([ys, c0[None]], axis=0)


def _sweep_up(c0, t, r, ttp=None, copy=None, other=None):
    """Sign-1 recurrence, left to right (dist_fn.fpp:5634-5645 / 2108-2114)::

        do ig = -ntgrid, ntgrid-1
           if (.not. is_ttp(ig)) then: g(ig+1) = t(ig) - g(ig)*r(ig)
           else: if (copy(ig)) g(ig) = other(ig)      ! g(ig+1) keeps its initial 0

    ``c0``: ``g(-ntgrid)``; ``t, r, ttp, copy, other``: ``(T-1, ...)`` for
    ``ig = -ntgrid..ntgrid-1`` (``ttp``/``copy`` broadcastable masks). The
    carry is the stored value of ``g(ig)`` when iteration ``ig`` starts.
    """
    if ttp is None:
        def step(c, x):
            ti, ri = x
            nxt = ti - cmul(c, ri)
            return nxt, c

        last, ys = lax.scan(step, c0, (t, r))
    else:
        if copy is None:
            def step(c, x):
                ti, ri, tp = x
                nxt = ti - cmul(c, ri)
                return jnp.where(tp, jnp.zeros_like(nxt), nxt), c

            last, ys = lax.scan(step, c0, (t, r, ttp))
        else:
            def step(c, x):
                ti, ri, tp, cp, oi = x
                nxt = ti - cmul(c, ri)
                out = jnp.where(cp, oi, c)
                return jnp.where(tp, jnp.zeros_like(nxt), nxt), out

            last, ys = lax.scan(step, c0, (t, r, ttp, copy, other))
    return jnp.concatenate([ys, last[None]], axis=0)


def _split_signs(arr7):
    """G7 -> two theta-first arrays ``(T, S, E, L, Y, X)`` (sign 1, sign 2)."""
    arr7 = jnp.asarray(arr7)
    return _tf(arr7[..., 0, :]), _tf(arr7[..., 1, :])


def _join_signs(g1, g2):
    """Two theta-first arrays -> G7."""
    return jnp.stack([_tl(g1), _tl(g2)], axis=-2)


# ============================================================================ kernels
def init_homogeneous_g(coeffs, data: InvertData, params: InvertParams):
    """Mirrors ``init_homogeneous_g`` (src/dist_fn.fpp:1937-2125), ``force_passing`` absent.

    Returns ``g_h`` (G7 complex) from ``coeffs.r``. The ``links_h``/``wfb_h``
    fills into ``g_adj`` (2120-2123) are not needed: ``linked`` gathers the
    endpoint values of ``g_h`` directly. Recompute on every dt change (G-19, D3).
    """
    p, d = params, data
    r1, r2 = _split_signs(coeffs.r)
    batch = r1.shape[1:]
    cdt = r1.dtype
    bc1 = jnp.broadcast_to(_ly(d.gh_bc1), batch).astype(cdt)     # g_h(-ntgrid, 1) initial
    bc2 = jnp.broadcast_to(_ly(d.gh_bc2), batch).astype(cdt)     # g_h(ntgrid, 2) initial
    g2up = _lt(d.gh_g2up).astype(cdt)                             # g2(ig, 2), (T,1,1,L,1,1)
    # sign 2: g_h(ig,2) = -g_h(ig+1,2)*r(ig,2) + g2(ig,2)
    # GS2 quirk: g2(ntgrid,2) is set at an upper bounce point at ig = ntgrid but never
    # used by the sweep (dist_fn.fpp:2054-2073; 03 §6.9)
    gh2 = _sweep_down(bc2, g2up[:-1], r2[:-1])
    # g2(ig+1,1) = g_h(ig+1,2) at lower bounce points (2080-2096)
    if p.has_trapped:
        g2low = jnp.where(_lt(d.lbsrc)[:-1], gh2[1:], jnp.zeros_like(gh2[1:]))
    else:
        g2low = jnp.zeros_like(gh2[1:])
    c0 = bc1
    if p.trapped_wfb and p.has_wfb:
        c0 = jnp.where(_l(d.wfb_trapped), gh2[0], c0)              # 2102-2105
    # sign 1: g_h(ig+1,1) = -g_h(ig,1)*r(ig,1) + g2(ig+1,1), skipped after ttp points
    if p.has_ttp:
        gh1 = _sweep_up(c0, g2low, r1[:-1], ttp=_lt(d.ttp)[:-1])
    else:
        gh1 = _sweep_up(c0, g2low, r1[:-1])
    skip = _ly(d.gh_skip)                                          # il < ilmin: g_h == 0
    zero = jnp.zeros((), cdt)
    gh1 = jnp.where(skip, zero, gh1)
    gh2 = jnp.where(skip, zero, gh2)
    return _join_signs(gh1, gh2)


def invert_rhs_1(source, phinew, bparnew, coeffs, data: InvertData, params: InvertParams,
                 *, gnew_in=None, kwork_filter=None):
    """Mirrors ``invert_rhs_1`` (src/dist_fn.fpp:5440-5707) including its final
    ``enforce_parity`` (5706). Returns ``gnew`` (G7 complex).

    ``source``: G7 ``g_work`` from the source term (its in-place lower-bounce
    modification, 5619, stays local); ``phinew``, ``bparnew``: F3 fields of this
    ``timeadv`` call; ``coeffs``: ``types.Coeffs`` (``ainv``, ``r``, ``g_h``;
    ``g_h`` is computed here with :func:`init_homogeneous_g` if None).

    ``gnew_in`` (keyword): ``gnew`` as it is before the call. GS2 leaves the
    ``kwork_filter``/``ieqzip`` cells untouched (5475-5477; box runs always have
    ``kwork_filter(1,1)``), so it is required whenever such a cell exists.
    ``kwork_filter`` (keyword): optional ``(Y, X)`` bool override of the exported
    filter (e.g. all False to mirror the initial response-matrix build).

    Per line, in Fortran order: pass-through cells; passing hybrid electrons
    (``set_h_zero``); ``gnew = 0``; incoming boundary values; sign-2 sweep;
    trapped lower-bounce source; trapped-wfb copy at ``-ntgrid``; sign-1 sweep
    with the ttp branch; ``self_periodic``; trapped homogeneous correction per
    well; trapped-wfb correction; ``zero_forbid``; then ``enforce_parity``.
    """
    p, d = params, data
    passthrough = p.any_passthrough or kwork_filter is not None
    if passthrough and gnew_in is None:
        raise ValueError("invert_rhs_1: gnew_in is required (kwork_filter/ieqzip cells keep "
                         "their incoming gnew, dist_fn.fpp:5475-5477)")
    nth = 2 * p.ntgrid + 1
    s1, s2 = _split_signs(source)
    a1, a2 = _split_signs(coeffs.ainv)
    r1, r2 = _split_signs(coeffs.r)
    g_h = coeffs.g_h
    if g_h is None:
        g_h = init_homogeneous_g(coeffs, data, params)
    gh1, gh2 = _split_signs(g_h)
    cdt = s1.dtype
    zero = jnp.zeros(s1.shape[1:], cdt)

    # ---- incoming boundary values: apply_parallel_boundary_conditions (5580-5583, 5833-5847)
    if p.nonad_zero:
        g1_0 = jnp.where(_lyx(d.bc_left), -get_adjust_end(0, phinew, bparnew, d, p), zero)
        g2_0 = jnp.where(_lyx(d.bc_right), -get_adjust_end(nth - 1, phinew, bparnew, d, p), zero)
    else:
        g1_0 = zero                       # nonad_zero = F leaves the initial gnew = 0
        g2_0 = zero

    # ---- sign 2, right to left (5596-5598): -gnew(ig+1,2)*r(ig,2) + ainv(ig,2)*source(ig,2)
    g2 = _sweep_down(g2_0, cmul(a2[:-1], s2[:-1]), r2[:-1])

    # ---- trapped: source(ig,1) = gnew(ig+1,2) left of a lower bounce point (5604-5622)
    if p.has_trapped:
        s1 = jnp.concatenate([jnp.where(_lt(d.lbsrc)[:-1], g2[1:], s1[:-1]), s1[-1:]], axis=0)

    # ---- trapped wfb: gnew(-ntgrid,1) = gnew(-ntgrid,2) (5628-5631)
    if p.trapped_wfb and p.has_wfb:
        g1_0 = jnp.where(_l(d.wfb_trapped), g2[0], g1_0)

    # ---- sign 1, left to right (5634-5645)
    # GS2 quirk: at a ttp point the sweep writes gnew(ig,1) = gnew(ig,2), not gnew(ig+1,1),
    # which keeps its initial 0 (dist_fn.fpp:5643; 03 §6.10)
    t1 = cmul(a1[:-1], s1[:-1])
    if p.has_ttp:
        g1 = _sweep_up(g1_0, t1, r1[:-1], ttp=_lt(d.ttp)[:-1], copy=_lt(d.ttp_copy)[:-1],
                       other=g2[:-1])
    else:
        g1 = _sweep_up(g1_0, t1, r1[:-1])

    # ---- self_periodic (5651, 5864-5906)
    if p.any_self_periodic:
        sp = _lyx(d.self_periodic)
        if p.nonad_zero and p.has_wfb and not p.passing_wfb:      # static form of sp_dadj
            adjl = _self_periodic_adj(0, phinew, bparnew, d, p)
            adjr = _self_periodic_adj(nth - 1, phinew, bparnew, d, p)
            dadj = jnp.where(_l(d.sp_dadj), adjl - adjr, zero)
        else:
            dadj = zero                    # dadj = cmplx(0.0, 0.0)
        h1e = gh1[-1]                      # g_h(ntgrid, 1)
        c1 = sp & ((h1e.real != 1.0) | (h1e.imag != 0.0))
        beta = cdiv((g1[-1] - g1[0]) - dadj, rsub(1.0, h1e))
        g1 = jnp.where(c1, g1 + cmul(beta, gh1), g1)
        h2e = gh2[0]                       # g_h(-ntgrid, 2)
        c2 = sp & ((h2e.real != 1.0) | (h2e.imag != 0.0))
        beta = cdiv((g2[0] - g2[-1]) + dadj, rsub(1.0, h2e))
        g2 = jnp.where(c2, g2 + cmul(beta, gh2), g2)

    # ---- trapped: homogeneous solution well by well (5658-5671), as a static gather
    if p.has_trapped:
        beta_all = cdiv(g1 - g2, rsub(1.0, gh1))
        beta = jnp.take_along_axis(beta_all, _lt(d.beta_top), axis=0)
        beta = jnp.where(_lt(d.beta_has_top), beta, jnp.zeros((), cdt))   # beta1 = 0.0
        upd = _lt(d.beta_upd)
        g1 = jnp.where(upd, g1 + cmul(beta, gh1), g1)
        g2 = jnp.where(upd, g2 + cmul(beta, gh2), g2)

    # ---- trapped wfb: bounce at the right end (5679-5688)
    if p.trapped_wfb and p.has_wfb:
        m = _l(d.wfb_trapped)
        beta = cdiv(g1[-1] - g2[-1], rsub(1.0, gh1[-1]))
        g1n = g1 + cmul(beta, gh1)
        g2n = g2 + cmul(beta, gh2)
        g1 = jnp.where(m, g1n, g1)
        g2 = jnp.where(m, g2n, g2)

    # ---- zero_forbid (5696-5701)
    if p.zero_forbid:
        fb = _lt(d.forbid)
        g1 = jnp.where(fb, jnp.zeros((), cdt), g1)
        g2 = jnp.where(fb, jnp.zeros((), cdt), g2)

    out = _join_signs(g1, g2)

    # ---- passing hybrid electrons: h = 0 (5487-5490)
    if p.any_hybrid:
        hyb = jnp.asarray(d.hybrid).transpose(0, 2, 1)[:, None, :, :, None, None, None]
        out = jnp.where(hyb, set_h_zero(phinew, bparnew, d, p)[..., None, :], out)

    # ---- kwork_filter / ieqzip cells keep their gnew (5475-5477)
    if passthrough:
        kw = jnp.asarray(d.kwork_filter if kwork_filter is None else kwork_filter, bool)
        pt = (kw | jnp.asarray(d.ieqzip, bool))[None, None, None, :, :, None, None]
        out = jnp.where(pt, jnp.asarray(gnew_in), out)

    # ---- enforce_parity (5706)
    if p.def_parity:
        out = _parity.enforce_parity(out, data, params)
    return out


def invert_rhs(source, phinew, bparnew, coeffs, data: InvertData, params: InvertParams,
               *, gnew_in=None, kwork_filter=None):
    """``invert_rhs`` after the source (src/dist_fn.fpp:6184-6191): ``invert_rhs_1``
    followed by the linked boundary conditions when ``boundary_option = 'linked'``."""
    if coeffs.g_h is None:
        coeffs = coeffs._replace(g_h=init_homogeneous_g(coeffs, data, params))
    gnew = invert_rhs_1(source, phinew, bparnew, coeffs, data, params, gnew_in=gnew_in,
                        kwork_filter=kwork_filter)
    if params.linked:
        gnew = _linked.apply_linked_boundary_conditions(gnew, phinew, bparnew, coeffs, data,
                                                        params)
    return gnew
