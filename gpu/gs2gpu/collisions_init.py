"""Collision operator: static data, parameters and the dt-dependent coefficients.

Spec 06 (sections 2-5) with the resolutions of spec 10 (A.11, B.5, G-01, G-15, G-24..G-27).
Only the default ``use_le_layout = .true.`` path is supported (G-24).

Containers
----------
``CollParams``   static sizes and switches (frozen dataclass, hashable: a ``jit`` static argument).
``CollData``     dt-independent arrays: the runtime arrays used by ``collisions.vspace_derivatives``
                 and, in ``init``, the host arrays the dt rebuild needs.
``CollCoeffs``   the dt/vnmult-dependent matrices ``c1le, betaale, qle`` (Lorentz LU),
                 ``ec1le, ebetaale, eqle`` (energy-diffusion LU) and the conservation vectors
                 ``z0le, s0le, w0le, bz0le, bs0le, bw0le``, plus the ``code_dt`` and ``vnmult``
                 they were built with.

LE-canonical arrays (spec 06 section 2.4) have NumPy shape ``(S, Y, X, T, NE, NX)``: Fortran
``x(ixi, ie, ile)`` with ``ile`` expanded to ``(ig, it, ik, is)``. Coefficients keep the export
width ``nxi_lim``; the runtime ``gle`` has ``nxi + 1`` rows (equal to ``nxi_lim`` on a trapped grid).

``compute_coll_coeffs`` rebuilds every ``CollCoeffs`` array on the host in NumPy float64 with the
Fortran association order and the ``integrate_moment_r34`` summation order of spec 06 section 5.5.
NumPy never contracts ``a*b + c`` into an FMA, so with the ``-ffp-contract=off`` reference build
the rebuild is expected to be bitwise identical to GS2.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import NamedTuple, Optional

import numpy as np

from .gadjust import GAdjustData
from .io.export import Export

#: Fortran ``epsilon(0.0)``: GS2 is built with ``-fdefault-real-8`` (spec 06 section 4).
EPS = float(np.finfo(np.float64).eps)

# collision_model_switch values (src/collisions.fpp:60-64)
MODEL_LORENTZ = 1
MODEL_NONE = 3
MODEL_LORENTZ_TEST = 5
MODEL_FULL = 6
MODEL_EDIFFUSE = 7
# lorentz_switch / ediff_switch values (src/collisions.fpp:66-70)
SCHEME_DEFAULT = 1
SCHEME_OLD = 2
# species types (src/species.f90:37-43); electron-like = electron or hybrid
ELECTRON_LIKE = (2, 4)
HYBRID_ELECTRON = 4

LORENTZ_MATRICES = ("c1le", "betaale", "qle")
EDIFFUSE_MATRICES = ("ec1le", "ebetaale", "eqle")
LORENTZ_VECTORS = ("z0le", "s0le", "w0le")
EDIFFUSE_VECTORS = ("bz0le", "bs0le", "bw0le")
COEFF_NAMES = LORENTZ_MATRICES + EDIFFUSE_MATRICES + LORENTZ_VECTORS + EDIFFUSE_VECTORS

# Scalars of /collisions read by CollParams (spec 09 section 3.10).
_FLAG_DEFAULTS = {
    "colls": False, "collision_model_switch": MODEL_NONE, "has_lorentz": False,
    "has_diffuse": False, "drag": False, "use_le_layout": True, "conserve_moments": True,
    "conservative": True, "resistivity": True, "special_wfb_lorentz": False,
    "vpar_zero_mean": True, "conserve_forbid_zero": True, "force_collisions": False,
    "adjust": True, "heating": False, "split_collisions": False, "hypermult": False,
    "hyper_colls": False, "lorentz_switch": SCHEME_DEFAULT, "ediff_switch": SCHEME_DEFAULT,
    "cfac": 1.0, "nxi_lim": None, "vary_vnew": False,
}
# Arrays of /collisions (spec 09 section 3.10); absent ones become None.
_COLL_ARRAYS = ("vnew", "vnew_s", "vnew_D", "vnew_E", "delvnew", "vnewh", "pitch_weights",
                "vpdiff", "vpdiffle", "aj0le", "vperp_aj1le", "vpa_aj0_le", "e15", "capg_half",
                "vnmult")


def _scalar(value, kind):
    if value is None:
        return None
    arr = np.asarray(value)
    return kind(arr.reshape(()).item() if arr.size == 1 else arr)


def _layout_of(ex: Export) -> str:
    layout = ex.attrs.get("layout")
    if layout is None and "gs2_layouts" in ex.root.groups:
        layout = ex["/gs2_layouts"].attrs.get("layout")
    return str(layout).strip() if layout is not None else "lxyes"


def raw_from_export(ex: Export) -> dict:
    """Every export array/scalar the collision operator reads, keyed by Fortran name.

    Shapes are the loader's canonical NumPy shapes (spec 06 section 2.4, spec 10 G-29). The
    ``/collisions`` group exists only when GS2 initialised collisions; without it ``colls`` is False.
    """
    le, th, kt = ex["/le_grids"], ex["/theta_grid"], ex["/kt_grids"]
    rp, sp, dfa, tm = ex["/run_parameters"], ex["/species"], ex["/dist_fn_arrays"], ex["/gs2_time"]
    raw = dict(
        ntgrid=ex.ntgrid, nth=ex.nth, ntheta0=ex.ntheta0, naky=ex.naky, nlambda=ex.nlambda,
        negrid=ex.negrid, nspec=ex.nspec, ng2=int(le["ng2"]), nxi=int(le["nxi"]),
        trapped=_scalar(le["grid_has_trapped_particles"], bool), layout=_layout_of(ex),
        fphi=float(rp["fphi"]), fbpar=float(rp["fbpar"]), beta=float(rp["beta"]),
        bmag=th["bmag"], al=le["al"], wl=le["wl"], jend=le["jend"], forbid=le["forbid"],
        xi=le["xi"], wxi=le["wxi"], ixi_to_il=le["ixi_to_il"], ixi_to_isgn=le["ixi_to_isgn"],
        w=le["w"], energy_maxw=le["energy_maxw"], speed_maxw=le["speed_maxw"],
        w_maxw=le["w_maxw"], kperp2=kt["kperp2"], kwork_filter=kt["kwork_filter"],
        ieqzip=rp["ieqzip"], tunits=tm["tunits"],
        spec_type=sp["spec_type"], spec_vnewk=sp["spec_vnewk"], spec_z=sp["spec_z"],
        spec_dens=sp["spec_dens"], spec_stm=sp["spec_stm"], spec_zstm=sp["spec_zstm"],
        spec_smz=sp["spec_smz"], spec_zt=sp["spec_zt"], nonmaxw_corr=sp["nonmaxw_corr"],
        vpa=dfa["vpa"], aj0=dfa["aj0"], aj1=dfa["aj1"], vperp2=dfa["vperp2"],
    )
    coll = ex["/collisions"] if "collisions" in ex.root.groups else None
    for name, default in _FLAG_DEFAULTS.items():
        raw[name] = coll.get(name, default) if coll is not None else default
    for name in _COLL_ARRAYS:
        raw[name] = coll.get(name) if coll is not None else None
    if coll is None:
        raw["colls"] = False
    return raw


# --------------------------------------------------------------------------- parameters
@dataclass(frozen=True)
class CollParams:
    """Static sizes and switches of the collision operator (hashable).

    ``nx = nxi + 1`` is the number of xi rows of the runtime ``gle`` (le_derivatives.fpp:79);
    ``nxi_lim`` is the width of the coefficient arrays (collisions.fpp:557-561). ``ie_outer`` is
    the velocity order of ``integrate_moment_r34`` (spec 06 section 5.5). ``exact_sums`` selects
    the sequential Fortran order (ie outer, ixi inner) for the runtime conservation moments
    (``False``: a tree reduction, round-off level differences only, spec 06 section 7).
    """

    nspec: int
    naky: int
    ntheta0: int
    nth: int
    ntgrid: int
    nlambda: int
    negrid: int
    ng2: int
    nxi: int
    nxi_lim: int
    nx: int
    trapped: bool
    colls: bool
    collision_model_switch: int
    has_lorentz: bool
    has_diffuse: bool
    drag: bool
    conserve_moments: bool
    conservative: bool
    special_wfb_lorentz: bool
    vpar_zero_mean: bool
    conserve_forbid_zero: bool
    force_collisions: bool
    adjust: bool
    heating: bool
    split_collisions: bool
    hypermult: bool
    hyper_colls: bool
    cfac: float
    fphi: float
    fbpar: float
    layout: str
    ie_outer: bool
    any_wfb_average: bool
    any_special_shift: bool
    exact_sums: bool = True

    @property
    def lorentz_test(self) -> bool:
        return self.collision_model_switch == MODEL_LORENTZ_TEST

    @classmethod
    def from_export(cls, ex: Export, **overrides) -> "CollParams":
        """Parameters from the export file (gates checked, spec 10 B.5)."""
        return cls.from_raw(raw_from_export(ex), **overrides)

    @classmethod
    def from_raw(cls, raw: dict, **overrides) -> "CollParams":
        """Parameters from a ``raw_from_export``-style dict (also used by synthetic tests).

        Mirrors read_parameters (src/collisions.fpp:466-562) for the derived flags.
        Raises ``NotImplementedError`` for gated options (spec 10 B.5, G-24, G-25; spec 06 section 8.1).
        """
        colls = bool(raw["colls"])
        model = int(raw["collision_model_switch"])
        if model == MODEL_NONE:
            colls = False
        has_lorentz = bool(raw["has_lorentz"]) and colls
        has_diffuse = bool(raw["has_diffuse"]) and colls
        ng2, nxi, nlambda = int(raw["ng2"]), int(raw["nxi"]), int(raw["nlambda"])
        trapped = bool(raw["trapped"])
        nxi_lim = raw.get("nxi_lim")
        nxi_lim = int(nxi_lim) if nxi_lim is not None else (nxi + 1 if nxi > 2 * ng2 else nxi)
        layout = str(raw.get("layout", "lxyes")).strip()
        if sorted(layout) != sorted("lxyes"):
            raise ValueError(f"unknown g_lo layout {layout!r}")
        conserve = bool(raw["conserve_moments"])
        if colls:
            if not bool(raw["use_le_layout"]):
                raise NotImplementedError(
                    "collisions_knobs:use_le_layout = .false. is GATED (spec 10 G-24)")
            if not bool(raw["conservative"]):
                raise NotImplementedError(
                    "collisions_knobs:conservative = .false. is GATED (vnew_E/delvnew end points "
                    "uninitialised, spec 06 section 0)")
            if has_lorentz and int(raw["lorentz_switch"]) == SCHEME_OLD:
                raise NotImplementedError("collisions_knobs:lorentz_scheme = 'old' is GATED")
            if has_diffuse and int(raw["ediff_switch"]) == SCHEME_OLD:
                raise NotImplementedError("collisions_knobs:ediff_scheme = 'old' is GATED")
            if np.any(np.isin(np.asarray(raw["spec_type"]), [HYBRID_ELECTRON])):
                raise NotImplementedError(
                    "hybrid electron species together with collisions are GATED (spec 10 G-25)")
            if (not trapped) and has_lorentz and (conserve or not bool(raw["use_le_layout"])):
                raise NotImplementedError(
                    "Lorentz collisions with conserve_moments on a grid without trapped particles "
                    "are GATED (out-of-bounds bug, spec 06 section 8.1)")
        jend = np.asarray(raw["jend"])
        is_wfb = jend == ng2 + 1
        special = bool(raw["special_wfb_lorentz"])
        vzm = bool(raw["vpar_zero_mean"])
        kw = dict(
            nspec=int(raw["nspec"]), naky=int(raw["naky"]), ntheta0=int(raw["ntheta0"]),
            nth=int(raw["nth"]), ntgrid=int(raw["ntgrid"]), nlambda=nlambda,
            negrid=int(raw["negrid"]), ng2=ng2, nxi=nxi, nxi_lim=nxi_lim, nx=nxi + 1,
            trapped=trapped, colls=colls, collision_model_switch=model,
            has_lorentz=has_lorentz, has_diffuse=has_diffuse,
            drag=bool(raw["drag"]) and has_lorentz, conserve_moments=conserve,
            conservative=bool(raw["conservative"]), special_wfb_lorentz=special,
            vpar_zero_mean=vzm, conserve_forbid_zero=bool(raw["conserve_forbid_zero"]),
            force_collisions=bool(raw["force_collisions"]), adjust=bool(raw["adjust"]),
            heating=bool(raw["heating"]), split_collisions=bool(raw["split_collisions"]),
            hypermult=bool(raw["hypermult"]), hyper_colls=bool(raw["hyper_colls"]),
            cfac=float(raw["cfac"]), fphi=float(raw["fphi"]), fbpar=float(raw["fbpar"]),
            layout=layout, ie_outer=layout.index("l") < layout.index("e"),
            any_wfb_average=bool(trapped and vzm and np.any(is_wfb)),
            any_special_shift=bool(special and np.any(is_wfb)),
        )
        kw.update(overrides)
        return cls(**kw)


# --------------------------------------------------------------------------- data
class CollInitData(NamedTuple):
    """Host (NumPy) inputs of the dt rebuild (spec 06 sections 5.1-5.5).

    Shapes: ``bmag (T,)``, ``al (L,)``, ``pitch_weights (T, L)``, ``xi (T, 2L)``,
    ``jend (T,)`` Fortran values, ``forbid, wl (L, T)``, ``w (S, E)``, ``energy, speed, w_maxw,
    e15 (E,)``, ``capg_half (E-1,)``, ``vnew, vnew_s, vnew_D, vnew_E, delvnew (S, E, Y)``,
    ``vnewh (S, Y, X, T)`` or None, ``tunits (Y,)``, ``kperp2 (Y, X, T)``, ``vnewk, zstm, smz (S,)``,
    ``is_electron (S,)``, ``aj0, aj1`` G6, ``vpa`` G7, ``vpdiff (L, 2, T)`` or None,
    ``kwork_filter`` (the exported, post-``init_fields`` value), ``ieqzip (Y, X)``, ``skip (S, Y)``.
    """

    bmag: np.ndarray
    al: np.ndarray
    pitch_weights: Optional[np.ndarray]
    xi: np.ndarray
    ixi_to_il: np.ndarray
    ixi_to_isgn: np.ndarray
    jend: np.ndarray
    forbid: np.ndarray
    wl: np.ndarray
    w: np.ndarray
    energy: np.ndarray
    speed: np.ndarray
    w_maxw: np.ndarray
    e15: np.ndarray
    capg_half: Optional[np.ndarray]
    vnew: np.ndarray
    vnew_s: Optional[np.ndarray]
    vnew_D: Optional[np.ndarray]
    vnew_E: Optional[np.ndarray]
    delvnew: Optional[np.ndarray]
    vnewh: Optional[np.ndarray]
    tunits: np.ndarray
    kperp2: np.ndarray
    vnewk: np.ndarray
    zstm: np.ndarray
    smz: np.ndarray
    is_electron: np.ndarray
    aj0: np.ndarray
    aj1: np.ndarray
    vpa: np.ndarray
    vpdiff: Optional[np.ndarray]
    kwork_filter: np.ndarray
    ieqzip: np.ndarray
    skip: np.ndarray


class CollData(NamedTuple):
    """dt-independent arrays of the collision operator (spec 06 section 3).

    Runtime arrays (device):

    * gather/scatter maps of ``g2le`` (le_grids.f90:3379-3599): ``le_cidx (T, NX)`` = ``2*il + isgn``
      (0-based) of each xi row, ``le_valid (T, NX)`` (row is an allowed point), ``g_ixi (L, 2, T)``
      (row of each allowed ``(il, isgn)``), ``g_allowed (L, T)`` (= ``.not. forbid``);
    * ``wxi (T, NX)``, ``w (S, E)``, ``speed (E,)``, ``vnew_D, vnew_E, delvnew (S, E, Y)``;
    * ``aj0le, vperp_aj1le, vpa_aj0_le`` LE6 ``(S, Y, X, T, E+1, nxi_lim)``, ``vpdiffle (T, nxi_lim)``;
    * drag: ``vnewk (S,)``, ``drag_den (S, E)`` = ``(((-(z*dens))*beta)*stm)*e15``, ``kperp2 (Y, X, T)``,
      ``is_electron (S,)``;
    * masks: ``apply (S, Y, X)`` (= not skip and active: the two solves), ``active (Y, X)``
      (= not kwork_filter and not ieqzip: every stage, spec 10 G-01), ``ed_rows (T, NX)`` (rows of the
      energy solve, ``ixi <= max(2*jend, 2*ng2)``);
    * Lorentz per-theta structure (solfp_lorentz_le_layout, collisions.fpp:3897-3980), 0-based:
      ``lz_je (T,)`` = ``JE-1``, ``lz_cj (T,)`` = ``jend-1``, ``lz_cj2 (T,)`` = ``2*jend-1``,
      ``lz_bounce (T,)`` (``jend /= 0``), ``lz_avg (T,)`` (wfb vpar=0 averaging), ``lz_special (T,)``,
      ``lz_shift, lz_shift_zero, lz_unshift (T, NX)`` (special_wfb_lorentz shifts);
    * ``gadj``: ``GAdjustData`` for the g <-> h transform.

    ``init``: ``CollInitData`` (host) for ``compute_coll_coeffs``; ``device()`` drops it.
    """

    le_cidx: np.ndarray
    le_valid: np.ndarray
    g_ixi: np.ndarray
    g_allowed: np.ndarray
    wxi: np.ndarray
    w: np.ndarray
    speed: np.ndarray
    vnew_D: Optional[np.ndarray]
    vnew_E: Optional[np.ndarray]
    delvnew: Optional[np.ndarray]
    aj0le: Optional[np.ndarray]
    vperp_aj1le: Optional[np.ndarray]
    vpa_aj0_le: Optional[np.ndarray]
    vpdiffle: Optional[np.ndarray]
    vnewk: np.ndarray
    drag_den: np.ndarray
    kperp2: np.ndarray
    is_electron: np.ndarray
    apply: np.ndarray
    active: np.ndarray
    ed_rows: np.ndarray
    lz_je: np.ndarray
    lz_cj: np.ndarray
    lz_cj2: np.ndarray
    lz_bounce: np.ndarray
    lz_avg: np.ndarray
    lz_special: np.ndarray
    lz_shift: np.ndarray
    lz_shift_zero: np.ndarray
    lz_unshift: np.ndarray
    gadj: GAdjustData
    init: Optional[CollInitData]

    @classmethod
    def from_export(cls, ex: Export, params: Optional[CollParams] = None) -> "CollData":
        raw = raw_from_export(ex)
        return cls.from_raw(raw, params if params is not None else CollParams.from_raw(raw))

    @classmethod
    def from_raw(cls, raw: dict, params: CollParams) -> "CollData":
        """Build the data from a ``raw_from_export``-style dict (see ``raw_from_export``)."""
        p = params
        S, E, L, Y, X, T = p.nspec, p.negrid, p.nlambda, p.naky, p.ntheta0, p.nth
        NX, XI, L2 = p.nx, p.nxi_lim, 2 * p.nlambda
        f64 = lambda a: None if a is None else np.asarray(a, np.float64)  # noqa: E731

        jend = np.asarray(raw["jend"]).astype(np.int64)
        forbid = np.asarray(raw["forbid"], bool)
        maps = le_maps(raw["ixi_to_il"], raw["ixi_to_isgn"], forbid, NX)
        wxi = np.zeros((T, NX))
        cols = min(NX, L2)
        wxi[:, :cols] = np.asarray(raw["wxi"], np.float64)[:, :cols]

        # vnew & co (S, E, Y); skip test of the two solves (collisions.fpp:3907, 4061)
        vnew = f64(raw["vnew"]) if raw.get("vnew") is not None else np.zeros((S, E, Y))
        skip = (np.abs(vnew[:, 0, :]) < 2.0 * EPS) & (not p.force_collisions)
        kwork = np.asarray(raw["kwork_filter"], bool)
        ieqzip = np.asarray(raw["ieqzip"], bool)
        active = ~kwork & ~ieqzip
        apply = (~skip)[:, :, None] & active[None, :, :]

        spec_type = np.asarray(raw["spec_type"])
        is_electron = np.isin(spec_type, ELECTRON_LIKE)
        energy = f64(raw["energy_maxw"])
        e15 = f64(raw.get("e15")) if raw.get("e15") is not None else energy ** 1.5
        z, dens, stm = f64(raw["spec_z"]), f64(raw["spec_dens"]), f64(raw["spec_stm"])
        # (((-(z*dens))*beta)*stm)*energy**1.5   (collisions.fpp:3074-3075)
        drag_den = (((-(z * dens)) * float(raw["beta"])) * stm)[:, None] * e15[None, :]

        # Lorentz per-theta structure (collisions.fpp:3917-3978)
        cj = jend
        JE = np.maximum(2 * cj, 2 * p.ng2 + 1)
        is_wfb = cj == p.ng2 + 1
        k = np.arange(NX)[None, :]
        lz_special = is_wfb & p.special_wfb_lorentz
        lz = dict(
            lz_je=(JE - 1).astype(np.int32),
            lz_cj=np.maximum(cj - 1, 0).astype(np.int32),
            lz_cj2=np.maximum(2 * cj - 1, 0).astype(np.int32),
            lz_bounce=cj != 0,
            lz_avg=(JE == 2 * cj) & is_wfb & p.vpar_zero_mean,
            lz_special=lz_special,
            lz_shift=lz_special[:, None] & (k >= p.ng2) & (k <= JE[:, None] - 3),
            lz_shift_zero=lz_special[:, None] & (k == JE[:, None] - 2),
            lz_unshift=lz_special[:, None] & (k >= p.ng2 + 1) & (k <= JE[:, None] - 2),
        )
        ed_rows = k < np.maximum(2 * cj, 2 * p.ng2)[:, None]

        # le-layout Bessel arrays: exported, or the init_le_bessel gather (collisions.fpp:634-687)
        aj0le, vperp_aj1le, vpa_aj0_le = (f64(raw.get(n)) for n in
                                          ("aj0le", "vperp_aj1le", "vpa_aj0_le"))
        if p.colls and (aj0le is None or vperp_aj1le is None or vpa_aj0_le is None):
            aj0le, vperp_aj1le, vpa_aj0_le = init_le_bessel(
                raw["aj0"], raw["aj1"], raw["vpa"], energy, f64(raw["al"]), maps, p)
        pitch_weights = f64(raw.get("pitch_weights"))
        vpdiff = f64(raw.get("vpdiff"))
        if vpdiff is None and p.has_lorentz and p.conservative and pitch_weights is not None:
            vpdiff = init_vpdiff(f64(raw["bmag"]), f64(raw["al"]), jend, pitch_weights, p)
        vpdiffle = f64(raw.get("vpdiffle"))
        if vpdiffle is None and vpdiff is not None:
            vpdiffle = vpdiff_le(vpdiff, raw["ixi_to_il"], raw["ixi_to_isgn"], XI)

        trapped = bool(p.trapped)
        il = np.arange(1, L + 1)[:, None]
        gadj = GAdjustData(
            vperp2=f64(raw["vperp2"]), aj0=f64(raw["aj0"]), aj1=f64(raw["aj1"]),
            zt=f64(raw["spec_zt"]), nonmaxw_corr=f64(raw["nonmaxw_corr"]),
            adjust_mask=~(trapped & (il > jend[None, :])))

        init = CollInitData(
            bmag=f64(raw["bmag"]), al=f64(raw["al"]), pitch_weights=pitch_weights,
            xi=f64(raw["xi"]), ixi_to_il=np.asarray(raw["ixi_to_il"]),
            ixi_to_isgn=np.asarray(raw["ixi_to_isgn"]), jend=jend, forbid=forbid,
            wl=f64(raw["wl"]), w=f64(raw["w"]),
            energy=energy, speed=f64(raw["speed_maxw"]), w_maxw=f64(raw["w_maxw"]), e15=e15,
            capg_half=f64(raw.get("capg_half")), vnew=vnew, vnew_s=f64(raw.get("vnew_s")),
            vnew_D=f64(raw.get("vnew_D")), vnew_E=f64(raw.get("vnew_E")),
            delvnew=f64(raw.get("delvnew")), vnewh=f64(raw.get("vnewh")),
            tunits=f64(raw["tunits"]), kperp2=f64(raw["kperp2"]), vnewk=f64(raw["spec_vnewk"]),
            zstm=f64(raw["spec_zstm"]), smz=f64(raw["spec_smz"]), is_electron=is_electron,
            aj0=f64(raw["aj0"]), aj1=f64(raw["aj1"]), vpa=f64(raw["vpa"]), vpdiff=vpdiff,
            kwork_filter=kwork, ieqzip=ieqzip, skip=skip)

        return cls(
            le_cidx=(2 * maps.le_il + maps.le_sg).astype(np.int32), le_valid=maps.le_valid,
            g_ixi=maps.g_ixi.astype(np.int32), g_allowed=maps.allowed, wxi=wxi, w=f64(raw["w"]),
            speed=f64(raw["speed_maxw"]), vnew_D=init.vnew_D, vnew_E=init.vnew_E,
            delvnew=init.delvnew, aj0le=aj0le, vperp_aj1le=vperp_aj1le, vpa_aj0_le=vpa_aj0_le,
            vpdiffle=vpdiffle, vnewk=init.vnewk, drag_den=drag_den, kperp2=init.kperp2,
            is_electron=is_electron, apply=apply, active=active, ed_rows=ed_rows,
            gadj=gadj, init=init, **lz)

    def device(self) -> "CollData":
        """Copy of the runtime arrays on the default JAX device, without ``init``."""
        from .types import to_jax

        return to_jax(self._replace(init=None))


# --------------------------------------------------------------------------- coefficients
class CollCoeffs(NamedTuple):
    """dt/vnmult-dependent collision arrays, LE6 ``(S, Y, X, T, NE, nxi_lim)`` (spec 06 section 5).

    ``NE = negrid`` for the LU arrays, ``negrid + 1`` for the conservation vectors. Arrays that GS2
    does not allocate for the chosen model are None. ``code_dt`` and ``vnmult`` (shape (2,)) are the
    values they were built with; the runtime conservation and drag terms read them from here.
    """

    code_dt: float
    vnmult: np.ndarray
    c1le: Optional[np.ndarray] = None
    betaale: Optional[np.ndarray] = None
    qle: Optional[np.ndarray] = None
    ec1le: Optional[np.ndarray] = None
    ebetaale: Optional[np.ndarray] = None
    eqle: Optional[np.ndarray] = None
    z0le: Optional[np.ndarray] = None
    s0le: Optional[np.ndarray] = None
    w0le: Optional[np.ndarray] = None
    bz0le: Optional[np.ndarray] = None
    bs0le: Optional[np.ndarray] = None
    bw0le: Optional[np.ndarray] = None

    @classmethod
    def from_export(cls, ex: Export, rung: Optional[int] = None) -> "CollCoeffs":
        """Coefficients at the initial dt (``/collisions``), or at dt-ladder rung ``rung``
        (``/dt_ladder/k<rung>/collisions``, 1-based, spec 09 section 6.3)."""
        if "collisions" not in ex.root.groups:
            return cls(code_dt=float(ex["/state_t0/code_dt"]), vnmult=np.ones(2))
        coll = ex["/collisions"]
        if rung is None:
            src = coll
            code_dt = float(ex["/state_t0/code_dt"])
        else:
            base = ex[f"/dt_ladder/k{rung}"]
            src = base["collisions"]
            code_dt = float(base.attrs["code_dt"])
        vnmult = np.asarray(coll["vnmult"], np.float64).reshape(2)
        arrays = {name: (None if src.get(name) is None else np.asarray(src.get(name), np.float64))
                  for name in COEFF_NAMES}
        return cls(code_dt=code_dt, vnmult=vnmult, **arrays)


# --------------------------------------------------------------------------- layout maps
class LeMaps(NamedTuple):
    """0-based index maps of the ``g2le`` redistribute (host NumPy)."""

    le_il: np.ndarray      # (T, NX) il of each xi row (0 on invalid rows)
    le_sg: np.ndarray      # (T, NX) isgn of each xi row
    le_valid: np.ndarray   # (T, NX) row holds an allowed (not forbidden) point
    g_ixi: np.ndarray      # (L, 2, T) xi row of each allowed (il, isgn, ig) (0 elsewhere)
    allowed: np.ndarray    # (L, T) = .not. forbid


def le_maps(ixi_to_il, ixi_to_isgn, forbid, nx: int) -> LeMaps:
    """Index maps of the ``g2le`` redistribute (setup_g2le_redistribute_local, le_grids.f90:3379-3599).

    Row ``ixi`` of ``gle(:, ie, ile)`` holds ``g(ig, ixi_to_isgn, ixi_to_il)``; rows beyond ``2*nlambda``
    (no-trapped grids, where ``nxi + 1 = 2*ng2 + 1``) are padding. Forbidden points are skipped in both
    directions (``skip_forbidden_region``, le_grids.f90:3403). Checks that every allowed point has
    exactly one row.
    """
    forbid = np.asarray(forbid, bool)
    L, T = forbid.shape
    il_map = np.asarray(ixi_to_il).astype(np.int64) - 1
    sg_map = np.asarray(ixi_to_isgn).astype(np.int64) - 1
    cols = min(nx, il_map.shape[1])
    le_il = np.zeros((T, nx), np.int64)
    le_sg = np.zeros((T, nx), np.int64)
    le_il[:, :cols] = il_map[:, :cols]
    le_sg[:, :cols] = sg_map[:, :cols]
    le_valid = np.zeros((T, nx), bool)
    tt = np.arange(T)[:, None]
    le_valid[:, :cols] = ~forbid[le_il[:, :cols], tt]
    g_ixi = np.full((L, 2, T), -1, np.int64)
    for t in range(T):
        for ixi in np.nonzero(le_valid[t])[0]:
            il, sg = le_il[t, ixi], le_sg[t, ixi]
            if g_ixi[il, sg, t] >= 0:
                raise ValueError(f"g2le map: (il={il + 1}, isgn={sg + 1}, t={t}) mapped twice")
            g_ixi[il, sg, t] = ixi
    allowed = ~forbid
    if np.any(g_ixi[:, 0][allowed] < 0) or np.any(g_ixi[:, 1][allowed] < 0):
        raise ValueError("g2le map: an allowed (il, isgn, theta) point has no xi row")
    g_ixi[g_ixi < 0] = 0
    le_il[~le_valid] = 0
    le_sg[~le_valid] = 0
    return LeMaps(le_il, le_sg, le_valid, g_ixi, allowed)


def gather_le(v, maps: LeMaps, nx: int):
    """``gather(g2le, v, out)`` into a zeroed buffer: G7 ``(S,E,L,Y,X,2,T)`` -> LE6 ``(S,Y,X,T,E+1,nx)``.

    Host version (NumPy) of a pure element copy (redistribute.f90:6425-6480).
    """
    v = np.asarray(v)
    S, E, L, Y, X, _, T = v.shape
    vt = v.transpose(0, 3, 4, 6, 1, 2, 5)                  # (S, Y, X, T, E, L, 2)
    tt = np.broadcast_to(np.arange(T)[:, None], maps.le_il.shape)
    core = vt[:, :, :, tt, :, maps.le_il, maps.le_sg]       # (T, NX, S, Y, X, E)
    core = np.where(maps.le_valid[:, :, None, None, None, None], core, 0)
    out = np.zeros((S, Y, X, T, E + 1, nx), v.dtype)
    out[:, :, :, :, :E, :] = core.transpose(2, 3, 4, 0, 5, 1)
    return out


def scatter_le(le, v, maps: LeMaps):
    """``scatter(g2le, le, v)``: allowed points of ``v`` (G7) take the LE6 values, forbidden points keep
    their values (host version)."""
    v = np.asarray(v)
    T = v.shape[-1]
    tt = np.broadcast_to(np.arange(T)[None, None, :], maps.g_ixi.shape)
    E = v.shape[1]
    got = np.asarray(le)[:, :, :, tt, :E, maps.g_ixi]      # (L, 2, T, S, Y, X, E)
    got = got.transpose(3, 6, 0, 4, 5, 1, 2)                 # (S, E, L, Y, X, 2, T)
    return np.where(maps.allowed[None, None, :, None, None, None, :], got, v)


def init_le_bessel(aj0, aj1, vpa, energy, al, maps, p: CollParams):
    """``init_le_bessel`` (src/collisions.fpp:634-687): ``aj0le``, ``vperp_aj1le = (aj1*E)*al``,
    ``vpa_aj0_le = vpa*aj0le``, LE6 ``(S, Y, X, T, E+1, nxi_lim)``; forbidden/padding entries 0."""
    aj0 = np.asarray(aj0, np.float64)
    aj1 = np.asarray(aj1, np.float64)
    XI = p.nxi_lim
    g0 = np.broadcast_to(aj0[..., None, :], aj0.shape[:-1] + (2, aj0.shape[-1]))
    g1 = np.broadcast_to(aj1[..., None, :], aj1.shape[:-1] + (2, aj1.shape[-1]))
    aj0le = gather_le(g0, maps, p.nx)[..., :XI]
    aj1le = gather_le(g1, maps, p.nx)[..., :XI]
    le_il = maps.le_il[:, :XI]                               # (T, XI) 0-based il
    E = energy.shape[0]
    vperp = np.zeros_like(aj1le)
    # vperp_aj1le(ixi,ie,ile) * energy(ie) * al(ixi_to_il(ixi,ig))   (collisions.fpp:675)
    vperp[..., :E, :] = (aj1le[..., :E, :] * energy[None, None, None, None, :, None]) \
        * al[le_il][None, None, None, :, None, :]
    vpale = gather_le(np.asarray(vpa, np.float64), maps, p.nx)[..., :XI]
    return aj0le, vperp, vpale * aj0le


def init_vpdiff(bmag, al, jend, pitch_weights, p: CollParams):
    """``init_vpdiff`` (src/collisions.fpp:4079-4140): ``vpdiff`` NumPy ``(L, 2, T)``."""
    L, T = p.nlambda, p.nth
    vp = np.zeros((L, 2, T))
    for t in range(T):
        je = int(jend[t])
        te = p.ng2 if (je <= p.ng2 or (je == p.ng2 + 1 and p.special_wfb_lorentz)) else je
        s = np.sqrt(np.maximum(0.0, 1.0 - bmag[t] * al))
        pw = pitch_weights[t]
        for il in range(2, te):
            slb0, slb1, slb2 = s[il - 2], s[il - 1], s[il]
            slbl = (slb1 + slb0) * 0.5
            slbr = (slb1 + slb2) * 0.5
            vp[il - 1, 0, t] = (slbl * slbl - slbr * slbr) / pw[il - 1]
        slb1, slb2 = s[0], s[1]
        slbr = 0.5 * (slb1 + slb2)
        vp[0, 0, t] = (1.0 - slbr * slbr) / pw[0]
        il = te
        slb0 = s[il - 2]
        if te == p.ng2:
            slb1 = s[il - 1]
            slb2 = -slb1
        else:
            slb1 = 0.0
            slb2 = -slb0
        slbl = (slb1 + slb0) * 0.5
        slbr = (slb1 + slb2) * 0.5
        if abs(pw[il - 1]) <= EPS:
            vp[il - 1, 0, t] = 0.0
        else:
            vp[il - 1, 0, t] = (slbl * slbl - slbr * slbr) / pw[il - 1]
        vp[:, 1, t] = -vp[:, 0, t]
    return vp


def vpdiff_le(vpdiff, ixi_to_il, ixi_to_isgn, nxi_lim: int):
    """``vpdiffle(ixi, ig) = vpdiff(ig, ixi_to_isgn, ixi_to_il)`` (collisions.fpp:4142-4152): ``(T, nxi_lim)``."""
    il = np.asarray(ixi_to_il)[:, :nxi_lim] - 1
    sg = np.asarray(ixi_to_isgn)[:, :nxi_lim] - 1
    T = il.shape[0]
    tt = np.broadcast_to(np.arange(T)[:, None], il.shape)
    return np.asarray(vpdiff)[il, sg, tt]


# --------------------------------------------------------------------------- the dt rebuild
def compute_coll_coeffs(code_dt: float, vnmult, data: CollData, params: CollParams, *,
                        kwork_filter=None, max_g_elements: int = 1 << 26) -> CollCoeffs:
    """Rebuild every dt/vnmult-dependent collision array (spec 06 section 5, spec 10 E.3).

    Mirrors init_lorentz (src/collisions.fpp:2185-2358) + get_lorentz_matrix (2392-2639),
    init_ediffuse (1945-2038) + get_ediffuse_matrix (2041-2182), init_lorentz_conserve (691-1194)
    and init_diffuse_conserve (1198-1679), le layout, host NumPy float64.

    ``vnmult``: the two multipliers (post-init values, >= 1; spec 06 section 8.7).

    ``kwork_filter``: the filter in force while GS2 builds the arrays. The conservation vectors run
    the runtime solves, which skip filtered modes. At the first initialisation GS2 builds them before
    ``init_fields`` sets ``kwork_filter(1,1)`` in box runs (fields.f90:376), i.e. with an all-False
    filter (pass ``kwork_filter=False``); on every dt reset the box filter persists (spec 10 G-01), which
    is the default (``None``: the exported value). Only the filtered mode's vectors differ, and the runtime
    never reads them (it skips that mode).

    The build is done in ky chunks of at most ``max_g_elements`` G7 elements (identical numbers: every
    integral, solve and matrix is per ky).
    """
    if data.init is None:
        raise ValueError("compute_coll_coeffs needs CollData.init (host arrays); "
                         "use CollData.from_export, not CollData.device()")
    d, p = data.init, params
    vnmult = np.asarray(vnmult, np.float64).reshape(2)
    code_dt = float(code_dt)
    out = dict(code_dt=code_dt, vnmult=vnmult.copy())
    if not p.colls:
        return CollCoeffs(**out)
    if kwork_filter is None:
        kwork = d.kwork_filter
    elif isinstance(kwork_filter, (bool, np.bool_)):
        kwork = np.full(d.kwork_filter.shape, bool(kwork_filter))
    else:
        kwork = np.asarray(kwork_filter, bool)
    maps = le_maps(d.ixi_to_il, d.ixi_to_isgn, d.forbid, p.nx)
    # Every array is independent per ky: build in ky chunks to bound the host G7 temporaries.
    per_ky = p.nspec * p.negrid * p.nlambda * p.ntheta0 * 2 * p.nth
    chunk = max(1, min(p.naky, int(max_g_elements) // max(per_ky, 1)))
    parts = []
    for k0 in range(0, p.naky, chunk):
        ky = slice(k0, min(k0 + chunk, p.naky))
        pk = replace(p, naky=ky.stop - ky.start)
        parts.append(_coeffs_chunk(code_dt, vnmult, _slice_ky(d, ky), pk, kwork[ky], maps))
    for name in COEFF_NAMES:
        vals = [part.get(name) for part in parts]
        if vals[0] is not None:
            out[name] = vals[0] if len(vals) == 1 else np.concatenate(vals, axis=1)
    return CollCoeffs(**out)


def _slice_ky(d: CollInitData, ky: slice) -> CollInitData:
    """``CollInitData`` restricted to the ky indices ``ky``."""
    def sl(a, axis):
        if a is None:
            return None
        idx = [slice(None)] * a.ndim
        idx[axis] = ky
        return a[tuple(idx)]

    return d._replace(
        vnew=sl(d.vnew, 2), vnew_s=sl(d.vnew_s, 2), vnew_D=sl(d.vnew_D, 2), vnew_E=sl(d.vnew_E, 2),
        delvnew=sl(d.delvnew, 2), vnewh=sl(d.vnewh, 1), tunits=sl(d.tunits, 0),
        kperp2=sl(d.kperp2, 0), aj0=sl(d.aj0, 3), aj1=sl(d.aj1, 3), vpa=sl(d.vpa, 3),
        kwork_filter=sl(d.kwork_filter, 0), ieqzip=sl(d.ieqzip, 0), skip=sl(d.skip, 1))


def _coeffs_chunk(code_dt, vnmult, d: CollInitData, p: CollParams, kwork, maps) -> dict:
    """All ``CollCoeffs`` arrays for the (ky-sliced) data ``d``."""
    out = {}
    # Build-time skip mask of the solves: skip(is,ik), kwork_filter(it,ik), ieqzip(it,ik).
    # GS2 quirk: the conservation vectors of a kwork_filter mode are built without its solve, and the
    # filter state differs between the first init and a dt reset (collisions.fpp:3908-3909, 4058-4059;
    # fields.f90:376)
    apply = (~d.skip)[:, :, None] & (~kwork & ~d.ieqzip)[None, :, :]
    if p.has_lorentz:
        c1, beta, q = lorentz_matrices(code_dt, vnmult[0], d, p)
        out.update(c1le=c1, betaale=beta, qle=q)
        if p.conserve_moments:
            z0, s0, w0 = lorentz_conserve_vectors(code_dt, vnmult[0], (c1, beta, q), apply,
                                                  maps, d, p)
            out.update(z0le=z0, s0le=s0, w0le=w0)
    if p.has_diffuse:
        ec1, eb, eq = ediffuse_matrices(code_dt, vnmult[1], d, p, maps)
        out.update(ec1le=ec1, ebetaale=eb, eqle=eq)
        if p.conserve_moments:
            bz0, bs0, bw0 = diffuse_conserve_vectors(code_dt, vnmult[1], (ec1, eb, eq), apply,
                                                     maps, d, p)
            out.update(bz0le=bz0, bs0le=bs0, bw0le=bw0)
    return out


def lorentz_matrices(code_dt: float, vnm1: float, d: CollInitData, p: CollParams):
    """Lorentz tridiagonal LU arrays ``c1le, betaale, qle`` LE6 ``(S, Y, X, T, E, nxi_lim)``.

    get_lorentz_matrix (src/collisions.fpp:2392-2639, lorentz_scheme='default') and the le-layout
    forward elimination (2240-2285); spec 06 section 5.1.
    """
    S, E, Y, X, T, XI = p.nspec, p.negrid, p.naky, p.ntheta0, p.nth, p.nxi_lim
    ng2 = p.ng2
    c1 = np.zeros((S, Y, X, T, E, XI))
    beta = np.zeros_like(c1)
    q = np.zeros_like(c1)
    if p.lorentz_test:
        # vn = vnmult(1) * abs(spec(is)%vnewk) * tunits(ik)
        vn_sy = (vnm1 * np.abs(d.vnewk))[:, None] * d.tunits[None, :]        # (S, Y)
        vn_base = np.broadcast_to(vn_sy[:, :, None, None], (S, Y, X, E))
    else:
        vnc = np.transpose(vnm1 * d.vnew, (0, 2, 1))[:, :, None, :]         # (S, Y, 1, E)
    E_ = d.energy[None, None, None, :]
    for t in range(T):
        bm = d.bmag[t]
        je = int(d.jend[t])
        if je <= ng2 or (je == ng2 + 1 and p.special_wfb_lorentz):
            te, te2, teh = ng2, 2 * ng2, ng2
        else:
            te, te2, teh = je, 2 * je - 1, je - 1
        if p.lorentz_test:
            vn = vn_base
        else:
            vhyp = d.vnewh[:, :, :, t][..., None] if p.hyper_colls else 0.0   # (S, Y, X, 1)
            if p.hypermult:
                vn = vnc * (1.0 + vhyp)                    # vnmult(1)*vnew*(1 + vhyp)
            else:
                vn = vnc + vhyp                            # vnmult(1)*vnew + vhyp
            vn = np.broadcast_to(vn, (S, Y, X, E))
        s = np.sqrt(np.maximum(0.0, 1.0 - bm * d.al))   # safe_sqrt(1 - bmag*al(il)), (L,)
        kp = d.kperp2[:, :, t][None, :, :, None]         # (1, Y, X, 1)
        zb = bm * d.zstm
        den = (zb * zb)[:, None, None, None]             # (bmag*zstm)**2, (S, 1, 1, 1)
        pw = d.pitch_weights[t]

        def ee_of(slb1):
            # 0.5*energy(ie)*(1 + slb1**2)*kperp2*cfac / (bmag*zstm)**2
            return ((((0.5 * E_) * (1.0 + slb1 * slb1)) * kp) * p.cfac) / den

        two_vn_dt = (2.0 * vn) * code_dt
        aa = np.zeros((S, Y, X, E, XI))
        bb = np.zeros_like(aa)
        cc = np.zeros_like(aa)
        for il in range(2, te):                          # interior, Fortran il = 2..te-1
            slb0, slb1, slb2 = s[il - 2], s[il - 1], s[il]
            slbl = (slb1 + slb0) * 0.5
            slbr = (slb1 + slb2) * 0.5
            ee = ee_of(slb1)
            cc[..., il - 1] = (two_vn_dt * (1.0 - slbr * slbr)) / (pw[il - 1] * (slb2 - slb1))
            aa[..., il - 1] = (two_vn_dt * (1.0 - slbl * slbl)) / (pw[il - 1] * (slb1 - slb0))
            bb[..., il - 1] = (1.0 - (aa[..., il - 1] + cc[..., il - 1])) + ((ee * vn) * code_dt)
        # boundary at xi = 1
        slb1, slb2 = s[0], s[1]
        slbr = (slb1 + slb2) * 0.5
        ee = ee_of(slb1)
        cc[..., 0] = (two_vn_dt * (1.0 - slbr * slbr)) / (pw[0] * (slb2 - slb1))
        aa[..., 0] = 0.0
        bb[..., 0] = (1.0 - cc[..., 0]) + ((ee * vn) * code_dt)
        # boundary at xi = 0, il = te
        il = te
        slb0 = s[il - 2]
        if te == ng2:
            slb1 = s[il - 1]
            slb2 = -slb1
        else:
            slb1 = 0.0
            slb2 = -slb0
        slbl = (slb1 + slb0) * 0.5
        slbr = (slb1 + slb2) * 0.5
        ee = ee_of(slb1)
        deltaxi = pw[il - 1]
        if te == je:
            # GS2 quirk: doubled weight at the vpar=0 end point in the matrix but not in vpdiff
            # (collisions.fpp:2538 vs 4132-4136)
            deltaxi = 2 * deltaxi
        if (not p.special_wfb_lorentz) and deltaxi == 0 and il == ng2 + 1:
            deltaxi = d.xi[t, ng2 - 1] - d.xi[t, ng2 + 1]  # xi(ng2,ig) - xi(ng2+2,ig), by ixi
        cc[..., il - 1] = (two_vn_dt * (1.0 - slbr * slbr)) / (deltaxi * (slb2 - slb1))
        aa[..., il - 1] = (two_vn_dt * (1.0 - slbl * slbl)) / (deltaxi * (slb1 - slb0))
        bb[..., il - 1] = (1.0 - (aa[..., il - 1] + cc[..., il - 1])) + ((ee * vn) * code_dt)
        # symmetric fill (collisions.fpp:2630-2632): x(te+1:te2) = y(teh:1:-1)
        src = np.arange(teh - 1, -1, -1)
        dst = np.arange(te, te2)
        aa0, bb0, cc0 = aa.copy(), bb.copy(), cc.copy()
        aa[..., dst] = cc0[..., src]
        bb[..., dst] = bb0[..., src]
        cc[..., dst] = aa0[..., src]
        # le-layout elimination (collisions.fpp:2268-2282)
        c1[:, :, :, t] = cc
        qt = np.zeros_like(cc)
        bt = np.zeros_like(cc)
        bt[..., 0] = 1.0 / bb[..., 0]
        for k in range(1, te2):
            qt[..., k] = aa[..., k] * bt[..., k - 1]
            bt[..., k] = 1.0 / (bb[..., k] - qt[..., k] * cc[..., k - 1])
        q[:, :, :, t] = qt
        beta[:, :, :, t] = bt
    return c1, beta, q


def ediffuse_matrices(code_dt: float, vnm2: float, d: CollInitData, p: CollParams, maps):
    """Energy-diffusion LU arrays ``ec1le, ebetaale, eqle`` LE6 ``(S, Y, X, T, E, nxi_lim)``.

    get_ediffuse_matrix (src/collisions.fpp:2041-2182, ediff_scheme='default') and the le-layout
    elimination (1981-1999); forbidden rows stay 0. Spec 06 section 5.2.
    """
    S, E, Y, X, T, XI = p.nspec, p.negrid, p.naky, p.ntheta0, p.nth, p.nxi_lim
    xe, wM, capg = d.speed, d.w_maxw, d.capg_half
    le_il = maps.le_il[:, :XI]
    rows = maps.le_valid[:, :XI]                               # allowed rows (not forbid)
    ec1 = np.zeros((S, Y, X, T, E, XI))
    eb = np.zeros_like(ec1)
    eq = np.zeros_like(ec1)
    # vn = vnmult(2) * spec(is)%vnewk * tunits(ik)
    vn = ((vnm2 * d.vnewk)[:, None] * d.tunits[None, :])[:, :, None, None]     # (S, Y, 1, 1)
    c025 = (0.25 * vn) * code_dt
    vnew_s = np.transpose(d.vnew_s, (0, 2, 1))[:, :, None, None, :]           # (S, Y, 1, 1, E)
    for t in range(T):
        bm = d.bmag[t]
        slb1 = np.sqrt(np.maximum(0.0, 1.0 - bm * d.al[le_il[t]]))          # (XI,)
        kp = d.kperp2[:, :, t][None, :, :, None, None]                       # (1, Y, X, 1, 1)
        zb = bm * d.zstm
        den = (zb * zb)[:, None, None, None, None]
        # ee = 0.125*(1 - slb1**2)*vnew_s*kperp2*cfac / (bmag*zstm)**2   -> (S, Y, X, XI, E)
        # GS2 quirk: the classical term ee*code_dt is not scaled by vnmult(2), unlike Lorentz's
        # ee*vn*code_dt (collisions.fpp:2080 vs 2473)
        eeE = ((((0.125 * (1.0 - slb1 * slb1))[None, None, None, :, None] * vnew_s) * kp)
               * p.cfac) / den
        aa = np.zeros((S, Y, X, XI, E))
        bb = np.zeros_like(aa)
        cc = np.zeros_like(aa)
        for ie in range(1, E - 1):                           # Fortran ie = 2..negrid-1
            cc[..., ie] = -((c025 * capg[ie]) / (wM[ie] * (xe[ie + 1] - xe[ie])))
            aa[..., ie] = -((c025 * capg[ie - 1]) / (wM[ie] * (xe[ie] - xe[ie - 1])))
            bb[..., ie] = (1.0 - (aa[..., ie] + cc[..., ie])) + (eeE[..., ie] * code_dt)
        cc[..., 0] = -((c025 * capg[0]) / (wM[0] * (xe[1] - xe[0])))
        aa[..., 0] = 0.0
        bb[..., 0] = (1.0 - cc[..., 0]) + (eeE[..., 0] * code_dt)
        cc[..., E - 1] = 0.0
        aa[..., E - 1] = -((c025 * capg[E - 2]) / (wM[E - 1] * (xe[E - 1] - xe[E - 2])))
        bb[..., E - 1] = (1.0 - aa[..., E - 1]) + (eeE[..., E - 1] * code_dt)
        bt = np.zeros_like(aa)
        qt = np.zeros_like(aa)
        bt[..., 0] = 1.0 / bb[..., 0]
        for ie in range(E - 1):
            qt[..., ie + 1] = aa[..., ie + 1] * bt[..., ie]
            bt[..., ie + 1] = 1.0 / (bb[..., ie + 1] - qt[..., ie + 1] * cc[..., ie])
        m = rows[t][None, None, None, :, None]
        ec1[:, :, :, t] = np.where(m, cc, 0.0).transpose(0, 1, 2, 4, 3)
        eb[:, :, :, t] = np.where(m, bt, 0.0).transpose(0, 1, 2, 4, 3)
        eq[:, :, :, t] = np.where(m, qt, 0.0).transpose(0, 1, 2, 4, 3)
    return ec1, eb, eq


# ----- host versions of the runtime solves (used by the conservation rebuild, spec 06 section 7)
def _pad_xi_np(a, nx):
    if a.shape[-1] >= nx:
        return a
    pad = [(0, 0)] * (a.ndim - 1) + [(0, nx - a.shape[-1])]
    return np.pad(a, pad)


def np_solfp_lorentz(gle, lu, apply, d: CollInitData, p: CollParams):
    """Host transcription of solfp_lorentz_le_layout (src/collisions.fpp:3897-3980) on LE6 ``gle``
    (real or complex); lines with ``apply == False`` are skipped (skip / kwork_filter / ieqzip)."""
    c1, beta, q = (_pad_xi_np(np.asarray(a), gle.shape[-1]) for a in lu)
    out = np.array(gle, copy=True)
    E, ng2 = p.negrid, p.ng2
    for t in range(p.nth):
        cj = int(d.jend[t])
        je = max(2 * cj, 2 * ng2 + 1)
        nscatt = je - 1
        is_wfb = cj == ng2 + 1
        if is_wfb and p.special_wfb_lorentz:
            nscatt -= 1
        g = out[:, :, :, t, :E, :].copy()                    # (S, Y, X, E, NX)
        if je == 2 * cj and is_wfb and p.vpar_zero_mean:
            g[..., cj - 1] = (g[..., cj - 1] + g[..., 2 * cj - 1]) / 2.0
        g[..., je - 1] = 0.0
        if is_wfb and p.special_wfb_lorentz:
            gwfb = g[..., ng2].copy()
            g[..., ng2:je - 2] = g[..., ng2 + 1:je - 1].copy()
            g[..., je - 2] = 0.0
        qt, bt, ct = q[:, :, :, t], beta[:, :, :, t], c1[:, :, :, t]
        delta = np.zeros_like(g)
        delta[..., 0] = g[..., 0]
        for il in range(1, nscatt + 1):
            delta[..., il] = g[..., il] - qt[..., il] * delta[..., il - 1]
        g[..., nscatt] = delta[..., nscatt] * bt[..., nscatt]
        for il in range(nscatt - 1, -1, -1):
            g[..., il] = (delta[..., il] - ct[..., il] * g[..., il + 1]) * bt[..., il]
        if is_wfb and p.special_wfb_lorentz:
            g[..., ng2 + 1:je - 1] = g[..., ng2:je - 2].copy()
            g[..., ng2] = gwfb
        if cj != 0:
            g[..., 2 * cj - 1] = g[..., cj - 1]
        out[:, :, :, t, :E, :] = np.where(apply[:, :, :, None, None], g, out[:, :, :, t, :E, :])
    return out


def np_solfp_ediffuse(gle, lu, apply, d: CollInitData, p: CollParams):
    """Host transcription of solfp_ediffuse_le_layout (src/collisions.fpp:4037-4076)."""
    ec1, eb, eq = (_pad_xi_np(np.asarray(a), gle.shape[-1]) for a in lu)
    out = np.array(gle, copy=True)
    E = p.negrid
    for t in range(p.nth):
        mx = max(2 * int(d.jend[t]), 2 * p.ng2)
        g = out[:, :, :, t, :, :mx]                           # (S, Y, X, E+1, mx)
        delta = np.zeros_like(g[..., :E, :])
        delta[..., 0, :] = g[..., 0, :]
        for ie in range(E - 1):
            delta[..., ie + 1, :] = g[..., ie + 1, :] - eq[:, :, :, t, ie + 1, :mx] * delta[..., ie, :]
        new = np.zeros_like(g)
        for ie in range(E - 1, -1, -1):
            new[..., ie, :] = (delta[..., ie, :] - ec1[:, :, :, t, ie, :mx] * new[..., ie + 1, :]) \
                * eb[:, :, :, t, ie, :mx]
        out[:, :, :, t, :, :mx] = np.where(apply[:, :, :, None, None], new, g)
    return out


# ----- G7 broadcasting helpers for the conservation rebuild (G7 = (S, E, L, Y, X, 2, T))
def _g_sey(a):            # (S, E, Y)
    return a[:, :, None, :, None, None, None]


def _g_e(a):              # (E,)
    return a[None, :, None, None, None, None, None]


def _g_l(a):              # (L,)
    return a[None, None, :, None, None, None, None]


def _g_s(a):              # (S,)
    return a[:, None, None, None, None, None, None]


def _g_6(a):              # G6 (S, E, L, Y, X, T)
    return a[:, :, :, :, :, None, :]


def _g_vpdiff(a):         # (L, 2, T)
    return a[None, None, :, None, None, :, :]


def _g_kperp2(a):         # (Y, X, T)
    return a[None, None, None, :, :, None, :]


def _g_t(a):              # (T,)
    return a[None, None, None, None, None, None, :]


def _g_syxt(a):           # (S, Y, X, T)
    return a[:, None, None, :, :, None, :]


def _g_lt(a):             # (L, T)
    return a[None, None, :, None, None, None, :]


def integrate_moment_r34(gtmp, d: CollInitData, p: CollParams):
    """integrate_moment_r34 (src/le_grids.f90:1431-1527), 1 rank, no OpenMP: G7 real -> (S, Y, X, T).

    ``total(:,it,ik,is) = total + (w(ie,is)*wl(:,il))*(g(:,1,iglo) + g(:,2,iglo))`` in ascending
    iglo order: for fixed (it, ik, is) ie outer / il inner, or il outer / ie inner for 'yxels'
    (spec 06 section 5.5).
    """
    S, E, L = p.nspec, p.negrid, p.nlambda
    gsum = gtmp[:, :, :, :, :, 0, :] + gtmp[:, :, :, :, :, 1, :]       # (S, E, L, Y, X, T)
    total = np.zeros((S, p.naky, p.ntheta0, p.nth))
    order = [(ie, il) for ie in range(E) for il in range(L)] if p.ie_outer else \
        [(ie, il) for il in range(L) for ie in range(E)]
    for ie, il in order:
        fac = d.w[:, ie][:, None, None, None] * d.wl[il][None, None, None, :]   # (S, 1, 1, T)
        total = total + fac * gsum[:, ie, il]
    return total


def _inv_where(du):
    # where (abs(duinv) > epsilon(0.0)) duinv = 1./duinv   (collisions.fpp:856-859, 1263-1266)
    with np.errstate(divide="ignore"):
        return np.where(np.abs(du) > EPS, 1.0 / np.where(du == 0.0, 1.0, du), du)


def _solve_g(v, solve, lu, apply, maps, d, p):
    """``gather(g2le, v, ctmp); solve(ctmp); scatter(g2le, ctmp, v)`` on a real G7 vector.

    GS2 reuses one zeroed ``ctmp`` for the three solves; the solves write only allowed rows (and the
    zero padding row), so a fresh zeroed buffer per solve is identical.
    """
    le = gather_le(v, maps, p.nx)
    le = solve(le, lu, apply, d, p)
    return scatter_le(le, v, maps)


def _to_le(v, maps, p):
    return gather_le(v, maps, p.nx)[..., :p.nxi_lim]


def lorentz_conserve_vectors(code_dt, vnm1, lu, apply, maps, d: CollInitData, p: CollParams):
    """``z0le, s0le, w0le`` LE6 ``(S, Y, X, T, E+1, nxi_lim)``: init_lorentz_conserve
    (src/collisions.fpp:691-1194), conservative branch, le layout; spec 06 section 5.3."""
    E15 = d.e15
    vns1 = vnm1 * d.vnew_D                                          # vns(:,:,:,1), (S, E, Y)
    vpa, aj0, aj1 = d.vpa, _g_6(d.aj0), _g_6(d.aj1)
    vpdiff, V, En = _g_vpdiff(d.vpdiff), _g_e(d.speed), _g_e(d.energy)
    al, bmag, kp2 = _g_l(d.al), _g_t(d.bmag), _g_kperp2(d.kperp2)
    IM = lambda g: integrate_moment_r34(g, d, p)                  # noqa: E731
    solve = lambda v: _solve_g(v, np_solfp_lorentz, lu, apply, maps, d, p)  # noqa: E731
    if p.drag:
        # vns(ik,:,is,3) = vnmult(1)*spec(is)%vnewk*tunits(ik)/energy**1.5, electron-like only
        # GS2 quirk: this z0 build includes tunits, the runtime drag term does not
        # (collisions.fpp:740 vs 3073-3075)
        vns3 = np.zeros_like(vns1)
        ele = d.is_electron
        vns3[ele] = (((vnm1 * d.vnewk[ele])[:, None] * d.tunits[None, :])[:, None, :]
                     / E15[None, :, None])
        # z0 = -2.0*code_dt*vns3*vpdiff*speed*aj0
        z0 = -(((((2.0 * code_dt) * _g_sey(vns3)) * vpdiff) * V) * aj0)
        z0 = np.broadcast_to(z0, vpa.shape).copy()
        z0 = solve(z0)
        D = IM((vpa * aj0) * z0)                                       # v0z0
        z0 = z0 / (1.0 + _g_syxt(D))
    else:
        z0 = np.zeros(vpa.shape)
    du = IM(np.broadcast_to(((_g_sey(vns1) * vpa) * vpdiff) * V, vpa.shape))
    duinv = _g_syxt(_inv_where(du))
    s0 = -(((((_g_sey(vns1) * vpdiff) * V) * aj0) * code_dt) * duinv)
    s0 = solve(np.broadcast_to(s0, vpa.shape).copy())
    D = IM((vpa * aj0) * s0)                                           # v0s0
    s0 = s0 - _g_syxt(D) * z0
    D = IM(((((_g_sey(vns1) * V) * vpdiff) * aj0) * s0))              # v1s0
    s0 = s0 / (1.0 + _g_syxt(D))
    smz2 = _g_s(d.smz * d.smz)
    w0 = -((((((((_g_sey(vns1) * En) * al) * aj1) * code_dt) * smz2) * kp2) * duinv) / bmag)
    w0 = np.broadcast_to(w0, vpa.shape)
    if p.conserve_forbid_zero:
        w0 = np.where(_g_lt(d.forbid), 0.0, w0)
    w0 = solve(np.array(w0, copy=True))
    D = IM((vpa * aj0) * w0)                                           # v0w0
    w0 = w0 - z0 * _g_syxt(D)
    D = IM(((((_g_sey(vns1) * V) * vpdiff) * aj0) * w0))              # v1w1
    w0 = w0 - s0 * _g_syxt(D)
    D = IM((((_g_sey(vns1) * En) * al) * aj1) * w0)                   # v2w2
    w0 = w0 / (1.0 + _g_syxt(D))
    return _to_le(z0, maps, p), _to_le(s0, maps, p), _to_le(w0, maps, p)


def diffuse_conserve_vectors(code_dt, vnm2, lu, apply, maps, d: CollInitData, p: CollParams):
    """``bz0le, bs0le, bw0le`` LE6 ``(S, Y, X, T, E+1, nxi_lim)``: init_diffuse_conserve
    (src/collisions.fpp:1198-1679), le layout; spec 06 section 5.4."""
    vns1 = vnm2 * d.delvnew                                         # vns(:,:,:,1), (S, E, Y)
    vpa, aj0, aj1 = d.vpa, _g_6(d.aj0), _g_6(d.aj1)
    En, al, bmag, kp2 = _g_e(d.energy), _g_l(d.al), _g_t(d.bmag), _g_kperp2(d.kperp2)
    vnE = _g_sey(d.vnew_E)
    IM = lambda g: integrate_moment_r34(g, d, p)                  # noqa: E731
    solve = lambda v: _solve_g(v, np_solfp_ediffuse, lu, apply, maps, d, p)  # noqa: E731
    shape = vpa.shape
    # du = int(energy*vnmult(2)*vnew_E)   (collisions.fpp:1251)
    du = IM(np.broadcast_to((En * vnm2) * vnE, shape))
    duinv = _g_syxt(_inv_where(du))
    bz0 = np.broadcast_to(-((((code_dt * vnm2) * vnE) * aj0) * duinv), shape)
    if p.conserve_forbid_zero:
        bz0 = np.where(_g_lt(d.forbid), 0.0, bz0)
    bz0 = solve(np.array(bz0, copy=True))
    D = IM(((vnm2 * vnE) * aj0) * bz0)                                 # v0z0
    bz0 = bz0 / (1.0 + _g_syxt(D))
    du = IM(np.broadcast_to(_g_sey(vns1) * (vpa * vpa), shape))
    duinv = _g_syxt(_inv_where(du))
    bs0 = -((((_g_sey(vns1) * vpa) * aj0) * code_dt) * duinv)
    bs0 = solve(np.broadcast_to(bs0, shape).copy())
    D = IM(((vnm2 * vnE) * aj0) * bs0)                                 # v0s0
    bs0 = bs0 - _g_syxt(D) * bz0
    D = IM(((_g_sey(vns1) * vpa) * aj0) * bs0)                         # v1s0
    bs0 = bs0 / (1.0 + _g_syxt(D))
    smz2 = _g_s(d.smz * d.smz)
    bw0 = -((((((((_g_sey(vns1) * En) * al) * aj1) * code_dt) * smz2) * kp2) * duinv) / bmag)
    bw0 = np.broadcast_to(bw0, shape)
    if p.conserve_forbid_zero:
        bw0 = np.where(_g_lt(d.forbid), 0.0, bw0)
    bw0 = solve(np.array(bw0, copy=True))
    D = IM(((vnm2 * vnE) * aj0) * bw0)                                 # v0w0
    bw0 = bw0 - bz0 * _g_syxt(D)
    D = IM(((_g_sey(vns1) * vpa) * aj0) * bw0)                         # v1w1
    bw0 = bw0 - bs0 * _g_syxt(D)
    D = IM((((_g_sey(vns1) * En) * al) * aj1) * bw0)                   # v2w2
    bw0 = bw0 / (1.0 + _g_syxt(D))
    return _to_le(bz0, maps, p), _to_le(bs0, maps, p), _to_le(bw0, maps, p)
