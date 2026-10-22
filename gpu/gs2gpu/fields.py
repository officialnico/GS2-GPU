"""Field equations and the implicit field solve (spec 04; spec 10 G-10, G-17, G-23, G-27, G-28).

Mirrors, with GS2's association and summation order:

* ``integrate_species_original`` (src/le_grids.f90:723-786), ``getan_from_dfn`` / ``getan``
  (src/dist_fn.fpp:6273-6380), ``ensure_single_val_fields_pass`` (6201-6233);
* ``calculate_flux_surface_average`` (dist_fn.fpp:8006-8034) with the Klein-compensated
  trapezoid (externals/utils/integration.f90:40-71, summation.f90:71-103);
* ``getfieldeq1`` / ``getfieldeq`` (dist_fn.fpp:6575-6680);
* ``get_field_vector`` / ``getfield`` / ``get_field_solution`` (src/fields_implicit.f90:139-312)
  and the field update of ``advance_implicit`` (fields_implicit.f90:354-356);
* ``get_fields_direct_from_dfn`` / ``get_init_field`` (dist_fn.fpp:7457-7799);
* ``calculate_potentials_from_nonadiabatic_dfn`` (dist_fn.fpp:7322-7455).

Layouts (11_api.md section 3): g-like arrays are G7 ``(S, E, L, Y, X, 2, T)``, fields F3
``(Y, X, T)``, ``(it, ik)`` arrays K2 ``(Y, X)``. The supercell / class structure of the
implicit solve (spec 04 section 3.11) lives in :class:`Topology`: per class ``c`` the member
tables ``ik[c]`` ``(M_c,)`` and ``it[c]`` ``(M_c, N_c)`` (0-based), and the response matrices
(``gs2gpu.response.Aminv``) are ``(M_c, R_c, R_c)`` in GS2's runtime "B" order
``B(ig, f, n) = igp + nth*(f + nf*n)`` for both axes (``[m, output, input]``).

Exact mode (``FieldParams.exact``, default) keeps every sum in GS2's order: the velocity
integral is a sequential scan over velocity points in g_lo layout order, the trapezoids use
the Klein sum, and ``u = -aminv.fl`` accumulates ``u = ((0 - S_1) - S_2) ...`` with each
``S_n`` a sequential sum over the ``nidx`` entries of cell ``n``. ``exact=False`` uses tree
reductions (``jnp.sum``/``einsum``), which differ at the 1e-16 level. In numerics strict mode
(no FMA contraction) ``getfield`` reproduces GS2's reference build bitwise.

Usage (integration layer; all kernels jit with ``params`` static)::

    params = FieldParams.from_export(ex)
    data = to_jax(FieldData.from_export(ex, params))
    aminv = to_jax(response.Aminv.from_export(ex))          # implicit or local (B layout)
    old = Fields(phi, apar, bpar)                            # fields_implicit.f90:338
    u = getfield(old, gnew_pass1, aminv, data, params)       # increment (after_getfield)
    new = field_update(u, old, params)                       # after_field_update
    init = get_init_field(gnew, data, params)                # t = 0 and after a dt reset
    aminv = response.build_aminv(timeadv0, data, params, code_dt=dt)   # after a dt change
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import NamedTuple, Optional

import jax
import jax.numpy as jnp
import numpy as np
from jax import lax

from .io.export import Export
from .numerics import cmul, ob, rmul
from .types import Fields

#: Fortran ``epsilon(0.0)`` under ``-fdefault-real-8``.
EPS = float(np.finfo(np.float64).eps)

# Switch values (src/dist_fn.fpp:124-138, src/fields.f90:35).
ADIABATIC_DEFAULT, ADIABATIC_ZERO, ADIABATIC_FIELDLINEAVG, ADIABATIC_YAVG = 1, 2, 3, 4
BOUNDARY_DEFAULT, BOUNDARY_ZERO, BOUNDARY_SELF_PERIODIC, BOUNDARY_ALTERNATE_ZERO, \
    BOUNDARY_LINKED = 0, 1, 2, 3, 4
FIELDOPT_IMPLICIT, FIELDOPT_TEST, FIELDOPT_LOCAL, FIELDOPT_GF_LOCAL = 1, 2, 3, 4
TRACER_SPECIES = 3


# --------------------------------------------------------------------------- export helpers
def _opt(ex: Export, path: str):
    """Value at ``path`` or None if GS2 did not allocate it / the exporter did not write it."""
    try:
        return ex[path]
    except KeyError:
        return None


def _first(ex: Export, *paths, default=None):
    for path in paths:
        value = _opt(ex, path)
        if value is not None:
            return value
    return default


def _scalar(value, kind=float):
    if value is None:
        return None
    arr = np.asarray(value)
    return kind(arr.reshape(-1)[0]) if arr.size else None


def _num(value, default, kind=float):
    """Scalar ``value`` or ``default`` if it is missing (0 is a valid value)."""
    v = _scalar(value, kind)
    return default if v is None else v


def _bool(value, default=False) -> bool:
    if value is None:
        return default
    return bool(np.asarray(value).reshape(-1)[0])


def _req(ex: Export, *paths):
    value = _first(ex, *paths)
    if value is None:
        raise KeyError(f"none of {paths} found in {ex.path}")
    return value


def species_group(ex: Export):
    """``/species`` (spec 09), ``/species_`` (loader alias) or ``/species_data`` (the exporter's
    name: a netCDF-4 group may not share the name of the root dimension ``species``)."""
    for path in ("/species", "/species_", "/species_data"):
        grp = _opt(ex, path)
        if grp is not None and hasattr(grp, "arrays"):
            return grp
    raise KeyError(f"{ex.path}: neither /species nor /species_data")


def _attr(ex: Export, name: str, groups=("/gs2_layouts",)):
    if name in ex.attrs:
        return ex.attrs[name]
    for group in groups:
        grp = _opt(ex, group)
        if grp is not None and name in getattr(grp, "attrs", {}):
            return grp.attrs[name]
    return None


# --------------------------------------------------------------------------- topology (host)
class Topology(NamedTuple):
    """Supercell classes of the implicit solve (fields_implicit.f90:1074-1169).

    ``ik[c]``: ``(M_c,)`` int32 0-based ky index of member ``m`` of class ``c``;
    ``it[c]``: ``(M_c, N_c)`` int32 0-based kx index of cell ``n`` (``n = 0`` leftmost);
    ``zonal[c]``: ``(M_c, N_c)`` bool, ``aky == 0 .and. akx == 0`` exactly (the members whose
    inverse rows GS2 zeroes, fields_implicit.f90:834-857);
    ``empty[c]``: ``(M_c, N_c)`` bool, ``|aky| < eps .and. |akx| < eps`` (cells field_option
    'local' ignores, fields_local.f90:580);
    ``perm``: ``(Y, X)`` int32, position of ``(ik, it)`` in the concatenation over
    ``c, m, n`` (C order) of all cells.
    """

    ik: tuple
    it: tuple
    zonal: tuple
    empty: tuple
    perm: np.ndarray

    @property
    def nclass(self) -> int:
        return len(self.ik)

    def sizes(self, c: int) -> tuple[int, int]:
        """``(M_c, N_c)`` (static)."""
        return tuple(int(s) for s in np.shape(self.it[c]))


def class_tables_from_links(naky: int, ntheta0: int, linked: bool, N_class, l_links=None,
                            r_links=None, itright=None):
    """Member tables ``f_lo(c)%ik(m,n)``, ``f_lo(c)%it(m,n)`` (0-based) per class.

    Mirrors ``finish_fields_layouts`` + ``kt2ki`` (src/fields_implicit.f90:1077-1169): members
    are found looping ``it`` outer, ``ik`` inner, at the leftmost cell (``n == 1``), and the
    chain is followed with ``itright``. Non-linked: one class with ``N = 1`` and
    ``m = ik + naky*(it-1)``. ``l_links``, ``r_links``, ``itright`` are NumPy ``(Y, X)``
    with Fortran values.
    """
    N_class = [int(n) for n in np.asarray(N_class).reshape(-1)]
    out_ik, out_it = [], []
    if not linked:
        ik_list, it_list = [], []
        for it in range(ntheta0):
            for ik in range(naky):
                ik_list.append(ik)
                it_list.append([it])
        return [np.asarray(ik_list, np.int32)], [np.asarray(it_list, np.int32)]
    l_links = np.asarray(l_links)
    r_links = np.asarray(r_links)
    itright = np.asarray(itright)
    for c, nc in enumerate(N_class):
        ik_list, it_list = [], []
        for it in range(ntheta0):
            for ik in range(naky):
                nn = 1 + int(l_links[ik, it]) + int(r_links[ik, it])
                n = 1 + int(l_links[ik, it])
                if nn not in N_class:
                    raise ValueError(f"kt2ki: supercell size {nn} of (ik={ik + 1}, it={it + 1}) "
                                     f"not in N_class {N_class}")
                if N_class.index(nn) != c or n != 1:
                    continue
                chain = [it]
                itr = it
                for _ in range(1, nc):
                    itr = int(itright[ik, itr]) - 1
                    chain.append(itr)
                ik_list.append(ik)
                it_list.append(chain)
        out_ik.append(np.asarray(ik_list, np.int32).reshape(-1))
        out_it.append(np.asarray(it_list, np.int32).reshape(-1, nc))
    return out_ik, out_it


def make_topology(ik_tables, it_tables, aky, akx) -> Topology:
    """Build :class:`Topology` from 0-based member tables (host, NumPy)."""
    aky = np.asarray(aky, np.float64)
    akx = np.asarray(akx, np.float64)
    naky, ntheta0 = aky.size, akx.size
    perm = np.full((naky, ntheta0), -1, np.int64)
    zonal, empty = [], []
    pos = 0
    for ikc, itc in zip(ik_tables, it_tables):
        ikc = np.asarray(ikc, np.int32)
        itc = np.asarray(itc, np.int32)
        M, N = itc.shape
        kyv = aky[ikc][:, None]
        kxv = akx[itc]
        zonal.append((kyv == 0.0) & (kxv == 0.0))
        # GS2 quirk: fields_local tests |aky| < eps .and. |akx| < eps (fields_local.f90:580),
        # fields_implicit aky == 0 .and. akx == 0 exactly (fields_implicit.f90:834)
        empty.append((np.abs(kyv) < EPS) & (np.abs(kxv) < EPS))
        for m in range(M):
            for n in range(N):
                if perm[ikc[m], itc[m, n]] != -1:
                    raise ValueError(f"cell (ik={ikc[m]}, it={itc[m, n]}) in two supercells")
                perm[ikc[m], itc[m, n]] = pos
                pos += 1
    if np.any(perm < 0):
        raise ValueError("supercell tables do not cover every (ik, it)")
    return Topology(ik=tuple(np.asarray(a, np.int32) for a in ik_tables),
                    it=tuple(np.asarray(a, np.int32) for a in it_tables),
                    zonal=tuple(zonal), empty=tuple(empty), perm=perm.astype(np.int32))


def topology_from_export(ex: Export, *, prefer_export: bool = True) -> Topology:
    """Class structure from ``/fields_implicit/cls_ik|cls_it`` (spec 09 section 7.4), or
    recomputed from ``/dist_fn`` links (``N_class``, ``l_links``, ``r_links``, ``itright``)."""
    aky = np.asarray(ex["/kt_grids/aky"], np.float64).reshape(-1)
    akx = np.asarray(ex["/kt_grids/akx"], np.float64).reshape(-1)
    N_class = np.asarray(_req(ex, "/dist_fn/N_class")).reshape(-1)
    M_class = np.asarray(_req(ex, "/dist_fn/M_class")).reshape(-1)
    cls_ik = _opt(ex, "/fields_implicit/cls_ik") if prefer_export else None
    cls_it = _opt(ex, "/fields_implicit/cls_it") if prefer_export else None
    if cls_ik is not None and cls_it is not None:
        cls_ik = np.asarray(cls_ik)
        cls_it = np.asarray(cls_it)
        ik_tables, it_tables = [], []
        for c, (nc, mc) in enumerate(zip(N_class, M_class)):
            nc, mc = int(nc), int(mc)
            ik_tables.append(cls_ik[c, 0, :mc].astype(np.int32) - 1)
            it_tables.append(cls_it[c, :nc, :mc].T.astype(np.int32) - 1)
        return make_topology(ik_tables, it_tables, aky, akx)
    bsw = _scalar(_first(ex, "/dist_fn/boundary_option_switch"), int)
    linked = bsw == BOUNDARY_LINKED
    if linked:
        ik_tables, it_tables = class_tables_from_links(
            aky.size, akx.size, True, N_class, ex["/dist_fn/l_links"], ex["/dist_fn/r_links"],
            ex["/dist_fn/itright"])
    else:
        ik_tables, it_tables = class_tables_from_links(aky.size, akx.size, False, N_class)
    for c, mc in enumerate(M_class):
        if ik_tables[c].size != int(mc):
            raise ValueError(f"class {c + 1}: {ik_tables[c].size} members found, M_class={mc}")
    return make_topology(ik_tables, it_tables, aky, akx)


# --------------------------------------------------------------------------- params / data
def _velocity_order(ex: Export) -> np.ndarray:
    """Flat ``(is, ie, il)`` indices (C order of ``(S, E, L)``) in the order
    ``integrate_species_original`` adds them for a fixed ``(it, ik)`` (le_grids.f90:747-760):
    increasing ``iglo`` of the layout (gs2_layouts.fpp:917-948)."""
    S, E, L = ex.nspec, ex.negrid, ex.nlambda
    iglo_of = _opt(ex, "/gs2_layouts/iglo_of")
    if iglo_of is not None:
        key = np.asarray(iglo_of)[:, :, :, 0, 0].reshape(-1)
        return np.argsort(key, kind="stable").astype(np.int32)
    layout = str(_attr(ex, "layout") or _first(ex, "/gs2_layouts/layout") or "lxyes").strip()
    return velocity_order_from_layout(layout, S, E, L)


def velocity_order_from_layout(layout: str, S: int, E: int, L: int) -> np.ndarray:
    """Velocity summation order from a layout string (first character fastest)."""
    idx = np.indices((S, E, L)).reshape(3, -1)
    key = {"s": idx[0], "e": idx[1], "l": idx[2]}
    fast_to_slow = [ch for ch in layout if ch in "sel"]
    if sorted(fast_to_slow) != ["e", "l", "s"]:
        raise ValueError(f"bad layout {layout!r}")
    return np.lexsort(tuple(key[ch] for ch in fast_to_slow)).astype(np.int32)


def _amp_dependent_hyper(ex: Export) -> bool:
    """True if ``hyper_diff`` inside ``timeadv(..., 0)`` couples modes through |phi|^2
    (iso_shear with const_amp=F, or aniso_shear; hyper.f90:658-784; spec 04 P9)."""
    hyp = _opt(ex, "/hyper")
    if hyp is None:
        return False
    if not _bool(hyp.get("hyper_on"), False):
        return False
    if not _bool(_first(ex, "/dist_fn/hyper_in_initialisation"), True):
        return False
    sw = _scalar(hyp.get("hyper_option_switch"), int)
    if sw is None:
        return True  # unknown: be conservative
    # hyper_option_{none,visc,res,both,simple3D} = 1..5 (hyper.f90:33-37)
    if sw in (1, 5):
        return False
    dvisc = _scalar(hyp.get("D_hypervisc"), float)
    if dvisc is not None and dvisc < 0.0:
        return False
    iso = _bool(hyp.get("isotropic_shear"), True)
    const = _bool(hyp.get("const_amp"), False)
    return (not iso) or (not const)


@dataclass(frozen=True)
class FieldParams:
    """Static switches and scalars of the field solve (hashable; jit static)."""

    ntgrid: int
    nth: int
    naky: int
    ntheta0: int
    nspec: int
    negrid: int
    nlambda: int
    has_phi: bool
    has_apar: bool
    has_bpar: bool
    beta: float
    tite: float
    apfac: float
    esv: bool
    linked: bool
    adiabatic_option: int
    has_electron_species: bool
    has_ion_species: bool
    phi_weight: float                 # calculate_potentials_from_nonadiabatic_dfn (7360-7366)
    field_option: str = "implicit"    # 'implicit' or 'local'
    amp_dependent_hyper: bool = False
    exact: bool = True
    unroll: int = 8

    @property
    def nfield(self) -> int:
        return int(self.has_phi) + int(self.has_apar) + int(self.has_bpar)

    @property
    def nidx(self) -> int:
        return self.nth * self.nfield

    @property
    def fl_avg_on(self) -> bool:
        """``.not. has_electron_species .and. adiabatic_option_switch == fieldlineavg``
        (dist_fn.fpp:6605-6609, 7637)."""
        # GS2 quirk: fl_avg needs "no electron-like species" only, while the +tite adiabatic
        # term of gamtot also applies when there are no ions (dist_fn.fpp:6557 vs 6605)
        return (not self.has_electron_species) and self.adiabatic_option == ADIABATIC_FIELDLINEAVG

    @property
    def slots(self) -> tuple[int, ...]:
        """Indices (0 phi, 1 apar, 2 bpar) of the present fields in solve order."""
        return tuple(i for i, on in enumerate((self.has_phi, self.has_apar, self.has_bpar)) if on)

    @classmethod
    def from_export(cls, ex: Export, *, exact: bool = True, unroll: int = 8) -> "FieldParams":
        rp = ex["/run_parameters"]
        has_phi = _bool(rp["has_phi"])
        has_apar = _bool(rp["has_apar"])
        has_bpar = _bool(rp["has_bpar"])
        afilter = _num(_first(ex, "/dist_fn/afilter"), 0.0)
        if afilter > EPS:
            raise NotImplementedError("afilter > 0 (gated: spec 04 section 4)")
        if _bool(_first(ex, "/dist_fn/gf_lo_integrate"), False):
            raise NotImplementedError("gf_lo_integrate = .true. (gated: spec 04 section 4)")
        fsw = _scalar(_first(ex, "/fields/fieldopt_switch"), int)
        if fsw is None:
            fsw = FIELDOPT_LOCAL if _opt(ex, "/fields_local") is not None \
                and _opt(ex, "/fields_implicit/aminv_c1") is None else FIELDOPT_IMPLICIT
        if fsw == FIELDOPT_IMPLICIT:
            field_option = "implicit"
        elif fsw == FIELDOPT_LOCAL:
            field_option = "local"
        else:
            raise NotImplementedError(f"field_option switch {fsw} (gf_local/test are gated)")
        bsw = _scalar(_first(ex, "/dist_fn/boundary_option_switch"), int)
        linked = bsw == BOUNDARY_LINKED
        esv = _bool(_first(ex, "/dist_fn/esv"), False) and linked   # dist_fn.fpp:1327
        spec = species_group(ex)
        stype = np.asarray(spec["spec_type"]).reshape(-1)
        has_e = spec.get("has_electron_species")
        has_e = _bool(has_e) if has_e is not None else bool(np.any((stype == 2) | (stype == 4)))
        has_i = spec.get("has_ion_species")
        has_i = _bool(has_i) if has_i is not None else bool(np.any(stype == 1))
        tite = float(_scalar(rp["tite"], float))
        dens = np.asarray(spec["spec_dens"], np.float64).reshape(-1)
        z = np.asarray(spec["spec_z"], np.float64).reshape(-1)
        temp = np.asarray(spec["spec_temp"], np.float64).reshape(-1)
        total = 0.0
        for d_, z_, t_ in zip(dens, z, temp):        # sum(spec%dens*spec%z*spec%z/spec%temp)
            total = total + ((float(d_) * float(z_)) * float(z_)) / float(t_)
        phi_weight = 1.0 / (total + tite) if ((not has_e) or (not has_i)) else 1.0 / total
        adiab = _num(_first(ex, "/dist_fn/adiabatic_option_switch"), ADIABATIC_DEFAULT, int)
        return cls(
            ntgrid=ex.ntgrid, nth=ex.nth, naky=ex.naky, ntheta0=ex.ntheta0, nspec=ex.nspec,
            negrid=ex.negrid, nlambda=ex.nlambda, has_phi=has_phi, has_apar=has_apar,
            has_bpar=has_bpar, beta=float(_scalar(rp["beta"], float)), tite=tite,
            apfac=float(_num(_first(ex, "/dist_fn/apfac"), 1.0)),
            esv=esv, linked=linked, adiabatic_option=adiab, has_electron_species=has_e,
            has_ion_species=has_i, phi_weight=float(phi_weight), field_option=field_option,
            amp_dependent_hyper=_amp_dependent_hyper(ex), exact=exact, unroll=unroll)


class FieldData(NamedTuple):
    """Arrays of the field solve (NumPy from ``from_export``; ``to_jax`` before use).

    Products that GS2 forms from static arrays inside the field routines are precomputed here
    with the same IEEE operations (host float64): ``fac2_*[v, ig] = (wgt(is)*w(ie,is))*wl(ig,il)``
    (le_grids.f90:759-760, ``wgt`` with tracer species zeroed, 741), ``gamtot_gf =
    gamtot*gridfac1``, ``kperp2_gf = kperp2*gridfac1`` (dist_fn.fpp:6621, 6639), ``bfac =
    (beta*apfac)/(bmag*bmag)`` (6658), the get_fields_direct_from_dfn factors (7655, 7746) and
    ``bpar_weight = -beta/bmag**2`` (7422).
    """

    aj0: np.ndarray            # G6
    aj1: np.ndarray            # G6
    vperp2: np.ndarray         # G6
    vpa: np.ndarray            # G7
    fac2_phi: np.ndarray       # (V, T), V = S*E*L flattened (is, ie, il)
    fac2_apar: np.ndarray
    fac2_bpar: np.ndarray
    vel_order: np.ndarray      # (V,) int32 summation order
    kwork_filter: np.ndarray   # K2 bool (value after init_fields)
    gamtot: np.ndarray         # F3
    gamtot1: np.ndarray
    gamtot2: np.ndarray
    gridfac1: np.ndarray
    gamtot_gf: np.ndarray
    kperp2_gf: np.ndarray
    inv_kperp2: np.ndarray
    bfac: np.ndarray           # (T,)
    bpar_weight: np.ndarray    # (T,)
    inv_phi_denom: np.ndarray  # F3
    inv_bpar_denom: np.ndarray
    phi_num_a: np.ndarray      # beta*gamtot2 + bmag**2
    phi_num_b: np.ndarray      # beta*gamtot1
    bpar_num_a: np.ndarray     # beta*gamtot
    bpar_num_b: np.ndarray     # (beta/2.0)*gamtot1
    dtheta: np.ndarray         # (T-1,) theta(2:) - theta(:n-1)
    jacob: np.ndarray          # (T,)
    awgt: np.ndarray           # K2
    esv_ik: np.ndarray         # (P,) int32 pairs of ensure_single_val_fields_pass
    esv_src: np.ndarray
    esv_dst: np.ndarray
    topo: Topology

    @classmethod
    def from_export(cls, ex: Export, params: Optional[FieldParams] = None,
                    topo: Optional[Topology] = None) -> "FieldData":
        p = params or FieldParams.from_export(ex)
        S, E, L, Y, X, T = ex.nspec, ex.negrid, ex.nlambda, ex.naky, ex.ntheta0, ex.nth
        arrays = ex["/dist_fn_arrays"]
        dfn = ex["/dist_fn"]
        kt = ex["/kt_grids"]
        th = ex["/theta_grid"]
        spec = species_group(ex)
        f64 = lambda a: np.asarray(a, np.float64)  # noqa: E731
        aj0, aj1, vperp2, vpa = (f64(arrays[n]) for n in ("aj0", "aj1", "vperp2", "vpa"))
        w = f64(ex["/le_grids/w"]).reshape(S, E)
        wl = f64(ex["/le_grids/wl"]).reshape(L, T)
        z = f64(spec["spec_z"]).reshape(-1)
        dens = f64(spec["spec_dens"]).reshape(-1)
        temp = f64(spec["spec_temp"]).reshape(-1)
        mass = f64(spec["spec_mass"]).reshape(-1)
        tracer = np.asarray(spec["spec_type"]).reshape(-1) == TRACER_SPECIES

        def fac2(wgt):
            wgt = np.where(tracer, 0.0, wgt)                     # le_grids.f90:741
            fac = wgt[:, None] * w                               # weights(is)*w(ie,is)
            return (fac[:, :, None, None] * wl[None, None, :, :]).reshape(S * E * L, T)

        fac2_phi = fac2(z * dens)                                             # 6313
        fac2_apar = fac2((((2.0 * p.beta) * z) * dens) * np.sqrt(temp / mass))  # 6335
        fac2_bpar = fac2(temp * dens)                                         # 6355
        bmag = f64(th["bmag"]).reshape(T)
        theta = f64(th["theta"]).reshape(T)
        jacob = f64(th["jacob"]).reshape(T)
        gamtot, gamtot1, gamtot2 = (f64(dfn[n]) for n in ("gamtot", "gamtot1", "gamtot2"))
        gridfac1 = f64(dfn["gridfac1"])
        kperp2 = f64(kt["kperp2"])
        inv_kperp2 = f64(kt["inv_kperp2"])
        beta = p.beta
        bmag2 = bmag * bmag
        phi_num_a = beta * gamtot2 + bmag2
        phi_num_b = beta * gamtot1
        bpar_num_a = beta * gamtot
        bpar_num_b = (beta / 2.0) * gamtot1
        inv_phi = _opt(ex, "/dist_fn/inv_phi_denominator_g")
        inv_bpar = _opt(ex, "/dist_fn/inv_bpar_denominator_g")
        if inv_phi is None or inv_bpar is None:
            cphi, cbpar = denominators(gamtot, gamtot1, gamtot2, bmag, beta,
                                       p.has_phi, p.has_bpar)
            inv_phi = cphi if inv_phi is None else inv_phi
            inv_bpar = cbpar if inv_bpar is None else inv_bpar
        # GS2 quirk: awgt is computed at the first calculate_flux_surface_average call and then
        # frozen, with the kwork_filter state of that call (dist_fn.fpp:8006-8020)
        awgt = _opt(ex, "/dist_fn/awgt")
        if awgt is None:
            if p.fl_avg_on:
                gamtot3 = _opt(ex, "/dist_fn/gamtot3")
                if gamtot3 is None:
                    gamtot3 = gamtot3_from_gamtot(gamtot, p.tite)
                awgt = awgt_host(theta, jacob, f64(gamtot3), f64(kt["aky"]).reshape(-1))
            else:
                awgt = np.zeros((Y, X))
        esv_ik, esv_src, esv_dst = esv_pairs(_opt(ex, "/dist_fn/itright") if p.esv else None,
                                             Y, X)
        if topo is None:
            topo = topology_from_export(ex)
        return cls(
            aj0=aj0, aj1=aj1, vperp2=vperp2, vpa=vpa,
            fac2_phi=fac2_phi, fac2_apar=fac2_apar, fac2_bpar=fac2_bpar,
            vel_order=_velocity_order(ex),
            kwork_filter=np.asarray(kt["kwork_filter"], bool).reshape(Y, X),
            gamtot=gamtot, gamtot1=gamtot1, gamtot2=gamtot2, gridfac1=gridfac1,
            gamtot_gf=gamtot * gridfac1, kperp2_gf=kperp2 * gridfac1, inv_kperp2=inv_kperp2,
            bfac=(beta * p.apfac) / bmag2, bpar_weight=-(beta / bmag2),
            inv_phi_denom=f64(inv_phi), inv_bpar_denom=f64(inv_bpar),
            phi_num_a=phi_num_a, phi_num_b=phi_num_b, bpar_num_a=bpar_num_a,
            bpar_num_b=bpar_num_b, dtheta=theta[1:] - theta[:-1], jacob=jacob,
            awgt=f64(awgt).reshape(Y, X), esv_ik=esv_ik, esv_src=esv_src, esv_dst=esv_dst,
            topo=topo)


def denominators(gamtot, gamtot1, gamtot2, bmag, beta, has_phi, has_bpar):
    """``setup_get_fields_direct_from_dfn_denominators`` (dist_fn.fpp:7457-7585), host."""
    bmag2 = np.asarray(bmag) * np.asarray(bmag)

    def inv(d):
        d = np.asarray(d, np.float64)
        out = np.zeros_like(d)
        ok = np.abs(d) > EPS
        np.divide(1.0, d, out=out, where=ok)
        return out

    shape = np.shape(gamtot)
    phi_d = np.zeros(shape)
    bpar_d = np.zeros(shape)
    if has_phi:
        if has_bpar:
            phi_d = inv((beta * gamtot2 + bmag2) * gamtot + ((beta / 2.0) * gamtot1) * gamtot1)
        else:
            phi_d = inv(np.array(gamtot, np.float64, copy=True))
    if has_bpar:
        if has_phi:
            bpar_d = inv(gamtot * (beta * gamtot2 + bmag2) + ((beta / 2.0) * gamtot1) * gamtot1)
        else:
            bpar_d = inv(beta * gamtot2 + bmag2)
    return phi_d, bpar_d


def gamtot3_from_gamtot(gamtot, tite):
    """``gamtot3 = (gamtot - tite)/gamtot``, ``< 2*eps -> 1`` (dist_fn.fpp:6564-6565),
    from the final (tite-including) gamtot."""
    g3 = (np.asarray(gamtot) - tite) / np.asarray(gamtot)
    return np.where(g3 < 2.0 * EPS, 1.0, g3)


def klein_sum_np(values) -> np.ndarray:
    """``klein_sum`` (externals/utils/summation.f90:71-103) over the last axis, NumPy."""
    a = np.moveaxis(np.asarray(values, np.float64), -1, 0)
    s = np.zeros(a.shape[1:])
    c1 = np.zeros(a.shape[1:])
    c2 = np.zeros(a.shape[1:])
    for x in a:
        t = s + x
        lc1 = np.where(np.abs(s) >= np.abs(x), (s - t) + x, (x - t) + s)
        s = t
        t = c1 + lc1
        lc2 = np.where(np.abs(c1) >= np.abs(lc1), (c1 - t) + lc1, (lc1 - t) + c1)
        c1 = t
        c2 = c2 + lc2
    return (s + c1) + c2


def trapz_real_np(x, y) -> np.ndarray:
    """``real_trapezoidal_integration`` over the last axis (integration.f90:55-61), NumPy."""
    x = np.asarray(x, np.float64)
    y = np.asarray(y, np.float64)
    v = (x[1:] - x[:-1]) * (y[..., :-1] + y[..., 1:])
    return klein_sum_np(v) * 0.5


def awgt_host(theta, jacob, gamtot3, aky) -> np.ndarray:
    """First-call ``awgt`` (dist_fn.fpp:8006-8020) with an all-false kwork_filter (spec 07 P16)."""
    Y, X, _ = np.shape(gamtot3)
    awgt = np.zeros((Y, X))
    for ik in range(Y):
        if aky[ik] > EPS:
            continue
        awgt[ik] = 1.0 / trapz_real_np(theta, jacob[None, :] * np.asarray(gamtot3)[ik])
    return awgt


def esv_pairs(itright, naky, ntheta0):
    """Copy list of ``ensure_single_val_fields_pass`` (dist_fn.fpp:6201-6233):
    ``arr(-ntgrid, itright(it,ik), ik) = arr(ntgrid, it, ik)``, loop ik outer, it inner;
    duplicate destinations keep the last writer."""
    if itright is None:
        e = np.zeros(0, np.int32)
        return e, e, e
    itright = np.asarray(itright)
    dest = {}
    for ik in range(naky):
        for it in range(ntheta0):
            link = int(itright[ik, it])
            if link < 0:
                continue
            dest[(ik, link - 1)] = it
    keys = sorted(dest)
    return (np.asarray([k[0] for k in keys], np.int32), np.asarray([dest[k] for k in keys], np.int32),
            np.asarray([k[1] for k in keys], np.int32))


# --------------------------------------------------------------------------- kernels
def _zeros_field(shape):
    return jnp.zeros(shape, jnp.complex128)


def _ensure_single_val(arr, data: FieldData):
    """``ensure_single_val_fields_pass(arr, force_all=.true.)`` (dist_fn.fpp:6201-6233)."""
    if np.size(data.esv_ik) == 0:
        return arr
    return arr.at[data.esv_ik, data.esv_dst, 0].set(arr[data.esv_ik, data.esv_src, -1])


def _velocity_integrals(g, data: FieldData, params: FieldParams, want, kwork_filter=None):
    """``integrate_species_original`` (le_grids.f90:723-786) for the g_work of getan_from_dfn
    (dist_fn.fpp:6290-6358), all requested fields in one pass over velocity space.

    ``want`` = (phi, apar, bpar) booleans. Per velocity point (is, ie, il) and field:
    ``s = g_work(:,1) + g_work(:,2)`` with ``g_work = aj0*g`` | ``(aj0*vpa)*g`` |
    ``(aj1*vperp2)*g``; ``term = fac2*s`` (componentwise); ``total = total + term``,
    sequential in g_lo order. Filtered modes are exactly ``+0`` (751-753).
    """
    S, E, L, Y, X, _, T = g.shape
    V = S * E * L
    gv = g.reshape((V, Y, X, 2, T))
    aj0 = jnp.reshape(data.aj0, (V, Y, X, T))
    aj1 = jnp.reshape(data.aj1, (V, Y, X, T))
    vperp2 = jnp.reshape(data.vperp2, (V, Y, X, T))
    vpa = jnp.reshape(data.vpa, (V, Y, X, 2, T))
    f_phi = jnp.asarray(data.fac2_phi)
    f_apar = jnp.asarray(data.fac2_apar)
    f_bpar = jnp.asarray(data.fac2_bpar)
    want_phi, want_apar, want_bpar = want

    def terms(v):
        gg = gv[v]
        g1, g2 = gg[..., 0, :], gg[..., 1, :]
        out = []
        if want_phi:
            a0 = aj0[v]
            s = rmul(g1, a0) + rmul(g2, a0)
            out.append(rmul(s, f_phi[v]))
        if want_apar:
            a0 = aj0[v]
            vp = vpa[v]
            s = rmul(g1, a0 * vp[..., 0, :]) + rmul(g2, a0 * vp[..., 1, :])
            out.append(rmul(s, f_apar[v]))
        if want_bpar:
            wb = aj1[v] * vperp2[v]
            s = rmul(g1, wb) + rmul(g2, wb)
            out.append(rmul(s, f_bpar[v]))
        return tuple(out)

    nwant = sum(bool(x) for x in want)
    if nwant == 0:
        return ()
    order = jnp.asarray(data.vel_order)
    if params.exact:
        def body(acc, v):
            return tuple(a + t for a, t in zip(acc, terms(v))), None

        init = tuple(_zeros_field((Y, X, T)) for _ in range(nwant))
        acc, _ = lax.scan(body, init, order, unroll=max(1, int(params.unroll)))
    else:
        acc = tuple(jnp.sum(t, axis=0) for t in jax.vmap(terms)(order))
    filt = data.kwork_filter if kwork_filter is None else kwork_filter
    filt = jnp.asarray(filt)[..., None]
    return tuple(jnp.where(filt, jnp.zeros((), jnp.complex128), a) for a in acc)


def getan(gnew, data: FieldData, params: FieldParams):
    """``getan_from_dfn(g_in, antot, antota, antotp)`` (src/dist_fn.fpp:6273-6360); ``getan``
    (6372-6380) is this with ``g_in = gnew``.

    Returns ``(antot, antota, antotp)``, F3 complex; the integral of an absent field is
    returned as zeros (GS2 leaves the stale array untouched; it is never read).
    """
    want = (params.has_phi, params.has_apar, params.has_bpar)
    res = list(_velocity_integrals(gnew, data, params, want))
    shape = gnew.shape[3:5] + gnew.shape[-1:]
    out = []
    for on in want:
        if on:
            arr = res.pop(0)
            if params.esv:
                arr = _ensure_single_val(arr, data)
            out.append(arr)
        else:
            out.append(_zeros_field(shape))
    return tuple(out)


def klein_sum(values, axis: int = -1):
    """``klein_sum`` (externals/utils/summation.f90:71-103): second-order compensated
    sequential sum over ``axis`` (real arrays)."""
    a = jnp.moveaxis(values, axis, 0)
    zero = jnp.zeros(a.shape[1:], a.dtype)

    def body(carry, x):
        s, c1, c2 = carry
        t = s + x
        lc1 = jnp.where(jnp.abs(s) >= jnp.abs(x), (s - t) + x, (x - t) + s)
        s = t
        t = c1 + lc1
        lc2 = jnp.where(jnp.abs(c1) >= jnp.abs(lc1), (c1 - t) + lc1, (lc1 - t) + c1)
        return (s, t, c2 + lc2), None

    (s, c1, c2), _ = lax.scan(body, (zero, zero, zero), a)
    return (s + c1) + c2


def trapz_complex(dtheta, y):
    """``complex_trapezoidal_integration(x, y)`` over the last axis (integration.f90:63-71):
    ``cmplx(klein(real(v)), klein(aimag(v))) * 0.5`` with ``v(i) = (x(i+1)-x(i))*(y(i)+y(i+1))``."""
    s = y[..., :-1] + y[..., 1:]
    v = rmul(s, dtheta)
    re = klein_sum(v.real, axis=-1)
    im = klein_sum(v.imag, axis=-1)
    return lax.complex(re * 0.5, im * 0.5)


def flux_surface_average(antot, data: FieldData, params: FieldParams, kwork_filter=None):
    """``calculate_flux_surface_average(fl_avg, antot)`` (src/dist_fn.fpp:8022-8034) with the
    frozen ``awgt`` (8006-8020): ``fl_avg = (tite*trapz(theta, (jacob*antot)/gamtot))*awgt``,
    0 for filtered modes. Returns K2 complex."""
    jac = jnp.asarray(data.jacob)
    gam = jnp.asarray(data.gamtot)
    y = lax.complex((jac * antot.real) / gam, (jac * antot.imag) / gam)
    tint = trapz_complex(jnp.asarray(data.dtheta), y)
    fl = rmul(rmul(tint, params.tite), jnp.asarray(data.awgt))
    filt = data.kwork_filter if kwork_filter is None else kwork_filter
    return jnp.where(jnp.asarray(filt), jnp.zeros((), jnp.complex128), fl)


def _fl_avg(antot, data, params, kwork_filter=None):
    if params.fl_avg_on:
        return flux_surface_average(antot, data, params, kwork_filter)
    return jnp.zeros(antot.shape[:2], jnp.complex128)


def getfieldeq(fields_new: Fields, antot, antota, antotp, data: FieldData,
               params: FieldParams, kwork_filter=None):
    """``getfieldeq1`` (src/dist_fn.fpp:6575-6665): field-equation residuals.

    ``fieldeq  = ((antot + bpar*gamtot1) - (gamtot*gridfac1)*phi) + fl_avg``      (6620-6623)
    ``fieldeqa = antota - (kperp2*gridfac1)*apar``                                (6638-6639)
    ``fieldeqp = ((antotp + bpar*gamtot2) + (0.5*phi)*gamtot1)*((beta*apfac)/bmag**2)
    + bpar*gridfac1``                                                             (6656-6659)

    ``fl_avg`` only with no electron-like species and adiabatic_option='iphi00=2'
    (6602-6609). Absent fields' residuals are returned as zeros.
    """
    phi, apar, bpar = fields_new
    shape = antot.shape
    fl = _fl_avg(antot, data, params, kwork_filter)
    fe = fea = fep = _zeros_field(shape)
    if params.has_phi:
        fe = ((antot + rmul(bpar, data.gamtot1)) - rmul(phi, data.gamtot_gf)) + fl[..., None]
    if params.has_apar:
        fea = antota - rmul(apar, data.kperp2_gf)
    if params.has_bpar:
        fep = rmul((antotp + rmul(bpar, data.gamtot2)) + rmul(rmul(phi, 0.5), data.gamtot1),
                   data.bfac) + rmul(bpar, data.gridfac1)
    return fe, fea, fep


getfieldeq1 = getfieldeq


def getfieldeq_from_g(fields_new: Fields, gnew, data: FieldData, params: FieldParams,
                      kwork_filter=None):
    """Fortran ``getfieldeq`` (dist_fn.fpp:6667-6680) = ``getan`` + ``getfieldeq1``."""
    ant = getan(gnew, data, params) if kwork_filter is None else \
        _getan_filtered(gnew, data, params, kwork_filter)
    return getfieldeq(fields_new, *ant, data, params, kwork_filter=kwork_filter)


def _getan_filtered(g, data, params, kwork_filter):
    want = (params.has_phi, params.has_apar, params.has_bpar)
    res = list(_velocity_integrals(g, data, params, want, kwork_filter=kwork_filter))
    shape = g.shape[3:5] + g.shape[-1:]
    out = []
    for on in want:
        if on:
            arr = res.pop(0)
            out.append(_ensure_single_val(arr, data) if params.esv else arr)
        else:
            out.append(_zeros_field(shape))
    return tuple(out)


def stack_present(fe3, params: FieldParams):
    """``fl`` of ``get_field_vector`` (fields_implicit.f90:139-171): present fields stacked
    as ``(nf, Y, X, T)`` in phi, apar, bpar order."""
    return jnp.stack([fe3[i] for i in params.slots], axis=0)


def gather_cells(fl, topo: Topology, c: int):
    """``fl(:, it_c(m,n), ik_c(m))`` for class ``c``: ``(M, N, nf*T)`` in B order within a
    cell (field outer, theta inner)."""
    ikc = jnp.asarray(topo.ik[c])
    itc = jnp.asarray(topo.it[c])
    M, N = topo.sizes(c)
    nf, T = fl.shape[0], fl.shape[-1]
    flc = fl[:, ikc[:, None], itc, :]                    # (nf, M, N, T)
    return jnp.transpose(flc, (1, 2, 0, 3)).reshape(M, N, nf * T)


def scatter_cells(u_list, topo: Topology, nf: int, nth: int):
    """Inverse of :func:`gather_cells` for the outputs: ``u_list[c]`` is ``(M_c, R_c)`` in B
    order; returns ``(nf, Y, X, T)`` (``get_field_solution``, fields_implicit.f90:174-219)."""
    parts = [u.reshape(u.shape[0] * (u.shape[1] // (nf * nth)), nf, nth) for u in u_list]
    U = jnp.concatenate(parts, axis=0)                   # (Ncell, nf, T)
    uf = U[jnp.asarray(topo.perm)]                       # (Y, X, nf, T)
    return jnp.moveaxis(uf, 2, 0)


#: Longest product sum that is unrolled in the HLO; longer ones use a fori_loop.
UNROLL_MAX = 64


def _seq_prod_sum(a, x):
    """``S = (((0 + a_1*x_1) + a_2*x_2) + ...)`` over the last axis (complex products
    ``numerics.cmul``, one rounded add per term, GS2's ``sum(a*b)`` and ``matmul`` order,
    spec 04 section 0.2.2). ``a`` and ``x`` broadcast against each other."""
    K = a.shape[-1]
    shape = jnp.broadcast_shapes(a.shape, x.shape)[:-1]
    S = jnp.zeros(shape, jnp.complex128)
    if K <= UNROLL_MAX:
        for k in range(K):
            S = S + cmul(a[..., k], x[..., k])
        return S

    def body(k, S):
        ak = lax.dynamic_index_in_dim(a, k, axis=a.ndim - 1, keepdims=False)
        xk = lax.dynamic_index_in_dim(x, k, axis=x.ndim - 1, keepdims=False)
        return S + cmul(ak, xk)

    return lax.fori_loop(0, K, body, S)


def apply_aminv(fl, mats, topo: Topology, params: FieldParams):
    """``u = -aminv.fl`` of ``getfield`` (fields_implicit.f90:261-289).

    For every output row ``r`` of member ``m``: ``u_r = ((0 - S_1) - S_2) - ... - S_N`` with
    ``S_n = (((0 + p_1) + p_2) + ...) + p_nidx`` and ``p_k = aminv(r, n, k)*fl_k(it_n, ik)``
    (k in B order). Returns ``(nf, Y, X, T)``.
    """
    nf, nth = params.nfield, params.nth
    nidx = nf * nth
    outs = []
    for c, mat in enumerate(mats):
        M, N = topo.sizes(c)
        R = N * nidx
        flc = gather_cells(fl, topo, c)                    # (M, N, nidx)
        mat4 = jnp.reshape(mat, (M, R, N, nidx))
        if params.exact:
            S = _seq_prod_sum(mat4, flc[:, None, :, :])          # (M, R, N)
            u = jnp.zeros((M, R), jnp.complex128)
            for n in range(N):
                u = u - S[:, :, n]
        else:
            u = -jnp.einsum("mrnk,mnk->mr", mat4, flc)
        outs.append(u)
    return scatter_cells(outs, topo, nf, nth)


def _fields_from_stack(u, fallback: Fields, params: FieldParams) -> Fields:
    out = list(fallback)
    for j, i in enumerate(params.slots):
        out[i] = u[j]
    return Fields(*out)


def getfield(fields_new: Fields, gnew, aminv, data: FieldData, params: FieldParams) -> Fields:
    """``getfield(phinew, aparnew, bparnew)`` (src/fields_implicit.f90:222-312).

    ``fields_new`` are the fields ``getfieldeq`` sees (in a time step ``phinew == phi`` of the
    first timeadv, fields_implicit.f90:338-352); ``gnew`` the predictor ``g``. Returns the
    increment ``u`` as Fields: the present fields hold ``u`` (``get_field_solution``), absent
    fields are returned unchanged from ``fields_new`` (GS2 does not write them). With
    ``params.field_option == 'local'`` this dispatches to
    :func:`gs2gpu.fields_local.getfield_local`.
    """
    if params.field_option == "local":
        from .fields_local import getfield_local
        return getfield_local(fields_new, gnew, aminv, data, params)
    fe = getfieldeq_from_g(fields_new, gnew, data, params)
    fl = stack_present(fe, params)
    mats = aminv.mats if hasattr(aminv, "mats") else aminv
    u = apply_aminv(fl, mats, data.topo, params)
    return _fields_from_stack(u, fields_new, params)


def field_update(u: Fields, fields_old: Fields, params: Optional[FieldParams] = None) -> Fields:
    """``phinew = phinew + phi ; aparnew = aparnew + apar ; bparnew = bparnew + bpar``
    (src/fields_implicit.f90:354-356), all three fields unconditionally (absent fields are 0).
    With ``params.field_option == 'local'`` only present fields are updated
    (fields_local.f90:4316-4320)."""
    if params is not None and params.field_option == "local":
        from .fields_local import field_update_local
        return field_update_local(u, fields_old, params)
    return Fields(u.phi + fields_old.phi, u.apar + fields_old.apar, u.bpar + fields_old.bpar)


def get_fields_direct_from_dfn(g_in, data: FieldData, params: FieldParams) -> Fields:
    """``get_fields_direct_from_dfn`` (src/dist_fn.fpp:7603-7790).

    ``phi  = [(beta*gamtot2 + bmag**2)*(antot + fl_avg) - (beta*gamtot1)*antotp] * inv_phi_denom``
    (``has_bpar``; else ``(antot + fl_avg)*inv_phi_denom``),
    ``apar = antota*inv_kperp2``,
    ``bpar = [-(beta*gamtot)*antotp - ((beta/2)*gamtot1)*antot] * inv_bpar_denom`` (``has_phi``;
    else ``(-beta*antotp)*inv_bpar_denom``). Absent fields are 0.
    """
    antot, antota, antotp = getan(g_in, data, params)
    shape = antot.shape
    zero = _zeros_field(shape)
    phi = apar = bpar = zero
    # GS2 quirk: gridfac1 does not enter the direct solve (consistent only for gridfac = 1;
    # dist_fn.fpp:7655-7772, spec 04 P27)
    if params.has_phi:
        a = antot + _fl_avg(antot, data, params)[..., None]
        if params.has_bpar:
            num = rmul(a, data.phi_num_a) - rmul(antotp, data.phi_num_b)
        else:
            num = a
        phi = rmul(num, data.inv_phi_denom)
    if params.has_apar:
        apar = rmul(antota, data.inv_kperp2)
    if params.has_bpar:
        if params.has_phi:
            num = (-rmul(antotp, data.bpar_num_a)) - rmul(antot, data.bpar_num_b)
        else:
            num = -rmul(antotp, params.beta)
        bpar = rmul(num, data.inv_bpar_denom)
    return Fields(phi, apar, bpar)


def get_init_field(gnew, data: FieldData, params: FieldParams) -> Fields:
    """``get_init_field`` (src/dist_fn.fpp:7792-7799): fields consistent with ``gnew``
    (t = 0, and after every dt reset with force_maxwell_reinit, gs2_init.fpp:347)."""
    return get_fields_direct_from_dfn(gnew, data, params)


def calculate_potentials_from_nonadiabatic_dfn(h, data: FieldData,
                                               params: FieldParams) -> Fields:
    """``calculate_potentials_from_nonadiabatic_dfn`` (src/dist_fn.fpp:7335-7455), used by
    ``collisions_advance`` with split_collisions (4130-4167, spec 10 G-10).

    ``phi = antot*phi_weight`` with ``phi_weight = 1/(sum((dens*z*z)/temp) [+ tite])``,
    ``apar = antota*inv_kperp2``, ``bpar = antotp*(-beta/bmag**2)``; no ``fl_avg``.
    """
    antot, antota, antotp = getan(h, data, params)
    zero = _zeros_field(antot.shape)
    # GS2 quirk: phi_weight ignores poisfac, yavg and fl_avg (dist_fn.fpp:7360-7366)
    phi = rmul(antot, params.phi_weight) if params.has_phi else zero
    apar = rmul(antota, data.inv_kperp2) if params.has_apar else zero
    bpar = rmul(antotp, data.bpar_weight) if params.has_bpar else zero
    return Fields(phi, apar, bpar)


# --------------------------------------------------------------------------- validation
def init_fieldeq_check(ex: Export, data: FieldData, params: FieldParams):
    """Recompute ``gamtot, gamtot1, gamtot2`` (and ``gamtot3``) as ``init_fieldeq`` does
    (src/dist_fn.fpp:6470-6572), for validating the integrator against the exported arrays.

    ``init_fieldeq`` runs at level 28 with an all-false kwork_filter (spec 04 P6).
    Returns a dict of NumPy arrays.
    """
    S, E, L, Y, X, T = ex.nspec, ex.negrid, ex.nlambda, ex.naky, ex.ntheta0, ex.nth
    spec = species_group(ex)
    f64 = lambda a: np.asarray(a, np.float64)  # noqa: E731
    z, dens, temp = (f64(spec[n]).reshape(-1) for n in ("spec_z", "spec_dens", "spec_temp"))
    nmc = f64(spec["nonmaxw_corr"]).reshape(S, E)[:, :, None, None, None, None]
    tracer = np.asarray(spec["spec_type"]).reshape(-1) == TRACER_SPECIES
    w = f64(ex["/le_grids/w"]).reshape(S, E)
    wl = f64(ex["/le_grids/wl"]).reshape(L, T)
    aj0, aj1, vperp2 = (jnp.asarray(f64(ex["/dist_fn_arrays"][n])) for n in ("aj0", "aj1", "vperp2"))
    kperp2 = f64(ex["/kt_grids/kperp2"])
    poisfac = _num(_first(ex, "/dist_fn/poisfac"), 0.0)
    aky = f64(ex["/kt_grids/aky"]).reshape(-1)

    def integ(gw_real, wgt):
        # integrate_species_original with both isgn copies equal (g(:,1)+g(:,2) = x+x) and an
        # all-false filter; only the real part of the complex Fortran sum is non-zero.
        wgt = np.where(tracer, 0.0, wgt)
        fac2 = ((wgt[:, None] * w)[:, :, None, None] * wl[None, None]).reshape(S * E * L, T)
        gw = jnp.asarray(gw_real).reshape(S * E * L, Y, X, T)
        order = jnp.asarray(data.vel_order)
        f2 = jnp.asarray(fac2)

        def body(acc, v):
            x = gw[v]
            return acc + ob(f2[v] * (x + x)), None

        tot, _ = lax.scan(body, jnp.zeros((Y, X, T)), order)
        return np.array(tot)

    g0 = (1.0 - ob(aj0 * aj0)) * nmc
    gamtot = integ(g0, ((z * z) * dens) / temp) + kperp2 * poisfac
    g1 = (((aj0 * aj1) * 2.0) * vperp2) * nmc
    gamtot1 = integ(g1, z * dens)
    g2 = (((aj1 * aj1) * 2.0) * (vperp2 * vperp2)) * nmc
    gamtot2 = integ(g2, temp * dens)
    if params.esv:
        for arr in (gamtot, gamtot1, gamtot2):
            arr[data.esv_ik, data.esv_dst, 0] = arr[data.esv_ik, data.esv_src, -1]
    gamtot3 = None
    if (not params.has_electron_species) or (not params.has_ion_species):
        if params.adiabatic_option == ADIABATIC_YAVG:
            for ik in range(Y):
                if aky[ik] > EPS:
                    gamtot[ik] = gamtot[ik] + params.tite
        elif params.adiabatic_option == ADIABATIC_FIELDLINEAVG:
            gamtot = gamtot + params.tite
            gamtot3 = gamtot3_from_gamtot(gamtot, params.tite)
        else:
            gamtot = gamtot + params.tite
    return {"gamtot": gamtot, "gamtot1": gamtot1, "gamtot2": gamtot2, "gamtot3": gamtot3}
