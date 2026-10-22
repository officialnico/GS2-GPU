"""Diagnostics of the new GS2 diagnostics module and the ``<run>.out.nc`` records (spec 08).

What this module reproduces (src/diagnostics/gs2_diagnostics_new.f90:339-519,
``run_diagnostics``) for the SUPPORTed options of spec 10 section B.7:

* every step (``istep > 0``): the instantaneous complex frequency, its ring buffer
  ``omegahist(0:navg-1)``, ``omega_average`` and the convergence / instability tests that can
  end the run (diagnostics_omega.f90:64-136);
* on write steps (``mod(istep, nwrite) == 0``, ``exit`` or the final forced call, never twice
  for the same istep): ``t``, omega, the field family (phi2 ... phi0, final fields), the QL
  metric, the fluxes (incl. the linear diffusive estimates, ``*_flux_vs_e`` and totals), the
  moments, ``phi_corr_2pi`` and the eigenfunctions;
* the ``g_adjust`` round trips of ``calculate_fluxes`` (gnew and g) and ``write_moments``
  (gnew) that make ``gnew <- (gnew + adj) - adj`` on every write step (spec 08 section 3.1.8,
  pitfall P1). These change the evolved state at round-off level and run on the device;
* ``finalize``: the forced final ``run_diagnostics`` call and ``run_old_final_routines``
  (final fields, epar, final moments with one more round trip, antot).

Split of the work. Everything that touches the distribution function (round trips,
``integrate_moment`` over velocity space, the per-``iglo`` theta integrals of
``get_flux_vs_e``) runs on the JAX device in :func:`diagnostics_device`. Its outputs are
small (``(S, Y, X, T)`` moments, ``(S, E, L, Y, X, 2)`` integrals); the theta integrals,
volume averages, normalisations and the scalar omega logic then run on the host in NumPy
with the Fortran operation order, one IEEE operation per NumPy call (NumPy's own complex
multiply/divide/abs are avoided: they contract to FMA or use other algorithms on this
platform; ``np.log`` of complex is the libm ``clog`` and ``np.hypot`` the libm ``hypot``
that gfortran calls for ``log`` and ``abs``). With ``numerics.set_strict(True)`` the device
part is FMA-free too, so the whole chain follows a ``-ffp-contract=off`` GS2 build.

Layouts are the canonical ones of spec 11 section 3: G7 ``(S, E, L, Y, X, 2, T)``, fields
F3 ``(Y, X, T)``, per-mode results ``(S, Y, X)`` (Fortran ``(it, ik, is)``), moments
``(S, Y, X, T)`` (Fortran ``(ig, it, ik, is)``).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, NamedTuple, Optional

import jax
import jax.numpy as jnp
import numpy as np
from jax import lax

from . import outnc
from .gadjust import FROM_G, TO_G, GAdjustData, g_adjust_direction
from .numerics import cdiv as _cdiv_jax
from .numerics import imul, ob, rmul
from .types import Fields

#: Fortran ``tiny(0.0)`` and ``epsilon(0.0)`` with ``-fdefault-real-8``.
TINY = float(np.finfo(np.float64).tiny)
EPS = float(np.finfo(np.float64).eps)

_MISSING = object()

# ============================================================================ export helpers


def _lookup(grp, names, default=_MISSING):
    """First of ``names`` present in export group ``grp`` (None = no group)."""
    if isinstance(names, str):
        names = (names,)
    if grp is not None:
        for name in names:
            if name in grp:
                return grp[name]
    if default is not _MISSING:
        return default
    raise KeyError(f"none of {names} in {getattr(grp, 'path', grp)}")


def _group(ex, path):
    try:
        return ex[path]
    except KeyError:
        return None


def _as_bool(v) -> bool:
    return bool(np.asarray(v).item())


def _as_int(v) -> int:
    return int(np.asarray(v).item())


def _as_float(v) -> float:
    return float(np.asarray(v).item())


# ============================================================================ gnostics flags

#: Defaults of ``diagnostics_base_config_type`` (src/diagnostics_base_configuration.f90:13-536)
#: after ``apply_new_defaults`` (src/diagnostics/diagnostics_config.f90:252-265).
GNOSTICS_DEFAULTS: dict[str, Any] = dict(
    write_any=True, nwrite=10, navg=10, omegatol=1e-3, omegatinst=1.0e6,
    exit_when_converged=True, igomega=0, nwrite_mult=10, ncheck=10, nc_sync_freq=10,
    print_line=False, print_flux_line=False, write_line=True, write_flux_line=False,
    write_fields=True, write_phi_over_time=False, write_apar_over_time=False,
    write_bpar_over_time=False, write_omega=True, write_ql_metric=True,
    write_fluxes=False, write_fluxes_by_mode=False, write_moments=True,
    write_ntot_over_time=False, write_density_over_time=False, write_upar_over_time=False,
    write_tperp_over_time=False, write_correlation=True, write_eigenfunc=False,
    write_final_fields=False, write_final_moments=False, write_final_epar=False,
    write_final_antot=False, write_ascii=True,
    # GATED (spec 10 B.7); kept so that a non-default value can be refused
    save_distfn=False, append_old=False, make_movie=False, dump_fields_periodically=False,
    write_heating=False, write_collisional=False, write_verr=False, write_max_verr=False,
    write_cerr=False, write_symmetry=False, write_nl_flux_dist=False, write_parity=False,
    write_cross_phase=False, write_correlation_extend=False, write_jext=False,
    write_lorentzian=False, write_zonal_transfer=False, write_kinetic_energy_transfer=False,
    write_full_moments_notgc=False, write_g=False, write_gyx=False, write_gs=False,
    write_kpar=False, write_final_db=False, use_nonlin_convergence=False,
)

#: Options GATED by spec 10 B.7 / 08 section 5 (refused when true).
GATED_FLAGS = ("save_distfn", "append_old", "make_movie", "dump_fields_periodically",
               "write_heating", "write_collisional", "write_verr", "write_max_verr", "write_cerr",
               "write_symmetry", "write_nl_flux_dist", "write_parity", "write_cross_phase",
               "write_correlation_extend", "write_jext", "write_lorentzian",
               "write_zonal_transfer", "write_kinetic_energy_transfer",
               "write_full_moments_notgc", "write_g", "write_gyx", "write_gs", "write_kpar",
               "write_final_db")


def _gnostics_from_namelist(text: Optional[str], nonlin: bool, use_le_layout: bool = True) -> dict:
    """``read_parameters`` + ``check_parameters`` from the input file text (fallback only).

    Mirrors diagnostics_config.f90:267-395 (smart defaults: ``write_fluxes`` and
    ``write_fluxes_by_mode`` are true for nonlinear runs; ``use_le_layout`` forces
    ``write_cerr = F``) and gs2_diagnostics_new.f90:188-208 (``print_line .or. write_line``
    forces ``write_fields`` and ``write_omega``; the flux lines force ``write_fields``).
    """
    flags = dict(GNOSTICS_DEFAULTS)
    if nonlin:
        flags["write_fluxes"] = True
        flags["write_fluxes_by_mode"] = True
    if text:
        import f90nml
        nml = f90nml.reads(text)
        knobs = nml.get("gs2_diagnostics_knobs", {})
        for key, value in knobs.items():
            if key.lower() in flags:
                flags[key.lower()] = value
    if use_le_layout:
        flags["write_cerr"] = False
    return _check_parameters(flags)


def _check_parameters(flags: dict) -> dict:
    """``check_parameters`` (gs2_diagnostics_new.f90:188-208)."""
    flags = dict(flags)
    if (flags["print_line"] or flags["write_line"]) and not (flags["write_fields"] and flags["write_omega"]):
        flags["write_fields"] = True
        flags["write_omega"] = True
    if (flags["print_flux_line"] or flags["write_flux_line"]) and not flags["write_fields"]:
        flags["write_fields"] = True
    return flags


# ============================================================================ layout orders


def layout_orders(layout: str, iglo_of=None) -> tuple[bool, tuple[str, str, str]]:
    """Summation orders fixed by the ``g_lo`` layout (gs2_layouts.fpp:917-930).

    Returns ``(le_e_outer, fve_order)``:

    * ``le_e_outer``: in ``integrate_moment`` (le_grids.f90, ``integrate_moment_c34``) the
      contributions to a fixed ``(it, ik, is)`` arrive in increasing ``iglo``; True means
      ``do ie: do il`` (il varies faster), False means ``do il: do ie``;
    * ``fve_order``: in ``get_flux_vs_e`` (diagnostics_fluxes.fpp:865-882) the order of
      ``(il, it, ik)`` for fixed ``(ie, is)``, outermost first, as names from
      ``{"il", "it", "ik"}``.

    The layout string lists the indices fastest first (``x``=it, ``y``=ik, ``l``=il, ``e``=ie,
    ``s``=is). If ``iglo_of`` (the exported G5 map ``(S, E, L, Y, X)`` of 0-based iglo) is
    given, the orders are taken from it and checked against the string.
    """
    layout = layout.strip().lower()
    if sorted(layout) != sorted("xyles"):
        raise ValueError(f"unknown layout {layout!r}")
    pos = {c: layout.index(c) for c in layout}
    le_e_outer = pos["l"] < pos["e"]
    names = {"l": "il", "x": "it", "y": "ik"}
    fve = tuple(names[c] for c in sorted("lxy", key=lambda c: -pos[c]))
    if iglo_of is not None:
        m = np.asarray(iglo_of)
        if m.ndim == 5:
            # strides of each index in iglo (S, E, L, Y, X axes)
            base = m[0, 0, 0, 0, 0]
            stride = {}
            for ax, key in ((1, "e"), (2, "l"), (3, "y"), (4, "x")):
                if m.shape[ax] > 1:
                    idx = [0] * 5
                    idx[ax] = 1
                    stride[key] = int(m[tuple(idx)] - base)
            if "e" in stride and "l" in stride:
                if (stride["l"] < stride["e"]) != le_e_outer:
                    raise ValueError(f"iglo_of disagrees with layout {layout!r} (e/l order)")
            have = [c for c in "lxy" if c in stride]
            order = tuple(names[c] for c in sorted(have, key=lambda c: -stride[c]))
            if order != tuple(o for o in fve if o in order):
                raise ValueError(f"iglo_of disagrees with layout {layout!r} (l/x/y order)")
    return le_e_outer, fve


# ============================================================================ parameters


class DiagParams(NamedTuple):
    """Static switches and sizes of the diagnostics (hashable; ``jit`` static argument).

    The ``gnostics`` flags are the values after ``read_parameters`` and ``check_parameters``
    as exported in ``/diagnostics`` (spec 08 section 2.3, 10 A.13). ``igomega`` keeps its
    Fortran value (``-ntgrid..ntgrid``).
    """

    ntgrid: int
    nth: int
    ntheta0: int
    naky: int
    nlambda: int
    negrid: int
    nspec: int
    nonlin: bool
    has_phi: bool
    has_apar: bool
    has_bpar: bool
    fphi: float
    fapar: float
    fbpar: float
    def_parity: bool
    even: bool
    layout: str
    le_e_outer: bool
    fve_order: tuple
    any_ieqzip: bool
    # gnostics
    write_any: bool
    nwrite: int
    navg: int
    omegatol: float
    omegatinst: float
    exit_when_converged: bool
    igomega: int
    nc_sync_freq: int
    ncheck: int
    write_omega: bool
    write_fields: bool
    write_phi_over_time: bool
    write_apar_over_time: bool
    write_bpar_over_time: bool
    write_ql_metric: bool
    write_fluxes: bool
    write_fluxes_by_mode: bool
    print_flux_line: bool
    write_flux_line: bool
    write_moments: bool
    write_ntot_over_time: bool
    write_density_over_time: bool
    write_upar_over_time: bool
    write_tperp_over_time: bool
    write_correlation: bool
    write_eigenfunc: bool
    write_final_fields: bool
    write_final_moments: bool
    write_final_epar: bool
    write_final_antot: bool
    print_line: bool
    write_line: bool
    write_ascii: bool

    @property
    def calculate_fluxes(self) -> bool:
        """gs2_diagnostics_new.f90:410-412."""
        return bool(self.write_fluxes or self.print_flux_line or self.write_flux_line)

    @property
    def omega_active(self) -> bool:
        """``calculate_omega`` returns immediately for ``nonlin .and. .not. any(ieqzip)``
        (diagnostics_omega.f90:80-82)."""
        return bool((not self.nonlin) or self.any_ieqzip)

    @property
    def ig(self) -> int:
        """NumPy theta index of ``igomega``."""
        return self.igomega + self.ntgrid

    @property
    def fields_present(self) -> tuple:
        return tuple(f for f, on in (("phi", self.has_phi), ("apar", self.has_apar),
                                     ("bpar", self.has_bpar)) if on)

    @classmethod
    def from_export(cls, ex) -> "DiagParams":
        """Sizes, field switches, layout orders and the parsed ``gnostics`` flags.

        Reads ``/diagnostics`` (exported post-``check_parameters``); if that group is
        missing, falls back to parsing ``gs2_diagnostics_knobs`` from the input text in the
        export's global attributes with GS2's defaults (``_gnostics_from_namelist``).
        """
        rp = _group(ex, "/run_parameters")
        nonlin = _as_bool(_lookup(_group(ex, "/nonlinear_terms"), "nonlin", False))
        fphi = _as_float(_lookup(rp, "fphi", 1.0))
        fapar = _as_float(_lookup(rp, "fapar", 0.0))
        fbpar = _as_float(_lookup(rp, "fbpar", 0.0))
        has_phi = _as_bool(_lookup(rp, "has_phi", fphi > EPS))
        has_apar = _as_bool(_lookup(rp, "has_apar", fapar > EPS))
        has_bpar = _as_bool(_lookup(rp, "has_bpar", fbpar > EPS))
        ieqzip = _lookup(rp, "ieqzip", None)
        any_ieqzip = bool(np.any(np.asarray(ieqzip))) if ieqzip is not None else False
        if _as_bool(_lookup(rp, "wstar_units", False)):
            raise NotImplementedError("knobs:wstar_units=T is gated (spec 10 G-16)")
        if _as_bool(_lookup(rp, "use_old_diagnostics", False)):
            raise NotImplementedError("knobs:use_old_diagnostics=T is gated (spec 10 G-09)")
        if any_ieqzip:
            raise NotImplementedError("knobs:eqzip_option is gated (spec 10 G-07)")
        df = _group(ex, "/dist_fn")
        def_parity = _as_bool(_lookup(df, "def_parity", False))
        even = _as_bool(_lookup(df, "even", True))

        layout = None
        lay = _group(ex, "/gs2_layouts")
        if lay is not None:
            layout = lay.attrs.get("layout") or _lookup(lay, "layout", None)
        if not layout:
            layout = ex.attrs.get("layout", "lxyes")
        if isinstance(layout, (bytes, np.ndarray)):
            layout = np.asarray(layout).tobytes().decode() if isinstance(layout, np.ndarray) \
                else layout.decode()
        iglo_of = _lookup(lay, "iglo_of", None)
        le_e_outer, fve_order = layout_orders(str(layout), iglo_of)

        dg = _group(ex, "/diagnostics")
        if dg is not None:
            flags = dict(GNOSTICS_DEFAULTS)
            if nonlin:
                flags["write_fluxes"] = True
                flags["write_fluxes_by_mode"] = True
            for key in list(flags):
                value = _lookup(dg, key, None)
                if value is None:
                    value = dg.attrs.get(key)
                if value is not None:
                    flags[key] = value
            flags = _check_parameters(flags)
        else:
            text = ex.attrs.get("input_file_text") or ex.attrs.get("used_inputs")
            coll = _group(ex, "/collisions")
            use_le = _as_bool(_lookup(coll, "use_le_layout", True))
            flags = _gnostics_from_namelist(text, nonlin, use_le)
        for key in GATED_FLAGS:
            if _as_bool(flags[key]):
                raise NotImplementedError(f"gs2_diagnostics_knobs:{key}=T is gated (spec 10 B.7)")
        if nonlin and _as_bool(flags["use_nonlin_convergence"]):
            raise NotImplementedError("use_nonlin_convergence=T is gated (spec 10 G-22)")
        if not _as_bool(flags["write_any"]):
            pass  # supported: run_diagnostics returns immediately
        ntgrid = ex.ntgrid
        igomega = _as_int(flags["igomega"])
        if not -ntgrid <= igomega <= ntgrid:
            raise ValueError(f"igomega={igomega} outside -ntgrid..ntgrid")
        b = _as_bool
        return cls(
            ntgrid=ntgrid, nth=ex.nth, ntheta0=ex.ntheta0, naky=ex.naky,
            nlambda=ex.nlambda, negrid=ex.negrid, nspec=ex.nspec,
            nonlin=nonlin, has_phi=has_phi, has_apar=has_apar, has_bpar=has_bpar,
            fphi=fphi, fapar=fapar, fbpar=fbpar, def_parity=def_parity, even=even,
            layout=str(layout).strip(), le_e_outer=le_e_outer, fve_order=fve_order,
            any_ieqzip=any_ieqzip,
            write_any=b(flags["write_any"]), nwrite=_as_int(flags["nwrite"]),
            navg=_as_int(flags["navg"]), omegatol=_as_float(flags["omegatol"]),
            omegatinst=_as_float(flags["omegatinst"]),
            exit_when_converged=b(flags["exit_when_converged"]), igomega=igomega,
            nc_sync_freq=_as_int(flags["nc_sync_freq"]), ncheck=_as_int(flags["ncheck"]),
            write_omega=b(flags["write_omega"]), write_fields=b(flags["write_fields"]),
            write_phi_over_time=b(flags["write_phi_over_time"]),
            write_apar_over_time=b(flags["write_apar_over_time"]),
            write_bpar_over_time=b(flags["write_bpar_over_time"]),
            write_ql_metric=b(flags["write_ql_metric"]), write_fluxes=b(flags["write_fluxes"]),
            write_fluxes_by_mode=b(flags["write_fluxes_by_mode"]),
            print_flux_line=b(flags["print_flux_line"]),
            write_flux_line=b(flags["write_flux_line"]), write_moments=b(flags["write_moments"]),
            write_ntot_over_time=b(flags["write_ntot_over_time"]),
            write_density_over_time=b(flags["write_density_over_time"]),
            write_upar_over_time=b(flags["write_upar_over_time"]),
            write_tperp_over_time=b(flags["write_tperp_over_time"]),
            write_correlation=b(flags["write_correlation"]),
            write_eigenfunc=b(flags["write_eigenfunc"]),
            write_final_fields=b(flags["write_final_fields"]),
            write_final_moments=b(flags["write_final_moments"]),
            write_final_epar=b(flags["write_final_epar"]),
            write_final_antot=b(flags["write_final_antot"]),
            print_line=b(flags["print_line"]), write_line=b(flags["write_line"]),
            write_ascii=b(flags["write_ascii"]),
        )


# ============================================================================ data


class DiagData(NamedTuple):
    """Arrays read by the diagnostics (NumPy from ``from_export``; ``to_jax`` before use).

    G7/G6 arrays are only used by the device kernel; the small arrays are also pulled to the
    host for the NumPy part. ``le_ie, le_il``: the ``(ie, il)`` sequence of
    ``integrate_moment`` (0-based), see :func:`layout_orders`. Derived geometry
    (``sq_bpol``, ``vperp_den``, ``gds_geo``, ``flux_wgt``) is evaluated once here, in
    NumPy, with the Fortran expression order.
    """

    gadj: GAdjustData
    vpa: Any            # G7 real
    vpac: Any           # G7 real
    w: Any              # (S, E)    le_grids w(ie,is)
    wl: Any             # (L, T)    le_grids wl(ig,il)
    energy: Any         # (S, E)
    nonmaxw_corr: Any   # (S, E)
    z: Any              # (S,) species
    dens: Any
    temp: Any
    mass: Any
    tprim: Any
    fprim: Any
    stm: Any
    tz: Any
    zt: Any
    smz: Any
    theta: Any          # (T,)
    jacob: Any
    grho: Any
    gradpar: Any
    delthet: Any
    rplot: Any
    sq_bpol: Any        # sqrt(1.0 - Bpol**2/bmag**2)       (diagnostics_fluxes.fpp:602)
    vperp_den: Any      # qval*shat*bmag**2                  (diagnostics_fluxes.fpp:614)
    gds_geo: Any        # (Y, X, T) gds21 + theta0*gds22     (diagnostics_fluxes.fpp:614)
    aky: Any            # (Y,)
    akx: Any            # (X,)
    theta0: Any         # (Y, X)
    kperp2: Any         # (Y, X, T)
    woutunits: Any      # (Y,)
    tunits: Any         # (Y,)
    le_ie: Any          # (E*L,) int32
    le_il: Any          # (E*L,) int32
    rhoc: float
    field_line_average_weight: float
    flux_wgt: float

    @classmethod
    def from_export(cls, ex) -> "DiagData":
        tg = ex["/theta_grid"]
        kt = ex["/kt_grids"]
        le = ex["/le_grids"]
        sp = ex["/species"]
        da = ex["/dist_fn_arrays"]
        gt = _group(ex, "/gs2_time")
        S, E, L = ex.nspec, ex.negrid, ex.nlambda
        Y, X, T = ex.naky, ex.ntheta0, ex.nth

        theta = np.asarray(tg["theta"], np.float64)
        jacob = np.asarray(tg["jacob"], np.float64)
        grho = np.asarray(tg["grho"], np.float64)
        bmag = np.asarray(tg["bmag"], np.float64)
        bpol = np.asarray(_lookup(tg, ("Bpol", "bpol")), np.float64)
        rplot = np.asarray(_lookup(tg, ("Rplot", "rplot")), np.float64)
        gds21 = np.asarray(tg["gds21"], np.float64)
        gds22 = np.asarray(tg["gds22"], np.float64)
        qval = _as_float(tg["qval"])
        shat = _as_float(tg["shat"])
        rhoc = _as_float(_lookup(tg, "rhoc"))
        theta0 = np.asarray(kt["theta0"], np.float64)

        # diagnostics_fluxes.fpp:602  sqrt(1.0-Bpol(ig)**2/bmag(ig)**2)
        sq_bpol = np.sqrt(1.0 - (bpol * bpol) / (bmag * bmag))
        # diagnostics_fluxes.fpp:614  (gds21+theta0(it,ik)*gds22) and (qval*shat*bmag**2)
        gds_geo = gds21[None, None, :] + theta0[:, :, None] * gds22[None, None, :]
        vperp_den = (qval * shat) * (bmag * bmag)
        # theta_grid.f90:520 (exported; recomputed here, the test asserts bitwise equality)
        # and diagnostics_fluxes.fpp:825 (recomputed: it is a local of get_flux)
        flw_re = float(1.0 / host_trapz(theta, jacob))
        flw = _as_float(_lookup(tg, "field_line_average_weight", flw_re))
        flux_wgt = float(1.0 / host_trapz(theta, grho * jacob))

        if gt is not None and "woutunits" in gt:
            woutunits = np.asarray(gt["woutunits"], np.float64)
            tunits = np.asarray(gt["tunits"], np.float64)
        else:   # gs2_time.f90:121-125 with wstar_units = F (gated otherwise)
            tunits = np.ones(Y)
            woutunits = 1.0 / tunits

        layout_str = None
        lay = _group(ex, "/gs2_layouts")
        if lay is not None:
            layout_str = lay.attrs.get("layout") or _lookup(lay, "layout", None)
        layout_str = layout_str or ex.attrs.get("layout", "lxyes")
        if isinstance(layout_str, np.ndarray):
            layout_str = layout_str.tobytes().decode()
        le_e_outer, _ = layout_orders(str(layout_str), _lookup(lay, "iglo_of", None))
        if le_e_outer:
            ie_seq, il_seq = np.meshgrid(np.arange(E), np.arange(L), indexing="ij")
        else:
            il_seq, ie_seq = np.meshgrid(np.arange(L), np.arange(E), indexing="ij")

        return cls(
            gadj=GAdjustData.from_export(ex),
            vpa=np.asarray(da["vpa"], np.float64), vpac=np.asarray(da["vpac"], np.float64),
            w=np.asarray(le["w"], np.float64), wl=np.asarray(le["wl"], np.float64),
            energy=np.asarray(le["energy"], np.float64),
            nonmaxw_corr=np.asarray(sp["nonmaxw_corr"], np.float64),
            z=np.asarray(sp["spec_z"], np.float64), dens=np.asarray(sp["spec_dens"], np.float64),
            temp=np.asarray(sp["spec_temp"], np.float64),
            mass=np.asarray(sp["spec_mass"], np.float64),
            tprim=np.asarray(sp["spec_tprim"], np.float64),
            fprim=np.asarray(sp["spec_fprim"], np.float64),
            stm=np.asarray(sp["spec_stm"], np.float64), tz=np.asarray(sp["spec_tz"], np.float64),
            zt=np.asarray(sp["spec_zt"], np.float64), smz=np.asarray(sp["spec_smz"], np.float64),
            theta=theta, jacob=jacob, grho=grho,
            gradpar=np.asarray(tg["gradpar"], np.float64),
            delthet=np.asarray(tg["delthet"], np.float64), rplot=rplot, sq_bpol=sq_bpol,
            vperp_den=vperp_den, gds_geo=gds_geo,
            aky=np.asarray(kt["aky"], np.float64), akx=np.asarray(kt["akx"], np.float64),
            theta0=theta0, kperp2=np.asarray(kt["kperp2"], np.float64),
            woutunits=woutunits, tunits=tunits,
            le_ie=ie_seq.reshape(-1).astype(np.int32), le_il=il_seq.reshape(-1).astype(np.int32),
            rhoc=rhoc, field_line_average_weight=flw, flux_wgt=flux_wgt,
        )


class HostData(NamedTuple):
    """NumPy copies of the small arrays of :class:`DiagData` used on the host."""

    w: np.ndarray
    energy: np.ndarray
    z: np.ndarray
    dens: np.ndarray
    temp: np.ndarray
    mass: np.ndarray
    tprim: np.ndarray
    fprim: np.ndarray
    stm: np.ndarray
    tz: np.ndarray
    zt: np.ndarray
    smz: np.ndarray
    theta: np.ndarray
    jacob: np.ndarray
    grho: np.ndarray
    gradpar: np.ndarray
    delthet: np.ndarray
    aky: np.ndarray
    akx: np.ndarray
    theta0: np.ndarray
    kperp2: np.ndarray
    woutunits: np.ndarray
    tunits: np.ndarray
    W: float            # field_line_average_weight
    flux_wgt: float

    @classmethod
    def from_data(cls, data: DiagData) -> "HostData":
        a = lambda x: np.asarray(x, dtype=np.float64)  # noqa: E731
        return cls(w=a(data.w), energy=a(data.energy), z=a(data.z), dens=a(data.dens),
                   temp=a(data.temp), mass=a(data.mass), tprim=a(data.tprim),
                   fprim=a(data.fprim), stm=a(data.stm), tz=a(data.tz), zt=a(data.zt),
                   smz=a(data.smz), theta=a(data.theta), jacob=a(data.jacob),
                   grho=a(data.grho), gradpar=a(data.gradpar), delthet=a(data.delthet),
                   aky=a(data.aky), akx=a(data.akx), theta0=a(data.theta0),
                   kperp2=a(data.kperp2), woutunits=a(data.woutunits), tunits=a(data.tunits),
                   W=float(np.asarray(data.field_line_average_weight)),
                   flux_wgt=float(np.asarray(data.flux_wgt)))


# ============================================================================ host numerics
#
# One IEEE operation per NumPy ufunc call, in Fortran association order. Complex values are
# complex128 arrays, but products/quotients/moduli are formed from their real parts.


def _cplx(re, im):
    re = np.asarray(re, np.float64)
    im = np.asarray(im, np.float64)
    out = np.empty(np.broadcast_shapes(re.shape, im.shape), np.complex128)
    out.real = re
    out.imag = im
    return out


def _cmul(a, b):
    """gfortran complex ``a*b``: ``(ar*br - ai*bi, ar*bi + ai*br)``, no FMA."""
    ar, ai, br, bi = a.real, a.imag, b.real, b.imag
    return _cplx(ar * br - ai * bi, ar * bi + ai * br)


def _cmul_conj(a, b):
    """``a*conjg(b)`` evaluated as the complex product with ``(br, -bi)``."""
    ar, ai, br, bi = a.real, a.imag, b.real, -b.imag
    return _cplx(ar * br - ai * bi, ar * bi + ai * br)


def _crmul(a, r):
    """complex times real (gfortran treats the real as having a zero imaginary part)."""
    return _cplx(a.real * r, a.imag * r)


def _crdiv(a, r):
    """complex / real: componentwise."""
    return _cplx(a.real / r, a.imag / r)


def _cdiv(a, b):
    """gfortran complex ``a/b`` (GCC's Smith method, ``expand_complex_div_wide``)."""
    ar, ai, br, bi = a.real, a.imag, b.real, b.imag
    t = np.abs(br) < np.abs(bi)
    with np.errstate(all="ignore"):
        r1 = br / bi
        d1 = br * r1 + bi
        x1 = ar * r1 + ai
        y1 = ai * r1 - ar
        r2 = bi / br
        d2 = bi * r2 + br
        x2 = ai * r2 + ar
        y2 = ai - ar * r2
        re = np.where(t, x1 / d1, x2 / d2)
        im = np.where(t, y1 / d1, y2 / d2)
    return _cplx(re, im)


def _cabs(a):
    """gfortran ``abs`` of complex = libm ``cabs`` = ``hypot``."""
    return np.hypot(a.real, a.imag)


def _clog(a):
    """gfortran ``log`` of complex = libm ``clog`` (``np.log`` calls it for complex128)."""
    return np.log(np.asarray(a, np.complex128))


def _seq_sum(x, axis):
    """Fortran ``SUM`` / an explicit ``s = s + x`` loop: sequential, left to right."""
    x = np.asarray(x)
    if x.shape[axis] == 0:
        return np.zeros(np.delete(x.shape, axis), x.dtype)
    if np.iscomplexobj(x):
        return _cplx(np.cumsum(x.real, axis=axis).take(-1, axis=axis),
                     np.cumsum(x.imag, axis=axis).take(-1, axis=axis))
    return np.cumsum(x, axis=axis).take(-1, axis=axis)


def host_klein_sum(a):
    """``klein_sum`` (externals/utils/summation.f90:71-103) along axis 0, vectorised."""
    a = np.asarray(a, np.float64)
    res = np.zeros(a.shape[1:])
    c1 = np.zeros(a.shape[1:])
    c2 = np.zeros(a.shape[1:])
    for x in a:
        t = res + x
        lc1 = np.where(np.abs(res) >= np.abs(x), (res - t) + x, (x - t) + res)
        res = t
        t = c1 + lc1
        lc2 = np.where(np.abs(c1) >= np.abs(lc1), (c1 - t) + lc1, (lc1 - t) + c1)
        c1 = t
        c2 = c2 + lc2
    return (res + c1) + c2


def host_trapz(x, y):
    """Real ``trapezoidal_integration`` over the last axis (integration.f90:55-62)."""
    x = np.asarray(x, np.float64)
    y = np.asarray(y, np.float64)
    d = (x[1:] - x[:-1]) * (y[..., :-1] + y[..., 1:])
    return host_klein_sum(np.moveaxis(d, -1, 0)) * 0.5


def host_trapz_c(x, z):
    """Complex ``trapezoidal_integration`` (integration.f90:64-71): Klein sums of the real
    and imaginary parts of ``dx*(z(:n-1)+z(2:))``, times 0.5."""
    x = np.asarray(x, np.float64)
    dx = x[1:] - x[:-1]
    s = z[..., :-1] + z[..., 1:]
    re = host_klein_sum(np.moveaxis(dx * s.real, -1, 0)) * 0.5
    im = host_klein_sum(np.moveaxis(dx * s.imag, -1, 0)) * 0.5
    return _cplx(re, im)


def host_flavg(q, h: HostData):
    """``field_line_average`` of a real quantity (theta_grid.f90:657-662), last axis theta."""
    return host_trapz(h.theta, q * h.jacob) * h.W


def host_flavg_c(q, h: HostData):
    """``field_line_average`` of a complex quantity (theta_grid.f90:665-670)."""
    return _crmul(host_trapz_c(h.theta, _crmul(q, h.jacob)), h.W)


def _abs2_conj(a):
    """``real(conjg(a)*a)`` = ``ar*ar - (-ai)*ai``."""
    return a.real * a.real - (-a.imag) * a.imag


def _ky_fac(aky):
    """``fac = 0.5; if (aky(ik) == 0.) fac = 1.0`` (volume_averages.f90:325-327)."""
    return np.where(np.asarray(aky) == 0.0, 1.0, 0.5)


def average_ky(f, aky):
    """``average_ky_xy`` / ``average_ky_xys``: ``f[..., Y, X] -> [..., X]`` (volume_averages.f90:373-441)."""
    fac = _ky_fac(aky)
    return _seq_sum(f * fac[:, None], axis=-2)


def average_kx(f, aky):
    """``average_kx_xy`` / ``average_kx_xys``: ``f[..., Y, X] -> [..., Y]`` (volume_averages.f90:444-489)."""
    fac = _ky_fac(aky)
    return _seq_sum(f * fac[:, None], axis=-1)


# ============================================================================ device numerics


def div_exact(x, d):
    """Correctly rounded ``x / d`` on the device.

    XLA's algebraic simplifier rewrites a division by a broadcast operand (a scalar, a
    constant or a lower-rank array) into a multiplication by the broadcast reciprocal,
    which is not correctly rounded (about a third of the results are 1 ulp off on XLA:CPU).
    Materialising the broadcast divisor behind an optimisation barrier keeps a true IEEE
    division; the cost is one extra divisor-sized array, only used on write steps here.
    """
    d = jnp.broadcast_to(jnp.asarray(d, x.dtype), jnp.shape(x))
    return x / lax.optimization_barrier(d)


def klein_sum_jax(a, axis: int = -1):
    """``klein_sum`` along ``axis`` with ``lax.scan`` (strictly sequential, as summation.f90:71-103)."""
    a = jnp.moveaxis(a, axis, 0)
    zero = jnp.zeros(a.shape[1:], a.dtype)

    def step(carry, x):
        res, c1, c2 = carry
        t = res + x
        lc1 = jnp.where(jnp.abs(res) >= jnp.abs(x), (res - t) + x, (x - t) + res)
        res = t
        t = c1 + lc1
        lc2 = jnp.where(jnp.abs(c1) >= jnp.abs(lc1), (c1 - t) + lc1, (lc1 - t) + c1)
        return (res, t, c2 + lc2), None

    (res, c1, c2), _ = lax.scan(step, (zero, zero, zero), a)
    return (res + c1) + c2


def trapz_jax(x, y):
    """Real ``trapezoidal_integration`` over the last axis on the device (integration.f90:55-62)."""
    dx = x[1:] - x[:-1]
    d = ob(dx * (y[..., :-1] + y[..., 1:]))
    return klein_sum_jax(d, -1) * 0.5


# ============================================================================ omega


def omega_quotient_np(num, den, code_dt: float, thresh: float):
    """``log(num/den)*zi/code_dt``, zero where ``|num|`` or ``|den| < thresh``.

    Mirrors the ``where`` blocks of diagnostics_omega.f90:99-106 and
    gs2_diagnostics.fpp:2952-2956: complex division (Smith), libm ``clog``, the product
    with ``zi = (0, 1)`` (``(-Im, Re)``, exact) and the division by the real ``code_dt``.
    """
    small = (_cabs(num) < thresh) | (_cabs(den) < thresh)
    with np.errstate(all="ignore"):
        q = _cdiv(num, den)
        lq = _clog(np.where(small, 1.0 + 0j, q))
    om = _cplx((-lq.imag) / code_dt, lq.real / code_dt)
    return np.where(small, 0.0 + 0.0j, om)


def omega_threshold(omegatol: float, kind: str) -> float:
    """``tiny(0.0)*fac`` with ``fac = 1000/abs(omegatol)`` (diagnostics_omega.f90:87-91, ``kind='omega'``)
    or ``tiny(0.0)*1000./abs(tol)`` (gs2_diagnostics.fpp:2944-2948, ``kind='inst'``)."""
    if omegatol == 0:
        return TINY
    if kind == "omega":
        return TINY * (1000.0 / abs(omegatol))
    return (TINY * 1000.0) / abs(omegatol)


def calculate_omega(state: "DiagState", phi_ig_old: np.ndarray, phi_ig_new: np.ndarray,
                    code_dt: float, istep: int, params: DiagParams) -> None:
    """``calculate_omega`` (src/diagnostics/diagnostics_omega.f90:64-136), host NumPy.

    ``phi_ig_old`` / ``phi_ig_new``: ``phi + apar + bpar`` at ``ig = igomega`` for the old and
    new time levels, shape ``(Y, X)`` complex (``((phi+apar)+bpar)``, all three fields
    unconditionally, as the Fortran does). Updates ``state.omegahist`` slot
    ``mod(istep, navg)``, ``state.omega_average`` and ``state.exit`` (convergence sets
    ``exit = exit_when_converged`` -- an assignment; instability sets ``exit = .true.``).
    """
    if not params.omega_active:
        return
    navg = params.navg
    k = istep % navg
    thresh = omega_threshold(params.omegatol, "omega")
    state.omegahist[k] = omega_quotient_np(phi_ig_new, phi_ig_old, code_dt, thresh)
    # omega_average = sum(omegahist/real(navg), dim=1)
    state.omega_average = _seq_sum(_crdiv(state.omegahist, float(navg)), axis=0)
    if istep > navg:
        dom = state.omega_average[None] - state.omegahist        # spread(...) - omegahist
        a = _cabs(dom)
        s = _seq_sum((a * a) / float(navg), axis=0)              # abs(domega)**2/real(navg)
        conv = np.all(np.sqrt(s) <= np.minimum(_cabs(state.omega_average), 1.0) * params.omegatol)
        if conv:
            state.exit = bool(params.exit_when_converged)
            state.converged_message = True
        if np.any(_cabs(state.omega_average) * code_dt > params.omegatinst):
            state.exit = True
            state.unstable = True


def calculate_omega_jax(omegahist, num, den, code_dt, istep, params: DiagParams):
    """Device variant of :func:`calculate_omega` (jittable, ``params`` static).

    ``omegahist``: ``(navg, Y, X)`` complex; ``num``/``den``: ``(Y, X)`` complex field sums
    at ``igomega``. Returns ``(omegahist, omega_average, converged, unstable)`` (booleans as
    0-d device arrays; ``converged`` already includes the ``istep > navg`` condition). Uses
    XLA's complex ``log`` and ``abs``, which differ from libm's at the ulp level (omega to
    about 1e-14 relative); the host path is bitwise.
    """
    navg = params.navg
    thresh = omega_threshold(params.omegatol, "omega")
    small = (jnp.abs(num) < thresh) | (jnp.abs(den) < thresh)
    q = _cdiv_jax(num, den)
    lq = jnp.log(jnp.where(small, 1.0 + 0j, q))
    om = lax.complex(div_exact(-lq.imag, code_dt), div_exact(lq.real, code_dt))
    om = jnp.where(small, jnp.zeros_like(om), om)
    k = jnp.mod(istep, navg)
    omegahist = lax.dynamic_update_index_in_dim(omegahist, om, k, 0)
    scaled = lax.complex(div_exact(omegahist.real, float(navg)),
                         div_exact(omegahist.imag, float(navg)))

    def add(c, x):
        return c + x, None
    avg, _ = lax.scan(add, jnp.zeros_like(om), scaled)
    dom = avg[None] - omegahist
    a = jnp.abs(dom)
    s, _ = lax.scan(add, jnp.zeros(a.shape[1:], a.dtype), (a * a) / float(navg))
    conv = jnp.all(jnp.sqrt(s) <= jnp.minimum(jnp.abs(avg), 1.0) * params.omegatol)
    unstable = jnp.any(jnp.abs(avg) * code_dt > params.omegatinst)
    late = istep > navg
    return omegahist, avg, jnp.logical_and(late, conv), jnp.logical_and(late, unstable)


# ============================================================================ device kernel


def _gadj(g, phi, bpar, direction, data: DiagData, params: DiagParams):
    """``g_adjust(g, phi, bpar, direction)`` (dist_fn_arrays.fpp:101-170)."""
    return g_adjust_direction(g, phi, bpar, direction, params.fphi, params.fbpar, data.gadj)


def _le_slice(x, ie, il):
    """``x[:, ie, il]`` of an ``(S, E, L, ...)`` array with traced indices."""
    x = lax.dynamic_index_in_dim(x, ie, axis=1, keepdims=False)
    return lax.dynamic_index_in_dim(x, il, axis=1, keepdims=False)


def _integrate_moments(make_terms, nterm: int, data: DiagData, params: DiagParams):
    """``integrate_moment`` of ``nterm`` integrands at once (le_grids.f90, integrate_moment_c34).

    ``make_terms(ie, il)`` returns the ``g_work(:, isgn, iglo)`` values of every term for one
    ``(ie, il)`` as ``(nterm, S, Y, X, 2, T)`` complex. The accumulation
    ``total(:,it,ik,is) = total + (w(ie,is)*wl(:,il))*(g_work(:,1)+g_work(:,2))`` runs over
    ``(ie, il)`` in the layout order of ``data.le_ie / le_il`` (spec 08 section 3.1.6).
    Returns ``(nterm, S, Y, X, T)``.
    """
    S, Y, X, T = params.nspec, params.naky, params.ntheta0, params.nth
    n = params.negrid * params.nlambda
    acc0 = jnp.zeros((nterm, S, Y, X, T), jnp.complex128)

    def body(k, acc):
        ie = data.le_ie[k]
        il = data.le_il[k]
        terms = make_terms(ie, il)
        s = terms[..., 0, :] + terms[..., 1, :]
        w_ie = lax.dynamic_index_in_dim(data.w, ie, axis=1, keepdims=False)       # (S,)
        wl_il = lax.dynamic_index_in_dim(data.wl, il, axis=0, keepdims=False)     # (T,)
        ww = ob(w_ie[:, None] * wl_il[None, :])                                    # (S, T)
        return acc + rmul(s, ww[None, :, None, None, :])

    return lax.fori_loop(0, n, body, acc0)


def _bc(x):
    """``(S,)`` species vector broadcast against ``(S, Y, X, 2, T)``."""
    return x[:, None, None, None, None]


def _flux_terms(hn, ho, ie, il, data: DiagData, params: DiagParams):
    """The ``g_work`` of ``flux`` (diagnostics_fluxes.fpp:574-738) and of the symmetric
    energy exchange (:1167-1185) for one ``(ie, il)``, with ``gnew``/``g`` holding h.

    Order of the returned terms is :data:`FLUX_TERMS` (filtered by the present fields).
    """
    h = _le_slice(hn, ie, il)                          # (S, Y, X, 2, T)
    aj0 = _le_slice(data.gadj.aj0, ie, il)[..., None, :]
    aj1 = _le_slice(data.gadj.aj1, ie, il)[..., None, :]
    vperp2 = _le_slice(data.gadj.vperp2, ie, il)[..., None, :]
    vpa = _le_slice(data.vpa, ie, il)
    vpac = _le_slice(data.vpac, ie, il)
    en = _bc(lax.dynamic_index_in_dim(data.energy, ie, axis=1, keepdims=False))
    out = []
    if params.has_phi:
        gw = rmul(h, aj0)                                            # gnew*aj0
        out.append(gw)                                               # pflux
        out.append(rmul(gw, en))                                     # *energy -> qflux(1)
        out.append(rmul(rmul(rmul(h, 2.0), ob(vpa * vpa)), aj0))     # gnew*2.*vpa**2*aj0
        out.append(rmul(rmul(h, vperp2), aj0))                       # gnew*vperp2*aj0
        # gnew*aj0*vpac*Rplot*sqrt(1.0-Bpol**2/bmag**2)
        out.append(rmul(rmul(rmul(rmul(h, aj0), vpac), data.rplot), data.sq_bpol))
        # -zi*aky*gnew*aj1*rhoc*(gds21+theta0*gds22)*vperp2*smz/(qval*shat*bmag**2)
        aky = data.aky[None, :, None, None, None]
        t = imul(aky, h)
        t = rmul(t, aj1)
        t = rmul(t, data.rhoc)
        t = rmul(t, data.gds_geo[None, :, :, None, :])
        t = rmul(t, vperp2)
        t = rmul(t, _bc(data.smz))
        t = lax.complex(div_exact(t.real, data.vperp_den), div_exact(t.imag, data.vperp_den))
        out.append(-t)                                               # vflux_perp
    if params.has_apar:
        base = rmul(rmul(rmul(h, aj0), _bc(data.stm)), vpa)          # gnew*aj0*stm*vpa
        gw = -base
        out.append(gw)                                               # pmflux
        out.append(rmul(gw, en))                                     # qmflux(1)
        out.append(-rmul(rmul(base, 2.0), ob(vpa * vpa)))            # ...*2.*vpa**2
        out.append(-rmul(base, vperp2))                              # ...*vperp2
        out.append(-rmul(base, vpac))                                # ...*vpac
    if params.has_bpar:
        base = rmul(rmul(rmul(rmul(h, aj1), 2.0), vperp2), _bc(data.tz))   # gnew*aj1*2.0*vperp2*tz
        out.append(base)                                             # pbflux
        out.append(rmul(base, en))                                   # qbflux(1)
        out.append(rmul(rmul(base, 2.0), ob(vpa * vpa)))             # *2.*vpa**2
        out.append(rmul(base, vperp2))                               # *vperp2
        out.append(rmul(base, vpac))                                 # *vpac
    if params.has_phi:
        o = _le_slice(ho, ie, il)
        a = ob(aj0 * 0.25)
        out.append(rmul(h + o, a))                                   # aj0*0.25*(gnew+g)
        out.append(rmul(h - o, a))                                   # aj0*0.25*(gnew-g)
    return jnp.stack(out)


def flux_term_names(params: DiagParams) -> list[str]:
    """Names of the :func:`_flux_terms` integrals, in order."""
    out = []
    if params.has_phi:
        out += ["pflux", "qflux1", "qflux2", "qflux3", "vflux_par", "vflux_perp"]
    if params.has_apar:
        out += ["pmflux", "qmflux1", "qmflux2", "qmflux3", "vmflux"]
    if params.has_bpar:
        out += ["pbflux", "qbflux1", "qbflux2", "qbflux3", "vbflux"]
    if params.has_phi:
        out += ["exch_total", "exch_total2"]
    return out


MOMENT_TERMS = ("density", "upar", "tpar", "tperp", "qparflux", "pperpj1", "qpperpj1", "ntot")


def _moment_terms(h, phinew, ie, il, data: DiagData, params: DiagParams):
    """The ``g_work`` of ``getmoms`` (diagnostics_moments.fpp:229-311) for one ``(ie, il)``,
    with ``gnew`` holding h. Order :data:`MOMENT_TERMS`."""
    hs = _le_slice(h, ie, il)
    aj0 = _le_slice(data.gadj.aj0, ie, il)[..., None, :]
    aj1 = _le_slice(data.gadj.aj1, ie, il)[..., None, :]
    vperp2 = _le_slice(data.gadj.vperp2, ie, il)[..., None, :]
    vpa = _le_slice(data.vpa, ie, il)
    en = _bc(lax.dynamic_index_in_dim(data.energy, ie, axis=1, keepdims=False))
    nm = _bc(lax.dynamic_index_in_dim(data.nonmaxw_corr, ie, axis=1, keepdims=False))
    g0 = rmul(hs, aj0)                                    # aj0*gnew
    g1 = rmul(g0, vpa)                                    # vpa*g_work
    g2 = rmul(g1, ob(2.0 * vpa))                          # 2.*vpa*g_work
    tp = rmul(rmul(hs, vperp2), aj0)                      # vperp2*gnew*aj0
    qp = rmul(rmul(rmul(hs, vpa), aj0), en)               # vpa*gnew*aj0*energy
    pp = rmul(rmul(rmul(rmul(hs, aj1), 2.0), vperp2), _bc(data.tz))   # gnew*aj1*2.0*vperp2*tz
    qq = rmul(pp, en)                                     # g_work*energy
    adb = rmul(rmul(phinew[None, :, :, None, :], _bc(data.zt)), nm)   # phinew*zt*nonmaxw_corr
    nt = g0 - adb                                         # aj0*gnew - ...
    return jnp.stack([g0, g1, g2, tp, qp, pp, qq, nt])


def _fve_integrals(hn, fields_new: Fields, data: DiagData, params: DiagParams):
    """Per-``(iglo, isgn)`` theta integrals of ``get_flux_vs_e`` (diagnostics_fluxes.fpp:874-881):
    ``trapz(theta, aimag(g_work*conjg(fld))*wl(:,il)*jacob)`` with the ``flux_vs_e`` g_work
    (:759-796). Returns a dict field -> ``(S, E, L, Y, X, 2)`` real."""
    aj0 = data.gadj.aj0[..., None, :]
    aj1 = data.gadj.aj1[..., None, :]
    vperp2 = data.gadj.vperp2[..., None, :]
    bc7 = lambda x: x[:, None, None, None, None, None, None]   # noqa: E731
    wl = data.wl[None, None, :, None, None, None, :]
    out = {}
    for name, on, fld in (("phi", params.has_phi, fields_new.phi),
                          ("apar", params.has_apar, fields_new.apar),
                          ("bpar", params.has_bpar, fields_new.bpar)):
        if not on:
            continue
        if name == "phi":
            gw = rmul(hn, aj0)
        elif name == "apar":
            gw = -rmul(rmul(rmul(hn, aj0), bc7(data.stm)), data.vpa)
        else:
            gw = rmul(rmul(rmul(rmul(hn, aj1), 2.0), vperp2), bc7(data.tz))
        f = fld[None, None, None, :, :, None, :]
        im = ob(gw.real * (-f.imag)) + ob(gw.imag * f.real)    # aimag(g*conjg(f))
        y = ob(ob(im * wl) * data.jacob)
        out[name] = trapz_jax(data.theta, y)
    return out


def _diagnostics_device_impl(gnew, g, fields_new: Fields, fields_old: Fields, data: DiagData,
                             params: DiagParams):
    """Device part of one write step: the ``g_adjust`` round trips and every velocity-space
    integral of the write block (spec 08 section 3.1.8).

    Mirrors, in this order, ``common_calculate_fluxes`` (diagnostics_fluxes.fpp:262-280:
    gnew and g from g to h, ``flux``, ``flux_vs_e``, ``eexchange``, back to g) when
    ``calculate_fluxes``, then ``getmoms`` (diagnostics_moments.fpp:223-329: gnew to h,
    moments, back) when ``write_moments``. Returns ``(gnew, g, results)``; ``results`` holds
    the raw ``integrate_moment`` outputs ``(S, Y, X, T)`` (before species normalisation) and
    the ``flux_vs_e`` theta integrals.
    """
    res = {}
    if params.calculate_fluxes:
        hn = _gadj(gnew, fields_new.phi, fields_new.bpar, FROM_G, data, params)
        ho = _gadj(g, fields_old.phi, fields_old.bpar, FROM_G, data, params)
        names = flux_term_names(params)
        if names:
            tot = _integrate_moments(lambda ie, il: _flux_terms(hn, ho, ie, il, data, params),
                                     len(names), data, params)
            for i, name in enumerate(names):
                res[name] = tot[i]
        for name, val in _fve_integrals(hn, fields_new, data, params).items():
            res[f"fve_{name}"] = val
        gnew = _gadj(hn, fields_new.phi, fields_new.bpar, TO_G, data, params)
        g = _gadj(ho, fields_old.phi, fields_old.bpar, TO_G, data, params)
    if params.write_moments:
        res.update(_getmoms_device(gnew, fields_new, data, params))
        gnew = _gadj(_gadj(gnew, fields_new.phi, fields_new.bpar, FROM_G, data, params),
                     fields_new.phi, fields_new.bpar, TO_G, data, params)
    return gnew, g, res


def _getmoms_device(gnew, fields_new: Fields, data: DiagData, params: DiagParams):
    """``integrate_moment`` calls of ``getmoms`` (diagnostics_moments.fpp:183-330) on device;
    the round trip on ``gnew`` is applied by the caller. Raw (unnormalised) results."""
    h = _gadj(gnew, fields_new.phi, fields_new.bpar, FROM_G, data, params)
    tot = _integrate_moments(lambda ie, il: _moment_terms(h, fields_new.phi, ie, il, data, params),
                             len(MOMENT_TERMS), data, params)
    return {f"mom_{name}": tot[i] for i, name in enumerate(MOMENT_TERMS)}


#: Jitted device part of a write step (``params`` static).
diagnostics_device = jax.jit(_diagnostics_device_impl, static_argnames=("params",))


def _final_moments_impl(gnew, fields_new: Fields, data: DiagData, params: DiagParams):
    """``getmoms`` at finalisation (diagnostics_final_routines.f90:330): moments and the
    round trip on gnew."""
    res = _getmoms_device(gnew, fields_new, data, params)
    gnew = _gadj(_gadj(gnew, fields_new.phi, fields_new.bpar, FROM_G, data, params),
                 fields_new.phi, fields_new.bpar, TO_G, data, params)
    return gnew, res


final_moments_device = jax.jit(_final_moments_impl, static_argnames=("params",))


# ============================================================================ host assembly


def _np_fields(fields: Fields) -> Fields:
    return Fields(*(np.asarray(x, np.complex128) for x in fields))


def field_record(fname: str, f: np.ndarray, params: DiagParams, h: HostData) -> tuple[dict, float]:
    """``nc_write_field`` (gs2_io.fpp:1318-1390) for one field: ``X2_by_mode`` (theta average of
    ``real(conjg(f)*f)``), ``X_t``, ``X2_by_kx``, ``X2_by_ky``, ``X0``,
    ``X_igomega_by_mode``, ``X2_by_mode``, ``X2``. Returns (record, X2)."""
    rec = {}
    f2 = host_flavg(_abs2_conj(f), h)                       # (Y, X)
    if getattr(params, f"write_{fname}_over_time"):
        rec[f"{fname}_t"] = f
    by_kx = average_ky(f2, h.aky)
    rec[f"{fname}2_by_kx"] = by_kx
    rec[f"{fname}2_by_ky"] = average_kx(f2, h.aky)
    x2 = float(_seq_sum(by_kx, axis=0))
    rec[f"{fname}0"] = f[..., params.ig]
    rec[f"{fname}_igomega_by_mode"] = f[..., params.ig]
    rec[f"{fname}2_by_mode"] = f2
    rec[f"{fname}2"] = np.float64(x2)
    return rec, x2


def instantaneous_omega(fold: Fields, fnew: Fields, code_dt: float, params: DiagParams):
    """``calculate_instantaneous_omega(ig=igomega, tolerance=omegatol)``
    (gs2_diagnostics.fpp:2897-2957): ``has_*``-guarded field sums."""
    ig = params.ig
    Y, X = params.naky, params.ntheta0
    fs = np.zeros((Y, X), np.complex128)
    fsn = np.zeros((Y, X), np.complex128)
    first = True
    for name, on in (("phi", params.has_phi), ("apar", params.has_apar), ("bpar", params.has_bpar)):
        if not on:
            continue
        o = getattr(fold, name)[..., ig]
        n = getattr(fnew, name)[..., ig]
        if first and name == "phi":
            fs, fsn = o.copy(), n.copy()
        else:
            fs = fs + o
            fsn = fsn + n
        first = False
    thresh = omega_threshold(params.omegatol, "inst")
    return omega_quotient_np(fsn, fs, code_dt, thresh)


def _weighted_kperp2(f, h: HostData):
    """``get_field_weighted_average_kperp2`` (gs2_diagnostics.fpp:2807-2824)."""
    a = _cabs(f)
    fsq = a * a
    num = host_flavg(fsq * h.kperp2, h)
    den = host_flavg(fsq, h)
    with np.errstate(all="ignore"):
        avg = num / den
    skip = np.max(np.abs(fsq), axis=-1) <= EPS
    return np.where(skip, 0.0, avg)


def ql_metric(fold: Fields, fnew: Fields, code_dt: float, params: DiagParams, h: HostData):
    """``write_ql_metric`` (gs2_diagnostics_new.f90:319-333): growth rate from the
    instantaneous omega and ``calculate_simple_quasilinear_flux_metric_by_k``
    (gs2_diagnostics.fpp:2827-2893)."""
    growth = instantaneous_omega(fold, fnew, code_dt, params).imag
    Y, X = params.naky, params.ntheta0
    ql = np.zeros((Y, X))
    if params.has_phi:
        norm = np.max(_cabs(fnew.phi), axis=-1)
    elif params.has_apar:
        norm = np.max(_cabs(fnew.apar), axis=-1)
    elif params.has_bpar:
        norm = np.max(_cabs(fnew.bpar), axis=-1)
    else:
        return ql
    with np.errstate(all="ignore"):
        norm = np.where(norm > 0, 1.0 / norm, norm)
    for name in ("phi", "apar", "bpar"):
        if not getattr(params, f"has_{name}"):
            continue
        f = getattr(fnew, name)
        avg = _weighted_kperp2(f, h)
        with np.errstate(all="ignore"):
            term = (np.max(_cabs(f), axis=-1) * norm) / avg
        if name == "phi":
            ql = np.where(avg > 0.0, term, ql)
        else:
            ql = np.where(avg > 0.0, ql + term, ql)
    lim = np.where(growth > 0, growth, 0.0)
    ql = np.where((h.aky == 0)[:, None], 0.0, ql)
    return ql * lim


def _get_flux(total, fld, h: HostData):
    """``get_flux`` (diagnostics_fluxes.fpp:801-838): ``total`` (S, Y, X, T) -> (S, Y, X)."""
    im = total.real * (-fld.imag)[None] + total.imag * fld.real[None]   # aimag(total*conjg(fld))
    flx = (host_trapz(h.theta, im * h.jacob) * h.flux_wgt) * h.aky[None, :, None]
    return flx * 0.5


def _fve_total(T, h: HostData, params: DiagParams):
    """``get_flux_vs_e`` accumulation (diagnostics_fluxes.fpp:865-887): ``T`` (S, E, L, Y, X, 2)
    -> ``0.5*total`` (S, E), contributions in ``iglo`` order, isgn 1 then 2."""
    fac = _ky_fac(h.aky)[None, None, None, :, None, None]
    aky = h.aky[None, None, None, :, None, None]
    term = ((fac * T) * aky) * h.flux_wgt
    axis = {"il": 2, "ik": 3, "it": 4}
    perm = (0, 1) + tuple(axis[o] for o in params.fve_order) + (5,)
    seq = np.transpose(term, perm).reshape(term.shape[0], term.shape[1], -1)
    return 0.5 * _seq_sum(seq, axis=-1)


def standard_flux_properties(name: str, F: np.ndarray, params: DiagParams, h: HostData) -> dict:
    """``calculate_standard_flux_properties`` (diagnostics_fluxes.fpp:434-531) for ``F`` (S, Y, X).

    Returns the written variables (if ``write_fluxes``) plus ``'_by_species'`` for the
    current-results bookkeeping."""
    by_species = _seq_sum(average_ky(F, h.aky), axis=-1)                  # average_all_real_xys
    by_mode = _seq_sum(F, axis=0)                                          # 0 + F(.,.,1) + ...
    by_ky = average_kx(by_mode, h.aky)
    by_kx = average_ky(by_mode, h.aky)
    rec = {"_by_species": by_species}
    if params.write_fluxes:
        rec[name] = by_species
        rec[f"total_{name}_by_ky"] = by_ky
        rec[f"total_{name}_by_kx"] = by_kx
        rec[f"total_{name}"] = np.float64(_seq_sum(by_kx, axis=0))
        if params.write_fluxes_by_mode:
            rec[f"{name}_by_mode"] = F
            rec[f"total_{name}_by_mode"] = by_mode
    return rec


def diffusive_estimates(omega_average, params: DiagParams, h: HostData) -> tuple[dict, float]:
    """``write_diffusive_estimates`` (diagnostics_fluxes.fpp:334-430), linear runs only."""
    ig = params.ig
    kp = h.kperp2[..., ig]                                           # (Y, X)
    gr = h.grho[ig]
    with np.errstate(all="ignore"):
        D = (np.maximum(omega_average.imag, 0.0) / kp) * 2.0
    D = np.where(kp == 0.0, 0.0, D)
    dens, temp = h.dens[:, None, None], h.temp[:, None, None]
    Q = (((D[None] * dens) * temp) * h.tprim[:, None, None]) * gr
    P = ((D[None] * (dens * dens)) * h.fprim[:, None, None]) * gr
    Q = np.where((kp == 0.0)[None], 0.0, Q)
    P = np.where((kp == 0.0)[None], 0.0, P)
    diffusivity = float(np.max(D)) * gr
    rec = {}
    for name, val in (("heat_flux_diff", Q), ("part_flux_diff", P)):
        r = standard_flux_properties(name, val, params, h)
        r.pop("_by_species")
        rec.update(r)
    hmax = np.max(np.max(Q, axis=2), axis=1)       # maxval over it, then over ik
    pmax = np.max(np.max(P, axis=2), axis=1)
    if params.write_fluxes:
        rec["es_heat_flux_diff_max"] = hmax
        rec["heat_flux_diff_max"] = np.float64(_seq_sum(hmax, axis=0))
        rec["es_particle_flux_diff_max"] = pmax
        rec["particle_flux_diff_max"] = np.float64(_seq_sum(pmax, axis=0))
        rec["diffusivity"] = np.float64(diffusivity)
    return rec, diffusivity


def fluxes_record(dev: dict, fold: Fields, fnew: Fields, omega_average, code_dt: float,
                  user_time: float, state: "DiagState", params: DiagParams, h: HostData) -> dict:
    """``calculate_fluxes`` (diagnostics_fluxes.fpp:101-242) from the device integrals."""
    rec = {}
    res = state.results
    res.update(total_heat_flux=0.0, total_momentum_flux=0.0, total_particle_flux=0.0)
    S = params.nspec
    species_heat = np.zeros(S)
    species_mom = np.zeros(S)
    species_part = np.zeros(S)
    if not params.nonlin:
        r, diffusivity = diffusive_estimates(omega_average, params, h)
        rec.update(r)
        res["diffusivity"] = diffusivity
    # common_calculate_fluxes: per-mode fluxes, then species normalisation (:284-329)
    dens = h.dens[:, None, None]
    DT = (h.dens * h.temp)[:, None, None]
    MOM = (h.dens * np.sqrt(h.mass * h.temp))[:, None, None]
    MOMP = ((h.dens * h.mass) * h.stm)[:, None, None]
    fl = {}
    fve = {}
    if params.has_phi:
        phi = fnew.phi
        q1 = _get_flux(dev["qflux1"], phi, h)
        q2 = _get_flux(dev["qflux2"], phi, h)
        q3 = _get_flux(dev["qflux3"], phi, h)
        pf = _get_flux(dev["pflux"], phi, h)
        vpar = _get_flux(dev["vflux_par"], phi, h)
        vperp = _get_flux(dev["vflux_perp"], phi, h)
        vfl = vpar + vperp                                   # vflux = vflux_par + vflux_perp
        # eexchange: exchange_sym (diagnostics_fluxes.fpp:1195-1206)
        t1, t2 = dev["exch_total"], dev["exch_total2"]
        dphi = fnew.phi - fold.phi
        sphi = fnew.phi + fold.phi
        a_re = t1.real * dphi.real[None] - t1.imag * (-dphi.imag)[None]   # real(total*conjg(dphi))
        b_re = sphi.real[None] * t2.real - sphi.imag[None] * (-t2.imag)   # real(sphi*conjg(total2))
        exch = host_flavg(a_re - b_re, h) / code_dt
        if params.nonlin:
            exch[:, 0, 0] = 0.0
        fl["es_heat_flux"] = q1 * DT
        fl["es_heat_flux_par"] = q2 * DT
        fl["es_heat_flux_perp"] = q3 * DT
        fl["es_part_flux"] = pf * dens
        fl["es_part_tormom_flux"] = np.zeros_like(pf) * dens      # no LOWFLOW: zeros
        fl["es_mom_flux"] = vfl * MOM
        fl["es_mom_flux_par"] = vpar * MOM
        fl["es_energy_exchange"] = (exch * dens) * h.z[:, None, None]
        fl["es_mom_flux_perp"] = vperp * MOMP
        fve["phi"] = _fve_total(dev["fve_phi"], h, params) * h.dens[:, None]
    if params.has_apar:
        apar = fnew.apar
        fl["apar_heat_flux"] = _get_flux(dev["qmflux1"], apar, h) * DT
        fl["apar_heat_flux_par"] = _get_flux(dev["qmflux2"], apar, h) * DT
        fl["apar_heat_flux_perp"] = _get_flux(dev["qmflux3"], apar, h) * DT
        fl["apar_part_flux"] = _get_flux(dev["pmflux"], apar, h) * dens
        fl["apar_mom_flux"] = _get_flux(dev["vmflux"], apar, h) * MOMP
        fve["apar"] = _fve_total(dev["fve_apar"], h, params) * h.dens[:, None]
    if params.has_bpar:
        bpar = fnew.bpar
        fl["bpar_heat_flux"] = _get_flux(dev["qbflux1"], bpar, h) * DT
        fl["bpar_heat_flux_par"] = _get_flux(dev["qbflux2"], bpar, h) * DT
        fl["bpar_heat_flux_perp"] = _get_flux(dev["qbflux3"], bpar, h) * DT
        fl["bpar_part_flux"] = _get_flux(dev["pbflux"], bpar, h) * dens
        fl["bpar_mom_flux"] = _get_flux(dev["vbflux"], bpar, h) * MOMP
        fve["bpar"] = _fve_total(dev["fve_bpar"], h, params) * h.dens[:, None]
    tot_h = tot_m = tot_p = 0.0
    for name, _desc, _units, fld in outnc.FLUX_PROPERTIES:
        if name not in fl:
            continue
        r = standard_flux_properties(name, fl[name], params, h)
        sp = r.pop("_by_species")
        rec.update(r)
        if name in ("es_heat_flux", "apar_heat_flux", "bpar_heat_flux"):
            tot_h = tot_h + float(_seq_sum(sp, axis=0))
            species_heat = species_heat + sp
        elif name in ("es_mom_flux", "apar_mom_flux", "bpar_mom_flux"):
            tot_m = tot_m + float(_seq_sum(sp, axis=0))
            species_mom = species_mom + sp
        elif name in ("es_part_flux", "apar_part_flux", "bpar_part_flux"):
            tot_p = tot_p + float(_seq_sum(sp, axis=0))
            species_part = species_part + sp
        if name == "es_heat_flux":
            res["species_es_heat_flux"] = sp
        if name == "es_energy_exchange":
            res["species_energy_exchange"] = sp
        if name in ("apar_heat_flux", "bpar_heat_flux"):
            res[f"species_{name}"] = sp
        if name == "es_mom_flux_perp" or (name == "apar_mom_flux") or (name == "bpar_mom_flux"):
            vname = outnc.FLUX_VS_E[fld][0]
            rec[vname] = fve[fld]
    res.update(total_heat_flux=tot_h, total_momentum_flux=tot_m, total_particle_flux=tot_p,
               species_heat_flux=species_heat, species_momentum_flux=species_mom,
               species_particle_flux=species_part)
    dt = user_time - state.user_time_old
    for key in ("heat", "particle", "momentum"):   # Trinity-only averages (:212-217)
        cur = res[f"species_{key}_flux"]
        res[f"species_{key}_flux_avg"] = res.get(f"species_{key}_flux_avg", np.zeros(S)) + cur * dt
    if params.write_fluxes:
        rec["heat_flux_tot"] = np.float64(tot_h)
        rec["hflux_tot"] = np.float64(tot_h)
        rec["mom_flux_tot"] = np.float64(tot_m)
        rec["vflux_tot"] = np.float64(tot_m)
        rec["part_flux_tot"] = np.float64(tot_p)
        rec["zflux_tot"] = np.float64(tot_p)
    return rec


def normalised_moments(dev: dict, h: HostData) -> dict:
    """Species normalisation of ``getmoms`` (diagnostics_moments.fpp:257-326), from the raw
    device integrals. Returns name -> (S, Y, X, T) complex."""
    m = {k: dev[f"mom_{k}"] for k in MOMENT_TERMS}
    tpar = m["tpar"] - m["density"]
    tperp = m["tperp"] - m["density"]
    b = lambda x: x[:, None, None, None]   # noqa: E731
    return {
        "ntot": _crmul(m["ntot"], b(h.dens)),
        "density": _crmul(m["density"], b(h.dens)),
        "upar": _crmul(m["upar"], b(h.stm)),
        "tpar": _crmul(tpar, b(h.temp)),
        "tperp": _crmul(tperp, b(h.temp)),
        "qparflux": _crmul(_crmul(_crmul(m["qparflux"], b(h.dens)), b(h.temp)), b(h.stm)),
        "pperpj1": _crmul(_crmul(m["pperpj1"], b(h.dens)), b(h.tz)),
        "qpperpj1": _crmul(_crmul(_crmul(m["qpperpj1"], b(h.dens)), b(h.temp)), b(h.tz)),
    }


def moments_record(moms: dict, params: DiagParams, h: HostData) -> dict:
    """``write_standard_moment_properties`` for the 8 moments (diagnostics_moments.fpp:83-180)."""
    rec = {}
    ig = params.ig
    for name, _desc, _units in outnc.MOMENTS:
        M = moms[name]
        m2 = host_flavg(_abs2_conj(M), h)                       # average_theta_txys (S, Y, X)
        rec[f"{name}2_by_ky"] = average_kx(m2, h.aky)            # (S, Y)
        by_kx = average_ky(m2, h.aky)                            # (S, X)
        rec[f"{name}2_by_kx"] = by_kx
        rec[f"{name}2"] = _seq_sum(by_kx, axis=-1)
        mbm = host_flavg_c(M, h)                                 # (S, Y, X) complex
        rec[f"{name}_flxsurf_avg"] = mbm[:, 0, :]
        rec[name] = M
        rec[f"{name}_igomega_by_mode"] = M[..., ig]
        rec[f"{name}_by_mode"] = mbm
        rec[f"{name}2_by_mode"] = m2
        flag = outnc.MOMENT_OVER_TIME.get(name)
        if flag and getattr(params, flag):
            rec[f"{name}_t"] = M
    return rec


def correlation(phinew, params: DiagParams):
    """``correlation`` (gs2_diagnostics.fpp:3389-3413): ``cf[ik, ig] = sum_it
    phinew(0,it,ik)*conjg(phinew(ig,it,ik))*fac``, ``fac = 1`` for ``ik == 1`` else 0.5."""
    p0 = phinew[..., params.ntgrid][..., None]                # theta index 0 hard-coded
    terms = _cmul_conj(np.broadcast_to(p0, phinew.shape), phinew)   # (Y, X, T)
    fac = np.where(np.arange(params.naky) == 0, 1.0, 0.5)[:, None, None]
    terms = _crmul(terms, fac)
    return _seq_sum(terms, axis=1)                            # (Y, T)


def get_phi0(fold: Fields, params: DiagParams, h: HostData):
    """``get_phi0`` (diagnostics_fields.f90:57-85), from the OLD fields."""
    i0 = params.ntgrid
    phi0 = fold.phi[..., i0].copy()
    if params.def_parity and params.fapar > 0 and not params.even:
        phi0 = fold.apar[..., i0].copy()
    alt = _crdiv(fold.phi[..., i0 + 1], h.theta[i0 + 1] - h.theta[i0])
    phi0 = np.where(_cabs(phi0) < 10.0 * EPS, alt, phi0)
    phi0 = np.where(_cabs(phi0) < 10.0 * EPS, 1.0 + 0.0j, phi0)
    return phi0


def eigenfunc_record(fold: Fields, params: DiagParams, h: HostData) -> dict:
    """``write_eigenfunc`` -> ``nc_eigenfunc`` (diagnostics_fields.f90:88-123, gs2_io.fpp:560-603)."""
    phase = get_phi0(fold, params, h)
    rec = {"phase": phase}
    for name in params.fields_present:
        f = getattr(fold, name)
        rec[f"{name}_norm"] = _cdiv(f, np.broadcast_to(phase[..., None], f.shape))
    return rec


def compute_record(dev: dict, fold: Fields, fnew: Fields, state: "DiagState", code_dt: float,
                   user_time: float, istep: int, params: DiagParams, h: HostData) -> dict:
    """Every variable of one write block (gs2_diagnostics_new.f90:434-475), host part.

    ``dev`` is the (host-converted) output of :func:`diagnostics_device`; ``fold``/``fnew``
    are NumPy ``Fields`` at the old/new time level. Returns ``{netCDF name: value}`` in
    GS2's write order (complex values complex128, without the ``ri`` axis).
    """
    rec = {"t": np.float64(user_time)}
    res = state.results
    if params.write_omega:      # nc_write_omega (gs2_io.fpp:1259-1294)
        wu = h.woutunits[:, None]
        k = istep % params.navg
        rec["omega"] = _crmul(state.omegahist[k], wu)
        avg = _crmul(state.omega_average, wu)
        rec["omegaavg"] = avg
        rec["omega_average"] = avg
    if params.write_fields:
        finals_done = False
        for name in params.fields_present:
            r, x2 = field_record(name, getattr(fnew, name), params, h)
            rec.update(r)
            res[f"{name}2"] = x2
            if not finals_done:
                for f in params.fields_present:
                    rec[f] = getattr(fnew, f)
                finals_done = True
    if params.write_ql_metric:
        rec["ql_metric_by_mode"] = ql_metric(fold, fnew, code_dt, params, h)
    if params.calculate_fluxes:
        rec.update(fluxes_record(dev, fold, fnew, state.omega_average, code_dt, user_time,
                                 state, params, h))
    if params.write_moments:
        rec.update(moments_record(normalised_moments(dev, h), params, h))
    if params.write_correlation and params.has_phi:
        rec["phi_corr_2pi"] = correlation(fnew.phi, params)
    if params.write_eigenfunc:
        rec.update(eigenfunc_record(fold, params, h))
    return rec


# ============================================================================ final routines


def get_epar(fold: Fields, fnew: Fields, code_dt: float, params: DiagParams, h: HostData):
    """``get_epar`` (diagnostics_final_routines.f90:201-230); ``epar(ntgrid)`` is left 0
    (uninitialised in Fortran)."""
    phi, apar, phin, aparn = fold.phi, fold.apar, fnew.phi, fnew.apar
    # phi_m = 0.5*(phi(ig+1)-phi(ig) + phinew(ig+1)-phinew(ig))*fphi
    phi_m = _crmul(_crmul((((phi[..., 1:] - phi[..., :-1]) + phin[..., 1:]) - phin[..., :-1]),
                          0.5), params.fphi)
    # apar_m = 0.5*(aparnew(ig+1)+aparnew(ig) - apar(ig+1)-apar(ig))*fapar
    apar_m = _crmul(_crmul((((aparn[..., 1:] + aparn[..., :-1]) - apar[..., 1:]) - apar[..., :-1]),
                           0.5), params.fapar)
    # epar = -phi_m/delthet(ig)*0.5*(abs(gradpar(ig)) + abs(gradpar(ig+1))) - apar_m/tunits(ik)/code_dt
    t1 = _crmul(_crmul(_crdiv(phi_m, h.delthet[:-1]), 0.5),
                np.abs(h.gradpar[:-1]) + np.abs(h.gradpar[1:]))
    t2 = _crdiv(_crdiv(apar_m, h.tunits[:, None, None]), code_dt)
    epar = np.zeros(phi.shape, np.complex128)
    epar[..., :-1] = (-t1) - t2
    return epar


# ============================================================================ state & control


@dataclass
class DiagState:
    """Host bookkeeping of the diagnostics between steps (spec 10 A.1).

    ``omegahist`` ``(navg, Y, X)`` complex and ``omega_average`` ``(Y, X)``: zero at t0 and
    not in restart files (P14). ``istep_last`` / ``nout``: 0 / 2 after the exporter wrote
    record 1 (G-05); use ``record1_written=False`` to let JAX write record 1 itself.
    """

    omegahist: np.ndarray
    omega_average: np.ndarray
    istep_last: int = 0
    nout: int = 2
    user_time_old: float = 0.0
    start_time: float = 0.0
    exit: bool = False
    converged_message: bool = False
    unstable: bool = False
    results: dict = field(default_factory=dict)

    @classmethod
    def initial(cls, params: DiagParams, *, record1_written: bool = True,
                start_time: float = 0.0) -> "DiagState":
        """State after ``init_gs2_diagnostics_new`` (``user_time_old = 0``,
        gs2_diagnostics_new.f90:103) and, if ``record1_written``, after the istep-0 write
        block of the exporter (``istep_last = 0``, ``nout = 2``,
        ``user_time_old = start_time = tstart``)."""
        Y, X = params.naky, params.ntheta0
        return cls(omegahist=np.zeros((params.navg, Y, X), np.complex128),
                   omega_average=np.zeros((Y, X), np.complex128),
                   istep_last=0 if record1_written else -1, nout=2 if record1_written else 1,
                   user_time_old=float(start_time) if record1_written else 0.0,
                   start_time=float(start_time))

    @classmethod
    def from_export(cls, ex, params: DiagParams) -> "DiagState":
        """State matching the exported ``/state_t0`` (post ``initialize_diagnostics``, G-05):
        ``nout``, ``user_time_old`` and ``start_time`` from ``/diagnostics``; record 1 written
        when the export says ``diagnostics_run_at_t0 = 1``."""
        dg = _group(ex, "/diagnostics")
        ran = bool(int(ex.attrs.get("diagnostics_run_at_t0", 1)))
        st = cls.initial(params, record1_written=ran,
                         start_time=_as_float(_lookup(dg, "start_time", 0.0)))
        if dg is not None and ran:
            st.nout = _as_int(_lookup(dg, "nout", st.nout))
            st.user_time_old = _as_float(_lookup(dg, "user_time_old", st.user_time_old))
            st.istep_last = 0
        return st


def _ig_sum(fields: Fields, ig: int) -> np.ndarray:
    """``phi(ig) + apar(ig) + bpar(ig)`` as ``((phi+apar)+bpar)``, host (Y, X) complex."""
    p = np.asarray(fields.phi[..., ig], np.complex128)
    a = np.asarray(fields.apar[..., ig], np.complex128)
    b = np.asarray(fields.bpar[..., ig], np.complex128)
    return (p + a) + b


def write_block(state: DiagState, gnew, g, fields_old: Fields, fields_new: Fields, *,
                istep: int, code_dt: float, user_time: float, data: DiagData,
                params: DiagParams, out: Optional[outnc.OutNC] = None,
                host: Optional[HostData] = None):
    """The write block of ``run_diagnostics`` (gs2_diagnostics_new.f90:432-508), unconditionally.

    Applies the round trips to ``gnew``/``g`` on the device, assembles the record on the host,
    appends it to ``out`` at ``state.nout`` and advances ``nout``. Returns
    ``(gnew, g, record)``.
    """
    host = host if host is not None else HostData.from_data(data)
    gnew, g, dev = diagnostics_device(gnew, g, fields_new, fields_old, data, params)
    dev = {k: np.asarray(v) for k, v in dev.items()}
    fold, fnew = _np_fields(fields_old), _np_fields(fields_new)
    rec = compute_record(dev, fold, fnew, state, code_dt, user_time, istep, params, host)
    if out is not None:
        out.add_specs(outnc.record_specs(params))
        out.write_record(state.nout, rec)
    state.user_time_old = user_time
    state.nout += 1
    return gnew, g, rec


def run_diagnostics(state: DiagState, gnew, g, fields_old: Fields, fields_new: Fields, *,
                    istep: int, code_dt: float, user_time: float, exit: bool,
                    data: DiagData, params: DiagParams, out: Optional[outnc.OutNC] = None,
                    host: Optional[HostData] = None, force: bool = False):
    """``run_diagnostics(istep, exit, force)`` (src/diagnostics/gs2_diagnostics_new.f90:339-519).

    Call after every accepted step (``istep >= 1``, after ``update_time``), with the state at
    diagnostic time (spec 08 section 1.2): ``gnew``/``fields_new`` at t(n+1), ``g``/
    ``fields_old`` at t(n), ``code_dt`` the dt of the step, ``user_time`` = t(n+1). Runs
    ``calculate_omega`` (every step), then the write block if
    ``istep /= istep_last .and. (mod(istep,nwrite)==0 .or. exit .or. force)``.

    Returns ``(gnew, g, exit, record)``; ``gnew`` and ``g`` carry the round-trip
    perturbation on write steps (and are the inputs unchanged otherwise); ``record`` is the
    written dict or None.
    """
    if not params.write_any:
        return gnew, g, exit, None
    state.exit = bool(exit)
    if istep == 0:
        state.start_time = user_time
    if istep > 0 and params.omega_active:
        ig = params.ig
        calculate_omega(state, _ig_sum(fields_old, ig), _ig_sum(fields_new, ig), code_dt,
                        istep, params)
    rec = None
    if istep != state.istep_last and (istep % params.nwrite == 0 or state.exit or force):
        state.istep_last = istep
        gnew, g, rec = write_block(state, gnew, g, fields_old, fields_new, istep=istep,
                                   code_dt=code_dt, user_time=user_time, data=data,
                                   params=params, out=out, host=host)
    # else: mod(istep, ncheck) == 0 only matters for write_verr (gated)
    return gnew, g, state.exit, rec


def finalize(state: DiagState, gnew, g, fields_old: Fields, fields_new: Fields, *,
             istep: int, code_dt: float, user_time: float, exit: bool, data: DiagData,
             params: DiagParams, out: Optional[outnc.OutNC] = None,
             host: Optional[HostData] = None,
             getan: Optional[Callable[[Any], tuple]] = None, close: bool = True):
    """End of run: the forced ``run_diagnostics`` and ``finish_gs2_diagnostics_new``.

    ``istep`` is the value GS2 passes to the forced call (gs2_main.fpp:870-881): ``nstep``
    after a normal loop end, the step of an explicit exit otherwise. Then
    ``run_old_final_routines`` (gs2_diagnostics_new.f90:544-602) in order: final fields
    (``phinew``...), ``epar``, final moments (one more ``g_adjust`` round trip on gnew),
    ``antot/antota/antotp`` from ``getan(gnew)`` -- pass ``getan`` (e.g. the fields module's
    ``getan`` bound to its data) when ``write_final_antot``. Restart files
    (``save_restart_dist_fn``) are the checkpoint module's job.

    Returns ``(gnew, g, exit, final_record)``.
    """
    host = host if host is not None else HostData.from_data(data)
    gnew, g, exit, _ = run_diagnostics(state, gnew, g, fields_old, fields_new, istep=istep,
                                       code_dt=code_dt, user_time=user_time, exit=exit,
                                       data=data, params=params, out=out, host=host, force=True)
    if not params.write_any:
        return gnew, g, exit, {}
    fold, fnew = _np_fields(fields_old), _np_fields(fields_new)
    final = {}
    if params.write_final_fields:
        for f in params.fields_present:
            final[f] = getattr(fnew, f)
    if params.write_final_epar:
        final["epar"] = get_epar(fold, fnew, code_dt, params, host)
    if params.write_final_moments:
        gnew, dev = final_moments_device(gnew, fields_new, data, params)
        moms = normalised_moments({k: np.asarray(v) for k, v in dev.items()}, host)
        for name, _d, _u in outnc.MOMENTS:
            final[name] = moms[name]
    if params.write_final_antot:
        if getan is None:
            raise ValueError("write_final_antot needs a getan callable (fields module)")
        antot, antota, antotp = getan(gnew)
        for name, f, val in (("antot", "phi", antot), ("antota", "apar", antota),
                             ("antotp", "bpar", antotp)):
            if f in params.fields_present:
                final[name] = np.asarray(val, np.complex128)
    if out is not None:
        out.add_specs(outnc.record_specs(params))
        specs = {s.name: s for s in outnc.final_specs(params)}
        for name in final:
            if name in out.ds.variables:   # keep the attributes of the existing variable
                specs.pop(name, None)
        out.add_specs(specs.values())
        out.write_static(final)
        if close:
            out.close()
    return gnew, g, exit, final


def _bind_fields_getan(ex) -> Callable:
    """``getan`` of the fields module (11_api section 8: ``getan(gnew, data, params)``) bound
    to the export, for ``write_final_antot`` (do_write_final_antot,
    diagnostics_final_routines.f90:419-456)."""
    from . import fields as _fields
    from .types import to_jax
    fdata = to_jax(_fields.FieldData.from_export(ex))
    fparams = _fields.FieldParams.from_export(ex)
    return lambda gnew: _fields.getan(gnew, fdata, fparams)


class Diagnostics:
    """Convenience runner holding data, params, state and the output file.

    >>> diag = Diagnostics.from_export(ex, out_path)              # record 1 from the exporter
    >>> gnew, g, exit = diag.step(istep, gnew, g, fields_old, fields_new, code_dt, user_time)
    >>> gnew, g = diag.finalize(istep, gnew, g, fields_old, fields_new, code_dt, user_time)
    """

    def __init__(self, data: DiagData, params: DiagParams, out_path: Optional[str | Path] = None,
                 *, record1_written: bool = True, start_time: float = 0.0,
                 getan: Optional[Callable] = None):
        from .types import to_jax
        self.data = to_jax(data)
        self.params = params
        self.host = HostData.from_data(data)
        self.state = DiagState.initial(params, record1_written=record1_written,
                                       start_time=start_time)
        self.out = None
        if out_path is not None and params.write_any:
            self.out = outnc.OutNC(out_path, outnc.record_specs(params),
                                   nc_sync_freq=params.nc_sync_freq)
        self.getan = getan

    @classmethod
    def from_export(cls, ex, out_path=None, **kw) -> "Diagnostics":
        """Runner for the export's t0 state: ``nout``/``user_time_old``/``start_time`` from
        ``/diagnostics`` (:meth:`DiagState.from_export`)."""
        params = DiagParams.from_export(ex)
        if params.write_final_antot and kw.get("getan") is None:
            kw["getan"] = _bind_fields_getan(ex)
        diag = cls(DiagData.from_export(ex), params, out_path, **kw)
        if "record1_written" not in kw and "start_time" not in kw:
            diag.state = DiagState.from_export(ex, params)
        return diag

    def step(self, istep, gnew, g, fields_old, fields_new, code_dt, user_time, exit=False):
        gnew, g, exit, _ = run_diagnostics(self.state, gnew, g, fields_old, fields_new,
                                           istep=istep, code_dt=code_dt, user_time=user_time,
                                           exit=exit, data=self.data, params=self.params,
                                           out=self.out, host=self.host)
        return gnew, g, exit

    def finalize(self, istep, gnew, g, fields_old, fields_new, code_dt, user_time, exit=False):
        gnew, g, exit, _ = finalize(self.state, gnew, g, fields_old, fields_new, istep=istep,
                                    code_dt=code_dt, user_time=user_time, exit=exit,
                                    data=self.data, params=self.params, out=self.out,
                                    host=self.host, getan=self.getan)
        self.out = None
        return gnew, g

    def close(self):
        if self.out is not None:
            self.out.close()
            self.out = None
