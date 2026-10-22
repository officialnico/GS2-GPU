"""Writer for GS2's ``<run>.out.nc`` (spec 08 section 4; 10 A.13 and G-05).

The Fortran exporter creates ``<run>.out.nc`` with the static part (metadata,
dimensions, grids, species, geometry, input file) and record 1 (istep 0),
using GS2's own ``gs2_io`` routines. This module opens that file in append
mode and writes records ``nout = 2, 3, ...`` and the "final" variables with
GS2's variable names, dimension names, dimension order and attributes.

Conventions (spec 08 section 0 and 4):

* netCDF dimension tuples are given in C order (the reverse of the Fortran
  ``dim_names`` lists in ``src/gs2_io.fpp``), e.g. Fortran
  ``[ri, kx, ky, t]`` is ``("t", "ky", "kx", "ri")`` here.
* complex data have a trailing ``ri`` axis (real, imaginary).
* every variable is ``double``; attributes are ``units`` (if any) then
  ``long_name`` (if any), the order in which ``neasyf_write`` creates them
  (externals/neasyf/src/neasyf.in.f90:640-648). Trailing blanks of the
  Fortran literals are dropped (``'Q_gB = '`` is stored as ``"Q_gB ="``),
  as observed in GS2's own output.
* variables are created without compression or an explicit fill value, as
  ``neasyf_write`` does with ``neasyf_default_compression = 0`` (the new
  diagnostics never call ``init_gs2_io``, which is the only place that turns
  compression on).
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional

import netCDF4
import numpy as np

# ----------------------------------------------------------------------------- descriptions
#: Field long names used by ``nc_write_field`` (diagnostics_fields.f90:23-39).
FIELD_LONG_NAME = {
    "phi": "Electrostatic potential",
    "apar": "Parallel magnetic potential",
    "bpar": "Parallel magnetic field",
}

#: (units, long_name) of the final fields written by ``nc_final_fields`` (gs2_io.fpp:680-706).
FINAL_FIELD_ATTRS = {
    "phi": ("T_ref/q_ref rho_ref/a_ref", "Electrostatic Potential"),
    "apar": (None, "Parallel Magnetic Potential"),
    "bpar": (None, "delta B Parallel"),
}

_APAR = " resulting from  perpendicular magnetic fluctuations"
_BPAR = " resulting from  parallel magnetic fluctuations"

#: (name, description, units, field) of every ``calculate_standard_flux_properties`` call,
#: in call order (diagnostics_fluxes.fpp:126-208). Units are the Fortran literals, trimmed.
FLUX_PROPERTIES = (
    ("es_heat_flux", "Turbulent flux of heat", "Q_gB =", "phi"),
    ("es_heat_flux_par", "Turbulent flux of parallel heat", "Q_gB =", "phi"),
    ("es_heat_flux_perp", "Turbulent flux of perpendicular heat", "Q_gB =", "phi"),
    ("es_part_flux", "Turbulent flux of particles", "n_r?", "phi"),
    ("es_part_tormom_flux", "Ask Jung-Pyo Lee...", "?", "phi"),
    ("es_mom_flux", "Flux of toroidal angular momentum", "Pi_gB =", "phi"),
    ("es_mom_flux_par", "Flux of the parallel component of the toroidal angular momentum",
     "Pi_gB =", "phi"),
    ("es_energy_exchange", "??", "Pi_gB =", "phi"),
    ("es_mom_flux_perp", "Flux of the perpendicular component of the toroidal angular momentum",
     "Pi_gB =", "phi"),
    ("apar_heat_flux", "Turbulent flux of heat" + _APAR, "Q_gB =", "apar"),
    ("apar_heat_flux_par", "Turbulent flux of parallel heat" + _APAR, "Q_gB =", "apar"),
    ("apar_heat_flux_perp", "Turbulent flux of perpendicular heat" + _APAR, "Q_gB =", "apar"),
    ("apar_part_flux", "Turbulent flux of particles" + _APAR, "TBC", "apar"),
    ("apar_mom_flux", "Turbulent flux of momentum" + _APAR, "Pi_gB =", "apar"),
    ("bpar_heat_flux", "Turbulent flux of heat" + _BPAR, "Q_gB =", "bpar"),
    ("bpar_heat_flux_par", "Turbulent flux of parallel heat" + _BPAR, "Q_gB =", "bpar"),
    ("bpar_heat_flux_perp", "Turbulent flux of perpendicular heat" + _BPAR, "Q_gB =", "bpar"),
    ("bpar_part_flux", "Turbulent flux of particles" + _BPAR, "TBC", "bpar"),
    ("bpar_mom_flux", "Turbulent flux of momentum" + _BPAR, "Pi_gB =", "bpar"),
)

#: Linear-only diffusive estimates (diagnostics_fluxes.fpp:394-400).
DIFFUSIVE_PROPERTIES = (
    ("heat_flux_diff", "Diffusive estimate of turbulent flux of heat", "Q_gB ="),
    ("part_flux_diff", "Diffusive estimate of turbulent flux of particles", "n_r?"),
)

#: ``*_flux_vs_e`` long names (diagnostics_fluxes.fpp:150-207), written right after the last
#: flux property of the field.
FLUX_VS_E = {
    "phi": ("es_flux_vs_e", "Electrostatic flux"),
    "apar": ("apar_flux_vs_e", "Perpendicular magnetic flux, as function of energy"),
    "bpar": ("bpar_flux_vs_e", "Parallel magnetic flux, as function of energy"),
}

#: (name, description, units) of the moments in ``write_moments`` order
#: (diagnostics_moments.fpp:32-49).
MOMENTS = (
    ("ntot", "The total perturbed species density", "n_ref"),
    ("density", "The non-adiabatic part of the perturbed species density", "n_ref"),
    ("upar", "The perturbed parallel flow", "v_ths"),
    ("tpar", "The perturbed parallel temperature", "T_ref"),
    ("tperp", "The perturbed perpendicular temperature", "T_ref"),
    ("qparflux", "The parallel heat flux", "n_ref T_ref v_thr"),
    ("pperpj1", "A modified perpendicular pressure which give the particle flux from bpar",
     "n_ref T_ref / q_ref"),
    ("qpperpj1",
     "A modified perpendicular pressure * energy which give the heat flux from bpar",
     "n_ref T_ref^2 / q_ref"),
)

#: Moments that can be written over time, with their switch (diagnostics_moments.fpp:164-168).
MOMENT_OVER_TIME = {"ntot": "write_ntot_over_time", "density": "write_density_over_time",
                    "upar": "write_upar_over_time", "tperp": "write_tperp_over_time"}

#: Long names of the final antot variables (gs2_io.fpp:800-826).
ANTOT_LONG_NAME = {
    "antot": "Velocity moment of g used in calculating phi",
    "antota": "Velocity space moment of g used in calculating apar",
    "antotp": "Velocity space moment of g used in calculating bpar",
}


@dataclass(frozen=True)
class VarSpec:
    """One ``<run>.out.nc`` variable: name, C-order dims, attributes, record variable or not."""

    name: str
    dims: tuple[str, ...]
    units: Optional[str] = None
    long_name: Optional[str] = None

    @property
    def time(self) -> bool:
        return bool(self.dims) and self.dims[0] == "t"

    @property
    def is_complex(self) -> bool:
        return bool(self.dims) and self.dims[-1] == "ri"


def _flux_property_specs(name, desc, units, by_mode) -> list[VarSpec]:
    """Variables of ``calculate_standard_flux_properties`` (diagnostics_fluxes.fpp:474-503)."""
    out = [
        VarSpec(name, ("t", "species"), units, desc + " averaged over kx and ky"),
        VarSpec(f"total_{name}_by_ky", ("t", "ky"), units,
                desc + " summed over species and averaged over kx"),
        VarSpec(f"total_{name}_by_kx", ("t", "kx"), units,
                desc + " summed over species and averaged over ky"),
        VarSpec(f"total_{name}", ("t",), units,
                desc + " summed over species and averaged over kx and ky"),
    ]
    if by_mode:
        out += [
            VarSpec(f"{name}_by_mode", ("t", "species", "ky", "kx"), units,
                    desc + " as a function of species, kx and ky"),
            VarSpec(f"total_{name}_by_mode", ("t", "ky", "kx"), units,
                    desc + " summed over species"),
        ]
    return out


def field_specs(field: str, over_time: bool) -> list[VarSpec]:
    """Variables of ``nc_write_field`` for one field (gs2_io.fpp:1318-1390), in write order."""
    ln = FIELD_LONG_NAME[field]
    out = []
    if over_time:
        out.append(VarSpec(f"{field}_t", ("t", "ky", "kx", "theta", "ri"), None,
                           ln + " over time"))
    out += [
        VarSpec(f"{field}2_by_kx", ("t", "kx"), None, ln + " squared and averaged over theta and ky"),
        VarSpec(f"{field}2_by_ky", ("t", "ky"), None, ln + " squared and averaged over theta and kx"),
        VarSpec(f"{field}0", ("t", "ky", "kx", "ri"), None, ln + " at ig=igomega"),
        VarSpec(f"{field}_igomega_by_mode", ("t", "ky", "kx", "ri"), None, ln + " at ig=igomega"),
        VarSpec(f"{field}2_by_mode", ("t", "ky", "kx"), None, ln + " per Fourier mode"),
        VarSpec(f"{field}2", ("t",), "(T_ref/q_ref rho_ref/a_ref)**2", ln + " volume average"),
    ]
    return out


def final_field_spec(field: str) -> VarSpec:
    units, ln = FINAL_FIELD_ATTRS[field]
    return VarSpec(field, ("ky", "kx", "theta", "ri"), units, ln)


def moment_specs(name: str, desc: str, units: str, over_time: bool) -> list[VarSpec]:
    """Variables of ``write_standard_moment_properties`` (diagnostics_moments.fpp:113-175)."""
    sq = "(" + units + ")^2"
    out = [
        VarSpec(f"{name}2_by_ky", ("t", "species", "ky"), sq,
                desc + " squared and averaged over theta and kx"),
        VarSpec(f"{name}2_by_kx", ("t", "species", "kx"), sq,
                desc + " squared and averaged over theta and ky"),
        VarSpec(f"{name}2", ("t", "species"), sq,
                desc + " squared and averaged over theta, kx and ky"),
        VarSpec(f"{name}_flxsurf_avg", ("t", "species", "kx", "ri"), units,
                desc + " flux surface averaged over theta, at ky=0 (actally ik==1)"),
        VarSpec(name, ("species", "ky", "kx", "theta", "ri"), units, desc),
        VarSpec(f"{name}_igomega_by_mode", ("t", "species", "ky", "kx", "ri"), units,
                desc + " at ig=igomega"),
        VarSpec(f"{name}_by_mode", ("t", "species", "ky", "kx", "ri"), units,
                desc + "  averaged over theta"),
        VarSpec(f"{name}2_by_mode", ("t", "species", "ky", "kx"), units,
                desc + " squared and averaged over theta"),
    ]
    if over_time:
        out.append(VarSpec(f"{name}_t", ("t", "species", "ky", "kx", "theta", "ri"), units,
                           desc + ": the whole moment"))
    return out


def record_specs(params) -> list[VarSpec]:
    """Every variable written by one ``run_diagnostics`` write block, in GS2's order.

    ``params`` is a ``diagnostics.DiagParams`` (any object with the same flag attributes).
    Mirrors the write sequence of ``run_diagnostics``
    (src/diagnostics/gs2_diagnostics_new.f90:434-475).
    """
    p = params
    fields = [f for f, on in (("phi", p.has_phi), ("apar", p.has_apar), ("bpar", p.has_bpar)) if on]
    out = [VarSpec("t", ("t",))]
    if p.write_omega:
        out += [
            VarSpec("omega", ("t", "ky", "kx", "ri"), "aref/vtref", "Complex frequency"),
            VarSpec("omegaavg", ("t", "ky", "kx", "ri"), "aref/vtref",
                    "Complex frequency averaged over navg timesteps"),
            VarSpec("omega_average", ("t", "ky", "kx", "ri"), "aref/vtref",
                    "Complex frequency averaged over navg timesteps"),
        ]
    if p.write_fields:
        over = {"phi": p.write_phi_over_time, "apar": p.write_apar_over_time,
                "bpar": p.write_bpar_over_time}
        finals_done = False
        for f in fields:
            out += field_specs(f, over[f])
            if not finals_done:   # nc_final_fields writes every present field (gs2_io.fpp:1384)
                out += [final_field_spec(x) for x in fields]
                finals_done = True
    if p.write_ql_metric:
        out.append(VarSpec("ql_metric_by_mode", ("t", "ky", "kx"), None,
                           "QL flux metric by mode number"))
    if p.calculate_fluxes:
        if not p.nonlin and p.write_fluxes:
            for name, desc, units in DIFFUSIVE_PROPERTIES:
                out += _flux_property_specs(name, desc, units, p.write_fluxes_by_mode)
            out += [
                VarSpec("es_heat_flux_diff_max", ("t", "species"), "Q_gB",
                        "Linear estimate of the heat flux"),
                VarSpec("heat_flux_diff_max", ("t",), "Q_gB", "Linear estimate of the heat flux"),
                VarSpec("es_particle_flux_diff_max", ("t", "species"), "Q_gB",
                        "Linear estimate of the particle flux"),
                VarSpec("particle_flux_diff_max", ("t",), "Q_gB",
                        "Linear estimate of the particle flux"),
                VarSpec("diffusivity", ("t",), None, "Linear estimate of the diffusivity"),
            ]
        for f in fields:
            if p.write_fluxes:
                for name, desc, units, fld in FLUX_PROPERTIES:
                    if fld == f:
                        out += _flux_property_specs(name, desc, units, p.write_fluxes_by_mode)
            vname, ln = FLUX_VS_E[f]
            out.append(VarSpec(vname, ("t", "species", "egrid"), None, ln))
        if p.write_fluxes:
            out += [
                VarSpec("heat_flux_tot", ("t",), "Q_gB", "Total heat flux"),
                VarSpec("hflux_tot", ("t",), "Q_gB",
                        "Total heat flux, same as heat_flux_tot, included for backwards compatiblity"),
                VarSpec("mom_flux_tot", ("t",), "Pi_gB", "Total momentum flux"),
                VarSpec("vflux_tot", ("t",), "Pi_gB",
                        "Total momentum flux, same as mom_flux_tot, included for backwards compatiblity"),
                VarSpec("part_flux_tot", ("t",), None, "Total particle flux"),
                VarSpec("zflux_tot", ("t",), None,
                        "Total particle flux, same as part_flux_tot, included for backwards compatiblity"),
            ]
    if p.write_moments:
        for name, desc, units in MOMENTS:
            flag = MOMENT_OVER_TIME.get(name)
            out += moment_specs(name, desc, units, bool(flag and getattr(p, flag)))
    if p.write_correlation and p.has_phi:
        out.append(VarSpec("phi_corr_2pi", ("t", "ky", "theta", "ri"), None,
                           "2 point correlation function calculated from the electric potential"))
    if p.write_eigenfunc:
        out.append(VarSpec("phase", ("ky", "kx", "ri"), None, "Normalising phase"))
        out += [VarSpec(f"{f}_norm", ("ky", "kx", "theta", "ri")) for f in fields]
    return out


def final_specs(params) -> list[VarSpec]:
    """Variables of ``run_old_final_routines`` (gs2_diagnostics_new.f90:544-602)."""
    p = params
    fields = [f for f, on in (("phi", p.has_phi), ("apar", p.has_apar), ("bpar", p.has_bpar)) if on]
    out = []
    if p.write_final_fields:
        out += [final_field_spec(f) for f in fields]
    if p.write_final_epar:
        out.append(VarSpec("epar", ("ky", "kx", "theta", "ri")))          # gs2_io.fpp:709-721
    if p.write_final_moments:   # nc_final_moments writes no attributes (gs2_io.fpp:725-746)
        out += [VarSpec(name, ("species", "ky", "kx", "theta", "ri")) for name, _, _ in MOMENTS]
    if p.write_final_antot:
        for name, f in (("antot", "phi"), ("antota", "apar"), ("antotp", "bpar")):
            if f in fields:
                out.append(VarSpec(name, ("ky", "kx", "theta", "ri"), None, ANTOT_LONG_NAME[name]))
    return out


# ----------------------------------------------------------------------------- writer
def _to_disk(value, spec: VarSpec):
    arr = np.asarray(value)
    if spec.is_complex:
        if not np.iscomplexobj(arr):
            raise ValueError(f"{spec.name}: complex variable needs a complex array")
        arr = np.stack([arr.real, arr.imag], axis=-1).astype(np.float64)
    else:
        if np.iscomplexobj(arr):
            raise ValueError(f"{spec.name}: real variable got a complex array")
        arr = arr.astype(np.float64)
    return arr


class OutNC:
    """Append records to a ``<run>.out.nc`` (opened in ``"a"`` mode).

    ``specs`` maps variable name to :class:`VarSpec`; variables missing from the file are
    created with GS2's dimensions and attributes (the exporter normally creates them all in
    record 1, so in practice only the final-only variables are created here).
    """

    def __init__(self, path: str | Path, specs: Iterable[VarSpec] = (), *, nc_sync_freq: int = 10):
        self.path = Path(path)
        if not self.path.exists():
            raise FileNotFoundError(
                f"{self.path}: the exporter must create <run>.out.nc (static part and record 1) "
                "before JAX appends to it; see OutNC.create_minimal for synthetic tests")
        self.ds = netCDF4.Dataset(self.path, "a")
        self.ds.set_auto_mask(False)
        self.specs: dict[str, VarSpec] = {s.name: s for s in specs}
        self.nc_sync_freq = int(nc_sync_freq)

    # -- context manager
    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()

    def close(self):
        if self.ds is not None and self.ds.isopen():
            self.ds.close()
        self.ds = None

    def sync(self):
        self.ds.sync()

    @property
    def nrecords(self) -> int:
        return len(self.ds.dimensions["t"]) if "t" in self.ds.dimensions else 0

    def add_specs(self, specs: Iterable[VarSpec]):
        for s in specs:
            self.specs[s.name] = s

    def _variable(self, spec: VarSpec):
        var = self.ds.variables.get(spec.name)
        if var is not None:
            if tuple(var.dimensions) != spec.dims:
                raise ValueError(f"{self.path}: variable {spec.name} has dims {var.dimensions}, "
                                 f"expected {spec.dims}")
            return var
        missing = [d for d in spec.dims if d not in self.ds.dimensions]
        if missing:
            raise ValueError(f"{self.path}: cannot create {spec.name}: missing dims {missing}")
        var = self.ds.createVariable(spec.name, "f8", spec.dims)
        if spec.units is not None:
            var.setncattr("units", spec.units)
        if spec.long_name is not None:
            var.setncattr("long_name", spec.long_name)
        return var

    def write(self, name: str, value, nout: Optional[int] = None):
        """Write one variable; record variables go to record ``nout`` (1-based, Fortran)."""
        spec = self.specs.get(name)
        if spec is None:
            raise KeyError(f"no VarSpec for {name}")
        var = self._variable(spec)
        data = _to_disk(value, spec)
        if spec.time:
            if name == "t":
                var[nout - 1] = data
            else:
                var[nout - 1, ...] = data
        else:
            var[...] = data

    def write_record(self, nout: int, values: dict):
        """Write every entry of ``values`` (record variables at ``nout``), in dict order, and
        sync every ``nc_sync_freq`` records as ``nc_sync`` does (gs2_io.fpp:1558-1598)."""
        if "t" in values:
            self.write("t", values["t"], nout)
        for name, value in values.items():
            if name != "t":
                self.write(name, value, nout)
        if self.nc_sync_freq > 0 and nout % self.nc_sync_freq == 0:
            self.sync()

    def write_static(self, values: dict):
        for name, value in values.items():
            self.write(name, value)

    # -- file creation helpers (tests / standalone use only; production files come from the
    #    exporter)
    @staticmethod
    def create_minimal(path: str | Path, *, akx, aky, theta, nspec: int, negrid: int,
                       al=None) -> Path:
        """Create a bare ``<run>.out.nc`` with the dimensions ``define_dims`` makes
        (gs2_io.fpp:297-315) and their coordinate variables, no records. For synthetic tests."""
        path = Path(path)
        with netCDF4.Dataset(path, "w") as ds:
            ds.title = "GS2 Simulation Data"
            for name, values, units, ln in (
                    ("kx", akx, "1/rho_r", "Wavenumber perpendicular to flux surface"),
                    ("ky", aky, "1/rho_r", "Wavenumber in direction of grad alpha"),
                    ("theta", theta, "rad", "Poloidal angle")):
                ds.createDimension(name, len(values))
                v = ds.createVariable(name, "f8", (name,))
                v.units = units
                v.long_name = ln
                v[:] = np.asarray(values)
            if al is not None:
                ds.createDimension("lambda", len(al))
                v = ds.createVariable("lambda", "f8", ("lambda",))
                v.units = "1/(2 B_a)"
                v.long_name = "Energy/magnetic moment"
                v[:] = np.asarray(al)
            ds.createDimension("species", nspec)
            ds.createVariable("species", "i4", ("species",))[:] = np.arange(1, nspec + 1)
            ds.createDimension("egrid", negrid)
            v = ds.createVariable("egrid", "i4", ("egrid",))
            v.long_name = "See 'energy' for energy grid values"
            v[:] = np.arange(1, negrid + 1)
            ds.createDimension("t", None)
            ds.createDimension("ri", 2)
            v = ds.createVariable("ri", "i4", ("ri",))
            v.units = "(real, imaginary)"
            v.long_name = "Complex components"
            v[:] = [1, 2]
            ds.createVariable("t", "f8", ("t",))
        return path


def copy_records(src: str | Path, dst: str | Path, nrec: int = 1,
                 exclude: Iterable[str] = ()) -> Path:
    """Copy ``src`` (a GS2 ``<run>.out.nc``) to ``dst`` keeping only the first ``nrec`` records.

    Reproduces what the exporter writes (static part plus record 1) from GS2's own output, so
    the JAX writer can be tested before the exporter exists. Dimensions, variables (dims,
    dtype, attributes in order) and global attributes are copied; variables named in
    ``exclude`` are left out (so that the writer has to create them).
    """
    src, dst = Path(src), Path(dst)
    exclude = set(exclude)
    with netCDF4.Dataset(src) as s, netCDF4.Dataset(dst, "w", format=s.data_model) as d:
        s.set_auto_mask(False)
        d.set_auto_mask(False)
        d.setncatts({k: s.getncattr(k) for k in s.ncattrs()})
        for name, dim in s.dimensions.items():
            d.createDimension(name, None if dim.isunlimited() else len(dim))
        for name, var in s.variables.items():
            if name in exclude:
                continue
            kw = {}
            if "_FillValue" in var.ncattrs():
                kw["fill_value"] = var.getncattr("_FillValue")
            v = d.createVariable(name, var.dtype, var.dimensions, **kw)
            v.setncatts({k: var.getncattr(k) for k in var.ncattrs() if k != "_FillValue"})
            if var.dimensions and var.dimensions[0] == "t":
                n = min(nrec, var.shape[0])
                if n > 0:
                    v[:n, ...] = var[:n, ...]
            else:
                v[...] = var[...]
    return dst


# ----------------------------------------------------------------------------- format check
_IGNORED_ATTRS = {"_FillValue", "_Storage", "_ChunkSizes", "_Endianness", "_DeflateLevel",
                  "_Shuffle", "_NoFill"}


def check_format(got_path: str | Path, ref_path: str | Path,
                 names: Optional[Iterable[str]] = None, *, attr_order: bool = False) -> list[str]:
    """Problems that would stop GS2 analysis scripts from reading ``got`` like ``ref``.

    For every variable of ``names`` (default: every variable of ``ref``) checks presence,
    dimension names and order, non-time dimension sizes, dtype and the attributes (names and
    values; their order too if ``attr_order``: HDF5 occasionally reports the creation order
    of a variable created in append mode reversed, which no reader depends on). Returns a
    list of messages (empty if OK).
    """
    problems = []
    with netCDF4.Dataset(ref_path) as ref, netCDF4.Dataset(got_path) as got:
        names = list(ref.variables) if names is None else list(names)
        for name in names:
            if name not in ref.variables:
                problems.append(f"{name}: not in reference file")
                continue
            rv = ref.variables[name]
            if name not in got.variables:
                problems.append(f"{name}: missing")
                continue
            gv = got.variables[name]
            if tuple(gv.dimensions) != tuple(rv.dimensions):
                problems.append(f"{name}: dims {gv.dimensions} != {rv.dimensions}")
                continue
            for dname, rlen, glen in zip(rv.dimensions, rv.shape, gv.shape):
                if dname != "t" and rlen != glen:
                    problems.append(f"{name}: size of {dname} {glen} != {rlen}")
            if gv.dtype != rv.dtype:
                problems.append(f"{name}: dtype {gv.dtype} != {rv.dtype}")
            ra = [(k, rv.getncattr(k)) for k in rv.ncattrs() if k not in _IGNORED_ATTRS]
            ga = [(k, gv.getncattr(k)) for k in gv.ncattrs() if k not in _IGNORED_ATTRS]
            if sorted(k for k, _ in ra) != sorted(k for k, _ in ga):
                problems.append(f"{name}: attributes {[k for k, _ in ga]} != {[k for k, _ in ra]}")
            else:
                if attr_order and [k for k, _ in ra] != [k for k, _ in ga]:
                    problems.append(f"{name}: attribute order {[k for k, _ in ga]} != "
                                    f"{[k for k, _ in ra]}")
                gd = dict(ga)
                for k, rval in ra:
                    if str(rval) != str(gd[k]):
                        problems.append(f"{name}: attribute {k} {gd[k]!r} != {rval!r}")
        for dname, rdim in ref.dimensions.items():
            if dname not in got.dimensions:
                problems.append(f"dimension {dname}: missing")
            elif rdim.isunlimited() != got.dimensions[dname].isunlimited():
                problems.append(f"dimension {dname}: unlimited mismatch")
    return problems
