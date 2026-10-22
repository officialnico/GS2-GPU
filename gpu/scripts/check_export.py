#!/usr/bin/env python3
"""Check a GS2 GPU export file (<case>.gs2export.nc) and, optionally, its trace file.

Checks (gpu/spec/09_export_infra.md, gpu/spec/10_critique.md §A, gpu/spec/09a_exporter_notes.md):

1. Schema: every variable of the master export list (10 §A) is present with the expected NumPy
   shape and dtype, or listed in its group's `absent` attribute (arrays not allocated in GS2).
   Every variable written carries the standard attributes (fortran_type, canonical).
2. Numerical spot checks: recompute exported arrays from other exported inputs with the GS2
   expressions (vperp2, vpa, vpac, vpar, wstar, tunits/wunits, kperp2, field_line_average_weight,
   has_*, the layout maps iglo_of/ile_of, the dt-ladder rungs) and compare geometry, grids and the
   t = 0 fields with GS2's own <run>.out.nc. Expressions without an FMA pattern must agree bitwise;
   FMA-prone ones bitwise for a reference build (-ffp-contract=off), within 4 ulp otherwise.
3. Trace invariants (09 §6.2): step_end(k) == advance_begin(k+1, attempt 1); pass-1 phinew == phi;
   after_timeadv pass 2 (last attempt) == advance_end; after_field_update == after_getfield + phi;
   the nonlinear history shuffle; update_time; advance_begin(1) == /state_t0; ...

Usage:
    check_export.py CASE.gs2export.nc [--trace CASE.gs2trace.nc] [--outnc EXPORT_RUN.out.nc]
                    [--ref-outnc GS2_RUN.out.nc] [--trace-outnc TRACE_RUN.out.nc]
                    [--compare-export OTHER.gs2export.nc] [--quiet]
Exit status 0 if every check passes.
"""
import argparse
import os
import sys

import numpy as np
import netCDF4

# ----------------------------------------------------------------------------------- helpers


class Report:
    def __init__(self, quiet=False):
        self.quiet = quiet
        self.n_pass = 0
        self.fails = []
        self.notes = []

    def ok(self, cond, what, detail=""):
        if cond:
            self.n_pass += 1
            if not self.quiet:
                print(f"  PASS  {what}" + (f"  ({detail})" if detail else ""))
        else:
            self.fails.append(f"{what}" + (f": {detail}" if detail else ""))
            print(f"  FAIL  {what}" + (f": {detail}" if detail else ""))
        return cond

    def note(self, msg):
        self.notes.append(msg)
        if not self.quiet:
            print(f"  NOTE  {msg}")


def bits_equal(a, b):
    a = np.ascontiguousarray(np.asarray(a))
    b = np.ascontiguousarray(np.asarray(b))
    return a.shape == b.shape and a.dtype == b.dtype and a.tobytes() == b.tobytes()


def max_ulp(a, b):
    """Max distance in units of the last place between two float64 arrays."""
    a = np.asarray(a, dtype=np.float64).ravel()
    b = np.asarray(b, dtype=np.float64).ravel()
    if a.size == 0:
        return 0
    ia = a.view(np.int64)
    ib = b.view(np.int64)
    ia = np.where(ia < 0, np.int64(-0x8000000000000000) - ia, ia)
    ib = np.where(ib < 0, np.int64(-0x8000000000000000) - ib, ib)
    return int(np.max(np.abs(ia - ib)))


def describe_diff(a, b):
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    if a.shape != b.shape:
        return f"shape {a.shape} vs {b.shape}"
    d = np.abs(a - b)
    return f"max|diff|={np.max(d):.3e} max|ref|={np.max(np.abs(b)):.3e} ulp={max_ulp(a, b)}"


def rd(grp, name):
    v = grp.variables[name]
    a = v[...]
    a = np.asarray(a)
    return a


def cplx(a):
    return a[..., 0] + 1j * a[..., 1]


def walk(grp, prefix=""):
    yield prefix or "/", grp
    for name, sub in grp.groups.items():
        yield from walk(sub, f"{prefix}/{name}")


def absent_set(grp):
    if "absent" in grp.ncattrs():
        return set(str(grp.getncattr("absent")).split())
    return set()


def klein_sum(arr):
    """Klein (second-order Neumaier) sum, summation.f90:71-103, strictly sequential."""
    res = 0.0
    c1 = 0.0
    c2 = 0.0
    for x in np.asarray(arr, dtype=np.float64).tolist():
        t = res + x
        if abs(res) >= abs(x):
            lc1 = (res - t) + x
        else:
            lc1 = (x - t) + res
        res = t
        t = c1 + lc1
        if abs(c1) >= abs(lc1):
            lc2 = (c1 - t) + lc1
        else:
            lc2 = (lc1 - t) + c1
        c1 = t
        c2 = c2 + lc2
    return (res + c1) + c2


def trapz_gs2(x, y):
    """integration.f90:55-62: klein_sum((x(2:)-x(:n-1))*(y(:n-1)+y(2:))) * 0.5"""
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    return klein_sum((x[1:] - x[:-1]) * (y[:-1] + y[1:])) * 0.5


LAYOUT_ORDS = {  # gs2_layouts.fpp:917-930: (it_ord, ik_ord, il_ord, ie_ord, is_ord)
    "lxyes": (2, 3, 1, 4, 5),
    "xyles": (1, 2, 3, 4, 5),
    "yxles": (2, 1, 3, 4, 5),
    "yxels": (2, 1, 4, 3, 5),
    "lexys": (3, 4, 1, 2, 5),
    "lyxes": (3, 2, 1, 4, 5),
}


# ----------------------------------------------------------------------------------- schema

def schema(d, ex):
    """Expected variables (group -> list of (name, shape tuple, dtype char)), 10 §A."""
    S, E, L, Y, X, T = d["S"], d["E"], d["L"], d["Y"], d["X"], d["T"]
    E1, Em1, L2, XI, XP, C = d["E1"], d["Em1"], d["L2"], d["XI"], d["XP"], d["C"]
    f, i, b = "f8", "i4", "i1"
    G7 = (S, E, L, Y, X, 2, T)
    G7c = G7 + (2,)
    G6 = (S, E, L, Y, X, T)
    G5 = (S, E, L, Y, X)
    F3 = (Y, X, T)
    F3c = F3 + (2,)
    K2 = (Y, X)
    sc = ()
    tg = [(n, (T,), f) for n in ("theta", "delthet", "bmag", "gradpar", "gbdrift", "gbdrift0",
                                 "cvdrift", "cvdrift0", "cdrift", "cdrift0", "itor_over_B", "gds2",
                                 "gds21", "gds22", "grho", "jacob", "Rplot", "Bpol")]
    tg += [(n, sc, f) for n in ("field_line_average_weight", "shat", "kxfac", "qval", "rhoc", "drhodpsi")]
    tg += [(n, sc, i) for n in ("ntheta", "ntgrid", "nperiod")] + [("is_effectively_zero_shear", sc, b)]
    kt = [("aky", (Y,), f), ("akx", (X,), f), ("ikx", (X,), i), ("theta0", K2, f), ("kperp2", F3, f),
          ("inv_kperp2", F3, f), ("nx", sc, i), ("ny", sc, i), ("box", sc, b), ("reality", sc, b),
          ("gridopt_switch", sc, i), ("kwork_filter", K2, b)]
    le = [("energy", (S, E), f), ("w", (S, E), f), ("energy_maxw", (E,), f), ("speed_maxw", (E,), f),
          ("w_maxw", (E,), f), ("al", (L,), f), ("wl", (L, T), f), ("forbid", (L, T), b),
          ("is_ttp", (L, T), b), ("is_bounce_point", (L, T), b), ("can_be_ttp", (L,), b),
          ("is_lower_bounce_point", (L, T), b), ("is_upper_bounce_point", (L, T), b), ("jend", (T,), i),
          ("xi", (T, L2), f), ("wxi", (T, L2), f), ("ixi_to_il", (T, L2), i), ("ixi_to_isgn", (T, L2), i),
          ("ng2", sc, i), ("nlambda", sc, i), ("lmax", sc, i), ("nxi", sc, i), ("negrid", sc, i),
          ("passing_wfb", sc, b), ("trapped_wfb", sc, b), ("mixed_wfb", sc, b), ("new_trap_int", sc, b),
          ("bouncefuzz", sc, f), ("il_is_wfb", (L,), b), ("il_is_passing", (L,), b),
          ("il_is_trapped", (L,), b), ("is_passing_hybrid_electron", (S, Y, L), b)]
    sp = [(f"spec_{c}", (S,), f) for c in ("z", "mass", "dens", "temp", "tprim", "fprim", "uprim",
                                            "uprim2", "vnewk", "nu_h", "bess_fac", "stm", "zstm", "tz",
                                            "zt", "smz")]
    sp += [("spec_type", (S,), i), ("spec_f0type", (S,), i), ("nonmaxw_corr", (S, E), f),
           ("dlnf0drho", (S, E), f), ("has_electron_species", sc, b), ("has_ion_species", sc, b),
           ("has_hybrid_electron_species", sc, b), ("nspec", sc, i)]
    rp = [(n, sc, f) for n in ("fphi", "fapar", "fbpar", "beta", "tite", "zeff", "max_sim_time")]
    rp += [(n, sc, b) for n in ("has_phi", "has_apar", "has_bpar", "immediate_reset", "wstar_units", "reset")]
    rp += [("ieqzip", K2, b), ("nstep", sc, i), ("ncheck_stop", sc, i), ("delt_option_switch", sc, i)]
    tm = [(n, (Y,), f) for n in ("tunits", "wunits", "woutunits")]
    timeblock = [(n, sc, f) for n in ("code_dt", "user_dt", "code_dt_prev1", "code_dt_prev2",
                                      "code_dt_old", "code_time", "user_time", "code_dt_cfl",
                                      "user_dt_cfl", "code_dt_min", "code_dt_max")]
    timeblock += [("istep", sc, i), ("reset", sc, b), ("vnmult", (2,), f), ("nl_istep_last", sc, i),
                  ("nl_initializing", sc, b), ("exb_first", sc, b)]
    tm += timeblock
    ri = [(n, sc, f) for n in ("delt_adj", "delt_cushion", "delt_minimum", "dt0")]
    ri += [("abort_rapid_time_step_change", sc, b), ("in_memory", sc, b)]
    ig = [("tstart", sc, f), ("new_field_init", sc, b), ("scale", sc, f), ("ginitopt_switch", sc, i)]
    fl = [("fieldopt_switch", sc, i), ("remove_zonal_flows_switch", sc, b), ("force_maxwell_reinit", sc, b)]
    nl = [(n, sc, b) for n in ("nonlin", "split_nonlinear", "accelerated", "use_cfl_limit", "use_2d_cfl",
                               "use_order_based_error", "include_phi", "include_apar", "include_bpar",
                               "zip", "nl_forbid_force_zero")]
    nl += [("cfl", sc, f), ("error_target", sc, f), ("istep_error_start", sc, i), ("cflx", sc, f), ("cfly", sc, f)]
    lay = [("iglo_of", G5, i), ("ile_of", (S, Y, X, T), i)]
    lay += [(f"g_lo_{x}_ord", sc, i) for x in ("ik", "it", "il", "ie", "is")]
    lay += [(f"le_lo_{x}_ord", sc, i) for x in ("ig", "it", "ik", "is")]
    dfa = [("vpa", G7, f), ("vpac", G7, f), ("vperp2", G6, f), ("aj0", G6, f), ("aj1", G6, f),
           ("wstar", (S, E, Y), f)]
    df = [("ufac", (S, E), f), ("fexp", (S,), f), ("bkdiff", (S,), f)]
    df += [(n, sc, f) for n in ("driftknob", "tpdriftknob", "mach", "g_exb", "g_exbfac", "g_exb_start_time",
                                "omprimfac", "phi_ext", "t0", "omega0", "gamma0", "source0", "wfb", "gridfac",
                                "apfac", "poisfac", "afilter")]
    df += [(n, sc, i) for n in ("g_exb_start_timestep", "source_option_switch", "boundary_option_switch",
                                "adiabatic_option_switch", "i_class")]
    df += [(n, sc, b) for n in ("esv", "nonad_zero", "zero_forbid", "def_parity", "even", "hyper_in_initialisation",
                                "start_from_previous_solution", "boundary_off_grid", "exponential_boundary",
                                "opt_source", "mult_imp", "gf_lo_integrate", "no_connections", "exb_first")]
    df += [(n, F3, f) for n in ("gamtot", "gamtot1", "gamtot2", "gamtot3", "gridfac1",
                                "inv_phi_denominator_g", "inv_bpar_denominator_g")]
    df += [("awgt", K2, f)]
    df += [(n, G7, f) for n in ("vpar", "wdrift", "wdriftttp")]
    df += [(n, G7c, f) for n in ("a", "b", "r", "ainv", "g_h")]
    df += [(n, (S, E, L, Y, X, 2, T - 1, 2), f) for n in ("source_coeffs_phim", "source_coeffs_phip",
                                                          "source_coeffs_aparm", "source_coeffs_aparp")]
    df += [(n, K2, i) for n in ("itleft", "itright", "l_links", "r_links", "itleftmost", "itrightmost",
                                "supercell_labels", "n_cells")]
    df += [("n_links", (Y, X, 2), i), ("n_supercells", (Y,), i), ("connections_neighbor", G5, b),
           ("save_h", G5 + (2,), b), ("M_class", (C,), i), ("N_class", (C,), i)]
    df += [(n, (S, Y, X, E), f) for n in ("wfb_aj0l", "wfb_aj1l", "wfb_aj0r", "wfb_aj1r")]
    df += [(n, (S, E), f) for n in ("wfb_vperp2left", "wfb_vperp2right")]
    hy = [("D_res", K2, f), ("hypervisc_filter", F3, f)]
    hy += [(n, sc, f) for n in ("akx4_max", "aky_max", "aky4_max", "akperp4_max", "kperp2_max",
                                "kperp2_max_zonal", "D_hypervisc", "D_hyperres", "omega_osc", "D_hyper3D",
                                "P_hyper3D", "ky_cut", "kx_cut")]
    hy += [(n, sc, b) for n in ("hyper_on", "const_amp", "isotropic_shear", "damp_zonal_only", "gridnorm",
                                "isotropic_model")]
    hy += [("nexp", sc, i), ("hyper_option_switch", sc, i)]
    st = timeblock + [("g", G7c, f), ("gnew", G7c, f)] + [(n, F3c, f) for n in (
        "phi", "apar", "bpar", "phinew", "aparnew", "bparnew")]
    st += [(n, G7c, f) for n in ("gexp_1", "gexp_2", "gexp_3")] + [("kx_shift", (Y,), f), ("theta0_shift", (Y,), f)]
    out = {
        "theta_grid": tg, "kt_grids": kt, "le_grids": le, "species_": sp, "run_parameters": rp,
        "gs2_time": tm, "gs2_reinit": ri, "init_g": ig, "fields": fl, "nonlinear_terms": nl,
        "gs2_layouts": lay, "dist_fn_arrays": dfa, "dist_fn": df, "hyper": hy,
        "state_t0_prediag": st, "state_t0": st,
    }
    if ex.fieldopt == 1:
        fi = [("nfield", sc, i), ("nidx", sc, i), ("cls_ik", (C, d["Nm"], d["Mm"]), i),
              ("cls_it", (C, d["Nm"], d["Mm"]), i), ("sc_class", K2, i), ("sc_member", K2, i), ("sc_cell", K2, i)]
        for c in range(1, C + 1):
            fi.append((f"aminv_c{c}", (d[f"M_c{c}"], d[f"R_c{c}"], d[f"R_c{c}"], 2), f))
        out["fields_implicit"] = fi
    # always present, also when the collision operator is inactive (split_collisions still re-solves the
    # fields in collisions_advance)
    out["collisions"] = [("collision_model_switch", sc, i), ("colls", sc, b), ("split_collisions", sc, b),
                         ("vnmult", (2,), f)]
    if ex.colls_active:
        co = [(n, (S, E, Y), f) for n in ("vnew", "vnew_s", "vnew_D", "vnew_E", "delvnew")]
        co += [("vnewh", (S, Y, X, T), f), ("pitch_weights", (T, L), f), ("vpdiff", (L, 2, T), f),
               ("vpdiffle", (T, XI), f), ("e15", (E,), f), ("capg_half", (Em1,), f), ("vnmult", (2,), f)]
        co += [(n, (S, Y, X, T, E, XI), f) for n in ("c1le", "betaale", "qle", "ec1le", "ebetaale", "eqle")]
        co += [(n, (S, Y, X, T, E1, XI), f) for n in ("z0le", "s0le", "w0le", "bz0le", "bs0le", "bw0le",
                                                      "aj0le", "vperp_aj1le", "vpa_aj0_le")]
        co += [(n, sc, b) for n in ("colls", "has_lorentz", "has_diffuse", "drag", "use_le_layout",
                                    "conserve_moments", "conservative", "resistivity", "const_v", "ei_coll_only",
                                    "hypermult", "hyper_colls", "special_wfb_lorentz", "vpar_zero_mean",
                                    "conserve_forbid_zero", "force_collisions", "adjust", "heating",
                                    "split_collisions", "vary_vnew")]
        co += [(n, sc, i) for n in ("collision_model_switch", "lorentz_switch", "ediff_switch", "nxi_lim")]
        co += [("cfac", sc, f)]
        out["collisions"] = co
    if ex.diag:
        dg = [(n, sc, i) for n in ("nwrite", "navg", "igomega", "nout")]
        dg += [(n, sc, f) for n in ("omegatol", "omegatinst")]
        dg += [(n, sc, b) for n in ("exit_when_converged", "calculate_fluxes", "write_omega", "write_fields",
                                    "write_fluxes", "write_moments", "write_correlation", "write_eigenfunc",
                                    "write_line", "write_ascii", "write_phi_over_time", "write_ql_metric")]
        out["diagnostics"] = dg
    return out


class Export:
    def __init__(self, path):
        self.path = path
        self.ds = netCDF4.Dataset(path)
        ds = self.ds
        for _, g in walk(ds):
            g.set_auto_mask(False)
        D = {k: len(v) for k, v in ds.dimensions.items()}
        self.D = D
        self.d = {"S": D["species"], "E": D["egrid"], "L": D["lambda"], "Y": D["ky"], "X": D["kx"],
                  "T": D["theta"], "E1": D["egrid_p1"], "Em1": D["egrid_m1"], "L2": D["xi2"],
                  "XI": D.get("xi_lim", -1), "XP": D["xi_p1"], "C": D.get("n_class", 0),
                  "Mm": D.get("sc_Mmax", 0), "Nm": D.get("sc_Nmax", 0)}
        for c in range(1, self.d["C"] + 1):
            if f"sc_R_c{c}" in D:
                self.d[f"R_c{c}"] = D[f"sc_R_c{c}"]
                self.d[f"M_c{c}"] = D[f"sc_M_c{c}"]
        self.fieldopt = int(ds["/fields/fieldopt_switch"][...])
        co = ds["/collisions"]
        self.colls_active = bool(int(co["colls"][...])) and int(co["collision_model_switch"][...]) != 3 \
            and "vnew" in co.variables
        self.diag = int(ds.getncattr("diagnostics_run_at_t0")) == 1
        self.reference_build = "-ffp-contract=off" in str(ds.getncattr("build_config"))
        self.ntgrid = int(ds["/theta_grid/ntgrid"][...])

    def var(self, path):
        return np.asarray(self.ds[path][...])

    def has(self, path):
        grp, name = path.rsplit("/", 1)
        try:
            return name in self.ds[grp or "/"].variables
        except (KeyError, IndexError):
            return False


# ----------------------------------------------------------------------------------- checks

def check_schema(ex, rep):
    print("[schema] expected variables (10 §A), shapes, dtypes, absent lists")
    ds = ex.ds
    rep.ok(str(ds.getncattr("export_format_version")).split(".")[0] == "1", "export_format_version major == 1",
           str(ds.getncattr("export_format_version")))
    rep.ok(int(ds.getncattr("real_kind_bytes")) == 8, "real_kind_bytes == 8")
    for att in ("git_hash", "build_config", "layout", "index_convention", "gate_status", "used_inputs",
                "input_file_text", "run_name"):
        rep.ok(att in ds.ncattrs(), f"global attribute {att} present")
    sch = schema(ex.d, ex)
    nmiss = 0
    nabs = 0
    nchecked = 0
    for gname, entries in sch.items():
        if gname not in ds.groups:
            rep.ok(False, f"group /{gname} present")
            continue
        grp = ds[gname]
        absent = absent_set(grp)
        for name, shape, dt in entries:
            if name in grp.variables:
                v = grp.variables[name]
                nchecked += 1
                good = tuple(v.shape) == tuple(shape) and v.dtype == np.dtype(dt)
                if not good:
                    rep.ok(False, f"/{gname}/{name} shape/dtype", f"{v.shape}/{v.dtype} expected {shape}/{dt}")
                    nmiss += 1
                ftype = v.getncattr("fortran_type") if "fortran_type" in v.ncattrs() else None
                if ftype is None or "canonical" not in v.ncattrs():
                    rep.ok(False, f"/{gname}/{name} has fortran_type and canonical attributes")
                    nmiss += 1
                elif dt == "i1" and ftype != "logical":
                    rep.ok(False, f"/{gname}/{name} int8 must be fortran_type=logical", ftype)
                    nmiss += 1
            elif name in absent or any(a.endswith("*") and name.startswith(a[:-1]) for a in absent):
                nabs += 1
            else:
                rep.ok(False, f"/{gname}/{name} present or listed in absent")
                nmiss += 1
    rep.ok(nmiss == 0, f"schema: {nchecked} variables with correct shape/dtype/attributes, "
                       f"{nabs} listed as absent, {nmiss} problems")
    # complex variables: trailing ri axis of length 2
    bad = []
    for gpath, grp in walk(ds):
        for name, v in grp.variables.items():
            if "fortran_type" in v.ncattrs() and v.getncattr("fortran_type") == "complex":
                if v.dimensions[-1] != "ri":
                    bad.append(f"{gpath}/{name}")
    rep.ok(not bad, "complex variables have a trailing ri axis", ", ".join(bad[:5]))
    # gates
    gates = ds["/gates"]
    viol = [n for n, v in gates.variables.items() if int(v[...]) != 0 and not n.startswith("hard_")]
    hard = [n for n, v in gates.variables.items() if int(v[...]) != 0 and n.startswith("hard_")]
    status = str(ds.getncattr("gate_status"))
    rep.ok(not hard, "no hard gate violated", ",".join(hard))
    rep.ok(status == ("ok" if not viol else ";".join(viol)), "gate_status consistent with /gates", status)
    if viol:
        rep.note(f"soft gates violated (JAX loader will refuse unless overridden): {status}")


def check_numerics(ex, rep, outnc=None, ref_outnc=None):
    print("[numerics] recomputation from exported inputs (GS2 expressions)")
    v = ex.var
    exact_fma = ex.reference_build  # FMA-prone expressions must be bitwise for reference builds

    def fma_cmp(name, mine, theirs):
        if exact_fma:
            return rep.ok(bits_equal(mine, theirs), f"{name} (bitwise; reference build)", describe_diff(mine, theirs))
        # stock build: FMA contraction changes the rounding (several ulp where terms cancel)
        dmax = float(np.max(np.abs(np.asarray(mine) - np.asarray(theirs))))
        scale = float(np.max(np.abs(theirs))) or 1.0
        return rep.ok(dmax <= 1e-14 * scale, f"{name} (stock build: <= 1e-14 relative)",
                      f"max|diff|/max|ref|={dmax / scale:.2e}, ulp={max_ulp(mine, theirs)}")

    ntg = ex.ntgrid
    bmag = v("/theta_grid/bmag")
    al = v("/le_grids/al")
    energy = v("/le_grids/energy")  # (S,E)
    forbid = v("/le_grids/forbid").astype(bool)  # (L,T)
    S, E, L, Y, X, T = (ex.d[k] for k in ("S", "E", "L", "Y", "X", "T"))

    # 1. vperp2 = bmag*al*e (dist_fn.fpp:1584), 0 where forbidden
    vp2 = (bmag[None, None, :] * al[None, :, None])[None, :, :, :] * energy[:, :, None, None]  # (S,E,L,T)
    vp2 = np.where(forbid[None, None, :, :], 0.0, vp2)
    mine = np.broadcast_to(vp2[:, :, :, None, None, :], (S, E, L, Y, X, T))
    rep.ok(bits_equal(np.ascontiguousarray(mine), v("/dist_fn_arrays/vperp2")), "vperp2 == (bmag*al)*energy (bitwise)")

    # 2. vpa = sqrt(e*max(0, 1-al*bmag))*vparknob (dist_fn.fpp:1641), -vpa for sign 2, 0 where forbidden
    vparknob = float(v("/dist_fn/vparknob")) if ex.has("/dist_fn/vparknob") else 1.0
    one_m = 1.0 - al[:, None] * bmag[None, :]  # (L,T)
    vpa1 = np.sqrt(energy[:, :, None, None] * np.maximum(0.0, one_m)[None, None]) * vparknob  # (S,E,L,T)
    vpa1 = np.where(forbid[None, None], 0.0, vpa1)
    vpa_mine = np.zeros((S, E, L, Y, X, 2, T))
    vpa_mine[:, :, :, :, :, 0, :] = vpa1[:, :, :, None, None, :]
    vpa_mine[:, :, :, :, :, 1, :] = -vpa1[:, :, :, None, None, :]
    vpa_mine = np.where(forbid[None, None, :, None, None, None, :], 0.0, vpa_mine)
    vpa = v("/dist_fn_arrays/vpa")
    if exact_fma:
        fma_cmp("vpa == sqrt(e*max(0,1-al*bmag))*vparknob", vpa_mine, vpa)
    else:
        # A stock build contracts 1.0 - al*bmag into an FMA: 1 ulp in general and sqrt(e*O(eps))
        # instead of 0 at bounce points (this is what the reference build removes).
        dmax = float(np.max(np.abs(vpa_mine - vpa)))
        rep.ok(dmax <= 1e-6 * float(np.max(np.abs(vpa))), "vpa == sqrt(e*max(0,1-al*bmag))*vparknob "
               "(stock build: FMA-level differences allowed)", f"max|diff|={dmax:.3e}, ulp={max_ulp(vpa_mine, vpa)}")

    # 3. vpac from the exported vpa (dist_fn.fpp:1666-1675)
    fb = forbid[:, :-1] | forbid[:, 1:]  # (L,T-1)
    vpac_mine = np.zeros_like(vpa)
    centred = 0.5 * (vpa[..., :-1] + vpa[..., 1:])
    vpac_mine[..., :-1] = np.where(fb[None, None, :, None, None, None, :], 0.0, centred)
    rep.ok(bits_equal(vpac_mine, v("/dist_fn_arrays/vpac")), "vpac == 0.5*(vpa(ig)+vpa(ig+1)) (bitwise)")

    # 4. vpar (dist_fn.fpp:1715-1725) from vpac, zstm, tunits, code_dt, delthet, gradpar
    def vpar_of(code_dt, tunits):
        zstm = v("/species_/spec_zstm")
        delthet = v("/theta_grid/delthet")
        gp = v("/theta_grid/gradpar")
        vpac = v("/dist_fn_arrays/vpac")
        out = np.zeros_like(vpac)
        fac = ((zstm[:, None] * tunits[None, :]) * code_dt) * 0.5  # (S,Y)
        gsum = np.abs(gp[:-1]) + np.abs(gp[1:])
        out[..., :-1] = ((fac[:, None, None, :, None, None, None] / delthet[:-1]) * gsum) * vpac[..., :-1]
        return out
    code_dt = float(v("/gs2_time/code_dt"))
    tunits = v("/gs2_time/tunits")
    rep.ok(bits_equal(vpar_of(code_dt, tunits), v("/dist_fn/vpar")), "vpar recomputed (02 §3.1 association) (bitwise)")

    # 5. wstar = -code_dt*wunits(ik)*dlnf0drho(ie,is) (dist_fn.fpp:1746)
    def wstar_of(code_dt, wunits):
        dl = v("/species_/dlnf0drho")  # (S,E)
        return -((code_dt * wunits[None, None, :]) * dl[:, :, None])
    wunits = v("/gs2_time/wunits")
    rep.ok(bits_equal(wstar_of(code_dt, wunits), v("/dist_fn_arrays/wstar")), "wstar == -(code_dt*wunits)*dlnf0drho (bitwise)")

    # 6. tunits / wunits / woutunits (gs2_time.f90:107-125)
    aky = v("/kt_grids/aky")
    if int(v("/run_parameters/wstar_units")) == 0:
        rep.ok(bits_equal(tunits, np.ones(Y)) and bits_equal(wunits, aky / 2.0)
               and bits_equal(v("/gs2_time/woutunits"), 1.0 / tunits), "tunits=1, wunits=aky/2, woutunits=1/tunits (bitwise)")

    # 7. has_* = f* > epsilon (run_parameters.f90:351-353)
    eps = np.finfo(np.float64).eps
    ok = all(int(v(f"/run_parameters/has_{x}")) == int(float(v(f"/run_parameters/f{x}")) > eps)
             for x in ("phi", "apar", "bpar"))
    rep.ok(ok, "has_phi/apar/bpar == fphi/fapar/fbpar > epsilon(0.0)")

    # 8. field_line_average_weight = 1/trapz(theta, jacob) with the Klein sum (theta_grid.f90:520)
    w = float(v("/theta_grid/field_line_average_weight"))
    wpy = 1.0 / trapz_gs2(v("/theta_grid/theta"), v("/theta_grid/jacob"))
    rep.ok(np.float64(w).tobytes() == np.float64(wpy).tobytes(), "field_line_average_weight == 1/trapz_klein(theta, jacob) (bitwise)",
           f"{w!r} vs {wpy!r}")
    if ex.has("/theta_grid/field_line_average_weight_check"):
        rep.ok(np.float64(w).tobytes() == np.float64(v("/theta_grid/field_line_average_weight_check")).tobytes(),
               "field_line_average_weight == exporter recompute (bitwise)")

    # 9. kperp2 (kt_grids.f90:313-333)
    gds2, gds21, gds22 = v("/theta_grid/gds2"), v("/theta_grid/gds21"), v("/theta_grid/gds22")
    shat = float(v("/theta_grid/shat"))
    akx = v("/kt_grids/akx")
    th0 = v("/kt_grids/theta0")  # (Y,X)
    kp = np.zeros((Y, X, T))
    for ik in range(Y):
        for it in range(X):
            if aky[ik] == 0.0:
                kp[ik, it] = ((akx[it] * akx[it]) * gds22) / (shat * shat)
            else:
                kp[ik, it] = (aky[ik] * aky[ik]) * ((gds2 + (2.0 * th0[ik, it]) * gds21) + (th0[ik, it] * th0[ik, it]) * gds22)
    fma_cmp("kperp2 recomputed from gds2/gds21/gds22/theta0", kp, v("/kt_grids/kperp2"))
    ikp = np.where(np.abs(v("/kt_grids/kperp2")) > eps, 1.0 / np.where(v("/kt_grids/kperp2") == 0, 1, v("/kt_grids/kperp2")), 0.0)
    rep.ok(bits_equal(ikp, v("/kt_grids/inv_kperp2")), "inv_kperp2 == 1/kperp2 where |kperp2| > eps (bitwise)")

    # 10. layout maps (gs2_layouts.fpp:917-948, 1606-1617, 2160-2170)
    layout = str(ex.ds.getncattr("layout")).strip()
    ords = LAYOUT_ORDS[layout]
    ok = tuple(int(v(f"/gs2_layouts/g_lo_{x}_ord")) for x in ("it", "ik", "il", "ie", "is")) == ords
    rep.ok(ok, f"g_lo orders match layout '{layout}'")
    sizes = {"it": X, "ik": Y, "il": L, "ie": E, "is": S}
    names = ("it", "ik", "il", "ie", "is")
    dim_size = [0] * 5
    for nme, o in zip(names, ords):
        dim_size[o - 1] = sizes[nme]
    idx = np.zeros((S, E, L, Y, X), dtype=np.int64)
    grids = np.meshgrid(np.arange(X), np.arange(Y), np.arange(L), np.arange(E), np.arange(S), indexing="ij")
    i5 = [None] * 5
    for nme, gidx, o in zip(names, grids, ords):
        i5[o - 1] = gidx
    flat = i5[0] + dim_size[0] * (i5[1] + dim_size[1] * (i5[2] + dim_size[2] * (i5[3] + dim_size[3] * i5[4])))
    idx[...] = np.transpose(flat, (4, 3, 2, 1, 0))
    rep.ok(np.array_equal(idx, v("/gs2_layouts/iglo_of")), "iglo_of == idx(g_lo, ik,it,il,ie,is) recomputed from layout")
    rep.ok(np.array_equal(np.sort(v("/gs2_layouts/iglo_of").ravel()), np.arange(S * E * L * Y * X)),
           "iglo_of is a permutation of 0..N_g-1")
    if not ex.has("/gs2_layouts/ile_of"):
        rep.note("le_lo not initialised (no collisions): ile_of and le_lo orders absent")
        return check_numerics_rest(ex, rep, outnc, ref_outnc, locals())
    itik = (2, 3) if layout in ("lexys", "lxyes", "xyles") else (3, 2)
    lords = {"ig": 1, "it": itik[0], "ik": itik[1], "is": 4}
    ok = all(int(v(f"/gs2_layouts/le_lo_{x}_ord")) == lords[x] for x in lords)
    rep.ok(ok, f"le_lo orders match layout '{layout}'")
    lsizes = {"ig": T, "it": X, "ik": Y, "is": S}
    lds = [0] * 4
    for nme, o in lords.items():
        lds[o - 1] = lsizes[nme]
    g4 = np.meshgrid(np.arange(T), np.arange(X), np.arange(Y), np.arange(S), indexing="ij")
    i4 = [None] * 4
    for nme, gidx in zip(("ig", "it", "ik", "is"), g4):
        i4[lords[nme] - 1] = gidx
    lflat = i4[0] + lds[0] * (i4[1] + lds[1] * (i4[2] + lds[2] * i4[3]))
    rep.ok(np.array_equal(np.transpose(lflat, (3, 2, 1, 0)), v("/gs2_layouts/ile_of")),
           "ile_of == idx(le_lo, ig,ik,it,is) recomputed from layout")
    return check_numerics_rest(ex, rep, outnc, ref_outnc, locals())


def check_numerics_rest(ex, rep, outnc, ref_outnc, env):
    v = ex.var
    ntg = ex.ntgrid
    S, E, L, Y, X, T = (ex.d[k] for k in ("S", "E", "L", "Y", "X", "T"))
    aky, akx = env["aky"], env["akx"]
    code_dt = env["code_dt"]
    vpar_of, wstar_of = env["vpar_of"], env["wstar_of"]
    forbid, al = env["forbid"], env["al"]

    # 11. kwork_filter: box -> only (1,1)
    kwf = v("/kt_grids/kwork_filter").astype(bool)
    expect = np.zeros((Y, X), dtype=bool)
    if int(v("/kt_grids/gridopt_switch")) == 4:
        expect[0, 0] = True
    rep.ok(np.array_equal(kwf, expect), "kwork_filter after init_fields: (1,1) only for box grids (G-01)")

    # 12. aminv: rows of aky = akx = 0 members are zero (fields_implicit.f90:855-857)
    if ex.fieldopt == 1 and ex.d["C"] > 0:
        cls_ik, cls_it = v("/fields_implicit/cls_ik"), v("/fields_implicit/cls_it")
        nidx = int(v("/fields_implicit/nidx"))
        ok = True
        for c in range(1, ex.d["C"] + 1):
            p = f"/fields_implicit/aminv_c{c}"
            if not ex.has(p):
                continue
            a = v(p)
            M = a.shape[0]
            for m in range(M):
                ik, it = cls_ik[c - 1, 0, m], cls_it[c - 1, 0, m]
                N = a.shape[1] // nidx
                for n in range(N):
                    ikn, itn = cls_ik[c - 1, n, m], cls_it[c - 1, n, m]
                    if aky[ikn - 1] == 0.0 and akx[itn - 1] == 0.0:
                        rows = a[m, n * nidx:(n + 1) * nidx, :, :]
                        ok = ok and not np.any(rows)
        rep.ok(ok, "aminv rows of the aky=akx=0 member are exactly zero")
        tot = sum(np.count_nonzero(cls_it[c - 1] > 0) for c in range(1, ex.d["C"] + 1))
        rep.ok(tot == Y * X, "class tables cover every (it,ik) exactly once", f"{tot} vs {Y * X}")

    # 12b. LE6 canonicaliser: init_le_bessel (collisions.fpp:634-687) gathers aj0, aj1, vpa into the
    # le layout with g2le, a pure copy that skips forbidden points (06 §2.3)
    if ex.colls_active and ex.has("/collisions/aj0le"):
        ixi_il = v("/le_grids/ixi_to_il")      # (T, L2) 1-based
        ixi_sg = v("/le_grids/ixi_to_isgn")    # (T, L2)
        aj0, aj1 = v("/dist_fn_arrays/aj0"), v("/dist_fn_arrays/aj1")  # (S,E,L,Y,X,T)
        vpa_ = v("/dist_fn_arrays/vpa")         # (S,E,L,Y,X,2,T)
        emax = v("/le_grids/energy_maxw")
        aj0le = v("/collisions/aj0le")          # (S,Y,X,T,E+1,XI)
        XI = aj0le.shape[-1]
        m_aj0 = np.zeros_like(aj0le)
        m_aj1 = np.zeros_like(aj0le)
        m_vpa = np.zeros_like(aj0le)
        for t in range(T):
            for ixi in range(XI):
                il = int(ixi_il[t, ixi]) - 1
                isg = int(ixi_sg[t, ixi]) - 1
                if il < 0 or forbid[il, t]:
                    continue
                # aj0[s,e,il,y,x,t] -> [s,y,x,t,e,ixi]
                m_aj0[:, :, :, t, :E, ixi] = np.transpose(aj0[:, :, il, :, :, t], (0, 2, 3, 1))
                m_aj1[:, :, :, t, :E, ixi] = (np.transpose(aj1[:, :, il, :, :, t], (0, 2, 3, 1))
                                              * emax[None, None, None, :]) * al[il]
                m_vpa[:, :, :, t, :E, ixi] = np.transpose(vpa_[:, :, il, :, :, isg, t], (0, 2, 3, 1))
        rep.ok(bits_equal(m_aj0, aj0le), "aj0le == g2le gather of aj0 (LE6 canonical map, bitwise)")
        rep.ok(bits_equal(m_aj1, v("/collisions/vperp_aj1le")), "vperp_aj1le == (aj1*energy_maxw)*al gathered (bitwise)")
        rep.ok(bits_equal(m_vpa * m_aj0, v("/collisions/vpa_aj0_le")), "vpa_aj0_le == vpa*aj0le gathered (bitwise)")
        if ex.has("/collisions/e15"):
            rep.ok(bits_equal(emax ** 1.5, v("/collisions/e15")), "e15 == energy_maxw**1.5 (bitwise)")

    # 13. dt ladder rungs
    if "dt_ladder" in ex.ds.groups:
        lad = ex.ds["dt_ladder"]
        nr = int(lad.getncattr("n_rungs")) if "n_rungs" in lad.ncattrs() else 0
        dtk = code_dt
        delt_adj = float(v("/gs2_reinit/delt_adj"))
        for k in range(1, nr + 1):
            base = f"/dt_ladder/k{k}"
            dtk = dtk / delt_adj
            rk = float(v(f"{base}/code_dt"))
            rep.ok(np.float64(rk).tobytes() == np.float64(dtk).tobytes(), f"ladder k{k}: code_dt == dt_init/delt_adj**{k}",
                   f"{rk!r}")
            rep.ok(bits_equal(vpar_of(rk, v(f"{base}/gs2_time/tunits")), v(f"{base}/dist_fn/vpar")),
                   f"ladder k{k}: vpar recomputed with the rung code_dt (bitwise)")
            rep.ok(bits_equal(wstar_of(rk, v(f"{base}/gs2_time/wunits")), v(f"{base}/dist_fn_arrays/wstar")),
                   f"ladder k{k}: wstar recomputed with the rung code_dt (bitwise)")
            if ex.has(f"{base}/hyper/D_res") and ex.has("/hyper/D_res"):
                rep.ok(bits_equal(v(f"{base}/hyper/D_res"), v("/hyper/D_res")), f"ladder k{k}: D_res frozen at the first-init dt")
            if ex.fieldopt == 1 and ex.has(f"{base}/fields_implicit/aminv_c1"):
                same = bits_equal(v(f"{base}/fields_implicit/aminv_c1"), v("/fields_implicit/aminv_c1"))
                rep.ok(not same, f"ladder k{k}: aminv_c1 rebuilt (differs from the initial-dt matrix)")
        if nr == 0:
            rep.note("dt ladder group present with 0 rungs")

    # 14. state_t0 vs trace-independent invariants
    st = "/state_t0"
    ok = bits_equal(v(f"{st}/phinew"), v(f"{st}/phi")) and bits_equal(v(f"{st}/aparnew"), v(f"{st}/apar")) \
        and bits_equal(v(f"{st}/bparnew"), v(f"{st}/bpar"))
    rep.ok(ok, "state_t0: phi == phinew, apar == aparnew, bpar == bparnew (bitwise)")
    for x in ("phinew", "aparnew", "bparnew"):
        if not int(v(f"/run_parameters/has_{x[:-3]}")):
            rep.ok(not np.any(v(f"{st}/{x}")), f"state_t0: absent field {x} is exactly zero")
    rep.ok(bits_equal(v(f"{st}/phinew"), v("/state_t0_prediag/phinew")), "istep-0 diagnostics do not change phinew")

    # 15. GS2's own out.nc
    for label, path in (("export-run out.nc", outnc), ("bin/gs2 out.nc", ref_outnc)):
        if not path or not os.path.exists(path):
            continue
        o = netCDF4.Dataset(path)
        o.set_auto_mask(False)
        pairs = [("theta", "/theta_grid/theta"), ("bmag", "/theta_grid/bmag"), ("gradpar", "/theta_grid/gradpar"),
                 ("gbdrift", "/theta_grid/gbdrift"), ("gbdrift0", "/theta_grid/gbdrift0"),
                 ("cvdrift", "/theta_grid/cvdrift"), ("cvdrift0", "/theta_grid/cvdrift0"),
                 ("gds2", "/theta_grid/gds2"), ("gds21", "/theta_grid/gds21"), ("gds22", "/theta_grid/gds22"),
                 ("grho", "/theta_grid/grho"), ("jacob", "/theta_grid/jacob"), ("cdrift", "/theta_grid/cdrift"),
                 ("bpol", "/theta_grid/Bpol"), ("shat", "/theta_grid/shat"), ("drhodpsi", "/theta_grid/drhodpsi"),
                 ("qval", "/theta_grid/qval"), ("kx", "/kt_grids/akx"), ("ky", "/kt_grids/aky"),
                 ("theta0", "/kt_grids/theta0"), ("kperp2", "/kt_grids/kperp2"), ("lambda", "/le_grids/al"),
                 ("energy", "/le_grids/energy")]
        bad = []
        n = 0
        for on, en in pairs:
            if on in o.variables and ex.has(en):
                n += 1
                if not bits_equal(np.asarray(o[on][...], dtype=np.float64), v(en)):
                    bad.append(on)
        rep.ok(not bad and n >= 15, f"{label}: {n} geometry/grid variables bitwise equal to the export", ",".join(bad))
        if "phi_igomega_by_mode" in o.variables and int(v("/run_parameters/has_phi")):
            igomega = int(v("/diagnostics/igomega")) if ex.has("/diagnostics/igomega") else 0
            ph0 = np.asarray(o["phi_igomega_by_mode"][0])
            mine = v(f"{st}/phinew")[:, :, ntg + igomega, :]
            rep.ok(bits_equal(ph0, mine), f"{label}: phi_igomega_by_mode(t=0) == state_t0 phinew(igomega) (bitwise)")
        finalised_at_t0 = "t" in o.dimensions and len(o.dimensions["t"]) == 1
        if label.startswith("export") and finalised_at_t0 and "phi" in o.variables \
                and int(v("/run_parameters/has_phi")):
            rep.ok(bits_equal(np.asarray(o["phi"][...]), v(f"{st}/phinew")),
                   f"{label}: final phi (finalised at t=0) == state_t0 phinew (bitwise)")
        o.close()


# ----------------------------------------------------------------------------------- trace

def check_trace(ex, path, rep):
    print(f"[trace] invariants of {os.path.basename(path)} (09 §6.2)")
    tr = netCDF4.Dataset(path)
    for _, g in walk(tr):
        g.set_auto_mask(False)
    steps = sorted(k for k in tr.groups if k.startswith("s"))
    rep.ok(len(steps) > 0, "trace has step groups", f"{len(steps)} steps, {tr.getncattr('n_records')} records")
    recs = {}  # (istep, attempt) -> list of (seq, tag, group)
    for s in steps:
        istep = int(s[1:])
        for a, ga in tr[s].groups.items():
            att = int(a[1:])
            lst = []
            for pname, gp in ga.groups.items():
                lst.append((int(gp.getncattr("seq")), str(gp.getncattr("tag")), gp))
            recs[(istep, att)] = sorted(lst, key=lambda x: x[0])

    def arr(g, n):
        return np.asarray(g.variables[n][...])

    def find(istep, attempt, tag, which=0):
        lst = [r for r in recs.get((istep, attempt), []) if r[1] == tag]
        if not lst:
            return None
        return lst[which][2] if which >= 0 else lst[-1][2]

    state_vars = ("gnew", "phinew", "aparnew", "bparnew", "gexp_1", "gexp_2", "gexp_3", "kx_shift", "theta0_shift")
    tb = ("code_dt", "user_dt", "code_dt_prev1", "code_dt_prev2", "code_dt_old", "code_time", "user_time",
          "code_dt_cfl", "code_dt_min", "code_dt_max", "reset", "vnmult", "nl_istep_last", "exb_first")

    # T0: advance_begin of the first step == /state_t0 (the JAX initial state)
    first = min(k[0] for k in recs) if recs else None
    if first is not None:
        ab = find(first, 1, "advance_begin")
        if ab is not None:
            ok = True
            bad = []
            for n in state_vars + tb:
                if n in ab.variables and ex.has(f"/state_t0/{n}"):
                    if not bits_equal(arr(ab, n), ex.var(f"/state_t0/{n}")):
                        ok = False
                        bad.append(n)
            rep.ok(ok, f"advance_begin(step {first}) == /state_t0 (bitwise)", ",".join(bad))

    immediate_reset = bool(int(ex.var("/run_parameters/immediate_reset")))
    # g2le index map over the allowed points (le_grids.f90:2972-3011, 3379-3599; 06 §2.3)
    gmap = None
    if ex.colls_active and int(ex.var("/collisions/use_le_layout")):
        ixi_il = ex.var("/le_grids/ixi_to_il")
        ixi_sg = ex.var("/le_grids/ixi_to_isgn")
        forbid = ex.var("/le_grids/forbid").astype(bool)
        XP = ex.d["XP"]
        gmap = [(t, ixi, int(ixi_il[t, ixi]) - 1, int(ixi_sg[t, ixi]) - 1)
                for t in range(ex.d["T"]) for ixi in range(min(XP, ixi_il.shape[1]))
                if int(ixi_il[t, ixi]) >= 1 and not forbid[int(ixi_il[t, ixi]) - 1, t]]
        E_ = ex.d["E"]

        def gather(h):  # h (S,E,L,Y,X,2,T) -> gle (S,Y,X,T,E+1,nxi+1)
            S_, _, _, Y_, X_, _, T_ = h.shape
            out = np.zeros((S_, Y_, X_, T_, E_ + 1, XP), dtype=h.dtype)
            for t, ixi, il, isg in gmap:
                out[:, :, :, t, :E_, ixi] = np.transpose(h[:, :, il, :, :, isg, t], (0, 2, 3, 1))
            return out

        def scatter(gle, h_old):
            out = h_old.copy()
            for t, ixi, il, isg in gmap:
                out[:, :, il, :, :, isg, t] = np.transpose(gle[:, :, :, t, :E_, ixi], (0, 3, 1, 2))
            return out
    attempts = {}
    for (istep, att) in recs:
        attempts[istep] = max(attempts.get(istep, 0), att)
    istep_list = sorted(attempts)
    nfail = {}

    def tally(name, cond):
        nfail.setdefault(name, [0, 0])
        nfail[name][0] += 1
        if not cond:
            nfail[name][1] += 1

    for k in istep_list:
        last = attempts[k]
        se = find(k, last, "step_end")
        # (a) step_end(k) == advance_begin(k+1, attempt 1)
        if se is not None and (k + 1) in attempts:
            ab = find(k + 1, 1, "advance_begin")
            se_seq = [x[0] for x in recs[(k, last)] if x[1] == "step_end"][0]
            reset_after = [r for r in recs[(k, last)] if r[1] == "after_reset_time_step" and r[0] > se_seq]
            if ab is not None and reset_after:
                # deferred reset (gs2_main.fpp:800-816): the state is restored at init level 34
                # (restart-file round trip unless in_memory: values identical, signed zeros may
                # change) and the fields are recomputed by get_init_field (round-off)
                rr = reset_after[-1][2]
                tally("deferred reset: after_reset_time_step time block == next advance_begin (except reset)",
                      all(bits_equal(arr(rr, n), arr(ab, n)) for n in tb if n in rr.variables and n != "reset"))
                ok = all(np.array_equal(arr(se, n), arr(ab, n)) for n in ("gnew", "gexp_1", "gexp_2", "gexp_3")
                         if n in se.variables and n in ab.variables)
                tally("deferred reset: gnew/gexp restored (numerically equal)", ok)
                ok = True
                for n in ("phinew", "aparnew", "bparnew"):
                    if n in se.variables:
                        a_, b_ = arr(se, n), arr(ab, n)
                        ok = ok and float(np.max(np.abs(a_ - b_))) <= 1e-12 * max(float(np.max(np.abs(a_))), 1e-300)
                tally("deferred reset: fields recomputed from gnew (round-off)", ok)
            elif ab is not None:
                names = [n for n in state_vars if n in se.variables and n in ab.variables]
                # with immediate_reset = .false. the CFL limit (code_dt_cfl) is only completed in
                # check_time_step, after the step_end probe (nonlinear_terms.fpp:571-575)
                skip = ("code_dt_cfl", "user_dt_cfl") if not immediate_reset else ()
                ok = all(bits_equal(arr(se, n), arr(ab, n)) for n in names)
                ok = ok and all(bits_equal(arr(se, n), arr(ab, n)) for n in tb if n not in skip)
                tally("step_end(k) == advance_begin(k+1, attempt 1)", ok)
        for att in range(1, last + 1):
            lst = recs[(k, att)]
            ta = [r[2] for r in lst if r[1] == "after_timeadv"]
            # (b) pass-1 phinew == phi
            if ta:
                g1 = ta[0]
                ok = int(g1.getncattr("pass")) == 1
                for f in ("phi", "apar", "bpar"):
                    if f in g1.variables:
                        ok = ok and bits_equal(arr(g1, f + "new"), arr(g1, f))
                tally("after_timeadv pass 1: phinew == phi (bitwise)", ok)
            # (d) after_field_update == after_getfield + phi (implicit fields: one complex add)
            gf = find(k, att, "after_getfield")
            fu = find(k, att, "after_field_update")
            if gf is not None and fu is not None and ta and ex.fieldopt == 1:
                ok = True
                for f in ("phi", "apar", "bpar"):
                    if f in ta[0].variables:
                        ok = ok and bits_equal(arr(gf, f + "new") + arr(ta[0], f), arr(fu, f + "new"))
                tally("after_field_update == after_getfield + phi (bitwise)", ok)
            # (c) last after_timeadv (pass 2) == advance_end when collisions are not split
            ae = find(k, att, "advance_end")
            if ae is not None and len(ta) >= 2:
                ca = find(k, att, "after_collisions_advance")
                src = ca if ca is not None else ta[-1]
                tally("advance_end gnew == after_timeadv pass 2 (or after_collisions_advance)",
                      bits_equal(arr(src, "gnew"), arr(ae, "gnew")))
            # (e) source precedes after_invert_rhs_1 in each pass; invert_rhs_1 == invert_rhs without linked BCs
            tags = [r[1] for r in lst]
            for i, t in enumerate(tags):
                if t == "after_invert_rhs_1":
                    tally("source record precedes after_invert_rhs_1", i > 0 and tags[i - 1] == "source")
            # (f) nonlinear history shuffle and re-use (05 §4.1, Pitfall 11.18)
            aet = [r[2] for r in lst if r[1] == "after_add_explicit_terms"]
            ab = find(k, att, "advance_begin")
            if aet and ab is not None and "gexp_1" in aet[0].variables and "gexp_1" in ab.variables:
                nlc = find(k, att, "nl_after_center")
                if att == 1 and nlc is not None:
                    ok = bits_equal(arr(aet[0], "gexp_3"), arr(ab, "gexp_2")) and \
                        bits_equal(arr(aet[0], "gexp_2"), arr(ab, "gexp_1")) and \
                        bits_equal(arr(aet[0], "gexp_1"), arr(nlc, "g1"))
                    tally("NL shuffle: gexp_3<-gexp_2, gexp_2<-gexp_1, gexp_1 = nl_after_center g1", ok)
                if len(aet) >= 2:
                    ok = all(bits_equal(arr(aet[0], n), arr(aet[1], n)) for n in ("gexp_1", "gexp_2", "gexp_3"))
                    tally("pass 2 reuses the pass-1 explicit terms", ok)
            if att > 1:
                tally("attempt > 1 has no nl_* records", not any(t.startswith("nl_") for t in tags))
            # (h) collisions: g2le gather / scatter are pure copies over the allowed points (06 §2.3)
            if gmap is not None:
                hin = [r[2] for r in lst if r[1] == "coll_h_in"]
                gin = [r[2] for r in lst if r[1] == "coll_gle_in"]
                gout = [r[2] for r in lst if r[1] == "coll_gle_out"]
                hout = [r[2] for r in lst if r[1] == "coll_h_out"]
                for a_, b_ in zip(hin, gin):
                    tally("coll_gle_in == g2le gather of coll_h_in (bitwise)", bits_equal(gather(cplx(arr(a_, "h"))),
                                                                                         cplx(arr(b_, "gle"))))
                for a_, b_, c_ in zip(gout, hout, hin):
                    tally("coll_h_out == g2le scatter of coll_gle_out, forbidden points kept (bitwise)",
                          bits_equal(scatter(cplx(arr(a_, "gle")), cplx(arr(c_, "h"))), cplx(arr(b_, "h"))))
        # (g) update_time: step_end_pre_diag times
        ae = find(k, last, "advance_end")
        sp = find(k, last, "step_end_pre_diag")
        if ae is not None and sp is not None:
            ct = float(arr(ae, "code_time")) + float(arr(ae, "code_dt"))
            ok = np.float64(arr(sp, "code_time")).tobytes() == np.float64(ct).tobytes() and \
                bits_equal(arr(sp, "code_dt_prev1"), arr(ae, "code_dt")) and \
                bits_equal(arr(sp, "code_dt_prev2"), arr(ae, "code_dt_prev1")) and \
                bits_equal(arr(sp, "code_dt_old"), arr(ae, "code_dt"))
            tally("update_time: code_time += code_dt, prev2<-prev1<-code_dt, code_dt_old = code_dt", ok)
            tally("advance_end gnew == step_end_pre_diag gnew", bits_equal(arr(ae, "gnew"), arr(sp, "gnew")))
    for name, (n, nf) in nfail.items():
        rep.ok(nf == 0, f"{name}", f"{n - nf}/{n} records")
    nresets = sum(1 for k in attempts if attempts[k] > 1)
    if nresets:
        rep.note(f"{nresets} step(s) with an immediate reset (attempt 2)")
    tr.close()


def check_restart(ex, path, rep):
    """09 §12 T2: the canonical map vs GS2's own restart writer (write_restart_t0 = .true.)."""
    print(f"[restart] canonical map vs {os.path.basename(path)} (09 §12 T2)")
    r = netCDF4.Dataset(path)
    r.set_auto_mask(False)
    iglo_of = ex.var("/gs2_layouts/iglo_of")  # (S,E,L,Y,X)
    if "iglo_kykxsel" in r.variables:
        kyk = np.asarray(r["iglo_kykxsel"][...])  # Fortran (aky,akx,species,egrid,lambda) -> NumPy reversed
        rep.ok(np.array_equal(np.transpose(iglo_of, (2, 1, 0, 4, 3)), kyk), "iglo_of == restart iglo_kykxsel")
    g = np.asarray(r["gr"][...]) + 1j * np.asarray(r["gi"][...])  # (glo, sign, theta)
    gn = cplx(ex.var("/state_t0/gnew"))
    rep.ok(bits_equal(g[iglo_of], gn), "state_t0 gnew == restart gr/gi re-indexed with iglo_of (bitwise)")
    if "vpa" in r.variables:
        rep.ok(bits_equal(np.asarray(r["vpa"][...])[iglo_of], ex.var("/dist_fn_arrays/vpa")),
               "vpa == restart vpa re-indexed with iglo_of (bitwise)")
    if "vperp2" in r.variables:
        rep.ok(bits_equal(np.asarray(r["vperp2"][...])[iglo_of], ex.var("/dist_fn_arrays/vperp2")),
               "vperp2 == restart vperp2 re-indexed with iglo_of (bitwise)")
    r.close()


def check_same_export(ex, other_path, rep):
    print(f"[determinism] export file vs {os.path.basename(other_path)}")
    o = netCDF4.Dataset(other_path)
    for _, g in walk(o):
        g.set_auto_mask(False)
    ga, gb = dict(walk(ex.ds)), dict(walk(o))
    n = 0
    bad = []
    for gpath in ga:
        if gpath.startswith("/dt_ladder") or gpath not in gb:
            continue
        for name, va in ga[gpath].variables.items():
            if name not in gb[gpath].variables:
                bad.append(f"{gpath}/{name} missing")
                continue
            n += 1
            if not bits_equal(va[...], gb[gpath].variables[name][...]):
                bad.append(f"{gpath}/{name}")
    missing_groups = [g for g in ga if not g.startswith("/dt_ladder") and g not in gb]
    rep.ok(not bad and not missing_groups, f"{n} variables bitwise identical in both export files",
           ", ".join((bad + missing_groups)[:8]))
    o.close()


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("export")
    p.add_argument("--trace")
    p.add_argument("--outnc", help="out.nc written by the gs2_export run (static part + istep-0 record)")
    p.add_argument("--ref-outnc", help="out.nc of the full bin/gs2 run")
    p.add_argument("--trace-outnc", help="out.nc of the traced gs2_export run (compared with --ref-outnc)")
    p.add_argument("--compare-export", help="second export file that must be identical (except /dt_ladder)")
    p.add_argument("--restart", help="GS2 restart file written with write_restart_t0 = .true. (T2 check)")
    p.add_argument("--quiet", action="store_true")
    a = p.parse_args(argv)
    rep = Report(a.quiet)
    ex = Export(a.export)
    print(f"check_export: {a.export}  (layout {ex.ds.getncattr('layout')}, "
          f"{'reference' if ex.reference_build else 'stock'} build)")
    check_schema(ex, rep)
    outnc = a.outnc
    if outnc is None:
        cand = os.path.join(os.path.dirname(a.export) or ".", f"{ex.ds.getncattr('run_name')}.out.nc")
        if os.path.exists(cand):
            outnc = cand
    check_numerics(ex, rep, outnc, a.ref_outnc)
    if a.trace:
        check_trace(ex, a.trace, rep)
    if a.compare_export:
        check_same_export(ex, a.compare_export, rep)
    if a.restart:
        check_restart(ex, a.restart, rep)
    if a.trace_outnc and a.ref_outnc:
        print("[trajectory] traced gs2_export run vs bin/gs2 (time records with equal t)")
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        import compare_netcdf
        ok = compare_netcdf.compare(a.ref_outnc, a.trace_outnc, records_of_b=True, time_only=True, quiet=True)
        rep.ok(ok, "traced run reproduces bin/gs2's time-dependent output bitwise")
    print(f"check_export: {rep.n_pass} passed, {len(rep.fails)} failed"
          + (f", notes: {len(rep.notes)}" if rep.notes else ""))
    return 0 if not rep.fails else 1


if __name__ == "__main__":
    sys.exit(main())
