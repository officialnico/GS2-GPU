"""Tests of gs2gpu.diagnostics and gs2gpu.outnc (spec 08; 10 B.7, G-05).

Three groups:

1. Synthetic (always run): the vectorised device/host implementation against a literal,
   loop-by-loop transcription of the Fortran (``_lit_*`` below: ``iglo`` loops in layout
   order, scalar Klein sums, explicit real/imaginary arithmetic), for several layouts and
   field sets. In strict mode (``numerics.set_strict(True)``, no FMA) the two must agree
   bitwise; this checks summation orders, index conventions and association order.
2. Format: the writer's variables, dimensions and attributes against GS2's own
   ``<case>.out.nc`` (runs whenever ``<REFDATA>/<case>/<case>.out.nc`` exists; with no
   export file the arrays come from a synthetic export with GS2's grids, which is enough
   for the format).
3. Reference data (skip when missing): the ``g_adjust`` round trips from the trace
   (``step_end_pre_diag`` -> ``step_end``) and from ``/state_t0_prediag`` -> ``/state_t0``,
   and every written variable against the same record of GS2's ``<case>.out.nc``.
"""

from __future__ import annotations

import math
from functools import lru_cache
from pathlib import Path

import jax
import jax.numpy as jnp
import netCDF4
import numpy as np
import pytest

from conftest import (CASES_DIR, REFDATA, RTOL_KERNEL, all_cases, assert_close, export_of, have,
                      out_of, trace_of)
from gs2gpu import diagnostics as dg
from gs2gpu import numerics, outnc
from gs2gpu.io.export import Export, ExportError, Group
from gs2gpu.types import Fields, fields_from_group, to_jax

# ============================================================================ synthetic export


def _grp(path, arrays, attrs=None):
    return Group(name=path.rsplit("/", 1)[-1], path=path, arrays=dict(arrays),
                 attrs=dict(attrs or {}))


def synthetic_export(*, S=2, E=3, L=4, Y=3, X=4, ntgrid=4, layout="lxyes", fphi=1.0, fapar=0.0,
                     fbpar=0.0, nonlin=False, seed=0, flags=None, zonal=True, grids=None) -> Export:
    """A self-consistent export with random but physical-looking arrays.

    ``grids`` (optional dict) overrides grid arrays (used to take GS2's grids from an out.nc).
    """
    rng = np.random.default_rng(seed)
    g = dict(grids or {})
    T = 2 * ntgrid + 1
    theta = g.get("theta", np.linspace(-np.pi, np.pi, T))
    T = len(theta)
    ntgrid = (T - 1) // 2
    aky = g.get("aky", np.concatenate([[0.0] if zonal else [], 0.2 * np.arange(1, Y + 1)])[:Y])
    akx = g.get("akx", 0.3 * np.fft.fftfreq(X, 1.0 / X))
    Y, X = len(aky), len(akx)
    bmag = g.get("bmag", 1.0 / (1.0 + 0.18 * np.cos(theta)))
    gradpar = g.get("gradpar", 0.7 + 0.05 * np.cos(theta))
    jacob = g.get("jacob", 1.0 / (1.1 * gradpar * bmag))
    grho = g.get("grho", 1.0 + 0.1 * np.sin(theta) ** 2)
    bpol = g.get("bpol", 0.12 * bmag)
    rplot = g.get("rplot", 3.0 + 0.18 * np.cos(theta))
    gds21 = g.get("gds21", -0.8 * theta * 0.3)
    gds22 = g.get("gds22", 0.64 * (1 + 0 * theta))
    shat = g.get("shat", 0.8)
    qval = g.get("qval", 1.4)
    theta0 = g.get("theta0", np.outer(np.where(aky == 0, 0.0, 1.0 / np.where(aky == 0, 1, aky)),
                                      akx) / shat)
    kperp2 = g.get("kperp2", (aky[:, None, None] ** 2 + (akx[None, :, None]
                                                           + shat * aky[:, None, None] * theta) ** 2)
                   * (1 + 0.1 * np.cos(theta)))
    energy = g.get("energy", np.tile(np.linspace(0.3, 4.0, E), (S, 1)))
    S, E = energy.shape
    L = len(g.get("al", np.arange(L)))
    z = g.get("charge", np.array([1.0, -1.0, 2.0])[:S])
    mass = g.get("mass", np.array([1.0, 0.01, 2.0])[:S])
    dens = g.get("dens", np.array([1.0, 1.1, 0.3])[:S])
    temp = g.get("temp", np.array([1.0, 0.9, 1.2])[:S])
    tprim = g.get("tprim", np.array([6.9, 5.0, 2.0])[:S])
    fprim = g.get("fprim", np.array([2.2, 2.0, 1.0])[:S])
    stm, zt, tz = np.sqrt(temp / mass), z / temp, temp / z
    smz = np.abs(np.sqrt(temp * mass) / z)
    delthet = np.zeros(T)
    delthet[:-1] = theta[1:] - theta[:-1]
    w = rng.uniform(0.1, 1.0, (S, E))
    jend = np.full(T, L, np.int32)
    if L > 2:
        jend[: max(1, ntgrid - 1)] = L - 1
        jend[-max(1, ntgrid - 1):] = L - 1
    forbid = np.zeros((L, T), bool)
    for ig in range(T):
        forbid[jend[ig]:, ig] = True
    wl = np.where(forbid, 0.0, rng.uniform(0.0, 1.0, (L, T)))
    vpa1 = rng.uniform(0.1, 2.0, (S, E, L, Y, X, T))
    vpa = np.stack([vpa1, -vpa1], axis=5)
    vpac = 0.5 * (vpa + np.roll(vpa, -1, axis=-1))
    vpac[..., -1] = 0.0
    vperp2 = rng.uniform(0.0, 3.0, (S, E, L, Y, X, T))
    aj0 = rng.uniform(0.2, 1.0, (S, E, L, Y, X, T))
    aj1 = rng.uniform(0.1, 0.5, (S, E, L, Y, X, T))
    nmc = rng.uniform(0.9, 1.1, (S, E))
    root = Group(name="/", path="/")

    def add(path, arrays, attrs=None):
        root.groups[path.strip("/")] = _grp(path, arrays, attrs)

    add("/theta_grid", dict(theta=theta, delthet=delthet, bmag=bmag, gradpar=gradpar, jacob=jacob,
                            grho=grho, Bpol=bpol, Rplot=rplot, gds21=gds21, gds22=gds22,
                            shat=shat, qval=qval, rhoc=0.5))
    add("/kt_grids", dict(aky=aky, akx=akx, theta0=theta0, kperp2=kperp2))
    add("/le_grids", dict(w=w, wl=wl, energy=energy, jend=jend, forbid=forbid,
                          grid_has_trapped_particles=True, ng2=max(L - 2, 1), nxi=2 * L))
    add("/species", dict(spec_z=z, spec_mass=mass, spec_dens=dens, spec_temp=temp,
                         spec_tprim=tprim, spec_fprim=fprim, spec_stm=stm, spec_zt=zt,
                         spec_tz=tz, spec_smz=smz, spec_type=np.ones(S, np.int32),
                         nonmaxw_corr=nmc))
    add("/dist_fn_arrays", dict(vpa=vpa, vpac=vpac, vperp2=vperp2, aj0=aj0, aj1=aj1))
    eps = np.finfo(float).eps
    add("/run_parameters", dict(fphi=fphi, fapar=fapar, fbpar=fbpar, has_phi=fphi > eps,
                                has_apar=fapar > eps, has_bpar=fbpar > eps))
    add("/nonlinear_terms", dict(nonlin=nonlin))
    add("/dist_fn", dict(def_parity=False, even=True))
    add("/gs2_layouts", {}, {"layout": layout})
    diag = dict(dg.GNOSTICS_DEFAULTS)
    diag.update(write_ascii=False, write_fluxes=True, write_fluxes_by_mode=True,
                write_eigenfunc=True, write_final_fields=True, write_final_epar=True,
                write_final_moments=True, write_phi_over_time=True, write_apar_over_time=True,
                write_bpar_over_time=True, write_ntot_over_time=True, write_upar_over_time=True,
                write_density_over_time=True, write_tperp_over_time=True, omegatol=-1e-3,
                nwrite=2, navg=3)
    diag.update(flags or {})
    add("/diagnostics", {k: v for k, v in diag.items()})
    dims = {"theta": T, "kx": X, "ky": Y, "lambda": L, "egrid": E, "species": S, "ri": 2}
    return Export(path=Path("synthetic"), root=root, dims=dims, attrs={"layout": layout})


def random_state(ex: Export, seed=1, old_scale=0.93):
    rng = np.random.default_rng(seed)

    def c(shape):
        return rng.standard_normal(shape) + 1j * rng.standard_normal(shape)
    gnew, g = c(ex.g_shape), c(ex.g_shape)
    phin, aparn, bparn = c(ex.field_shape), 0.1 * c(ex.field_shape), 0.01 * c(ex.field_shape)
    rot = np.exp(-0.07j)
    fold = Fields(phin * old_scale * rot, aparn * old_scale * rot, bparn * old_scale * rot)
    fnew = Fields(phin, aparn, bparn)
    # absent fields are exactly zero
    rp = ex["/run_parameters"]
    fnew = Fields(*(f if bool(rp[f"has_{n}"]) else np.zeros_like(f)
                    for n, f in zip(("phi", "apar", "bpar"), fnew)))
    fold = Fields(*(f if bool(rp[f"has_{n}"]) else np.zeros_like(f)
                    for n, f in zip(("phi", "apar", "bpar"), fold)))
    # g is zero in the forbidden region (GS2 never sets it there)
    le = ex["/le_grids"]
    mask = ~(le["grid_has_trapped_particles"]
             & (np.arange(1, ex.nlambda + 1)[:, None] > le["jend"][None, :]))
    m7 = mask[None, None, :, None, None, None, :]
    return np.where(m7, gnew, 0), np.where(m7, g, 0), fold, fnew


# ============================================================================ literal reference
# Straight transcriptions of the Fortran loops (Python scalar Klein sums, iglo loops in layout
# order, theta handled as elementwise NumPy array syntax with explicit re/im arithmetic).


def _decode(layout):
    """iglo -> (it, ik, il, ie, is) (0-based) for gs2_layouts' idx with this layout."""
    def gen(sizes):
        names = {"x": "it", "y": "ik", "l": "il", "e": "ie", "s": "is"}
        order = [names[c] for c in layout]
        n = int(np.prod([sizes[o] for o in order]))
        for iglo in range(n):
            rem = iglo
            idx = {}
            for o in order:
                idx[o] = rem % sizes[o]
                rem //= sizes[o]
            yield iglo, idx["it"], idx["ik"], idx["il"], idx["ie"], idx["is"]
    return gen


def _lit_klein(a):
    res = c1 = c2 = 0.0
    for x in a:
        x = float(x)
        t = res + x
        lc1 = (res - t) + x if abs(res) >= abs(x) else (x - t) + res
        res = t
        t = c1 + lc1
        lc2 = (c1 - t) + lc1 if abs(c1) >= abs(lc1) else (lc1 - t) + c1
        c1 = t
        c2 = c2 + lc2
    return res + c1 + c2


def _lit_trapz(x, y):
    return _lit_klein([(x[i + 1] - x[i]) * (y[i] + y[i + 1]) for i in range(len(x) - 1)]) * 0.5


class Lit:
    """Literal Fortran transcription of the diagnostics on a (synthetic) export."""

    def __init__(self, ex: Export, params: dg.DiagParams):
        self.ex, self.p = ex, params
        tg, kt, le, sp, da = (ex[k] for k in ("/theta_grid", "/kt_grids", "/le_grids",
                                               "/species", "/dist_fn_arrays"))
        self.theta = np.asarray(tg["theta"])
        self.jacob = np.asarray(tg["jacob"])
        self.grho = np.asarray(tg["grho"])
        self.bmag = np.asarray(tg["bmag"])
        self.gradpar = np.asarray(tg["gradpar"])
        self.delthet = np.asarray(tg["delthet"])
        self.gds21, self.gds22 = np.asarray(tg["gds21"]), np.asarray(tg["gds22"])
        self.Rplot, self.Bpol = np.asarray(tg["Rplot"]), np.asarray(tg["Bpol"])
        self.rhoc, self.qval, self.shat = float(tg["rhoc"]), float(tg["qval"]), float(tg["shat"])
        self.aky, self.akx = np.asarray(kt["aky"]), np.asarray(kt["akx"])
        self.theta0, self.kperp2 = np.asarray(kt["theta0"]), np.asarray(kt["kperp2"])
        self.w, self.wl, self.energy = np.asarray(le["w"]), np.asarray(le["wl"]), np.asarray(le["energy"])
        self.jend, self.forbid = np.asarray(le["jend"]), np.asarray(le["forbid"])
        self.trapped = bool(le["grid_has_trapped_particles"])
        self.sp = {k[5:]: np.asarray(v) for k, v in sp.arrays.items() if k.startswith("spec_")}
        self.nm = np.asarray(sp["nonmaxw_corr"])
        self.vpa, self.vpac = np.asarray(da["vpa"]), np.asarray(da["vpac"])
        self.vperp2, self.aj0, self.aj1 = (np.asarray(da[k]) for k in ("vperp2", "aj0", "aj1"))
        self.iglos = list(_decode(params.layout)({"it": ex.ntheta0, "ik": ex.naky, "il": ex.nlambda,
                                                   "ie": ex.negrid, "is": ex.nspec}))
        # theta_grid.f90:520, diagnostics_fluxes.fpp:825
        self.W = 1.0 / _lit_trapz(self.theta, self.jacob)
        self.wgt = 1.0 / _lit_trapz(self.theta, self.grho * self.jacob)

    # -- g_adjust (dist_fn_arrays.fpp:118-168)
    def g_adjust(self, g, phi, bpar, facphi, facbpar):
        g = g.copy()
        with_bpar, with_phi = facbpar != 0, facphi != 0
        bpar_factor = 2.0 * facbpar if with_bpar else 0.0
        phi_factor = 0.0
        for iglo, it, ik, il, ie, is_ in self.iglos:
            if with_phi:
                phi_factor = self.sp["zt"][is_] * self.nm[is_, ie] * facphi
            for ig in range(len(self.theta)):
                if self.trapped and il + 1 > self.jend[ig]:
                    continue
                cb = bpar_factor * self.vperp2[is_, ie, il, ik, it, ig] * self.aj1[is_, ie, il, ik, it, ig]
                b = bpar[ik, it, ig]
                p = phi[ik, it, ig]
                a0 = self.aj0[is_, ie, il, ik, it, ig]
                adj_r = cb * b.real + (phi_factor * p.real) * a0
                adj_i = cb * b.imag + (phi_factor * p.imag) * a0
                for isgn in range(2):
                    v = g[is_, ie, il, ik, it, isgn, ig]
                    g[is_, ie, il, ik, it, isgn, ig] = complex(v.real + adj_r, v.imag + adj_i)
        return g

    # -- integrate_moment (le_grids.f90 integrate_moment_c34)
    def integrate_moment(self, gw):
        S, Y, X, T = self.ex.nspec, self.ex.naky, self.ex.ntheta0, self.ex.nth
        tr = np.zeros((S, Y, X, T))
        ti = np.zeros((S, Y, X, T))
        for iglo, it, ik, il, ie, is_ in self.iglos:
            ww = self.w[is_, ie] * self.wl[il]
            s = gw[is_, ie, il, ik, it, 0] + gw[is_, ie, il, ik, it, 1]
            tr[is_, ik, it] = tr[is_, ik, it] + ww * s.real
            ti[is_, ik, it] = ti[is_, ik, it] + ww * s.imag
        return tr + 1j * ti

    # -- per-iglo g_work builders (diagnostics_fluxes.fpp flux / flux_vs_e / eexchange,
    #    diagnostics_moments.fpp getmoms); h, ho G7 complex
    def gwork(self, kind, h, ho=None, phinew=None):
        out = np.zeros(h.shape, complex)
        for iglo, it, ik, il, ie, is_ in self.iglos:
            aj0 = self.aj0[is_, ie, il, ik, it]
            aj1 = self.aj1[is_, ie, il, ik, it]
            vp2 = self.vperp2[is_, ie, il, ik, it]
            en = self.energy[is_, ie]
            stm, tz, smz, zt = (self.sp[k][is_] for k in ("stm", "tz", "smz", "zt"))
            for isgn in range(2):
                g = h[is_, ie, il, ik, it, isgn]
                gr, gi = g.real, g.imag
                vpa = self.vpa[is_, ie, il, ik, it, isgn]
                vpac = self.vpac[is_, ie, il, ik, it, isgn]
                if kind == "pflux":
                    r, i = gr * aj0, gi * aj0
                elif kind == "qflux1":
                    r, i = (gr * aj0) * en, (gi * aj0) * en
                elif kind == "qflux2":
                    v2 = vpa * vpa
                    r, i = ((gr * 2.0) * v2) * aj0, ((gi * 2.0) * v2) * aj0
                elif kind == "qflux3":
                    r, i = (gr * vp2) * aj0, (gi * vp2) * aj0
                elif kind == "vflux_par":
                    sq = np.sqrt(1.0 - (self.Bpol * self.Bpol) / (self.bmag * self.bmag))
                    r = (((gr * aj0) * vpac) * self.Rplot) * sq
                    i = (((gi * aj0) * vpac) * self.Rplot) * sq
                elif kind == "vflux_perp":
                    a = self.aky[ik]
                    r, i = -(a * gi), a * gr                      # (zi*aky)*gnew
                    geo = self.gds21 + self.theta0[ik, it] * self.gds22
                    den = (self.qval * self.shat) * (self.bmag * self.bmag)
                    for f in (aj1, self.rhoc, geo, vp2, smz):
                        r, i = r * f, i * f
                    r, i = -(r / den), -(i / den)
                elif kind in ("pmflux", "qmflux1", "qmflux2", "qmflux3", "vmflux"):
                    br, bi = ((gr * aj0) * stm) * vpa, ((gi * aj0) * stm) * vpa
                    if kind == "pmflux":
                        r, i = -br, -bi
                    elif kind == "qmflux1":
                        r, i = (-br) * en, (-bi) * en
                    elif kind == "qmflux2":
                        r, i = -((br * 2.0) * (vpa * vpa)), -((bi * 2.0) * (vpa * vpa))
                    elif kind == "qmflux3":
                        r, i = -(br * vp2), -(bi * vp2)
                    else:
                        r, i = -(br * vpac), -(bi * vpac)
                elif kind in ("pbflux", "qbflux1", "qbflux2", "qbflux3", "vbflux"):
                    br = (((gr * aj1) * 2.0) * vp2) * tz
                    bi = (((gi * aj1) * 2.0) * vp2) * tz
                    if kind == "pbflux":
                        r, i = br, bi
                    elif kind == "qbflux1":
                        r, i = br * en, bi * en
                    elif kind == "qbflux2":
                        r, i = (br * 2.0) * (vpa * vpa), (bi * 2.0) * (vpa * vpa)
                    elif kind == "qbflux3":
                        r, i = br * vp2, bi * vp2
                    else:
                        r, i = br * vpac, bi * vpac
                elif kind in ("exch_total", "exch_total2"):
                    o = ho[is_, ie, il, ik, it, isgn]
                    sgn = 1.0 if kind == "exch_total" else -1.0
                    s_r = gr + o.real if sgn > 0 else gr - o.real
                    s_i = gi + o.imag if sgn > 0 else gi - o.imag
                    a = aj0 * 0.25
                    r, i = a * s_r, a * s_i
                elif kind == "density":
                    r, i = aj0 * gr, aj0 * gi
                elif kind == "upar":
                    r, i = vpa * (aj0 * gr), vpa * (aj0 * gi)
                elif kind == "tpar":
                    r, i = (2.0 * vpa) * (vpa * (aj0 * gr)), (2.0 * vpa) * (vpa * (aj0 * gi))
                elif kind == "tperp":
                    r, i = (vp2 * gr) * aj0, (vp2 * gi) * aj0
                elif kind == "qparflux":
                    r, i = ((vpa * gr) * aj0) * en, ((vpa * gi) * aj0) * en
                elif kind in ("pperpj1", "qpperpj1"):
                    r = (((gr * aj1) * 2.0) * vp2) * tz
                    i = (((gi * aj1) * 2.0) * vp2) * tz
                    if kind == "qpperpj1":
                        r, i = r * en, i * en
                elif kind == "ntot":
                    p = phinew[ik, it]
                    zr = (p.real * zt) * self.nm[is_, ie]
                    zi = (p.imag * zt) * self.nm[is_, ie]
                    r, i = aj0 * gr - zr, aj0 * gi - zi
                else:
                    raise KeyError(kind)
                out[is_, ie, il, ik, it, isgn] = r + 1j * i
        return out

    # -- theta averages
    def flavg(self, q):
        return _lit_trapz(self.theta, q * self.jacob) * self.W

    def get_flux(self, total, fld):
        S, Y, X = total.shape[:3]
        flx = np.zeros((S, Y, X))
        for is_ in range(S):
            for ik in range(Y):
                for it in range(X):
                    t, f = total[is_, ik, it], fld[ik, it]
                    im = t.real * (-f.imag) + t.imag * f.real
                    flx[is_, ik, it] = (_lit_trapz(self.theta, im * self.jacob) * self.wgt) * self.aky[ik]
        return flx * 0.5

    def get_flux_vs_e(self, gw, fld):
        S, E = self.ex.nspec, self.ex.negrid
        total = np.zeros((E, S))
        for iglo, it, ik, il, ie, is_ in self.iglos:
            fac = 1.0 if self.aky[ik] == 0.0 else 0.5
            for isgn in range(2):
                g, f = gw[is_, ie, il, ik, it, isgn], fld[ik, it]
                im = g.real * (-f.imag) + g.imag * f.real
                y = (im * self.wl[il]) * self.jacob
                total[ie, is_] = total[ie, is_] + ((fac * _lit_trapz(self.theta, y)) * self.aky[ik]) * self.wgt
        return (0.5 * total).T                                  # -> (S, E)

    # volume averages (volume_averages.f90)
    def avg_ky(self, f):            # f[it, ik] Fortran order given as numpy (Y, X)
        Y, X = f.shape
        out = np.zeros(X)
        for ik in range(Y):
            fac = 1.0 if self.aky[ik] == 0.0 else 0.5
            for it in range(X):
                out[it] = out[it] + f[ik, it] * fac
        return out

    def avg_kx(self, f):
        Y, X = f.shape
        out = np.zeros(Y)
        for ik in range(Y):
            fac = 1.0 if self.aky[ik] == 0.0 else 0.5
            for it in range(X):
                out[ik] = out[ik] + f[ik, it] * fac
        return out

    @staticmethod
    def fsum(v):
        s = 0.0
        for x in v:
            s = s + x
        return s

    def standard(self, name, F):
        S = F.shape[0]
        by_s = np.array([self.fsum(self.avg_ky(F[is_])) for is_ in range(S)])
        by_mode = np.zeros(F.shape[1:])
        for ik in range(F.shape[1]):
            for it in range(F.shape[2]):
                for is_ in range(S):
                    by_mode[ik, it] = by_mode[ik, it] + F[is_, ik, it]
        by_kx = self.avg_ky(by_mode)
        return {name: by_s, f"total_{name}_by_ky": self.avg_kx(by_mode),
                f"total_{name}_by_kx": by_kx, f"total_{name}": self.fsum(by_kx),
                f"{name}_by_mode": F, f"total_{name}_by_mode": by_mode}


def _lit_record(ex, params, gnew, g, fold, fnew, omega_average, code_dt):
    """Fluxes, moments and field families of one write block, literally; returns
    (record, gnew_after, g_after)."""
    L = Lit(ex, params)
    p = params
    rec = {}
    hn = L.g_adjust(gnew, fnew.phi, fnew.bpar, p.fphi, p.fbpar)
    ho = L.g_adjust(g, fold.phi, fold.bpar, p.fphi, p.fbpar)
    sp = L.sp
    DT = sp["dens"] * sp["temp"]
    MOM = sp["dens"] * np.sqrt(sp["mass"] * sp["temp"])
    MOMP = (sp["dens"] * sp["mass"]) * sp["stm"]
    b = lambda v: v[:, None, None]   # noqa: E731
    fl, fve = {}, {}
    if p.has_phi:
        F = {k: L.get_flux(L.integrate_moment(L.gwork(k, hn)), fnew.phi)
             for k in ("pflux", "qflux1", "qflux2", "qflux3", "vflux_par", "vflux_perp")}
        t1 = L.integrate_moment(L.gwork("exch_total", hn, ho))
        t2 = L.integrate_moment(L.gwork("exch_total2", hn, ho))
        S, Y, X = t1.shape[:3]
        ex_ = np.zeros((S, Y, X))
        for is_ in range(S):
            for ik in range(Y):
                for it in range(X):
                    if p.nonlin and it == 0 and ik == 0:
                        continue
                    a, c = t1[is_, ik, it], t2[is_, ik, it]
                    pn, po = fnew.phi[ik, it], fold.phi[ik, it]
                    dr, di = pn.real - po.real, pn.imag - po.imag
                    sr, si = pn.real + po.real, pn.imag + po.imag
                    q = (a.real * dr - a.imag * (-di)) - (sr * c.real - si * (-c.imag))
                    ex_[is_, ik, it] = L.flavg(q) / code_dt
        vfl = F["vflux_par"] + F["vflux_perp"]
        fl.update(es_heat_flux=F["qflux1"] * b(DT), es_heat_flux_par=F["qflux2"] * b(DT),
                  es_heat_flux_perp=F["qflux3"] * b(DT), es_part_flux=F["pflux"] * b(sp["dens"]),
                  es_part_tormom_flux=np.zeros_like(F["pflux"]), es_mom_flux=vfl * b(MOM),
                  es_mom_flux_par=F["vflux_par"] * b(MOM),
                  es_energy_exchange=(ex_ * b(sp["dens"])) * b(sp["z"]),
                  es_mom_flux_perp=F["vflux_perp"] * b(MOMP))
        fve["es_flux_vs_e"] = L.get_flux_vs_e(L.gwork("pflux", hn), fnew.phi) * sp["dens"][:, None]
    if p.has_apar:
        F = {k: L.get_flux(L.integrate_moment(L.gwork(k, hn)), fnew.apar)
             for k in ("pmflux", "qmflux1", "qmflux2", "qmflux3", "vmflux")}
        fl.update(apar_heat_flux=F["qmflux1"] * b(DT), apar_heat_flux_par=F["qmflux2"] * b(DT),
                  apar_heat_flux_perp=F["qmflux3"] * b(DT),
                  apar_part_flux=F["pmflux"] * b(sp["dens"]), apar_mom_flux=F["vmflux"] * b(MOMP))
        fve["apar_flux_vs_e"] = L.get_flux_vs_e(L.gwork("pmflux", hn), fnew.apar) * sp["dens"][:, None]
    if p.has_bpar:
        F = {k: L.get_flux(L.integrate_moment(L.gwork(k, hn)), fnew.bpar)
             for k in ("pbflux", "qbflux1", "qbflux2", "qbflux3", "vbflux")}
        fl.update(bpar_heat_flux=F["qbflux1"] * b(DT), bpar_heat_flux_par=F["qbflux2"] * b(DT),
                  bpar_heat_flux_perp=F["qbflux3"] * b(DT),
                  bpar_part_flux=F["pbflux"] * b(sp["dens"]), bpar_mom_flux=F["vbflux"] * b(MOMP))
        fve["bpar_flux_vs_e"] = L.get_flux_vs_e(L.gwork("pbflux", hn), fnew.bpar) * sp["dens"][:, None]
    tot = {"heat": 0.0, "mom": 0.0, "part": 0.0}
    for name, val in fl.items():
        r = L.standard(name, val)
        rec.update(r)
        for key, names in (("heat", ("es_heat_flux", "apar_heat_flux", "bpar_heat_flux")),
                           ("mom", ("es_mom_flux", "apar_mom_flux", "bpar_mom_flux")),
                           ("part", ("es_part_flux", "apar_part_flux", "bpar_part_flux"))):
            if name in names:
                tot[key] = tot[key] + L.fsum(r[name])
    rec.update(fve)
    rec.update(heat_flux_tot=tot["heat"], mom_flux_tot=tot["mom"], part_flux_tot=tot["part"])
    if not p.nonlin:
        ig = p.ig
        Y, X = omega_average.shape
        D = np.zeros((Y, X))
        S = len(sp["dens"])
        Q = np.zeros((S, Y, X))
        P = np.zeros((S, Y, X))
        for ik in range(Y):
            for it in range(X):
                if L.kperp2[ik, it, ig] == 0.0:
                    continue
                D[ik, it] = max(omega_average[ik, it].imag, 0.0) / L.kperp2[ik, it, ig] * 2.0
                for is_ in range(S):
                    Q[is_, ik, it] = D[ik, it] * sp["dens"][is_] * sp["temp"][is_] * sp["tprim"][is_] * L.grho[ig]
                    P[is_, ik, it] = D[ik, it] * (sp["dens"][is_] * sp["dens"][is_]) * sp["fprim"][is_] * L.grho[ig]
        rec.update(L.standard("heat_flux_diff", Q))
        rec.update(L.standard("part_flux_diff", P))
        rec["diffusivity"] = np.max(D) * L.grho[ig]
        rec["heat_flux_diff_max"] = L.fsum(Q.max(axis=2).max(axis=1))
    gnew1 = L.g_adjust(hn, fnew.phi, fnew.bpar, -p.fphi, -p.fbpar)
    g1 = L.g_adjust(ho, fold.phi, fold.bpar, -p.fphi, -p.fbpar)
    # getmoms
    h2 = L.g_adjust(gnew1, fnew.phi, fnew.bpar, p.fphi, p.fbpar)
    raw = {k: L.integrate_moment(L.gwork(k, h2, phinew=fnew.phi)) for k in dg.MOMENT_TERMS}
    raw["tpar"] = raw["tpar"] - raw["density"]
    raw["tperp"] = raw["tperp"] - raw["density"]
    bb = lambda v: v[:, None, None, None]   # noqa: E731
    mom = {"ntot": raw["ntot"] * bb(sp["dens"]), "density": raw["density"] * bb(sp["dens"]),
           "upar": raw["upar"] * bb(sp["stm"]), "tpar": raw["tpar"] * bb(sp["temp"]),
           "tperp": raw["tperp"] * bb(sp["temp"]),
           "qparflux": ((raw["qparflux"] * bb(sp["dens"])) * bb(sp["temp"])) * bb(sp["stm"]),
           "pperpj1": (raw["pperpj1"] * bb(sp["dens"])) * bb(sp["tz"]),
           "qpperpj1": ((raw["qpperpj1"] * bb(sp["dens"])) * bb(sp["temp"])) * bb(sp["tz"])}
    for name, M in mom.items():
        S, Y, X = M.shape[:3]
        m2 = np.zeros((S, Y, X))
        mbm = np.zeros((S, Y, X), complex)
        for is_ in range(S):
            for ik in range(Y):
                for it in range(X):
                    a = M[is_, ik, it]
                    m2[is_, ik, it] = L.flavg(a.real * a.real - (-a.imag) * a.imag)
                    mbm[is_, ik, it] = complex(_lit_trapz(L.theta, a.real * L.jacob) * L.W,
                                               _lit_trapz(L.theta, a.imag * L.jacob) * L.W)
        rec[name] = M
        rec[f"{name}2_by_mode"] = m2
        rec[f"{name}_by_mode"] = mbm
        rec[f"{name}2_by_ky"] = np.array([L.avg_kx(m2[s]) for s in range(S)])
        rec[f"{name}2_by_kx"] = np.array([L.avg_ky(m2[s]) for s in range(S)])
        rec[f"{name}2"] = np.array([L.fsum(L.avg_ky(m2[s])) for s in range(S)])
    gnew2 = L.g_adjust(h2, fnew.phi, fnew.bpar, -p.fphi, -p.fbpar)
    # fields
    for name in p.fields_present:
        f = getattr(fnew, name)
        Y, X = f.shape[:2]
        f2 = np.array([[L.flavg(f[ik, it].real * f[ik, it].real - (-f[ik, it].imag) * f[ik, it].imag)
                        for it in range(X)] for ik in range(Y)])
        rec[f"{name}2_by_mode"] = f2
        rec[f"{name}2_by_kx"] = L.avg_ky(f2)
        rec[f"{name}2_by_ky"] = L.avg_kx(f2)
        rec[f"{name}2"] = L.fsum(L.avg_ky(f2))
    if p.has_phi:
        Y, X, T = fnew.phi.shape
        cf = np.zeros((Y, T), complex)
        for ik in range(Y):
            fac = 1.0 if ik == 0 else 0.5
            for it in range(X):
                a = fnew.phi[ik, it, p.ntgrid]
                c = fnew.phi[ik, it]
                re = a.real * c.real - a.imag * (-c.imag)
                im = a.real * (-c.imag) + a.imag * c.real
                cf[ik] = (cf[ik].real + re * fac) + 1j * (cf[ik].imag + im * fac)
        rec["phi_corr_2pi"] = cf
    return rec, gnew2, g1


# ============================================================================ helpers


def _device(strict: bool):
    """The jitted device kernel for this mode (strict mode is read at trace time, so each mode
    has its own jit wrapper and every call sets the flag before tracing)."""
    numerics.set_strict(strict)
    if strict not in _KERNELS:
        # a distinct Python function per mode: jax.jit shares its trace cache between wrappers
        # of the same function, so jit(f) twice would reuse the first mode's trace
        def kernel(gnew, g, fields_new, fields_old, data, params):
            return dg._diagnostics_device_impl(gnew, g, fields_new, fields_old, data, params)
        _KERNELS[strict] = jax.jit(kernel, static_argnames=("params",))
    return _KERNELS[strict]


def _final_moments_kernel(strict):
    numerics.set_strict(strict)
    key = ("final", strict)
    if key not in _KERNELS:
        def kernel(gnew, fields_new, data, params):
            return dg._final_moments_impl(gnew, fields_new, data, params)
        _KERNELS[key] = jax.jit(kernel, static_argnames=("params",))
    return _KERNELS[key]


_KERNELS: dict = {}


def _run_block(ex, params, data, gnew, g, fold, fnew, *, istep, code_dt, user_time, state,
               strict=False, return_dev=False):
    kernel = _device(strict)
    try:
        fn = Fields(*map(jnp.asarray, fnew))
        fo = Fields(*map(jnp.asarray, fold))
        gn, go, dev = kernel(jnp.asarray(gnew), jnp.asarray(g), fn, fo, to_jax(data), params)
        dev = {k: np.asarray(v) for k, v in dev.items()}
        host = dg.HostData.from_data(data)
        rec = dg.compute_record(dev, dg._np_fields(fold), dg._np_fields(fnew), state, code_dt,
                                user_time, istep, params, host)
    finally:
        numerics.set_strict(False)
    if return_dev:
        return np.asarray(gn), np.asarray(go), rec, dev
    return np.asarray(gn), np.asarray(go), rec


def _maxrel(ref, got):
    ref, got = np.asarray(ref), np.asarray(got)
    return numerics.max_abs_rel(ref, got)


SYNTH_CONFIGS = [
    dict(layout="lxyes", fphi=1.0, fapar=0.0, fbpar=0.0, nonlin=False),
    dict(layout="yxels", fphi=1.0, fapar=1.0, fbpar=1.0, nonlin=False, S=1),
    dict(layout="xyles", fphi=1.0, fapar=1.0, fbpar=0.0, nonlin=True),
    dict(layout="lexys", fphi=1.0, fapar=0.0, fbpar=1.0, nonlin=True, S=2),
    dict(layout="yxles", fphi=1.0, fapar=1.0, fbpar=1.0, nonlin=False),
    dict(layout="lyxes", fphi=1.0, fapar=0.0, fbpar=0.0, nonlin=False, zonal=False),
]


# ============================================================================ 1. synthetic


def test_klein_and_trapz_agree():
    """Host vectorised Klein/trapz == device lax.scan == scalar loop, bitwise."""
    rng = np.random.default_rng(5)
    a = rng.standard_normal((40, 7)) * 10.0 ** rng.uniform(-8, 8, (40, 7))
    lit = np.array([_lit_klein(a[:, j]) for j in range(7)])
    host = dg.host_klein_sum(a)
    dev = np.asarray(jax.jit(lambda x: dg.klein_sum_jax(x, 0))(a))
    assert np.array_equal(lit, host)
    assert np.array_equal(lit, dev)
    x = np.sort(rng.uniform(-3, 3, 33))
    y = rng.standard_normal((5, 33))
    lt = np.array([_lit_trapz(x, y[j]) for j in range(5)])
    assert np.array_equal(lt, dg.host_trapz(x, y))
    assert np.array_equal(lt, np.asarray(jax.jit(dg.trapz_jax)(x, y)))


@pytest.mark.parametrize("layout", ["lxyes", "xyles", "yxles", "lexys", "lyxes", "yxels"])
def test_layout_orders(layout):
    """le / flux_vs_e orders from the layout string agree with an explicit iglo decode."""
    sizes = {"it": 3, "ik": 2, "il": 4, "ie": 5, "is": 2}
    seq = list(_decode(layout)(sizes))
    le_e_outer, fve = dg.layout_orders(layout)
    first = [(ie, il) for _, it, ik, il, ie, is_ in seq if (it, ik, is_) == (1, 1, 1)]
    if le_e_outer:
        expect = [(ie, il) for ie in range(5) for il in range(4)]
    else:
        expect = [(ie, il) for il in range(4) for ie in range(5)]
    assert first == expect
    got = [dict(il=il, it=it, ik=ik) for _, it, ik, il, ie, is_ in seq if (ie, is_) == (2, 1)]
    import itertools
    rng = {"il": range(4), "it": range(3), "ik": range(2)}
    expect = [dict(zip(fve, v)) for v in itertools.product(*(rng[o] for o in fve))]
    assert got == expect
    # iglo_of (S, E, L, Y, X) cross-check
    m = np.zeros((2, 5, 4, 2, 3), int)
    for iglo, it, ik, il, ie, is_ in seq:
        m[is_, ie, il, ik, it] = iglo
    assert dg.layout_orders(layout, m) == (le_e_outer, fve)


@pytest.mark.parametrize("cfg", SYNTH_CONFIGS, ids=lambda c: f"{c['layout']}-{c['fapar']}{c['fbpar']}-nl{int(c['nonlin'])}")
def test_write_block_matches_literal_fortran(cfg):
    """Device+host write block == literal Fortran transcription (strict mode: bitwise), and
    the round trips change gnew/g exactly as the literal g_adjust loops do."""
    ex = synthetic_export(**cfg)
    params = dg.DiagParams.from_export(ex)
    data = dg.DiagData.from_export(ex)
    gnew, g, fold, fnew = random_state(ex)
    state = dg.DiagState.initial(params)
    state.omega_average[:] = (0.3 + 0.1j) * np.arange(1, params.ntheta0 + 1)[None, :] \
        * np.linspace(-1, 1, params.naky)[:, None]
    code_dt = 0.05
    gn, go, rec = _run_block(ex, params, data, gnew, g, fold, fnew, istep=4, code_dt=code_dt,
                             user_time=0.2, state=state, strict=True)
    lit, lgnew, lg = _lit_record(ex, params, gnew, g, fold, fnew, state.omega_average, code_dt)
    assert np.array_equal(gn, lgnew), _maxrel(lgnew, gn)
    assert np.array_equal(go, lg), _maxrel(lg, go)
    assert not np.array_equal(gn, gnew)        # the round trip is not the identity
    bad = {}
    for name, ref in lit.items():
        if name not in rec:
            bad[name] = "missing"
            continue
        if not np.array_equal(np.asarray(ref), np.asarray(rec[name])):
            bad[name] = _maxrel(ref, rec[name])
    assert not bad, bad


def test_write_block_nonstrict_close():
    """Without strict mode (XLA may fuse a*b+c) the result stays within round-off."""
    ex = synthetic_export(**SYNTH_CONFIGS[1])
    params = dg.DiagParams.from_export(ex)
    data = dg.DiagData.from_export(ex)
    gnew, g, fold, fnew = random_state(ex)
    st1, st2 = dg.DiagState.initial(params), dg.DiagState.initial(params)
    a = _run_block(ex, params, data, gnew, g, fold, fnew, istep=2, code_dt=0.05, user_time=0.1,
                   state=st1, strict=True)
    b = _run_block(ex, params, data, gnew, g, fold, fnew, istep=2, code_dt=0.05, user_time=0.1,
                   state=st2, strict=False)
    assert _maxrel(a[0], b[0]) <= 1e-15
    for name, ref in a[2].items():
        assert _maxrel(ref, b[2][name]) <= 1e-12, name


def _lit_omega_step(hist, avg, num, den, code_dt, istep, params):
    """diagnostics_omega.f90:64-136 literally (scalar loops; Smith division, libm clog)."""
    navg = params.navg
    fac = 1.0 if params.omegatol == 0 else 1000 / abs(params.omegatol)
    k = istep % navg
    Y, X = num.shape
    for ik in range(Y):
        for it in range(X):
            n, d = complex(num[ik, it]), complex(den[ik, it])
            if math.hypot(n.real, n.imag) < dg.TINY * fac or math.hypot(d.real, d.imag) < dg.TINY * fac:
                hist[k, ik, it] = 0.0
                continue
            ar, ai, br, bi = n.real, n.imag, d.real, d.imag
            if abs(br) < abs(bi):
                r = br / bi
                dd = br * r + bi
                q = complex((ar * r + ai) / dd, (ai * r - ar) / dd)
            else:
                r = bi / br
                dd = bi * r + br
                q = complex((ai * r + ar) / dd, (ai - ar * r) / dd)
            lq = complex(np.log(np.complex128(q)))
            hist[k, ik, it] = complex(-lq.imag / code_dt, lq.real / code_dt)
    out = np.zeros((Y, X), complex)
    for ik in range(Y):
        for it in range(X):
            sr = si = 0.0
            for m in range(navg):
                sr = sr + hist[m, ik, it].real / navg
                si = si + hist[m, ik, it].imag / navg
            out[ik, it] = complex(sr, si)
    conv = unstable = False
    if istep > navg:
        ok = True
        for ik in range(Y):
            for it in range(X):
                s = 0.0
                for m in range(navg):
                    dd = out[ik, it] - hist[m, ik, it]
                    a = math.hypot(dd.real, dd.imag)
                    s = s + (a * a) / navg
                if not math.sqrt(s) <= min(math.hypot(out[ik, it].real, out[ik, it].imag), 1.0) * params.omegatol:
                    ok = False
        conv = ok
        unstable = any(math.hypot(v.real, v.imag) * code_dt > params.omegatinst for v in out.ravel())
    return out, conv, unstable


@pytest.mark.parametrize("omegatol", [-1e-3, 1e-3, 0.0, 10.0])
def test_omega_history_and_exit(omegatol):
    """calculate_omega (host, bitwise vs literal; device within 1e-13) over many steps,
    including the convergence and instability exits."""
    ex = synthetic_export(flags=dict(omegatol=omegatol, navg=4, omegatinst=50.0))
    params = dg.DiagParams.from_export(ex)
    rng = np.random.default_rng(3)
    Y, X = params.naky, params.ntheta0
    w0 = 0.4 - 0.05j
    state = dg.DiagState.initial(params)
    lhist = np.zeros((params.navg, Y, X), complex)
    jhist = jnp.zeros((params.navg, Y, X), jnp.complex128)
    den = rng.standard_normal((Y, X)) + 1j * rng.standard_normal((Y, X))
    den[0, 0] = 0.0                                   # below tiny: omega = 0 for this mode
    dt = 0.05
    calc = jax.jit(dg.calculate_omega_jax, static_argnames=("params",))
    for istep in range(1, 14):
        num = den * np.exp(-1j * w0 * dt) * (1 + 1e-6 * rng.standard_normal((Y, X)))
        num[0, 0] = 0.0
        state.exit = False
        dg.calculate_omega(state, den, num, dt, istep, params)
        lavg, lconv, lunst = _lit_omega_step(lhist, None, num, den, dt, istep, params)
        assert np.array_equal(state.omegahist, lhist)
        assert np.array_equal(state.omega_average, lavg)
        expect_exit = (params.exit_when_converged if lconv else False) or lunst
        assert state.exit == expect_exit, (istep, lconv, lunst)
        jhist, javg, jconv, junst = calc(jhist, jnp.asarray(num), jnp.asarray(den), dt, istep, params)
        assert _maxrel(lavg, javg) <= 1e-13
        assert bool(jconv) == lconv and bool(junst) == lunst
        den = num
    if omegatol > 0:
        assert state.exit                              # converged (10.0) or not (1e-3)?
    if omegatol <= 0:
        assert not state.exit


def test_run_diagnostics_schedule(tmp_path):
    """Write steps (mod(istep,nwrite)==0, exit, force, never twice), nout bookkeeping and the
    appended records, on a synthetic file."""
    ex = synthetic_export(flags=dict(nwrite=3, navg=2))
    params = dg.DiagParams.from_export(ex)
    data = dg.DiagData.from_export(ex)
    path = outnc.OutNC.create_minimal(tmp_path / "run.out.nc", akx=data.akx, aky=data.aky,
                                      theta=data.theta, nspec=params.nspec, negrid=params.negrid)
    diag = dg.Diagnostics(data, params, path, record1_written=False)
    gnew, g, fold, fnew = random_state(ex)
    fo, fn = Fields(*map(jnp.asarray, fold)), Fields(*map(jnp.asarray, fnew))
    gnew, g = jnp.asarray(gnew), jnp.asarray(g)
    written = []
    t = 0.0
    g0 = np.asarray(gnew)
    for istep in range(0, 8):
        nout = diag.state.nout
        gnew, g, exit_ = diag.step(istep, gnew, g, fo, fn, 0.1, t)
        if diag.state.nout != nout:
            written.append(istep)
        if istep == 0:     # the first round trip is not the identity
            assert not np.array_equal(g0, np.asarray(gnew))
        t += 0.1
    gnew, g = diag.finalize(7, gnew, g, fo, fn, 0.1, t - 0.1)
    assert written == [0, 3, 6]
    assert diag.state.istep_last == 7
    with netCDF4.Dataset(path) as ds:
        tt = np.asarray(ds["t"][:])
        assert len(tt) == 4                           # istep 0, 3, 6 and the forced 7
        np.testing.assert_allclose(tt, [0.0, 0.3, 0.6, 0.7], rtol=0, atol=1e-15)
        assert "epar" in ds.variables and "ntot" in ds.variables
        assert ds["phi_corr_2pi"].shape[0] == 4


def test_get_epar_literal():
    """get_epar (diagnostics_final_routines.f90:201-230) against a scalar transcription."""
    ex = synthetic_export(fphi=1.0, fapar=1.0, fbpar=0.0)
    params = dg.DiagParams.from_export(ex)
    data = dg.DiagData.from_export(ex)
    _, _, fold, fnew = random_state(ex)
    host = dg.HostData.from_data(data)
    got = dg.get_epar(dg._np_fields(fold), dg._np_fields(fnew), 0.05, params, host)
    gp, dth, tun = host.gradpar, host.delthet, host.tunits
    Y, X, T = fold.phi.shape
    for ik in range(Y):
        for it in range(X):
            for ig in range(T - 1):
                p0, p1 = fold.phi[ik, it, ig], fold.phi[ik, it, ig + 1]
                n0, n1 = fnew.phi[ik, it, ig], fnew.phi[ik, it, ig + 1]
                a0, a1 = fold.apar[ik, it, ig], fold.apar[ik, it, ig + 1]
                m0, m1 = fnew.apar[ik, it, ig], fnew.apar[ik, it, ig + 1]
                pr = ((p1.real - p0.real) + n1.real) - n0.real
                pi = ((p1.imag - p0.imag) + n1.imag) - n0.imag
                pmr, pmi = (0.5 * pr) * params.fphi, (0.5 * pi) * params.fphi
                ar = ((m1.real + m0.real) - a1.real) - a0.real
                ai = ((m1.imag + m0.imag) - a1.imag) - a0.imag
                amr, ami = (0.5 * ar) * params.fapar, (0.5 * ai) * params.fapar
                g2 = abs(gp[ig]) + abs(gp[ig + 1])
                er = (-(((pmr / dth[ig]) * 0.5) * g2)) - ((amr / tun[ik]) / 0.05)
                ei = (-(((pmi / dth[ig]) * 0.5) * g2)) - ((ami / tun[ik]) / 0.05)
                assert got[ik, it, ig] == complex(er, ei)
    assert np.all(got[..., -1] == 0)


def test_finalize_final_routines(tmp_path):
    """finalize: forced write only for a new istep, final fields/epar/moments/antot written,
    one more round trip on gnew (diagnostics_final_routines.f90:330)."""
    ex = synthetic_export(fphi=1.0, fapar=1.0, fbpar=1.0, flags=dict(nwrite=5, write_final_antot=True))
    params = dg.DiagParams.from_export(ex)
    data = dg.DiagData.from_export(ex)
    path = outnc.OutNC.create_minimal(tmp_path / "f.out.nc", akx=data.akx, aky=data.aky,
                                      theta=data.theta, nspec=params.nspec, negrid=params.negrid)
    seen = {}

    def getan(gnew):
        seen["gnew"] = np.asarray(gnew)
        f = np.ones(ex.field_shape, complex)
        return f, 2 * f, 3 * f
    diag = dg.Diagnostics(data, params, path, record1_written=False, getan=getan)
    gnew, g, fold, fnew = random_state(ex)
    fo, fn = Fields(*map(jnp.asarray, fold)), Fields(*map(jnp.asarray, fnew))
    gnew, g = diag.step(0, jnp.asarray(gnew), jnp.asarray(g), fo, fn, 0.05, 0.0)[:2]
    gnew, g = diag.step(5, gnew, g, fo, fn, 0.05, 0.25)[:2]
    before = np.asarray(gnew)
    gnew, g = diag.finalize(5, gnew, g, fo, fn, 0.05, 0.25)
    # same istep as the last write: no extra record; the final moments round trip ran
    assert np.array_equal(seen["gnew"], np.asarray(gnew))
    expect = dg._gadj(dg._gadj(jnp.asarray(before), fn.phi, fn.bpar, dg.FROM_G, to_jax(data), params),
                      fn.phi, fn.bpar, dg.TO_G, to_jax(data), params)
    # eager (op by op) vs jitted (XLA may contract a*b+c): equal up to round-off
    assert _maxrel(np.asarray(expect), np.asarray(gnew)) <= 1e-15
    with netCDF4.Dataset(path) as ds:
        assert len(ds["t"]) == 2
        for name in ("epar", "antot", "antota", "antotp", "ntot", "phi", "apar", "bpar"):
            assert name in ds.variables, name
        np.testing.assert_array_equal(ds["antota"][..., 0], 2.0)
        assert ds["antot"].long_name == outnc.ANTOT_LONG_NAME["antot"]


def test_gnostics_flags_from_namelist():
    """Smart defaults (nonlinear: write_fluxes, write_fluxes_by_mode) and check_parameters
    (write_line forces write_fields and write_omega) when /diagnostics is not exported."""
    text = "&gs2_diagnostics_knobs\n write_fields = .false.\n write_omega = .false.\n nwrite = 7\n/\n"
    f = dg._gnostics_from_namelist(text, nonlin=True)
    assert f["write_fluxes"] and f["write_fluxes_by_mode"] and f["nwrite"] == 7
    assert f["write_fields"] and f["write_omega"]          # forced by write_line = T (default)
    f = dg._gnostics_from_namelist(text.replace("/\n", " write_line = .false.\n/\n"), nonlin=False)
    assert not f["write_fluxes"] and not f["write_fields"] and not f["write_omega"]
    with pytest.raises(NotImplementedError):
        dg.DiagParams.from_export(synthetic_export(flags=dict(write_verr=True)))


# ============================================================================ 2. format


def _out_cases():
    return [c for c in all_cases() if (REFDATA / c / f"{c}.out.nc").exists()] or ["<no-refdata>"]


def _grids_from_outnc(path):
    with netCDF4.Dataset(path) as ds:
        ds.set_auto_mask(False)
        v = lambda n: np.asarray(ds.variables[n][...])   # noqa: E731
        return dict(theta=v("theta"), akx=v("kx"), aky=v("ky"), al=v("lambda"), energy=v("energy"),
                    kperp2=v("kperp2"), theta0=v("theta0"), bmag=v("bmag"), gradpar=v("gradpar"),
                    jacob=v("jacob"), grho=v("grho"), bpol=v("bpol"), rplot=v("rplot"),
                    gds21=v("gds21"), gds22=v("gds22"), shat=float(v("shat")), qval=float(v("qval")),
                    charge=v("charge"), mass=v("mass"), dens=v("dens"), temp=v("temp"),
                    tprim=v("tprim"), fprim=v("fprim"))


def _export_for_format(case):
    """The real export if present, else a synthetic one with GS2's grids and the case flags."""
    if have(case):
        return export_of(case)
    import f90nml
    nml = f90nml.read(CASES_DIR / f"{case}.in")
    kn = nml.get("knobs", {})
    fphi, fapar, fbpar = (float(kn.get(k, d)) for k, d in (("fphi", 1.0), ("fapar", 0.0), ("fbpar", -1.0)))
    if float(kn.get("beta", 0.0)) == 0.0:
        fapar = fbpar = 0.0
    nonlin = nml.get("nonlinear_terms_knobs", {}).get("nonlinear_mode", "off") == "on"
    flags = dg._gnostics_from_namelist((CASES_DIR / f"{case}.in").read_text(), nonlin)
    grids = _grids_from_outnc(REFDATA / case / f"{case}.out.nc")
    lay = nml.get("layouts_knobs", {}).get("layout", "lxyes")
    return synthetic_export(layout=lay, fphi=fphi, fapar=max(fapar, 0.0), fbpar=max(fbpar, 0.0),
                            nonlin=nonlin, flags=flags, grids=grids)


@pytest.mark.parametrize("case", _out_cases())
def test_outnc_format(case, tmp_path):
    """Run the write sequence for every record of GS2's file into a copy holding only the
    static part, so the writer creates every diagnostic variable itself; then compare
    variables, dims and attributes with GS2's ``<case>.out.nc``."""
    if case == "<no-refdata>":
        pytest.skip("no GS2 out.nc in REFDATA")
    ref = REFDATA / case / f"{case}.out.nc"
    ex = _export_for_format(case)
    params = dg.DiagParams.from_export(ex)
    data = dg.DiagData.from_export(ex)
    excl = [s.name for s in outnc.record_specs(params) + outnc.final_specs(params) if s.name != "t"]
    path = outnc.copy_records(ref, tmp_path / f"{case}.out.nc", nrec=0, exclude=excl)
    with netCDF4.Dataset(ref) as ds:
        nrec = len(ds.dimensions["t"])
    diag = dg.Diagnostics(data, params, path, record1_written=False)
    gnew, g, fold, fnew = random_state(ex)
    fo, fn = Fields(*map(jnp.asarray, fold)), Fields(*map(jnp.asarray, fnew))
    gnew, g = jnp.asarray(gnew), jnp.asarray(g)
    for k in range(nrec):
        istep = k * params.nwrite
        gnew, g, _ = diag.step(istep, gnew, g, fo, fn, 0.05, 0.1 * k)
    diag.finalize((nrec - 1) * params.nwrite, gnew, g, fo, fn, 0.05, 0.1 * (nrec - 1))
    problems = outnc.check_format(path, ref)
    assert not problems, problems
    with netCDF4.Dataset(path) as ds:
        assert len(ds.dimensions["t"]) == nrec


# ============================================================================ 3. reference data


def _cases_with(kind):
    return [c for c in all_cases() if have(c, kind) and have(c)] or ["<no-refdata>"]


def _skip_placeholder(case):
    if case == "<no-refdata>":
        pytest.skip("no reference data (run gpu/scripts/make_refdata.sh)")


@lru_cache(maxsize=4)
def _setup(case):
    ex = export_of(case)
    return ex, dg.DiagParams.from_export(ex), dg.DiagData.from_export(ex)


def _fields(grp, dims, new=True):
    return fields_from_group(grp, dims, new=new)


def _dims(ex):
    from gs2gpu.types import Dims
    return Dims.from_export(ex)


def _load(tr, tag, istep):
    try:
        return tr.load(tr.record(tag, istep))
    except ExportError:
        return None


def _trace_inputs(tr, istep, dims):
    """State at diagnostic time of ``istep`` from the trace (spec 08 section 1.2), or None.

    gnew and the new fields: ``step_end_pre_diag`` (after update_time, gs2_main.fpp:771), or
    ``advance_end`` of the kept attempt (nothing changes them in between) with the time block
    of ``step_end`` (same code_dt and user_time: check_time_step runs after it). Old fields:
    ``advance_end`` (phi, apar, bpar). g: ``after_exb_shear`` gnew when flow shear is on
    (g = gnew right after the remap, fields_implicit.f90:335-337), else ``advance_begin`` gnew.
    """
    pre = _load(tr, "step_end_pre_diag", istep)
    end = _load(tr, "advance_end", istep)
    if pre is None:
        post = _load(tr, "step_end", istep)
        if end is None or post is None:
            return None
        pre, timeblock = end, post
    else:
        timeblock = pre
    fnew = _fields(pre, dims, new=True)
    src = _load(tr, "after_exb_shear", istep) or _load(tr, "advance_begin", istep)
    g = pre.get("g") if pre is not end else None
    if g is None:
        if src is None:
            return None
        g = src["gnew"]
    if end is not None and "phi" in end:
        fold = _fields(end, dims, new=False)
    elif src is not None:
        fold = _fields(src, dims, new=True)
    else:
        return None
    return dict(gnew=np.asarray(pre["gnew"]), g=np.asarray(g), fold=fold, fnew=fnew,
                code_dt=float(timeblock["code_dt"]), user_time=float(timeblock["user_time"]))


def _write_steps(tr, params):
    return [s for s in tr.steps if s > 0 and s % params.nwrite == 0]


@pytest.mark.parametrize("case", _cases_with("trace"))
def test_roundtrip_trace(case):
    """step_end_pre_diag -> (round trips) -> step_end, for gnew and g (spec 08 section 3.1.8)."""
    _skip_placeholder(case)
    ex, params, data = _setup(case)
    tr = trace_of(case)
    dims = _dims(ex)
    steps = _write_steps(tr, params)
    if not steps:
        pytest.skip(f"{case}: no write step (mod(istep,{params.nwrite})==0) in the trace")
    errs = {}
    for istep in steps:
        inp = _trace_inputs(tr, istep, dims)
        if inp is None:
            continue
        post = tr.load(tr.record("step_end", istep))
        state = dg.DiagState.initial(params)
        gn, go, _ = _run_block(ex, params, data, inp["gnew"], inp["g"], inp["fold"], inp["fnew"],
                               istep=istep, code_dt=inp["code_dt"], user_time=inp["user_time"],
                               state=state)
        errs[istep] = (assert_close(post["gnew"], gn, RTOL_KERNEL, f"{case} istep={istep} gnew"),
                       assert_close(post["g"], go, RTOL_KERNEL, f"{case} istep={istep} g")
                       if "g" in post else None)
    if not errs:
        pytest.skip(f"{case}: trace lacks the inputs of every write step")
    print(f"{case} round trip max errors (gnew, g): {errs}")


@pytest.mark.parametrize("case", _cases_with("export"))
def test_roundtrip_t0(case):
    """/state_t0_prediag -> istep-0 write block -> /state_t0 (gnew and g; spec 10 G-05)."""
    _skip_placeholder(case)
    ex, params, data = _setup(case)
    if "state_t0_prediag" not in ex.root.groups:
        pytest.skip("no /state_t0_prediag")
    pre, post = ex["/state_t0_prediag"], ex["/state_t0"]
    dims = _dims(ex)
    state = dg.DiagState.initial(params, record1_written=False)
    gn, go, _ = _run_block(ex, params, data, pre["gnew"], pre["g"], _fields(pre, dims, False),
                           _fields(pre, dims, True), istep=0, code_dt=float(pre["code_dt"]),
                           user_time=float(pre["user_time"]), state=state)
    assert_close(post["gnew"], gn, RTOL_KERNEL, f"{case} t0 gnew")
    assert_close(post["g"], go, RTOL_KERNEL, f"{case} t0 g")


@pytest.mark.parametrize("case", _cases_with("export"))
def test_field_line_average_weight(case):
    """The recomputed Klein trapezoid weight equals the exported one bitwise (spec 08 2.1)."""
    _skip_placeholder(case)
    ex, params, data = _setup(case)
    tg = ex["/theta_grid"]
    if "field_line_average_weight" not in tg:
        pytest.skip("weight not exported")
    w = 1.0 / dg.host_trapz(np.asarray(tg["theta"]), np.asarray(tg["jacob"]))
    assert w == float(tg["field_line_average_weight"])




#: Variables whose error is reported, not asserted at 1e-12 (they divide by code_dt).
OMEGA_VARS = ("omega", "omegaavg", "omega_average")

#: Flux families: name -> (raw device integrals, field, species normalisation).
_FLUX_SRC = {
    "es_heat_flux": (("qflux1",), "phi", "DT"), "es_heat_flux_par": (("qflux2",), "phi", "DT"),
    "es_heat_flux_perp": (("qflux3",), "phi", "DT"), "es_part_flux": (("pflux",), "phi", "dens"),
    "es_mom_flux": (("vflux_par", "vflux_perp"), "phi", "MOM"),
    "es_mom_flux_par": (("vflux_par",), "phi", "MOM"), "es_mom_flux_perp": (("vflux_perp",), "phi", "MOMP"),
    "apar_heat_flux": (("qmflux1",), "apar", "DT"), "apar_heat_flux_par": (("qmflux2",), "apar", "DT"),
    "apar_heat_flux_perp": (("qmflux3",), "apar", "DT"), "apar_part_flux": (("pmflux",), "apar", "dens"),
    "apar_mom_flux": (("vmflux",), "apar", "MOMP"),
    "bpar_heat_flux": (("qbflux1",), "bpar", "DT"), "bpar_heat_flux_par": (("qbflux2",), "bpar", "DT"),
    "bpar_heat_flux_perp": (("qbflux3",), "bpar", "DT"), "bpar_part_flux": (("pbflux",), "bpar", "dens"),
    "bpar_mom_flux": (("vbflux",), "bpar", "MOMP"),
}


def _magnitudes(dev, fold, fnew, code_dt, params, host):
    """Uncancelled magnitudes of the flux-type outputs (spec 08 P10).

    Fluxes are theta integrals of ``Im(total*conjg(fld))`` that cancel to round-off when the
    flux vanishes physically (e.g. the electrostatic particle/heat fluxes of the single-mode
    s-alpha cases are ~1e-25 at every record). The comparison floor is the same integral of
    ``|total|*|fld|`` (summed over modes and species), so round-off of the terms is accepted
    and nothing else.
    """
    h = host
    norm = {"DT": np.abs(h.dens * h.temp), "dens": np.abs(h.dens),
            "MOM": np.abs(h.dens * np.sqrt(h.mass * h.temp)), "MOMP": np.abs((h.dens * h.mass) * h.stm)}
    fields = dg._np_fields(fnew)
    mags = {}
    for name, (keys, fld, nrm) in _FLUX_SRC.items():
        if not all(k in dev for k in keys):
            continue
        f = np.abs(getattr(fields, fld))
        m = 0.0
        for k in keys:
            y = np.abs(dev[k]) * f[None] * h.jacob
            per = (dg.host_trapz(h.theta, y) * h.flux_wgt) * np.abs(h.aky)[None, :, None] * 0.5
            m = m + float(np.sum(per * norm[nrm][:, None, None]))
        mags[name] = m
    if "exch_total" in dev:
        fo = dg._np_fields(fold)
        y = np.abs(dev["exch_total"]) * np.abs(fields.phi - fo.phi)[None] \
            + np.abs(fields.phi + fo.phi)[None] * np.abs(dev["exch_total2"])
        mags["es_energy_exchange"] = float(np.sum(dg.host_flavg(y, h) / abs(code_dt)
                                                  * np.abs(h.dens * h.z)[:, None, None]))
    out = {}
    for name, m in mags.items():
        for v in (name, f"total_{name}_by_ky", f"total_{name}_by_kx", f"total_{name}",
                  f"{name}_by_mode", f"total_{name}_by_mode"):
            out[v] = m
    for tot, fam in (("heat", "heat_flux"), ("mom", "mom_flux"), ("part", "part_flux")):
        m = sum(mags.get(f"{p}_{fam}", 0.0) for p in ("es", "apar", "bpar"))
        for v in ({"heat": ("heat_flux_tot", "hflux_tot"), "mom": ("mom_flux_tot", "vflux_tot"),
                   "part": ("part_flux_tot", "zflux_tot")}[tot]):
            out[v] = m
    for fld, (vname, _ln) in outnc.FLUX_VS_E.items():
        prefix = {"phi": "es", "apar": "apar", "bpar": "bpar"}[fld]
        if f"{prefix}_part_flux" in mags:
            out[vname] = mags[f"{prefix}_part_flux"]
    # theta averages of the moments: <M> cancels for some modes (e.g. the zonal ntot of an
    # adiabatic-electron run vanishes by quasineutrality); floor = max over modes of <|M|>
    if "mom_density" in dev:
        for name, M in dg.normalised_moments(dev, h).items():
            m = float(np.max(dg.host_flavg(np.abs(M), h)))
            out[f"{name}_flxsurf_avg"] = m
            out[f"{name}_by_mode"] = m
    return out


def _ref_record(path, user_time, guess):
    with netCDF4.Dataset(path) as ds:
        t = np.asarray(ds["t"][:])
    hit = np.nonzero(t == user_time)[0]
    return int(hit[0]) if hit.size else (guess if guess < len(t) else None)


def _read_ref(ds, name, idx, spec):
    v = ds.variables[name]
    data = np.asarray(v[idx] if spec.time else v[...])
    if spec.is_complex:
        data = data[..., 0] + 1j * data[..., 1]
    return data


def _compare_record(params, rec, path, idx, *, static=False, scales=None):
    """``{name: (max-abs error / max(max|ref|, scale), bitwise)}`` for every variable of ``rec``
    present in ``path`` (record ``idx`` for record variables; static ones only if ``static``)."""
    specs = {s.name: s for s in outnc.record_specs(params) + outnc.final_specs(params)}
    scales = scales or {}
    out = {}
    with netCDF4.Dataset(path) as ds:
        ds.set_auto_mask(False)
        for name, val in rec.items():
            spec = specs.get(name)
            if spec is None or name not in ds.variables or (not spec.time and not static):
                continue
            ref = _read_ref(ds, name, idx, spec)
            got = np.asarray(val)
            if name == "epar":
                ref, got = ref[..., :-1], got[..., :-1]
            denom = max(float(np.max(np.abs(ref))) if ref.size else 0.0, scales.get(name, 0.0))
            diff = float(np.max(np.abs(got - ref))) if ref.size else 0.0
            out[name] = (diff / denom if denom > 0 else diff, bool(np.array_equal(ref, got)))
    return out


def _omega_history(tr, istep, params, dims, device=False):
    """omegahist/omega_average at ``istep`` replayed from every traced step 1..istep (or None);
    ``device=True`` replays with :func:`calculate_omega_jax` and returns ``(hist, avg)``."""
    state = dg.DiagState.initial(params)
    ig = params.ig
    hist = jnp.zeros_like(jnp.asarray(state.omegahist))
    avg = None
    calc = jax.jit(dg.calculate_omega_jax, static_argnames=("params",))
    for k in range(1, istep + 1):
        inp = _trace_inputs(tr, k, dims) if k in tr.steps else None
        if inp is None:
            return None
        num, den = dg._ig_sum(inp["fnew"], ig), dg._ig_sum(inp["fold"], ig)
        if device:
            if params.omega_active:
                hist, avg, _, _ = calc(hist, jnp.asarray(num), jnp.asarray(den), inp["code_dt"], k,
                                       params)
        else:
            dg.calculate_omega(state, den, num, inp["code_dt"], k, params)
    if device:
        return np.asarray(hist), (np.asarray(avg) if avg is not None else np.zeros_like(state.omega_average))
    return state


def _check(errs, strict, what):
    if strict:
        bad = {k: e for k, (e, bw) in errs.items() if not bw}
    else:
        bad = {k: e for k, (e, bw) in errs.items() if e > RTOL_KERNEL and k not in OMEGA_VARS}
    worst = sorted(((k, e) for k, (e, _) in errs.items()), key=lambda kv: -kv[1])[:5]
    nbit = sum(bw for _, bw in errs.values())
    print(f"{what}: {nbit}/{len(errs)} bitwise; worst {[(k, f'{e:.2e}') for k, e in worst]}")
    assert not bad, bad


MODES = pytest.mark.parametrize("strict", [True, False], ids=["strict", "fast"])


@MODES
@pytest.mark.parametrize("case", _cases_with("export"))
def test_values_record1(case, strict):
    """The istep-0 write block from ``/state_t0_prediag`` vs record 1 of GS2's ``<case>.out.nc``
    and, for the static variables (final fields, moments, eigenfunctions at t = 0), vs the
    exporter's ``export/<case>.out.nc``. Strict mode: every value bitwise; fast mode: 1e-12
    relative to max(|ref|, uncancelled magnitude)."""
    _skip_placeholder(case)
    ex, params, data = _setup(case)
    if "state_t0_prediag" not in ex.root.groups:
        pytest.skip("no /state_t0_prediag")
    pre, dims = ex["/state_t0_prediag"], _dims(ex)
    state = dg.DiagState.initial(params, record1_written=False)
    fold, fnew = _fields(pre, dims, False), _fields(pre, dims, True)
    code_dt = float(pre["code_dt"])
    _, _, rec, dev = _run_block(ex, params, data, pre["gnew"], pre["g"], fold, fnew, istep=0,
                                code_dt=code_dt, user_time=float(pre["user_time"]), state=state,
                                strict=strict, return_dev=True)
    scales = {} if strict else _magnitudes(dev, fold, fnew, code_dt, params, dg.HostData.from_data(data))
    errs = _compare_record(params, rec, out_of(case), 0, scales=scales)
    exp_out = REFDATA / case / "export" / f"{case}.out.nc"
    if exp_out.exists() and not params.write_final_moments:
        static = {k: v for k, v in rec.items() if k in {s.name for s in outnc.record_specs(params)
                                                       if not s.time}}
        errs.update({f"export:{k}": v for k, v in
                     _compare_record(params, static, exp_out, 0, static=True).items()})
    _check(errs, strict, f"{case} record 1 ({'strict' if strict else 'fast'})")


@MODES
@pytest.mark.parametrize("case", _cases_with("trace"))
def test_values_trace(case, strict):
    """Every variable of the traced write steps vs the same record of GS2's ``<case>.out.nc``,
    with the trace states as inputs (tests the diagnostics alone) and the omega history
    replayed from every traced step."""
    _skip_placeholder(case)
    ex, params, data = _setup(case)
    tr, ref, dims = trace_of(case), out_of(case), _dims(ex)
    steps = _write_steps(tr, params)
    if not steps:
        pytest.skip(f"{case}: no write step (mod(istep,{params.nwrite})==0) in the trace")
    errs = {}
    for istep in steps:
        inp = _trace_inputs(tr, istep, dims)
        state = _omega_history(tr, istep, params, dims)
        if inp is None or state is None:
            continue
        _, _, rec, dev = _run_block(ex, params, data, inp["gnew"], inp["g"], inp["fold"], inp["fnew"],
                                    istep=istep, code_dt=inp["code_dt"], user_time=inp["user_time"],
                                    state=state, strict=strict, return_dev=True)
        idx = _ref_record(ref, inp["user_time"], istep // params.nwrite)
        if idx is None:
            continue
        scales = {} if strict else _magnitudes(dev, inp["fold"], inp["fnew"], inp["code_dt"],
                                               params, dg.HostData.from_data(data))
        errs.update({f"{istep}:{k}": v
                     for k, v in _compare_record(params, rec, ref, idx, scales=scales).items()})
    if not errs:
        pytest.skip(f"{case}: trace lacks the inputs of every write step")
    _check(errs, strict, f"{case} traced write steps ({'strict' if strict else 'fast'})")


@pytest.mark.parametrize("case", _cases_with("trace"))
def test_omega_device_variant(case):
    """``calculate_omega_jax`` (XLA complex log/abs) over the replayed history vs GS2's omega and
    omega_average: reported; must stay within 1e-12 relative to the largest |omega|."""
    _skip_placeholder(case)
    ex, params, data = _setup(case)
    if not params.omega_active:
        pytest.skip("nonlinear run: omega is not computed (zeros written)")
    tr, ref, dims = trace_of(case), out_of(case), _dims(ex)
    steps = _write_steps(tr, params)
    report = {}
    for istep in steps:
        res = _omega_history(tr, istep, params, dims, device=True)
        inp = _trace_inputs(tr, istep, dims)
        if res is None or inp is None:
            continue
        hist, avg = res
        idx = _ref_record(ref, inp["user_time"], istep // params.nwrite)
        with netCDF4.Dataset(ref) as ds:
            r_om = ds["omega"][idx][..., 0] + 1j * ds["omega"][idx][..., 1]
            r_av = ds["omega_average"][idx][..., 0] + 1j * ds["omega_average"][idx][..., 1]
        wu = np.asarray(data.woutunits)[:, None]
        report[f"{istep}:omega"] = _maxrel(r_om, hist[istep % params.navg] * wu)
        report[f"{istep}:omega_average"] = _maxrel(r_av, avg * wu)
    if not report:
        pytest.skip("no traced write step")
    print(f"{case}: device omega errors {report}")
    assert max(report.values()) <= RTOL_KERNEL, report


@pytest.mark.parametrize("case", _cases_with("trace"))
def test_append_to_export_outnc(case, tmp_path, monkeypatch):
    """Production flow on the exporter's own file: open ``export/<case>.out.nc`` (static part +
    record 1) in append mode, run ``Diagnostics.step`` for every traced step from the trace
    states (nout from 2), then ``finalize``; the result must equal the Fortran traced run's
    ``trace/<case>.out.nc`` bitwise (strict mode) record by record, the round-tripped gnew must
    equal ``step_end`` at every step, and the format must match GS2's full-run file."""
    _skip_placeholder(case)
    ex, params, data = _setup(case)
    src = REFDATA / case / "export" / f"{case}.out.nc"
    fortran = REFDATA / case / "trace" / f"{case}.out.nc"
    if not src.exists() or not fortran.exists():
        pytest.skip("no export/ or trace/ out.nc")
    import shutil
    path = tmp_path / f"{case}.out.nc"
    shutil.copy(src, path)
    tr, dims = trace_of(case), _dims(ex)
    if not tr.find("step_end"):
        pytest.skip("GS2 aborted before completing a step; the exit path (zeroed final "
                    "state) is covered by test_end_to_end")
    # If GS2 changed dt at the last traced step, its final write sees the post-reset state
    # (fields recomputed from g), which this harness does not rebuild: compare the records
    # only; the final write of such a case is covered by test_end_to_end.
    reset_at_end = bool(tr.find("after_reset_time_step", istep=max(tr.steps)))
    diag = dg.Diagnostics.from_export(ex, path)
    assert diag.state.nout == 2 and diag.state.istep_last == 0
    last = None
    try:
        monkeypatch.setattr(dg, "final_moments_device", _final_moments_kernel(True))
        monkeypatch.setattr(dg, "diagnostics_device", _device(True))   # (sets strict mode)
        for istep in sorted(s for s in tr.steps if s > 0):
            inp = _trace_inputs(tr, istep, dims)
            assert inp is not None, istep
            gnew, g, exit_ = diag.step(istep, jnp.asarray(inp["gnew"]), jnp.asarray(inp["g"]),
                                       Fields(*map(jnp.asarray, inp["fold"])),
                                       Fields(*map(jnp.asarray, inp["fnew"])),
                                       inp["code_dt"], inp["user_time"])
            post = tr.load(tr.record("step_end", istep))
            assert np.array_equal(np.asarray(gnew), post["gnew"]), f"gnew istep {istep}"
            assert not exit_
            last = (istep, gnew, g, inp)
        istep, gnew, g, inp = last
        diag.finalize(istep, gnew, g, Fields(*map(jnp.asarray, inp["fold"])),
                      Fields(*map(jnp.asarray, inp["fnew"])), inp["code_dt"], inp["user_time"])
    finally:
        numerics.set_strict(False)
        diag.close()
    specs = outnc.record_specs(params) + ([] if reset_at_end else outnc.final_specs(params))
    with netCDF4.Dataset(path) as got, netCDF4.Dataset(fortran) as ref:
        assert len(got["t"]) == len(ref["t"])
        bad = {}
        for s in specs:
            if s.name not in ref.variables:
                continue
            if reset_at_end and ref.variables[s.name].dimensions[:1] != ("t",):
                continue    # overwritten by GS2's final write after the reset (see above)
            a, b = np.asarray(got[s.name][...]), np.asarray(ref[s.name][...])
            if s.name == "epar":
                a, b = a[..., :-1, :], b[..., :-1, :]
            if a.shape != b.shape or not np.array_equal(a, b):
                bad[s.name] = _maxrel(b, a) if a.shape == b.shape else "shape"
        n = len([s for s in specs if s.name in ref.variables
                 and not (reset_at_end and ref.variables[s.name].dimensions[:1] != ("t",))])
        print(f"{case}: {n - len(bad)}/{n} variables bitwise equal to trace/{case}.out.nc over "
              f"{len(ref['t'])} records")
        assert not bad, bad
    problems = outnc.check_format(path, out_of(case))
    assert not problems, problems
