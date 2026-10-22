"""Tests for B7: hyperviscosity (gs2gpu.hyper), ExB flow shear (gs2gpu.exb) and
zonal-flow removal (gs2gpu.zonal). Specs 07 sections 3.2, 3.3, 3.12; 01 sections 3.7-3.8.

Part A runs without reference data: host bookkeeping against hand-computed values,
and the kernels against loop-literal translations of the Fortran on synthetic
arrays. Part B compares with the GS2 trace (skips when gpu/refdata is missing):

* hyper_diff: input ``after_invert_rhs`` gnew of the same pass, phinew valid in
  that pass; output ``after_hyper`` (gnew, hypervisc_filter); every case.
* exb_shear: input ``advance_begin``; output ``after_exb_shear`` (gnew, fields,
  kx_shift/theta0_shift, jump, ikx_indexed); every case (identity when inactive).
* remove_zonal_flows: input ``after_field_update``; output ``after_remove_zonal_flows``.
"""

from __future__ import annotations

import dataclasses
import math
from pathlib import Path

import jax
import jax.numpy as jnp
import numpy as np
import pytest
from jax import lax

from conftest import RTOL_KERNEL, assert_close, cases_where, export_of, trace_of
from gs2gpu import exb, hyper, zonal
from gs2gpu.io.export import Export, Group
from gs2gpu.numerics import max_abs_rel
from gs2gpu.types import Dims, Fields, TimeBlock, fields_from_group, to_jax

EPS = float(np.finfo(np.float64).eps)
RNG = np.random.default_rng(20261021)


# =============================================================================
# helpers
# =============================================================================
def _crandn(*shape):
    return RNG.standard_normal(shape) + 1j * RNG.standard_normal(shape)


def _box_ikx(nt):
    """Box-grid ikx: 0, 1, .., (nt-1)/2, -(nt-1)/2, .., -1 (kt_grids_box.f90:382-389)."""
    half = (nt - 1) // 2
    return np.array(list(range(0, half + 1)) + list(range(-half, 0)))


def _time(code_dt, code_dt_old, code_time=0.0):
    return TimeBlock(code_dt=code_dt, user_dt=code_dt, code_dt_prev1=-1.0,
                     code_dt_prev2=-1.0, code_dt_old=code_dt_old, code_time=code_time,
                     user_time=code_time, code_dt_cfl=1e8, user_dt_cfl=1e8,
                     code_dt_min=1e-6, code_dt_max=code_dt)


def _group(path, arrays=None, groups=None, attrs=None):
    return Group(name=path.rsplit("/", 1)[-1] or "/", path=path, arrays=dict(arrays or {}),
                 groups=dict(groups or {}), attrs=dict(attrs or {}))


def fake_export(*, naky=4, ntheta0=7, ntheta=12, nperiod=1, box=True, shat=0.8,
                g_exb=0.1, g_exbfac=1.0, start_time=-1.0, start_timestep=-1,
                hyper_vars=None, has_apar=False, has_bpar=False, remove_zonal=False,
                tunits=None, nspec=1, negrid=2, nlambda=3):
    """A minimal in-memory Export with the groups B7 reads (names as in spec 09)."""
    ntgrid = ntheta // 2 + (nperiod - 1) * ntheta
    nth = 2 * ntgrid + 1
    theta = np.linspace(-np.pi * (2 * nperiod - 1), np.pi * (2 * nperiod - 1), nth)
    jacob = 1.0 + 0.3 * np.cos(theta) + 0.01 * RNG.standard_normal(nth)
    if box:
        ikx = _box_ikx(ntheta0)
        akx = 0.3 * ikx
        aky = 0.2 * np.arange(naky)
        theta0 = np.zeros((naky, ntheta0))
    else:
        ikx = np.arange(ntheta0)
        akx = np.linspace(-0.4, 0.4, ntheta0)
        aky = np.linspace(0.3, 0.5, naky)
        theta0 = np.tile(np.linspace(-0.3, 0.3, ntheta0), (naky, 1))
    gds2 = 1.0 + (shat * theta) ** 2
    gds21 = -shat * shat * theta
    gds22 = np.full(nth, shat * shat)
    kperp2 = (aky[:, None, None] ** 2 * gds2[None, None, :]
              + 2 * aky[:, None, None] * akx[None, :, None] * gds21 / shat
              + akx[None, :, None] ** 2 * gds22 / shat ** 2)
    trapz = np.sum(np.diff(theta) * (jacob[:-1] + jacob[1:])) * 0.5
    state = {"exb_first": True, "code_dt": 0.05}
    if abs(g_exb * g_exbfac) > EPS:
        state["kx_shift" if (box or shat == 0.0) else "theta0_shift"] = np.zeros(naky)
    hv = {"hyper_option_switch": 1, "nexp": 2, "D_hypervisc": -10.0, "D_hyperres": -10.0,
          "D_hyper": -10.0, "omega_osc": 0.4, "const_amp": False, "isotropic_shear": True,
          "damp_zonal_only": False, "gridnorm": True, "hyper_on": False,
          "D_hyper3D": -10.0, "P_hyper3D": 4.0, "ky_cut": -10.0, "kx_cut": -10.0,
          "isotropic_model": True, "D_res": np.zeros((naky, ntheta0))}
    hv.update(hyper_vars or {})
    groups = {
        "kt_grids": _group("/kt_grids", {"aky": aky, "akx": akx, "ikx": ikx.astype(np.int32),
                                         "theta0": theta0, "kperp2": kperp2, "box": box}),
        "theta_grid": _group("/theta_grid", {
            "theta": theta, "jacob": jacob, "shat": shat, "ntheta": ntheta, "ntgrid": ntgrid,
            "gds2": gds2, "gds21": gds21, "gds22": gds22,
            "field_line_average_weight": 1.0 / trapz}),
        "hyper": _group("/hyper", hv),
        "dist_fn": _group("/dist_fn", {"g_exb": g_exb, "g_exbfac": g_exbfac,
                                       "g_exb_start_time": start_time,
                                       "g_exb_start_timestep": start_timestep,
                                       "hyper_in_initialisation": True}),
        "run_parameters": _group("/run_parameters", {"has_phi": True, "has_apar": has_apar,
                                                     "has_bpar": has_bpar,
                                                     "wstar_units": False}),
        "gs2_time": _group("/gs2_time", {"tunits": np.ones(naky) if tunits is None
                                         else np.asarray(tunits, float)}),
        "state_t0": _group("/state_t0", state),
        "fields": _group("/fields", {"remove_zonal_flows_switch": remove_zonal}),
    }
    root = _group("/", groups=groups)
    root.path = "/"
    dims = {"theta": nth, "kx": ntheta0, "ky": naky, "lambda": nlambda, "egrid": negrid,
            "species": nspec}
    return Export(path=Path("fake.gs2export.nc"), root=root, dims=dims, attrs={})


# ---------------------------------------------------------------- Fortran-literal refs
def ref_powi(x: float, n: int) -> float:
    """libgcc ``__powidf2`` on a Python float."""
    u = abs(n)
    y = x if u % 2 else 1.0
    while True:
        u >>= 1
        if not u:
            break
        x = x * x
        if u % 2:
            y = y * x
    return 1.0 / y if n < 0 else y


def ref_klein_sum(values) -> float:
    """summation.f90:71-103 on Python floats."""
    s = c1 = c2 = 0.0
    for x in values:
        x = float(x)
        t = s + x
        lc1 = (s - t) + x if abs(s) >= abs(x) else (x - t) + s
        s = t
        t = c1 + lc1
        lc2 = (c1 - t) + lc1 if abs(c1) >= abs(lc1) else (lc1 - t) + c1
        c1 = t
        c2 = c2 + lc2
    return (s + c1) + c2


def ref_hyper_diff(g, phi, code_dt, akx, aky, p: hyper.HyperParams, filter_init):
    """hyper.f90:620-786, loop-literal (NumPy vectors along theta only)."""
    g = np.array(g, copy=True)
    naky, ntheta0, nth = phi.shape
    if not p.hyper_on:
        return g, filter_init
    filt = np.array(filter_init, copy=True)
    fac = np.full(naky, 0.5)
    fac[0] = 1.0

    def sq(z):
        return z.real * z.real + z.imag * z.imag

    if p.switch in (hyper.HYPER_BOTH, hyper.HYPER_RES, hyper.HYPER_VISC):
        if p.D_hypervisc < 0.0:
            return g, filter_init
        if p.isotropic_shear:
            if p.const_amp:
                s = np.ones(nth)
            else:
                s = np.zeros(nth)
                for ik in range(naky):
                    for it in range(ntheta0):
                        k2 = akx[it] * akx[it] + aky[ik] * aky[ik]
                        s = s + ((sq(phi[ik, it]) * (k2 * k2)) * fac[ik])
                s = s ** 0.5
            for ik in range(naky):
                for it in range(ntheta0):
                    if p.damp_zonal_only and not aky[ik] == EPS:
                        continue
                    q = aky[ik] * aky[ik] + akx[it] * akx[it]
                    filt[ik, it] = np.exp(-((p.D_hypervisc * code_dt)
                                            * ((s * ref_powi(q, p.nexp)) / p.akperp4_max)))
        else:
            snz = np.zeros(nth)
            for ik in range(1, naky):
                for it in range(ntheta0):
                    k2 = akx[it] * akx[it] + aky[ik] * aky[ik]
                    snz = snz + ((sq(phi[ik, it]) * (k2 * k2)) * fac[ik])
            snz = 0.5 * (-p.omega_osc + np.sqrt(p.omega_osc * p.omega_osc + 2 * snz))
            sz = np.zeros(nth)
            for it in range(ntheta0):
                k2 = akx[it] * akx[it] + aky[0] * aky[0]
                sz = sz + ((sq(phi[0, it]) * (k2 * k2)) * fac[0])
            sz = 0.5 * (-p.omega_osc + np.sqrt(p.omega_osc * p.omega_osc + 2 * sz))
            sznz = np.zeros(nth)
            for ik in range(1, naky):
                for it in range(ntheta0):
                    a2 = aky[ik] * aky[ik]
                    sznz = sznz + ((sq(phi[ik, it]) * (a2 * a2)) * fac[ik])
            sznz = 0.5 * (-p.omega_osc + np.sqrt(p.omega_osc * p.omega_osc + 2 * sznz))
            for ik in range(naky):
                for it in range(ntheta0):
                    x2 = akx[it] * akx[it]
                    x4 = x2 * x2
                    if aky[ik] == 0.0:
                        arg = (p.D_hypervisc * code_dt) * ((sznz * x4) / p.akx4_max)
                    else:
                        q = aky[ik] * aky[ik] + akx[it] * akx[it]
                        arg = (p.D_hypervisc * code_dt) * (
                            ((snz * (q * q)) / p.akperp4_max)
                            + ((((sz * x4) / p.akx4_max) * aky[ik]) / p.aky_max))
                    filt[ik, it] = np.exp(-arg)
    elif p.switch != hyper.HYPER_SIMPLE3D:
        return g, filter_init
    for ik in range(naky):
        for it in range(ntheta0):
            if (p.switch != hyper.HYPER_SIMPLE3D and p.isotropic_shear and p.damp_zonal_only
                    and not aky[ik] == EPS):
                continue
            f = filt[ik, it]
            sl = g[:, :, :, ik, it]
            g[:, :, :, ik, it] = (sl.real * f) + 1j * (sl.imag * f)
    return g, filt


def ref_exb_kx(g, fields, jump, ikx_indexed, has):
    """dist_fn.fpp:4490-4599, in-place loops (ikx_indexed is Fortran 1-based)."""
    g = np.array(g, copy=True)
    xs = [np.array(f, copy=True) for f in fields]
    naky, nt = xs[0].shape[:2]

    def ikxi(it):
        return int(ikx_indexed[it - 1]) - 1

    for ik in range(naky, 1, -1):
        J = int(jump[ik - 1])
        k = ik - 1
        targets = [x for x, flag in zip(xs, has) if flag]
        if J < 0:
            for x in targets:
                for it in range(1, nt + J + 1):
                    x[k, ikxi(it)] = x[k, ikxi(it - J)]
                for it in range(nt + J + 1, nt + 1):
                    x[k, ikxi(it)] = 0.0
            for it in range(1, nt + J + 1):
                g[:, :, :, k, ikxi(it)] = g[:, :, :, k, ikxi(it - J)]
            for it in range(nt + J + 1, nt + 1):
                g[:, :, :, k, ikxi(it)] = 0.0
        if J > 0:
            for x in targets:
                for it in range(nt, J, -1):
                    x[k, ikxi(it)] = x[k, ikxi(it - J)]
                for it in range(J, 0, -1):
                    x[k, ikxi(it)] = 0.0
            for it in range(nt, J, -1):
                g[:, :, :, k, ikxi(it)] = g[:, :, :, k, ikxi(it - J)]
            for it in range(J, 0, -1):
                g[:, :, :, k, ikxi(it)] = 0.0
    return g, xs


def ref_exb_theta0(g, fields, jump, ntgrid, ntheta, has):
    """dist_fn.fpp:4341-4489, in-place loops with temp/temp2 and the all-ky zeroing."""
    g = np.array(g, copy=True)
    xs = [np.array(f, copy=True) for f in fields]
    naky, nt, nth = xs[0].shape
    ileft = -ntgrid + ntheta
    iright = ntgrid - ntheta

    def r(ig):  # Fortran theta index -> NumPy index
        return ig + ntgrid

    for ik in range(1, naky + 1):
        j = int(jump[ik - 1])
        if j == 0:
            continue
        assert abs(j) < nt
        k = ik - 1
        iit = nt + 1 - abs(j)
        iib = abs(j)
        for x, flag in zip(xs, has):
            if not flag:
                continue
            if j < 0:
                temp2 = x[k, :iib, :].copy()                       # phi(:, :iib, ik)
                for i in range(1, iit):
                    x[k, i - 1, :] = x[k, i - j - 1, :]
                n = max(0, r(iright) + 1)
                x[k, iit - 1:, r(ileft):] = temp2[:, :n]           # phi(ileft:, iit:, ik)
                x[k, iit - 1:, :max(0, r(ileft))] = 0.0            # phi(:ileft-1, iit:, ik)
            else:
                temp2 = x[k, iit - 1:, :].copy()                   # phi(:, iit:, ik)
                for i in range(nt, iib, -1):
                    x[k, i - 1, :] = x[k, i - j - 1, :]
                n = max(0, r(iright) + 1)
                x[k, :iib, :n] = temp2[:, r(ileft):]               # phi(:iright, :iib, ik)
                x[:, :iib, n:] = 0.0                                # ALL ky (4370)
        # distribution function, all (is, ie, il), both signs
        if j < 0:
            temp = g[:, :, :, k, :iib].copy()                      # (S,E,L,iib,2,T)
            for it in range(1, iit):
                g[:, :, :, k, it - 1] = g[:, :, :, k, it - j - 1]
            for it in range(iit, nt + 1):
                n = max(0, r(iright) + 1)
                g[:, :, :, k, it - 1, :, r(ileft):] = temp[:, :, :, it - iit, :, :n]
                g[:, :, :, k, it - 1, :, :max(0, r(ileft))] = 0.0
        else:
            temp = np.stack([g[:, :, :, k, iit + it - 2] for it in range(1, j + 1)], axis=3)
            for it in range(nt, iib, -1):
                g[:, :, :, k, it - 1] = g[:, :, :, k, it - j - 1]
            for it in range(1, iib + 1):
                n = max(0, r(iright) + 1)
                g[:, :, :, k, it - 1, :, :n] = temp[:, :, :, it - 1, :, r(ileft):]
                g[:, :, :, k, it - 1, :, n:] = 0.0
    return g, xs


def ref_remove_zonal(phinew, theta, jacob, weight):
    """fields_implicit.f90:374-425 + theta_grid.f90:665-670 + integration.f90, literal."""
    out = np.array(phinew, copy=True)
    d = theta[1:] - theta[:-1]
    for it in range(phinew.shape[1]):
        q = phinew[0, it]
        qr, qi = q.real * jacob, q.imag * jacob
        yr = d * (qr[:-1] + qr[1:])
        yi = d * (qi[:-1] + qi[1:])
        ir, ii = ref_klein_sum(yr) * 0.5, ref_klein_sum(yi) * 0.5
        ar, ai = ir * weight, ii * weight
        out[0, it] = (phinew[0, it].real - ar) + 1j * (phinew[0, it].imag - ai)
    return out


# =============================================================================
# Part A: no reference data
# =============================================================================
def test_fortran_nint():
    cases = {0.5: 1, -0.5: -1, 1.5: 2, -1.5: -2, 2.5: 3, -2.5: -3, 0.0: 0, -0.0: 0,
             0.49999999999999994: 0, -0.49999999999999994: 0, 0.5000000000000001: 1,
             -0.5000000000000001: -1, 1.4999999999999998: 1, 2.0: 2, -3.0: -3,
             4503599627370495.5: 4503599627370496, -1e-300: 0, 7.3: 7, -7.7: -8}
    for x, want in cases.items():
        assert exb.fortran_nint(x) == want, x
    # jnp.round / Python round are half-to-even and would differ on ties
    assert round(2.5) == 2 and exb.fortran_nint(2.5) == 3


def test_powi_matches_libgcc_and_integer_pow():
    x = RNG.uniform(-3, 3, 200)
    for n in range(-3, 10):
        got = hyper.powi(x, n)
        want = np.array([ref_powi(float(v), n) for v in x])
        assert np.array_equal(got, want), n
        if n >= 0:
            jx = np.asarray(lax.integer_pow(jnp.asarray(x), n))
            assert np.array_equal(jx, want), n


def test_resolve_hyper_parameters():
    R = hyper.resolve_hyper_parameters
    base = dict(D_hypervisc=0.1, D_hyperres=0.2, D_hyper=-10.0, D_hyper3D=-10.0,
                P_hyper3D=4.0, nexp=2, isotropic_shear=True)
    e = R(hyper.HYPER_NONE, **base)
    assert e == dict(hyper_on=False, D_hypervisc=-10.0, D_hyperres=-10.0, nexp=2)
    e = R(hyper.HYPER_VISC, **base)
    assert e == dict(hyper_on=True, D_hypervisc=0.1, D_hyperres=-10.0, nexp=2)
    e = R(hyper.HYPER_VISC, **{**base, "D_hypervisc": -1.0})
    assert not e["hyper_on"]
    e = R(hyper.HYPER_RES, **base)
    assert e == dict(hyper_on=True, D_hypervisc=-10.0, D_hyperres=0.2, nexp=2)
    e = R(hyper.HYPER_RES, **{**base, "D_hyperres": -1.0})
    assert not e["hyper_on"]
    e = R(hyper.HYPER_BOTH, **{**base, "D_hyper": 0.0})        # D_hyper == 0 counts as set
    assert e == dict(hyper_on=True, D_hypervisc=0.0, D_hyperres=0.0, nexp=2)
    e = R(hyper.HYPER_BOTH, **{**base, "D_hypervisc": -1.0, "D_hyperres": -2.0})
    assert not e["hyper_on"]
    e = R(hyper.HYPER_BOTH, **{**base, "D_hypervisc": -1.0})
    assert e["hyper_on"] and e["D_hypervisc"] == -1.0
    assert not R(hyper.HYPER_SIMPLE3D, **{**base, "D_hyper3D": 0.1, "P_hyper3D": 3.0})["hyper_on"]
    assert not R(hyper.HYPER_SIMPLE3D, **{**base, "D_hyper3D": -0.1})["hyper_on"]
    assert R(hyper.HYPER_SIMPLE3D, **{**base, "D_hyper3D": 0.1})["hyper_on"]
    assert R(hyper.HYPER_VISC, **{**base, "nexp": 4, "isotropic_shear": False})["nexp"] == 2
    assert R(hyper.HYPER_VISC, **{**base, "nexp": 4})["nexp"] == 4
    # idempotent: applying it to the effective (exported) values changes nothing
    for sw in (1, 2, 3, 4, 5):
        for _ in range(20):
            raw = dict(D_hypervisc=float(RNG.choice([-10, -1, 0, 0.1])),
                       D_hyperres=float(RNG.choice([-10, -1, 0, 0.2])),
                       D_hyper=float(RNG.choice([-10, 0, 0.3])),
                       D_hyper3D=float(RNG.choice([-1, 0.1])),
                       P_hyper3D=float(RNG.choice([2.0, 4.0, 6.0])),
                       nexp=int(RNG.choice([1, 2, 3])), isotropic_shear=bool(RNG.integers(2)))
            once = R(sw, **raw)
            twice = R(sw, **{**raw, "D_hypervisc": once["D_hypervisc"],
                             "D_hyperres": once["D_hyperres"], "nexp": once["nexp"]})
            assert once == twice
    with pytest.raises(NotImplementedError):
        R(7, **base)


def test_normalisers_hand_values():
    akx = np.array([0.0, 0.5, 1.0, -1.0, -0.5])     # akx(ntheta0/2+1) = akx(3) = 1.0
    aky = np.array([0.0, 0.25, 0.5])
    akx4, kymax, aky4, akp4 = hyper.normalisers(akx, aky, 2, True)
    assert (akx4, kymax, aky4, akp4) == (1.0, 0.5, 0.0625, 1.5625)   # (1+0.25)^2
    assert hyper.normalisers(akx, aky, 3, True)[0] == 1.0
    assert hyper.normalisers(akx, aky, 2, False) == (1.0, 1.0, 1.0, 1.0)
    d = hyper.d_res_formula(akx, aky, 2, 0.5, 0.1, akp4)
    assert d[1, 1] == ((0.5 * 0.1) * (0.0625 + 0.25) ** 2) / 1.5625
    assert np.all(hyper.d_res_formula(akx, aky, 2, -10.0, 0.1, akp4) == 0.0)


def _hyper_setup(params_kw, *, naky=4, ntheta0=7, nth=13, box=True):
    akx = 0.3 * _box_ikx(ntheta0) if box else np.linspace(-0.4, 0.5, ntheta0)
    aky = 0.2 * np.arange(naky) if box else np.linspace(0.3, 0.5, naky)
    switch = params_kw.get("switch", hyper.HYPER_VISC)
    nexp = params_kw.get("nexp", 2)
    norms = hyper.normalisers(akx, aky, nexp, True)
    base = dict(switch=switch, hyper_on=True, D_hypervisc=0.1, D_hyperres=-10.0, nexp=nexp,
                const_amp=False, isotropic_shear=True, damp_zonal_only=False, omega_osc=0.4,
                gridnorm=True, akx4_max=norms[0], aky_max=norms[1], aky4_max=norms[2],
                akperp4_max=norms[3], D_hyper3D=-10.0, P_hyper3D=4.0)
    base.update(params_kw)
    p = hyper.HyperParams(**base)
    k2 = akx[None, :] * akx[None, :] + aky[:, None] * aky[:, None]
    q = aky[:, None] * aky[:, None] + akx[None, :] * akx[None, :]
    fac = np.full(naky, 0.5)
    fac[0] = 1.0
    filt = RNG.uniform(0.5, 1.0, (naky, ntheta0, nth)) if switch == hyper.HYPER_SIMPLE3D \
        else np.ones((naky, ntheta0, nth))
    d = hyper.HyperData(k4=k2 * k2, qn=hyper.powi(q, nexp), x4=(akx * akx) * (akx * akx),
                        a4=(aky * aky) * (aky * aky), fac=fac, aky=aky, aky_zero=aky == 0.0,
                        zonal_skip=np.full(naky, p.damp_zonal_only) & ~(aky == EPS),
                        D_res=np.zeros((naky, ntheta0)), filter_static=filt)
    return akx, aky, p, d


HYPER_VARIANTS = {
    "iso": {},
    "iso_nexp3": {"nexp": 3},
    "iso_const_amp": {"const_amp": True},
    "iso_damp_zonal_only": {"damp_zonal_only": True},
    "aniso": {"isotropic_shear": False},
    "aniso_const_amp": {"isotropic_shear": False, "const_amp": True},
    "both": {"switch": hyper.HYPER_BOTH, "D_hyperres": 0.1},
    "simple3d": {"switch": hyper.HYPER_SIMPLE3D, "D_hyper3D": 0.1},
    "simple3d_off": {"switch": hyper.HYPER_SIMPLE3D, "hyper_on": False},
    "res_only": {"switch": hyper.HYPER_RES, "D_hypervisc": -10.0, "D_hyperres": 0.1},
    "none": {"switch": hyper.HYPER_NONE, "hyper_on": False},
}


@pytest.mark.parametrize("variant", sorted(HYPER_VARIANTS))
@pytest.mark.parametrize("box", [True, False])
def test_hyper_diff_vs_fortran_loops(variant, box):
    akx, aky, p, d = _hyper_setup(HYPER_VARIANTS[variant], box=box)
    phi = _crandn(*d.filter_static.shape) * 0.3
    g = _crandn(2, 2, 3, *d.filter_static.shape[:2], 2, d.filter_static.shape[2])
    code_dt = 0.05
    want_g, want_f = ref_hyper_diff(g, phi, code_dt, akx, aky, p, d.filter_static)
    fn = jax.jit(hyper.hyper_diff, static_argnames=("params",))
    got_g, got_f = fn(jnp.asarray(g), jnp.asarray(phi), code_dt, to_jax(d), p)
    assert max_abs_rel(want_g, got_g) <= 1e-14, variant
    assert max_abs_rel(want_f, got_f) <= 1e-14, variant
    if not p.applies_filter or (p.damp_zonal_only and p.isotropic_shear):
        assert np.array_equal(np.asarray(got_g), g)            # exact identity
    # traced code_dt gives the same numbers
    got_g2, _ = fn(jnp.asarray(g), jnp.asarray(phi), jnp.float64(code_dt), to_jax(d), p)
    assert np.array_equal(np.asarray(got_g2), np.asarray(got_g))


def test_hyper_box_zonal_mode_is_identity():
    """(it, ik) = (1, 1) in box: q = 0 -> F = exp(-0) = 1 exactly (07 section 3.2.7)."""
    for kw in ({}, {"isotropic_shear": False}):
        akx, aky, p, d = _hyper_setup(kw)
        phi = _crandn(*d.filter_static.shape)
        g = _crandn(1, 1, 2, *d.filter_static.shape[:2], 2, d.filter_static.shape[2])
        got_g, got_f = hyper.hyper_diff(jnp.asarray(g), jnp.asarray(phi), 0.05, to_jax(d), p)
        assert np.all(np.asarray(got_f)[0, 0] == 1.0)
        assert np.array_equal(np.asarray(got_g)[..., 0, 0, :, :], g[..., 0, 0, :, :])


def test_sequential_sum_is_fortran_order():
    terms = RNG.standard_normal((300, 5)) * np.logspace(-8, 8, 300)[:, None]
    got = np.asarray(hyper._sequential_sum(jnp.asarray(terms), True))
    want = np.zeros(5)
    for t in terms:
        want = want + t
    assert np.array_equal(got, want)
    assert np.asarray(hyper._sequential_sum(jnp.zeros((0, 3)), True)).tolist() == [0.0] * 3


def test_hyper_from_fake_export():
    ex = fake_export(hyper_vars={"hyper_option_switch": 2, "D_hypervisc": 0.1,
                                 "hyper_on": True, "isotropic_shear": False})
    p = hyper.HyperParams.from_export(ex)
    akx, aky = ex["/kt_grids/akx"], ex["/kt_grids/aky"]
    assert (p.akx4_max, p.aky_max, p.aky4_max, p.akperp4_max) == \
        hyper.normalisers(akx, aky, 2, True)
    assert p.switch == hyper.HYPER_VISC and p.hyper_on and p.applies_filter
    d = hyper.HyperData.from_export(ex, p)
    assert d.filter_static.shape == ex.field_shape and np.all(d.filter_static == 1.0)
    # inconsistent exported hyper_on is caught
    bad = fake_export(hyper_vars={"hyper_option_switch": 2, "D_hypervisc": -1.0,
                                  "hyper_on": True})
    with pytest.raises(ValueError):
        hyper.HyperParams.from_export(bad)
    # simple3D without an exported filter falls back to the init_hyper formula
    ex3 = fake_export(hyper_vars={"hyper_option_switch": 5, "D_hyper3D": 0.1, "hyper_on": True})
    p3 = hyper.HyperParams.from_export(ex3)
    d3 = hyper.HyperData.from_export(ex3, p3)
    tg = ex3["/theta_grid"]
    want = hyper.simple3d_filter_formula(ex3["/kt_grids/kperp2"], tg["gds2"], tg["gds21"],
                                         tg["gds22"], tg["shat"], akx, aky, -10.0, -10.0,
                                         True, 0.1, 4.0, 0.05)
    assert np.array_equal(d3.filter_static, want[0])
    assert want[3] == aky[-1] and want[4] == akx[(akx.size + 1) // 2 - 1]
    assert np.all(want[0][0, 0] == 1.0)                       # kperp2 = 0 at (1,1)


# ------------------------------------------------------------------------------- ExB host
def test_exb_first_step_half_dt_and_jump_bookkeeping():
    """kx branch with exactly representable numbers: dkx = 0.25, aky = (0, .5, 1),
    g_exb = 0.5, dt = 0.125 -> full-step shift aky/16, first step half of it."""
    ex = fake_export(naky=3, ntheta0=5, g_exb=0.5)
    kt = ex["/kt_grids"]
    kt.arrays["akx"] = 0.25 * _box_ikx(5)
    kt.arrays["aky"] = np.array([0.0, 0.5, 1.0])
    params = exb.ExbParams.from_export(ex)
    state = exb.ExbState.from_export(ex, params)
    assert state.exb_first and state.istep_last == 0 and params.use_kx
    dt = 0.125
    # step 1: code_dt_old = 0 -> gdt = dt/2
    state, jump = exb.exb_update(state, 1, _time(dt, 0.0), params)
    assert jump is None and not state.exb_first
    assert state.kx_shift.tolist() == [0.0, -0.015625, -0.03125]
    assert state.ikx_indexed.tolist() == [4, 5, 1, 2, 3] and state.itmin == 4
    assert state.dkx == 0.25 and state.istep_last == 1
    # retry of the same istep (immediate reset): no second shift
    again, jump = exb.exb_update(state, 1, _time(dt / 2, 0.0), params)
    assert jump is None and again is state
    # step 2: -0.03125 - 0.0625 = -0.09375 -> nint(-0.375) = 0
    state, jump = exb.exb_update(state, 2, _time(dt, dt), params)
    assert jump is None and state.kx_shift[2] == -0.09375
    # step 3: -0.15625 -> nint(-0.625) = -1 -> kx_shift = 0.09375
    state, jump = exb.exb_update(state, 3, _time(dt, dt), params)
    assert jump.tolist() == [0, 0, -1]
    assert state.kx_shift.tolist() == [0.0, -0.078125, 0.09375]
    # a dt change: gdt = 0.5*(new + old)
    state, jump = exb.exb_update(state, 4, _time(dt / 2, dt), params)
    assert state.kx_shift[2] == 0.09375 - 0.75 * 0.0625


def test_exb_nint_tie_rounds_away_from_zero():
    ex = fake_export(naky=2, ntheta0=5, g_exb=0.5)
    ex["/kt_grids"].arrays["akx"] = 0.25 * _box_ikx(5)
    ex["/kt_grids"].arrays["aky"] = np.array([0.0, 1.0])
    params = exb.ExbParams.from_export(ex)
    state = exb.ExbState.from_export(ex, params)
    dt = 0.125
    state, jump = exb.exb_update(state, 1, _time(dt, dt), params)   # -0.0625 -> -0.25
    assert jump is None
    state, jump = exb.exb_update(state, 2, _time(dt, dt), params)   # -0.125 -> -0.5 -> -1
    assert jump.tolist() == [0, -1] and state.kx_shift[1] == 0.125


def test_exb_start_guards():
    ex = fake_export(naky=2, ntheta0=5, g_exb=0.5, start_timestep=3)
    params = exb.ExbParams.from_export(ex)
    state = exb.ExbState.from_export(ex, params)
    for istep in (1, 2):
        state, jump = exb.exb_update(state, istep, _time(0.1, 0.1 * (istep > 1)), params)
        assert jump is None and state.istep_last == 0 and np.all(state.kx_shift == 0.0)
        assert not state.exb_first                       # the setup block runs anyway
    state, _ = exb.exb_update(state, 3, _time(0.1, 0.1), params)
    assert state.istep_last == 3 and state.kx_shift[1] == -(((0.2 * 0.5) * 1.0) * 0.1)
    ex = fake_export(naky=2, ntheta0=5, g_exb=0.5, start_time=1.0)
    params = exb.ExbParams.from_export(ex)
    state = exb.ExbState.from_export(ex, params)
    state, jump = exb.exb_update(state, 1, _time(0.1, 0.0, code_time=0.5), params)
    assert jump is None and state.istep_last == 0
    state, _ = exb.exb_update(state, 2, _time(0.1, 0.1, code_time=1.0), params)
    assert state.istep_last == 2
    # inactive: nothing happens, exb_first stays true
    ex = fake_export(naky=2, ntheta0=5, g_exb=0.0)
    params = exb.ExbParams.from_export(ex)
    state = exb.ExbState.from_export(ex, params)
    out, jump = exb.exb_update(state, 1, _time(0.1, 0.0), params)
    assert jump is None and out is state and out.exb_first


def test_exb_theta0_bookkeeping_and_abort():
    ex = fake_export(naky=2, ntheta0=3, box=False, shat=0.8, g_exb=0.05,
                     tunits=[1.0, 2.0])
    params = exb.ExbParams.from_export(ex)
    assert not params.use_kx
    state = exb.ExbState.from_export(ex, params)
    dt = 0.05
    state, jump = exb.exb_update(state, 1, _time(dt, 0.0), params)
    th = ex["/kt_grids/theta0"]
    assert state.dtheta0 == th[0, 1] - th[0, 0]
    want = [0.0 - ((((0.05 * 1.0) * 0.025) / 0.8) * t) for t in (1.0, 2.0)]
    assert state.theta0_shift.tolist() == want
    # force a large jump: |jump| >= ntheta0 aborts (mp_abort, 4347-4351)
    big = exb.ExbState(theta0_shift=np.array([0.0, -0.3 * 2.6]), exb_first=True)
    with pytest.raises(RuntimeError):
        exb.exb_update(big, 5, _time(dt, dt), params)
    # g_exb /= 0 with ntheta0 == 1 is gated in both branches (critique B.6) ...
    for box in (True, False):
        with pytest.raises(NotImplementedError):
            exb.ExbParams.from_export(fake_export(naky=2, ntheta0=1, box=box, g_exb=0.1))
    # ... but the Fortran theta0 fallback dtheta0 = 2*pi (4264) is still mirrored
    p1 = dataclasses.replace(params, ntheta0=1, theta0_k1=(0.0,))
    s1, _ = exb.exb_update(exb.ExbState(theta0_shift=np.zeros(2)), 1, _time(dt, 0.0), p1)
    assert s1.dtheta0 == 2.0 * math.pi


# ----------------------------------------------------------------------------- ExB device
KX_JUMPS = [[0, 1, 0, 0], [0, -1, 0, 0], [3, -2, 2, -1], [0, 0, 7, -7], [0, 5, -3, 1],
            [-1, 0, 0, 0], [0, 0, 0, 0]]


@pytest.mark.parametrize("jump", KX_JUMPS)
@pytest.mark.parametrize("has", [(True, False, False), (True, True, True)])
def test_exb_kx_shift_vs_fortran_loops(jump, has):
    ex = fake_export(naky=4, ntheta0=7, has_apar=has[1], has_bpar=has[2])
    params = exb.ExbParams.from_export(ex)
    data = exb.ExbData.from_params(params)
    ikx_indexed, _ = exb.ikx_indexed_of(params.ikx)
    fields = [_crandn(4, 7, 13) for _ in range(3)]
    g = _crandn(1, 2, 3, 4, 7, 2, 13)
    want_g, want_f = ref_exb_kx(g, fields, jump, ikx_indexed, has)
    fn = jax.jit(exb.exb_shift, static_argnames=("params",))
    got_g, got_f = fn(jnp.asarray(g), Fields(*map(jnp.asarray, fields)),
                      jnp.asarray(jump, jnp.int32), to_jax(data), params)
    assert np.array_equal(np.asarray(got_g), want_g)
    for a, b in zip(got_f, want_f):
        assert np.array_equal(np.asarray(a), b)
    assert np.array_equal(np.asarray(got_g)[:, :, :, 0], g[:, :, :, 0])   # ik = 1 never shifted


THETA0_CASES = [  # (ntheta, nperiod, ntheta0, jumps)
    (12, 1, 3, [1, 0]), (12, 1, 3, [-1, -1]), (12, 1, 3, [2, -2]),
    (8, 2, 5, [1, 0, -2, 2]), (8, 2, 5, [0, 3, 1, -4]), (8, 2, 5, [-1, -1, -1, -1]),
    (8, 2, 5, [2, 2, 2, 2]), (8, 3, 4, [-3, 1, 0, 3]), (8, 2, 5, [0, 0, 0, 0]),
]


@pytest.mark.parametrize("ntheta,nperiod,ntheta0,jump", THETA0_CASES)
@pytest.mark.parametrize("has", [(True, False, False), (True, True, True)])
def test_exb_theta0_shift_vs_fortran_loops(ntheta, nperiod, ntheta0, jump, has):
    naky = len(jump)
    ex = fake_export(naky=naky, ntheta0=ntheta0, ntheta=ntheta, nperiod=nperiod, box=False,
                     has_apar=has[1], has_bpar=has[2], g_exb=0.05)
    params = exb.ExbParams.from_export(ex)
    data = exb.ExbData.from_params(params)
    nth = ex.nth
    fields = [_crandn(naky, ntheta0, nth) for _ in range(3)]
    g = _crandn(1, 2, 2, naky, ntheta0, 2, nth)
    want_g, want_f = ref_exb_theta0(g, fields, jump, ex.ntgrid, ntheta, has)
    fn = jax.jit(exb.exb_shift, static_argnames=("params",))
    got_g, got_f = fn(jnp.asarray(g), Fields(*map(jnp.asarray, fields)),
                      jnp.asarray(jump, jnp.int32), to_jax(data), params)
    assert np.array_equal(np.asarray(got_g), want_g)
    for a, b in zip(got_f, want_f):
        assert np.array_equal(np.asarray(a), b)


def test_exb_shear_end_to_end_synthetic():
    """exb_shear = exb_update + exb_shift; repeated steps against the literal refs."""
    ex = fake_export(naky=3, ntheta0=5, g_exb=4.0)
    params = exb.ExbParams.from_export(ex)
    data = to_jax(exb.ExbData.from_params(params))
    state = exb.ExbState.from_export(ex, params)
    g = _crandn(1, 1, 2, 3, 5, 2, 13)
    f = [_crandn(3, 5, 13), np.zeros((3, 5, 13), complex), np.zeros((3, 5, 13), complex)]
    gj, fj = jnp.asarray(g), Fields(*map(jnp.asarray, f))
    shifted = 0
    for istep in range(1, 12):
        t = _time(0.05, 0.0 if istep == 1 else 0.05)
        new_state, jump = exb.exb_update(state, istep, t, params)
        gj, fj, state = exb.exb_shear(gj, fj, state, istep, t, data, params)
        assert np.array_equal(state.kx_shift, new_state.kx_shift)
        if jump is not None:
            shifted += 1
            g, f = ref_exb_kx(g, f, jump, state.ikx_indexed, (True, False, False))
        assert np.array_equal(np.asarray(gj), g)
        assert np.array_equal(np.asarray(fj.phi), f[0])
    assert shifted > 0


# --------------------------------------------------------------------------------- zonal
def test_klein_sum_bitwise():
    for n in (1, 2, 7, 24, 101):
        a = RNG.standard_normal((n, 3)) * np.logspace(-12, 12, n)[:, None]
        a[::3] *= -1
        got = np.asarray(zonal.klein_sum(jnp.asarray(a), axis=0))
        want = [ref_klein_sum(a[:, k]) for k in range(3)]
        assert got.tolist() == want, n
    tricky = np.array([1.0, 1e100, 1.0, -1e100])
    assert float(zonal.klein_sum(jnp.asarray(tricky))) == ref_klein_sum(tricky) == 2.0


@pytest.mark.parametrize("nperiod", [1, 2])
def test_remove_zonal_flows_vs_fortran(nperiod):
    ex = fake_export(naky=4, ntheta0=7, ntheta=12, nperiod=nperiod, remove_zonal=True)
    data = zonal.ZonalData.from_export(ex)
    params = zonal.ZonalParams.from_export(ex)
    assert params.remove_zonal_flows_switch
    tg = ex["/theta_grid"]
    phi = _crandn(*ex.field_shape)
    want = ref_remove_zonal(phi, tg["theta"], tg["jacob"], tg["field_line_average_weight"])
    got_eager = np.asarray(zonal.remove_zonal_flows(jnp.asarray(phi), to_jax(data), params))
    assert np.array_equal(got_eager, want)                     # op-by-op: no FMA, bitwise
    fn = jax.jit(zonal.remove_zonal_flows, static_argnames=("params",))
    got = np.asarray(fn(jnp.asarray(phi), to_jax(data), params))
    assert max_abs_rel(want, got) <= 1e-15
    assert np.array_equal(got[1:], phi[1:])                    # other ky untouched
    off = zonal.ZonalParams(False)
    assert np.array_equal(np.asarray(fn(jnp.asarray(phi), to_jax(data), off)), phi)
    # after the removal the field-line average of the zonal row is ~0
    resid = np.asarray(zonal.field_line_average(jnp.asarray(got[0]), to_jax(data)))
    assert np.max(np.abs(resid)) <= 1e-14 * np.max(np.abs(phi[0]))


def test_field_line_average_weight_formula():
    ex = fake_export(nperiod=2)
    data = zonal.ZonalData.from_export(ex)
    tg = ex["/theta_grid"]
    d = tg["theta"][1:] - tg["theta"][:-1]
    want = 1.0 / (ref_klein_sum(d * (tg["jacob"][:-1] + tg["jacob"][1:])) * 0.5)
    assert float(zonal.field_line_average_weight(data)) == want


def test_kernels_jit_trace():
    """All public kernels trace under jit with static params (no host syncs)."""
    ex = fake_export(hyper_vars={"hyper_option_switch": 2, "D_hypervisc": 0.1,
                                 "hyper_on": True}, remove_zonal=True)
    hp = hyper.HyperParams.from_export(ex)
    hd = to_jax(hyper.HyperData.from_export(ex, hp))
    g = jnp.asarray(_crandn(*ex.g_shape))
    phi = jnp.asarray(_crandn(*ex.field_shape))
    jax.jit(hyper.hyper_diff, static_argnames=("params",)).lower(g, phi, 0.05, hd, hp)
    ep = exb.ExbParams.from_export(ex)
    ed = to_jax(exb.ExbData.from_params(ep))
    jax.jit(exb.exb_shift, static_argnames=("params",)).lower(
        g, Fields(phi, phi, phi), jnp.zeros(ex.naky, jnp.int32), ed, ep)
    zp, zd = zonal.ZonalParams.from_export(ex), to_jax(zonal.ZonalData.from_export(ex))
    jax.jit(zonal.remove_zonal_flows, static_argnames=("params",)).lower(phi, zd, zp)


# =============================================================================
# Part B: GS2 reference data (skips when gpu/refdata is missing)
# =============================================================================
def _hyper_on(ex) -> bool:
    return hyper.HyperParams.from_export(ex).switch != hyper.HYPER_NONE


def _zonal_on(ex) -> bool:
    return zonal.ZonalParams.from_export(ex).remove_zonal_flows_switch


ALL_CASES = cases_where(lambda ex: True)
HYPER_CASES = cases_where(_hyper_on)
ZONAL_CASES = cases_where(_zonal_on)


def _latest(tr, rec, tags):
    """The latest record before ``rec`` (same istep/attempt) among ``tags``."""
    best = None
    for tag in tags:
        r = tr.previous(rec, tag)
        if r is not None and (best is None or r.seq > best.seq):
            best = r
    return best


def _report(name, errs):
    worst = max(errs.values()) if errs else float("nan")
    print(f"[B7] {name}: max-abs-relative errors {errs} (worst {worst:.3e})")


@pytest.mark.parametrize("case", HYPER_CASES)
def test_hyper_static_arrays(case):
    """Exported effective normalisers, frozen D_res and simple3D filter vs init_hyper."""
    ex = export_of(case)
    p = hyper.HyperParams.from_export(ex)
    d = hyper.HyperData.from_export(ex, p)
    akx, aky = ex["/kt_grids/akx"], ex["/kt_grids/aky"]
    errs = {}
    if p.switch in (hyper.HYPER_VISC, hyper.HYPER_RES, hyper.HYPER_BOTH):
        rec = hyper.normalisers(akx, aky, p.nexp, p.gridnorm)
        assert rec == (p.akx4_max, p.aky_max, p.aky4_max, p.akperp4_max)
        g = ex["/hyper"]
        dt0 = float(g.meta.get("D_res", {}).get("code_dt", ex["/state_t0/code_dt"]))
        want = hyper.d_res_formula(akx, aky, p.nexp, p.D_hyperres, dt0, p.akperp4_max)
        errs["D_res"] = assert_close(d.D_res, want, RTOL_KERNEL, f"{case} D_res") \
            if np.any(d.D_res) else float(np.max(np.abs(want)))
        assert np.all(want == 0.0) or np.any(d.D_res)
    if p.switch == hyper.HYPER_SIMPLE3D:
        g, tg = ex["/hyper"], ex["/theta_grid"]
        dt0 = float(g.meta.get("hypervisc_filter", {}).get("code_dt",
                                                           ex["/state_t0/code_dt"]))
        # effective cut-offs if exported (the formula re-applies the `.not. > 0` defaults)
        want = hyper.simple3d_filter_formula(
            ex["/kt_grids/kperp2"], tg["gds2"], tg["gds21"], tg["gds22"], float(tg["shat"]),
            akx, aky, float(hyper._lookup(g, "ky_cut", -10.0)),
            float(hyper._lookup(g, "kx_cut", -10.0)),
            hyper._as_bool(hyper._lookup(g, "isotropic_model")), p.D_hyper3D, p.P_hyper3D, dt0)
        errs["filter_s3d"] = assert_close(d.filter_static, want[0], RTOL_KERNEL,
                                          f"{case} simple3D filter")
        for name, value in zip(("kperp2_max", "kperp2_max_zonal", "ky_cut", "kx_cut"),
                               want[1:]):
            got = hyper._lookup(g, name, None)
            if got is not None:
                errs[name] = abs(float(got) - value) / abs(value)
                assert errs[name] <= RTOL_KERNEL, name
    _report(f"{case} hyper static", errs)


def _gnew_untouched_after_hyper(ex) -> bool:
    """True if no stage between hyper_diff and the end of timeadv can modify gnew
    (no collisions inside timeadv, no parity; dist_fn.fpp:4083-4091)."""
    coll = ex.root.groups.get("collisions")
    colls = coll is not None and bool(np.asarray(coll.get("colls", True)).item()) and \
        not bool(np.asarray(coll.get("split_collisions", False)).item())
    parity = bool(np.asarray(ex["/dist_fn"].get("def_parity", False)).item())
    return not colls and not parity


@pytest.mark.parametrize("case", ALL_CASES)
def test_hyper_diff_trace(case):
    """hyper_diff on after_invert_rhs vs after_hyper, every step/attempt/pass.

    The exporter records after_hyper only when hyper_on; for the other cases the
    identity is checked on every after_invert_rhs record (and against GS2's own
    after_timeadv when nothing in between can modify gnew)."""
    ex, tr = export_of(case), trace_of(case)
    p = hyper.HyperParams.from_export(ex)
    d = to_jax(hyper.HyperData.from_export(ex, p))
    fn = jax.jit(hyper.hyper_diff, static_argnames=("params",))
    errs = {"gnew": 0.0, "filter": 0.0}
    n = 0
    if not p.hyper_on:
        untouched = _gnew_untouched_after_hyper(ex)
        nchecked = 0
        for rin in tr.find("after_invert_rhs"):
            gin = tr.load(rin)
            rph = _latest(tr, rin, ("advance_begin", "after_exb_shear", "after_field_update",
                                    "after_remove_zonal_flows"))
            phi = tr.load(rph)["phinew"] if rph is not None else np.zeros(ex.field_shape, complex)
            got_g, got_f = fn(jnp.asarray(gin["gnew"]), jnp.asarray(phi),
                              float(gin["code_dt"]), d, p)
            what = f"{case} istep={rin.istep} attempt={rin.attempt} pass={rin.pass_}"
            assert np.array_equal(np.asarray(got_g), gin["gnew"]), what
            assert np.array_equal(np.asarray(got_f), np.asarray(d.filter_static)), what
            nxt = tr.next(rin)
            if untouched and nxt is not None and nxt.tag == "after_timeadv":
                assert np.array_equal(tr.load(nxt)["gnew"], gin["gnew"]), what + " (GS2)"
                nchecked += 1
            n += 1
        print(f"[B7] {case} hyper off: identity on {n} records "
              f"({nchecked} also against GS2's after_timeadv)")
        return
    for rh in tr.find("after_hyper"):
        rin = tr.previous(rh, "after_invert_rhs")
        rph = _latest(tr, rh, ("advance_begin", "after_exb_shear", "after_field_update",
                               "after_remove_zonal_flows"))
        assert rin is not None and rph is not None, f"{case}: inputs of {rh.path}"
        gh, gin, gph = tr.load(rh), tr.load(rin), tr.load(rph)
        got_g, got_f = fn(jnp.asarray(gin["gnew"]), jnp.asarray(gph["phinew"]),
                          float(gh["code_dt"]), d, p)
        what = f"{case} istep={rh.istep} attempt={rh.attempt} pass={rh.pass_}"
        errs["gnew"] = max(errs["gnew"], assert_close(gh["gnew"], got_g, RTOL_KERNEL,
                                                      what + " gnew"))
        if "hypervisc_filter" in gh:
            errs["filter"] = max(errs["filter"], assert_close(
                gh["hypervisc_filter"], got_f, RTOL_KERNEL, what + " hypervisc_filter"))
        if not p.applies_filter:
            assert np.array_equal(np.asarray(got_g), gin["gnew"]), what
        n += 1
    assert n > 0, f"{case}: no after_hyper records in the trace"
    _report(f"{case} hyper_diff ({n} records)", errs)


@pytest.mark.parametrize("case", ALL_CASES)
def test_exb_shear_trace(case):
    """exb_shear on advance_begin vs after_exb_shear, carrying the host state over
    every traced step and attempt (identity when the shear is inactive)."""
    ex, tr = export_of(case), trace_of(case)
    dims = Dims.from_export(ex)
    params = exb.ExbParams.from_export(ex)
    data = to_jax(exb.ExbData.from_params(params))
    steps = tr.steps
    chained = steps == list(range(1, len(steps) + 1))
    state = exb.ExbState.from_export(ex, params)
    errs = {"gnew": 0.0, "fields": 0.0, "shift": 0.0}
    njumps = ncompared = 0
    jumps = []
    prev_istep = 0
    for istep in steps:
        for attempt in tr.attempts(istep):
            rb = tr.record("advance_begin", istep, attempt=attempt)
            found = tr.find("after_exb_shear", istep=istep, attempt=attempt)
            if not found:
                continue
            gb, ga = tr.load(rb), tr.load(found[0])
            ncompared += 1
            if params.active and not chained and attempt == 1:
                # gaps in the trace: re-seed from the record (a retry attempt keeps the
                # state of attempt 1, as GS2 does)
                state = exb.ExbState.from_group(gb, params, istep_last=prev_istep)
            if params.active and state.shift is not None:
                key = "kx_shift" if params.use_kx else "theta0_shift"
                assert_close(gb[key], state.shift, RTOL_KERNEL, f"{case} {istep}: {key} chain")
            if params.active and chained and "exb_first" in gb:
                assert bool(gb["exb_first"]) == state.exb_first, f"{case} {istep}: exb_first"
            time = TimeBlock.from_group(gb)
            fields = to_jax(fields_from_group(gb, dims))
            before = state
            got_g, got_f, state = exb.exb_shear(jnp.asarray(gb["gnew"]), fields, state,
                                                istep, time, data, params)
            what = f"{case} istep={istep} attempt={attempt}"
            errs["gnew"] = max(errs["gnew"], assert_close(ga["gnew"], got_g, RTOL_KERNEL,
                                                          what + " gnew"))
            want_f = fields_from_group(ga, dims)
            for name, a, b in zip(("phinew", "aparnew", "bparnew"), want_f, got_f):
                errs["fields"] = max(errs["fields"], assert_close(a, b, RTOL_KERNEL,
                                                                  f"{what} {name}"))
            if params.active:
                key = "kx_shift" if params.use_kx else "theta0_shift"
                if key in ga:
                    errs["shift"] = max(errs["shift"], assert_close(
                        ga[key], state.shift, RTOL_KERNEL, f"{what} {key}"))
                ran = state.istep_last == istep and before.istep_last != istep
                if "jump" in ga and state.jump is not None and (chained or ran):
                    assert np.array_equal(np.asarray(ga["jump"]), state.jump), what + " jump"
                if "ikx_indexed" in ga and state.ikx_indexed is not None:
                    assert np.array_equal(np.asarray(ga["ikx_indexed"]), state.ikx_indexed)
                if ran and np.any(state.jump != 0):
                    njumps += 1
                    jumps.append((istep, attempt, state.jump.tolist()))
            else:
                assert np.array_equal(np.asarray(got_g), gb["gnew"]), what
        prev_istep = istep
    assert ncompared > 0 or not params.active, f"{case}: no after_exb_shear records"
    print(f"[B7] {case} exb_shear: {ncompared} records, {njumps} steps with a shift, "
          f"chained={chained}; (istep, attempt, jump) = {jumps}")
    _report(f"{case} exb_shear", errs)


@pytest.mark.parametrize("case", ZONAL_CASES)
def test_remove_zonal_flows_trace(case):
    """remove_zonal_flows on after_field_update vs after_remove_zonal_flows."""
    ex, tr = export_of(case), trace_of(case)
    dims = Dims.from_export(ex)
    data, params = to_jax(zonal.ZonalData.from_export(ex)), zonal.ZonalParams.from_export(ex)
    fn = jax.jit(zonal.remove_zonal_flows, static_argnames=("params",))
    errs = {"phinew": 0.0}
    recs = tr.find("after_remove_zonal_flows")
    assert recs, f"{case}: no after_remove_zonal_flows records"
    for rz in recs:
        rf = tr.previous(rz, "after_field_update")
        assert rf is not None
        gz, gf = tr.load(rz), tr.load(rf)
        what = f"{case} istep={rz.istep} attempt={rz.attempt}"
        got = fn(jnp.asarray(gf["phinew"]), data, params)
        errs["phinew"] = max(errs["phinew"], assert_close(gz["phinew"], got, RTOL_KERNEL,
                                                          what + " phinew"))
        before, after = fields_from_group(gf, dims), fields_from_group(gz, dims)
        assert np.array_equal(before.apar, after.apar) and np.array_equal(before.bpar,
                                                                          after.bpar)
    _report(f"{case} remove_zonal_flows ({len(recs)} records)", errs)


@pytest.mark.parametrize("case", ALL_CASES)
def test_field_line_average_weight_export(case):
    """field_line_average_weight (theta_grid.f90:520) recomputed vs the export."""
    ex = export_of(case)
    data = zonal.ZonalData.from_export(ex)
    got = float(zonal.field_line_average_weight(data))
    want = float(ex["/theta_grid/field_line_average_weight"])
    assert abs(got - want) <= RTOL_KERNEL * abs(want)
