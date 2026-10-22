"""Tests of gs2gpu.source (source term g_work, AB coefficients, sourcefac).

* Host functions: ``ab_coefficients`` / ``adams_bashforth_coefficients`` against the
  literal Fortran formula (gs2_time.f90:250-312) for orders 1, 2, 3, changed dt,
  ``maximum_order`` and ``target_dt``; ``sourcefac`` on both sides of ``t0``.
* Synthetic (always run): ``get_source_term_on_local_domain`` against ``ref_source``, a
  literal scalar transcription of get_source_term_on_local_domain / get_source_term /
  set_source (src/dist_fn.fpp:4856-5420), bit for bit in strict mode, for full /
  phiext_full (phi_ext = 0) / homogeneous, linear and nonlinear (AB orders 1-3), ES and
  EM, with ttp points, ``il > lmax`` and a ``kwork_filter`` mode; jit with a traced
  ``istep``.  The trace-driven harness used for the reference data is also run on a
  synthetic trace file written in the layout of spec 09 section 6.2.
* Reference data (skip without ``gpu/refdata``): the trace probe ``source`` (``g_work``)
  for every traced step, attempt and pass, inputs from the earlier probes of the same
  step (spec 01 section 3.3), ``RTOL_KERNEL``; and a strict-mode run (expect bitwise).
"""

from __future__ import annotations

import math

import jax
import jax.numpy as jnp
import netCDF4
import numpy as np
import pytest

from conftest import RTOL_KERNEL, assert_close, export_of, have, all_cases, trace_of
from gs2gpu import coeffs as co
from gs2gpu import source as src
from gs2gpu.io.trace import Trace
from gs2gpu.types import Coeffs, Dims, Fields, TimeBlock, fields_from_group, to_jax
from test_coeffs import (bitwise_equal, ladder_rungs, make_synthetic_export, ref_coeffs,
                         strict_mode, strict_rtol)

EPS = float(np.finfo(np.float64).eps)


# ================================================================ AB coefficients
def fortran_ab(this_dt, p1, p2, order):
    """The Fortran text of get_adams_bashforth_coefficients, transcribed verbatim."""
    if order == 3:
        c1 = (((this_dt +
                p1)**2/p1)*(
            (this_dt+p1)/3. + p2/2.) -
            p1*(p1/3. + p2/2.)) / \
            (this_dt*(p1 + p2))
        c2 = - this_dt * (this_dt/3. +
                          (p1+p2)/2.) / \
            (p1*p2)
        c3 = this_dt * (this_dt/3. +
                        p1/2.) / (
            p2*(p1+p2))
        return (c1, c2, c3)
    if order == 2:
        return ((0.5*this_dt + p1) / p1, -0.5*this_dt / p1)
    return (1.0,)


def _tb(dt, p1, p2, time=0.0):
    return TimeBlock(code_dt=dt, user_dt=dt, code_dt_prev1=p1, code_dt_prev2=p2,
                     code_dt_old=0.0, code_time=time, user_time=time, code_dt_cfl=1e8,
                     user_dt_cfl=1e8, code_dt_min=1e-6, code_dt_max=dt)


def test_time_history_available():
    assert src.time_history_available(-1.0, -1.0) == 1
    assert src.time_history_available(0.1, -1.0) == 2
    assert src.time_history_available(0.1, 0.1) == 3
    assert src.time_history_available(-1.0, 0.1) == 3      # Fortran tests prev2 first


@pytest.mark.parametrize("dt,p1,p2", [
    (0.05, -1.0, -1.0), (0.05, 0.05, -1.0), (0.05, 0.05, 0.05),     # fresh start, equal dt
    (0.025, 0.05, 0.05), (0.025, 0.025, 0.05), (0.1, 0.05, 0.025),  # decrease / increase
    (0.0137, 0.0291, 0.00731), (1.0e-5, 3.0e-3, 2.0e-1),
])
def test_ab_coefficients_literal(dt, p1, p2):
    order = 3 if p2 != -1.0 else (2 if p1 != -1.0 else 1)
    want = fortran_ab(dt, p1, p2, order)
    got = src.ab_coefficients(_tb(dt, p1, p2), istep=7)
    assert len(got) == order
    assert all(isinstance(c, float) for c in got)
    assert got == want                                   # bitwise (Python float = Fortran real)
    assert src.adams_bashforth_coefficients(dt, p1, p2) == want


@pytest.mark.parametrize("dt", [0.1, 0.05, 0.03, 0.0137, 0.01, 1e-3])
def test_ab_coefficients_equal_steps(dt):
    """Equal steps: the variable-step formula, not the textbook rationals (02 pitfall 12);
    for e.g. dt = 0.03, 0.01, 1e-3 they differ in the last bit."""
    got = src.ab_coefficients(_tb(dt, dt, dt), 3)
    assert got == fortran_ab(dt, dt, dt, 3)
    assert np.allclose(got, (23 / 12, -4 / 3, 5 / 12), rtol=1e-15, atol=0)
    if dt in (0.03, 0.01, 1e-3):
        assert got != (23 / 12, -4 / 3, 5 / 12)


def test_ab_coefficients_options():
    dt, p1, p2 = 0.03, 0.05, 0.07
    assert src.adams_bashforth_coefficients(dt, p1, p2, maximum_order=2) == fortran_ab(dt, p1, p2, 2)
    assert src.adams_bashforth_coefficients(dt, p1, p2, maximum_order=1) == (1.0,)
    assert src.adams_bashforth_coefficients(dt, p1, -1.0, maximum_order=3) == fortran_ab(dt, p1, -1.0, 2)
    assert src.adams_bashforth_coefficients(dt, p1, p2, target_dt=0.01) == fortran_ab(0.01, p1, p2, 3)
    assert src.ab_coefficients(_tb(dt, p1, p2), istep=0) == ()        # include_explicit false
    assert src.ab_coefficients(_tb(dt, p1, p2)) == fortran_ab(dt, p1, p2, 3)


def test_sourcefac():
    params = src.SourceParams(1, False, False, False, False, False, False, False,
                              t0=2.0, omega0=0.7, gamma0=0.1, source0=1.5)
    for time in (0.0, 0.3, 2.0, 2.5, 10.0):
        e = complex(math.exp(0.1 * time) * math.cos(-(0.7 * time)),
                    math.exp(0.1 * time) * math.sin(-(0.7 * time)))
        fac = 1.5 if time > 2.0 else 0.5 - 0.5 * math.cos(math.pi * time / 2.0)
        got = src.sourcefac(time, params)
        assert abs(got - fac * e) <= 4 * EPS * abs(fac * e) + 1e-300
        traced = complex(jax.jit(lambda t: src.sourcefac(t, params))(time))
        assert abs(traced - got) <= 8 * EPS * max(abs(got), 1e-300)


# ===================================================== scalar Fortran transcription
def cmplx_list(x):
    return np.asarray(x).tolist()


def fcmla(c: complex, x: complex) -> complex:
    """Fortran ``c*x`` in a loop that gfortran vectorises for apple-m1 (FCMLA #0, #90):
    ``re = fma(-c.im, x.im, c.re*x.re)``, ``im = fma(c.im, x.re, c.re*x.im)``."""
    return complex(math.fma(-c.imag, x.imag, c.real * x.real),
                   math.fma(c.imag, x.real, c.real * x.imag))


def ref_source(ex, coeffs: dict, g, fields_old: Fields, fields_new: Fields, gexp, ab,
               istep: int, sourcefac=0j):
    """Literal scalar transcription of get_source_term_on_local_domain
    (src/dist_fn.fpp:5363-5420), get_source_term (4856-5012) and set_source
    (5017-5128), ``opt_source = .false.``, non-LOWFLOW branches.

    All operands are Python floats/complex (see ``test_coeffs.ref_coeffs``), so the
    arithmetic is IEEE without FMA, Fortran association order, and GCC's componentwise
    handling of real*complex -- except the matrix terms ``b*g + a*g``, whose loop gfortran
    vectorises with FCMLA fused complex products (:func:`fcmla`; verified against a
    gfortran build of the verbatim loop with the reference flags).  The ttp loop is
    scalar code: naive products.
    """
    da, sp, df, le, kt, rp = (ex[p] for p in ("/dist_fn_arrays", "/species", "/dist_fn",
                                              "/le_grids", "/kt_grids", "/run_parameters"))
    S, E, L, Y, X, _, T = np.shape(g)
    lst = lambda grp, n: np.asarray(grp[n]).tolist()          # noqa: E731
    aj0, aj1, vperp2, vpac = (lst(da, n) for n in ("aj0", "aj1", "vperp2", "vpac"))
    vpar, wdrift, wdriftttp = (cmplx_list(coeffs[n]) for n in ("vpar", "wdrift", "wdriftttp"))
    wstar = cmplx_list(coeffs["wstar"])
    a, b = cmplx_list(coeffs["a"]), cmplx_list(coeffs["b"])
    code_dt = float(coeffs["code_dt"])
    gl = cmplx_list(g)
    phi, apar, bpar = (cmplx_list(f) for f in fields_old)
    phinew, aparnew, bparnew = (cmplx_list(f) for f in fields_new)
    if gexp is not None:
        gexp_1, gexp_2, gexp_3 = (cmplx_list(x) for x in gexp)
    nonmaxw_corr, ufac = lst(sp, "nonmaxw_corr"), lst(df, "ufac")
    zstm, stm, tz = lst(sp, "spec_zstm"), lst(sp, "spec_stm"), lst(sp, "spec_tz")
    fexp, bkdiff = np.atleast_1d(df["fexp"]).tolist(), np.atleast_1d(df["bkdiff"]).tolist()
    D_res = np.asarray(ex["/hyper/D_res"]).tolist()
    itor_over_B = np.asarray(ex["/theta_grid/itor_over_B"]).tolist()
    wunits = np.asarray(ex["/gs2_time/wunits"]).tolist()
    aky = lst(kt, "aky")
    kwork_filter = lst(kt, "kwork_filter")
    is_ttp, can_be_ttp = lst(le, "is_ttp"), lst(le, "can_be_ttp")
    lmax = int(le["lmax"])
    trapped = L > int(le["ng2"])
    fphi, fapar, fbpar = float(rp["fphi"]), float(rp["fapar"]), float(rp["fbpar"])
    omprimfac, g_exb = float(df["omprimfac"]), float(df["g_exb"])
    phi_ext = float(df["phi_ext"])
    option = int(df["source_option_switch"])
    nonlin = bool(ex["/nonlinear_terms/nonlin"])
    split = bool(ex["/nonlinear_terms/split_nonlinear"])
    zi = 1j
    epsilon = EPS
    out = np.zeros(np.shape(g), complex)

    include_explicit = istep > 0
    include_explicit = include_explicit and nonlin and (not split)
    ab_coefficients = list(ab) if include_explicit else None

    for s in range(S):
        for e in range(E):
            for l in range(L):
                il = l + 1
                for y in range(Y):
                    for x in range(X):
                        if kwork_filter[y][x]:
                            continue
                        for isgn in (1, 2):
                            sg = isgn - 1
                            source = [0j] * T
                            phigavg = [0j] * T
                            apargavg = [0j] * T
                            for ig in range(T):
                                phigavg[ig] = (fexp[s]*phi[y][x][ig] + (1.0-fexp[s])*phinew[y][x][ig]) \
                                    * aj0[s][e][l][y][x][ig]*fphi \
                                    + (fexp[s]*bpar[y][x][ig] + (1.0-fexp[s])*bparnew[y][x][ig]) \
                                    * aj1[s][e][l][y][x][ig]*fbpar*2.0*vperp2[s][e][l][y][x][ig]*tz[s]
                                apargavg[ig] = (fexp[s]*apar[y][x][ig] + (1.0-fexp[s])*aparnew[y][x][ig]) \
                                    * aj0[s][e][l][y][x][ig]*fapar

                            def explicit_terms(ig):
                                k = len(ab_coefficients)
                                if k == 1:
                                    return 0.5*code_dt*(
                                        ab_coefficients[0]*gexp_1[s][e][l][y][x][sg][ig])
                                if k == 2:
                                    return 0.5*code_dt*(
                                        ab_coefficients[0]*gexp_1[s][e][l][y][x][sg][ig] +
                                        ab_coefficients[1]*gexp_2[s][e][l][y][x][sg][ig])
                                return 0.5*code_dt*(
                                    ab_coefficients[0]*gexp_1[s][e][l][y][x][sg][ig] +
                                    ab_coefficients[1]*gexp_2[s][e][l][y][x][sg][ig] +
                                    ab_coefficients[2]*gexp_3[s][e][l][y][x][sg][ig])

                            def set_source():
                                bd = bkdiff[s]
                                bdfac_p = 1.+bd*(3.-2.*float(isgn))
                                bdfac_m = 1.-bd*(3.-2.*float(isgn))
                                nm = nonmaxw_corr[s][e]
                                for ig in range(T - 1):
                                    phi_p = bdfac_p*phigavg[ig+1]+bdfac_m*phigavg[ig]
                                    phi_m = phigavg[ig+1]-phigavg[ig]
                                    apar_p = bdfac_p*apargavg[ig+1]+bdfac_m*apargavg[ig]
                                    apar_m = bdfac_p * (aparnew[y][x][ig+1] - apar[y][x][ig+1]) * fapar + \
                                        bdfac_m * (aparnew[y][x][ig] - apar[y][x][ig]) * fapar
                                    vpr = vpar[s][e][l][y][x][sg][ig]
                                    vpc = vpac[s][e][l][y][x][sg][ig]
                                    wdr = wdrift[s][e][l][y][x][sg][ig]
                                    source[ig] = (-(2.0*vpr*nm*phi_m)
                                                  - zstm[s]*vpc*nm
                                                  * ((aj0[s][e][l][y][x][ig+1] + aj0[s][e][l][y][x][ig])*0.5*apar_m
                                                     + D_res[y][x]*apar_p)
                                                  - zi*wdr*nm*phi_p) \
                                        + zi*(wstar[s][e][y]
                                              + vpc*code_dt*wunits[y]*ufac[s][e]
                                              - 2.0*omprimfac*vpc*nm*code_dt*wunits[y]*g_exb*itor_over_B[ig]/stm[s]) \
                                        * (phi_p - apar_p*stm[s]*vpc)
                                if include_explicit:
                                    for ig in range(T - 1):
                                        source[ig] = source[ig] + explicit_terms(ig)

                            if option == 1:
                                if il <= lmax:
                                    set_source()
                                else:
                                    source = [0j] * T
                            elif option == 5:
                                if il <= lmax:
                                    set_source()
                                    if aky[y] < epsilon:
                                        for ig in range(T - 1):
                                            source[ig] = source[ig] \
                                                - zi*wdrift[s][e][l][y][x][sg][ig]*nonmaxw_corr[s][e] \
                                                * 2.0*phi_ext*sourcefac*aj0[s][e][l][y][x][ig]
                                else:
                                    source = [0j] * T
                            elif option == 6:
                                source = [0j] * T
                            if il <= lmax:
                                # vectorised by gfortran: complex products are FCMLA pairs
                                if isgn == 1:
                                    for ig in range(T - 1):
                                        source[ig] = source[ig] \
                                            + fcmla(b[s][e][l][y][x][0][ig], gl[s][e][l][y][x][0][ig]) \
                                            + fcmla(a[s][e][l][y][x][0][ig], gl[s][e][l][y][x][0][ig+1])
                                else:
                                    for ig in range(T - 1):
                                        source[ig] = source[ig] \
                                            + fcmla(a[s][e][l][y][x][1][ig], gl[s][e][l][y][x][1][ig]) \
                                            + fcmla(b[s][e][l][y][x][1][ig], gl[s][e][l][y][x][1][ig+1])
                            if option in (1, 5):
                                if trapped and isgn == 2 and can_be_ttp[l]:
                                    for ig in range(T - 1):
                                        if not is_ttp[l][ig]:
                                            continue
                                        source[ig] = gl[s][e][l][y][x][1][ig]*a[s][e][l][y][x][1][ig] \
                                            - zi*(wdriftttp[s][e][l][y][x][1][ig])*nonmaxw_corr[s][e]*phigavg[ig] \
                                            + zi*wstar[s][e][y]*phigavg[ig]
                                    if option == 5 and aky[y] < epsilon:
                                        for ig in range(T - 1):
                                            if not is_ttp[l][ig]:
                                                continue
                                            source[ig] = source[ig] - zi * \
                                                wdriftttp[s][e][l][y][x][sg][ig]*nonmaxw_corr[s][e] \
                                                * 2.0*phi_ext*sourcefac*aj0[s][e][l][y][x][ig]
                                    if include_explicit:
                                        for ig in range(T - 1):
                                            if not is_ttp[l][ig]:
                                                continue
                                            source[ig] = source[ig] + explicit_terms(ig)
                            out[s, e, l, y, x, sg, :] = source
    return out


# ============================================================== synthetic state / runner
def random_state(ex, seed: int, *, nonlin: bool):
    rng = np.random.default_rng(seed)
    gshape = ex.g_shape
    fshape = ex.field_shape

    def c(shape, scale=1.0):
        return scale * (rng.normal(size=shape) + 1j * rng.normal(size=shape))

    has_apar = bool(ex["/run_parameters/has_apar"])
    has_bpar = bool(ex["/run_parameters/has_bpar"])
    zero = np.zeros(fshape, complex)
    old = Fields(c(fshape), c(fshape, 0.3) if has_apar else zero, c(fshape, 0.2) if has_bpar else zero)
    new = Fields(old.phi + c(fshape, 0.05), old.apar + c(fshape, 0.05) if has_apar else zero,
                 old.bpar + c(fshape, 0.05) if has_bpar else zero)
    g = c(gshape)
    gexp = tuple(c(gshape, 0.1) for _ in range(3)) if nonlin else None
    return g, old, new, gexp


def coeff_dict(c: Coeffs) -> dict:
    return {n: np.asarray(getattr(c, n)) for n in ("vpar", "wdrift", "wdriftttp", "wstar", "a", "b")} \
        | {"code_dt": float(c.code_dt)}


def jit_source(params, static_istep=True):
    if static_istep:
        return jax.jit(lambda g, fo, fn, gexp, ab, coeffs, data, istep, code_time:
                       src.get_source_term_on_local_domain(g, fo, fn, gexp, ab, istep, coeffs,
                                                           data, params, code_time=code_time),
                       static_argnums=(7,))
    return jax.jit(lambda g, fo, fn, gexp, ab, coeffs, data, istep, code_time:
                   src.get_source_term_on_local_domain(g, fo, fn, gexp, ab, istep, coeffs, data,
                                                       params, code_time=code_time))


SYNTH = {
    "linear_es_full": dict(ex=dict(nonlin=False, has_apar=False, has_bpar=False), istep=1, ab=None),
    "linear_em_full": dict(ex=dict(nonlin=False), istep=3, ab=None),
    "nl_em_ab1": dict(ex=dict(nonlin=True), istep=1, ab=(-1.0, -1.0)),
    "nl_em_ab2": dict(ex=dict(nonlin=True, seed=1), istep=2, ab=(0.05, -1.0)),
    "nl_es_ab3_dtchange": dict(ex=dict(nonlin=True, has_apar=False, has_bpar=False, seed=2),
                               istep=5, ab=(0.1, 0.05)),
    "nl_istep0_probe": dict(ex=dict(nonlin=True), istep=0, ab=(0.05, 0.05)),
    "phiext_full_zero": dict(ex=dict(source_option=5), istep=2, ab=None),
    "homogeneous": dict(ex=dict(source_option=6, nonlin=True), istep=2, ab=(0.05, 0.05)),
    "untrapped_maxwell": dict(ex=dict(trapped=False, nonmaxwell=False, box=False, zonal=False,
                                      g_exb=0.0, mach=0.0, seed=4), istep=1, ab=None),
    "one_species_bpar_only": dict(ex=dict(nspec=1, negrid=2, ntgrid=3, naky=2, ntheta0=3,
                                          has_apar=False, nonlin=True, seed=5),
                                  istep=4, ab=(0.05, 0.05)),
}


def synth_case(name):
    cfg = SYNTH[name]
    exk = dict(cfg["ex"])
    seed = exk.pop("seed", 0)
    ex = make_synthetic_export(seed, ladder=(), **exk)
    nonlin = exk.get("nonlin", False)
    g, old, new, gexp = random_state(ex, 100 + seed, nonlin=nonlin)
    coeffs = Coeffs.from_export(ex)
    if cfg["ab"] is None:
        ab = ()
    else:
        p1, p2 = cfg["ab"]
        ab = src.ab_coefficients(_tb(coeffs.code_dt, p1, p2), cfg["istep"])
    return ex, g, old, new, gexp, ab, cfg["istep"], coeffs


@pytest.fixture(scope="module", params=sorted(SYNTH))
def synth(request):
    name = request.param
    ex, g, old, new, gexp, ab, istep, coeffs = synth_case(name)
    sfac = src.sourcefac(0.0, src.SourceParams.from_export(ex))   # _run uses code_time = 0
    ref = ref_source(ex, coeff_dict(coeffs), g, old, new, gexp, ab, istep, sfac)
    # g = 0 (as in the response probes): the matrix terms (~2g) no longer swamp
    # last-bit differences in the set_source terms, so the strict test is sharper.
    ref0 = ref_source(ex, coeff_dict(coeffs), np.zeros_like(g), old, new, gexp, ab, istep, sfac)
    return name, ex, g, old, new, gexp, ab, istep, coeffs, ref, ref0


def _run(ex, g, old, new, gexp, ab, istep, coeffs, *, jit=True, static_istep=True, code_time=0.0):
    data, params = to_jax(src.SourceData.from_export(ex)), src.SourceParams.from_export(ex)
    args = (jnp.asarray(g), Fields(*map(jnp.asarray, old)), Fields(*map(jnp.asarray, new)),
            None if gexp is None else tuple(map(jnp.asarray, gexp)), tuple(ab), to_jax(coeffs),
            data)
    if not jit:
        return src.get_source_term_on_local_domain(args[0], args[1], args[2], args[3], args[4],
                                                   istep, args[5], data, params,
                                                   code_time=code_time)
    return jit_source(params, static_istep)(*args, istep, code_time)


def test_source_synthetic_strict_bitwise(synth):
    name, ex, g, old, new, gexp, ab, istep, coeffs, ref, ref0 = synth
    with strict_mode():
        got = _run(ex, g, old, new, gexp, ab, istep, coeffs)
        got0 = _run(ex, np.zeros_like(g), old, new, gexp, ab, istep, coeffs)
    for r, x, what in ((ref, got, "g"), (ref0, got0, "g=0")):
        assert bitwise_equal(r, x), \
            f"{name} ({what}): {np.sum(np.asarray(x) != r)} elements differ, max-abs-rel " \
            f"{np.max(np.abs(np.asarray(x) - r)) / np.max(np.abs(r)):.3e}"


def test_source_synthetic_default_mode(synth):
    name, ex, g, old, new, gexp, ab, istep, coeffs, ref, _ = synth
    got = _run(ex, g, old, new, gexp, ab, istep, coeffs)
    assert_close(ref, got, 1e-14, f"{name} jit")
    got = _run(ex, g, old, new, gexp, ab, istep, coeffs, jit=False)
    assert_close(ref, got, 1e-14, f"{name} eager")


def test_source_synthetic_traced_istep(synth):
    name, ex, g, old, new, gexp, ab, istep, coeffs, ref, _ = synth
    got = _run(ex, g, old, new, gexp, ab, istep, coeffs, static_istep=False)
    assert_close(ref, got, 1e-14, f"{name} traced istep")


def test_fcmla_mul_matches_fma_emulation():
    """XLA:CPU in strict mode reproduces the FCMLA pair exactly (barriered c.re products,
    fused c.im products); naive products differ in ~25 % of the cases."""
    rng = np.random.default_rng(11)
    c = rng.normal(size=(64, 50)) + 1j * rng.normal(size=(64, 50))
    x = rng.normal(size=(64, 50)) + 1j * rng.normal(size=(64, 50))
    want = np.vectorize(fcmla)(c, x)
    with strict_mode():
        got = jax.jit(lambda c, x: src.fcmla_mul(c, x))(c, x)
    assert bitwise_equal(want, got)
    assert np.sum(np.vectorize(lambda u, v: complex(u.real * v.real - u.imag * v.imag,
                                                    u.real * v.imag + u.imag * v.real))(c, x)
                  != want) > 0


def test_source_structure():
    """ig = ntgrid and kwork_filter entries are 0; il > lmax keeps only ttp points;
    pass 1 (phinew == phi) has apar_m == 0 exactly."""
    ex, g, old, new, gexp, ab, istep, coeffs = synth_case("nl_em_ab2")
    got = np.asarray(_run(ex, g, old, old, gexp, ab, istep, coeffs))
    assert np.all(got[..., -1] == 0)
    kw = ex["/kt_grids/kwork_filter"]
    assert kw.any() and np.all(got[:, :, :, kw] == 0)
    le = ex["/le_grids"]
    lmax = int(le["lmax"])
    beyond = got[:, :, lmax:]                                   # il > lmax (0-based l >= lmax)
    ttp = le["can_be_ttp"][lmax:, None] & le["is_ttp"][lmax:]
    mask = np.zeros(beyond.shape, bool)
    mask[..., 1, :] = ttp[None, None, :, None, None, :]
    mask[:, :, :, kw] = False
    assert np.all(beyond[~mask] == 0) and np.all(beyond[mask] != 0)


def test_source_option_gates():
    ex = make_synthetic_export(0, ladder=(), source_option=5)
    ex["/dist_fn"].arrays["phi_ext"] = 0.1
    with pytest.raises(NotImplementedError):
        src.SourceParams.from_export(ex)
    ex = make_synthetic_export(0, ladder=())
    ex["/dist_fn"].arrays["opt_source"] = True
    with pytest.raises(NotImplementedError):
        src.SourceParams.from_export(ex)


def test_source_with_recomputed_coeffs():
    """The kernel fed with compute_coeffs output (what the integration does)."""
    ex, g, old, new, gexp, ab, istep, _ = synth_case("nl_em_ab2")
    inputs, cparams = to_jax(co.CoeffInputs.from_export(ex)), co.CoeffParams.from_export(ex)
    coeffs = co.compute_coeffs(0.02, inputs, cparams)
    ref = ref_source(ex, coeff_dict(Coeffs(code_dt=0.02, g_h=None, **{
        k: v for k, v in ref_coeffs(ex, 0.02).items()})), g, old, new, gexp, ab, istep)
    got = _run(ex, g, old, new, gexp, ab, istep, coeffs)
    assert_close(ref, got, 1e-14, "recomputed coeffs")


# ============================================================ trace-driven harness
def _coeffs_for(ex, code_dt, cache, inputs, cparams):
    """Exported coefficients for ``code_dt`` (initial dt or a ladder rung, matched by bit
    pattern), else recomputed with ``compute_coeffs`` (cached per dt)."""
    if code_dt in cache:
        return cache[code_dt]
    base = Coeffs.from_export(ex)
    found = base if base.code_dt == code_dt else None
    for k in ladder_rungs(ex):
        rung = Coeffs.from_export(ex, rung=k)
        if found is None and rung.code_dt == code_dt:
            found = rung
    if found is None:
        found = co.compute_coeffs(code_dt, inputs, cparams)
    cache[code_dt] = to_jax(found)
    return cache[code_dt]


def source_probe_inputs(tr: Trace, rec, dims: Dims):
    """Inputs of the ``source`` probe ``rec`` from the earlier probes of the same step
    and attempt (spec 01 section 3.3, 09 section 5.3):

    * ``g`` and ``fields_old``: ``gnew``/``phinew...`` at ``after_exb_shear`` (else
      ``advance_begin``) -- the copies ``g = gnew``, ``phi = phinew`` of
      fields_implicit.f90:337-338;
    * ``fields_new``: pass 1 = ``fields_old``; pass 2 = the latest of
      ``after_field_update`` / ``after_remove_zonal_flows`` before ``rec``;
    * ``gexp``: ``gexp_1..3`` at ``after_add_explicit_terms`` (if allocated);
    * the time block of ``rec`` (``code_dt``, ``code_dt_prev1/2``, ``code_time``).
    """
    grp = tr.load(rec)
    start = tr.previous(rec, "after_exb_shear") or tr.previous(rec, "advance_begin")
    if start is None:
        raise AssertionError(f"no advance_begin before {rec}")
    sgrp = tr.load(start)
    fields_old = fields_from_group(sgrp, dims)
    cands = [r for r in (tr.previous(rec, "after_field_update"),
                         tr.previous(rec, "after_remove_zonal_flows")) if r is not None]
    fields_new = fields_old
    if cands:
        fields_new = fields_from_group(tr.load(max(cands, key=lambda r: r.seq)), dims)
    gexp = None
    erec = tr.previous(rec, "after_add_explicit_terms")
    if erec is not None:
        eg = tr.load(erec)
        if "gexp_1" in eg:
            gexp = (eg["gexp_1"], eg["gexp_2"], eg["gexp_3"])
    return grp["g_work"], sgrp["gnew"], fields_old, fields_new, gexp, TimeBlock.from_group(grp)


def check_source_trace(ex, tr: Trace, *, rtol=RTOL_KERNEL, strict=False, label=""):
    """Run the source kernel for every ``source`` probe of the trace; returns the errors."""
    dims = Dims.from_export(ex)
    data, params = to_jax(src.SourceData.from_export(ex)), src.SourceParams.from_export(ex)
    inputs, cparams = to_jax(co.CoeffInputs.from_export(ex)), co.CoeffParams.from_export(ex)
    cache: dict = {}
    recs = tr.find("source")
    if not recs and not tr.find("step_end"):
        pytest.skip(f"{label}: GS2 aborted before computing a source term (no completed step; "
                    "the exit path is covered by test_end_to_end)")
    assert recs, f"{label}: trace has no 'source' probes"
    errs = []
    with strict_mode(strict):
        fn = jit_source(params)
        for rec in recs:
            ref, g, fo, fnew, gexp, time = source_probe_inputs(tr, rec, dims)
            coeffs = _coeffs_for(ex, time.code_dt, cache, inputs, cparams)
            ab = src.ab_coefficients(time, rec.istep)
            if params.nonlin and rec.istep > 0:
                assert gexp is not None, f"{label}: no gexp before {rec.path}"
            got = fn(jnp.asarray(g), Fields(*map(jnp.asarray, fo)), Fields(*map(jnp.asarray, fnew)),
                     None if gexp is None else tuple(map(jnp.asarray, gexp)), tuple(ab), coeffs,
                     data, rec.istep, time.code_time)
            what = f"{label} istep={rec.istep} attempt={rec.attempt} pass={rec.pass_}"
            errs.append((rec.istep, rec.attempt, rec.pass_, assert_close(ref, got, rtol, what)))
    return errs


# --------------------------------------------------------------- synthetic trace file
def _write_var(grp, name, value):
    value = np.asarray(value)
    if np.iscomplexobj(value):
        raw = np.stack([value.real, value.imag], axis=-1)
        ftype = "complex"
    elif value.dtype == bool:
        raw = value.astype(np.int8)
        ftype = "logical"
    else:
        raw = value.astype(np.float64)
        ftype = "real"
    dims = []
    for k, n in enumerate(raw.shape):
        dname = f"d{n}"
        if dname not in grp.dimensions and dname not in _root_dims(grp):
            _root(grp).createDimension(dname, n)
        dims.append(dname)
    var = grp.createVariable(name, raw.dtype, tuple(dims))
    var.setncattr("fortran_type", ftype)
    var[...] = raw


def _root(grp):
    while grp.parent is not None:
        grp = grp.parent
    return grp


def _root_dims(grp):
    return _root(grp).dimensions


def write_synthetic_trace(path, records):
    """records: list of (istep, attempt, seq, tag, pass, arrays dict)."""
    with netCDF4.Dataset(path, "w") as ds:
        for istep, attempt, seq, tag, pass_, arrays in records:
            sname, aname = f"s{istep:06d}", f"a{attempt}"
            sg = ds.groups.get(sname) or ds.createGroup(sname)
            ag = sg.groups.get(aname) or sg.createGroup(aname)
            pg = ag.createGroup(f"p{seq:03d}_{tag}")
            for k, v in dict(tag=tag, istep=istep, attempt=attempt, seq=seq, pass_=pass_,
                             in_advance=1).items():
                pg.setncattr("pass" if k == "pass_" else k, v)
            for name, value in arrays.items():
                _write_var(pg, name, value)


def make_synthetic_trace(path):
    """A 3-step nonlinear EM synthetic run (AB orders 1, 2, 3 and a dt halving at step 3
    = ladder rung k1, plus an immediate-reset attempt at step 2), with ``g_work`` from
    ``ref_source``.  Probe contents follow spec 09 section 5.3."""
    ex = make_synthetic_export(7, nonlin=True, ladder=(0.025,))
    dts = [0.05, 0.05, 0.025]
    p1, p2, time = -1.0, -1.0, 0.0
    records = []
    for istep, dt in enumerate(dts, start=1):
        attempts = [1, 2] if istep == 2 else [1]
        for attempt in attempts:
            seq = 0
            g, old, new, gexp = random_state(ex, 1000 + 10 * istep + attempt, nonlin=True)
            tb = _tb(dt, p1, p2, time)._asdict()

            def rec(tag, pass_, arrays):
                nonlocal seq
                seq += 1
                records.append((istep, attempt, seq, tag, pass_, tb | dict(arrays)))

            fdict = lambda f, suffix="new": {f"phi{suffix}": f.phi, f"apar{suffix}": f.apar,  # noqa: E731
                                             f"bpar{suffix}": f.bpar}
            rec("advance_begin", 0, {"gnew": g + 0.5, **fdict(old)})      # pre-shear state
            rec("after_exb_shear", 0, {"gnew": g, **fdict(old)})
            rec("after_add_explicit_terms", 1, {"gexp_1": gexp[0], "gexp_2": gexp[1],
                                                "gexp_3": gexp[2]})
            if istep == 2 and attempt == 1:
                continue                                            # reset: return before invert
            coeffs = Coeffs.from_export(ex) if dt == 0.05 else Coeffs.from_export(ex, rung=1)
            ab = src.ab_coefficients(TimeBlock(**tb), istep)
            src1 = ref_source(ex, coeff_dict(coeffs), g, old, old, gexp, ab, istep)
            rec("source", 1, {"g_work": src1})
            rec("after_invert_rhs_1", 1, {"gnew": src1 * 0.5})
            rec("after_getfield", 0, fdict(Fields(new.phi - old.phi, new.apar - old.apar,
                                                  new.bpar - old.bpar)))
            rec("after_field_update", 0, fdict(new))
            src2 = ref_source(ex, coeff_dict(coeffs), g, old, new, gexp, ab, istep)
            rec("source", 2, {"g_work": src2})
        p2, p1, time = p1, dt, time + dt
    write_synthetic_trace(path, records)
    return ex


def test_source_trace_harness_synthetic(tmp_path):
    path = tmp_path / "synthetic.gs2trace.nc"
    ex = make_synthetic_trace(path)
    tr = Trace(path)
    try:
        assert tr.steps == [1, 2, 3] and tr.attempts(2) == [1, 2]
        errs = check_source_trace(ex, tr, rtol=1e-14, label="synthetic")
        assert len(errs) == 6
        errs_strict = check_source_trace(ex, tr, rtol=0.0, strict=True, label="synthetic strict")
        assert all(e == 0.0 for *_, e in errs_strict)
    finally:
        tr.close()


# ================================================================== refdata tests
TRACED = [c for c in all_cases() if have(c, "export") and have(c, "trace")] or ["<no-refdata>"]


@pytest.mark.parametrize("case", TRACED)
def test_source_vs_trace(case):
    ex, tr = export_of(case), trace_of(case)
    errs = check_source_trace(ex, tr, label=case)
    print(f"{case}: max error {max(e for *_, e in errs):.2e} over {len(errs)} source probes: "
          + ", ".join(f"s{i}a{a}p{p}={e:.1e}" for i, a, p, e in errs))


@pytest.mark.parametrize("case", TRACED)
def test_source_vs_trace_strict(case):
    """Strict mode (no FMA, FCMLA-mirrored matrix terms) against the -ffp-contract=off
    reference build: bitwise."""
    ex, tr = export_of(case), trace_of(case)
    errs = check_source_trace(ex, tr, rtol=strict_rtol(ex), strict=True, label=f"{case} strict")
    print(f"{case} strict: max error {max(e for *_, e in errs):.2e}")
