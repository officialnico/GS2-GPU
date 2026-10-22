"""Tests of the field equations and the field solve (B4: fields, fields_local, response).

Two groups:

* synthetic tests (always run): a small fake export with linked supercells of 1, 2 and 3
  cells, compared with literal NumPy/Python transcriptions of the Fortran loops (bitwise in
  numerics strict mode), plus algebraic checks of the response construction;
* reference-data tests (skip without gpu/refdata): every kernel against the GS2 export and the
  stage probes of the trace (spec 09 section 5).
"""

from __future__ import annotations

import math
from functools import lru_cache, partial
from pathlib import Path

import jax
import jax.numpy as jnp
import numpy as np
import pytest

from conftest import RTOL_KERNEL, assert_close, cases_where, export_of, trace_of
from gs2gpu import fields as F
from gs2gpu import fields_local as FL
from gs2gpu import numerics
from gs2gpu import response as RS
from gs2gpu.io.export import Export, ExportError, Group
from gs2gpu.types import Dims, Fields, fields_from_group, to_jax

EPS = F.EPS


# =========================================================================== helpers
class strict_mode:
    """numerics strict mode (no FMA contraction) for functions traced inside the block."""

    def __enter__(self):
        numerics.set_strict(True)
        jax.clear_caches()

    def __exit__(self, *exc):
        numerics.set_strict(False)
        jax.clear_caches()


def maxrel(ref, got):
    return numerics.max_abs_rel(np.asarray(ref), np.asarray(got))


def cm(a, b):
    """Complex product with separately rounded real operations (no contraction)."""
    ar, ai, br, bi = np.real(a), np.imag(a), np.real(b), np.imag(b)
    return np.subtract(np.multiply(ar, br), np.multiply(ai, bi)) + \
        1j * np.add(np.multiply(ar, bi), np.multiply(ai, br))


def rm(z, r):
    return np.multiply(np.real(z), r) + 1j * np.multiply(np.imag(z), r)


def sdiv(a, b):
    """GCC Smith division (spec 04 P2), separately rounded operations, scalar divisor."""
    ar, ai = np.real(a), np.imag(a)
    br, bi = float(np.real(b)), float(np.imag(b))
    if abs(br) < abs(bi):
        r = br / bi
        den = br * r + bi
        return np.divide(np.add(np.multiply(ar, r), ai), den) + \
            1j * np.divide(np.subtract(np.multiply(ai, r), ar), den)
    r = bi / br
    den = bi * r + br
    return np.divide(np.add(np.multiply(ai, r), ar), den) + \
        1j * np.divide(np.subtract(ai, np.multiply(ar, r)), den)


def _chains_to_links(chains, Y, X):
    itright = np.full((Y, X), -1, np.int32)
    itleft = np.full((Y, X), -1, np.int32)
    l_links = np.zeros((Y, X), np.int32)
    r_links = np.zeros((Y, X), np.int32)
    for ik, chs in chains.items():
        for ch in chs:
            if ik == 0:      # ky = 0 in a box: self-linked, N = 1 (dist_fn.fpp:2441-2515)
                itright[ik, ch[0]] = ch[0] + 1
                itleft[ik, ch[0]] = ch[0] + 1
                continue
            for pos, it in enumerate(ch):
                l_links[ik, it] = pos
                r_links[ik, it] = len(ch) - 1 - pos
            for a, b in zip(ch, ch[1:]):
                itright[ik, a] = b + 1
                itleft[ik, b] = a + 1
    return itleft, itright, l_links, r_links


def make_synth(*, S=2, E=3, L=3, ntgrid=2, linked=True, fields=(True, True, True),
               species_types=(1, 2), adiabatic=F.ADIABATIC_DEFAULT, esv=False, layout="lxyes",
               zero_wl=True, seed=0, tite=1.0, beta=0.04, apfac=1.0, field_option=1,
               with_denoms=True, with_aminv=True, filter11=True):
    """A small consistent fake export (names of spec 09; see fields.FieldData.from_export)."""
    rng = np.random.default_rng(seed)
    T = 2 * ntgrid + 1
    if linked:
        aky = np.array([0.0, 0.25, 0.5])
        akx = np.array([0.0, 0.1, 0.2, -0.2, -0.1])
        Y, X = 3, 5
        chains = {0: [[it] for it in range(X)], 1: [[0, 1], [2], [3, 4]],
                  2: [[4, 0, 2], [1], [3]]}
    else:
        aky = np.array([0.0, 0.3])
        akx = np.array([0.0, 0.15])
        Y, X = 2, 2
        chains = {0: [[0], [1]], 1: [[0], [1]]}
    itleft, itright, l_links, r_links = _chains_to_links(chains, Y, X)
    lens = sorted({len(ch) for chs in chains.values() for ch in chs}) if linked else [1]
    N_class = np.asarray(lens, np.int32)
    M_class = np.asarray([sum(len(ch) == n for chs in chains.values() for ch in chs)
                          if linked else Y * X for n in lens], np.int32)
    has_phi, has_apar, has_bpar = fields
    nf = int(has_phi) + int(has_apar) + int(has_bpar)
    S = len(species_types)
    g6 = (S, E, L, Y, X, T)
    aj0 = rng.uniform(0.3, 1.0, g6)
    aj1 = rng.uniform(0.1, 0.5, g6)
    vperp2 = rng.uniform(0.1, 2.0, g6)
    vp = rng.uniform(0.1, 2.0, g6)
    vpa = np.stack([vp, -vp], axis=-2)
    w = rng.uniform(0.05, 0.3, (S, E))
    wl = rng.uniform(0.05, 0.3, (L, T))
    forbid = np.zeros((L, T), bool)
    if zero_wl:
        forbid[-1, 0] = forbid[-1, -1] = True
        wl[forbid] = 0.0
    theta = np.linspace(-np.pi, np.pi, T)
    bmag = rng.uniform(0.8, 1.2, T)
    jacob = rng.uniform(0.5, 1.5, T)
    kperp2 = rng.uniform(0.05, 1.0, (Y, X, T))
    kperp2[0, 0] = 0.0
    inv_kperp2 = np.where(np.abs(kperp2) > EPS, 1.0 / np.where(kperp2 == 0, 1, kperp2), 0.0)
    gamtot = rng.uniform(1.0, 2.0, (Y, X, T))
    gamtot1 = rng.uniform(0.0, 0.3, (Y, X, T))
    gamtot2 = rng.uniform(0.0, 0.3, (Y, X, T))
    gridfac1 = np.ones((Y, X, T))
    gridfac1[1:, :, 0] = 1.3
    kwork = np.zeros((Y, X), bool)
    if filter11 and linked:
        kwork[0, 0] = True
    dist_fn = dict(gamtot=gamtot, gamtot1=gamtot1, gamtot2=gamtot2, gridfac1=gridfac1,
                   N_class=N_class, M_class=M_class, i_class=len(lens), l_links=l_links,
                   r_links=r_links, itright=itright, itleft=itleft,
                   boundary_option_switch=F.BOUNDARY_LINKED if linked else F.BOUNDARY_ZERO,
                   adiabatic_option_switch=adiabatic, apfac=apfac, poisfac=0.0, esv=esv,
                   afilter=0.0, hyper_in_initialisation=True)
    absent = []
    if adiabatic == F.ADIABATIC_FIELDLINEAVG:
        dist_fn["gamtot3"] = F.gamtot3_from_gamtot(gamtot, tite)
        aw = F.awgt_host(theta, jacob, dist_fn["gamtot3"], aky)
        dist_fn["awgt"] = aw
    if with_denoms:
        a, b = F.denominators(gamtot, gamtot1, gamtot2, bmag, beta, has_phi, has_bpar)
        dist_fn["inv_phi_denominator_g"] = a
        dist_fn["inv_bpar_denominator_g"] = b
    else:
        absent += ["inv_phi_denominator_g", "inv_bpar_denominator_g"]
    sz = np.array([1.0, -1.0, 1.0][:S])
    sd = rng.uniform(0.5, 1.5, S)
    st = rng.uniform(0.5, 2.0, S)
    sm = np.array([1.0, 0.01, 2.0][:S])
    groups = {
        "theta_grid": dict(theta=theta, bmag=bmag, jacob=jacob),
        "kt_grids": dict(aky=aky, akx=akx, kperp2=kperp2, inv_kperp2=inv_kperp2,
                         kwork_filter=kwork, box=linked),
        "le_grids": dict(w=w, wl=wl, forbid=forbid, ng2=L - 1, nxi=2 * L),
        "nonlinear_terms": dict(nonlin=False),
        "species": dict(spec_z=sz, spec_dens=sd, spec_temp=st, spec_mass=sm,
                        spec_stm=np.sqrt(st / sm), spec_type=np.asarray(species_types, np.int32),
                        nonmaxw_corr=np.ones((S, E))),
        "run_parameters": dict(has_phi=has_phi, has_apar=has_apar, has_bpar=has_bpar,
                               beta=beta, tite=tite, fphi=1.0, fapar=1.0, fbpar=1.0),
        "dist_fn_arrays": dict(vpa=vpa, aj0=aj0, aj1=aj1, vperp2=vperp2),
        "dist_fn": dist_fn,
        "fields": dict(fieldopt_switch=field_option),
        "state_t0": dict(code_dt=0.05),
    }
    root = Group(name="/", path="/")
    for gname, arrays in groups.items():
        grp = Group(name=gname, path=f"/{gname}", arrays=dict(arrays))
        if gname == "dist_fn":
            grp.absent = frozenset(absent + [n for n in ("gamtot3", "awgt") if n not in arrays])
        root.groups[gname] = grp
    ex = Export(path=Path("synthetic.gs2export.nc"), root=root,
                dims=dict(theta=T, kx=X, ky=Y, egrid=E, species=S, **{"lambda": L}),
                attrs=dict(layout=layout))
    topo = F.topology_from_export(ex)
    fi = Group(name="fields_implicit", path="/fields_implicit")
    Mmax, Nmax = int(M_class.max()), int(N_class.max())
    cls_ik = -np.ones((len(lens), Nmax, Mmax), np.int32)
    cls_it = -np.ones((len(lens), Nmax, Mmax), np.int32)
    for c in range(len(lens)):
        M, N = topo.sizes(c)
        cls_ik[c, :N, :M] = topo.ik[c][None, :] + 1
        cls_it[c, :N, :M] = topo.it[c].T + 1
    fi.arrays.update(cls_ik=cls_ik, cls_it=cls_it, nfield=nf, nidx=nf * T)
    if with_aminv:
        for c in range(len(lens)):
            M, N = topo.sizes(c)
            R = N * nf * T
            mat = rng.standard_normal((M, R, R)) + 1j * rng.standard_normal((M, R, R))
            fi.arrays[f"aminv_c{c + 1}"] = mat / R
    root.groups["fields_implicit"] = fi
    return ex


def synth_state(ex, seed=1):
    rng = np.random.default_rng(seed)
    T, X, Y = ex.nth, ex.ntheta0, ex.naky
    g = rng.standard_normal(ex.g_shape) + 1j * rng.standard_normal(ex.g_shape)
    fl = [rng.standard_normal((Y, X, T)) + 1j * rng.standard_normal((Y, X, T)) for _ in range(3)]
    return g, Fields(*fl)


@lru_cache(maxsize=None)
def _synth(key):
    kw = dict(key)
    ex = make_synth(**kw)
    params = F.FieldParams.from_export(ex)
    data = F.FieldData.from_export(ex, params)
    return ex, params, data


def synth(**kw):
    return _synth(tuple(sorted(kw.items())))


# --------------------------------------------------------------------------- references
def ref_integrate(gw, wgt, ex, filt):
    """Literal integrate_species_original: per iglo in layout order, total += fac2*(g1+g2)."""
    S, E, L, Y, X, _, T = gw.shape
    w = np.asarray(ex["/le_grids/w"])
    wl = np.asarray(ex["/le_grids/wl"])
    tracer = np.asarray(ex["/species/spec_type"]) == F.TRACER_SPECIES
    wgt = np.where(tracer, 0.0, wgt)
    order = F._velocity_order(ex)
    tot = np.zeros((Y, X, T), complex)
    for v in order:
        is_, ie, il = np.unravel_index(v, (S, E, L))
        fac2 = (wgt[is_] * w[is_, ie]) * wl[il]
        s = gw[is_, ie, il, :, :, 0, :] + gw[is_, ie, il, :, :, 1, :]
        tot = tot + rm(s, fac2)
    tot[filt] = 0.0
    return tot


def ref_getan(g, ex, params):
    a = ex["/dist_fn_arrays"]
    sp = ex["/species"]
    z, dens, temp, mass = (np.asarray(sp[n]) for n in ("spec_z", "spec_dens", "spec_temp",
                                                        "spec_mass"))
    filt = np.asarray(ex["/kt_grids/kwork_filter"])
    aj0, aj1, vperp2, vpa = a["aj0"][..., None, :], a["aj1"][..., None, :], \
        a["vperp2"][..., None, :], a["vpa"]
    out = []
    shape = g.shape[3:5] + g.shape[-1:]
    if params.has_phi:
        out.append(ref_integrate(rm(g, aj0), z * dens, ex, filt))
    else:
        out.append(np.zeros(shape, complex))
    if params.has_apar:
        wgt = (((2.0 * params.beta) * z) * dens) * np.sqrt(temp / mass)
        out.append(ref_integrate(rm(g, np.multiply(aj0, vpa)), wgt, ex, filt))
    else:
        out.append(np.zeros(shape, complex))
    if params.has_bpar:
        out.append(ref_integrate(rm(g, np.multiply(aj1, vperp2)), temp * dens, ex, filt))
    else:
        out.append(np.zeros(shape, complex))
    if params.esv:
        itright = np.asarray(ex["/dist_fn/itright"])
        for arr in out:
            for ik in range(itright.shape[0]):
                for it in range(itright.shape[1]):
                    if itright[ik, it] >= 1:
                        arr[ik, itright[ik, it] - 1, 0] = arr[ik, it, -1]
    return tuple(out)


def ref_klein(a):
    s = c1 = c2 = 0.0
    for x in a:
        t = s + x
        lc1 = (s - t) + x if abs(s) >= abs(x) else (x - t) + s
        s = t
        t = c1 + lc1
        lc2 = (c1 - t) + lc1 if abs(c1) >= abs(lc1) else (lc1 - t) + c1
        c1 = t
        c2 = c2 + lc2
    return (s + c1) + c2


def ref_fl_avg(antot, ex, params):
    th = ex["/theta_grid"]
    theta, jacob = np.asarray(th["theta"]), np.asarray(th["jacob"])
    gam = np.asarray(ex["/dist_fn/gamtot"])
    awgt = np.asarray(ex["/dist_fn/awgt"])
    filt = np.asarray(ex["/kt_grids/kwork_filter"])
    Y, X, T = antot.shape
    out = np.zeros((Y, X), complex)
    for ik in range(Y):
        for it in range(X):
            if filt[ik, it]:
                continue
            y = (jacob * antot[ik, it].real) / gam[ik, it] + 1j * ((jacob * antot[ik, it].imag)
                                                                    / gam[ik, it])
            v = rm(y[:-1] + y[1:], theta[1:] - theta[:-1])
            tint = complex(ref_klein(v.real), ref_klein(v.imag))
            tint = complex(tint.real * 0.5, tint.imag * 0.5)
            out[ik, it] = complex((params.tite * tint.real) * awgt[ik, it],
                                  (params.tite * tint.imag) * awgt[ik, it])
    return out


def ref_getfieldeq(fl3, ant, ex, params):
    phi, apar, bpar = (np.asarray(x) for x in fl3)
    antot, antota, antotp = ant
    d = ex["/dist_fn"]
    gamtot, gamtot1, gamtot2, gridfac1 = (np.asarray(d[n]) for n in
                                          ("gamtot", "gamtot1", "gamtot2", "gridfac1"))
    kperp2 = np.asarray(ex["/kt_grids/kperp2"])
    bmag = np.asarray(ex["/theta_grid/bmag"])
    fl = ref_fl_avg(antot, ex, params) if params.fl_avg_on else np.zeros(antot.shape[:2], complex)
    fe = (antot + rm(bpar, gamtot1)) - rm(phi, gamtot * gridfac1) + fl[..., None]
    fea = antota - rm(apar, kperp2 * gridfac1)
    fep = rm((antotp + rm(bpar, gamtot2)) + rm(rm(phi, 0.5), gamtot1),
             (params.beta * params.apfac) / (bmag * bmag)) + rm(bpar, gridfac1)
    return fe, fea, fep


def ref_getfield(fl, mats, topo, nf, nth):
    """Literal getfield loops (fields_implicit.f90:261-289)."""
    nidx = nf * nth
    out = np.zeros(fl.shape, complex)
    for c, mat in enumerate(mats):
        M, N = topo.sizes(c)
        for m in range(M):
            ik = topo.ik[c][m]
            vecs = [np.concatenate([fl[f, ik, topo.it[c][m, n]] for f in range(nf)])
                    for n in range(N)]
            for r in range(N * nidx):
                u = 0.0 + 0.0j
                for n in range(N):
                    s = 0.0 + 0.0j
                    row = mat[m, r, n * nidx:(n + 1) * nidx]
                    for k in range(nidx):
                        p = cm(row[k], vecs[n][k])
                        s = complex(s.real + p.real, s.imag + p.imag)
                    u = complex(u.real - s.real, u.imag - s.imag)
                n_o, f_o, igp = r // nidx, (r // nth) % nf, r % nth
                out[f_o, ik, topo.it[c][m, n_o], igp] = u
    return out


def ref_gauss_jordan(AMT, zonal_col, nf, nth, N):
    """Literal init_inverse_matrix loops (fields_implicit.f90:760-908) in K order."""
    AMT = np.array(AMT, complex)
    M, n, _ = AMT.shape
    AIT = np.zeros((M, n, n), complex)
    for m in range(M):
        AIT[m] = np.eye(n)
    skip = RS.skip_set(N, nf, nth)
    for i in range(n):
        if skip[i]:
            continue
        L = AMT[:, i, :].copy()
        Rr = AIT[:, i, :].copy()
        for m in range(M):
            for j in range(n):
                if skip[j]:
                    continue
                if not zonal_col[m, j]:
                    fac = sdiv(AMT[m, j, i], L[m, i])
                    AMT[m, j, i] = fac
                    others = np.arange(n) != i
                    AMT[m, j, others] = AMT[m, j, others] - cm(L[m, others], fac)
                    if j == i:
                        AIT[m, j, :] = sdiv(AIT[m, j, :], L[m, i])
                    else:
                        AIT[m, j, :] = AIT[m, j, :] - sdiv(cm(Rr[m, :], L[m, j]), L[m, i])
                else:
                    AIT[m, j, :] = 0.0
    AIT[:, :, skip] = 0.0
    for j in np.nonzero(skip)[0]:
        AIT[:, j, :] = AIT[:, j + nf, :]
    return AIT


# =========================================================================== synthetic tests
def test_velocity_order_layouts():
    S, E, L = 2, 3, 4
    lx = F.velocity_order_from_layout("lxyes", S, E, L)
    np.testing.assert_array_equal(lx, np.arange(S * E * L))        # il innermost, is outer
    for lay in ("lyxes", "lexys", "xyles", "yxles"):
        np.testing.assert_array_equal(F.velocity_order_from_layout(lay, S, E, L), lx)
    yx = F.velocity_order_from_layout("yxels", S, E, L)
    expect = [np.ravel_multi_index((s, e, l), (S, E, L))
              for s in range(S) for l in range(L) for e in range(E)]
    np.testing.assert_array_equal(yx, expect)


def test_topology_tables():
    ex, params, data = synth()
    topo = data.topo
    assert [topo.sizes(c) for c in range(topo.nclass)] == [(8, 1), (2, 2), (1, 3)]
    # chains follow itright; every cell exactly once
    itright = np.asarray(ex["/dist_fn/itright"])
    seen = set()
    for c in range(topo.nclass):
        M, N = topo.sizes(c)
        for m in range(M):
            ik = topo.ik[c][m]
            for n in range(N):
                seen.add((ik, topo.it[c][m, n]))
                if n + 1 < N:
                    assert itright[ik, topo.it[c][m, n]] - 1 == topo.it[c][m, n + 1]
    assert len(seen) == ex.naky * ex.ntheta0
    # recompute from links == exported cls tables; member order: it outer, ik inner
    t2 = F.topology_from_export(ex, prefer_export=False)
    for a, b in zip(topo.it + topo.ik, t2.it + t2.ik):
        np.testing.assert_array_equal(a, b)
    assert topo.zonal[0][np.where((topo.ik[0] == 0) & (topo.it[0][:, 0] == 0))[0][0], 0]
    assert sum(int(np.asarray(z).sum()) for z in topo.zonal) == 1


@pytest.mark.parametrize("fields", [(True, False, False), (True, True, True),
                                    (True, False, True), (False, True, False)])
@pytest.mark.parametrize("layout", ["lxyes", "yxels"])
def test_getan_matches_fortran_loops(fields, layout):
    ex, params, data = synth(fields=fields, layout=layout, species_types=(1, 2, 3))
    g, _ = synth_state(ex)
    with strict_mode():
        got = jax.jit(F.getan, static_argnames="params")(jnp.asarray(g), to_jax(data), params)
    ref = ref_getan(g, ex, params)
    for r, x in zip(ref, got):
        assert np.array_equal(np.asarray(x), r), maxrel(r, x)
    # filtered mode is exactly zero; tracer species contributes nothing
    assert np.all(np.asarray(got[0])[0, 0] == 0)


def test_getan_fast_mode_close():
    ex, params, data = synth()
    g, _ = synth_state(ex)
    fast = F.FieldParams(**{**params.__dict__, "exact": False})
    a = F.getan(jnp.asarray(g), to_jax(data), params)
    b = F.getan(jnp.asarray(g), to_jax(data), fast)
    for x, y in zip(a, b):
        assert maxrel(x, y) < 1e-14


def test_esv_pass():
    ex, params, data = synth(esv=True)
    assert params.esv
    g, _ = synth_state(ex)
    with strict_mode():
        got = F.getan(jnp.asarray(g), to_jax(data), params)
    ref = ref_getan(g, ex, params)
    for r, x in zip(ref, got):
        assert np.array_equal(np.asarray(x), r)
    # zonal (self-linked) modes become periodic
    assert np.all(np.asarray(got[0])[0, 1:, 0] == np.asarray(got[0])[0, 1:, -1])


def test_klein_and_trapezoid():
    rng = np.random.default_rng(3)
    a = rng.standard_normal((4, 37)) * 10.0 ** rng.integers(-8, 8, (4, 37))
    got = np.asarray(F.klein_sum(jnp.asarray(a)))
    ref = np.array([ref_klein(row) for row in a])
    assert np.array_equal(got, ref)
    assert np.array_equal(F.klein_sum_np(a), ref)


@pytest.mark.parametrize("fields", [(True, False, False), (True, True, True)])
@pytest.mark.parametrize("adiabatic", [F.ADIABATIC_DEFAULT, F.ADIABATIC_FIELDLINEAVG])
def test_getfieldeq_matches_fortran(fields, adiabatic):
    ex, params, data = synth(fields=fields, adiabatic=adiabatic, species_types=(1,))
    assert params.fl_avg_on == (adiabatic == F.ADIABATIC_FIELDLINEAVG)
    g, fl3 = synth_state(ex)
    fl3 = Fields(*[x if on else np.zeros_like(x) for x, on in zip(fl3, fields)])
    with strict_mode():
        jd = to_jax(data)
        ant = F.getan(jnp.asarray(g), jd, params)
        got = jax.jit(F.getfieldeq, static_argnames="params")(
            Fields(*map(jnp.asarray, fl3)), *ant, jd, params)
    ref = ref_getfieldeq(fl3, tuple(map(np.asarray, ant)), ex, params)
    for on, r, x in zip(fields, ref, got):
        if on:
            assert np.array_equal(np.asarray(x), r), maxrel(r, x)


@pytest.mark.parametrize("unroll_max", [64, 0])
@pytest.mark.parametrize("fields", [(True, False, False), (True, True, True), (False, True, False)])
def test_getfield_matches_fortran_loops(fields, unroll_max, monkeypatch):
    monkeypatch.setattr(F, "UNROLL_MAX", unroll_max)     # unrolled and fori_loop sums
    ex, params, data = synth(fields=fields)
    g, fl3 = synth_state(ex)
    aminv = RS.Aminv.from_export(ex)
    jd = to_jax(data)
    with strict_mode():
        fn = jax.jit(F.getfield, static_argnames="params")
        u = fn(Fields(*map(jnp.asarray, fl3)), jnp.asarray(g), to_jax(aminv), jd, params)
        fe = F.getfieldeq_from_g(Fields(*map(jnp.asarray, fl3)), jnp.asarray(g), jd, params)
    fl = np.stack([np.asarray(fe[i]) for i in params.slots])
    ref = ref_getfield(fl, aminv.mats, data.topo, params.nfield, params.nth)
    for j, i in enumerate(params.slots):
        assert np.array_equal(np.asarray(u[i]), ref[j]), maxrel(ref[j], u[i])
    for i in range(3):
        if i not in params.slots:          # absent fields untouched
            assert np.array_equal(np.asarray(u[i]), fl3[i])
    new = F.field_update(u, Fields(*map(jnp.asarray, fl3)), params)
    for i in range(3):
        assert np.array_equal(np.asarray(new[i]), np.asarray(u[i]) + fl3[i])


def test_get_init_field_and_potentials():
    for fields in [(True, False, False), (True, True, True), (False, False, True),
                   (True, False, True)]:
        ex, params, data = synth(fields=fields, adiabatic=F.ADIABATIC_FIELDLINEAVG,
                                 species_types=(1,), with_denoms=False)
        g, _ = synth_state(ex)
        jd = to_jax(data)
        with strict_mode():
            got = jax.jit(F.get_init_field, static_argnames="params")(jnp.asarray(g), jd, params)
            pot = jax.jit(F.calculate_potentials_from_nonadiabatic_dfn,
                          static_argnames="params")(jnp.asarray(g), jd, params)
            ant = tuple(np.asarray(a) for a in F.getan(jnp.asarray(g), jd, params))
        d = ex["/dist_fn"]
        gamtot, gamtot1, gamtot2 = (np.asarray(d[n]) for n in ("gamtot", "gamtot1", "gamtot2"))
        bmag = np.asarray(ex["/theta_grid/bmag"])
        beta = params.beta
        ip, ib = F.denominators(gamtot, gamtot1, gamtot2, bmag, beta, params.has_phi,
                                params.has_bpar)
        fl = ref_fl_avg(ant[0], ex, params)[..., None]
        phi = apar = bpar = np.zeros_like(ant[0])
        if params.has_phi:
            a = ant[0] + fl
            num = (rm(a, beta * gamtot2 + bmag * bmag) - rm(ant[2], beta * gamtot1)) \
                if params.has_bpar else a
            phi = rm(num, ip)
        if params.has_apar:
            apar = rm(ant[1], np.asarray(ex["/kt_grids/inv_kperp2"]))
        if params.has_bpar:
            num = (-rm(ant[2], beta * gamtot) - rm(ant[0], (beta / 2.0) * gamtot1)) \
                if params.has_phi else -rm(ant[2], beta)
            bpar = rm(num, ib)
        for r, x in zip((phi, apar, bpar), got):
            assert np.array_equal(np.asarray(x), r), maxrel(r, x)
        # calculate_potentials_from_nonadiabatic_dfn
        sp = ex["/species"]
        tot = 0.0
        for dn, z, t in zip(sp["spec_dens"], sp["spec_z"], sp["spec_temp"]):
            tot = tot + ((dn * z) * z) / t
        pw = 1.0 / (tot + params.tite)          # no electrons -> + tite
        assert params.phi_weight == pw
        if params.has_phi:
            assert np.array_equal(np.asarray(pot.phi), rm(ant[0], pw))
        if params.has_bpar:
            assert np.array_equal(np.asarray(pot.bpar), rm(ant[2], -(beta / (bmag * bmag))))


def test_smith_hoisted_equals_cdiv():
    rng = np.random.default_rng(5)
    a = rng.standard_normal(4000) + 1j * rng.standard_normal(4000)
    b = rng.standard_normal(4000) * 10.0 ** rng.integers(-3, 3, 4000) + \
        1j * rng.standard_normal(4000)
    b[:5] = [1 + 1j, -2 + 2j, 3.0 + 0j, 1j, 0.5 - 0.5j]
    with strict_mode():
        ref = np.asarray(jax.jit(numerics.cdiv)(a, b))
        got = np.asarray(jax.jit(lambda a, b: RS._smith_apply(a, *RS._smith_factors(b)))(a, b))
    assert np.array_equal(ref, got)
    for k in range(50):
        assert sdiv(a[k], b[k]) == got[k]


@pytest.mark.parametrize("nf,N", [(1, 1), (1, 3), (2, 2), (3, 1), (3, 3)])
def test_gauss_jordan_matches_fortran_loops(nf, N):
    nth = 5
    R = nf * nth * N
    M = 3
    rng = np.random.default_rng(nf * 10 + N)
    AMT = rng.standard_normal((M, R, R)) + 1j * rng.standard_normal((M, R, R)) + 4 * np.eye(R)
    zonal_col = np.zeros((M, R), bool)
    zonal_col[1] = True                         # a zonal member
    with strict_mode():
        AIT = np.asarray(jax.jit(RS.gauss_jordan, static_argnums=(2, 3, 4))(
            AMT, zonal_col, nf, nth, N))
    ref = ref_gauss_jordan(AMT, zonal_col, nf, nth, N)
    assert np.array_equal(AIT, ref), maxrel(ref, AIT)
    # the result inverts the reduced system (spec 04 section 3.14)
    skip = RS.skip_set(N, nf, nth)
    keep = ~skip
    for m in (0, 2):
        A = AMT[m].T                                     # A[k, j] = am(k, j)
        X = AIT[m]                                       # X[j, k] = (A^-1)[j, k]
        prod = X[np.ix_(keep, keep)] @ A[np.ix_(keep, keep)]
        assert np.max(np.abs(prod - np.eye(keep.sum()))) < 1e-12
        assert np.all(X[:, skip] == 0)
        for j in np.nonzero(skip)[0]:
            assert np.array_equal(X[j], X[j + nf])
    assert np.all(AIT[1] == 0)                           # zonal member zeroed


def test_reorder_k_b_roundtrip():
    for nf, N, nth in [(1, 2, 5), (3, 3, 7)]:
        R = nf * N * nth
        a = np.arange(2 * R * R).reshape(2, R, R)
        b = RS.reorder_k_to_b(a, nf, nth, N)
        np.testing.assert_array_equal(RS.reorder_b_to_k(b, nf, nth, N), a)
        kofb = RS.k_of_b(N, nf, nth)
        # B(ig, f, n) = igp + nth*(f + nf*n) <-> K = f + nf*(igp + nth*n)
        for bi in (0, 1, nth, R - 1):
            igp, f, n = bi % nth, (bi // nth) % nf, bi // (nf * nth)
            assert kofb[bi] == f + nf * (igp + nth * n)


def _synthetic_timeadv0(ex, seed=7):
    """A linear stand-in for timeadv(.., istep=0): per mode a dense theta coupling plus a
    coupling to the right-linked cell (so supercells are connected), no cross-supercell terms."""
    rng = np.random.default_rng(seed)
    S, E, L, Y, X, _, T = ex.g_shape
    K = rng.standard_normal((3, S, E, L, Y, X, 2, T, T)) * 0.3
    C = rng.standard_normal((3, S, E, L, Y, X, 2, T, T)) * 0.1
    itright = np.asarray(ex["/dist_fn/itright"])
    right = np.where(itright > 0, itright - 1, np.arange(X)[None, :])
    has_right = (itright > 0) & (np.arange(Y)[:, None] > 0)
    Kj, Cj = jnp.asarray(K + 0j), jnp.asarray(C + 0j)
    rj = jnp.asarray(right)
    hr = jnp.asarray(has_right)
    iky = jnp.arange(Y)[:, None]

    def timeadv0(fields):
        g = 0
        for f, x in enumerate(fields):
            xr = jnp.where(hr[..., None], x[iky, rj], 0)
            g = g + jnp.einsum("seluvxtk,uvk->seluvxt", Kj[f], x) + \
                jnp.einsum("seluvxtk,uvk->seluvxt", Cj[f], xr)
        return g

    return timeadv0


@pytest.mark.parametrize("fields", [(True, False, False), (True, True, True)])
def test_build_aminv_synthetic(fields):
    ex, params, data = synth(fields=fields)
    jd = to_jax(data)
    timeadv0 = _synthetic_timeadv0(ex)
    nf, nth = params.nfield, params.nth
    am = RS.build_aminv(timeadv0, jd, params, batch_size=4)
    topo = data.topo
    for c in range(topo.nclass):
        M, N = topo.sizes(c)
        R = N * nf * nth
        AMT = np.asarray(jax.jit(partial(RS.response_rows, timeadv0, params=params, c=c))(jd))
        # brute force: perturb one member at a time, same probe semantics
        js, ns, igps, fs = RS.probe_table(N, nf, nth)
        for m in range(M):
            sub_ik, sub_it = topo.ik[c][m:m + 1], topo.it[c][m:m + 1]
            for j, n, igp, f in zip(js[::5], ns[::5], igps[::5], fs[::5]):
                probe = RS.probe_fields(n, igp, f, sub_ik, sub_it,
                                        (ex.naky, ex.ntheta0, nth), params.slots)
                gg = timeadv0(probe)
                fe = F.getfieldeq(probe, *F.getan(gg, jd, params), jd, params)
                row = RS._class_rows(F.stack_present(fe, params), sub_ik, sub_it)[0]
                assert np.allclose(np.asarray(row), AMT[m, j], rtol=0, atol=1e-13)
        # aminv inverts the reduced response of non-zonal members
        skip = RS.skip_set(N, nf, nth)
        kofb = RS.k_of_b(N, nf, nth)
        AIT = RS.reorder_b_to_k(np.asarray(am.mats[c]), nf, nth, N)
        keep = ~skip
        for m in range(M):
            if np.asarray(topo.zonal[c])[m].any():
                assert np.all(AIT[m] == 0)
                continue
            A = AMT[m].T
            prod = AIT[m][np.ix_(keep, keep)] @ A[np.ix_(keep, keep)]
            assert np.max(np.abs(prod - np.eye(keep.sum()))) < 1e-11
        assert kofb.shape == (R,)


def test_build_aminv_timeadv_args():
    """Coefficients passed as traced timeadv_args give the same matrices as a closure."""
    ex, params, data = synth(fields=(True, False, True))
    jd = to_jax(data)
    t0 = _synthetic_timeadv0(ex, seed=3)
    scale = jnp.asarray(1.0)
    a = RS.build_aminv(t0, jd, params)
    b = RS.build_aminv(lambda f, s: t0(Fields(*(x * s for x in f))), jd, params,
                       timeadv_args=scale)
    for x, y in zip(a.mats, b.mats):
        assert np.array_equal(np.asarray(x), np.asarray(y))


def test_field_solve_residual_vanishes():
    """End-to-end on the synthetic system: fields + getfield increment satisfy the reduced
    field equations, gnew(fields) = G0 + timeadv0(fields) (linear response)."""
    ex, params, data = synth(fields=(True, True, True), adiabatic=F.ADIABATIC_DEFAULT)
    jd = to_jax(data)
    timeadv0 = _synthetic_timeadv0(ex)
    am = RS.build_aminv(timeadv0, jd, params)
    g0, old = synth_state(ex, seed=11)
    old = Fields(*map(jnp.asarray, old))
    gpred = jnp.asarray(g0) + timeadv0(old)
    u = F.getfield(old, gpred, am, jd, params)
    new = F.field_update(u, old, params)
    du = Fields(*(a - b for a, b in zip(new, old)))
    gcorr = gpred + timeadv0(du)
    res = F.getfieldeq_from_g(new, gcorr, jd, params)
    res0 = F.getfieldeq_from_g(old, gpred, jd, params)
    topo = data.topo
    nth = params.nth
    for c in range(topo.nclass):
        M, N = topo.sizes(c)
        for m in range(M):
            if np.asarray(topo.zonal[c])[m].any():
                continue
            for n in range(N):
                ik, it = topo.ik[c][m], topo.it[c][m, n]
                top = nth - 1 if n < N - 1 else nth      # duplicate equations are dropped
                for i in params.slots:
                    r = np.asarray(res[i])[ik, it, :top]
                    scale = np.max(np.abs(np.asarray(res0[i])))
                    assert np.max(np.abs(r)) < 1e-11 * scale


@pytest.mark.parametrize("unroll_max", [64, 0])
def test_local_layout_roundtrip_and_matvec(unroll_max, monkeypatch):
    monkeypatch.setattr(F, "UNROLL_MAX", unroll_max)
    ex, params, data = synth(fields=(True, True, True))
    nf, nth, ntgrid = params.nfield, params.nth, params.ntgrid
    topo = data.topo
    rng = np.random.default_rng(9)
    # local matrices in extended order -> B layout -> back
    loc = {}
    mats = []
    for c in range(topo.nclass):
        M, N = topo.sizes(c)
        nextend = 1 + 2 * ntgrid * N
        mb = np.zeros((M, N * nf * nth, N * nf * nth), complex)
        for m in range(M):
            r = rng.standard_normal((nf * nextend,) * 2) + 1j * rng.standard_normal((nf * nextend,) * 2)
            loc[(c, m)] = r
            mb[m] = FL.local_to_b(r, N, nf, ntgrid)
            np.testing.assert_array_equal(FL.b_to_local(mb[m], N, nf, ntgrid), r)
        mats.append(mb)
    g, fl3 = synth_state(ex)
    jd = to_jax(data)
    with strict_mode():
        fe = FL.getfieldeq_local(Fields(*map(jnp.asarray, fl3)),
                                 *F.getan(jnp.asarray(g), jd, params), jd, params)
        fl = F.stack_present(fe, params)
        got = np.asarray(jax.jit(FL.apply_local, static_argnames="params")(
            fl, tuple(map(jnp.asarray, mats)), jd.topo, params))
    fl = np.asarray(fl)
    # literal Fortran: rowblocks per cell/field, matmul(vect, data), cell and supercell sums
    ref = np.zeros_like(got)
    for c in range(topo.nclass):
        M, N = topo.sizes(c)
        nextend = 1 + 2 * ntgrid * N
        for m in range(M):
            ik = topo.ik[c][m]
            its = topo.it[c][m]
            resp = loc[(c, m)]
            empty = np.asarray(topo.empty[c])[m]
            nrow = nf * nextend
            sc = np.zeros(nrow, complex)
            for ic in range(N):
                if empty[ic]:
                    continue
                ncol = 2 * ntgrid + (1 if ic == N - 1 else 0)
                tmp = np.zeros(nrow, complex)
                for f in range(nf):
                    cl = ic * 2 * ntgrid + f * nextend
                    vect = fl[f, ik, its[ic], :ncol]
                    S = np.zeros(nrow, complex)
                    for k in range(ncol):
                        S = S + cm(vect[k], resp[:, cl + k])
                    tmp = tmp - S
                sc = sc + tmp
            iex = 0
            for f in range(nf):
                for ic in range(N):
                    bnd = 0 if ic == N - 1 else 1
                    for igp in range(nth - bnd):
                        ref[f, ik, its[ic], igp] = sc[iex]
                        iex += 1
                for ic in range(N - 1):
                    ref[f, ik, its[ic], nth - 1] = ref[f, ik, its[ic + 1], 0]
    assert np.array_equal(got, ref), maxrel(ref, got)
    assert np.all(got[:, 0, 0] == 0)                      # empty ky = kx = 0 supercell


def test_local_residual_association():
    ex, params, data = synth(fields=(True, False, True))
    g, fl3 = synth_state(ex)
    jd = to_jax(data)
    with strict_mode():
        ant = F.getan(jnp.asarray(g), jd, params)
        a = FL.getfieldeq_local(Fields(*map(jnp.asarray, fl3)), *ant, jd, params)
        b = F.getfieldeq(Fields(*map(jnp.asarray, fl3)), *ant, jd, params)
    ant = tuple(map(np.asarray, ant))
    d = ex["/dist_fn"]
    gg = np.asarray(d["gamtot"]) * np.asarray(d["gridfac1"])
    ref = (ant[0] - rm(fl3[0], gg)) + rm(fl3[2], np.asarray(d["gamtot1"]))
    assert np.array_equal(np.asarray(a[0]), ref)
    assert maxrel(b[0], a[0]) < 1e-14                    # same equation, other rounding
    assert np.array_equal(np.asarray(a[2]), np.asarray(b[2]))   # commutative sum


def test_local_aminv_from_export_synthetic():
    ex, params, data = synth(fields=(True, False, False), field_option=3, with_aminv=False)
    topo = data.topo
    nf, ntgrid, nth = 1, ex.ntgrid, ex.nth
    rng = np.random.default_rng(4)
    leftmost = FL._leftmost_from_export(ex, topo)
    order = FL.local_supercell_order(leftmost, ex.naky, ex.ntheta0)
    grp = Group(name="fields_local", path="/fields_local")
    expect = {}
    for (ik, is_), lm in order.items():
        c, m = next((c, m) for c in range(topo.nclass) for m in range(topo.sizes(c)[0])
                    if topo.ik[c][m] == ik - 1 and topo.it[c][m, 0] == lm - 1)
        if np.asarray(topo.empty[c])[m].all():
            continue                                   # zonal supercell is not dumped
        N = topo.sizes(c)[1]
        nextend = 1 + 2 * ntgrid * N
        r = rng.standard_normal((nf * nextend,) * 2) + 0j
        grp.arrays[f"resp_ik{ik}_is{is_}"] = r
        expect[(c, m)] = FL.local_to_b(r, N, nf, ntgrid)
    ex.root.groups["fields_local"] = grp
    am = RS.Aminv.from_export(ex)
    for (c, m), mat in expect.items():
        np.testing.assert_array_equal(am.mats[c][m], mat)
    c0 = 0
    zon = [m for m in range(topo.sizes(c0)[0]) if np.asarray(topo.empty[c0])[m].all()]
    assert all(np.all(am.mats[c0][m] == 0) for m in zon)


def test_jit_tracing_all_entry_points():
    ex, params, data = synth()
    jd = to_jax(data)
    g, fl3 = synth_state(ex)
    am = to_jax(RS.Aminv.from_export(ex))
    flj = Fields(*map(jnp.asarray, fl3))
    for fn, args in [(F.getan, (jnp.asarray(g),)), (F.get_init_field, (jnp.asarray(g),)),
                     (F.calculate_potentials_from_nonadiabatic_dfn, (jnp.asarray(g),)),
                     (F.getfield, (flj, jnp.asarray(g), am)),
                     (FL.getfield_local, (flj, jnp.asarray(g), am))]:
        out = jax.jit(fn, static_argnames="params")(*args, jd, params)
        jax.block_until_ready(out)
        assert all(np.all(np.isfinite(np.asarray(x))) for x in jax.tree_util.tree_leaves(out))


# =========================================================================== reference data
def _fopt(ex):
    try:
        return F.FieldParams.from_export(ex).field_option
    except Exception:
        return "implicit"


ALL = cases_where(lambda ex: True)
IMPLICIT = cases_where(lambda ex: _fopt(ex) == "implicit")
LOCAL = cases_where(lambda ex: _fopt(ex) == "local")


@lru_cache(maxsize=4)
def _setup(case):
    ex = export_of(case)
    params = F.FieldParams.from_export(ex)
    data = F.FieldData.from_export(ex, params)
    return ex, params, data, to_jax(data), Dims.from_export(ex)


def _jit(fn):
    return jax.jit(fn, static_argnames="params")


def _time(grp):
    try:
        return float(grp["code_dt"])
    except Exception:
        return math.nan


def _aminv_for_dt(ex, code_dt, cache):
    """Exported aminv for this code_dt: initial, or a dt-ladder rung with the same bits."""
    if "init" not in cache:
        try:
            cache["init"] = RS.Aminv.from_export(ex)
        except KeyError as exc:
            pytest.skip(f"no exported response matrix (write_response off?): {exc}")
        cache["ladder"] = {}
        ladder = ex.root.groups.get("dt_ladder")
        if ladder is not None:
            for name, grp in ladder.groups.items():
                if name.startswith("k") and "code_dt" in grp.attrs:
                    cache["ladder"][float(grp.attrs["code_dt"])] = int(name[1:])
    init = cache["init"]
    if not math.isfinite(code_dt) or code_dt == init.code_dt or math.isnan(init.code_dt):
        return init
    rung = cache["ladder"].get(code_dt)
    if rung is None:
        return None
    key = f"k{rung}"
    if key not in cache:
        cache[key] = RS.Aminv.from_export(ex, rung=rung)
    return cache[key]


def _old_fields(tr, rec_t, dims):
    """phi/apar/bpar of the step (copied at fields_implicit.f90:338): dumped by the
    after_timeadv record itself (gs2_export_trace.fpp tr_point), else by advance_end."""
    grp = tr.load(rec_t)
    if "phi" in grp:
        return fields_from_group(grp, dims, new=False)
    end = tr.next(rec_t, "advance_end")
    if end is not None and "phi" in tr.load(end):
        return fields_from_group(tr.load(end), dims, new=False)
    return None


@pytest.mark.parametrize("case", ALL)
def test_class_tables_vs_export(case):
    ex = export_of(case)
    if F._opt(ex, "/fields_implicit/cls_ik") is None:
        pytest.skip("no exported cls tables")
    a = F.topology_from_export(ex, prefer_export=True)
    b = F.topology_from_export(ex, prefer_export=False)
    for x, y in zip(a.ik + a.it, b.ik + b.it):
        np.testing.assert_array_equal(x, y)


@pytest.mark.parametrize("case", ALL)
def test_fieldeq_coefficients_vs_export(case):
    """init_fieldeq recomputed with our integrator (gamtot*, denominators, awgt) vs export."""
    ex, params, data, jd, dims = _setup(case)
    chk = F.init_fieldeq_check(ex, data, params)
    for name in ("gamtot", "gamtot1", "gamtot2"):
        assert_close(ex[f"/dist_fn/{name}"], chk[name], RTOL_KERNEL, f"{case} {name}")
    g3 = F._opt(ex, "/dist_fn/gamtot3")
    if g3 is not None and chk["gamtot3"] is not None:
        assert_close(g3, chk["gamtot3"], RTOL_KERNEL, f"{case} gamtot3")
    d = ex["/dist_fn"]
    ip, ib = F.denominators(d["gamtot"], d["gamtot1"], d["gamtot2"], ex["/theta_grid/bmag"],
                            params.beta, params.has_phi, params.has_bpar)
    for name, val in (("inv_phi_denominator_g", ip), ("inv_bpar_denominator_g", ib)):
        ref = F._opt(ex, f"/dist_fn/{name}")
        if ref is not None:
            assert_close(ref, val, 0.0, f"{case} {name}")
    aw = F._opt(ex, "/dist_fn/awgt")
    if aw is not None and params.fl_avg_on:
        g3 = g3 if g3 is not None else F.gamtot3_from_gamtot(d["gamtot"], params.tite)
        assert_close(aw, F.awgt_host(ex["/theta_grid/theta"], ex["/theta_grid/jacob"], g3,
                                     ex["/kt_grids/aky"]), 0.0, f"{case} awgt")


@pytest.mark.parametrize("case", ALL)
def test_get_init_field_t0(case):
    """get_init_field(gnew) vs the t = 0 fields (computed before initialize_diagnostics, so the
    pre-diagnostics gnew is the exact input when exported)."""
    ex, params, data, jd, dims = _setup(case)
    grp = F._opt(ex, "/state_t0_prediag") or ex["/state_t0"]
    got = _jit(F.get_init_field)(jnp.asarray(grp["gnew"]), jd, params)
    ref = fields_from_group(grp, dims)
    for name, i in zip(("phi", "apar", "bpar"), range(3)):
        if i in params.slots:
            err = assert_close(ref[i], got[i], RTOL_KERNEL, f"{case} t0 {name}")
            print(f"{case} t0 {name}: {err:.2e}")


@pytest.mark.parametrize("case", ALL)
def test_get_init_field_after_reset(case):
    tr = trace_of(case)
    ex, params, data, jd, dims = _setup(case)
    resets = tr.find("after_reset_time_step")
    if not resets:
        pytest.skip("no dt reset in the trace")
    checked = 0
    for rec in resets:
        nxt = [r for r in tr.records if r.tag == "advance_begin" and
               ((r.istep == rec.istep and r.attempt > rec.attempt) or r.istep == rec.istep + 1)]
        if not nxt:
            continue
        grp = tr.load(nxt[0])
        got = _jit(F.get_init_field)(jnp.asarray(grp["gnew"]), jd, params)
        ref = fields_from_group(grp, dims)
        for i in params.slots:
            err = assert_close(ref[i], got[i], RTOL_KERNEL, f"{case} reset@{rec.istep} f{i}")
            print(f"{case} reset@{rec.istep} field {i}: {err:.2e}")
        checked += 1
    if not checked:
        pytest.skip("reset steps not traced")


def _field_solve_trace(case, local):
    tr = trace_of(case)
    ex, params, data, jd, dims = _setup(case)
    cache = {}
    fn = _jit(FL.getfield_local if local else F.getfield)
    upd = FL.field_update_local if local else F.field_update
    nsteps = 0
    worst = 0.0
    for istep in tr.steps:
        try:
            rec_g = tr.record("after_getfield", istep)
        except ExportError:
            continue
        rec_t = tr.previous(rec_g, "after_timeadv")
        if rec_t is None:
            continue
        old = _old_fields(tr, rec_t, dims)
        if old is None:
            continue
        aminv = _aminv_for_dt(ex, _time(tr.load(rec_t)), cache)
        if aminv is None:
            print(f"{case} istep={istep}: no exported aminv for this dt; skipped")
            continue
        gnew = jnp.asarray(tr.load(rec_t)["gnew"])
        oldj = Fields(*map(jnp.asarray, old))
        got = fn(oldj, gnew, to_jax(aminv), jd, params)
        ref = fields_from_group(tr.load(rec_g), dims)
        for i in params.slots:
            err = assert_close(ref[i], got[i], RTOL_KERNEL, f"{case} istep={istep} u{i}")
            worst = max(worst, err)
        rec_u = tr.next(rec_g, "after_field_update")
        if rec_u is not None:
            refu = fields_from_group(tr.load(rec_u), dims)
            new_ref = upd(Fields(*map(jnp.asarray, ref)), oldj, params)
            new_got = upd(got, oldj, params)
            for i in params.slots:
                assert_close(refu[i], new_ref[i], 0.0, f"{case} istep={istep} update{i}")
                err = assert_close(refu[i], new_got[i], RTOL_KERNEL,
                                   f"{case} istep={istep} new{i}")
                worst = max(worst, err)
        nsteps += 1
    if not nsteps:
        pytest.skip("no after_getfield records")
    print(f"{case}: {nsteps} steps, max error {worst:.2e}")


@pytest.mark.parametrize("case", IMPLICIT)
def test_getfield_trace(case):
    """getfield increment vs after_getfield (input: gnew of after_timeadv pass 1 and the old
    fields), and the field update vs after_field_update, every traced step."""
    _field_solve_trace(case, local=False)


@pytest.mark.parametrize("case", LOCAL)
def test_getfield_local_trace(case):
    _field_solve_trace(case, local=True)


@pytest.mark.parametrize("case", IMPLICIT)
def test_gauss_jordan_vs_exported_am(case):
    """init_inverse_matrix on the exported am_c<c> (probe P-R1) vs the exported aminv_c<c>.

    The current exporter cannot provide am_c<c> (a probe inside init_response_matrix changes
    the inlined Gauss-Jordan code generation; gs2_export_static.fpp export_fields_implicit),
    so this skips; build_aminv is checked against aminv_c<c> by test_build_aminv_vs_export."""
    ex, params, data, jd, dims = _setup(case)
    topo = data.topo
    if all(RS.amt_from_export(ex, c) is None for c in range(topo.nclass)):
        pytest.skip("am_c<c> not exported (P-R1 unavailable in the exporter)")
    aminv = RS.Aminv.from_export(ex)
    nf, nth = params.nfield, params.nth
    found = False
    for c in range(topo.nclass):
        AMT = RS.amt_from_export(ex, c)
        if AMT is None:
            continue
        found = True
        with strict_mode():
            mat = np.asarray(jax.jit(RS.invert_response, static_argnums=(2, 3, 4))(
                AMT, jd.topo, c, nf, nth))
        err = maxrel(aminv.mats[c], mat)
        cond = RS.condition_numbers(aminv.mats[c], topo, c, nf, nth)
        print(f"{case} class {c + 1}: max-abs-rel {err:.2e}, cond max {np.nanmax(cond):.2e}")
        assert err <= max(RTOL_KERNEL, 100 * EPS * np.nanmax(cond))
    if not found:
        pytest.skip("am_c<c> not exported (trace_response_build off)")


@pytest.mark.parametrize("case", cases_where(lambda ex: bool(np.asarray(
    F._first(ex, "/collisions/split_collisions", default=False)).reshape(-1)[0])))
def test_potentials_split_collisions(case):
    """collisions_advance (dist_fn.fpp:4130-4167): g_adjust(from) of the post-collision gnew,
    calculate_potentials_from_nonadiabatic_dfn, vs the fields at advance_end."""
    from gs2gpu.gadjust import FROM_G, TO_G, GAdjustData, g_adjust_direction
    tr = trace_of(case)
    ex, params, data, jd, dims = _setup(case)
    ga = to_jax(GAdjustData.from_export(ex))
    rp = ex["/run_parameters"]
    fphi, fbpar = float(rp["fphi"]), float(rp["fbpar"])
    adjust = bool(np.asarray(F._first(ex, "/collisions/adjust", default=True)).reshape(-1)[0])
    n = 0
    for istep in tr.steps:
        hrecs = [r for r in tr.find("coll_h_out", istep=istep, attempt=tr.last_attempt(istep))
                 if r.pass_ == 2]
        end = tr.find("advance_end", istep=istep, attempt=tr.last_attempt(istep))
        upd = tr.find("after_field_update", istep=istep, attempt=tr.last_attempt(istep))
        if not (hrecs and end and upd):
            continue
        cur = fields_from_group(tr.load(upd[-1]), dims)
        zr = tr.find("after_remove_zonal_flows", istep=istep, attempt=tr.last_attempt(istep))
        if zr:
            cur = fields_from_group(tr.load(zr[-1]), dims)
        h = jnp.asarray(tr.load(hrecs[-1])["h"])       # tr_g writes coll_h_* as 'h'
        g1 = g_adjust_direction(h, cur.phi, cur.bpar, TO_G, fphi, fbpar, ga) if adjust else h
        g2 = g_adjust_direction(g1, cur.phi, cur.bpar, FROM_G, fphi, fbpar, ga)
        got = _jit(F.calculate_potentials_from_nonadiabatic_dfn)(g2, jd, params)
        ref = fields_from_group(tr.load(end[-1]), dims)
        for i in params.slots:
            err = assert_close(ref[i], got[i], RTOL_KERNEL, f"{case} istep={istep} f{i}")
            print(f"{case} istep={istep} field {i}: {err:.2e}")
        n += 1
    if not n:
        pytest.skip("needed probes (coll_h_out pass 2, after_field_update, advance_end) absent")


# --------------------------------------------------------------------------- rebuild vs GS2
def _gs2_timeadv0(ex, *, rung=None, coeffs=None, coll_coeffs=None):
    """``timeadv(phi=0, apar=0, bpar=0, phinew, aparnew, bparnew, istep=0)`` from ``g = 0``
    (init_response_row, fields_implicit.f90:553-557, 669), assembled from the B2 source, B3
    invert_rhs/linked/parity, B6 collisions and B7 hyper kernels in the order of spec 10
    section C (response-matrix paragraph). A test stand-in for the integration layer's
    timeadv: no explicit/NL terms, hyper iff hyper_in_initialisation, collisions iff not
    split_collisions with the probe aparnew in the drag (G-27), parity after collisions."""
    from gs2gpu import hyper as HY
    from gs2gpu import invert_rhs as IR
    from gs2gpu import source as SRC
    from gs2gpu.types import Coeffs
    coeffs = to_jax(Coeffs.from_export(ex, rung=rung) if coeffs is None else coeffs)
    sparams = SRC.SourceParams.from_export(ex)
    sdata = to_jax(SRC.SourceData.from_export(ex))
    idata, iparams = IR.from_export(ex)
    idata = to_jax(idata)
    hparams = HY.HyperParams.from_export(ex)
    hdata = to_jax(HY.HyperData.from_export(ex, hparams))
    use_hyper = HY.hyper_applies(0, hparams)
    coll = None
    if "collisions" in ex.root.groups:
        from gs2gpu import collisions as CO
        from gs2gpu import collisions_init as CI
        cparams = CI.CollParams.from_export(ex)
        if cparams.colls and not cparams.split_collisions:
            cdata = CI.CollData.from_export(ex, cparams)
            cdata = to_jax(cdata._replace(init=None))
            ccoeffs = to_jax(CI.CollCoeffs.from_export(ex, rung=rung) if coll_coeffs is None
                             else coll_coeffs)
            coll = (CO, cparams, cdata, ccoeffs)
    zero_g = jnp.zeros(ex.g_shape, jnp.complex128)
    zf = jnp.zeros(ex.field_shape, jnp.complex128)
    zero_f = Fields(zf, zf, zf)

    def timeadv0(probe):
        src = SRC.get_source_term_on_local_domain(zero_g, zero_f, probe, None, (), 0, coeffs,
                                                  sdata, sparams, code_time=0.0)
        g = IR.invert_rhs(src, probe.phi, probe.bpar, coeffs, idata, iparams, gnew_in=zero_g)
        if use_hyper:
            g, _ = HY.hyper_diff(g, probe.phi, coeffs.code_dt, hdata, hparams)
        if coll is not None:
            CO, cp, cd, cc = coll
            g = CO.vspace_derivatives(g, zero_g, zero_f, probe, probe.apar, None, cc, cd, cp)
        if iparams.def_parity:
            from gs2gpu import parity as PA
            g = PA.enforce_parity(g, idata, iparams)
        return g

    return timeadv0


@pytest.mark.parametrize("case", ALL)
def test_build_aminv_vs_export(case):
    """response.build_aminv with a timeadv(istep=0) stand-in vs the exported inverse (implicit:
    the same Gauss-Jordan, local: invert_block; 00_decisions D5). Rebuilding inverts the same
    matrix, so agreement is ~cond*eps: measured error and condition numbers are reported."""
    ex, params, data, jd, dims = _setup(case)
    try:
        ref = RS.Aminv.from_export(ex)
    except KeyError as exc:
        pytest.skip(f"no exported response matrix: {exc}")
    try:
        timeadv0 = _gs2_timeadv0(ex)
    except (ImportError, AttributeError, NotImplementedError, KeyError, TypeError,
            ValueError) as exc:
        pytest.skip(f"timeadv stand-in unavailable: {type(exc).__name__}: {exc}")
    with strict_mode():
        got = RS.build_aminv(timeadv0, jd, params, batch_size=8)
    topo = data.topo
    worst = 0.0
    for c in range(topo.nclass):
        err = maxrel(ref.mats[c], got.mats[c])
        cond = RS.condition_numbers(ref.mats[c], topo, c, params.nfield, params.nth)
        cmax = float(np.nanmax(cond)) if np.isfinite(cond).any() else 1.0
        print(f"{case} class {c + 1} (M={topo.sizes(c)[0]}, N={topo.sizes(c)[1]}): "
              f"max-abs-rel {err:.2e}, cond max {cmax:.2e}, err/(cond*eps) {err / (cmax * EPS):.1f}")
        worst = max(worst, err / max(1.0, cmax * EPS * 1e3))
        assert err <= max(1e-10, 1e3 * EPS * cmax), f"class {c + 1}: {err:.2e} (cond {cmax:.2e})"


# --------------------------------------------------------------------------- harness self-test
def _fake_trace(records):
    """An in-memory ``io.trace.Trace``: ``records`` is a list of (tag, istep, attempt, pass,
    {name: array}) in execution order."""
    from gs2gpu.io.trace import Record, Trace
    tr = Trace.__new__(Trace)
    tr.records, groups, seq = [], {}, {}
    for tag, istep, attempt, pass_, arrays in records:
        k = (istep, attempt)
        seq[k] = seq.get(k, 0) + 1
        path = f"/s{istep:06d}/a{attempt}/p{seq[k]:03d}_{tag}"
        tr.records.append(Record(istep, attempt, seq[k], tag, pass_, True, path))
        groups[path] = Group(name=path, path=path, arrays=dict(arrays))
    tr._load = groups.__getitem__
    tr.path = Path("fake.gs2trace.nc")
    return tr


def test_refdata_harness_on_synthetic_trace(monkeypatch):
    """Drives the refdata tests with a synthetic export + trace generated by this module's own
    kernels, to exercise their record plumbing (not the numerics) before refdata exists."""
    import sys
    mod = sys.modules[__name__]
    ex, params, data = synth(fields=(True, True, True), adiabatic=F.ADIABATIC_DEFAULT)
    jd = to_jax(data)
    timeadv0 = _synthetic_timeadv0(ex)
    aminv = RS.build_aminv(timeadv0, jd, params, code_dt=0.05)
    for c, mat in enumerate(aminv.mats):
        ex["/fields_implicit"].arrays[f"aminv_c{c + 1}"] = np.asarray(mat)
    g0, old = synth_state(ex, seed=21)
    old = Fields(*map(jnp.asarray, old))
    zero = np.zeros(ex.field_shape, complex)
    recs = []
    for istep in (1, 2):
        gpred = np.asarray(jnp.asarray(g0) * istep + timeadv0(old))
        u = F.getfield(old, jnp.asarray(gpred), aminv, jd, params)
        new = F.field_update(u, old, params)
        t = dict(code_dt=0.05)
        oldd = dict(phi=np.asarray(old.phi), apar=np.asarray(old.apar), bpar=np.asarray(old.bpar))
        recs += [("advance_begin", istep, 1, 0, dict(gnew=gpred, phinew=oldd["phi"], **t)),
                 ("after_timeadv", istep, 1, 1, dict(gnew=gpred, phinew=oldd["phi"],
                                                     aparnew=oldd["apar"], bparnew=oldd["bpar"],
                                                     **oldd, **t)),
                 ("after_getfield", istep, 1, 1, dict(phinew=np.asarray(u.phi),
                                                      aparnew=np.asarray(u.apar),
                                                      bparnew=np.asarray(u.bpar), **t)),
                 ("after_field_update", istep, 1, 1, dict(phinew=np.asarray(new.phi),
                                                          aparnew=np.asarray(new.apar),
                                                          bparnew=np.asarray(new.bpar), **t)),
                 ("advance_end", istep, 1, 2, dict(gnew=gpred, **oldd, **t))]
    # an immediate reset at istep 2: attempt 2 restarts from get_init_field(gnew)
    recs.append(("after_reset_time_step", 2, 1, 0, dict(code_dt=0.025)))
    fi = F.get_init_field(jnp.asarray(g0), jd, params)
    recs.append(("advance_begin", 2, 2, 0, dict(gnew=np.asarray(g0), phinew=np.asarray(fi.phi),
                                                aparnew=np.asarray(fi.apar),
                                                bparnew=np.asarray(fi.bpar), code_dt=0.025)))
    tr = _fake_trace(recs)
    st0 = Group(name="state_t0", path="/state_t0",
                arrays=dict(code_dt=0.05, gnew=np.asarray(g0), phinew=np.asarray(fi.phi),
                            aparnew=np.asarray(fi.apar), bparnew=np.asarray(fi.bpar)))
    ex.root.groups["state_t0"] = st0
    monkeypatch.setattr(mod, "export_of", lambda case: ex)
    monkeypatch.setattr(mod, "trace_of", lambda case: tr)
    mod._setup.cache_clear()
    try:
        test_get_init_field_t0("synthetic")
        test_get_init_field_after_reset("synthetic")
        _field_solve_trace("synthetic", local=False)       # steps 1 (and 2 only via attempt 1)
        test_class_tables_vs_export("synthetic")
    finally:
        mod._setup.cache_clear()


@pytest.mark.parametrize("case", IMPLICIT)
def test_getfield_after_dt_change_rebuild(case):
    """Steps whose dt has no exported matrix (dt increases: the ladder only goes down, G-19):
    coefficients recomputed at the new dt (B1 coeffs; B6 for collisions), aminv rebuilt with
    response.build_aminv and the timeadv stand-in (00_decisions D3), then getfield vs
    after_getfield. The rebuild inverts the same matrix: agreement ~cond*eps (reported)."""
    tr = trace_of(case)
    ex, params, data, jd, dims = _setup(case)
    cache = {}
    todo = []
    for istep in tr.steps:
        try:
            rec_g = tr.record("after_getfield", istep)
        except ExportError:
            continue
        rec_t = tr.previous(rec_g, "after_timeadv")
        if rec_t is None:
            continue
        dt = _time(tr.load(rec_t))
        if _aminv_for_dt(ex, dt, cache) is None:
            todo.append((istep, dt, rec_t, rec_g))
    if not todo:
        pytest.skip("every traced step has an exported response matrix")
    from gs2gpu import coeffs as CF
    inputs, cparams = CF.CoeffInputs.from_export(ex), CF.CoeffParams.from_export(ex)
    rebuilt = {}
    fn = _jit(F.getfield)
    worst = 0.0
    for istep, dt, rec_t, rec_g in todo:
        if dt not in rebuilt:
            co = CF.compute_coeffs(dt, inputs, cparams)
            coll = None
            if "collisions" in ex.root.groups:
                from gs2gpu import collisions_init as CI
                cp = CI.CollParams.from_export(ex)
                if cp.colls and not cp.split_collisions:
                    cd = CI.CollData.from_export(ex, cp)
                    coll = CI.compute_coll_coeffs(dt, np.asarray(ex["/collisions/vnmult"]), cd, cp)
            t0 = _gs2_timeadv0(ex, coeffs=co, coll_coeffs=coll)
            rebuilt[dt] = RS.build_aminv(t0, jd, params, batch_size=8, code_dt=dt)
        grp = tr.load(rec_t)
        old = Fields(*map(jnp.asarray, fields_from_group(grp, dims, new=False)))
        got = fn(old, jnp.asarray(grp["gnew"]), to_jax(rebuilt[dt]), jd, params)
        ref = fields_from_group(tr.load(rec_g), dims)
        for i in params.slots:
            err = maxrel(ref[i], got[i])
            worst = max(worst, err)
            print(f"{case} istep={istep} dt={dt}: rebuilt-aminv getfield field {i}: {err:.2e}")
    assert worst < 1e-9, f"rebuild deviation {worst:.2e}"


STRICT_CASES = [c for c in ("lin_single_ae", "lin_single_em_ke_miller", "lin_box_ae_linked",
                            "lin_box_ae_fieldlineavg_lexys", "lin_single_ae_local")
                if c in ALL]


@pytest.mark.parametrize("case", STRICT_CASES or ["<no-refdata>"])
def test_strict_mode_agreement(case):
    """In numerics strict mode (no FMA contraction, no reciprocal rewrite) get_init_field and the
    implicit getfield reproduce GS2's reference build bitwise; the local matvec differs by the
    FMA-fused products of libgfortran's matmul (<= 1 ulp). Reports the number of differing
    elements; asserts a 1e-15 bound."""
    tr = trace_of(case)
    ex, params, data, jd, dims = _setup(case)
    with strict_mode():
        grp = ex["/state_t0_prediag"]
        got = _jit(F.get_init_field)(jnp.asarray(grp["gnew"]), jd, params)
        ref = fields_from_group(grp, dims)
        nd = [int(np.count_nonzero(np.asarray(got[i]) != ref[i])) for i in params.slots]
        print(f"{case} strict get_init_field: differing elements {nd}")
        for i in params.slots:
            assert maxrel(ref[i], got[i]) <= 1e-15
        aminv = to_jax(RS.Aminv.from_export(ex))
        fn = _jit(FL.getfield_local if params.field_option == "local" else F.getfield)
        for istep in tr.steps[:2]:
            rec_g = tr.record("after_getfield", istep)
            g = tr.load(tr.previous(rec_g, "after_timeadv"))
            old = Fields(*map(jnp.asarray, fields_from_group(g, dims, new=False)))
            u = fn(old, jnp.asarray(g["gnew"]), aminv, jd, params)
            ref = fields_from_group(tr.load(rec_g), dims)
            nd = [int(np.count_nonzero(np.asarray(u[i]) != ref[i])) for i in params.slots]
            print(f"{case} strict getfield istep={istep}: differing elements {nd}")
            for i in params.slots:
                assert maxrel(ref[i], u[i]) <= 1e-15


# --------------------------------------------------------------------------- local rebuild
LOCAL_HYPER_CASE = "nl_box_ae_local_hyper_dt"


@lru_cache(maxsize=4)
def _local_rebuild(case, rung):
    """Local response matrices rebuilt with GS2's fm_populate grouping at dt-ladder rung
    ``rung`` (exported rung coefficients), and the exported rung matrices."""
    ex, params, data, jd, dims = _setup(case)
    ref = FL.local_aminv_from_export(ex, rung=rung)
    t0 = _gs2_timeadv0(ex, rung=rung)
    with strict_mode():
        got = RS.build_aminv(t0, jd, params, batch_size=8, code_dt=ref.code_dt, local=True)
    return ex, params, data, jd, dims, ref, got, t0


@pytest.mark.parametrize("rung", [1, 2])
def test_local_rebuild_vs_ladder(rung):
    """(a) fm_populate grouping: the rebuilt local matrices equal GS2's dt-ladder rung (dt/2,
    dt/4 of delt = 0.4) to round-off, with amplitude-dependent hyperviscosity (iso, const_amp=F),
    where the implicit grouping would not."""
    if LOCAL_HYPER_CASE not in ALL:
        pytest.skip(f"no refdata for {LOCAL_HYPER_CASE}")
    ex, params, data, jd, dims, ref, got, t0 = _local_rebuild(LOCAL_HYPER_CASE, rung)
    assert params.field_option == "local" and params.amp_dependent_hyper
    topo = data.topo
    worst = 0.0
    for c in range(topo.nclass):
        err = maxrel(ref.mats[c], got.mats[c])
        cond = RS.condition_numbers(ref.mats[c], topo, c, params.nfield, params.nth)
        cmax = float(np.nanmax(cond))
        worst = max(worst, err)
        print(f"rung {rung} (dt={ref.code_dt}) class {c + 1} (M={topo.sizes(c)[0]}, "
              f"N={topo.sizes(c)[1]}): max-abs-rel {err:.2e}, cond max {cmax:.2e}")
        assert err <= max(1e-13, 1e3 * EPS * cmax)
    # the implicit grouping (previous behaviour) is far off with this hyperviscosity
    with strict_mode():
        old = [np.asarray(RS.invert_response(
            RS.response_rows(t0, jd, params, c, batch_size=8, local=True), jd.topo, c,
            params.nfield, params.nth)) for c in range(topo.nclass)]
    dev = max(maxrel(ref.mats[c], old[c]) for c in range(topo.nclass))
    print(f"rung {rung}: implicit grouping would deviate by {dev:.2e}")
    assert dev > 1e3 * max(worst, EPS)


def test_local_rebuild_getfield_trace():
    """(b) getfield_local with the fm_populate rebuild at dt = 0.1 vs after_getfield and
    after_field_update for the traced steps (all at dt = 0.1 after the step-1 reset)."""
    if LOCAL_HYPER_CASE not in ALL:
        pytest.skip(f"no refdata for {LOCAL_HYPER_CASE}")
    tr = trace_of(LOCAL_HYPER_CASE)
    ex = export_of(LOCAL_HYPER_CASE)
    rungs = {float(g.attrs["code_dt"]): int(n[1:]) for n, g in ex.root.groups["dt_ladder"]
             .groups.items() if n.startswith("k")}
    fn = _jit(FL.getfield_local)
    n = 0
    worst = 0.0
    for istep in tr.steps:
        try:
            rec_g = tr.record("after_getfield", istep)
        except ExportError:
            continue
        rec_t = tr.previous(rec_g, "after_timeadv")
        dt = _time(tr.load(rec_t))
        if dt not in rungs:
            print(f"istep={istep}: dt={dt} is not a ladder rung; skipped")
            continue
        _, params, data, jd, dims, ref_am, got_am, _ = _local_rebuild(LOCAL_HYPER_CASE, rungs[dt])
        grp = tr.load(rec_t)
        old = Fields(*map(jnp.asarray, fields_from_group(grp, dims, new=False)))
        u = fn(old, jnp.asarray(grp["gnew"]), to_jax(got_am), jd, params)
        new = FL.field_update_local(u, old, params)
        ref_u = fields_from_group(tr.load(rec_g), dims)
        ref_new = fields_from_group(tr.load(tr.next(rec_g, "after_field_update")), dims)
        for i in params.slots:
            e1 = assert_close(ref_u[i], u[i], RTOL_KERNEL, f"istep={istep} u{i}")
            e2 = assert_close(ref_new[i], new[i], RTOL_KERNEL, f"istep={istep} new{i}")
            worst = max(worst, e1, e2)
            print(f"istep={istep} dt={dt}: increment {e1:.2e}, updated field {e2:.2e}")
        n += 1
    assert n >= 10, f"only {n} traced steps checked"
    print(f"{LOCAL_HYPER_CASE}: {n} steps, max error {worst:.2e}")


def _fortran_populate_impulses(ex, topo, nfield):
    """Literal transcription of fm_populate / fm_init_next_field_points' bookkeeping
    (fields_local.f90:1973-2167): per field, calls until pts_remain == 0; returns, per call,
    the set of (ik, it, ig) impulses and the (supercell, iex) rows stored."""
    ntgrid = ex.ntgrid
    lm = FL._leftmost_from_export(ex, topo)
    order = FL.local_supercell_order(lm, ex.naky, ex.ntheta0)
    chains = {}
    for c in range(topo.nclass):
        for m in range(topo.sizes(c)[0]):
            chains[(int(topo.ik[c][m]) + 1, int(topo.it[c][m, 0]) + 1)] = \
                [int(x) + 1 for x in topo.it[c][m]]
    scs = [(ik, is_, chains[(ik, lmost)]) for (ik, is_), lmost in sorted(order.items())]
    calls = []
    for ifl in range(1, nfield + 1):
        nextend = {k: 1 + 2 * ntgrid * len(ch) for k, (_, _, ch) in enumerate(scs)}
        initialised = {k: [False] * nextend[k] for k in nextend}
        initdone = {k: False for k in nextend}
        pts_remain = sum(nextend.values())
        while pts_remain > 0:
            imp, rows, set_here = set(), [], 0
            for k, (ik, is_, ch) in enumerate(scs):
                if initdone[k]:
                    continue
                iex = initialised[k].index(False) + 1
                if iex == nextend[k]:
                    initdone[k] = True
                initialised[k][iex - 1] = True
                set_here += 1
                ig = ntgrid if iex == nextend[k] else -ntgrid + (iex - 1) % (2 * ntgrid)
                ic = min(len(ch), 1 + (iex - 1) // (2 * ntgrid))
                it = ch[ic - 1]
                imp.add((ik, it, ig))
                if ig == -ntgrid and it != ch[0]:
                    imp.add((ik, ch[ic - 2], ntgrid))
                rows.append((k, iex))
            pts_remain -= set_here
            calls.append((ifl, imp, rows))
    return calls


@pytest.mark.parametrize("fields", [(True, False, False), (True, True, True)])
def test_local_probe_pattern_matches_fm_populate(fields):
    ex, params, data = synth(fields=fields)
    topo = data.topo
    calls = _fortran_populate_impulses(ex, topo, params.nfield)
    pmax = max(1 + 2 * ex.ntgrid * topo.sizes(c)[1] for c in range(topo.nclass))
    assert len(calls) == params.nfield * pmax
    for k, (ifl, imp, rows) in enumerate(calls):
        f, p = ifl - 1, k % pmax + 1
        pf = RS.local_probe_fields(jnp.asarray(f), jnp.asarray(p), topo, params)
        x = np.asarray(pf[params.slots[f]])
        got = {(ik + 1, it + 1, igp - ex.ntgrid) for ik, it, igp in zip(*np.nonzero(x))}
        assert got == imp, (ifl, p)
        assert np.all(x[x != 0] == 1.0)
        for i in range(3):
            if i != params.slots[f]:
                assert not np.any(np.asarray(pf[i]))
        # stored rows: iex = p of every active supercell <-> K index j(f, p)
        assert all(iex == p for _, iex in rows)
        for c in range(topo.nclass):
            N = topo.sizes(c)[1]
            active, n, igp, j = RS.local_probe_index(f, p, N, params.nfield, params.nth)
            assert bool(active) == (p <= 1 + 2 * ex.ntgrid * N)
            if bool(active):
                e = (FL.ext_index_of_b(N, params.nfield, ex.ntgrid)[0])[RS.reorder_k_to_b(
                    np.arange((N * params.nfield * params.nth) ** 2).reshape(
                        1, N * params.nfield * params.nth, -1), params.nfield, params.nth,
                    N)[0, :, 0] // (N * params.nfield * params.nth) == int(j)]
                assert int(e[0]) == f * (1 + 2 * ex.ntgrid * N) + p - 1


def test_build_aminv_local_grouping_synthetic():
    """With a mode-coupling (amplitude-dependent) stand-in, build_aminv(local=True) equals a
    brute-force construction that evaluates each fm_populate call literally."""
    ex, params, data = synth(fields=(True, False, True), field_option=3)
    assert params.field_option == "local"
    jd = to_jax(data)
    base = _synthetic_timeadv0(ex, seed=5)
    kx = jnp.asarray(np.asarray(ex["/kt_grids/akx"]))[None, :, None]
    ky = jnp.asarray(np.asarray(ex["/kt_grids/aky"]))[:, None, None]

    def coupled(f):                     # like iso_shear: global factor from sum |phi|^2 k^4
        s = jnp.sum(jnp.abs(f.phi) ** 2 * (kx ** 2 + ky ** 2) ** 2, axis=(0, 1))   # (T,)
        return base(f) * jnp.exp(-0.3 * jnp.sqrt(s))[None, None, None, None, None, None, :]

    am = RS.build_aminv(coupled, jd, params, batch_size=3)
    topo = data.topo
    nf, nth = params.nfield, params.nth
    calls = _fortran_populate_impulses(ex, topo, nf)
    pmax = max(1 + 2 * ex.ntgrid * topo.sizes(c)[1] for c in range(topo.nclass))
    AMT = [np.zeros((M, N * nf * nth, N * nf * nth), complex)
           for M, N in (topo.sizes(c) for c in range(topo.nclass))]
    for k, (ifl, imp, rows) in enumerate(calls):
        x = np.zeros(ex.field_shape, complex)
        for ik, it, ig in imp:
            x[ik - 1, it - 1, ig + ex.ntgrid] = 1.0
        fl3 = [np.zeros(ex.field_shape, complex) for _ in range(3)]
        fl3[params.slots[ifl - 1]] = x
        probe = Fields(*map(jnp.asarray, fl3))
        g = coupled(probe)
        fe = FL.getfieldeq_local(probe, *F.getan(g, jd, params), jd, params)
        fl = F.stack_present(fe, params)
        p = k % pmax + 1
        for c in range(topo.nclass):
            M, N = topo.sizes(c)
            active, _, _, j = RS.local_probe_index(ifl - 1, p, N, nf, nth)
            if bool(active):
                AMT[c][:, int(j), :] = np.asarray(RS._class_rows(fl, topo.ik[c], topo.it[c]))
    for c in range(topo.nclass):
        ref = np.asarray(RS.invert_response(jnp.asarray(AMT[c]), jd.topo, c, nf, nth))
        assert maxrel(ref, am.mats[c]) < 1e-13
    # and the implicit grouping differs for this coupled operator
    imp = RS.build_aminv(coupled, jd, params, local=False)
    assert max(maxrel(a, b) for a, b in zip(am.mats, imp.mats)) > 1e-6
