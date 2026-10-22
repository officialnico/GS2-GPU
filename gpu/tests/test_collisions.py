"""Tests of the collision operator (spec 06; owner B6).

Two families:

* **Synthetic** (always run): a small but consistent trapped-particle grid built with the le_grids
  formulas (jend/forbid, xigridset), checked against *literal* scalar transcriptions of the Fortran
  loops written here independently of the vectorised code: the g2le gather/scatter formula, the Lorentz
  and energy-diffusion matrices and eliminations, the two conservation-vector builds with the
  ``integrate_moment_r34`` iglo order, and every runtime stage of ``solfp1_le_layout``. These pin the
  batched ``lax.scan`` solves to the Fortran forward/back substitution order.
* **Reference data** (skip without ``gpu/refdata``): coefficients vs the export at the initial dt and
  every dt-ladder rung, the static le arrays, every collision probe of the trace, the whole
  ``vspace_derivatives`` vs ``after_vspace_derivatives`` and ``collisions_advance`` vs
  ``after_collisions_advance`` (split case). Tolerance ``RTOL_KERNEL`` (00_decisions D1).
"""

from __future__ import annotations

import math

import jax
import jax.numpy as jnp
import numpy as np
import pytest

from conftest import RTOL_KERNEL, assert_close, cases_where, export_of, have, trace_of
from gs2gpu import collisions as coll
from gs2gpu import collisions_init as ci
from gs2gpu.numerics import max_abs_rel
from gs2gpu.types import Fields, to_jax

EPS = np.finfo(np.float64).eps
COLL_CASES_ALL = ["lin_single_coll_ke_full", "lin_single_coll_ke_lorentz_local",
                  "lin_single_coll_ae_ediffuse", "lin_single_coll_ke_split", "nl_box_em_ke_coll"]


# ======================================================================== synthetic grid
def _capg(x):
    # capg / capg_kernel, collisions.fpp:2165-2176
    sqrt_pi = math.sqrt(math.pi)
    kern = math.erf(x) - 2 * x * math.exp(-x ** 2) / sqrt_pi
    return 2 * math.exp(-x ** 2) * kern / (x * sqrt_pi)


def synth_raw(*, ntgrid=2, ng2=2, negrid=3, nspec=2, naky=2, ntheta0=2, trapped=True,
              model=ci.MODEL_FULL, conserve_moments=True, drag=True, special_wfb=False,
              vpar_zero_mean=True, hypermult=False, hyper_colls=False, pw_wfb_zero=False,
              layout="lxyes", force_collisions=False, zero_vnew_species=None, kwork00=True,
              vnmult=(1.3, 1.1), cfac=1.0, seed=1, literal_static=True):
    """A ``raw_from_export``-style dict for a small, internally consistent grid."""
    rng = np.random.default_rng(seed)
    T = 2 * ntgrid + 1
    theta = np.linspace(-np.pi, np.pi, T)
    bmag = 1.0 / (1.0 + 0.3 * np.cos(theta))
    bmax = bmag.max()
    xg = np.linspace(0.92, 0.35, ng2)
    al = (1.0 - xg ** 2) / bmax                              # passing, increasing
    if trapped:
        al = np.concatenate([al, 1.0 / bmag[:ntgrid + 1]])  # wfb = 1/bmax, then 1/bset
    L = al.size
    arg = 1.0 - al[:, None] * bmag[None, :]                  # (L, T)
    if trapped:
        forbid = arg < -1e-8
        jend = (arg > -1e-8).sum(axis=0).astype(np.int32)
    else:
        forbid = np.zeros((L, T), bool)
        jend = np.zeros(T, np.int32)
    xi = np.zeros((T, 2 * L))
    il_map = np.zeros((T, 2 * L), np.int32)
    sg_map = np.zeros((T, 2 * L), np.int32)
    for t in range(T):
        je = int(jend[t])
        if je > 0:
            xi[t, :je] = np.sqrt(np.maximum(1.0 - al[:je] * bmag[t], 0.0))
            if je >= 2:
                xi[t, je:2 * je - 1] = -xi[t, je - 2::-1]
        for ixi in range(1, 2 * L + 1):                      # xigridset, le_grids.f90:2972-3011
            if je == 0:
                sg, il = (2, 2 * L - ixi + 1) if ixi > L else (1, ixi)
            elif ixi >= L + je + 1:
                sg, il = 2, 2 * L + je + 1 - ixi
            elif ixi >= 2 * je:
                sg, il = 1, ixi - je
            elif ixi >= je:
                sg, il = 2, 2 * je - ixi
            else:
                sg, il = 1, ixi
            il_map[t, ixi - 1], sg_map[t, ixi - 1] = il, sg
    allowed = ~forbid
    wl = np.where(allowed, rng.uniform(0.05, 0.3, (L, T)), 0.0)
    pw = np.where(allowed.T, rng.uniform(0.05, 0.3, (T, L)), 0.0)
    if trapped and pw_wfb_zero:
        for t in range(T):
            if jend[t] == ng2 + 1:
                pw[t, ng2] = 0.0                             # wfb weight 0 at its bounce point
    wxi = np.zeros((T, 2 * L))
    for t in range(T):
        wxi[t] = wl[il_map[t] - 1, t]
    nxi = max(2 * L - 1, 2 * ng2)
    nxi_lim = nxi + 1 if nxi > 2 * ng2 else nxi
    S, E, Y, X = nspec, negrid, naky, ntheta0
    energy = np.sort(rng.uniform(0.2, 4.0, E))
    speed = np.sqrt(energy)
    w_maxw = rng.uniform(0.1, 0.5, E)
    w = rng.uniform(0.1, 0.5, (S, E))
    capg_half = np.array([_capg((speed[j] + speed[j + 1]) * 0.5) for j in range(E - 1)])
    spec_type = np.array([1, 2, 1, 1][:S], np.int32)
    vnew = rng.uniform(0.01, 0.2, (S, E, Y))
    if zero_vnew_species is not None:
        vnew[zero_vnew_species] = 0.0
    kperp2 = rng.uniform(0.05, 2.0, (Y, X, T))
    kwork = np.zeros((Y, X), bool)
    kwork[0, 0] = kwork00
    aj0 = np.where(allowed[None, None, :, None, None, :], rng.uniform(0.6, 1.0, (S, E, L, Y, X, T)),
                   0.0)
    aj1 = np.where(allowed[None, None, :, None, None, :], rng.uniform(0.3, 0.5, (S, E, L, Y, X, T)),
                   0.0)
    vpa1 = np.sqrt(energy[None, :, None, None, None, None]
                   * np.maximum(arg, 0.0)[None, None, :, None, None, :])
    vpa1 = np.broadcast_to(vpa1, (S, E, L, Y, X, T))
    vpa = np.stack([vpa1, -vpa1], axis=5)
    vperp2 = np.broadcast_to(energy[None, :, None, None, None, None] * al[None, None, :, None, None, None]
                             * bmag[None, None, None, None, None, :], (S, E, L, Y, X, T)).copy()
    has_lorentz = model in (ci.MODEL_FULL, ci.MODEL_LORENTZ, ci.MODEL_LORENTZ_TEST)
    has_diffuse = model in (ci.MODEL_FULL, ci.MODEL_EDIFFUSE)
    raw = dict(
        ntgrid=ntgrid, nth=T, ntheta0=X, naky=Y, nlambda=L, negrid=E, nspec=S, ng2=ng2, nxi=nxi,
        trapped=trapped, layout=layout, fphi=1.0, fbpar=1.0, beta=0.01,
        bmag=bmag, al=al, wl=wl, jend=jend, forbid=forbid, xi=xi, wxi=wxi, ixi_to_il=il_map,
        ixi_to_isgn=sg_map, w=w, energy_maxw=energy, speed_maxw=speed, w_maxw=w_maxw,
        kperp2=kperp2, kwork_filter=kwork, ieqzip=np.zeros((Y, X), bool),
        tunits=rng.uniform(0.8, 1.2, Y),
        spec_type=spec_type, spec_vnewk=np.array([0.01, 0.1, 0.02, 0.03][:S]),
        spec_z=np.array([1.0, -1.0, 2.0, 1.0][:S]), spec_dens=np.array([1.0, 1.0, 0.3, 0.2][:S]),
        spec_stm=np.array([1.0, 60.0, 0.7, 0.5][:S]), spec_zstm=np.array([1.0, -60.0, 1.4, 0.5][:S]),
        spec_smz=np.array([1.0, -0.0165, 0.6, 1.2][:S]), spec_zt=np.array([1.0, -1.0, 2.0, 1.0][:S]),
        nonmaxw_corr=np.ones((S, E)), vpa=vpa, aj0=aj0, aj1=aj1, vperp2=vperp2,
        colls=True, collision_model_switch=model, has_lorentz=has_lorentz,
        has_diffuse=has_diffuse, drag=drag and has_lorentz, use_le_layout=True,
        conserve_moments=conserve_moments, conservative=True, resistivity=True,
        special_wfb_lorentz=special_wfb, vpar_zero_mean=vpar_zero_mean, conserve_forbid_zero=True,
        force_collisions=force_collisions, adjust=True, heating=False, split_collisions=False,
        hypermult=hypermult, hyper_colls=hyper_colls, lorentz_switch=1, ediff_switch=1, cfac=cfac,
        nxi_lim=nxi_lim, vary_vnew=False,
        vnew=vnew, vnew_s=rng.uniform(0.01, 0.2, (S, E, Y)), vnew_D=rng.uniform(0.01, 0.2, (S, E, Y)),
        vnew_E=rng.uniform(0.01, 0.2, (S, E, Y)), delvnew=rng.uniform(-0.1, 0.1, (S, E, Y)),
        vnewh=rng.uniform(0.0, 0.3, (S, Y, X, T)) if hyper_colls else None,
        pitch_weights=pw, e15=energy ** 1.5, capg_half=capg_half, vnmult=np.array(vnmult),
    )
    raw["vpdiff"] = raw["vpdiffle"] = raw["aj0le"] = raw["vperp_aj1le"] = raw["vpa_aj0_le"] = None
    if not literal_static:            # large grids: CollData builds them with the vectorised code
        return raw
    lit = Lit(raw)
    if has_lorentz:
        raw["vpdiff"] = lit.init_vpdiff()
        raw["vpdiffle"] = lit.vpdiffle(raw["vpdiff"], nxi_lim)
    raw["aj0le"], raw["vperp_aj1le"], raw["vpa_aj0_le"] = lit.init_le_bessel(nxi_lim)
    return raw


# ======================================================================== literal transcriptions
class Lit:
    """Scalar, loop-by-loop transcriptions of the Fortran (1-based indices, ig in -ntgrid..ntgrid).

    Written from src/collisions.fpp, src/le_derivatives.fpp, src/le_grids.f90 and
    src/dist_fn_arrays.fpp without reusing gs2gpu code.
    """

    def __init__(self, raw):
        self.r = raw
        self.nt = raw["ntgrid"]
        self.S, self.E, self.L = raw["nspec"], raw["negrid"], raw["nlambda"]
        self.Y, self.X, self.T = raw["naky"], raw["ntheta0"], raw["nth"]
        self.ng2 = raw["ng2"]
        self.nxi = raw["nxi"]
        self.NX = self.nxi + 1
        self.XI = raw["nxi_lim"]

    # ----- accessors
    def jend(self, ig):
        return int(self.r["jend"][ig + self.nt])

    def forbid(self, ig, il):
        return bool(self.r["forbid"][il - 1, ig + self.nt])

    def bmag(self, ig):
        return float(self.r["bmag"][ig + self.nt])

    def al(self, il):
        return float(self.r["al"][il - 1])

    def pw(self, il, ig):
        return float(self.r["pitch_weights"][ig + self.nt, il - 1])

    def xi(self, ixi, ig):
        return float(self.r["xi"][ig + self.nt, ixi - 1])

    def wxi(self, ixi, ig):
        return float(self.r["wxi"][ig + self.nt, ixi - 1])

    def ixi_to_il(self, ixi, ig):
        return int(self.r["ixi_to_il"][ig + self.nt, ixi - 1])

    def ixi_to_isgn(self, ixi, ig):
        return int(self.r["ixi_to_isgn"][ig + self.nt, ixi - 1])

    def kperp2(self, ig, it, ik):
        return float(self.r["kperp2"][ik - 1, it - 1, ig + self.nt])

    def v(self, name, ik, ie, is_):
        return float(self.r[name][is_ - 1, ie - 1, ik - 1])

    def sp(self, name, is_):
        return float(self.r[name][is_ - 1])

    def skip(self, is_, ik):
        return (abs(self.v("vnew", ik, 1, is_)) < 2.0 * EPS) and not self.r["force_collisions"]

    def filt(self, it, ik, kwork=None):
        kw = self.r["kwork_filter"] if kwork is None else kwork
        return bool(kw[ik - 1, it - 1]) or bool(self.r["ieqzip"][ik - 1, it - 1])

    def electron(self, is_):
        return int(self.r["spec_type"][is_ - 1]) in (2, 4)

    def iles(self):
        for is_ in range(1, self.S + 1):
            for ik in range(1, self.Y + 1):
                for it in range(1, self.X + 1):
                    for ig in range(-self.nt, self.nt + 1):
                        yield is_, ik, it, ig

    def iglos(self):
        """(it, ik, il, ie, is) in ascending iglo order of the layout string (first letter fastest)."""
        sizes = {"x": self.X, "y": self.Y, "l": self.L, "e": self.E, "s": self.S}
        lay = self.r["layout"]
        n = self.X * self.Y * self.L * self.E * self.S
        for iglo in range(n):
            idx, rem = {}, iglo
            for c in lay:
                idx[c] = rem % sizes[c] + 1
                rem //= sizes[c]
            yield idx["x"], idx["y"], idx["l"], idx["e"], idx["s"]

    # ----- g2le (le_grids.f90:3519-3566)
    def ixi_of(self, ig, il0, isign):
        je, nl = self.jend(ig), self.L
        il = il0
        if je == 0:
            if isign == 2:
                il = 2 * nl + 1 - il
        else:
            if il == je:
                if isign == 1:
                    il = 2 * je
            elif il > je:
                il = il + je if isign == 1 else 2 * nl + 1 - il + je
            else:
                if isign == 2:
                    il = 2 * je - il
        return il

    def gather(self, h, nrows=None):
        nrows = self.NX if nrows is None else nrows
        out = np.zeros((self.S, self.Y, self.X, self.T, self.E + 1, nrows), h.dtype)
        for is_, ik, it, ig in self.iles():
            for ie in range(1, self.E + 1):
                for isign in (1, 2):
                    for il0 in range(1, self.L + 1):
                        if self.forbid(ig, il0):
                            continue
                        ixi = self.ixi_of(ig, il0, isign)
                        out[is_ - 1, ik - 1, it - 1, ig + self.nt, ie - 1, ixi - 1] = \
                            h[is_ - 1, ie - 1, il0 - 1, ik - 1, it - 1, isign - 1, ig + self.nt]
        return out

    def scatter(self, gle, h):
        h = np.array(h, copy=True)
        for is_, ik, it, ig in self.iles():
            for ie in range(1, self.E + 1):
                for isign in (1, 2):
                    for il0 in range(1, self.L + 1):
                        if self.forbid(ig, il0):
                            continue
                        ixi = self.ixi_of(ig, il0, isign)
                        h[is_ - 1, ie - 1, il0 - 1, ik - 1, it - 1, isign - 1, ig + self.nt] = \
                            gle[is_ - 1, ik - 1, it - 1, ig + self.nt, ie - 1, ixi - 1]
        return h

    # ----- g_adjust (dist_fn_arrays.fpp:117-170)
    def g_adjust(self, g, phi, bpar, facphi, facbpar):
        g = np.array(g, copy=True)
        trapped = self.L > self.ng2
        bpar_factor = 2.0 * facbpar if facbpar != 0 else 0.0
        for is_ in range(1, self.S + 1):
            for ie in range(1, self.E + 1):
                phi_factor = (self.sp("spec_zt", is_) * self.r["nonmaxw_corr"][is_ - 1, ie - 1]) \
                    * facphi if facphi != 0 else 0.0
                for il in range(1, self.L + 1):
                    for ik in range(1, self.Y + 1):
                        for it in range(1, self.X + 1):
                            for ig in range(-self.nt, self.nt + 1):
                                if trapped and il > self.jend(ig):
                                    continue
                                t = ig + self.nt
                                i6 = (is_ - 1, ie - 1, il - 1, ik - 1, it - 1, t)
                                adj = bpar_factor * self.r["vperp2"][i6] * self.r["aj1"][i6] \
                                    * bpar[ik - 1, it - 1, t] \
                                    + phi_factor * phi[ik - 1, it - 1, t] * self.r["aj0"][i6]
                                for s in (0, 1):
                                    g[is_ - 1, ie - 1, il - 1, ik - 1, it - 1, s, t] += adj
        return g

    # ----- static le arrays
    def init_vpdiff(self):
        vp = np.zeros((self.L, 2, self.T))
        ng2 = self.ng2
        sq = lambda x: math.sqrt(max(0.0, x))  # noqa: E731
        for ig in range(-self.nt, self.nt + 1):
            je = self.jend(ig)
            te = ng2 if (je <= ng2 or (je == ng2 + 1 and self.r["special_wfb_lorentz"])) else je
            b = self.bmag(ig)
            for il in range(2, te):
                slb0, slb1, slb2 = sq(1.0 - b * self.al(il - 1)), sq(1.0 - b * self.al(il)), \
                    sq(1.0 - b * self.al(il + 1))
                slbl = (slb1 + slb0) * 0.5
                slbr = (slb1 + slb2) * 0.5
                vp[il - 1, 0, ig + self.nt] = (slbl * slbl - slbr * slbr) / self.pw(il, ig)
            slb1, slb2 = sq(1.0 - b * self.al(1)), sq(1.0 - b * self.al(2))
            slbr = 0.5 * (slb1 + slb2)
            vp[0, 0, ig + self.nt] = (1.0 - slbr * slbr) / self.pw(1, ig)
            il = te
            slb0 = sq(1.0 - b * self.al(il - 1))
            if te == ng2:
                slb1 = sq(1.0 - b * self.al(il))
                slb2 = -slb1
            else:
                slb1 = 0.0
                slb2 = -slb0
            slbl = (slb1 + slb0) * 0.5
            slbr = (slb1 + slb2) * 0.5
            if abs(self.pw(il, ig)) <= EPS:
                vp[il - 1, 0, ig + self.nt] = 0.0
            else:
                vp[il - 1, 0, ig + self.nt] = (slbl * slbl - slbr * slbr) / self.pw(il, ig)
            vp[:, 1, ig + self.nt] = -vp[:, 0, ig + self.nt]
        return vp

    def vpdiffle(self, vpdiff, nxi_lim):
        out = np.zeros((self.T, nxi_lim))
        for ig in range(-self.nt, self.nt + 1):
            for ixi in range(1, nxi_lim + 1):
                il, isgn = self.ixi_to_il(ixi, ig), self.ixi_to_isgn(ixi, ig)
                out[ig + self.nt, ixi - 1] = vpdiff[il - 1, isgn - 1, ig + self.nt]
        return out

    def init_le_bessel(self, XI):
        z = self.r["aj0"] + 1j * self.r["aj1"]
        zb = np.stack([z, z], axis=5)
        ctmp = self.gather(zb, nrows=XI)
        aj0le = ctmp.real.copy()
        vperp = ctmp.imag.copy()
        for is_, ik, it, ig in self.iles():
            for ie in range(1, self.E + 1):
                for ixi in range(1, XI + 1):
                    il = self.ixi_to_il(ixi, ig)
                    i = (is_ - 1, ik - 1, it - 1, ig + self.nt, ie - 1, ixi - 1)
                    vperp[i] = vperp[i] * self.r["energy_maxw"][ie - 1] * self.al(il)
        vpale = self.gather(self.r["vpa"].astype(complex), nrows=XI).real
        return aj0le, vperp, vpale * aj0le

    # ----- matrices
    def lorentz_matrix(self, ig, ik, it, ie, is_, code_dt, vnm1):
        """get_lorentz_matrix (collisions.fpp:2392-2639): aa, bb, cc (1-based lists, length nxi_lim+1)."""
        r, ng2 = self.r, self.ng2
        je = self.jend(ig)
        special = r["special_wfb_lorentz"]
        if je <= ng2 or (je == ng2 + 1 and special):
            te, te2, teh = ng2, 2 * ng2, ng2
        else:
            te, te2, teh = je, 2 * je - 1, je - 1
        if r["collision_model_switch"] == ci.MODEL_LORENTZ_TEST:
            vn = vnm1 * abs(self.sp("spec_vnewk", is_)) * float(r["tunits"][ik - 1])
        else:
            vhyp = float(r["vnewh"][is_ - 1, ik - 1, it - 1, ig + self.nt]) if r["hyper_colls"] \
                else 0.0
            if r["hypermult"]:
                vn = vnm1 * self.v("vnew", ik, ie, is_) * (1 + vhyp)
            else:
                vn = vnm1 * self.v("vnew", ik, ie, is_) + vhyp
        n = self.XI
        aa, bb, cc = [0.0] * (n + 1), [0.0] * (n + 1), [0.0] * (n + 1)
        b = self.bmag(ig)
        sq = lambda x: math.sqrt(max(0.0, x))  # noqa: E731
        E_ = float(r["energy_maxw"][ie - 1])
        kp = self.kperp2(ig, it, ik)
        zstm = self.sp("spec_zstm", is_)
        cfac = r["cfac"]

        def eef(slb1):
            x = b * zstm
            return 0.5 * E_ * (1 + slb1 * slb1) * kp * cfac / (x * x)

        for il in range(2, te):
            slb0, slb1, slb2 = sq(1.0 - b * self.al(il - 1)), sq(1.0 - b * self.al(il)), \
                sq(1.0 - b * self.al(il + 1))
            slbl = (slb1 + slb0) * 0.5
            slbr = (slb1 + slb2) * 0.5
            ee = eef(slb1)
            cc[il] = 2.0 * vn * code_dt * (1 - slbr * slbr) / (self.pw(il, ig) * (slb2 - slb1))
            aa[il] = 2.0 * vn * code_dt * (1 - slbl * slbl) / (self.pw(il, ig) * (slb1 - slb0))
            bb[il] = 1.0 - (aa[il] + cc[il]) + ee * vn * code_dt
        slb1, slb2 = sq(1.0 - b * self.al(1)), sq(1.0 - b * self.al(2))
        slbr = (slb1 + slb2) * 0.5
        ee = eef(slb1)
        cc[1] = 2.0 * vn * code_dt * (1.0 - slbr * slbr) / (self.pw(1, ig) * (slb2 - slb1))
        aa[1] = 0.0
        bb[1] = 1.0 - cc[1] + ee * vn * code_dt
        il = te
        slb0 = sq(1.0 - b * self.al(il - 1))
        if te == ng2:
            slb1 = sq(1.0 - b * self.al(il))
            slb2 = -slb1
        else:
            slb1 = 0.0
            slb2 = -slb0
        slbl = (slb1 + slb0) * 0.5
        slbr = (slb1 + slb2) * 0.5
        ee = eef(slb1)
        deltaxi = self.pw(il, ig)
        if te == je:
            deltaxi = 2 * deltaxi
        if (not special) and deltaxi == 0 and il == ng2 + 1:
            deltaxi = self.xi(ng2, ig) - self.xi(ng2 + 2, ig)
        cc[il] = 2.0 * vn * code_dt * (1 - slbr * slbr) / (deltaxi * (slb2 - slb1))
        aa[il] = 2.0 * vn * code_dt * (1 - slbl * slbl) / (deltaxi * (slb1 - slb0))
        bb[il] = 1.0 - (aa[il] + cc[il]) + ee * vn * code_dt
        a0, b0, c0 = list(aa), list(bb), list(cc)
        for k in range(1, te2 - te + 1):        # aa(te+1:te2) = cc(teh:1:-1), etc.
            aa[te + k] = c0[teh + 1 - k]
            bb[te + k] = b0[teh + 1 - k]
            cc[te + k] = a0[teh + 1 - k]
        return aa, bb, cc, te2

    def lorentz_lu(self, code_dt, vnm1):
        S, Y, X, T, E, n = self.S, self.Y, self.X, self.T, self.E, self.XI
        c1, beta, q = (np.zeros((S, Y, X, T, E, n)) for _ in range(3))
        for is_, ik, it, ig in self.iles():
            for ie in range(1, E + 1):
                aa, bb, cc, te2 = self.lorentz_matrix(ig, ik, it, ie, is_, code_dt, vnm1)
                i = (is_ - 1, ik - 1, it - 1, ig + self.nt, ie - 1)
                for ixi in range(1, n + 1):
                    c1[i + (ixi - 1,)] = cc[ixi]
                q[i + (0,)] = 0.0
                beta[i + (0,)] = 1.0 / bb[1]
                for ixi in range(1, te2):
                    q[i + (ixi,)] = aa[ixi + 1] * beta[i + (ixi - 1,)]
                    beta[i + (ixi,)] = 1.0 / (bb[ixi + 1] - q[i + (ixi,)] * c1[i + (ixi - 1,)])
        return c1, beta, q

    def ediffuse_lu(self, code_dt, vnm2):
        r = self.r
        S, Y, X, T, E, n = self.S, self.Y, self.X, self.T, self.E, self.XI
        ec1, eb, eq = (np.zeros((S, Y, X, T, E, n)) for _ in range(3))
        xe, wM, cap = r["speed_maxw"], r["w_maxw"], r["capg_half"]
        sq = lambda x: math.sqrt(max(0.0, x))  # noqa: E731
        for is_, ik, it, ig in self.iles():
            for ixi in range(1, n + 1):
                il = self.ixi_to_il(ixi, ig)
                if self.forbid(ig, il):
                    continue
                vn = vnm2 * self.sp("spec_vnewk", is_) * float(r["tunits"][ik - 1])
                slb1 = sq(1.0 - self.bmag(ig) * self.al(il))
                x = self.bmag(ig) * self.sp("spec_zstm", is_)
                kp = self.kperp2(ig, it, ik)
                ee = lambda ie: 0.125 * (1 - slb1 * slb1) * self.v("vnew_s", ik, ie, is_) * kp \
                    * r["cfac"] / (x * x)  # noqa: E731
                aa, bb, cc = [0.0] * (E + 1), [0.0] * (E + 1), [0.0] * (E + 1)
                for ie in range(2, E):
                    capgr, capgl = cap[ie - 1], cap[ie - 2]
                    cc[ie] = -0.25 * vn * code_dt * capgr / (wM[ie - 1] * (xe[ie] - xe[ie - 1]))
                    aa[ie] = -0.25 * vn * code_dt * capgl / (wM[ie - 1] * (xe[ie - 1] - xe[ie - 2]))
                    bb[ie] = 1.0 - (aa[ie] + cc[ie]) + ee(ie) * code_dt
                cc[1] = -0.25 * vn * code_dt * cap[0] / (wM[0] * (xe[1] - xe[0]))
                aa[1] = 0.0
                bb[1] = 1.0 - cc[1] + ee(1) * code_dt
                cc[E] = 0.0
                aa[E] = -0.25 * vn * code_dt * cap[E - 2] / (wM[E - 1] * (xe[E - 1] - xe[E - 2]))
                bb[E] = 1.0 - aa[E] + ee(E) * code_dt
                i = (is_ - 1, ik - 1, it - 1, ig + self.nt)
                for ie in range(1, E + 1):
                    ec1[i + (ie - 1, ixi - 1)] = cc[ie]
                eb[i + (0, ixi - 1)] = 1.0 / bb[1]
                for ie in range(1, E):
                    eq[i + (ie, ixi - 1)] = aa[ie + 1] * eb[i + (ie - 1, ixi - 1)]
                    eb[i + (ie, ixi - 1)] = 1.0 / (bb[ie + 1] - eq[i + (ie, ixi - 1)]
                                                    * ec1[i + (ie - 1, ixi - 1)])
        return ec1, eb, eq

    # ----- runtime stages on one ile (gle view (E+1, NX): Fortran gle(ixi, ie) -> view[ie-1, ixi-1])
    @staticmethod
    def _m(arr, i, ie, ixi):
        return 0.0 if ixi > arr.shape[-1] else float(arr[i + (ie - 1, ixi - 1)])

    def ediffuse_ile(self, gle, i, ile, cc, kwork=None):
        is_, ik, it, ig = ile
        if self.skip(is_, ik) or self.filt(it, ik, kwork):
            return
        E = self.E
        max_xi = max(2 * self.jend(ig), 2 * self.ng2)
        for ixi in range(1, max_xi + 1):
            delta = [0j] * (E + 2)
            delta[1] = gle[0, ixi - 1]
            for ie in range(1, E):
                delta[ie + 1] = gle[ie, ixi - 1] - self._m(cc.eqle, i, ie + 1, ixi) * delta[ie]
            gle[E, ixi - 1] = 0.0
            for ie in range(E, 0, -1):
                gle[ie - 1, ixi - 1] = (delta[ie] - self._m(cc.ec1le, i, ie, ixi) * gle[ie, ixi - 1]) \
                    * self._m(cc.ebetaale, i, ie, ixi)

    def lorentz_ile(self, gle, i, ile, cc, kwork=None):
        is_, ik, it, ig = ile
        if self.skip(is_, ik) or self.filt(it, ik, kwork):
            return
        r, ng2, E = self.r, self.ng2, self.E
        cur = self.jend(ig)
        je = max(2 * cur, 2 * ng2 + 1)
        nscatt = je - 1
        is_wfb = cur == ng2 + 1
        special = r["special_wfb_lorentz"]
        if is_wfb and special:
            nscatt -= 1
        for ie in range(1, E + 1):
            col = ie - 1
            if je == 2 * cur and is_wfb and r["vpar_zero_mean"]:
                gle[col, cur - 1] = (gle[col, cur - 1] + gle[col, 2 * cur - 1]) / 2.0
            gle[col, je - 1] = 0.0
            if is_wfb and special:
                gwfb = gle[col, ng2]
                for k in range(ng2 + 1, je - 1):            # gle(ng2+1:je-2) = gle(ng2+2:je-1)
                    gle[col, k - 1] = gle[col, k]
                gle[col, je - 2] = 0.0
            delta = [0j] * (self.NX + 2)
            delta[1] = gle[col, 0]
            for il in range(1, nscatt + 1):
                delta[il + 1] = gle[col, il] - self._m(cc.qle, i, ie, il + 1) * delta[il]
            gle[col, nscatt] = delta[nscatt + 1] * self._m(cc.betaale, i, ie, nscatt + 1)
            for il in range(nscatt, 0, -1):
                gle[col, il - 1] = (delta[il] - self._m(cc.c1le, i, ie, il) * gle[col, il]) \
                    * self._m(cc.betaale, i, ie, il)
            if is_wfb and special:
                for k in range(je - 1, ng2 + 1, -1):        # gle(ng2+2:je-1) = gle(ng2+1:je-2)
                    gle[col, k - 1] = gle[col, k - 2]
                gle[col, ng2] = gwfb
            if cur != 0:
                gle[col, 2 * cur - 1] = gle[col, cur - 1]

    def drag_ile(self, gle, i, ile, cc, apar):
        is_, ik, it, ig = ile
        if not self.electron(is_) or self.filt(it, ik):
            return
        r = self.r
        for ie in range(1, self.E + 1):
            den = (-self.sp("spec_z", is_) * self.sp("spec_dens", is_)) * r["beta"] \
                * self.sp("spec_stm", is_) * float(r["e15"][ie - 1])
            za = cc.vnmult[0] * self.sp("spec_vnewk", is_) * cc.code_dt \
                * self.kperp2(ig, it, ik) * complex(apar[ik - 1, it - 1, ig + self.nt])
            tmp = complex(za.real / den, za.imag / den)      # gfortran: complex/real componentwise
            for ixi in range(1, self.XI + 1):
                gle[ie - 1, ixi - 1] += tmp * float(r["vpa_aj0_le"][i + (ie - 1, ixi - 1)])

    def _sum(self, gle, fac, i, ig):
        v = 0j
        for ie in range(1, self.E + 1):
            for ixi in range(1, self.XI + 1):
                v = v + fac(ie, ixi) * gle[ie - 1, ixi - 1] * self.wxi(ixi, ig)
        return v

    def _sub(self, gle, vec, i, v):
        for ie in range(1, self.E + 2):
            for ixi in range(1, self.XI + 1):
                gle[ie - 1, ixi - 1] = gle[ie - 1, ixi - 1] - float(vec[i + (ie - 1, ixi - 1)]) * v

    def conserve_diffuse_ile(self, gle, i, ile, cc):
        is_, ik, it, ig = ile
        if self.filt(it, ik):
            return
        r, v2 = self.r, cc.vnmult[1]
        W = lambda ie: float(r["w"][is_ - 1, ie - 1])  # noqa: E731
        le = lambda name, ie, ixi: float(r[name][i + (ie - 1, ixi - 1)])  # noqa: E731
        vnE = lambda ie: v2 * self.v("vnew_E", ik, ie, is_) * W(ie)  # noqa: E731
        v = self._sum(gle, lambda ie, ixi: vnE(ie) * le("aj0le", ie, ixi), i, ig)
        self._sub(gle, cc.bz0le, i, v)
        dn = lambda ie: v2 * self.v("delvnew", ik, ie, is_) * W(ie)  # noqa: E731
        v = self._sum(gle, lambda ie, ixi: le("vpa_aj0_le", ie, ixi) * dn(ie), i, ig)
        self._sub(gle, cc.bs0le, i, v)
        v = self._sum(gle, lambda ie, ixi: dn(ie) * le("vperp_aj1le", ie, ixi), i, ig)
        self._sub(gle, cc.bw0le, i, v)

    def conserve_lorentz_ile(self, gle, i, ile, cc):
        is_, ik, it, ig = ile
        if self.filt(it, ik):
            return
        r, v1 = self.r, cc.vnmult[0]
        W = lambda ie: float(r["w"][is_ - 1, ie - 1])  # noqa: E731
        le = lambda name, ie, ixi: float(r[name][i + (ie - 1, ixi - 1)])  # noqa: E731
        if r["drag"]:
            v = 0j
            for ie in range(1, self.E + 1):
                for ixi in range(1, self.XI + 1):
                    v = v + le("vpa_aj0_le", ie, ixi) * gle[ie - 1, ixi - 1] * W(ie) * self.wxi(ixi, ig)
            self._sub(gle, cc.z0le, i, v)
        vpdle = lambda ixi: float(r["vpdiffle"][ig + self.nt, ixi - 1])  # noqa: E731
        nud = lambda ie: float(r["speed_maxw"][ie - 1]) * v1 * self.v("vnew_D", ik, ie, is_) * W(ie)  # noqa: E731,E501
        v = self._sum(gle, lambda ie, ixi: vpdle(ixi) * nud(ie) * le("aj0le", ie, ixi), i, ig)
        self._sub(gle, cc.s0le, i, v)
        nud2 = lambda ie: v1 * self.v("vnew_D", ik, ie, is_) * W(ie)  # noqa: E731
        v = self._sum(gle, lambda ie, ixi: nud2(ie) * le("vperp_aj1le", ie, ixi), i, ig)
        self._sub(gle, cc.w0le, i, v)

    def solfp1(self, gle_all, cc, apar):
        """solfp1_le_layout on a whole LE6 array; returns (gle, stages) with the probe states."""
        gle_all = np.array(gle_all, copy=True)
        r = self.r
        stages = {}

        def each(fn, *args):
            for ile in self.iles():
                is_, ik, it, ig = ile
                i = (is_ - 1, ik - 1, it - 1, ig + self.nt)
                fn(gle_all[i], i, ile, cc, *args)

        has_l = r["collision_model_switch"] in (ci.MODEL_FULL, ci.MODEL_LORENTZ, ci.MODEL_LORENTZ_TEST)
        has_d = r["collision_model_switch"] in (ci.MODEL_FULL, ci.MODEL_EDIFFUSE)
        if has_d:
            each(self.ediffuse_ile)
            stages["coll_after_ediffuse"] = gle_all.copy()
            if r["conserve_moments"]:
                each(self.conserve_diffuse_ile)
            stages["coll_after_ediffuse_conserve"] = gle_all.copy()
        if has_l:
            if r["drag"]:
                each(self.drag_ile, apar)
            stages["coll_after_drag"] = gle_all.copy()
            each(self.lorentz_ile)
            stages["coll_after_lorentz"] = gle_all.copy()
            if r["conserve_moments"]:
                each(self.conserve_lorentz_ile)
            stages["coll_after_lorentz_conserve"] = gle_all.copy()
        return gle_all, stages

    def vspace(self, gnew, phi, apar, bpar, cc):
        """vspace_derivatives, le path (le_derivatives.fpp:67-112)."""
        f = self.r["fphi"], self.r["fbpar"]
        h = self.g_adjust(gnew, phi, bpar, f[0], f[1])
        gle = self.gather(h)
        gle_out, stages = self.solfp1(gle, cc, apar)
        h_out = self.scatter(gle_out, h)
        out = self.g_adjust(h_out, phi, bpar, -f[0], -f[1])
        return out, {"coll_h_in": h, "coll_gle_in": gle, **stages, "coll_gle_out": gle_out,
                     "coll_h_out": h_out}

    # ----- conservation vectors (init_lorentz_conserve / init_diffuse_conserve), real G7 arrays
    def integrate(self, g):
        tot = np.zeros((self.S, self.Y, self.X, self.T))
        for it, ik, il, ie, is_ in self.iglos():
            for ig in range(-self.nt, self.nt + 1):
                t = ig + self.nt
                tot[is_ - 1, ik - 1, it - 1, t] = tot[is_ - 1, ik - 1, it - 1, t] + \
                    self.r["w"][is_ - 1, ie - 1] * self.r["wl"][il - 1, t] * \
                    (g[is_ - 1, ie - 1, il - 1, ik - 1, it - 1, 0, t]
                     + g[is_ - 1, ie - 1, il - 1, ik - 1, it - 1, 1, t])
        return tot

    def _solve_g(self, v, which, cc, kwork):
        le = self.gather(v.astype(complex), nrows=self.XI)
        for ile in self.iles():
            is_, ik, it, ig = ile
            i = (is_ - 1, ik - 1, it - 1, ig + self.nt)
            if which == "lorentz":
                self.lorentz_ile(le[i], i, ile, cc, kwork)
            else:
                self.ediffuse_ile(le[i], i, ile, cc, kwork)
        return self.scatter(le, v.astype(complex)).real

    def _loop(self, fn):
        """for iglo (layout order): for isgn: for ig: fn(idx7, it, ik, il, ie, is, isgn, ig)."""
        for it, ik, il, ie, is_ in self.iglos():
            for isgn in (1, 2):
                for ig in range(-self.nt, self.nt + 1):
                    fn((is_ - 1, ie - 1, il - 1, ik - 1, it - 1, isgn - 1, ig + self.nt),
                       it, ik, il, ie, is_, isgn, ig)

    def lorentz_conserve(self, code_dt, cc, kwork):
        r = self.r
        v1 = cc.vnmult[0]
        shape = r["vpa"].shape
        vpa, aj0, aj1, vpd = r["vpa"], r["aj0"], r["aj1"], r["vpdiff"]
        V, En = r["speed_maxw"], r["energy_maxw"]
        a6 = lambda i: i[:5] + i[6:]  # noqa: E731  G7 index -> G6 index
        vns1 = lambda ik, ie, is_: v1 * self.v("vnew_D", ik, ie, is_)  # noqa: E731
        z0, s0, w0, gtmp = (np.zeros(shape) for _ in range(4))
        syx = lambda D, i: D[i[0], i[3], i[4], i[6]]  # noqa: E731

        def setg(expr):
            def fn(i, it, ik, il, ie, is_, isgn, ig):
                gtmp[i] = expr(i, it, ik, il, ie, is_, isgn, ig)
            self._loop(fn)
            return self.integrate(gtmp)

        if r["drag"]:
            def vns3(ik, ie, is_):
                if not self.electron(is_):
                    return 0.0
                return v1 * self.sp("spec_vnewk", is_) * float(r["tunits"][ik - 1]) / float(r["e15"][ie - 1])

            def fz(i, it, ik, il, ie, is_, isgn, ig):
                z0[i] = -2.0 * code_dt * vns3(ik, ie, is_) * vpd[il - 1, isgn - 1, ig + self.nt] \
                    * V[ie - 1] * aj0[a6(i)]
            self._loop(fz)
            z0 = self._solve_g(z0, "lorentz", cc, kwork)
            D = setg(lambda i, *a: vpa[i] * aj0[a6(i)] * z0[i])
            z0n = z0.copy()
            self._loop(lambda i, *a: z0n.__setitem__(i, z0[i] / (1.0 + syx(D, i))))
            z0 = z0n
        du = setg(lambda i, it, ik, il, ie, is_, isgn, ig: vns1(ik, ie, is_) * vpa[i]
                  * vpd[il - 1, isgn - 1, ig + self.nt] * V[ie - 1])
        duinv = np.where(np.abs(du) > EPS, 1.0 / np.where(du == 0, 1.0, du), du)

        def fs(i, it, ik, il, ie, is_, isgn, ig):
            s0[i] = -vns1(ik, ie, is_) * vpd[il - 1, isgn - 1, ig + self.nt] * V[ie - 1] \
                * aj0[a6(i)] * code_dt * syx(duinv, i)
        self._loop(fs)
        s0 = self._solve_g(s0, "lorentz", cc, kwork)
        D = setg(lambda i, *a: vpa[i] * aj0[a6(i)] * s0[i])
        s0 = self._apply(s0, lambda i: s0[i] - syx(D, i) * z0[i])
        D = setg(lambda i, it, ik, il, ie, is_, isgn, ig: vns1(ik, ie, is_) * V[ie - 1]
                 * vpd[il - 1, isgn - 1, ig + self.nt] * aj0[a6(i)] * s0[i])
        s0 = self._apply(s0, lambda i: s0[i] / (1.0 + syx(D, i)))
        smz = lambda is_: self.sp("spec_smz", is_)  # noqa: E731

        def fw(i, it, ik, il, ie, is_, isgn, ig):
            if not (self.forbid(ig, il) and r["conserve_forbid_zero"]):
                w0[i] = -vns1(ik, ie, is_) * En[ie - 1] * self.al(il) * aj1[a6(i)] * code_dt \
                    * (smz(is_) * smz(is_)) * self.kperp2(ig, it, ik) * syx(duinv, i) / self.bmag(ig)
            else:
                w0[i] = 0.0
        self._loop(fw)
        w0 = self._solve_g(w0, "lorentz", cc, kwork)
        D = setg(lambda i, *a: vpa[i] * aj0[a6(i)] * w0[i])
        w0 = self._apply(w0, lambda i: w0[i] - z0[i] * syx(D, i))
        D = setg(lambda i, it, ik, il, ie, is_, isgn, ig: vns1(ik, ie, is_) * V[ie - 1]
                 * vpd[il - 1, isgn - 1, ig + self.nt] * aj0[a6(i)] * w0[i])
        w0 = self._apply(w0, lambda i: w0[i] - s0[i] * syx(D, i))
        D = setg(lambda i, it, ik, il, ie, is_, isgn, ig: vns1(ik, ie, is_) * En[ie - 1]
                 * self.al(il) * aj1[a6(i)] * w0[i])
        w0 = self._apply(w0, lambda i: w0[i] / (1.0 + syx(D, i)))
        g = lambda v: self.gather(v.astype(complex), nrows=self.XI).real  # noqa: E731
        return g(z0), g(s0), g(w0)

    def diffuse_conserve(self, code_dt, cc, kwork):
        r = self.r
        v2 = cc.vnmult[1]
        shape = r["vpa"].shape
        vpa, aj0, aj1 = r["vpa"], r["aj0"], r["aj1"]
        En = r["energy_maxw"]
        a6 = lambda i: i[:5] + i[6:]  # noqa: E731
        syx = lambda D, i: D[i[0], i[3], i[4], i[6]]  # noqa: E731
        vns1 = lambda ik, ie, is_: v2 * self.v("delvnew", ik, ie, is_)  # noqa: E731
        vnE = lambda ik, ie, is_: self.v("vnew_E", ik, ie, is_)  # noqa: E731
        gtmp, bz0, bs0, bw0 = (np.zeros(shape) for _ in range(4))

        def setg(expr):
            def fn(i, it, ik, il, ie, is_, isgn, ig):
                gtmp[i] = expr(i, it, ik, il, ie, is_, isgn, ig)
            self._loop(fn)
            return self.integrate(gtmp)

        du = setg(lambda i, it, ik, il, ie, is_, isgn, ig: En[ie - 1] * v2 * vnE(ik, ie, is_))
        duinv = np.where(np.abs(du) > EPS, 1 / np.where(du == 0, 1.0, du), du)

        def fz(i, it, ik, il, ie, is_, isgn, ig):
            if not (self.forbid(ig, il) and r["conserve_forbid_zero"]):
                bz0[i] = -code_dt * v2 * vnE(ik, ie, is_) * aj0[a6(i)] * syx(duinv, i)
            else:
                bz0[i] = 0.0
        self._loop(fz)
        bz0 = self._solve_g(bz0, "ediffuse", cc, kwork)
        D = setg(lambda i, it, ik, il, ie, is_, isgn, ig: v2 * vnE(ik, ie, is_) * aj0[a6(i)] * bz0[i])
        bz0 = self._apply(bz0, lambda i: bz0[i] / (1.0 + syx(D, i)))
        du = setg(lambda i, it, ik, il, ie, is_, isgn, ig: vns1(ik, ie, is_) * (vpa[i] * vpa[i]))
        duinv = np.where(np.abs(du) > EPS, 1 / np.where(du == 0, 1.0, du), du)

        def fs(i, it, ik, il, ie, is_, isgn, ig):
            bs0[i] = -vns1(ik, ie, is_) * vpa[i] * aj0[a6(i)] * code_dt * syx(duinv, i)
        self._loop(fs)
        bs0 = self._solve_g(bs0, "ediffuse", cc, kwork)
        D = setg(lambda i, it, ik, il, ie, is_, isgn, ig: v2 * vnE(ik, ie, is_) * aj0[a6(i)] * bs0[i])
        bs0 = self._apply(bs0, lambda i: bs0[i] - syx(D, i) * bz0[i])
        D = setg(lambda i, it, ik, il, ie, is_, isgn, ig: vns1(ik, ie, is_) * vpa[i] * aj0[a6(i)]
                 * bs0[i])
        bs0 = self._apply(bs0, lambda i: bs0[i] / (1.0 + syx(D, i)))
        smz = lambda is_: self.sp("spec_smz", is_)  # noqa: E731

        def fw(i, it, ik, il, ie, is_, isgn, ig):
            if not (self.forbid(ig, il) and r["conserve_forbid_zero"]):
                bw0[i] = -vns1(ik, ie, is_) * En[ie - 1] * self.al(il) * aj1[a6(i)] * code_dt \
                    * (smz(is_) * smz(is_)) * self.kperp2(ig, it, ik) * syx(duinv, i) / self.bmag(ig)
            else:
                bw0[i] = 0.0
        self._loop(fw)
        bw0 = self._solve_g(bw0, "ediffuse", cc, kwork)
        D = setg(lambda i, it, ik, il, ie, is_, isgn, ig: v2 * vnE(ik, ie, is_) * aj0[a6(i)] * bw0[i])
        bw0 = self._apply(bw0, lambda i: bw0[i] - bz0[i] * syx(D, i))
        D = setg(lambda i, it, ik, il, ie, is_, isgn, ig: vns1(ik, ie, is_) * vpa[i] * aj0[a6(i)]
                 * bw0[i])
        bw0 = self._apply(bw0, lambda i: bw0[i] - bs0[i] * syx(D, i))
        D = setg(lambda i, it, ik, il, ie, is_, isgn, ig: vns1(ik, ie, is_) * En[ie - 1]
                 * self.al(il) * aj1[a6(i)] * bw0[i])
        bw0 = self._apply(bw0, lambda i: bw0[i] / (1.0 + syx(D, i)))
        g = lambda v: self.gather(v.astype(complex), nrows=self.XI).real  # noqa: E731
        return g(bz0), g(bs0), g(bw0)

    def _apply(self, v, fn):
        out = np.array(v, copy=True)
        for i in np.ndindex(v.shape):
            out[i] = fn(i)
        return out


# ======================================================================== helpers
def _build(raw, **param_overrides):
    params = ci.CollParams.from_raw(raw, **param_overrides)
    data = ci.CollData.from_raw(raw, params)
    return params, data


def _rand_state(raw, seed=7):
    rng = np.random.default_rng(seed)
    S, E, L, Y, X, T = (raw[k] for k in ("nspec", "negrid", "nlambda", "naky", "ntheta0", "nth"))
    g = rng.normal(size=(S, E, L, Y, X, 2, T)) + 1j * rng.normal(size=(S, E, L, Y, X, 2, T))
    f = [rng.normal(size=(Y, X, T)) + 1j * rng.normal(size=(Y, X, T)) for _ in range(3)]
    return g, Fields(*f)


def _err(ref, got):
    return max_abs_rel(np.asarray(ref), np.asarray(got))


SYNTH_VARIANTS = {
    "full": dict(),
    "full_yxels_hyper": dict(layout="yxels", hyper_colls=True, hypermult=True, seed=3),
    "full_special_wfb": dict(special_wfb=True, seed=4),
    "full_pwzero_novzm": dict(pw_wfb_zero=True, vpar_zero_mean=False, seed=5),
    "lorentz_nodrag": dict(model=ci.MODEL_LORENTZ, drag=False, seed=6),
    "lorentz_test": dict(model=ci.MODEL_LORENTZ_TEST, seed=8),
    "ediffuse_skip": dict(model=ci.MODEL_EDIFFUSE, zero_vnew_species=0, seed=9),
    "full_noconserve": dict(conserve_moments=False, seed=10),
    "notrapped_ediffuse": dict(model=ci.MODEL_EDIFFUSE, trapped=False, seed=11),
    "notrapped_lorentz_nocons": dict(model=ci.MODEL_LORENTZ, trapped=False,
                                     conserve_moments=False, seed=12),
}


@pytest.fixture(scope="module", params=sorted(SYNTH_VARIANTS))
def synth(request):
    raw = synth_raw(**SYNTH_VARIANTS[request.param])
    params, data = _build(raw)
    lit = Lit(raw)
    dt = 0.037
    coeffs = ci.compute_coll_coeffs(dt, raw["vnmult"], data, params, kwork_filter=False)
    return request.param, raw, params, data, lit, dt, coeffs


# ======================================================================== synthetic tests
def test_synthetic_grid_is_consistent():
    """The synthetic grid satisfies the le_grids invariants used by the port."""
    raw = synth_raw()
    lit = Lit(raw)
    T = raw["nth"]
    for t in range(T):
        ig = t - raw["ntgrid"]
        je = lit.jend(ig)
        assert je >= raw["ng2"] + 1
        assert all(not lit.forbid(ig, il) for il in range(1, je + 1))
        assert all(lit.forbid(ig, il) for il in range(je + 1, raw["nlambda"] + 1))
        for ixi in range(1, 2 * raw["nlambda"] + 1):   # g2le formula == ixi_to_* maps
            il, sg = lit.ixi_to_il(ixi, ig), lit.ixi_to_isgn(ixi, ig)
            if not lit.forbid(ig, il):
                assert lit.ixi_of(ig, il, sg) == ixi


def test_gather_scatter(synth):
    name, raw, params, data, lit, dt, coeffs = synth
    g, _ = _rand_state(raw)
    gle = coll.g2le(jnp.asarray(g), data, params)
    ref = lit.gather(g)
    assert gle.shape == ref.shape
    assert _err(ref, gle) == 0.0
    # host version used by the rebuild
    maps = ci.le_maps(raw["ixi_to_il"], raw["ixi_to_isgn"], raw["forbid"], params.nx)
    assert _err(ref, ci.gather_le(g, maps, params.nx)) == 0.0
    # scatter of a different array onto g: allowed points replaced, forbidden kept
    other = ref * (1.5 - 0.5j)
    got = coll.le2g(jnp.asarray(other), jnp.asarray(g), data, params)
    assert _err(lit.scatter(other, g), got) == 0.0
    assert _err(lit.scatter(other, g), ci.scatter_le(other, g, maps)) == 0.0
    # round trip
    assert _err(g, coll.le2g(gle, jnp.asarray(g), data, params)) == 0.0


def test_static_le_arrays(synth):
    """init_le_bessel / init_vpdiff transcriptions reproduce the (literal) exported arrays."""
    name, raw, params, data, lit, dt, coeffs = synth
    maps = ci.le_maps(raw["ixi_to_il"], raw["ixi_to_isgn"], raw["forbid"], params.nx)
    got = ci.init_le_bessel(raw["aj0"], raw["aj1"], raw["vpa"], raw["energy_maxw"], raw["al"],
                            maps, params)
    for ref, val in zip((raw["aj0le"], raw["vperp_aj1le"], raw["vpa_aj0_le"]), got):
        assert ref.shape == val.shape
        assert _err(ref, val) == 0.0
    if raw["vpdiff"] is not None:
        vp = ci.init_vpdiff(raw["bmag"], raw["al"], raw["jend"], raw["pitch_weights"], params)
        assert _err(raw["vpdiff"], vp) == 0.0
        vle = ci.vpdiff_le(vp, raw["ixi_to_il"], raw["ixi_to_isgn"], params.nxi_lim)
        assert _err(raw["vpdiffle"], vle) == 0.0


def test_matrices_literal(synth):
    """Lorentz / energy-diffusion LU arrays == literal get_*_matrix + elimination (bitwise)."""
    name, raw, params, data, lit, dt, coeffs = synth
    v = raw["vnmult"]
    if params.has_lorentz:
        for ref, got in zip(lit.lorentz_lu(dt, v[0]), (coeffs.c1le, coeffs.betaale, coeffs.qle)):
            assert got.shape == ref.shape
            assert _err(ref, got) <= 1e-15, name
    else:
        assert coeffs.c1le is None
    if params.has_diffuse:
        for ref, got in zip(lit.ediffuse_lu(dt, v[1]), (coeffs.ec1le, coeffs.ebetaale, coeffs.eqle)):
            assert _err(ref, got) <= 1e-15, name
    else:
        assert coeffs.ec1le is None


def test_lorentz_matrix_solves_tridiagonal_system():
    """The elimination + GS2 sweep solve the tridiagonal (I + C dt) system (orientation check)."""
    raw = synth_raw(model=ci.MODEL_LORENTZ, drag=False, conserve_moments=False, kwork00=False,
                    vpar_zero_mean=False)
    params, data = _build(raw)
    lit = Lit(raw)
    dt = 0.05
    coeffs = ci.compute_coll_coeffs(dt, raw["vnmult"], data, params)
    rng = np.random.default_rng(0)
    for ig in (-raw["ntgrid"], 0, 1):              # bmax (wfb bounce), bmin, intermediate
        aa, bb, cc, te2 = lit.lorentz_matrix(ig, 1, 1, 2, 1, dt, raw["vnmult"][0])
        A = np.zeros((te2, te2))
        for k in range(te2):
            A[k, k] = bb[k + 1]
            if k > 0:
                A[k, k - 1] = aa[k + 1]
            if k < te2 - 1:
                A[k, k + 1] = cc[k + 1]
        rhs = rng.normal(size=te2)
        le = np.zeros((params.nspec, params.naky, params.ntheta0, params.nth, params.negrid + 1,
                       params.nx), complex)
        t = ig + raw["ntgrid"]
        le[0, 0, 0, t, 1, :te2] = rhs
        cj = lit.jend(ig)
        out = np.asarray(coll.solfp_lorentz(jnp.asarray(le), coeffs, data, params))
        x = out[0, 0, 0, t, 1, :te2].real
        assert np.max(np.abs(A @ x - rhs)) <= 1e-12 * np.max(np.abs(rhs))
        assert out[0, 0, 0, t, 1, 2 * cj - 1] == out[0, 0, 0, t, 1, cj - 1]   # bounce copy


@pytest.mark.parametrize("kwork_build", [False, None])
def test_conservation_vectors_literal(synth, kwork_build):
    """init_lorentz_conserve / init_diffuse_conserve == literal transcriptions (iglo summation order),
    for the first-init filter (all False) and the dt-reset filter (box (1,1) skipped)."""
    name, raw, params, data, lit, dt, coeffs = synth
    if not params.conserve_moments:
        assert coeffs.z0le is None and coeffs.bz0le is None
        return
    cc = ci.compute_coll_coeffs(dt, raw["vnmult"], data, params, kwork_filter=kwork_build)
    kwork = np.zeros_like(raw["kwork_filter"]) if kwork_build is False else raw["kwork_filter"]
    if params.has_lorentz:
        for ref, got, nm in zip(lit.lorentz_conserve(dt, cc, kwork), (cc.z0le, cc.s0le, cc.w0le),
                                ("z0le", "s0le", "w0le")):
            assert ref.shape == got.shape, nm
            assert _err(ref, got) <= 1e-15, f"{name} {nm} {_err(ref, got):.3e}"
    if params.has_diffuse:
        for ref, got, nm in zip(lit.diffuse_conserve(dt, cc, kwork), (cc.bz0le, cc.bs0le, cc.bw0le),
                                ("bz0le", "bs0le", "bw0le")):
            assert _err(ref, got) <= 1e-15, f"{name} {nm} {_err(ref, got):.3e}"


def test_kwork_filter_build_state_matters_only_for_filtered_mode():
    raw = synth_raw(kwork00=True)
    params, data = _build(raw)
    a = ci.compute_coll_coeffs(0.04, raw["vnmult"], data, params, kwork_filter=False)
    b = ci.compute_coll_coeffs(0.04, raw["vnmult"], data, params)        # exported filter
    for nm in ci.LORENTZ_VECTORS + ci.EDIFFUSE_VECTORS:
        x, y = getattr(a, nm), getattr(b, nm)
        assert np.array_equal(x[:, 1:], y[:, 1:]) and np.array_equal(x[:, :, 1:], y[:, :, 1:])
    assert not np.array_equal(a.s0le[:, 0, 0], b.s0le[:, 0, 0])
    for nm in ci.LORENTZ_MATRICES + ci.EDIFFUSE_MATRICES:
        assert np.array_equal(getattr(a, nm), getattr(b, nm))


def test_stages_and_vspace_literal(synth):
    """Every solfp1 stage and the whole vspace_derivatives (jit and eager) == literal Fortran loops."""
    name, raw, params, data, lit, dt, coeffs = synth
    g, fields = _rand_state(raw)
    ref, ref_stages = lit.vspace(g, fields.phi, fields.apar, fields.bpar, coeffs)
    jd, jc = data.device(), to_jax(coeffs)
    out, stages = coll.vspace_derivatives(jnp.asarray(g), None, None, to_jax(fields),
                                          jnp.asarray(fields.apar), None, jc, jd, params,
                                          return_stages=True)
    assert set(stages) == set(ref_stages)
    for tag, val in ref_stages.items():
        e = _err(val, stages[tag])
        assert e <= 1e-14, f"{name} {tag}: {e:.3e}"
    assert _err(ref, out) <= 1e-14
    out_jit = coll.vspace_derivatives_jit(jnp.asarray(g), None, None, to_jax(fields),
                                          jnp.asarray(fields.apar), 1, jc, data, params)
    assert _err(ref, out_jit) <= 1e-14
    # forbidden points of h never change through the operator (g_adjust skips them too)
    forbid = raw["forbid"][None, None, :, None, None, None, :]
    assert np.array_equal(np.where(forbid, np.asarray(out), 0), np.where(forbid, g, 0))


def test_stage_kernels_individually(synth):
    """Each stage kernel on the literal previous-stage input (as the trace test does)."""
    name, raw, params, data, lit, dt, coeffs = synth
    g, fields = _rand_state(raw, seed=11)
    _, st = lit.vspace(g, fields.phi, fields.apar, fields.bpar, coeffs)
    order = [t for t in coll.STAGE_TAGS if t in st]
    kern = {
        "coll_after_ediffuse": lambda x: coll.solfp_ediffuse(x, coeffs, data, params),
        "coll_after_ediffuse_conserve": (lambda x: coll.conserve_diffuse(x, coeffs, data, params))
        if params.conserve_moments else (lambda x: x),
        "coll_after_drag": (lambda x: coll.add_drag(x, jnp.asarray(fields.apar), coeffs, data, params))
        if params.drag else (lambda x: x),
        "coll_after_lorentz": lambda x: coll.solfp_lorentz(x, coeffs, data, params),
        "coll_after_lorentz_conserve": (lambda x: coll.conserve_lorentz(x, coeffs, data, params))
        if params.conserve_moments else (lambda x: x),
    }
    for prev, tag in zip(order[1:-3], order[2:-2]):
        got = kern[tag](jnp.asarray(st[prev]))
        assert _err(st[tag], got) <= 1e-14, f"{name} {tag}"


def test_exact_vs_tree_sums():
    raw = synth_raw()
    params, data = _build(raw)
    fast = ci.CollParams.from_raw(raw, exact_sums=False)
    coeffs = ci.compute_coll_coeffs(0.05, raw["vnmult"], data, params)
    g, f = _rand_state(raw)
    a = coll.vspace_derivatives(jnp.asarray(g), None, None, f, f.apar, None, coeffs, data, params)
    b = coll.vspace_derivatives(jnp.asarray(g), None, None, f, f.apar, None, coeffs, data, fast)
    assert _err(a, b) <= 1e-13


def test_collisions_advance_synthetic():
    raw = synth_raw()
    raw["split_collisions"] = True
    params, data = _build(raw)
    lit = Lit(raw)
    coeffs = ci.compute_coll_coeffs(0.05, raw["vnmult"], data, params)
    g, f = _rand_state(raw)
    scale = np.linspace(0.5, 1.5, raw["nth"])

    def potentials(h):          # a stand-in for calculate_potentials_from_nonadiabatic_dfn
        s = jnp.sum(h, axis=(0, 1, 2, 5)) * scale
        return Fields(s, 2.0 * s, -0.5 * s)

    out, newf = coll.collisions_advance(jnp.asarray(g), None, None, f, potentials, coeffs, data,
                                        params)
    ref, _ = lit.vspace(g, f.phi, f.apar, f.bpar, coeffs)
    h = lit.g_adjust(ref, f.phi, f.bpar, 1.0, 1.0)
    nf = potentials(jnp.asarray(h))
    ref = lit.g_adjust(h, np.asarray(nf.phi), np.asarray(nf.bpar), -1.0, -1.0)
    assert _err(ref, out) <= 1e-13
    assert _err(nf.phi, newf.phi) <= 1e-13
    off = ci.CollParams.from_raw(raw, split_collisions=False)
    same, f2 = coll.collisions_advance(jnp.asarray(g), None, None, f, potentials, coeffs, data, off)
    assert np.array_equal(np.asarray(same), g) and f2 is f


def test_no_collisions_is_identity():
    raw = synth_raw()
    raw["colls"] = False
    params, data = _build(raw)
    coeffs = ci.compute_coll_coeffs(0.05, raw["vnmult"], data, params)
    assert coeffs.c1le is None
    g, f = _rand_state(raw)
    out = coll.vspace_derivatives(jnp.asarray(g), None, None, f, f.apar, None, coeffs, data, params)
    assert np.array_equal(np.asarray(out), g)


@pytest.mark.parametrize("bad", [dict(use_le_layout=False), dict(conservative=False),
                                 dict(lorentz_switch=2), dict(ediff_switch=2)])
def test_gates(bad):
    raw = synth_raw()
    raw.update(bad)
    with pytest.raises(NotImplementedError):
        ci.CollParams.from_raw(raw)


def test_gate_hybrid_and_notrapped():
    raw = synth_raw()
    raw["spec_type"] = np.array([1, 4], np.int32)
    with pytest.raises(NotImplementedError):
        ci.CollParams.from_raw(raw)
    raw = synth_raw(trapped=False, model=ci.MODEL_LORENTZ, conserve_moments=True)
    with pytest.raises(NotImplementedError):
        ci.CollParams.from_raw(raw)


def test_strict_mode_bitwise(synth):
    """With FMA-free products (numerics strict mode) the jitted kernel equals the literal loops
    bitwise, stage by stage (the reference build uses -ffp-contract=off, 00_decisions D2)."""
    from gs2gpu import numerics
    name, raw, params, data, lit, dt, coeffs = synth
    g, f = _rand_state(raw, seed=3)
    ref, ref_stages = lit.vspace(g, f.phi, f.apar, f.bpar, coeffs)
    numerics.set_strict(True)
    try:
        fn = jax.jit(coll.vspace_derivatives, static_argnums=(5, 8),
                     static_argnames=("return_stages",))
        out, stages = fn(jnp.asarray(g), None, None, to_jax(f), jnp.asarray(f.apar), None,
                         to_jax(coeffs), data.device(), params, return_stages=True)
    finally:
        numerics.set_strict(False)
    for tag, val in ref_stages.items():
        assert np.array_equal(np.asarray(stages[tag]), val), f"{name} {tag}: {_err(val, stages[tag]):.2e}"
    assert np.array_equal(np.asarray(out), ref)


# ======================================================================== reference-data tests
def _coll_cases():
    return cases_where(lambda ex: "collisions" in ex.root.groups and
                       bool(ex["/collisions"].get("colls", False)))


def _ladder_rungs(ex):
    if "dt_ladder" not in ex.root.groups:
        return []
    return sorted(int(k[1:]) for k in ex["/dt_ladder"].groups if k.startswith("k"))


def _setup(case):
    ex = export_of(case)
    params = ci.CollParams.from_export(ex)
    data = ci.CollData.from_export(ex, params)
    return ex, params, data


@pytest.mark.parametrize("case", _coll_cases())
def test_ref_coeffs_initial_dt(case):
    ex, params, data = _setup(case)
    ref = ci.CollCoeffs.from_export(ex)
    got = ci.compute_coll_coeffs(ref.code_dt, ref.vnmult, data, params, kwork_filter=False)
    n = 0
    for nm in ci.COEFF_NAMES:
        r, g = getattr(ref, nm), getattr(got, nm)
        assert (r is None) == (g is None), f"{case} {nm}: export {r is None} vs rebuild {g is None}"
        if r is not None:
            assert_close(r, g, RTOL_KERNEL, f"{case} {nm}")
            n += 1
    assert n > 0


@pytest.mark.parametrize("case", _coll_cases())
def test_ref_coeffs_dt_ladder(case):
    ex, params, data = _setup(case)
    rungs = _ladder_rungs(ex)
    if not rungs:
        pytest.skip(f"{case}: no dt ladder in the export")
    for k in rungs:
        ref = ci.CollCoeffs.from_export(ex, rung=k)
        got = ci.compute_coll_coeffs(ref.code_dt, ref.vnmult, data, params)
        for nm in ci.COEFF_NAMES:
            r = getattr(ref, nm)
            if r is not None:
                assert_close(r, getattr(got, nm), RTOL_KERNEL, f"{case} rung {k} {nm}")


@pytest.mark.parametrize("case", _coll_cases())
def test_ref_static_le_arrays(case):
    """Our g2le gather of aj0/aj1/vpa reproduces the exported le arrays (validates the index maps)."""
    ex, params, data = _setup(case)
    raw = ci.raw_from_export(ex)
    maps = ci.le_maps(raw["ixi_to_il"], raw["ixi_to_isgn"], raw["forbid"], params.nx)
    got = ci.init_le_bessel(raw["aj0"], raw["aj1"], raw["vpa"], np.asarray(raw["energy_maxw"]),
                            np.asarray(raw["al"]), maps, params)
    for nm, val in zip(("aj0le", "vperp_aj1le", "vpa_aj0_le"), got):
        if raw[nm] is not None:
            assert_close(raw[nm], val, 0.0 if nm == "aj0le" else RTOL_KERNEL, f"{case} {nm}")
    if raw["vpdiff"] is not None and raw["pitch_weights"] is not None:
        vp = ci.init_vpdiff(np.asarray(raw["bmag"]), np.asarray(raw["al"]), np.asarray(raw["jend"]),
                            np.asarray(raw["pitch_weights"]), params)
        assert_close(raw["vpdiff"], vp, RTOL_KERNEL, f"{case} vpdiff")


def _coeffs_provider(ex, data, params):
    """Coefficients valid at a trace record: the exported set at the initial dt, else (after a dt
    reset inside the traced steps) the rebuild at the record's ``code_dt`` (reset filter semantics)."""
    init = ci.CollCoeffs.from_export(ex)
    cache = {init.code_dt: to_jax(init)}

    def get(grp):
        dt = float(grp.get("code_dt", init.code_dt))
        if dt not in cache:
            vnm = np.asarray(grp.get("vnmult", init.vnmult), np.float64)
            cache[dt] = to_jax(ci.compute_coll_coeffs(dt, vnm, data, params))
        return cache[dt]

    return get


def _probe_array(grp, shape=None, prefer=("gle", "h", "arr", "gnew")):
    """The array a probe record dumped: ``gle`` (le probes) or ``h`` (``coll_h_*``), as written by
    src/gs2_export_trace.fpp (tr_le, tr_g); fall back to the unique array with the expected shape."""
    for nm in prefer:
        if nm in grp.arrays:
            a = grp.arrays[nm]
            if shape is None or np.shape(a) == tuple(shape):
                return a
    cands = [a for a in grp.arrays.values() if np.ndim(a) >= 6
             and (shape is None or np.shape(a) == tuple(shape))]
    if len(cands) == 1:
        return cands[0]
    raise KeyError(f"cannot identify the probe array among {sorted(grp.arrays)}")


def _fields_for(tr, rec, dims_shape):
    """phinew/aparnew/bparnew valid during ``rec``'s timeadv pass (pass 1: start of advance,
    after exb; pass 2: after the field update / zonal removal)."""
    tags = ("after_exb_shear", "advance_begin") if rec.pass_ == 1 else \
        ("after_remove_zonal_flows", "after_field_update")
    if rec.tag == "after_timeadv":          # records the fields valid at its pass itself
        tags = (rec.tag,)
    for tag in tags:
        prev = rec if tag == rec.tag else tr.previous(rec, tag)
        if prev is not None:
            grp = tr.load(prev)
            out = []
            for nm in ("phinew", "aparnew", "bparnew"):
                v = grp.get(nm)
                out.append(np.zeros(dims_shape, complex) if v is None else v)
            return Fields(*out)
    raise KeyError(f"no field record before {rec}")


def _vspace_records(tr):
    return [r for r in tr.records if r.tag == "after_vspace_derivatives" and r.in_advance]


@pytest.mark.parametrize("case", _coll_cases())
def test_ref_collision_stages(case):
    """Each collision probe of the trace from the preceding probe (spec 09 section 5.2)."""
    ex, params, data = _setup(case)
    if not have(case, "trace"):
        pytest.skip(f"{case}: no trace")
    tr = trace_of(case)
    coeffs_at = _coeffs_provider(ex, data, params)
    jd = data.device()
    recs = [r for r in tr.records if r.tag == "coll_h_in" and r.in_advance]
    if not recs:
        pytest.skip(f"{case}: no collision probes in the trace")
    le_shape = (params.nspec, params.naky, params.ntheta0, params.nth, params.negrid + 1, params.nx)
    for rec in recs:
        f = _fields_for(tr, rec, ex.field_shape)
        coeffs = coeffs_at(tr.load(rec))
        h = _probe_array(tr.load(rec), ex.g_shape)
        prev = h
        chain = {}
        r = rec
        while True:
            r = tr.next(r)
            if r is None or not r.tag.startswith("coll_"):
                break
            chain[r.tag] = r
        # gather
        if "coll_gle_in" in chain:
            gle_in = _probe_array(tr.load(chain["coll_gle_in"]), le_shape)
            assert_close(gle_in, coll.g2le(jnp.asarray(h), jd, params), RTOL_KERNEL,
                         f"{case} {rec.istep}/{rec.pass_} coll_gle_in")
            prev = gle_in
        kern = {
            "coll_after_ediffuse": lambda x: coll.solfp_ediffuse(x, coeffs, jd, params),
            "coll_after_ediffuse_conserve": lambda x: coll.conserve_diffuse(x, coeffs, jd, params)
            if params.conserve_moments else x,
            "coll_after_drag": lambda x: coll.add_drag(x, jnp.asarray(f.apar), coeffs, jd, params)
            if params.drag else x,
            "coll_after_lorentz": lambda x: coll.solfp_lorentz(x, coeffs, jd, params),
            "coll_after_lorentz_conserve": lambda x: coll.conserve_lorentz(x, coeffs, jd, params)
            if params.conserve_moments else x,
        }
        for tag in ("coll_after_ediffuse", "coll_after_ediffuse_conserve", "coll_after_drag",
                    "coll_after_lorentz", "coll_after_lorentz_conserve"):
            if tag not in chain:
                continue
            ref = _probe_array(tr.load(chain[tag]), le_shape)
            got = kern[tag](jnp.asarray(prev))
            assert_close(ref, got, RTOL_KERNEL, f"{case} {rec.istep}/{rec.pass_} {tag}")
            prev = ref
        if "coll_h_out" in chain and "coll_gle_out" in chain:
            gle_out = _probe_array(tr.load(chain["coll_gle_out"]), le_shape)
            ref = _probe_array(tr.load(chain["coll_h_out"]), ex.g_shape)
            assert_close(ref, coll.le2g(jnp.asarray(gle_out), jnp.asarray(h), jd, params),
                         RTOL_KERNEL, f"{case} {rec.istep}/{rec.pass_} coll_h_out")


@pytest.mark.parametrize("case", _coll_cases())
def test_ref_vspace_derivatives(case):
    """The whole vspace_derivatives: input gnew of after_hyper (or after_invert_rhs), output
    after_vspace_derivatives, for every traced step and both passes."""
    ex, params, data = _setup(case)
    if params.split_collisions:
        pytest.skip(f"{case}: split_collisions (collisions run in collisions_advance)")
    if not have(case, "trace"):
        pytest.skip(f"{case}: no trace")
    tr = trace_of(case)
    coeffs_at = _coeffs_provider(ex, data, params)
    recs = _vspace_records(tr)
    if not recs:
        pytest.skip(f"{case}: no after_vspace_derivatives records")
    for rec in recs:
        src = tr.previous(rec, "after_hyper") or tr.previous(rec, "after_invert_rhs")
        assert src is not None
        gin = tr.load(src)["gnew"]
        f = _fields_for(tr, rec, ex.field_shape)
        got = coll.vspace_derivatives_jit(jnp.asarray(gin), None, None, to_jax(f),
                                          jnp.asarray(f.apar), rec.pass_ - 1 or None,
                                          coeffs_at(tr.load(rec)), data, params)
        assert_close(tr.load(rec)["gnew"], got, RTOL_KERNEL,
                     f"{case} istep={rec.istep} attempt={rec.attempt} pass={rec.pass_}")


@pytest.mark.parametrize("case", _coll_cases())
def test_ref_vspace_strict_bitwise(case):
    """With FMA-free products (numerics strict mode, XLA:CPU) every collision stage and the whole
    vspace_derivatives are bitwise identical to the -ffp-contract=off reference build."""
    from gs2gpu import numerics
    if jax.default_backend() != "cpu":
        pytest.skip("strict bitwise check is for XLA:CPU")
    ex, params, data = _setup(case)
    if params.split_collisions or not have(case, "trace"):
        pytest.skip(f"{case}: split collisions or no trace")
    tr = trace_of(case)
    coeffs_at = _coeffs_provider(ex, data, params)
    numerics.set_strict(True)
    try:
        fn = jax.jit(coll.vspace_derivatives, static_argnums=(5, 8),
                     static_argnames=("return_stages",))
        for rec in _vspace_records(tr):
            src = tr.previous(rec, "after_hyper") or tr.previous(rec, "after_invert_rhs")
            f = _fields_for(tr, rec, ex.field_shape)
            out, st = fn(jnp.asarray(tr.load(src)["gnew"]), None, None, to_jax(f),
                         jnp.asarray(f.apar), None, coeffs_at(tr.load(rec)), data.device(), params,
                         return_stages=True)
            what = f"{case} istep={rec.istep} pass={rec.pass_}"
            assert np.array_equal(np.asarray(out), tr.load(rec)["gnew"]), what
            for tag in coll.STAGE_TAGS:
                prev = tr.previous(rec, tag)
                if prev is not None and tag in st:
                    assert np.array_equal(np.asarray(st[tag]), _probe_array(tr.load(prev))), \
                        f"{what} {tag}"
    finally:
        numerics.set_strict(False)


def _potentials_fn(ex):
    """``calculate_potentials_from_nonadiabatic_dfn`` from the fields module (B4) when it exists,
    else a test-local transcription of dist_fn.fpp:7322-7455 (spec 04 sections 3.1, 3.2, 3.9)."""
    try:
        from gs2gpu import fields as fl
        fdata, fparams = to_jax(fl.FieldData.from_export(ex)), fl.FieldParams.from_export(ex)
        return lambda h: fl.calculate_potentials_from_nonadiabatic_dfn(h, fdata, fparams)
    except Exception:  # noqa: BLE001 - B4 not available yet
        pass
    from gs2gpu.numerics import rmul
    sp, le, rp, dfa = ex["/species"], ex["/le_grids"], ex["/run_parameters"], ex["/dist_fn_arrays"]
    kt, th = ex["/kt_grids"], ex["/theta_grid"]
    z, dens, temp, mass = (np.asarray(sp[f"spec_{n}"]) for n in ("z", "dens", "temp", "mass"))
    types = np.asarray(sp["spec_type"])
    w, wl = np.asarray(le["w"]), np.asarray(le["wl"])
    aj0, aj1, vperp2, vpa = (np.asarray(dfa[n]) for n in ("aj0", "aj1", "vperp2", "vpa"))
    kwork = np.asarray(kt["kwork_filter"], bool)
    layout = ci._layout_of(ex)
    ie_outer = layout.index("l") < layout.index("e")
    has_phi, has_apar, has_bpar = (bool(rp[f"has_{n}"]) for n in ("phi", "apar", "bpar"))
    beta, tite = float(rp["beta"]), float(rp.get("tite", 1.0))
    S, E, L = ex.nspec, ex.negrid, ex.nlambda

    def integrate(gw, wgt):
        wgt = np.where(types == 3, 0.0, wgt)
        tot = jnp.zeros(ex.field_shape, complex)
        order = [(s, e, l) for s in range(S) for e in range(E) for l in range(L)] if ie_outer else \
            [(s, e, l) for s in range(S) for l in range(L) for e in range(E)]
        for s, e, l in order:
            fac2 = (wgt[s] * w[s, e]) * wl[l]
            tot = tot + rmul(gw[s, e, l, :, :, 0, :] + gw[s, e, l, :, :, 1, :], fac2)
        return jnp.where(kwork[:, :, None], 0.0 + 0.0j, tot)

    ele, ion = np.isin(types, (2, 4)).any(), (types == 1).any()
    acc = 0.0
    for s in range(S):
        acc = acc + ((dens[s] * z[s]) * z[s]) / temp[s]
    phi_weight = 1.0 / (acc + tite) if (not ele or not ion) else 1.0 / acc

    def pot(h):
        zero = jnp.zeros(ex.field_shape, complex)
        phi = apar = bpar = zero
        if has_phi:
            phi = rmul(integrate(rmul(h, aj0[..., None, :]), z * dens), phi_weight)
        if has_apar:
            gw = rmul(h, aj0[..., None, :] * vpa)
            wg = (((2.0 * beta) * z) * dens) * np.sqrt(temp / mass)
            apar = rmul(integrate(gw, wg), np.asarray(kt["inv_kperp2"]))
        if has_bpar:
            gw = rmul(h, (aj1 * vperp2)[..., None, :])
            bm = np.asarray(th["bmag"])
            bpar = rmul(integrate(gw, temp * dens), -beta / (bm * bm))
        return Fields(phi, apar, bpar)

    return pot


@pytest.mark.parametrize("case", [c for c in _coll_cases()
                                  if c in ("lin_single_coll_ke_split", "<no-refdata>")])
def test_ref_collisions_advance(case):
    ex, params, data = _setup(case)
    if not params.split_collisions:
        pytest.skip(f"{case}: not split_collisions")
    if not have(case, "trace"):
        pytest.skip(f"{case}: no trace")
    tr = trace_of(case)
    coeffs_at = _coeffs_provider(ex, data, params)
    pot = _potentials_fn(ex)
    recs = [r for r in tr.records if r.tag == "after_collisions_advance"]
    if not recs:
        pytest.skip(f"{case}: no after_collisions_advance records")
    for rec in recs:
        tin = tr.previous(rec, "after_timeadv")
        assert tin is not None and tin.pass_ == 2
        f = _fields_for(tr, tin, ex.field_shape)
        gin = tr.load(tin)["gnew"]
        got, newf = coll.collisions_advance(jnp.asarray(gin), None, None, to_jax(f), pot,
                                            coeffs_at(tr.load(rec)), data.device(), params)
        assert_close(tr.load(rec)["gnew"], got, RTOL_KERNEL,
                     f"{case} istep={rec.istep} attempt={rec.attempt}")
        end = tr.next(rec, "advance_end")
        if end is not None:
            grp = tr.load(end)
            for nm, val in zip(("phinew", "aparnew", "bparnew"), newf):
                if grp.get(nm) is not None:
                    assert_close(grp[nm], val, RTOL_KERNEL, f"{case} istep={rec.istep} {nm}")


def test_reference_cases_listed():
    """The five collision cases of gpu/cases exist (their refdata may not yet)."""
    from conftest import all_cases
    for c in COLL_CASES_ALL:
        assert c in all_cases()


def test_ky_chunked_rebuild_identical():
    raw = synth_raw(naky=3, seed=31)
    params, data = _build(raw)
    a = ci.compute_coll_coeffs(0.05, raw["vnmult"], data, params)
    b = ci.compute_coll_coeffs(0.05, raw["vnmult"], data, params, max_g_elements=1)
    for nm in ci.COEFF_NAMES:
        assert np.array_equal(getattr(a, nm), getattr(b, nm)), nm
