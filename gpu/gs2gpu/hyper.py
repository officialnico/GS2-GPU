"""Hyperviscosity and hyperresistivity (src/hyper.f90), spec 07 section 3.2.

GS2 8.2.1 has exactly the options ``default|none|visc_only|res_only|both|simple3D``
(hyper.f90:458-464). This module ports:

* ``read_parameters`` switch logic (hyper.f90:452-601) as ``resolve_hyper_parameters``;
* the normalisers of ``init_hyper`` (hyper.f90:386-396) and the frozen arrays
  ``D_res`` (400-411) and the simple3D ``hypervisc_filter`` (419-445) as NumPy
  formulas (``normalisers``, ``d_res_formula``, ``simple3d_filter_formula``), used
  for validation and as a fallback; the time step uses the **exported** arrays;
* ``hyper_diff`` (hyper.f90:620-786) with ``iso_shear`` (753-784), ``aniso_shear``
  (688-750) and ``simple3Dfilter`` (673-685) as a jit-compatible kernel.

GS2 quirks replicated (00_decisions D4):

* ``D_res`` and the simple3D filter contain ``code_dt`` at the first ``init_hyper``
  call and are never rebuilt (level 19 is not torn down by ``reset_time_step``;
  hyper.f90:370 ``initialized`` guard). iso/aniso filters use the current dt.
* ``damp_zonal_only`` tests ``aky(ik) == epsilon(0.0)`` (hyper.f90:776), so the
  iso filter is applied to (practically) no mode and ``hypervisc_filter`` is
  left at its allocation value 1.0 (hyper.f90:613).
* simple3D "zonal" means ``ik == 1`` by index (hyper.f90:433-437), and its
  ``kx_cut`` default uses ``akx((ntheta0+1)/2)`` while the visc/res normaliser uses
  ``akx(ntheta0/2+1)`` (spec 07 P13).

Arithmetic follows the Fortran association order; ``x**n`` with integer ``n``
uses the libgcc ``__powidf2`` multiplication sequence (``powi``), ``**0.5`` is
``jnp.power(x, 0.5)`` (libm ``pow``, not ``sqrt``: Apple libm ``pow(x, 0.5)`` differs from
the correctly rounded ``sqrt`` by 1 ulp for ~0.1% of inputs; XLA:CPU's eager ``power``
reproduces libm bitwise, the fused jit one to 1 ulp; spec 07 P3). ``exp`` is XLA's
(1 ulp from libm for ~2.5% of inputs). Line numbers refer to the
upstream sources at 4d8c94bcf (the exporter's ``public`` statements shift the
working-tree lines slightly).
"""

from __future__ import annotations

from typing import Any, NamedTuple, Optional

import jax.numpy as jnp
import numpy as np
from jax import lax

from .numerics import ob, rmul

# hyper.f90:27-31
HYPER_NONE = 1
HYPER_VISC = 2
HYPER_RES = 3
HYPER_BOTH = 4
HYPER_SIMPLE3D = 5

# hyper.f90:458-464 (text_option table of read_parameters)
OPTION_SWITCH = {"default": HYPER_NONE, "none": HYPER_NONE, "visc_only": HYPER_VISC,
                 "res_only": HYPER_RES, "both": HYPER_BOTH, "simple3D": HYPER_SIMPLE3D}

#: ``epsilon(0.0)`` with ``-fdefault-real-8``.
EPS = float(np.finfo(np.float64).eps)

_MISSING = object()


class _Empty:
    arrays: dict = {}
    attrs: dict = {}
    path = "/hyper (absent)"


_EMPTY_GROUP = _Empty()


# ----------------------------------------------------------------------------- helpers
def _lookup(grp, name: str, default: Any = _MISSING):
    """Variable or attribute ``name`` of an export group (exact, then case-insensitive)."""
    for table in (getattr(grp, "arrays", {}), getattr(grp, "attrs", {})):
        if name in table:
            return table[name]
    low = name.lower()
    for table in (getattr(grp, "arrays", {}), getattr(grp, "attrs", {})):
        for key, value in table.items():
            if key.lower() == low:
                return value
    if default is _MISSING:
        raise KeyError(f"{getattr(grp, 'path', '?')}/{name} not found in the export "
                       f"(arrays: {sorted(getattr(grp, 'arrays', {}))})")
    return default


def _as_bool(value) -> bool:
    if isinstance(value, (bytes, str)):
        return value.strip().lower() in ("1", "t", "true", ".true.")
    return bool(np.asarray(value).item())


def _as_float(value) -> float:
    return float(np.asarray(value).item())


def _as_int(value) -> int:
    return int(np.asarray(value).item())


def powi(x, n: int):
    """``x**n`` for integer ``n`` exactly as libgcc ``__powidf2`` (spec 07 section 3.1.5).

    ``y = (n odd ? x : 1); while (n >>= 1) { x = x*x; if (n odd) y = y*x; }``;
    ``1/y`` for negative ``n``. NumPy float64, no FMA. ``jax.lax.integer_pow``
    performs the same multiplications.
    """
    x = np.asarray(x, dtype=np.float64)
    m = int(n)
    u = -m if m < 0 else m
    y = x.copy() if u % 2 else np.ones_like(x)
    while True:
        u >>= 1
        if not u:
            break
        x = x * x
        if u % 2:
            y = y * x
    return 1.0 / y if m < 0 else y


def resolve_hyper_parameters(switch: int, *, D_hypervisc: float, D_hyperres: float,
                             D_hyper: float, D_hyper3D: float, P_hyper3D: float,
                             nexp: int, isotropic_shear: bool) -> dict:
    """Effective hyper parameters, mirroring ``read_parameters`` (hyper.f90:452-601).

    Returns ``dict(hyper_on, D_hypervisc, D_hyperres, nexp)``. The logic is
    idempotent, so it may be applied to raw namelist values or to the effective
    values GS2 exported. ``hyper_on`` starts ``.false.`` (hyper.f90:35).
    """
    hyper_on = False
    if not isotropic_shear and nexp != 2:          # 498-501
        nexp = 2
    if switch == HYPER_NONE:                        # 506-516
        if D_hyperres > 0.0:
            D_hyperres = -10.0
        if D_hypervisc > 0.0:
            D_hypervisc = -10.0
    elif switch == HYPER_VISC:                      # 518-534
        hyper_on = True
        if D_hyperres > 0.0:
            D_hyperres = -10.0
        if D_hypervisc < 0.0:
            hyper_on = False
    elif switch == HYPER_RES:                       # 536-552
        hyper_on = True
        if D_hyperres < 0.0:
            hyper_on = False
        if D_hypervisc > 0.0:
            D_hypervisc = -10.0
    elif switch == HYPER_BOTH:                      # 554-582
        hyper_on = True
        if not D_hyper < 0.0:                       # `if (D_hyper < 0.) ... else` (557-566)
            D_hyperres = D_hyper
            D_hypervisc = D_hyper
        if D_hypervisc < 0.0 and D_hyperres < 0.0:
            hyper_on = False
    elif switch == HYPER_SIMPLE3D:                  # 584-599
        hyper_on = True
        if P_hyper3D < 4.0:
            hyper_on = False
        if D_hyper3D < 0.0:
            hyper_on = False
    else:
        raise NotImplementedError(f"hyper_option_switch = {switch} is not a GS2 8.2.1 option")
    return dict(hyper_on=hyper_on, D_hypervisc=float(D_hypervisc),
                D_hyperres=float(D_hyperres), nexp=int(nexp))


def normalisers(akx, aky, nexp: int, gridnorm: bool) -> tuple[float, float, float, float]:
    """``(akx4_max, aky_max, aky4_max, akperp4_max)`` of ``init_hyper`` (hyper.f90:386-396)."""
    if not gridnorm:
        return 1.0, 1.0, 1.0, 1.0
    akx = np.asarray(akx, np.float64)
    aky = np.asarray(aky, np.float64)
    kxm = akx[akx.shape[0] // 2]                    # akx(ntheta0/2 + 1), integer division
    kym = aky[-1]                                    # aky(naky)
    akx4_max = float(powi(kxm, 2 * nexp))
    aky_max = float(kym)
    aky4_max = float(powi(kym, 2 * nexp))
    akperp4_max = float(powi(kxm * kxm + kym * kym, nexp))
    return akx4_max, aky_max, aky4_max, akperp4_max


def d_res_formula(akx, aky, nexp: int, D_hyperres: float, dt0: float,
                  akperp4_max: float) -> np.ndarray:
    """``D_res`` K2 ``(Y, X)`` of ``init_hyper`` (hyper.f90:400-411), frozen at ``dt0``.

    ``D_res(it,ik) = ((D_hyperres*dt0) * powi(aky*aky + akx*akx, nexp)) / akperp4_max``
    if ``D_hyperres > 0`` else 0 (only reached for visc/res/both; 0 otherwise).
    """
    akx = np.asarray(akx, np.float64)
    aky = np.asarray(aky, np.float64)
    if not D_hyperres > 0.0:
        return np.zeros((aky.size, akx.size))
    q = aky[:, None] * aky[:, None] + akx[None, :] * akx[None, :]
    return ((D_hyperres * dt0) * powi(q, nexp)) / akperp4_max


def simple3d_filter_formula(kperp2, gds2, gds21, gds22, shat: float, akx, aky,
                            ky_cut: float, kx_cut: float, isotropic_model: bool,
                            D_hyper3D: float, P_hyper3D: float, dt0: float):
    """simple3D ``hypervisc_filter`` F3 of ``init_hyper`` (hyper.f90:419-445) at ``dt0``.

    Returns ``(filter, kperp2_max, kperp2_max_zonal, ky_cut, kx_cut)`` (effective cuts).
    ``ik = 1`` (by index) uses ``kperp2_max_zonal`` (433-437).
    """
    kperp2 = np.asarray(kperp2, np.float64)
    akx = np.asarray(akx, np.float64)
    aky = np.asarray(aky, np.float64)
    if not ky_cut > 0.0:                                       # 421
        ky_cut = float(aky[-1])
    if not kx_cut > 0.0:                                       # 422: akx((ntheta0+1)/2)
        kx_cut = float(akx[(akx.size + 1) // 2 - 1])
    gds2 = np.asarray(gds2, np.float64)
    gds21 = np.asarray(gds21, np.float64)
    gds22 = np.asarray(gds22, np.float64)
    kperp2_max = float(np.max(((ky_cut * ky_cut) * gds2
                               + (((2.0 * ky_cut) * kx_cut) * gds21) / shat)
                              + ((kx_cut * kx_cut) * gds22) / (shat * shat)))  # 424-426
    if isotropic_model:                                        # 427-431
        kperp2_max_zonal = kperp2_max
    else:
        kperp2_max_zonal = float(np.max(((kx_cut * kx_cut) * gds22) / (shat * shat)))
    c = D_hyper3D * dt0
    p = P_hyper3D / 2.0
    filt = np.empty_like(kperp2)
    filt[0] = np.exp(-(c * np.power(kperp2[0] / kperp2_max_zonal, p)))     # 433-437
    filt[1:] = np.exp(-(c * np.power(kperp2[1:] / kperp2_max, p)))         # 439-444
    return filt, kperp2_max, kperp2_max_zonal, ky_cut, kx_cut


# ----------------------------------------------------------------------------- containers
class HyperParams(NamedTuple):
    """Static (hashable) effective hyper parameters (``/hyper``, ``/dist_fn``).

    ``faithful_sum`` selects the Fortran sequential order for the shear-rate
    sums over ``(ik, it)`` (``lax.scan``); ``False`` uses a tree ``jnp.sum``
    (differs by ~``naky*ntheta0*eps`` relative in S, spec 07 section 5).
    """

    switch: int
    hyper_on: bool
    D_hypervisc: float
    D_hyperres: float
    nexp: int
    const_amp: bool
    isotropic_shear: bool
    damp_zonal_only: bool
    omega_osc: float
    gridnorm: bool
    akx4_max: float
    aky_max: float
    aky4_max: float
    akperp4_max: float
    D_hyper3D: float
    P_hyper3D: float
    hyper_in_initialisation: bool = True
    faithful_sum: bool = True

    @property
    def applies_filter(self) -> bool:
        """True if ``hyper_diff`` modifies ``gnew`` (hyper.f90:652-667)."""
        if not self.hyper_on:
            return False
        if self.switch in (HYPER_BOTH, HYPER_RES, HYPER_VISC):
            return not self.D_hypervisc < 0.0
        return self.switch == HYPER_SIMPLE3D

    @classmethod
    def from_export(cls, ex, *, faithful_sum: bool = True) -> "HyperParams":
        # A missing /hyper group means the namelist defaults (hyper_option = 'default').
        g = ex["/hyper"] if "hyper" in ex.root.groups else _EMPTY_GROUP
        sw = _lookup(g, "hyper_option_switch", HYPER_NONE if g is _EMPTY_GROUP else None)
        if sw is None:
            option = str(_lookup(g, "hyper_option")).strip()
            if option not in OPTION_SWITCH:
                raise NotImplementedError(f"hyper_option = {option!r}")
            sw = OPTION_SWITCH[option]
        switch = _as_int(sw)
        visc_res = switch in (HYPER_BOTH, HYPER_RES, HYPER_VISC)
        s3d = switch == HYPER_SIMPLE3D
        # Parameters that matter for the active option must be exported; the others
        # fall back to the hyper_knobs defaults (hyper.f90:63-126).
        required = {"D_hypervisc": visc_res, "D_hyperres": visc_res,
                    "D_hyper": switch == HYPER_BOTH, "nexp": visc_res,
                    "isotropic_shear": visc_res, "const_amp": visc_res,
                    "damp_zonal_only": visc_res, "omega_osc": visc_res, "gridnorm": visc_res,
                    "D_hyper3D": s3d, "P_hyper3D": s3d}
        defaults = {"D_hypervisc": -10.0, "D_hyperres": -10.0, "D_hyper": -10.0, "nexp": 2,
                    "isotropic_shear": True, "const_amp": False, "damp_zonal_only": False,
                    "omega_osc": 0.4, "gridnorm": True, "D_hyper3D": -10.0, "P_hyper3D": 4.0}

        def raw(name):
            return _lookup(g, name, _MISSING if required.get(name) else defaults[name])

        def f(name, default=_MISSING):
            v = raw(name) if name in defaults else _lookup(g, name, default)
            return v if v is None else _as_float(v)

        def b(name, default=_MISSING):
            v = raw(name) if name in defaults else _lookup(g, name, default)
            return v if v is None else _as_bool(v)

        isotropic_shear = b("isotropic_shear")
        eff = resolve_hyper_parameters(
            switch, D_hypervisc=f("D_hypervisc"), D_hyperres=f("D_hyperres"),
            D_hyper=f("D_hyper"), D_hyper3D=f("D_hyper3D"), P_hyper3D=f("P_hyper3D"),
            nexp=_as_int(raw("nexp")), isotropic_shear=isotropic_shear)
        exported_on = b("hyper_on", None)
        if exported_on is not None and exported_on != eff["hyper_on"]:
            raise ValueError(f"/hyper/hyper_on = {exported_on} but read_parameters logic "
                             f"gives {eff['hyper_on']} (switch {switch})")
        gridnorm = b("gridnorm")
        akx = np.asarray(ex["/kt_grids/akx"], np.float64)
        aky = np.asarray(ex["/kt_grids/aky"], np.float64)
        rec = normalisers(akx, aky, eff["nexp"], gridnorm)
        norms = []
        for name, value in zip(("akx4_max", "aky_max", "aky4_max", "akperp4_max"), rec):
            # Only set by init_hyper for visc/res/both (hyper.f90:382-396).
            got = f(name, None) if visc_res else None
            # Exported effective value wins (it is what GS2 used); the recomputation
            # (BOTH, critique A.12) is the fallback and is checked in test_misc.
            norms.append(got if got is not None and np.isfinite(got) else value)
        akx4_max, aky_max, aky4_max, akperp4_max = norms
        hyper_in_init = _lookup(ex["/dist_fn"], "hyper_in_initialisation", True)
        params = cls(
            switch=switch, hyper_on=eff["hyper_on"], D_hypervisc=eff["D_hypervisc"],
            D_hyperres=eff["D_hyperres"], nexp=eff["nexp"], const_amp=b("const_amp"),
            isotropic_shear=isotropic_shear, damp_zonal_only=b("damp_zonal_only"),
            omega_osc=f("omega_osc"), gridnorm=gridnorm, akx4_max=akx4_max,
            aky_max=aky_max, aky4_max=aky4_max, akperp4_max=akperp4_max,
            D_hyper3D=f("D_hyper3D"), P_hyper3D=f("P_hyper3D"),
            hyper_in_initialisation=_as_bool(hyper_in_init), faithful_sum=faithful_sum)
        params.check_gates()
        return params

    def check_gates(self) -> None:
        """Division-by-zero normalisers are gated (spec 07 section 3.2.7, App. B)."""
        if not self.hyper_on:
            return
        if self.switch in (HYPER_BOTH, HYPER_RES, HYPER_VISC):
            if self.gridnorm and self.akperp4_max == 0.0:
                raise NotImplementedError("hyper: gridnorm with akperp4_max == 0 (gated)")
            if self.applies_filter and not self.isotropic_shear and (
                    self.akx4_max == 0.0 or self.aky_max == 0.0):
                raise NotImplementedError(
                    "hyper: isotropic_shear = F with akx4_max == 0 or aky_max == 0 (gated)")


class HyperData(NamedTuple):
    """Arrays for ``hyper_diff``; k-dependent factors are precomputed on the host
    in NumPy (exact IEEE float64, the Fortran operation order).

    * ``k4`` K2: ``(akx*akx + aky*aky) * (akx*akx + aky*aky)`` (shear-rate sums);
    * ``qn`` K2: ``powi(aky*aky + akx*akx, nexp)`` (filters);
    * ``x4`` ``(X,)``: ``(akx*akx)*(akx*akx)``; ``a4`` ``(Y,)``: ``(aky*aky)*(aky*aky)``;
    * ``fac`` ``(Y,)``: 1.0 for ik = 1, 0.5 otherwise; ``aky`` ``(Y,)``;
    * ``aky_zero`` ``(Y,)`` bool: ``aky == 0.0`` exactly (aniso branch, 737);
    * ``zonal_skip`` ``(Y,)`` bool: ``damp_zonal_only .and. .not. aky == epsilon(0.0)``;
    * ``D_res`` K2: frozen export (consumed by the source term, spec 02);
    * ``filter_static`` F3: the exported simple3D filter for ``switch == simple3D``,
      otherwise ones (``hypervisc_filter`` allocation value, hyper.f90:613). It is the
      filter GS2 holds whenever ``hyper_diff`` does not overwrite it.
    """

    k4: np.ndarray
    qn: np.ndarray
    x4: np.ndarray
    a4: np.ndarray
    fac: np.ndarray
    aky: np.ndarray
    aky_zero: np.ndarray
    zonal_skip: np.ndarray
    D_res: np.ndarray
    filter_static: np.ndarray

    @classmethod
    def from_export(cls, ex, params: Optional[HyperParams] = None) -> "HyperData":
        if params is None:
            params = HyperParams.from_export(ex)
        g = ex["/hyper"] if "hyper" in ex.root.groups else _EMPTY_GROUP
        akx = np.asarray(ex["/kt_grids/akx"], np.float64)
        aky = np.asarray(ex["/kt_grids/aky"], np.float64)
        naky, ntheta0, nth = aky.size, akx.size, ex.nth
        k2 = akx[None, :] * akx[None, :] + aky[:, None] * aky[:, None]
        q = aky[:, None] * aky[:, None] + akx[None, :] * akx[None, :]
        x2 = akx * akx
        a2 = aky * aky
        fac = np.full(naky, 0.5)
        fac[0] = 1.0
        d_res = _lookup(g, "D_res", None)
        if d_res is None:
            d_res = np.zeros((naky, ntheta0))
        filt = None
        if params.switch == HYPER_SIMPLE3D:
            filt = _lookup(g, "hypervisc_filter", None)
            if filt is None:   # fallback: recompute at the initial dt (07 section 3.2.2)
                tg = ex["/theta_grid"]
                filt = simple3d_filter_formula(
                    ex["/kt_grids/kperp2"], tg["gds2"], tg["gds21"], tg["gds22"],
                    _as_float(tg["shat"]), akx, aky, _as_float(_lookup(g, "ky_cut")),
                    _as_float(_lookup(g, "kx_cut")), _as_bool(_lookup(g, "isotropic_model")),
                    params.D_hyper3D, params.P_hyper3D,
                    _as_float(ex["/state_t0/code_dt"]))[0]
        if filt is None:
            filt = np.ones((naky, ntheta0, nth))
        return cls(
            k4=k2 * k2, qn=powi(q, params.nexp), x4=x2 * x2, a4=a2 * a2, fac=fac, aky=aky,
            aky_zero=aky == 0.0,
            zonal_skip=np.full(naky, bool(params.damp_zonal_only)) & ~(aky == EPS),
            D_res=np.asarray(d_res, np.float64), filter_static=np.asarray(filt, np.float64))


# ----------------------------------------------------------------------------- kernels
def hyper_applies(istep: int, params: HyperParams) -> bool:
    """``(istep == 0 .and. hyper_in_initialisation) .or. istep > 0`` (dist_fn.fpp:4073)."""
    return (istep == 0 and params.hyper_in_initialisation) or istep > 0


def _sequential_sum(terms, faithful: bool):
    """Sum over axis 0 starting from 0.0, left to right (a Fortran accumulation loop)."""
    if terms.shape[0] == 0:
        return jnp.zeros(terms.shape[1:], terms.dtype)
    if not faithful:
        return jnp.sum(terms, axis=0)

    def body(s, t):
        return s + t, None

    n = terms.shape[0]
    s, _ = lax.scan(body, jnp.zeros(terms.shape[1:], terms.dtype), terms,
                    unroll=max(1, min(8, n)))
    return s


def _phi_sq(phi):
    """``real(conjg(phi)*phi) = phi.re*phi.re + phi.im*phi.im`` (spec 07 section 0.4)."""
    return ob(phi.real * phi.real) + ob(phi.imag * phi.imag)


def _osc(s, params: HyperParams):
    """``0.5 * ( -omega_osc + (omega_osc**2 + 2*S)**0.5 )`` (hyper.f90:706, 717, 728)."""
    w = params.omega_osc
    return 0.5 * ((-w) + jnp.power((w * w) + 2.0 * s, 0.5))


def shear_rate_iso(phinew, data: HyperData, params: HyperParams):
    """Belli-Hammett isotropic shear rate ``S(theta)`` (hyper.f90:757-771), ``(T,)``.

    ``S = sum_{ik, it} ((|phi|^2 * (k2*k2)) * fac(ik))`` with ik outer, it inner,
    then ``S**0.5``; ``S = 1`` if ``const_amp``.
    """
    nth = phinew.shape[-1]
    if params.const_amp:
        return jnp.ones((nth,), jnp.float64)
    terms = ob(_phi_sq(phinew) * data.k4[:, :, None]) * data.fac[:, None, None]
    s = _sequential_sum(ob(terms).reshape(-1, nth), params.faithful_sum)
    return jnp.power(s, 0.5)


def shear_rates_aniso(phinew, data: HyperData, params: HyperParams):
    """``(shear_rate_nz, shear_rate_z, shear_rate_z_nz)`` of ``aniso_shear``
    (hyper.f90:694-728), each ``(T,)``."""
    nth = phinew.shape[-1]
    b = _phi_sq(phinew)
    t_k4 = ob(ob(b * data.k4[:, :, None]) * data.fac[:, None, None])
    t_a4 = ob(ob(b * data.a4[:, None, None]) * data.fac[:, None, None])
    s_nz = _sequential_sum(t_k4[1:].reshape(-1, nth), params.faithful_sum)    # ik = 2..naky
    s_z = _sequential_sum(t_k4[:1].reshape(-1, nth), params.faithful_sum)     # ik = 1
    s_znz = _sequential_sum(t_a4[1:].reshape(-1, nth), params.faithful_sum)   # ik = 2..naky
    return _osc(s_nz, params), _osc(s_z, params), _osc(s_znz, params)


def hyper_filter(phinew, code_dt, data: HyperData, params: HyperParams):
    """The ``hypervisc_filter`` F3 ``(Y, X, T)`` that ``hyper_diff`` applies, or
    ``None`` if it returns without touching ``gnew`` (hyper.f90:652-667)."""
    if not params.hyper_on:                                               # 652
        return None
    if params.switch in (HYPER_BOTH, HYPER_RES, HYPER_VISC):
        if params.D_hypervisc < 0.0:                                      # 657
            return None
        ddt = params.D_hypervisc * code_dt
        if params.isotropic_shear:                                        # iso_shear
            s = shear_rate_iso(phinew, data, params)
            arg = ddt * (ob(s[None, None, :] * data.qn[:, :, None]) / params.akperp4_max)
            filt = jnp.exp(-arg)                                          # 777-778
            # GS2 quirk: damp_zonal_only tests aky(ik) == epsilon(0.0) exactly, so the
            # filter is applied to (practically) no mode and hypervisc_filter keeps its
            # allocation value 1.0 there (hyper.f90:776, 613).
            filt = jnp.where(data.zonal_skip[:, None, None], 1.0, filt)
            return filt
        s_nz, s_z, s_znz = shear_rates_aniso(phinew, data, params)        # aniso_shear
        x4 = data.x4[:, None]                                             # (X, 1)
        arg_zero = ddt * (ob(s_znz[None, :] * x4) / params.akx4_max)      # 738-739, (X, T)
        term_nz = ob(s_nz[None, None, :] * data.qn[:, :, None]) / params.akperp4_max
        term_z = ob((ob(s_z[None, :] * x4) / params.akx4_max)[None, :, :]
                    * data.aky[:, None, None]) / params.aky_max
        arg_nz = ddt * (term_nz + term_z)                                 # 741-743
        arg = jnp.where(data.aky_zero[:, None, None], arg_zero[None, :, :], arg_nz)
        return jnp.exp(-arg)
    if params.switch == HYPER_SIMPLE3D:                                   # simple3Dfilter
        return jnp.asarray(data.filter_static)
    return None


def hyper_diff(gnew, phinew, code_dt, data: HyperData, params: HyperParams):
    """``hyper_diff(g0 = gnew, phi = phinew)`` (src/hyper.f90:620-786).

    ``gnew`` G7 complex ``(S, E, L, Y, X, 2, T)``; ``phinew`` F3 complex (the
    fields passed to ``timeadv``: the old fields in pass 1, the new ones in pass
    2); ``code_dt`` the current time step (Python float or traced scalar).
    Returns ``(gnew, hypervisc_filter)`` with ``gnew(:, isgn, iglo) *= filter(:, it, ik)``
    componentwise for both signs (iso 779-780, aniso 746-747, simple3D 681-682).
    The filter returned when nothing is applied is the one GS2 holds
    (``data.filter_static``). Applied to every cell, including kwork_filter and
    forbidden points (spec 07 section 3.2.7). jit with ``params`` static.
    """
    filt = hyper_filter(phinew, code_dt, data, params)
    if filt is None:
        return gnew, jnp.asarray(data.filter_static)
    return rmul(gnew, filt[None, None, None, :, :, None, :]), filt
