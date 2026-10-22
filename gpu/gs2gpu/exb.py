"""ExB flow shear: ``exb_shear`` (src/dist_fn.fpp:4183-4600), spec 07 section 3.3,
spec 01 section 3.7.

Called once per ``advance(istep)``, before ``g = gnew`` and the first ``timeadv``
(fields_implicit.f90:335, fields_local.f90:4280), on ``gnew, phinew, aparnew,
bparnew``. It is active iff ``kx_shift`` or ``theta0_shift`` is allocated, i.e.
``abs(g_exb*g_exbfac) > epsilon(0.)`` (dist_fn.fpp:3983-3991).

Split, as decided (00_decisions D9, critique G-34):

* host (plain Python float64, the Fortran operation order): ``ExbState`` holds
  ``kx_shift``/``theta0_shift``, ``jump``, ``ikx_indexed``, ``exb_first``,
  ``istep_last``, ``dkx``/``dtheta0``/``itmin``; ``exb_update`` mirrors the setup
  block, the guards and the shift bookkeeping (4235-4314);
* device (jit, ``params`` static, ``jump`` a traced int32 ``(naky,)`` array):
  ``exb_shift`` is the permutation + zero fill of ``gnew`` and the fields along kx
  (box or ``shat == 0``, 4490-4599) or theta0 (ballooning, 4341-4489), written as a
  simultaneous gather. ``jump == 0`` everywhere is the identity.

GS2 quirks replicated (00_decisions D4):

* ``code_dt_old`` starts at 0 (gs2_time.f90:31; not saved in restarts), so the
  first applied shear uses ``gdt = code_dt/2`` (4283);
* ``NINT`` rounds half away from zero (4302, 4308);
* ``istep_last`` (SAVE, 4217) blocks a second shift in the retry ``advance(istep)``
  after an immediate reset, and is only updated when the start conditions pass
  (4270-4286); the ``ExbState`` must be carried across resets unchanged;
* the kx branch never shifts ``ik = 1`` (``do ik = naky, 2, -1``, 4492), but
  ``jump(1)`` is still computed and enters ``all(jump == 0)``;
* the theta0 branch zeroes ``phi(iright+1:ntgrid, :iib, :)`` for **all** ky (4370,
  4387, 4404), in ``ik`` order, and aborts if ``|jump| >= ntheta0`` (4347-4351).

Line numbers refer to the upstream sources at 4d8c94bcf.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, replace
from typing import NamedTuple, Optional

import jax
import jax.numpy as jnp
import numpy as np
from jax import lax

from .types import Fields, TimeBlock

#: ``epsilon(0.)`` with ``-fdefault-real-8``.
EPS = float(np.finfo(np.float64).eps)
TWOPI = 2.0 * math.pi  # constants::twopi = 2.0*pi with pi = 3.1415926535897931


def _get(grp, name, default=None):
    """Variable/attribute of an export group, exact then case-insensitive."""
    for table in (getattr(grp, "arrays", {}), getattr(grp, "attrs", {})):
        if name in table:
            return table[name]
    low = name.lower()
    for table in (getattr(grp, "arrays", {}), getattr(grp, "attrs", {})):
        for key, value in table.items():
            if key.lower() == low:
                return value
    return default


def _bool(value) -> bool:
    if isinstance(value, (bytes, str)):
        return value.strip().lower() in ("1", "t", "true", ".true.")
    return bool(np.asarray(value).item())


def fortran_nint(x: float) -> int:
    """Fortran ``NINT`` for real(8): round half away from zero (spec 07 section 3.1.4).

    ``x - trunc(x)`` is exact for every finite double, so ties are detected exactly
    (``floor(|x| + 0.5)`` would misround 0.49999999999999994).
    """
    x = float(x)
    if not math.isfinite(x):
        raise ValueError(f"NINT of non-finite value {x!r}")
    t = math.trunc(x)
    frac = x - t
    if abs(frac) >= 0.5:
        return int(t) + (1 if x > 0 else -1)
    return int(t)


def ikx_indexed_of(ikx) -> tuple[np.ndarray, int]:
    """``ikx_indexed`` (Fortran 1-based it values) and ``itmin`` (dist_fn.fpp:4238-4249).

    ``itmin = minloc(ikx)`` (first occurrence); ``ikx_indexed(p)`` lists it in the
    order ``itmin, ..., ntheta0, 1, ..., itmin-1`` (increasing kx for box grids).
    """
    ikx = np.asarray(ikx)
    nt = ikx.size
    itmin = int(np.argmin(ikx)) + 1                 # minloc: first occurrence, 1-based
    out = np.empty(nt, np.int64)
    for it in range(itmin, nt + 1):                 # 4243-4245
        out[it + 1 - itmin - 1] = it
    for it in range(1, itmin):                      # 4247-4249
        out[nt - itmin + 1 + it - 1] = it
    return out, itmin


# ----------------------------------------------------------------------------- containers
@dataclass(frozen=True)
class ExbParams:
    """Static (hashable) inputs of ``exb_shear``; host copies of the small grid
    vectors are kept as tuples so the bookkeeping never touches the device."""

    active: bool                 # kx_shift/theta0_shift allocated
    use_kx: bool                 # should_use_kx = box .or. (shat .eq. 0.0)   (4228)
    naky: int
    ntheta0: int
    ntgrid: int
    nth: int
    ntheta: int
    has_phi: bool
    has_apar: bool
    has_bpar: bool
    g_exb: float
    g_exbfac: float
    g_exb_start_time: float
    g_exb_start_timestep: int
    shat: float
    aky: tuple
    akx: tuple
    ikx: tuple
    theta0_k1: tuple             # theta0(:, 1) (Fortran theta0(it, ik=1))
    tunits: tuple

    @property
    def ileft(self) -> int:
        """Fortran ``ileft = -ntgrid + ntheta`` (4354)."""
        return -self.ntgrid + self.ntheta

    @property
    def iright(self) -> int:
        """Fortran ``iright = ntgrid - ntheta`` (4354)."""
        return self.ntgrid - self.ntheta

    @classmethod
    def from_export(cls, ex) -> "ExbParams":
        dist_fn, kt, tg, rp = ex["/dist_fn"], ex["/kt_grids"], ex["/theta_grid"], \
            ex["/run_parameters"]
        g_exb = float(_get(dist_fn, "g_exb", 0.0))
        g_exbfac = float(_get(dist_fn, "g_exbfac", 1.0))
        state = ex["/state_t0"]
        allocated = ("kx_shift" in state) or ("theta0_shift" in state)
        active = abs(g_exb * g_exbfac) > EPS                         # 3983
        if allocated != active:
            raise ValueError(f"exb: kx_shift/theta0_shift allocated = {allocated} but "
                             f"abs(g_exb*g_exbfac) > epsilon is {active}")
        box = _bool(kt["box"])
        shat = float(tg["shat"])
        aky = np.asarray(kt["aky"], np.float64)
        akx = np.asarray(kt["akx"], np.float64)
        theta0 = np.asarray(kt["theta0"], np.float64)               # K2 (Y, X)
        tunits = _get(ex["/gs2_time"], "tunits") if "gs2_time" in ex.root.groups else None
        if tunits is None:                                           # gs2_time.f90:108-124
            wstar_units = _bool(_get(rp, "wstar_units", False))
            tunits = np.where(aky != 0.0, 2.0 / np.where(aky != 0.0, aky, 1.0), 0.0) \
                if wstar_units else np.ones_like(aky)
        ntheta = _get(tg, "ntheta")
        if ntheta is None:
            raise KeyError("/theta_grid/ntheta missing from the export")
        params = cls(
            active=active, use_kx=bool(box or shat == 0.0), naky=ex.naky,
            ntheta0=ex.ntheta0, ntgrid=ex.ntgrid, nth=ex.nth, ntheta=int(ntheta),
            has_phi=_bool(rp["has_phi"]), has_apar=_bool(rp["has_apar"]),
            has_bpar=_bool(rp["has_bpar"]), g_exb=g_exb, g_exbfac=g_exbfac,
            g_exb_start_time=float(_get(dist_fn, "g_exb_start_time", -1.0)),
            g_exb_start_timestep=int(_get(dist_fn, "g_exb_start_timestep", -1)),
            shat=shat, aky=tuple(float(v) for v in aky), akx=tuple(float(v) for v in akx),
            ikx=tuple(int(v) for v in np.asarray(kt["ikx"])),
            theta0_k1=tuple(float(v) for v in theta0[0]),
            tunits=tuple(float(v) for v in np.asarray(tunits, np.float64)))
        if params.active and params.ntheta0 == 1:
            # Gated (critique B.6, spec 01 section 6.24): the kx branch never sets dkx
            # (4251-4257); the theta0 branch (dtheta0 = 2*pi) aborts at |theta0_shift| ~ pi.
            raise NotImplementedError("g_exb /= 0 with ntheta0 == 1 (gated)")
        return params


class ExbData(NamedTuple):
    """Device gather tables of the kx branch: ``ikx_indexed0[p]`` = 0-based it at
    ascending-kx position ``p``; ``pos[it]`` = its position (inverse permutation)."""

    ikx_indexed0: np.ndarray
    pos: np.ndarray

    @classmethod
    def from_params(cls, params: ExbParams) -> "ExbData":
        idx, _ = ikx_indexed_of(params.ikx)
        idx0 = (idx - 1).astype(np.int32)
        pos = np.empty_like(idx0)
        pos[idx0] = np.arange(idx0.size, dtype=np.int32)
        return cls(ikx_indexed0=idx0, pos=pos)

    @classmethod
    def from_export(cls, ex) -> "ExbData":
        return cls.from_params(ExbParams.from_export(ex))


@dataclass(frozen=True)
class ExbState:
    """Host-side ExB bookkeeping (dist_fn.fpp:4212-4219, dist_fn_arrays kx_shift /
    theta0_shift). Immutable; ``exb_update`` returns a new state.

    ``jump``/``ikx_indexed``/``dkx``/``dtheta0``/``itmin`` are None until the
    ``exb_first`` block has run. ``ikx_indexed`` and ``itmin`` keep Fortran 1-based
    values. Carry the same state across a time-step reset (nothing here is
    re-initialised by ``reset_time_step``, spec 07 section 1.4).
    """

    kx_shift: Optional[np.ndarray] = None
    theta0_shift: Optional[np.ndarray] = None
    jump: Optional[np.ndarray] = None
    ikx_indexed: Optional[np.ndarray] = None
    exb_first: bool = True
    istep_last: int = 0
    dkx: Optional[float] = None
    dtheta0: Optional[float] = None
    itmin: Optional[int] = None

    @classmethod
    def from_group(cls, grp, params: ExbParams, *, istep_last: int = 0) -> "ExbState":
        """State from a group holding ``kx_shift``/``theta0_shift`` and ``exb_first``
        (``/state_t0`` or a trace record). If ``exb_first`` is already false the
        (deterministic) setup block is re-run to fill ``ikx_indexed``/``dkx``/``dtheta0``;
        ``jump`` is taken from the group when present."""
        kx = _get(grp, "kx_shift")
        th = _get(grp, "theta0_shift")
        first = _get(grp, "exb_first")
        state = cls(kx_shift=None if kx is None else np.array(kx, np.float64, ndmin=1),
                    theta0_shift=None if th is None else np.array(th, np.float64, ndmin=1),
                    exb_first=True if first is None else _bool(first),
                    istep_last=int(istep_last))
        if params.active and state.kx_shift is None and state.theta0_shift is None:
            raise KeyError(f"exb: neither kx_shift nor theta0_shift in {getattr(grp, 'path', '?')}")
        if not state.exb_first:
            state = _setup(state, params)
            jump = _get(grp, "jump")
            if jump is not None:
                state = replace(state, jump=np.array(jump, np.int64, ndmin=1))
        return state

    @classmethod
    def from_export(cls, ex, params: Optional[ExbParams] = None) -> "ExbState":
        """Initial state (``/state_t0``): ``istep_last = 0`` (SAVE initial value, 4217)."""
        if params is None:
            params = ExbParams.from_export(ex)
        return cls.from_group(ex["/state_t0"], params, istep_last=0)

    @property
    def shift(self) -> Optional[np.ndarray]:
        """``kx_shift`` (kx branch) or ``theta0_shift`` (theta0 branch)."""
        return self.kx_shift if self.kx_shift is not None else self.theta0_shift


# ----------------------------------------------------------------------------- host part
def _setup(state: ExbState, params: ExbParams) -> ExbState:
    """The ``exb_first`` block (dist_fn.fpp:4235-4267)."""
    jump = np.zeros(params.naky, np.int64)                           # 4237-4238
    if params.use_kx:
        idx, itmin = ikx_indexed_of(params.ikx)                      # 4239-4249
        dkx = params.akx[1] - params.akx[0] if params.ntheta0 > 1 else None   # 4251-4257
        return replace(state, exb_first=False, jump=jump, ikx_indexed=idx, itmin=itmin,
                       dkx=dkx)
    if params.ntheta0 > 1:                                           # 4261-4265
        dtheta0 = params.theta0_k1[1] - params.theta0_k1[0]          # theta0(2,1)-theta0(1,1)
    else:
        dtheta0 = TWOPI
    return replace(state, exb_first=False, jump=jump, dtheta0=dtheta0)


def exb_update(state: ExbState, istep: int, time: TimeBlock,
               params: ExbParams) -> tuple[ExbState, Optional[np.ndarray]]:
    """Host bookkeeping of ``exb_shear`` (dist_fn.fpp:4235-4314), literal float64.

    Returns the new state and the ``jump`` array (int, ``(naky,)``) if a shift must
    be applied, else ``None`` (inactive, guarded out, or ``all(jump == 0)``).
    ``time`` is the time block at the start of the step (``update_time`` runs after
    ``advance``). Raises ``RuntimeError`` where GS2 calls ``mp_abort`` (theta0
    branch with ``|jump| >= ntheta0``, 4347-4351).
    """
    if not params.active:            # exb_shear is not called (fields_implicit.f90:335)
        return state, None
    if state.exb_first:
        state = _setup(state, params)
    if istep == state.istep_last:                                    # 4270
        return state, None
    if params.g_exb_start_timestep > istep:                          # 4271
        return state, None
    if params.g_exb_start_time >= 0 and time.code_time < params.g_exb_start_time:  # 4272
        return state, None
    gdt = 0.5 * (float(time.code_dt) + float(time.code_dt_old))      # 4283
    jump = np.array(state.jump, np.int64, copy=True)
    g_exb, g_exbfac = params.g_exb, params.g_exbfac
    if params.use_kx:                                                # 4299-4304
        if state.dkx is None:
            raise NotImplementedError("exb_shear: dkx unset (ntheta0 == 1, gated)")
        dkx = state.dkx
        kx = np.array(state.kx_shift, np.float64, copy=True)
        for ik in range(params.naky):
            v = float(kx[ik]) - ((((params.aky[ik] * g_exb) * g_exbfac) * gdt)
                                 * params.tunits[ik])
            j = fortran_nint(v / dkx)
            jump[ik] = j
            kx[ik] = v - float(j) * dkx
        state = replace(state, kx_shift=kx, jump=jump, istep_last=istep)   # 4286
    else:                                                            # 4305-4310
        dth = state.dtheta0
        th = np.array(state.theta0_shift, np.float64, copy=True)
        for ik in range(params.naky):
            v = float(th[ik]) - ((((g_exb * g_exbfac) * gdt) / params.shat)
                                 * params.tunits[ik])
            j = fortran_nint(v / dth)
            jump[ik] = j
            th[ik] = v - dth * float(j)
        state = replace(state, theta0_shift=th, jump=jump, istep_last=istep)
    if np.all(jump == 0):                                            # 4314
        return state, None
    if not params.use_kx:
        for ik in range(params.naky):                                # 4344-4351
            if jump[ik] != 0 and abs(int(jump[ik])) >= params.ntheta0:
                raise RuntimeError(
                    f"in exb_shear: jump(ik)={int(jump[ik])} > ntheta0 ={params.ntheta0} "
                    f"for ik={ik + 1}. => reduce timestep or increase ntheta0")
    return state, jump.copy()


# ----------------------------------------------------------------------------- device part
def _kx_tables(jump, data: ExbData, params: ExbParams):
    """Source it (0-based) and validity of ``X_new(:, it, ik)`` in the kx branch.

    ``q = pos(it) - jump(ik)``; ``X_new(ikx_indexed(p)) = X_old(ikx_indexed(p - jump))``
    for ``1 <= p - jump <= ntheta0``, else 0 (4492-4598); ``ik = 1`` unchanged.
    """
    nt = params.ntheta0
    pos = jnp.asarray(data.pos)
    idx0 = jnp.asarray(data.ikx_indexed0)
    q = pos[None, :] - jump[:, None]                                  # (Y, X)
    valid = (q >= 0) & (q < nt)
    src = idx0[jnp.clip(q, 0, nt - 1)]
    first = (jnp.arange(params.naky) == 0)[:, None]                   # ik = 1 never shifted
    src = jnp.where(first, jnp.arange(nt, dtype=src.dtype)[None, :], src)
    return src, valid | first


def _theta0_tables(jump, params: ExbParams):
    """kx-axis source and the theta wrap selectors of the theta0 branch.

    With Fortran ``i`` and ``c = i - j`` (0-based ``c0 = i0 - j``):
    ``0 <= c0 < nt``: ``X_new(ig, i) = X_old(ig, c)``;
    ``c0 < 0`` (j > 0): ``X_old(ig + ntheta, c + nt)`` for ``ig <= iright``, else 0;
    ``c0 >= nt`` (j < 0): ``X_old(ig - ntheta, c - nt)`` for ``ig >= ileft``, else 0
    (dist_fn.fpp:4353-4484). Requires ``|j| < ntheta0`` (checked on the host).
    """
    nt = params.ntheta0
    c0 = jnp.arange(nt)[None, :] - jump[:, None]                      # (Y, X)
    lo = c0 < 0
    hi = c0 >= nt
    src = jnp.clip(jnp.where(lo, c0 + nt, jnp.where(hi, c0 - nt, c0)), 0, nt - 1)
    return src, lo, hi


def _shift_theta(x, k: int):
    """``out[..., t] = x[..., t + k]`` where ``0 <= t + k < nth``, else 0 (static ``k``)."""
    nth = x.shape[-1]
    if abs(k) >= nth:
        return jnp.zeros_like(x)
    pad = jnp.zeros(x.shape[:-1] + (abs(k),), x.dtype)
    if k >= 0:
        return jnp.concatenate([x[..., k:], pad], axis=-1)
    return jnp.concatenate([pad, x[..., :nth + k]], axis=-1)


def _gather_fields(x, src):
    """``x[ik, src[ik, it], :]`` for an F3 ``(Y, X, T)`` array."""
    iy = jnp.arange(x.shape[0])[:, None]
    return x[iy, src]


def _gather_g(g, src):
    """``g[..., ik, src[ik, it], :, :]`` for a G7 ``(S, E, L, Y, X, 2, T)`` array."""
    iy = jnp.arange(g.shape[3])[:, None]
    return g[:, :, :, iy, src]


def _theta0_apply(x1, lo, hi, params: ExbParams):
    """Theta wrap of the gathered array ``x1`` (theta last; lo/hi already broadcast)."""
    up = _shift_theta(x1, params.ntheta)       # source ig + ntheta (j > 0, c < 1)
    dn = _shift_theta(x1, -params.ntheta)      # source ig - ntheta (j < 0, c > ntheta0)
    return jnp.where(lo, up, jnp.where(hi, dn, x1))


def exb_shift(gnew, fields: Fields, jump, data: ExbData, params: ExbParams):
    """Device part of ``exb_shear`` (dist_fn.fpp:4341-4599): shift ``gnew`` and the
    fields along kx (``params.use_kx``) or theta0 by ``jump`` (int32 ``(naky,)``).

    Pure permutation + zero fill (no arithmetic on values). Fields are shifted only
    if ``has_phi/has_apar/has_bpar``; ``gnew`` always (all ``il, ie, is``, both signs,
    all theta). ``jump == 0`` gives the identity. jit with ``params`` static.
    Returns ``(gnew, fields)``.
    """
    jump = jnp.asarray(jump, jnp.int32)
    flags = (params.has_phi, params.has_apar, params.has_bpar)
    if params.use_kx:
        src, valid = _kx_tables(jump, data, params)
        out = []
        for x, flag in zip(fields, flags):
            out.append(jnp.where(valid[:, :, None], _gather_fields(x, src), 0)
                       if flag else x)
        gout = jnp.where(valid[None, None, None, :, :, None, None], _gather_g(gnew, src), 0)
        return gout, Fields(*out)

    src, lo, hi = _theta0_tables(jump, params)
    # GS2 quirk: for j > 0 the field statements zero phi(iright+1:ntgrid, :iib, :) for
    # ALL ky (dist_fn.fpp:4370, 4387, 4404), processed in ascending ik. Zeroing regions
    # are nested in it, so row ik sees the union of the zeroings of the ky processed
    # before it (prefix max of positive jumps), then its own gather, then those of
    # itself and later ky (suffix max). Equivalent to the sequential Fortran loop.
    jpos = jnp.where(jump > 0, jump, 0)
    m_after = lax.cummax(jpos, axis=0, reverse=True)                  # max over ik' >= ik
    m_before = jnp.concatenate([jnp.zeros((1,), jpos.dtype),
                                lax.cummax(jpos, axis=0)[:-1]])       # max over ik' < ik
    t = jnp.arange(params.nth)
    tail = t > (params.iright + params.ntgrid)                        # ig > iright
    it = jnp.arange(params.ntheta0)
    zero_before = (it[None, :, None] < m_before[:, None, None]) & tail[None, None, :]
    zero_after = (it[None, :, None] < m_after[:, None, None]) & tail[None, None, :]
    out = []
    for x, flag in zip(fields, flags):
        if not flag:
            out.append(x)
            continue
        x = jnp.where(zero_before, 0, x)
        x = _theta0_apply(_gather_fields(x, src), lo[:, :, None], hi[:, :, None], params)
        out.append(jnp.where(zero_after, 0, x))
    # gnew: per-ik shift with the theta wrap, no cross-ky effect (4406-4484).
    lo_g = lo[None, None, None, :, :, None, None]
    hi_g = hi[None, None, None, :, :, None, None]
    gout = _theta0_apply(_gather_g(gnew, src), lo_g, hi_g, params)
    return gout, Fields(*out)


_exb_shift_jit = jax.jit(exb_shift, static_argnames=("params",))


def exb_shear(gnew, fields: Fields, exb_state: ExbState, istep: int, time: TimeBlock,
              data: ExbData, params: ExbParams):
    """``exb_shear(g0 = gnew, phinew, aparnew, bparnew, istep)`` (src/dist_fn.fpp:4183-4600).

    Host bookkeeping (``exb_update``) then, only if some ``jump /= 0``, the jitted
    device shift (``exb_shift``). ``fields`` are ``phinew, aparnew, bparnew``;
    ``time`` the ``TimeBlock`` at the start of the step. Returns
    ``(gnew, fields, exb_state)``. To embed the shift in a larger jitted step, call
    ``exb_update`` on the host and ``exb_shift`` inside the jit instead.
    """
    exb_state, jump = exb_update(exb_state, istep, time, params)
    if jump is None:
        return gnew, fields, exb_state
    gnew, fields = _exb_shift_jit(gnew, Fields(*fields), jnp.asarray(jump, jnp.int32),
                                  data, params)
    return gnew, fields, exb_state
