"""Zonal-flow removal (``remove_zonal_flows``, src/fields_implicit.f90:374-387) and
the field-line average it uses; spec 01 section 3.8, spec 07 section 3.12,
spec 04 section 3.7.

Applied after ``phinew = phinew + phi`` and before the second ``timeadv`` iff
``fields_knobs:remove_zonal_flows_switch`` (fields_implicit.f90:362; the local
and gf_local solvers call the same routine)::

    phi_avg = 0 ; for it: phi_avg(:, it, 1) = field_line_average(phinew(:, it, 1))
    phinew = phinew - phi_avg

``ik = 1`` is treated as zonal **by index** (fields_implicit.f90:384; spec 04 P24).
Only ``phinew`` changes. ``field_line_average`` is the Klein-compensated
trapezoid of theta_grid.f90:665-670 / integration.f90:42-71 /
summation.f90:71-103, ported literally (sequential ``lax.scan``; a plain or
pairwise sum would differ by O(eps*nth), spec 07 P19).
"""

from __future__ import annotations

from typing import NamedTuple

import jax.numpy as jnp
import numpy as np
from jax import lax

from .numerics import rmul


def klein_sum(a, axis: int = -1):
    """``klein_sum`` (externals/utils/summation.f90:71-103) along ``axis``.

    Second-order compensated (Kahan-Babuska-Klein) sum, strictly sequential in
    index order, vectorised over the other axes; result
    ``(the_result + compensation_1) + compensation_2``. Additions only.
    """
    a = jnp.moveaxis(jnp.asarray(a), axis, 0)
    zero = jnp.zeros(a.shape[1:], a.dtype)
    if a.shape[0] == 0:
        return zero

    def body(carry, x):
        s, c1, c2 = carry
        t = s + x
        lc1 = jnp.where(jnp.abs(s) >= jnp.abs(x), (s - t) + x, (x - t) + s)
        s = t
        t = c1 + lc1
        lc2 = jnp.where(jnp.abs(c1) >= jnp.abs(lc1), (c1 - t) + lc1, (lc1 - t) + c1)
        c1 = t
        c2 = c2 + lc2
        return (s, c1, c2), None

    (s, c1, c2), _ = lax.scan(body, (zero, zero, zero), a)
    return (s + c1) + c2


def trapezoidal_integration_complex(dtheta, y):
    """``complex_trapezoidal_integration(x, y)`` (integration.f90:64-71) along the last axis.

    ``local_sum((x(2:) - x(:n-1)) * (y(:n-1) + y(2:))) * 0.5`` with
    ``local_sum = cmplx(klein_sum(real), klein_sum(aimag))`` (42-52) and every
    real*complex product componentwise. ``dtheta = x(2:) - x(:n-1)`` (precomputed).
    """
    v = rmul(y[..., :-1] + y[..., 1:], dtheta)
    total = lax.complex(klein_sum(v.real, axis=-1), klein_sum(v.imag, axis=-1))
    return rmul(total, 0.5)


def trapezoidal_integration_real(dtheta, y):
    """``real_trapezoidal_integration(x, y)`` (integration.f90:55-62), NumPy-free JAX."""
    v = dtheta * (y[..., :-1] + y[..., 1:])
    return klein_sum(v, axis=-1) * 0.5


class ZonalParams(NamedTuple):
    """``/fields/remove_zonal_flows_switch`` (static)."""

    remove_zonal_flows_switch: bool

    @classmethod
    def from_export(cls, ex) -> "ZonalParams":
        if "fields" not in ex.root.groups:                    # namelist default .false.
            return cls(remove_zonal_flows_switch=False)
        value = ex["/fields"].get("remove_zonal_flows_switch")
        if value is None:
            value = ex["/fields"].attrs.get("remove_zonal_flows_switch", False)
        return cls(remove_zonal_flows_switch=bool(np.asarray(value).item()))


class ZonalData(NamedTuple):
    """``jacob`` ``(T,)``, ``dtheta = theta(2:) - theta(:n-1)`` ``(T-1,)`` (host NumPy),
    and ``field_line_average_weight`` (0-d; exported, theta_grid.f90:520)."""

    jacob: np.ndarray
    dtheta: np.ndarray
    weight: np.ndarray

    @classmethod
    def from_export(cls, ex) -> "ZonalData":
        tg = ex["/theta_grid"]
        theta = np.asarray(tg["theta"], np.float64)
        return cls(jacob=np.asarray(tg["jacob"], np.float64), dtheta=theta[1:] - theta[:-1],
                   weight=np.asarray(float(tg["field_line_average_weight"]), np.float64))


def field_line_average_weight(data: ZonalData):
    """``1.0/trapezoidal_integration(theta, jacob)`` (theta_grid.f90:520), for validation."""
    return 1.0 / trapezoidal_integration_real(jnp.asarray(data.dtheta), jnp.asarray(data.jacob))


def field_line_average(q, data: ZonalData):
    """``field_line_average_complex`` (theta_grid.f90:665-670) along the last (theta) axis:
    ``trapezoidal_integration(theta, quantity*jacob) * field_line_average_weight``."""
    y = rmul(q, jnp.asarray(data.jacob))
    return rmul(trapezoidal_integration_complex(jnp.asarray(data.dtheta), y),
                jnp.asarray(data.weight))


def remove_zonal_flows(phinew, data: ZonalData, params: ZonalParams):
    """``remove_zonal_flows`` (src/fields_implicit.f90:374-387, ``fieldline_average_phi``
    398-425): subtract the field-line average of ``phinew(:, it, ik=1)`` from that
    row for every ``it``; other ky are unchanged (``x - 0 == x`` exactly).

    ``phinew`` F3 complex ``(Y, X, T)``. Returns ``phinew`` unchanged when
    ``params.remove_zonal_flows_switch`` is false (the caller may also skip the
    call, as GS2 does at fields_implicit.f90:362). jit with ``params`` static.
    """
    if not params.remove_zonal_flows_switch:
        return phinew
    avg = field_line_average(phinew[0], data)                      # (X,), ik = 1 by index
    return phinew.at[0].set(phinew[0] - avg[:, None])              # componentwise subtract
