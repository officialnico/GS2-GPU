"""Transform between GS2's ``g`` and the non-adiabatic ``h`` (``g_adjust``).

Mirrors ``g_adjust_floats`` and ``g_adjust_direction``
(src/dist_fn_arrays.fpp:101-170)::

    adj = ((bpar_factor*vperp2)*aj1)*bpar + ((phi_factor*phi)*aj0)
    g(:, 1) += adj ;  g(:, 2) += adj        (skipped where trapped .and. il > jend(ig))

with ``phi_factor = (zt(is)*nonmaxw_corr(ie,is))*facphi`` and
``bpar_factor = 2.0*facbpar``. ``direction = FROM_G`` (g -> h) uses
``(+fphi, +fbpar)``; ``TO_G`` (h -> g) uses ``(-fphi, -fbpar)``.
Used by the nonlinear term, collisions and diagnostics.
"""

from __future__ import annotations

from typing import NamedTuple

import jax.numpy as jnp
import numpy as np

from .io.export import Export
from .numerics import rmul

FROM_G = "from_g_gs2"
TO_G = "to_g_gs2"


class GAdjustData(NamedTuple):
    """Inputs of ``g_adjust``.

    ``vperp2, aj0, aj1``: G6 ``(S, E, L, Y, X, T)``; ``zt``: ``(S,)``;
    ``nonmaxw_corr``: ``(S, E)``; ``adjust_mask``: ``(L, T)`` bool, True where
    the point is adjusted (Fortran skips ``trapped .and. il > jend(ig)``).
    """

    vperp2: np.ndarray
    aj0: np.ndarray
    aj1: np.ndarray
    zt: np.ndarray
    nonmaxw_corr: np.ndarray
    adjust_mask: np.ndarray

    @classmethod
    def from_export(cls, ex: Export) -> "GAdjustData":
        arrays = ex["/dist_fn_arrays"]
        le = ex["/le_grids"]
        trapped = bool(np.asarray(le["grid_has_trapped_particles"]).item())
        jend = np.asarray(le["jend"])                       # (T,), Fortran values
        il = np.arange(1, ex.nlambda + 1)[:, None]          # Fortran il, (L, 1)
        skip = trapped & (il > jend[None, :])
        return cls(vperp2=arrays["vperp2"], aj0=arrays["aj0"], aj1=arrays["aj1"],
                   zt=ex["/species/spec_zt"], nonmaxw_corr=ex["/species/nonmaxw_corr"],
                   adjust_mask=~skip)


def adjustment(phi, bpar, facphi: float, facbpar: float, data: GAdjustData):
    """``adj`` of ``g_adjust_floats`` as a G6 complex array ``(S, E, L, Y, X, T)``."""
    phi_factor = (data.zt[:, None] * data.nonmaxw_corr) * facphi   # (S, E)
    bpar_factor = 2.0 * facbpar if facbpar != 0 else 0.0
    if facphi == 0:
        phi_factor = jnp.zeros_like(phi_factor)                    # Fortran: phi_factor = 0.0
    coef_b = (bpar_factor * data.vperp2) * data.aj1                # G6 real
    term_b = rmul(bpar[None, None, None], coef_b)
    term_p = rmul(rmul(phi[None, None, None], phi_factor[:, :, None, None, None, None]),
                  data.aj0)
    return term_b + term_p


def g_adjust(g, phi, bpar, facphi: float, facbpar: float, data: GAdjustData):
    """``g_adjust_floats``: returns ``g + adj`` on both signs (masked points unchanged).

    ``g`` is G7 ``(S, E, L, Y, X, 2, T)``; ``phi``, ``bpar`` are F3 ``(Y, X, T)``.
    """
    adj = adjustment(phi, bpar, facphi, facbpar, data)[..., None, :]
    mask = jnp.asarray(data.adjust_mask)[None, None, :, None, None, None, :]
    return jnp.where(mask, g + adj, g)


def g_adjust_direction(g, phi, bpar, direction: str, fphi: float, fbpar: float,
                       data: GAdjustData):
    """``g_adjust_direction``: ``FROM_G`` uses (+fphi, +fbpar), ``TO_G`` (-fphi, -fbpar)."""
    if direction == TO_G:
        return g_adjust(g, phi, bpar, -fphi, -fbpar, data)
    if direction == FROM_G:
        return g_adjust(g, phi, bpar, fphi, fbpar, data)
    raise ValueError(f"unknown g_adjust direction {direction!r}")
