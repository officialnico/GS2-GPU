"""Parity projection ``enforce_parity`` (spec 07 section 3.6, spec 03 section 3.12.2).

Mirrors ``enforce_parity`` (src/dist_fn.fpp:4796-4853), ``init_enforce_parity``
(4603-4738) and ``get_parity_conn`` (4742-4793). Line numbers refer to GS2
8.2.1, commit 4d8c94bcf (the specs' numbering).

``mult = +1`` if ``even`` else ``-1`` (an integer, so ``mult*z`` is exact).

* Linked boundaries (flux tube): the redistribute ``parity_redist`` is a
  scatter from the ``to`` list ``(ig, 2, iglo_conn)`` into the ``from`` list
  ``(-ig, 1, iglo)`` for **every** ``ig`` (including 0) and every ``iglo``
  (all pitch angles, kwork_filter cells included), i.e.
  ``gnew(-ig, 1, iglo) = gnew(ig, 2, conn(iglo))``; then, if ``mult /= 1``,
  ``gnew(:, 1, :) = mult*gnew(:, 1, :)``. ``conn`` is the mirror cell of the
  supercell chain (``get_parity_conn``). Source (sign 2) and destination
  (sign 1) never overlap, so it is a simultaneous gather.
* Otherwise (ballooning space): for ``j = 1..ntgrid``
  ``gnew(-j, 1) = mult*gnew(j, 2)`` and ``gnew(j, 1) = mult*gnew(-j, 2)``;
  ``gnew(0, 1)`` and sign 2 are untouched.

Parity can write non-zero values into forbidden points of an up-down
asymmetric geometry; this is replicated (no masking).
"""

from __future__ import annotations

from typing import NamedTuple

import jax.numpy as jnp
import numpy as np


class ParityData(NamedTuple):
    """``par_it``: ``(Y, X)`` int32, 0-based it index of the parity partner
    cell of ``(ik, it)`` (identity when the boundaries are not linked)."""

    par_it: np.ndarray

    @classmethod
    def from_export(cls, ex) -> "ParityData":
        """Built from the exported ``/dist_fn`` links (``l_links``, ``r_links``,
        ``itleft``, ``itright``) and ``boundary_option_switch``."""
        from .invert_rhs import BOUNDARY_LINKED, _opt_array  # local import: no cycle at load

        dfn = ex["/dist_fn"]
        linked = int(dfn["boundary_option_switch"]) == BOUNDARY_LINKED
        shape = (ex.naky, ex.ntheta0)
        l_links = _opt_array(dfn, "l_links", np.zeros(shape, np.int32))
        r_links = _opt_array(dfn, "r_links", np.zeros(shape, np.int32))
        itleft = _opt_array(dfn, "itleft", None)
        itright = _opt_array(dfn, "itright", None)
        return build_parity(l_links, r_links, itleft, itright, linked)


def parity_partner(l_links, r_links, itleft, itright, linked: bool) -> np.ndarray:
    """``it`` index (0-based) of ``get_parity_conn`` for every ``(ik, it)``.

    Literal walk of src/dist_fn.fpp:4742-4793 on the host: equal link counts
    keep the cell; more links on the left walk ``itleft`` until the reached
    cell has ``r_links == l_links(it)``; otherwise walk ``itright`` until
    ``l_links == r_links(it)``. Non-linked runs have no connections, so every
    cell is its own partner (and ``enforce_parity`` does not use it there).
    """
    l_links = np.asarray(l_links)
    r_links = np.asarray(r_links)
    naky, ntheta0 = l_links.shape
    par = np.tile(np.arange(ntheta0, dtype=np.int32), (naky, 1))
    if not linked:
        return par
    itleft = np.asarray(itleft)
    itright = np.asarray(itright)
    for ik in range(naky):
        for it in range(ntheta0):
            nl, nr = int(l_links[ik, it]), int(r_links[ik, it])
            if nl == nr:
                continue
            tmp = it
            if nl > nr:
                for _ in range(nl):
                    conn = int(itleft[ik, tmp]) - 1           # 1-based Fortran value
                    if conn < 0:
                        raise ValueError(f"get_parity_conn: broken left link at ik={ik + 1} "
                                         f"it={tmp + 1}")
                    tmp = conn
                    if int(r_links[ik, conn]) == nl:
                        break
            else:
                for _ in range(nr):
                    conn = int(itright[ik, tmp]) - 1
                    if conn < 0:
                        raise ValueError(f"get_parity_conn: broken right link at ik={ik + 1} "
                                         f"it={tmp + 1}")
                    tmp = conn
                    if int(l_links[ik, conn]) == nr:
                        break
            par[ik, it] = tmp
    return par


def build_parity(l_links, r_links, itleft, itright, linked: bool) -> ParityData:
    """Host-side construction of :class:`ParityData` (init_enforce_parity, 4603-4738)."""
    return ParityData(par_it=parity_partner(l_links, r_links, itleft, itright, linked))


def enforce_parity(gnew, data, params):
    """Mirrors ``enforce_parity`` (src/dist_fn.fpp:4796-4853).

    ``gnew``: G7 ``(S, E, L, Y, X, 2, T)``. ``data``: an ``InvertData`` (uses
    ``data.parity``) or a :class:`ParityData`. ``params``: ``InvertParams``
    (uses ``def_parity``, ``even``, ``linked``, ``ntgrid``). Returns ``gnew``
    unchanged when ``def_parity`` is false. Called at the end of
    ``invert_rhs_1`` (dist_fn.fpp:5706) and in ``timeadv`` after the
    collisions (dist_fn.fpp:4091).
    """
    if not params.def_parity:
        return gnew
    pd = data.parity if hasattr(data, "parity") else data
    gnew = jnp.asarray(gnew)
    # GS2 quirk: no forbid mask, so an up-down asymmetric geometry gets parity-mapped values
    # in forbidden points (dist_fn.fpp:4819-4851; 07 §3.6)
    g2 = gnew[..., 1, :]                                    # (S, E, L, Y, X, T)
    if params.linked:
        par_it = jnp.asarray(pd.par_it)                     # (Y, X)
        iy = jnp.arange(par_it.shape[0])[:, None]
        # gnew(-ig, 1, iglo) = gnew(ig, 2, conn(iglo)) for every ig (incl. 0)
        g1 = g2[:, :, :, iy, par_it, :][..., ::-1]
        if not params.even:
            g1 = -g1                                        # mult*gnew(:,1,:), mult = -1
    else:
        nt = params.ntgrid
        flipped = g2[..., ::-1]
        if not params.even:
            flipped = -flipped
        # gnew(0, 1) is not touched (dist_fn.fpp:4848-4849 loop over ig = 1..ntgrid)
        g1 = jnp.concatenate([flipped[..., :nt], gnew[..., 0, nt:nt + 1],
                              flipped[..., nt + 1:]], axis=-1)
    return jnp.stack([g1, g2], axis=-2)
