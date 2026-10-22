"""Arithmetic helpers that mirror how gfortran lowers complex arithmetic.

Spec 02 section 3.0 and spec 03 section 3.0. GS2 is double precision, so all
of these operate on float64/complex128.

* ``rmul(z, r)``: Fortran ``z*r`` or ``r*z`` with ``r`` real (componentwise).
* ``imul(w, z)``: Fortran ``(zi*w)*z`` with ``w`` real; ``zi*w == (0, w)``.
* ``cmul(a, b)``: general complex times complex.
* ``cdiv(a, b)``: complex division with GCC's Smith method (``-fcx-fortran-rules``
  is not in effect; gfortran uses the wide Smith expansion).

"Strict" mode wraps every product in ``lax.optimization_barrier`` so XLA
cannot contract ``a*b + c`` into a fused multiply-add. With the reference
Fortran build (``-ffp-contract=off``) this gives bitwise agreement on XLA:CPU
(spec 03 section A.2). Strict mode is for tests only; it stops XLA from
fusing kernels, so it is slow on a GPU. Switch it with ``set_strict``
*before* tracing/jitting the code that should use it.
"""

from __future__ import annotations

import jax.numpy as jnp
from jax import lax

_STRICT = False


def set_strict(flag: bool) -> None:
    """Enable/disable strict (FMA-free) products for functions traced afterwards."""
    global _STRICT
    _STRICT = bool(flag)


def strict() -> bool:
    return _STRICT


def ob(x):
    """Optimisation barrier in strict mode, identity otherwise."""
    return lax.optimization_barrier(x) if _STRICT else x


def rmul(z, r):
    """Fortran ``z*r`` / ``r*z`` with real ``r`` and complex ``z``."""
    return lax.complex(ob(z.real * r), ob(z.imag * r))


def imul(w, z):
    """Fortran ``(zi*w)*z`` with real ``w``: ``(0, w) * z``."""
    return lax.complex(-ob(w * z.imag), ob(w * z.real))


def cmul(a, b):
    """General Fortran complex*complex."""
    return lax.complex(ob(a.real * b.real) - ob(a.imag * b.imag),
                       ob(a.real * b.imag) + ob(a.imag * b.real))


def cdiv(a, b):
    """Fortran complex/complex, GCC's Smith method (ties take the second branch).

    Branches not taken may contain inf/nan (division by a zero component);
    ``jnp.where`` discards them.
    """
    ar, ai, br, bi = a.real, a.imag, b.real, b.imag
    t = jnp.abs(br) < jnp.abs(bi)
    r1 = br / bi
    d1 = ob(br * r1) + bi
    x1 = ob(ar * r1) + ai
    y1 = ob(ai * r1) - ar
    r2 = bi / br
    d2 = ob(bi * r2) + br
    x2 = ob(ai * r2) + ar
    y2 = ai - ob(ar * r2)
    return lax.complex(jnp.where(t, x1 / d1, x2 / d2), jnp.where(t, y1 / d1, y2 / d2))


def recip(z):
    """Fortran ``1.0/z`` with complex ``z`` (real numerator promoted to complex)."""
    one = lax.complex(jnp.ones_like(z.real), jnp.zeros_like(z.real))
    return cdiv(one, z)


def radd(r, z):
    """Fortran ``r + z`` with real ``r``: only the real part changes."""
    return lax.complex(r + z.real, z.imag)


def rsub(r, z):
    """Fortran ``r - z`` with real ``r``: ``(r - z.re, -z.im)``."""
    return lax.complex(r - z.real, -z.imag)


def ith(ig: int, ntgrid: int) -> int:
    """Theta array index (0-based, NumPy) of Fortran ``ig`` in ``-ntgrid:ntgrid``."""
    return ig + ntgrid


def max_abs_rel(ref, got) -> float:
    """``max|got - ref| / max|ref|`` (0 if both are identically 0)."""
    import numpy as np

    ref = np.asarray(ref)
    got = np.asarray(got)
    scale = float(np.max(np.abs(ref))) if ref.size else 0.0
    diff = float(np.max(np.abs(got - ref))) if ref.size else 0.0
    if scale == 0.0:
        return diff
    return diff / scale
