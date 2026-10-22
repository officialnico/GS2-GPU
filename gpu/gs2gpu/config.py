"""JAX configuration shared by every module.

GS2 is built with ``-fdefault-real-8``: every real is float64 and every
complex is complex128. JAX defaults to 32-bit, so x64 is enabled here, on
import of ``gs2gpu``, before any array is created.
"""

import jax

jax.config.update("jax_enable_x64", True)

REAL = "float64"
COMPLEX = "complex128"


def device(kind: str | None = None):
    """Return the first JAX device of ``kind`` ('gpu', 'cpu'), or the default."""
    if kind is None:
        return jax.devices()[0]
    return jax.devices(kind)[0]
