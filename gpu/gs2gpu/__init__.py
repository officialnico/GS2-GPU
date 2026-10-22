"""gs2gpu: a faithful JAX (float64) port of the GS2 time step.

Fortran GS2 (``bin/gs2_export``) is the front end: it parses the input, builds
geometry, grids and coefficients, and writes them to ``<run>.gs2export.nc``.
This package reads that file and advances the gyrokinetic system on a JAX
device with the same numerics as GS2 8.2.1. See ``gpu/spec/`` for the
contract every module implements.
"""

from . import config as _config  # noqa: F401  (enables float64 on import)

__version__ = "0.1.0"
