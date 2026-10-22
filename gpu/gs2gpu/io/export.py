"""Reader for ``<run>.gs2export.nc`` and the groups of ``<run>.gs2trace.nc``.

Implements the loader contract of spec 09 section 6.6:

* netCDF stores every array in Fortran dimension order, so NumPy sees the
  dimensions reversed (canonical g: NumPy ``(S, E, L, Y, X, 2, T)``). Axes are
  never reordered here (spec 10 G-29).
* ``fortran_type == "complex"``: the trailing NumPy ``ri`` axis becomes
  complex128 (``a[..., 0] + 1j*a[..., 1]``).
* ``fortran_type == "logical"``: int8 0/1 becomes bool.
* reals stay float64; integers stay int32 with their Fortran values (1-based
  indices stay 1-based, -1 sentinels stay -1).
* a variable that GS2 did not allocate is not written, and its group's
  ``absent`` attribute lists it; asking for a name that is neither present nor
  listed as absent is an error (it was forgotten by the exporter).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterator

import netCDF4
import numpy as np

SUPPORTED_FORMAT_MAJOR = 1

# netCDF-4/HDF5 cannot hold a root dimension and a root group of the same
# name, so groups whose module name collides with a root dimension are stored
# with a trailing underscore; lookups by the module name are redirected.
GROUP_ALIASES = {"species": "species_"}


class ExportError(RuntimeError):
    """Raised for malformed export files, missing variables or gate violations."""


@dataclass
class Group:
    """One netCDF group: arrays, subgroups and attributes.

    Arrays and subgroups are reachable as attributes (``g.ainv``) or items
    (``g["ainv"]``). ``meta[name]`` holds the netCDF attributes of each array.
    """

    name: str
    path: str
    arrays: dict[str, Any] = field(default_factory=dict)
    groups: dict[str, "Group"] = field(default_factory=dict)
    attrs: dict[str, Any] = field(default_factory=dict)
    meta: dict[str, dict[str, Any]] = field(default_factory=dict)
    absent: frozenset[str] = frozenset()

    def __getattr__(self, name: str):
        # Only called when normal attribute lookup fails.
        if name.startswith("__"):
            raise AttributeError(name)
        try:
            return self[name]
        except KeyError as exc:
            raise AttributeError(str(exc)) from None

    def _group(self, name: str) -> "Group | None":
        if name in self.groups:
            return self.groups[name]
        alias = GROUP_ALIASES.get(name)
        if alias is not None and alias in self.groups:
            return self.groups[alias]
        return None

    def __getitem__(self, name: str):
        stripped = name.strip("/")
        if "/" in stripped:
            head, rest = stripped.split("/", 1)
            group = self._group(head)
            if group is None:
                raise KeyError(f"{self.path}/{head} is not a group of the export file; "
                               f"groups: {sorted(self.groups)}")
            return group[rest]
        name = stripped
        # Groups win over same-named arrays (e.g. the root coordinate variable
        # 'species' vs the group '/species', which the file stores as '/species_').
        group = self._group(name)
        if group is not None:
            return group
        if name in self.arrays:
            return self.arrays[name]
        if name in self.absent:
            raise KeyError(f"{self.path}/{name} was not allocated in GS2 (listed as absent)")
        raise KeyError(
            f"{self.path}/{name} is not in the export file and not listed as absent; "
            f"available: {sorted(self.arrays) + sorted(self.groups)}")

    def __contains__(self, name: str) -> bool:
        try:
            self[name]
        except KeyError:
            return False
        return True

    def has(self, name: str) -> bool:
        """True if the array/group exists (allocated in GS2)."""
        return name in self

    def is_absent(self, name: str) -> bool:
        return name in self.absent

    def get(self, name: str, default=None):
        try:
            return self[name]
        except KeyError:
            return default

    def walk(self) -> Iterator[tuple[str, Any]]:
        """Yield (path, array) for every array in this group and below."""
        for name, value in self.arrays.items():
            yield f"{self.path}/{name}".replace("//", "/"), value
        for group in self.groups.values():
            yield from group.walk()

    def __repr__(self) -> str:  # keep reprs short; groups can hold GBs
        return (f"Group({self.path!r}, arrays={sorted(self.arrays)}, "
                f"groups={sorted(self.groups)})")


def _attrs(obj) -> dict[str, Any]:
    out = {}
    for key in obj.ncattrs():
        value = obj.getncattr(key)
        if isinstance(value, np.ndarray) and value.size == 1:
            value = value.reshape(()).item()
        elif isinstance(value, np.generic):
            value = value.item()
        out[key] = value
    return out


def convert(var: netCDF4.Variable, raw: np.ndarray, attrs: dict[str, Any]):
    """Apply the 09 section 6.6 type rules to one variable's raw data."""
    ftype = str(attrs.get("fortran_type", "")).lower()
    dims = var.dimensions
    if not ftype:
        # Fallback heuristics for variables written without attributes.
        if dims and dims[0] == "ri" and raw.shape[-1:] == (2,):
            ftype = "complex"
        elif raw.dtype == np.int8:
            ftype = "logical"
    if ftype == "complex":
        if raw.shape[-1:] != (2,):
            raise ExportError(f"{var.group().path}/{var.name}: complex without trailing ri axis, "
                              f"shape {raw.shape}")
        data = raw[..., 0].astype(np.float64) + 1j * raw[..., 1].astype(np.float64)
        data = data.astype(np.complex128)
    elif ftype == "logical":
        data = raw != 0
    elif raw.dtype.kind == "f":
        data = raw.astype(np.float64, copy=False)
    elif raw.dtype.kind in "iu":
        data = raw.astype(np.int32, copy=False) if raw.dtype != np.int32 else raw
    else:
        data = raw
    if isinstance(data, np.ndarray) and data.ndim == 0:
        return data.item()
    return data


def read_group(nc, path: str) -> Group:
    """Read one netCDF group (and its subgroups) eagerly."""
    attrs = _attrs(nc)
    absent_attr = attrs.get("absent", "")
    absent = frozenset(str(absent_attr).split()) if absent_attr else frozenset()
    group = Group(name=path.rsplit("/", 1)[-1] or "/", path=path or "/", attrs=attrs,
                  absent=absent)
    for name, var in nc.variables.items():
        vattrs = _attrs(var)
        group.meta[name] = vattrs | {"dimensions": var.dimensions}
        group.arrays[name] = convert(var, np.asarray(var[...]), vattrs)
    for name, sub in nc.groups.items():
        group.groups[name] = read_group(sub, f"{path}/{name}")
    return group


@dataclass
class Export:
    """A loaded export file: ``root`` group, root dimensions and global attributes."""

    path: Path
    root: Group
    dims: dict[str, int]
    attrs: dict[str, Any]

    def __getattr__(self, name: str):
        if name.startswith("__"):
            raise AttributeError(name)
        return getattr(self.root, name)

    def __getitem__(self, name: str):
        return self.root[name]

    def __contains__(self, name: str) -> bool:
        return name in self.root

    @property
    def gate_status(self) -> str:
        return str(self.attrs.get("gate_status", "ok"))

    # Convenience sizes used everywhere.
    @property
    def nth(self) -> int:
        return self.dims["theta"]

    @property
    def ntgrid(self) -> int:
        return (self.dims["theta"] - 1) // 2

    @property
    def ntheta0(self) -> int:
        return self.dims["kx"]

    @property
    def naky(self) -> int:
        return self.dims["ky"]

    @property
    def nlambda(self) -> int:
        return self.dims["lambda"]

    @property
    def negrid(self) -> int:
        return self.dims["egrid"]

    @property
    def nspec(self) -> int:
        return self.dims["species"]

    @property
    def g_shape(self) -> tuple[int, ...]:
        """NumPy shape of a canonical G7 array, (S, E, L, Y, X, 2, T)."""
        return (self.nspec, self.negrid, self.nlambda, self.naky, self.ntheta0, 2, self.nth)

    @property
    def field_shape(self) -> tuple[int, int, int]:
        """NumPy shape of a canonical F3 field array, (Y, X, T)."""
        return (self.naky, self.ntheta0, self.nth)


def _check_header(path, attrs, allow_gates: bool):
    version = str(attrs.get("export_format_version", "0"))
    major = int(version.split(".")[0]) if version.split(".")[0].isdigit() else 0
    if major != SUPPORTED_FORMAT_MAJOR:
        raise ExportError(f"{path}: export_format_version {version!r}; "
                          f"this loader supports major version {SUPPORTED_FORMAT_MAJOR}")
    status = str(attrs.get("gate_status", "ok")).strip() or "ok"
    if status != "ok" and not allow_gates:
        raise ExportError(
            f"{path}: the input uses options the JAX port does not support "
            f"(gate_status): {status}. Pass allow_gates=True only for debugging.")


def load_export(path: str | Path, *, allow_gates: bool = False) -> Export:
    """Load a ``.gs2export.nc`` file eagerly into NumPy arrays."""
    path = Path(path)
    with netCDF4.Dataset(path) as ds:
        ds.set_auto_mask(False)
        attrs = _attrs(ds)
        _check_header(path, attrs, allow_gates)
        dims = {name: len(dim) for name, dim in ds.dimensions.items()}
        root = read_group(ds, "")
    root.path = "/"
    return Export(path=path, root=root, dims=dims, attrs=attrs)
