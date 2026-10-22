"""Reader for ``<run>.gs2trace.nc`` (spec 09 sections 5.3 and 6.2).

The trace holds one group per recorded probe::

    /s<istep:06d>/a<attempt>/p<seq:03d>_<tag>   attrs: tag, istep, attempt, seq, pass, in_advance
    /final                                       snapshot after evolve_equations returned

``Trace`` indexes every record on open and reads a record's arrays only when
asked (a trace can hold many full distribution functions). Records are
ordered by (istep, attempt, seq), i.e. in execution order.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

import netCDF4

from .export import ExportError, Group, read_group

_STEP = re.compile(r"^s(\d+)$")
_ATTEMPT = re.compile(r"^a(\d+)$")
_PROBE = re.compile(r"^p(\d+)_(.+)$")


@dataclass(frozen=True)
class Record:
    """Index entry of one probe record."""

    istep: int
    attempt: int
    seq: int
    tag: str
    pass_: int
    in_advance: bool
    path: str

    @property
    def key(self) -> tuple[int, int, int]:
        return (self.istep, self.attempt, self.seq)


class Trace:
    """Random access to the probe records of a trace file."""

    def __init__(self, path: str | Path):
        self.path = Path(path)
        self._ds = netCDF4.Dataset(self.path)
        self._ds.set_auto_mask(False)
        self.attrs = {k: self._ds.getncattr(k) for k in self._ds.ncattrs()}
        self.records: list[Record] = []
        for sname, sgrp in self._ds.groups.items():
            ms = _STEP.match(sname)
            if not ms:
                continue
            for aname, agrp in sgrp.groups.items():
                ma = _ATTEMPT.match(aname)
                if not ma:
                    continue
                for pname, pgrp in agrp.groups.items():
                    mp = _PROBE.match(pname)
                    if not mp:
                        continue
                    attrs = {k: pgrp.getncattr(k) for k in pgrp.ncattrs()}
                    self.records.append(Record(
                        istep=int(attrs.get("istep", ms.group(1))),
                        attempt=int(attrs.get("attempt", ma.group(1))),
                        seq=int(attrs.get("seq", mp.group(1))),
                        tag=str(attrs.get("tag", mp.group(2))),
                        pass_=int(attrs.get("pass", 0)),
                        in_advance=bool(int(attrs.get("in_advance", 0))),
                        path=f"/{sname}/{aname}/{pname}",
                    ))
        self.records.sort(key=lambda r: r.key)
        self._load = lru_cache(maxsize=16)(self._load_uncached)

    # ------------------------------------------------------------------ index
    def close(self):
        self._ds.close()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()

    @property
    def steps(self) -> list[int]:
        return sorted({r.istep for r in self.records})

    def attempts(self, istep: int) -> list[int]:
        return sorted({r.attempt for r in self.records if r.istep == istep})

    def last_attempt(self, istep: int) -> int:
        attempts = self.attempts(istep)
        if not attempts:
            raise ExportError(f"{self.path}: no records for istep {istep}")
        return attempts[-1]

    def find(self, tag: str | None = None, *, istep: int | None = None,
             attempt: int | None = None, pass_: int | None = None) -> list[Record]:
        """All records matching the given fields (None = any), in execution order."""
        out = []
        for r in self.records:
            if tag is not None and r.tag != tag:
                continue
            if istep is not None and r.istep != istep:
                continue
            if attempt is not None and r.attempt != attempt:
                continue
            if pass_ is not None and r.pass_ != pass_:
                continue
            out.append(r)
        return out

    def record(self, tag: str, istep: int, *, attempt: int | None = None,
               pass_: int | None = None, occurrence: int = 0) -> Record:
        """The ``occurrence``-th record of ``tag`` in (istep, attempt[, pass]).

        ``attempt`` defaults to the last attempt of the step (the one that was
        kept after an immediate reset).
        """
        if attempt is None:
            attempt = self.last_attempt(istep)
        found = self.find(tag, istep=istep, attempt=attempt, pass_=pass_)
        if len(found) <= occurrence:
            raise ExportError(
                f"{self.path}: no record #{occurrence} of tag {tag!r} at istep={istep} "
                f"attempt={attempt} pass={pass_}; tags there: "
                f"{[r.tag for r in self.find(istep=istep, attempt=attempt)]}")
        return found[occurrence]

    def previous(self, rec: Record, tag: str | None = None) -> Record | None:
        """The closest earlier record of the same (istep, attempt), optionally of ``tag``."""
        same = [r for r in self.records
                if r.istep == rec.istep and r.attempt == rec.attempt and r.seq < rec.seq
                and (tag is None or r.tag == tag)]
        return same[-1] if same else None

    def next(self, rec: Record, tag: str | None = None) -> Record | None:
        """The closest later record of the same (istep, attempt), optionally of ``tag``."""
        same = [r for r in self.records
                if r.istep == rec.istep and r.attempt == rec.attempt and r.seq > rec.seq
                and (tag is None or r.tag == tag)]
        return same[0] if same else None

    # ------------------------------------------------------------------- data
    def _load_uncached(self, path: str) -> Group:
        return read_group(self._ds[path], path)

    def load(self, rec: Record | str) -> Group:
        """Arrays of one record (or of ``/final`` when given that path)."""
        path = rec.path if isinstance(rec, Record) else rec
        return self._load(path)

    def get(self, tag: str, istep: int, **kwargs) -> Group:
        """Shortcut: ``load(record(tag, istep, ...))``."""
        return self.load(self.record(tag, istep, **kwargs))

    @property
    def final(self) -> Group:
        return self.load("/final")
