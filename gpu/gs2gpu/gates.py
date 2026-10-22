"""Refuse inputs that use GS2 options the JAX port does not implement.

The exporter evaluates every gate of spec 10 section B and 00_decisions D7
on the Fortran side and writes one int8 flag per gate to ``/gates``
(1 = violated, with a ``description`` attribute) and the list of violated
gates to the global attribute ``gate_status`` (spec 09a section 5.4).
``check_gates`` turns any violation into one error that lists them all, so a
run can never silently produce a result that is not GS2's.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .io.export import Export


class GateError(RuntimeError):
    """The input uses options the JAX port does not support."""


@dataclass(frozen=True)
class Violation:
    name: str
    description: str


def violations(ex: Export) -> list[Violation]:
    """Every violated gate of an export (from ``/gates`` and ``gate_status``)."""
    found: dict[str, Violation] = {}
    if "gates" in ex.root.groups:
        grp = ex["/gates"]
        for name, value in grp.arrays.items():
            if bool(np.asarray(value).item()):
                desc = str(grp.meta.get(name, {}).get("description", ""))
                found[name] = Violation(name, desc)
    status = str(ex.attrs.get("gate_status", "ok")).strip()
    if status and status != "ok":
        for item in status.split(";"):
            item = item.strip()
            if item and item not in found:
                found[item] = Violation(item, "listed in gate_status")
    return sorted(found.values(), key=lambda v: v.name)


def check_gates(ex: Export) -> None:
    """Raise ``GateError`` if the export violates any gate."""
    bad = violations(ex)
    if bad:
        lines = [f"  - {v.name}: {v.description}" if v.description else f"  - {v.name}"
                 for v in bad]
        raise GateError(
            f"{ex.path}: this input uses GS2 options the JAX port does not implement, "
            "so it cannot reproduce GS2's numbers:\n" + "\n".join(lines))
