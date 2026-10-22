# 11 — Python API contract for the kernel modules

Binding for every module in `gpu/gs2gpu/`. Read `00_decisions.md` (tolerances, ownership) and
`10_critique.md` first; your owner spec(s) define the numerics, this file defines the code shape.

## 1. Environment

* Python: `/Users/nico/Workspace/gs2/gpu/.venv/bin/python` (JAX 0.12 CPU, numpy, scipy, netCDF4,
  pytest; package `gs2gpu` installed editable). Run tests with
  `cd /Users/nico/Workspace/gs2/gpu && .venv/bin/python -m pytest -q tests/test_<module>.py`.
* `import gs2gpu` enables float64 (`jax_enable_x64`). Every array is float64 / complex128.
* Reference data: `gpu/refdata/<case>/<case>.gs2export.nc`, `<case>.gs2trace.nc`, `<case>.out.nc`,
  produced by `gpu/scripts/make_refdata.sh` from `gpu/cases/*.in` (see `gpu/cases/README.md` for what
  each case covers). The exporter's implementation notes, including the exact variable and probe names it
  writes, are in `gpu/spec/09a_exporter_notes.md` once it exists; where they differ from spec 09, the notes win.

## 2. Reading reference data (owner A, already written)

```python
from gs2gpu.io.export import load_export          # eager NumPy load of an export file
from gs2gpu.io.trace import Trace                 # lazy access to probe records
ex = load_export(path)                            # Export
ex["/dist_fn/ainv"]                               # complex128, canonical NumPy shape (G7)
ex["/theta_grid"].bmag                            # attribute access works too
"gamtot3" in ex["/dist_fn"]                       # False if not allocated
ex["/dist_fn"].is_absent("gamtot3")               # True if GS2 did not allocate it
ex.dims, ex.ntgrid, ex.nth, ex.g_shape, ex.field_shape
ex["/dist_fn"].meta["ainv"]                       # netCDF attributes (fortran_name, source, dt_dependent...)

tr = Trace(path)
rec = tr.record("source", istep=1, pass_=1)       # attempt defaults to the last attempt of the step
grp = tr.load(rec)                                # Group with the arrays dumped at that probe
grp["g_work"]; grp["gnew"]; grp["phinew"]         # names as in spec 09 §5.3 (Fortran names)
tr.previous(rec, "advance_begin"); tr.next(rec, "after_invert_rhs_1")
tr.find("after_timeadv", istep=1)                 # both passes, in execution order
```

Conversions (spec 09 §6.6): complex → complex128 (the trailing `ri` axis is folded), logical → bool,
integer → int32 with Fortran values (1-based indices stay 1-based, `-1` stays `-1`). Axes are never
reordered: NumPy shape is the Fortran shape reversed.

## 3. Canonical layouts (spec 09 §0.2, spec 10 §0.1)

| tag | NumPy shape | Fortran |
|---|---|---|
| G7 | `(S, E, L, Y, X, 2, T)` | `g(ig, isgn, it, ik, il, ie, is)` |
| G6 | `(S, E, L, Y, X, T)` | `x(ig, it, ik, il, ie, is)` |
| F3 | `(Y, X, T)` | `phi(ig, it, ik)` |
| K2 | `(Y, X)` | `x(it, ik)` |
| `wl`, `forbid`, `is_ttp`, … | `(L, T)` | `x(ig, il)` |
| `wstar` | `(S, E, Y)` | `wstar(ik, ie, is)` |
| LE6 | `(S, Y, X, T, NE, NX)` | `x(ixi, ie, ile)` (spec 06 §2.4) |

S = nspec, E = negrid, L = nlambda, Y = naky, X = ntheta0, T = nth = 2·ntgrid+1. Theta index `t`
↔ Fortran `ig = t - ntgrid` (`gs2gpu.numerics.ith`). Sign index 0/1 ↔ Fortran `isgn` 1/2. Index-valued
integers read from the export keep Fortran values: subtract 1 before using them as NumPy indices.
Kernels keep these layouts at their boundaries; inside a kernel you may transpose (e.g. move theta to
axis 0 for `lax.scan`) but must transpose back.

## 4. Shared modules (owner A, already written — import, do not edit)

* `gs2gpu.numerics`: `rmul, imul, cmul, cdiv, recip, radd, rsub` (mirror gfortran complex lowering,
  spec 02 §3.0 / 03 §3.0), `set_strict` (FMA-free test mode), `ith`, `max_abs_rel`.
* `gs2gpu.types`: `Dims` (frozen dataclass of sizes/switches), and NamedTuple containers built with
  `from_export(ex)`: `Geometry`, `KGrids`, `VGrids`, `Species`, `PointArrays` (vpa, vpac, vperp2, aj0,
  aj1), `Coeffs` (dt-dependent: code_dt, vpar, wdrift, wdriftttp, wstar, a, b, r, ainv, g_h;
  `Coeffs.from_export(ex, rung=k)` reads a dt-ladder rung), `Fields(phi, apar, bpar)`, `TimeBlock`;
  `fields_from_group(grp, dims)`; `to_jax(tree)`.
* `gs2gpu.gadjust`: `GAdjustData.from_export(ex)`, `g_adjust(g, phi, bpar, facphi, facbpar, data)`,
  `g_adjust_direction(g, phi, bpar, FROM_G|TO_G, fphi, fbpar, data)`.

If a shared container lacks something you need, build it in your own module from the export rather
than editing the shared file, and say so in your report.

## 5. Shape of a kernel module

```python
class XParams(NamedTuple):            # or a frozen dataclass: Python scalars/bools/strings only (static)
    ...
    @classmethod
    def from_export(cls, ex): ...

class XData(NamedTuple):              # arrays (NumPy from the loader; to_jax() before use)
    ...
    @classmethod
    def from_export(cls, ex): ...

def kernel(state_arrays..., data: XData, params: XParams, ...) -> outputs:
    """Mirrors <routine> (src/<file>:<lines>). ..."""
```

Rules:
1. Kernels are pure functions of explicit arguments and must be `jax.jit`-compatible with `params`
   (and other Python scalars that select code paths) static. No module-level mutable state. No
   Python loops over grid points; Python loops over a handful of static items (species, classes,
   options) are fine. No host synchronisation inside a kernel.
2. Use broadcasting/`vmap` for independent lines and `lax.scan` (or `lax.associative_scan`) for theta
   recurrences. Masks (forbid, ttp, wfb, bounce points, kwork_filter) are boolean arrays applied with
   `jnp.where`. Irregular structures (linked chains, supercells) are padded to rectangular arrays built
   on the host in `from_export`.
3. Host-side scalar logic (dt resets, Adams–Bashforth coefficients, exb jump bookkeeping, counters)
   is plain Python/NumPy float64 that mirrors the Fortran literally.
4. Keep the Fortran association order of every expression and the summation order wherever your spec
   says it matters. Use the `numerics` helpers for mixed real/complex products.
5. Every public function's docstring cites the Fortran routine and `file:line` it mirrors. Replicated
   GS2 quirks get a comment `# GS2 quirk: <what> (file:line)` (00_decisions D4).
6. Gated options: if your module sees an option it does not implement, raise `NotImplementedError`
   naming the option (the loader's gate check should already have refused it).

## 6. Tests

* `gpu/tests/test_<module>.py`, pytest, using `tests/conftest.py`: `export_of(case)`,
  `trace_of(case)`, `out_of(case)`, `cases_where(predicate)` (parametrize over cases that exercise your
  code), `assert_close(ref, got, rtol, what)` (max-abs-relative), `RTOL_KERNEL = 1e-12`.
* The standard kernel test: take the kernel's inputs from the trace record **before** the kernel (or
  from the export for coefficients) and compare the output with the record **after** it, for every
  traced step and both passes, over every case that exercises the kernel. Coefficient recomputation is
  tested against the export at the initial dt and against every dt-ladder rung.
* A test that needs a looser tolerance than 00_decisions D1 is a bug to find, not a tolerance to change.
  If you cannot find it, mark the test `xfail(strict=True)` with the measured error and an explanation,
  and report it.
* Tests must skip (not fail) when reference data is missing.

## 7. Worked example (pattern for a kernel test)

```python
import pytest, jax.numpy as jnp
from conftest import export_of, trace_of, cases_where, assert_close, RTOL_KERNEL
from gs2gpu.types import Dims, to_jax
from gs2gpu import invert_rhs as ir

CASES = cases_where(lambda ex: True)

@pytest.mark.parametrize("case", CASES)
def test_invert_rhs_1(case):
    ex, tr = export_of(case), trace_of(case)
    data, params = to_jax(ir.InvertData.from_export(ex)), ir.InvertParams.from_export(ex)
    for istep in tr.steps:
        for pass_ in (1, 2):
            src = tr.load(tr.record("source", istep, pass_=pass_))["g_work"]
            out = tr.load(tr.record("after_invert_rhs_1", istep, pass_=pass_))["gnew"]
            fields = ...   # phinew/bparnew valid during this pass (see 01 §3.3)
            got = ir.invert_rhs_1(jnp.asarray(src), fields.phi, fields.bpar, data, params)
            assert_close(out, got, RTOL_KERNEL, f"{case} istep={istep} pass={pass_}")
```

## 8. Module interfaces the integration layer (owner C) will call

These signatures are the contract between B-modules and C. Extra keyword arguments are fine; record
any deviation in your report.

| module | must provide |
|---|---|
| `coeffs` (B1) | `CoeffInputs.from_export(ex)`, `CoeffParams.from_export(ex)`, `compute_coeffs(code_dt, inputs, params) -> types.Coeffs` with `g_h=None` (vpar, wdrift, wdriftttp, wstar, a, b, r, ainv per spec 02 §3.1 / 03 §3.3) |
| `source` (B2) | `SourceData/SourceParams.from_export`, `ab_coefficients(time: TimeBlock, istep) -> tuple` (host, spec 01 §3.10 / 02), `get_source_term_on_local_domain(g, fields_old: Fields, fields_new: Fields, gexp: tuple or None, ab: tuple, istep: int, coeffs: Coeffs, data, params, *, code_time, aparnew_module=None) -> g_work` (spec 02 §3.2-3.5; returns the whole G7 source incl. zeros at `ig = ntgrid` and at kwork_filter points) |
| `invert_rhs` + `linked` + `parity` (B3) | `InvertData/InvertParams.from_export`, `invert_rhs_1(source, phinew, bparnew, coeffs, data, params) -> gnew`, `linked.apply_linked_boundary_conditions(gnew, phinew, bparnew, coeffs, data, params) -> gnew`, `parity.enforce_parity(gnew, data, params) -> gnew`, `init_homogeneous_g(coeffs, data, params) -> g_h` (spec 03 §3.4) |
| `fields` + `fields_local` + `response` (B4) | `FieldData/FieldParams.from_export`, `getan(gnew, data, params) -> (antot, antota, antotp)`, `getfieldeq(fields_new: Fields, antot..., data, params) -> (fieldeq, fieldeqa, fieldeqp)`, `getfield(fields_new: Fields, gnew, aminv, data, params) -> Fields` (the increment u), `get_init_field(gnew, data, params) -> Fields` (get_fields_direct_from_dfn), `calculate_potentials_from_nonadiabatic_dfn(h, data, params) -> Fields`, `response.Aminv.from_export(ex)`, `response.build_aminv(timeadv0: Callable[[Fields], gnew], data, params) -> Aminv` (probe construction of spec 04 §5.4-5.5, given a callable that runs `timeadv(…, istep=0)` for unit field perturbations), and the field_option='local' equivalents |
| `nonlinear` + `explicit` (B5) | `NLData/NLParams.from_export`, `add_nl(g, fields: Fields, data, params) -> (nl_term G7, max_vel_components)`, `nl_center(g1, bd, data, params)`, `explicit.add_explicit_terms(...)` host+device bookkeeping of gexp_1..3 with istep_last, CFL and order-based-error checks returning the new code_dt_cfl and whether a reset is required (spec 05 §4.1-4.10, 02 §3.7) |
| `collisions` + `collisions_init` (B6) | `CollData/CollParams.from_export`, `collisions_init.compute_coll_coeffs(code_dt, vnmult, data, params) -> CollCoeffs` (dt-dependent matrices, spec 06 §5), `CollCoeffs.from_export(ex, rung=None)`, `vspace_derivatives(gnew, g, fields_old: Fields, fields_new: Fields, aparnew, diagnostics: int, coll_coeffs, data, params) -> gnew` (spec 06 §4, le layout), `collisions_advance(...)` for split_collisions (06 §4.10) given a callable for `calculate_potentials_from_nonadiabatic_dfn` |
| `hyper`, `exb`, `zonal` (B7) | `hyper_diff(gnew, phinew, code_dt, data, params) -> (gnew, hypervisc_filter)`; `exb.ExbState` (host bookkeeping) and `exb_shear(gnew, fields, exb_state, istep, time: TimeBlock, data, params) -> (gnew, fields, exb_state)`; `zonal.remove_zonal_flows(phinew, data, params)` (spec 07, 01 §3.7-3.8) |
| `diagnostics` + `outnc` (B8) | `DiagData/DiagParams.from_export`, per-step `calculate_omega` state, the write-step routine that computes and appends one record to `<run>.out.nc` (including the g_adjust round trips of spec 08 §3.1.8 that modify `gnew` and `g`), exit conditions (omega convergence), and `finalize` (spec 08) |
