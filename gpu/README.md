# gs2gpu: GS2's time step on JAX (GPU)

A faithful port of the GS2 8.2.1 time step to [JAX](https://github.com/jax-ml/jax) in double
precision. Fortran GS2 stays the front end: a small Fortran program, `bin/gs2_export`, runs GS2's
own initialisation (input parsing, geometry, grids, coefficients, response matrices, initial
condition) and writes everything the time step reads to a netCDF file. `gs2gpu` reads that file
and advances the run on a GPU (or CPU) with GS2's numerics: same discretisation, same implicit
scheme, same response-matrix field solve, same collision operator, same time-step control, and
the same diagnostics written to a GS2-format `<run>.out.nc`.

The goal is replication, not approximation: on every validation case the JAX run agrees with GS2
to within GS2's own round-off spread (see [Validation](#validation)).

## Workflow

```bash
# 1. Build GS2 and the exporter (once). Default = reference build (-ffp-contract=off);
#    --stock uses GS2's normal flags. Either works as a front end.
gpu/scripts/build_reference.sh

# 2. Export the case (CPU, one MPI rank): writes case.gs2export.nc and case.out.nc
#    (GS2's static output and the t = 0 record) next to the input.
mpirun -np 1 bin/gs2_export case.in

# 3. Run it with JAX; the output is a GS2-format out.nc.
gpu/.venv/bin/gs2gpu run case.gs2export.nc -o case.jax.out.nc --progress 100 --checkpoint-every 1000

# Resume after an interruption (same step counter, same output file):
gpu/.venv/bin/gs2gpu run case.gs2export.nc -o case.jax.out.nc --restart case.jax.out.nc.ckpt

# 4. Optional: compare with a GS2 run of the same input, variable by variable.
gpu/.venv/bin/gs2gpu compare case.out.nc case.jax.out.nc
```

The export file records every option the input uses. Inputs that use an option the port does not
implement are refused with a list of the offending options (the *gates*), so that a run does not
silently produce something that is not GS2's answer.

## GPU setup

The port needs a GPU with full-rate FP64: NVIDIA A100/H100/H200 or AMD MI250/MI300. Consumer cards
(RTX 4090/5090) and FP64-reduced data-centre parts (L40S, RTX PRO 6000, B300) run FP64 at 1/32 to
1/64 rate. In the project venv:

```bash
gpu/.venv/bin/pip install -U "jax[cuda12]"        # NVIDIA; see the JAX docs for ROCm
gpu/.venv/bin/python -c "import jax; print(jax.devices())"
```

`import gs2gpu` enables float64 (`jax_enable_x64`). Nothing else changes between CPU and GPU.

## What is supported

Linear and nonlinear flux-tube runs; electrostatic and electromagnetic (A∥, B∥); adiabatic,
kinetic and hybrid electrons; `single`, `range` and `box` k-grids; s-alpha, Miller and numerical
geometry (whatever GS2's front end builds); every boundary option (`linked`, `zero`,
`self-periodic`, …) and wfb treatment; `field_option = implicit` and `local`; collisions
(`default`, `lorentz`, `ediffuse`, with conservation, drag and `split_collisions`); hyperviscosity
(`visc_only`, `res_only`, `both`, `simple3D`, iso/anisotropic); ExB flow shear (both branches);
zonal-flow removal; parity; CFL-controlled time-step changes on both reset paths; the standard
diagnostics (`phi2`, `omega`, fluxes, moments, eigenfunctions, …).

Gated (refused): the antenna, `opt_source`, `split_nonlinear`, `zip`, `eqzip`, `wstar_units`,
`use_old_diagnostics`, `use_le_layout = F`, experimental boundary options, hybrid electrons *with*
collisions, and diagnostics that feed back into the run (`write_verr`, heating, …). The full table
with reasons is `spec/10_critique.md` §B (plus `spec/00_decisions.md` D7).

## Validation

`gpu/cases/` holds 31 small cases that together exercise every supported option family and the
run-control exit paths (`cases/README.md`). For each case `gpu/scripts/make_refdata.sh` produces GS2's full run, the
export, and a *trace*: GS2's state dumped between every stage of the first 10 steps.

* **Kernels.** Every stage of the step (source, parallel solve, linked boundaries, parity, field
  solve, nonlinear term, collisions, hyperviscosity, flow shear, zonal removal, diagnostics) is
  compared with GS2's dump of the same stage, every traced step and pass, every case. Worst
  errors: source 1.0e-15, parallel solve 1.05e-15, field solve 1.0e-14, nonlinear term 1.3e-15
  (all of it XLA's FFT versus FFTW), collisions 1.9e-15, hyperviscosity / flow shear / zonal
  removal 1.5e-16 (the shifts bitwise), diagnostics 5.6e-15. With XLA's FMA contraction disabled
  ("strict" mode) the collisions, field solve, source and diagnostics are bitwise identical to
  GS2.
* **One step.** One JAX step from GS2's state reproduces GS2's next state to <= 1e-12 in every
  case (`tests/test_step.py`), including steps with immediate and end-of-step time-step changes
  and the response-matrix rebuild at a new time step.
* **Whole runs.** Every case is run to the end with the JAX driver and its `out.nc` compared with
  GS2's, variable by variable (`tests/test_end_to_end.py`, 225–304 variables per case). The worst
  relative error is 3e-15–5e-13 in most cases and up to 4e-10 in a few nonlinear and flow-shear
  cases; in every case it is the same size as GS2's own spread when only its summation order
  changes (GS2 rerun with another data layout, `scripts/make_noise_refdata.sh`: e.g. 3.7e-10 for
  JAX against 5.2e-10 for GS2 against itself in the anisotropic-hyperviscosity case). Quantities
  that are zero up to round-off (e.g. the particle flux with adiabatic electrons, ~1e-24) agree
  to within 0.5–6x GS2's own round-off. Growth rates and frequencies agree to all printed digits,
  and time-step changes happen at the same steps with the same values.
* **Restart.** Stopping and resuming from a checkpoint gives a bitwise identical `out.nc`.

The suite is 1742 tests: 1590 pass and 152 skip because they do not apply to a case (e.g. parity checks when `def_parity` is off); none fail.

Results are not bitwise identical to GS2: XLA's FFTs differ from FFTW's in the last bit, XLA
contracts multiply-adds, and GS2's own compiler fuses vectorised complex products (`fcmla`) even
with `-ffp-contract=off`. These are 1-ulp effects per operation; for nonlinear runs they grow
like any round-off difference in a chaotic system, exactly as between two GS2 runs on different
core counts, so nonlinear results are compared statistically beyond the first steps.

Run the tests:

```bash
gpu/scripts/build_reference.sh && gpu/scripts/make_refdata.sh && gpu/scripts/make_noise_refdata.sh
cd gpu && .venv/bin/python -m pytest -q          # ~30 min on an M1 Max CPU
```

## Performance

`gs2gpu bench case.gs2export.nc` times steady-state steps on the current device. On an M1 Max CPU
(JAX, multi-threaded) against GS2 on one core (reference build), for the small test grids:

| case | JAX CPU | GS2 1 core |
|---|---|---|
| `lin_single_ae` | 0.20 ms/step | 0.40 ms/step |
| `lin_single_coll_ke_full` | 0.78 ms/step | 0.96 ms/step |
| `nl_box_ae` (16x16) | 5.9 ms/step | 12.4 ms/step |
| `nl_box_em_ke_coll` | 6.6 ms/step | 8.6 ms/step |

These grids are far too small to keep a GPU busy; GPU throughput has not been measured yet.

## Known limitations

* **Single device.** No multi-GPU sharding yet; a case must fit in one GPU's memory (several
  copies of the distribution function, e.g. ~5–10 GB each for a large kinetic-electron run).
* **Front end on one MPI rank.** `gs2_export` runs GS2's initialisation, including GS2's own
  response-matrix build, on one CPU core. For large cases that can take a while (it is done once
  per run; restarts reuse the export).
* **Not yet tuned for GPUs.** The step is jitted per stage, with one host round trip per step
  (the CFL check) and more on diagnostic write steps. Memory chunking for very large grids and
  batching of the response-matrix rebuild after a time-step change are not done yet. Performance
  has only been measured on CPU.
* **Diagnostics.** No ASCII outputs (`.out` text lines). `write_final_epar`,
  `write_final_moments`, `write_final_antot` and the `*_over_time` outputs are implemented but
  validated only on synthetic data.
* **Round-off-level GS2 quirk not replicated.** After its first time-step change GS2 re-reads
  `collisions_knobs`, which turns the `heating` flag back on and adds one more `g_adjust` round
  trip on `g` per step; this changes g-dependent write-step diagnostics in the last bit only.

## Code layout

| module | GS2 counterpart |
|---|---|
| `io/export.py`, `io/trace.py` | readers for the export and trace files |
| `types.py`, `numerics.py`, `gadjust.py`, `gates.py` | shared containers, gfortran-faithful complex arithmetic, `g_adjust`, option gates |
| `coeffs.py` | `init_vpar`, `init_wdrift`, `init_wstar`, `init_invert_rhs` |
| `source.py` | `get_source_term`, Adams–Bashforth coefficients |
| `invert_rhs.py`, `linked.py`, `parity.py` | `invert_rhs_1`, linked boundary conditions, `enforce_parity`, `init_homogeneous_g` |
| `fields.py`, `fields_local.py`, `response.py` | field equations, implicit/local field solve, response-matrix construction |
| `nonlinear.py`, `explicit.py` | nonlinear term, explicit history, CFL control |
| `collisions.py`, `collisions_init.py` | `vspace_derivatives` (le layout), collision matrices |
| `hyper.py`, `exb.py`, `zonal.py` | hyperviscosity, ExB shear, zonal-flow removal |
| `diagnostics.py`, `outnc.py` | `run_diagnostics`, the out.nc writer |
| `step.py`, `reset.py`, `driver.py`, `checkpoint.py`, `cli.py` | `advance`, `timeadv`, `reset_time_step`, `evolve_equations`, restart, command line |

`spec/` holds the numerical specifications every module implements (`00_decisions.md` first,
then `10_critique.md`), and `spec/09a_exporter_notes.md` documents the Fortran exporter.
