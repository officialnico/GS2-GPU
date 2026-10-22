# 00 — Project decisions (binding)

These resolve the open questions in `10_critique.md` §F and fix the project structure. Where this file
and any other spec disagree, this file wins; after it, `10_critique.md` wins over specs 01–09.

## D1. Goal and tolerance policy

- Faithful replication of GS2 8.2.1 numerics in JAX (float64 / complex128), single device.
- No bitwise chasing. Plain `jnp` arithmetic, but **keep the Fortran association order** of every
  expression (the specs give it), and the Fortran summation order for reductions where a spec says
  it matters (velocity integrals, matvecs). Use the specs' `cmul`/`rmul` helpers for mixed
  real/complex products.
- Do **not** wrap everything in `optimization_barrier`; that kills GPU performance. It may be used in a
  test-only "strict" mode if a kernel needs it to prove an exactness claim.
- Test tolerances (relative to the max-abs of the reference array, unless stated):
  - single kernel on GS2-dumped inputs: `rtol <= 1e-12` (target 1e-14);
  - one full time step: `<= 1e-12`;
  - linear runs over the whole run: growth rate and frequency agree to `<= 1e-8` relative,
    eigenfunctions to `<= 1e-8`;
  - nonlinear runs: per-step state agreement `<= 1e-10` for the first 20 steps (before chaos
    amplifies round-off); after that, time-averaged fluxes within statistical scatter.
- A test that needs a looser tolerance than the above is a bug report, not a tolerance change.

## D2. Reference build

- The reference binaries `bin/gs2` and `bin/gs2_export` are built by `gpu/scripts/build_reference.sh`
  with `-ffp-contract=off` appended to the Fortran flags (removes FMA contraction, so differences
  from JAX are pure association/rounding), single MPI rank, no OpenMP. The exact command is recorded
  in that script (this resolves Q-01 and Q-11).
- The script must also be able to build the stock flags so we can confirm that the public-edit
  changes to the Fortran do not change GS2's output (09 §4.4).

## D3. dt changes (Q-02)

- RECOMPUTE in JAX is mandatory for every dt-dependent array (critique §E). No ladder cache in v1.
- The exporter supports dumping coefficients at extra dt values (`dt_ladder_n`) **only as test
  oracles** for the recompute path.

## D4. Bug replication (Q-03)

- Replicate every GS2 behaviour that is not gated, including the quirks listed in Q-03 (theta0-branch
  all-ky zeroing, `damp_zonal_only`, frozen `D_res`/simple3D filter). Each replicated quirk gets a
  comment `# GS2 quirk: <description> (file:line)`.

## D5. field_option='local' (Q-04)

- Supported at the initial dt with the exported `fieldmat` and the 04/G-23 matvec order.
- On a dt change in local mode, rebuild with the JAX response builder using GS2's own local probe
  grouping (`fm_populate`: one impulse in every unfinished supercell of every ky per
  `timeadv(...,0)` call), which matters with amplitude-dependent hyperviscosity. The inversion is
  Gauss–Jordan rather than a port of `invert_block`: a round-off-level deviation (measured
  <= 7.2e-16 on the dt-ladder rungs of `nl_box_ae_local_hyper_dt`); it is not a gate.

## D6. Diagnostics (Q-05)

- Reference runs keep default diagnostics. JAX replicates the `g_adjust` round trips on write steps.

## D7. Gates added by this file

- Q-06: the exporter asserts that bmax is at `ig = ±ntgrid` when trapped particles exist; else GATE.
- Q-07: GATE `shat == 0` exactly when any `aky == 0` mode exists.
- Q-08: replicate (freeze fields of non-box `aky = akx = 0` modes).
- Q-09: the exporter asserts no `1 - al*bmag == -bouncefuzz` point; else GATE.
- Q-10: JAX writes its own checkpoints that include every SAVE counter; mid-run GS2 trace restarts are
  out of scope.
- Q-12: stays GATED.

## D8. Package layout and file ownership

All Python lives in `gpu/`. Each implementation task owns its files exclusively; do not edit files
owned by another task (if you need a change there, record it in your result instead).

| owner | files |
|---|---|
| A (foundation) | `gpu/pyproject.toml`, `gpu/gs2gpu/__init__.py`, `gpu/gs2gpu/config.py`, `gpu/gs2gpu/numerics.py`, `gpu/gs2gpu/types.py`, `gpu/gs2gpu/gadjust.py`, `gpu/gs2gpu/io/*`, `gpu/gs2gpu/gates.py`, `gpu/tests/conftest.py`, `gpu/tests/test_io.py`, `gpu/cases/*`, `gpu/scripts/*`, Fortran exporter files under `src/`, Makefile edits, `gpu/spec/11_api.md` |
| B1 coeffs | `gpu/gs2gpu/coeffs.py`, `gpu/tests/test_coeffs.py` |
| B2 source | `gpu/gs2gpu/source.py`, `gpu/tests/test_source.py` |
| B3 parallel solve | `gpu/gs2gpu/invert_rhs.py`, `gpu/gs2gpu/linked.py`, `gpu/gs2gpu/parity.py`, `gpu/tests/test_invert_rhs.py` |
| B4 fields | `gpu/gs2gpu/fields.py`, `gpu/gs2gpu/fields_local.py`, `gpu/gs2gpu/response.py`, `gpu/tests/test_fields.py` |
| B5 nonlinear | `gpu/gs2gpu/nonlinear.py`, `gpu/gs2gpu/explicit.py`, `gpu/tests/test_nonlinear.py` |
| B6 collisions | `gpu/gs2gpu/collisions.py`, `gpu/gs2gpu/collisions_init.py`, `gpu/tests/test_collisions.py` |
| B7 misc | `gpu/gs2gpu/hyper.py`, `gpu/gs2gpu/exb.py`, `gpu/gs2gpu/zonal.py`, `gpu/tests/test_misc.py` |
| B8 diagnostics | `gpu/gs2gpu/diagnostics.py`, `gpu/gs2gpu/outnc.py`, `gpu/tests/test_diagnostics.py` |
| C integration | `gpu/gs2gpu/timeadv.py`, `gpu/gs2gpu/step.py`, `gpu/gs2gpu/reset.py`, `gpu/gs2gpu/driver.py`, `gpu/gs2gpu/checkpoint.py`, `gpu/gs2gpu/cli.py`, `gpu/tests/test_step.py`, `gpu/tests/test_end_to_end.py`, `gpu/README.md` |

## D9. Code style

- JAX-first: kernels are pure functions of explicit arguments, `jax.jit`-compatible, no Python loops
  over grid points (Python loops over a handful of static options/classes are fine), no host syncs
  inside a step. Use `vmap`/broadcasting for independent lines and `lax.scan` for theta recurrences.
- Host-side scalar control logic (dt resets, AB coefficients, exb jump bookkeeping, counters) stays
  in plain Python/NumPy float64, mirroring the Fortran literally.
- Every function's docstring cites the Fortran routine and file:line it mirrors.
- Tests use pytest and read reference data from `gpu/refdata/` (generated, gitignored). Tests skip,
  not fail, when their reference case is missing, but `scripts/make_refdata.sh` must generate all.
- Python: `gpu/.venv/bin/python` (JAX 0.12, CPU). Enable x64 via `gs2gpu.config`.
