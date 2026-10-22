# 09a — GS2 exporter: implementation notes

Status: implementation record of spec 09 (`09_export_infra.md`) with the corrections of `10_critique.md`
and `00_decisions.md`. Read 09 for the design; this file lists what was actually built, every deviation
from 09 (and why), the exact names written, and how to run it. Where this file and 09 disagree about
what the files *contain*, this file describes the code.

---------------------------------------------------------------------------------------------------

## 1. What was built

| file | content |
|---|---|
| `src/gs2_export_probe.f90` | module `gs2_export_probe` (09 §5.1): `probe_active`, `probe_in_advance`, `probe_istep`, `probe_pass`, `probe_point/g/le/real/c2` and their hook pointers. No GS2 dependencies. |
| `src/gs2_export_writer.fpp` | module `gs2_export_writer`: netCDF-4 file/group/dimension handling (plain `netcdf` API), attributes, `absent` lists, scalar/V-array writers, canonicalisers `ex_g7c, ex_g7r, ex_g6r, ex_g7m1c, ex_glc, ex_g5i, ex_g5l, ex_g5sl, ex_le6r, ex_le6c` (09 §2.4, per-species slabs, `int64` offsets). |
| `src/gs2_export_static.fpp` | module `gs2_export_static`: knobs type, root dimensions, global attributes, every static group, the state groups, the dt-dependent set of a ladder rung, `/diagnostics`, `/gates`. |
| `src/gs2_export_trace.fpp` | module `gs2_export_trace`: probe hooks and the trace file (09 §5.3, §6.2, §7.9). |
| `src/programs/gs2_export.fpp` | program `gs2_export` (09 §1.4, §7.8, §8.1, §8.3). |
| public edits | `public ::` blocks of 09 §4.2 in `dist_fn.fpp`, `collisions.fpp`, `hyper.f90`, `nonlinear_terms.fpp`, `fields_implicit.f90`, `theta_grid.f90`, `le_grids.f90`, `run_parameters.f90`, plus `gs2_time.f90` (`user_dt_cfl`) and `init_g.f90` (`ginitopt_switch`) — see §3. |
| probe calls | one-line `if (probe_active) ...` calls (and a `use gs2_export_probe, only: ...`) in `fields.f90`, `fields_implicit.f90`, `fields_local.f90`, `dist_fn.fpp`, `nonlinear_terms.fpp`, `le_derivatives.fpp`, `collisions.fpp`, `gs2_main.fpp` (§4). All edits to existing sources are pure insertions (`git diff` shows no deleted or modified line). |
| `Makefile.target_gs2` | 09 §9.2: `gs2_export` in `needs_git_hash`, new `gpu_port_targets` (not part of `gs2_all`). |
| `gpu/scripts/build_reference.sh` | builds `bin/gs2`, `bin/gs2_export` (reference or `--stock` flags), records the mode, forces a full rebuild on a mode change. |
| `gpu/scripts/make_refdata.sh` | reference data for every case (§7). |
| `gpu/scripts/check_export.py` | schema, numerical spot checks, trace invariants, determinism, trajectory checks (§8). |
| `gpu/scripts/compare_netcdf.py` | bitwise comparison of two netCDF files (used for the numerics-unchanged proof). |
| `gpu/scripts/compare_layouts.py` | compares two exports of one case run with different `layout` strings. |

## 2. How to build and run

```
gpu/scripts/build_reference.sh            # reference build: stock flags + -ffp-contract=off
gpu/scripts/build_reference.sh --stock    # the original flags (numerics-unchanged proof)
gpu/scripts/make_refdata.sh               # all gpu/cases/*.in -> gpu/refdata/<case>/
gpu/scripts/make_refdata.sh nl_box_ae lin_single_ae      # selected cases
GS2_REFDATA_NSTEP_TRACE=20 gpu/scripts/make_refdata.sh   # 20 traced steps
```

The exact make commands are in the header of `build_reference.sh`. The stock command is the one that built
the original `bin/gs2`:

```
make -j8 -I Makefiles GK_SYSTEM=macosx-homebrew FFT_INC=-I/opt/homebrew/include \
     FFT_LIB="-L/opt/homebrew/lib -lfftw3" NETCDF_INC="-I/opt/homebrew/include" \
     NETCDF_LIB="-L/opt/homebrew/lib -lnetcdf -lnetcdff" gs2 gs2_export
```

preceded by the same command with target `depend`. The reference build adds the single command-line
variable `F90OPTFLAGS="-O3 -ffp-contract=off"` (`F90OPTFLAGS` is the Fortran optimisation variable:
`Makefile:149` defines it, `Makefiles/enable_features.make` appends `-O3`, `Makefile:547` appends it to
`F90FLAGS`, hence also to `LDFLAGS`). `bin/gs2 --build-config` then reports
`... -mcpu=apple-m1 ... -O3 ... -ffp-contract=off`; the reference objects contain no `fmadd/fmsub/fnmadd/
fnmsub/fmla/fmls` instruction (2094 of them in the stock `dist_fn.o`). The mode is recorded in
`objs/.gs2_gpu_build_mode`; a different requested mode runs `make clean` first (full rebuild, ~1-2 min).
A hand-typed `make` in between is not tracked (use `--clean`).

Running the exporter (single MPI rank only):

```
mpirun -np 1 bin/gs2_export case.in [knobs.nml]
```

`&gs2_export_knobs` is read from `knobs.nml` (command argument 2) or else from the input file; GS2 itself
ignores that group. `bin/gs2_export --build-config` / `--version` work (it calls `parse_command_line`).
Two runs are needed per case because `dt_ladder_n > 0` and `nstep_trace > 0` are mutually exclusive
(09 §8.3): `make_refdata.sh` does an export run (with the ladder) and a separate traced run.

`&gs2_export_knobs` (09 §8.1; defaults in brackets): `export_file` [`<run>.gs2export.nc`],
`write_response` [T], `write_collisions` [T], `write_g_adj` [F], `write_raw_glo` [F],
`write_restart_t0` [F], `run_diagnostics` [T], `nstep_trace` [0], `trace_file` [`<run>.gs2trace.nc`],
`trace_stages` [`'all'`; or `'step'`, or a space-separated tag list], `trace_first_step` [1], `trace_every`
[1], `trace_response_build` [F, no effect, §3.2], `dt_ladder_n` [0], `allow_multi_rank` [F, no effect,
§3.3], `strict` [T], `compression` [0; `make_refdata.sh` uses 1 — shuffle+deflate, verified lossless].

## 3. Deviations from spec 09

### 3.1 `/species` is `/species_`
netCDF-4/HDF5 cannot hold a root dimension `species` and a root group `/species` (HDF5 dimension-scale
name clash: the first write into the group fails with "NetCDF: HDF error"; reproduced with netCDF4-python).
The root dimension `species` is kept (every g-like array uses it) and the group is **`/species_`** (same
contents: `spec_<component>`, `nonmaxw_corr`, `f0_values`, `dlnf0drho`, `has_*` flags, `nspec`). The
root attribute `group_renames = "species:species_"` records it for the loader. No other root
dimension/group name collides.

### 3.2 Probe P-R1 (`response_am`, pre-inversion `am`) is not available
Any probe inside `fields_implicit::init_response_matrix` — a plain `probe_c2('response_am', i, am)` as in
09 §5.2, or one passing a copy `(am)` — changes how gfortran compiles the inlined Gauss-Jordan inversion
(`init_fields_implicit`: 29 → 50 resp. 44 FMA instructions; aliasing/inlining decisions) and changed the
Cyclone test output at the 1e-16 level. P-R1 is optional (10 §A.10), so it was removed: the call site is
not in the sources; `probe_c2` and its hook pointer remain in the probe module (unused). The knob
`trace_response_build` only prints a warning and adds `/fields_implicit:am_note`; `am_c*` is listed in
`absent`. The JAX probe builder can be validated against the exported inverse `aminv_c<c>` instead.

### 3.3 No multi-rank export
The gathers of 09 §2.5 are not implemented. `gs2_export` aborts if `nproc > 1` regardless of
`allow_multi_rank` (reference data must be single-rank anyway, 09 §11.1, 00 D2). The writer and the
canonicalisers read the rank-0 g_lo/le_lo arrays directly.

### 3.4 Extra `public` statements (not listed in 09 §4.2)
* `gs2_time.f90`: `public :: user_dt_cfl` (time block of 09 §7.6; it is always equal to `code_dt_cfl`).
* `init_g.f90`: `public :: ginitopt_switch` (10 §A.5 lists it).
`scale`, `ginit_option`, `restart_dir`, `read_many`, `include_explicit_source_in_restart` come from
`get_init_g_config()`; the `gs2_reinit` knobs from `get_reinit_config()`; the fields knobs from
`get_fields_config()`; the layout attributes from `get_layouts_config()` — all read before any ladder
rung (G-20).

### 3.5 Probe placement details (re-anchored on the Fortran)
* `coll_h_in` (le path) is placed after the `end if` of the `if (adjust)` block (le_derivatives.fpp:73),
  not inside it, so it fires whenever the collision branch runs; standard path: after line 135.
* `coll_gle_out` is after the `end if` of `if (colls)` (le_derivatives.fpp:91); `coll_h_out` after the
  `scatter` (94) and, standard path, after `solfp1` (140).
* `coll_after_drag` is after the `end if` of `if (drag)`, inside `if (has_lorentz)`; the ediffuse probes
  inside `if (has_diffuse)`.
* `nl_cfl` fires after both blocking calls of `set_cfl_time_step_and_check_if_too_large`
  (nonlinear_terms.fpp:569, 582). With `immediate_reset = .false.` the max-velocity reduction is
  non-blocking and completed in `gs2_reinit::check_time_step` after `step_end`; there is no probe there,
  so `code_dt_cfl` of that step appears first in the next `advance_begin`.
* `after_parity` is after `if (def_parity) call enforce_parity(...)` (dist_fn.fpp timeadv).
* `step_end` is after `if (state%exit) state%converged = .true.` (before `check_time_step`), as in 09.
* All other probes are at the positions of 09 §5.2 (line numbers matched the Fortran).

### 3.6 Hook behaviour beyond 09 §5.3
* Stages that do nothing in a run are not recorded: `after_exb_shear` only if `kx_shift` or
  `theta0_shift` is allocated, `after_hyper` only if `hyper_on`, `after_vspace_derivatives` only if
  `colls`, `after_parity` only if `def_parity`, `after_remove_zonal_flows` only if
  `remove_zonal_flows_switch`, `after_collisions_advance` only if `split_collisions`.
* `after_timeadv` records also contain `phi, apar, bpar, phinew, aparnew, bparnew` (so that the pass-1
  invariant `phinew == phi` and `after_field_update = after_getfield + phi` can be checked).
* Names of the explicit-array probe variables: `h` (`coll_h_in`, `coll_h_out`), `g1` (`nl_before_center`,
  `nl_after_center`), `gle` (all LE probes, LE6 complex `(ri, xi_p1, egrid_p1, theta, kx, ky, species)`),
  `max_vel_components` (`nl_cfl`, dimension `three`).
* `step_end` records call `nf90_sync` so a partial trace is readable.

### 3.7 Other differences
* The export-only run (`nstep_trace = 0`) calls `finalize_diagnostics` right after `/state_t0` is
  written and **before** the dt ladder, so `export/<case>.out.nc` holds GS2's static part, the istep-0
  record and the "final" variables (e.g. `phi`) of the t = 0 state (10 G-05); the ladder (which
  re-initialises the equations) cannot pollute it.
* `write_used_inputs_file` is called in both diagnostics branches.
* `/fields_local`: the supercell matrices come from `fields::dump_response_to_file` with a temporary
  `response_file = <run>.gs2export_tmp` (files deleted after reading); sizes from
  `fieldmat%kyb(ik)%supercells(is)%nrow` (components of the public variable are readable). Root
  dimensions `lr_ax_ik<ik>_is<is>`; attributes `ik, is, it_leftmost, condition, dt, layout`; scalars
  `it_leftmost_ik<ik>_is<is>`.
* `/collisions` always holds `collision_model_switch`, `colls`, `split_collisions` and `vnmult`, also when
  the collision operator is inactive (`collisions_advance` re-solves the fields from h whenever
  `split_collisions`, even with collisions off); the other collision variables only if collisions are
  active (`colls .and. collision_model_switch /= none`), otherwise the group has a `note` attribute.
* `use_le_layout = .false.` (gated, G-24): the standard-layout g-like vectors `s0, w0, z0, bs0, bw0, bz0`
  are exported if allocated; the LZ6/E7 matrices `c1, betaa, ql, ec1, ebetaa, eql` are not (listed in
  `absent`).
* `le_lo` is only initialised when collisions are set up; otherwise `ile_of` and `le_lo_*_ord` are listed
  in `/gs2_layouts:absent`. The root dimension `xi_lim` exists whenever `nxi_lim > 0`.
* `absent` lists may contain wildcard entries ending in `*` (`am_c*`, `aminv_c*`, `resp_*`).
* `fortran_decl` is written by the canonical (G7/G6/G7m1/GL/LE6) writers only. Every variable has
  `fortran_type`, `canonical`, and where applicable `fortran_lbounds`, `fortran_name`, `source`,
  `owner_spec`, `dt_dependent` + `code_dt`, `index_base`.
* The ladder zeroes `gnew` before going down to `override_timestep` (as `reset_time_step` does,
  gs2_reinit.f90:163); it does not use the initial-value overrides (the state is discarded).
* Scalars of the time block are also written to `/gs2_time`; `exb_first` also to `/dist_fn`;
  `nl_istep_last`, `nl_initializing` also to `/nonlinear_terms`; `vnmult` also to `/collisions`. The
  authoritative JAX initial state is `/state_t0`.
* `/raw` is not a separate group: with `write_raw_glo = .true.` the state groups get `gnew_glo`
  `(ri, theta, sign, glo)`.
* Hard gates: abort on single precision, SHMEM, LOWFLOW, OpenMP > 1 thread (`strict`), `nproc > 1`,
  list/ensemble runs, `dt_ladder_n > 0 .and. nstep_trace > 0`, `do_eigsolve`; `nstep_trace > nstep` is
  clamped. Soft gates (10 §B, 00 D7) are written to `/gates` (§5.4); `gate_status` lists the violated ones.

## 4. Probe sites (tags) actually in the sources

| tag | file : position | recorded data (besides the time block) |
|---|---|---|
| `advance_begin` | fields.f90 `advance`, before `select case` | `gnew, phinew, aparnew, bparnew, gexp_1..3, kx_shift, theta0_shift` |
| `after_exb_shear` | fields_implicit.f90 / fields_local.f90 after the `exb_shear` line | `gnew, phinew, aparnew, bparnew, kx_shift/theta0_shift, jump, ikx_indexed` |
| `nl_before_center`, `nl_after_center` | nonlinear_terms.fpp `add_explicit`, after `add_nl` / `nl_center` | `g1` (the dummy associated with `gexp_1`) |
| `nl_cfl` | nonlinear_terms.fpp after `set_cfl_time_step_and_check_if_too_large` (2 sites) | `max_vel_components` |
| `after_add_explicit_terms` | dist_fn.fpp `timeadv`, before the reset return | `gexp_1..3` |
| `source` | dist_fn.fpp `invert_rhs`, after `get_source_term_on_local_domain` | `g_work` |
| `after_invert_rhs_1` | dist_fn.fpp `invert_rhs`, after `invert_rhs_1` | `gnew` |
| `after_invert_rhs` | dist_fn.fpp `timeadv`, after `invert_rhs` (incl. linked BCs) | `gnew` |
| `after_hyper` | dist_fn.fpp `timeadv`, after `hyper_diff` | `gnew, hypervisc_filter` |
| `coll_h_in`, `coll_h_out` | le_derivatives.fpp (le and standard paths) | `h` |
| `coll_gle_in`, `coll_gle_out` | le_derivatives.fpp, after `gather` / after the `solfp1` block | `gle` |
| `coll_after_ediffuse`, `coll_after_ediffuse_conserve`, `coll_after_drag`, `coll_after_lorentz`, `coll_after_lorentz_conserve` | collisions.fpp `solfp1_le_layout` | `gle` |
| `after_vspace_derivatives` | dist_fn.fpp `timeadv` | `gnew` |
| `after_parity` | dist_fn.fpp `timeadv` | `gnew` |
| `after_timeadv` (pass 1 and 2) | fields_implicit.f90 / fields_local.f90 after each `timeadv` | `gnew, phi, apar, bpar, phinew, aparnew, bparnew` |
| `after_getfield`, `after_field_update`, `after_remove_zonal_flows` | fields_implicit.f90 / fields_local.f90 | `phinew, aparnew, bparnew` |
| `after_collisions_advance` | after `collisions_advance` | `gnew` |
| `advance_end` | fields.f90 `advance`, after `end select` | as `advance_begin` + `phi, apar, bpar` |
| `after_reset_time_step` | gs2_main.fpp after both `reset_time_step` calls | time block only |
| `step_end_pre_diag` | gs2_main.fpp after `update_time` | as `advance_begin` |
| `step_end` | gs2_main.fpp after `state%converged` | as `advance_begin` + `g` |

`probe_pass` is set to 1/2 before the two `timeadv` calls of `advance_implicit`/`advance_local` and reset
to 0 by `advance_end`. Inner tags are recorded only while `probe_in_advance` (the `istep = 0` calls of
the response-matrix build are filtered out).

## 5. Files and names actually written

### 5.1 Export file `<run>.gs2export.nc` (NETCDF4)

Global attributes: GS2 metadata (`title = "GS2 GPU export"`, `software_name`, `software_version`,
`netcdf_version`, `date_created`, `id`, `Conventions`, `references`) and `export_format_version = "1.0"`,
`git_hash`, `git_state`, `build_config`, `run_name`, `input_file`, `nproc`, `nproc_g_lo`, `nproc_le_lo`,
`layout`, `index_convention`, `real_kind_bytes = 8`, `group_renames = "species:species_"`, `gate_status`,
`diagnostics_run_at_t0`, `used_inputs` (text of `<run>.used_inputs.in`), `input_file_text`.

Root dimensions (coordinate variable in brackets): `ri` (1,2), `theta` (theta), `theta_m1`, `sign`, `kx`
(akx), `ky` (aky), `lambda` (al), `egrid`, `egrid_p1`, `egrid_m1`, `species` (1..nspec), `xi2`, `xi_p1`,
`xi_lim` (if `nxi_lim > 0`), `bset` (bset, if `nbset > 0`), `two`, `three`, `n_class`, `sc_Mmax`,
`sc_Nmax`, `sc_R_c<c>`, `sc_M_c<c>`, `sc_RM_c<c>` (implicit fields), `lr_ax_ik<ik>_is<is>` (local
fields), `glo` (`write_raw_glo`), `nlinks_max` (`write_g_adj`).

Groups: `/theta_grid`, `/kt_grids`, `/le_grids`, `/species_`, `/run_parameters`, `/gs2_time`,
`/gs2_reinit`, `/init_g`, `/fields`, `/gs2_layouts`, `/dist_fn_arrays`, `/dist_fn`, `/fields_implicit`
or `/fields_local`, `/collisions`, `/hyper`, `/nonlinear_terms`, `/state_t0_prediag`, `/state_t0`,
`/diagnostics`, `/gates`, `/dt_ladder/k<k>/{dist_fn, dist_fn_arrays, collisions, fields_implicit|
fields_local, hyper, nonlinear_terms, gs2_time}` (+ time block and attribute `code_dt` in `/dt_ladder/k<k>`;
`/dt_ladder` attributes `n_rungs`, `delt_adj`). The variables of each group (from real export files) are
listed in §9.

### 5.2 Trace file `<run>.gs2trace.nc`

Same global attributes and root dimensions, plus `nstep_trace`, `trace_stages`, `trace_first_step`,
`trace_every`, `n_records`. Records `/s<istep:06d>/a<attempt>/p<seq:03d>_<tag>` with attributes `tag`,
`istep`, `attempt`, `seq`, `pass`, `in_advance`, the time block (`code_dt, user_dt, code_dt_prev1,
code_dt_prev2, code_dt_old, code_dt_cfl, user_dt_cfl, code_dt_min, code_dt_max, code_time, user_time,
istep, reset, vnmult(two), nl_istep_last, nl_initializing, exb_first`) and the arrays of §4. `/final`:
state after `evolve_equations` returned (after the forced final diagnostics call).

### 5.3 Conventions (unchanged from 09 §0.2)

Fortran dimension order (NumPy reversed); complex = trailing NumPy `ri` axis; logical = int8 with
`fortran_type = "logical"`; integers keep Fortran values (1-based indices, −1 sentinels); `iglo_of`,
`ile_of` are 0-based (`index_base = 0`). `aminv_c<c>`: NumPy `[m, row, col, ri]`, rows/cols in B order
(`(ig+ntgrid) + nth*(ifield-1) + nidx*(n-1)`), raw runtime rows (no duplicate removal; zonal rows 0).

### 5.4 Gates written to `/gates`
`hard_single_precision, hard_lowflow, hard_multi_rank, hard_do_eigsolve` (always 0 in a written file),
`antenna_on, eqzip, afilter, split_nonlinear, zip, field_option_test_or_gf_local, gf_lo_integrate,
gf_local_fields, fft_measure_plan, new_field_init_false, wstar_units, use_old_diagnostics,
phiext_full_with_phi_ext, opt_source, linked_non_box, experimental_boundary, nonlinear_non_box,
exb_single_kx, hyper_zero_normaliser, scale_not_one, hybrid_with_collisions,
nonmaxwellian_with_collisions, use_le_layout_false, no_trapped_lorentz_conserve, conservative_false,
old_collision_scheme, Q06_bmax_not_at_ends, Q07_zero_shear_with_zonal, Q09_bouncefuzz_edge` and, when
diagnostics are initialised, `diag_*` for every gated diagnostics flag of 10 §B.7. Each has a
`description` attribute; 1 = violated.

## 6. Numerics-unchanged proof (09 §4.4)

See §10 for the results of the runs. Method:
1. Object code (step 2 of 09 §4.4): the `public`-only edits were built with stock flags and compared
   with the objects of the original build (`otool -tv`). 8 of 10 edited modules (`collisions`, `hyper`,
   `nonlinear_terms`, `fields_implicit`, `le_grids`, `run_parameters`, `gs2_time`, `init_g`) are
   instruction-for-instruction identical except the source-line-number immediates passed to the
   runtime library (they shift by the number of inserted lines). `dist_fn.o` gains one out-of-line copy
   of `get_rightmost_it` (24 integer instructions; it was private and fully inlined before; all callers
   unchanged). `theta_grid.o` differs only in the scheduling of the argument set-up of one I/O call in
   `check_theta_grid`. No floating-point code changed.
2. Outputs (step 3): a pristine worktree of `HEAD` (`git worktree add --detach`, submodule contents
   copied) and the edited tree were both built with the stock command; GS2 test inputs were run with
   `mpirun -np 1` in both and every variable of `<run>.out.nc` compared bitwise
   (`gpu/scripts/compare_netcdf.py`).
3. Trajectory with active probes: for every refdata case `check_export.py` compares the traced
   `gs2_export` run's `<run>.out.nc` with the full `bin/gs2` run at the time records both contain.

## 7. Reference data layout (`gpu/refdata/<case>/`)

| path | content |
|---|---|
| `gs2/` | full `bin/gs2` run (`<case>.out.nc`, restart files, ascii outputs, `gs2.log`) |
| `export/` | `gs2_export` export run with `dt_ladder_n = 2`, `write_restart_t0 = .true.`: `<case>.gs2export.nc`, `<case>.out.nc` (static part + istep-0 record + final variables at t = 0), the GS2 restart file of the t = 0 state `<restart_file>.gs2export_t0.0`, `<case>.used_inputs.in`, `gs2_export.log` |
| `trace/` | traced `gs2_export` run (`nstep_trace = N`, `trace_stages = 'all'`): `<case>.gs2trace.nc`, the traced run's GS2 outputs (`<case>.out.nc` with records up to step N), `gs2_export.log` |
| `<case>.gs2export.nc`, `<case>.gs2trace.nc`, `<case>.out.nc` | symlinks to `export/`, `trace/`, `gs2/` |
| `check_export.log`, `STATUS` | output of `check_export.py`; one-line PASS/FAIL summary |

Export and trace files are written with `compression = 1` (shuffle + deflate, lossless: bitwise identical
to uncompressed files). Environment knobs: `GS2_REFDATA_NSTEP_TRACE` (N, default 10),
`GS2_REFDATA_DT_LADDER` (2), `GS2_REFDATA_COMPRESSION` (1), `GS2_REFDATA_DIR`,
`GS2_REFDATA_KEEP_TRACE_EXPORT`, `GS2_REFDATA_BUILD`, `GS2_REFDATA_ALLOW_STOCK`,
`GS2_REFDATA_ONLY=export` (regenerate only `export/` of the listed cases, keep `gs2/` and `trace/`, re-run
the checks against them), `MPIRUN`. The script
refuses non-reference binaries unless `GS2_REFDATA_ALLOW_STOCK=1`, unsets `GK_FFTW3_WISDOM`, sets
`OMP_NUM_THREADS=1`, uses fresh run directories (no stale FFTW wisdom), copies `!include`d files, and
picks up every `gpu/cases/*.in` present when it starts.

## 8. What `check_export.py` checks

* Schema: every name of the 10 §A master list (per group, with NumPy shape and dtype) present or listed
  in the group's `absent`; `fortran_type`/`canonical` attributes; complex variables end with `ri`;
  gates consistent with `gate_status`; no hard gate.
* Spot checks (recomputed in NumPy with the GS2 association order): `vperp2`, `vpa`, `vpac`, `vpar`,
  `wstar`, `tunits/wunits/woutunits`, `has_*`, `field_line_average_weight` (Klein-summed trapezoid,
  bitwise) and its exporter recompute, `kperp2`, `inv_kperp2`, `iglo_of` and `ile_of` from the layout
  string, `kwork_filter`, zero `aminv` rows of the zonal member and class-table coverage, dt-ladder rungs
  (`code_dt = dt/delt_adj^k` bitwise, `vpar`/`wstar` recomputed with the rung dt, `D_res` frozen,
  `aminv` rebuilt), `state_t0` field identities, 23 geometry/grid variables vs GS2's `.out.nc`,
  `phi_igomega_by_mode(t=0)` and the final `phi` vs `/state_t0/phinew`. FMA-prone expressions (`vpa`,
  `kperp2`) must be bitwise for a reference build (the build flags are read from the file).
* Trace: `advance_begin(1) == /state_t0`; `step_end(k) == advance_begin(k+1)` (time block too; with
  `immediate_reset = .false.` without `code_dt_cfl`; after a deferred reset against the
  `after_reset_time_step` record, gnew/gexp numerically, fields within round-off); pass-1
  `phinew == phi`; `after_field_update == after_getfield + phi` (bitwise); `advance_end == after_timeadv`
  pass 2 (or `after_collisions_advance`); record order `source` → `after_invert_rhs_1`; NL history shuffle
  and pass-2 reuse; no `nl_*` records in attempt 2; `update_time` arithmetic; `advance_end ==
  step_end_pre_diag`.
* LE6 map (collisional cases): `aj0le`, `vperp_aj1le`, `vpa_aj0_le` rebuilt bitwise from `aj0`, `aj1`,
  `vpa` with `ixi_to_il/ixi_to_isgn` (g2le skips forbidden points); `e15 == energy_maxw**1.5`; in the
  trace `coll_gle_in` = gather(`coll_h_in`) and `coll_h_out` = scatter(`coll_gle_out`) with forbidden
  points of `h` unchanged.
* Restart (`--restart`, 09 §12 T2): `iglo_of == iglo_kykxsel`, and `gr/gi`, `vpa`, `vperp2` of GS2's own
  restart file re-indexed with `iglo_of` equal `/state_t0/gnew`, `vpa`, `vperp2` bitwise.
* Determinism: the trace run's own export file is bitwise identical to the export run's (except
  `/dt_ladder`).
* Trajectory: the traced run's time-dependent `.out.nc` variables equal bin/gs2's at the common records.

## 9. Variable inventory (generated from real export files)

Union over all 28 `gpu/refdata` export files (symbols: `S` nspec, `E` negrid, `L` nlambda, `Y` naky,
`X` ntheta0, `T` 2*ntgrid+1, `XI` nxi_lim, `C` i_class; NumPy order; `c (ri)` = complex with trailing
`ri`; `i4`/`i1` = int32 / int8-logical, unmarked = float64). `/dt_ladder/k<k>` repeats the `k1` layout.

| group | variables |
|---|---|
| `/collisions` | `collision_model_switch`() i4; `colls`() i1; `split_collisions`() i1 (always); `vnmult`(2); `ec1le`(S,Y,X,T,E,XI); `ebetaale`(S,Y,X,T,E,XI); `eqle`(S,Y,X,T,E,XI); `bz0le`(S,Y,X,T,E+1,XI); `bs0le`(S,Y,X,T,E+1,XI); `bw0le`(S,Y,X,T,E+1,XI); `vnew`(S,E,Y); `vnew_s`(S,E,Y); `vnew_D`(S,E,Y); `vnew_E`(S,E,Y); `delvnew`(S,E,Y); `aj0le`(S,Y,X,T,E+1,XI); `vperp_aj1le`(S,Y,X,T,E+1,XI); `vpa_aj0_le`(S,Y,X,T,E+1,XI); `e15`(E); `capg_half`(E-1); `nxi_lim`() i4; `cfac`(); `drag`() i1; `has_lorentz`() i1; `has_diffuse`() i1; `const_v`() i1; `conserve_moments`() i1; `conservative`() i1; `resistivity`() i1; `lorentz_switch`() i4; `ediff_switch`() i4; `ei_coll_only`() i1; `special_wfb_lorentz`() i1; `vpar_zero_mean`() i1; `conserve_forbid_zero`() i1; `force_collisions`() i1; `hypermult`() i1; `vnfac`(); `vnslow`(); `vnm_init`(2); `timesteps_between_collisions`() i4; `use_le_layout`() i1; `adjust`() i1; `heating`() i1; `hyper_colls`() i1; `vary_vnew`() i1; `ncheck`() i4; `c1le`(S,Y,X,T,E,XI); `betaale`(S,Y,X,T,E,XI); `qle`(S,Y,X,T,E,XI); `z0le`(S,Y,X,T,E+1,XI); `s0le`(S,Y,X,T,E+1,XI); `w0le`(S,Y,X,T,E+1,XI); `pitch_weights`(T,L); `vpdiff`(L,2,T); `vpdiffle`(T,XI) |
| `/diagnostics` | `nout`() i4; `istep`() i4; `nwrite`() i4; `nwrite_mult`() i4; `navg`() i4; `igomega`() i4; `omegatinst`(); `omegatol`(); `exit_when_converged`() i1; `calculate_fluxes`() i1; `user_time_old`(); `start_time`(); `ncheck`() i4; `nsave`() i4; `nc_sync_freq`() i4; `nmovie`() i4; `conv_nstep_av`() i4; `conv_test_multiplier`(); `conv_min_step`() i4; `conv_max_step`() i4; `conv_nsteps_converged`() i4; `write_any`() i1; `append_old`() i1; `enable_parallel`() i1; `serial_netcdf4`() i1; `print_line`() i1; `print_flux_line`() i1; `write_line`() i1; `write_flux_line`() i1; `write_fields`() i1; `write_phi_over_time`() i1; `write_apar_over_time`() i1; `write_bpar_over_time`() i1; `make_movie`() i1; `dump_fields_periodically`() i1; `write_moments`() i1; `write_full_moments_notgc`() i1; `write_ntot_over_time`() i1; `write_density_over_time`() i1; `write_upar_over_time`() i1; `write_tperp_over_time`() i1; `write_fluxes`() i1; `write_fluxes_by_mode`() i1; `write_symmetry`() i1; `write_nl_flux_dist`() i1; `write_parity`() i1; `write_omega`() i1; `write_ql_metric`() i1; `write_verr`() i1; `write_cerr`() i1; `write_max_verr`() i1; `write_heating`() i1; `write_ascii`() i1; `write_gyx`() i1; `write_g`() i1; `use_nonlin_convergence`() i1; `write_cross_phase`() i1; `write_correlation`() i1; `write_correlation_extend`() i1; `write_jext`() i1; `write_lorentzian`() i1; `write_eigenfunc`() i1; `write_final_fields`() i1; `write_kpar`() i1; `write_final_epar`() i1; `write_final_db`() i1; `write_final_moments`() i1; `write_final_antot`() i1; `write_gs`() i1; `save_for_restart`() i1; `save_many`() i1; `file_safety_check`() i1; `save_distfn`() i1; `save_glo_info_and_grids`() i1; `save_velocities`() i1; `write_zonal_transfer`() i1; `write_kinetic_energy_transfer`() i1; `write_collisional`() i1; `write_omavg`() i1; `ob_midplane`() i1; `write_avg_moments`() i1 |
| `/dist_fn` | `vpar`(S,E,L,Y,X,2,T); `wdrift`(S,E,L,Y,X,2,T); `wdriftttp`(S,E,L,Y,X,2,T); `a`(S,E,L,Y,X,2,T,ri) c (ri); `b`(S,E,L,Y,X,2,T,ri) c (ri); `r`(S,E,L,Y,X,2,T,ri) c (ri); `ainv`(S,E,L,Y,X,2,T,ri) c (ri); `g_h`(S,E,L,Y,X,2,T,ri) c (ri); `ufac`(S,E); `gamtot`(Y,X,T); `gamtot1`(Y,X,T); `gamtot2`(Y,X,T); `gamtot3`(Y,X,T); `gridfac1`(Y,X,T); `inv_phi_denominator_g`(Y,X,T); `awgt`(Y,X); `itleft`(Y,X) i4; `itright`(Y,X) i4; `l_links`(Y,X) i4; `r_links`(Y,X) i4; `n_links`(Y,X,2) i4; `M_class`(C) i4; `N_class`(C) i4; `i_class`() i4; `supercell_labels`(Y,X) i4; `n_cells`(Y,X) i4; `n_supercells`(Y) i4; `itleftmost`(Y,X) i4; `itrightmost`(Y,X) i4; `connections_neighbor`(S,E,L,Y,X) i1; `connections_it_left`(S,E,L,Y,X) i4; `connections_it_right`(S,E,L,Y,X) i4; `save_h`(S,E,L,Y,X,2) i1; `no_connections`() i1; `exb_first`() i1; `wfb_aj0l`(S,Y,X,E); `wfb_aj1l`(S,Y,X,E); `wfb_aj0r`(S,Y,X,E); `wfb_aj1r`(S,Y,X,E); `wfb_vperp2left`(S,E); `wfb_vperp2right`(S,E); `fexp`(S); `bkdiff`(S); `gridfac`(); `apfac`(); `driftknob`(); `tpdriftknob`(); `poisfac`(); `vparknob`(); `t0`(); `omega0`(); `gamma0`(); `source0`(); `phi_ext`(); `afilter`(); `exponential_boundary_factor`(); `wfb`(); `g_exb`(); `g_exbfac`(); `omprimfac`(); `btor_slab`(); `mach`(); `g_exb_start_time`(); `g_exb_error_limit`(); `g_exb_start_timestep`() i4; `hyper_in_initialisation`() i1; `start_from_previous_solution`() i1; `boundary_off_grid`() i1; `exponential_boundary`() i1; `nonad_zero`() i1; `lf_default`() i1; `lf_decompose`() i1; `esv`() i1; `opt_source`() i1; `source_option_switch`() i4; `boundary_option_switch`() i4; `adiabatic_option_switch`() i4; `def_parity`() i1; `even`() i1; `zero_forbid`() i1; `gf_lo_integrate`() i1; `mult_imp`() i1; `inv_bpar_denominator_g`(Y,X,T) |
| `/dist_fn_arrays` | `vpa`(S,E,L,Y,X,2,T); `vpac`(S,E,L,Y,X,2,T); `vperp2`(S,E,L,Y,X,T); `aj0`(S,E,L,Y,X,T); `aj1`(S,E,L,Y,X,T); `wstar`(S,E,Y) |
| `/dt_ladder/k1` | `code_dt`(); `user_dt`(); `code_dt_prev1`(); `code_dt_prev2`(); `code_dt_old`(); `code_dt_cfl`(); `user_dt_cfl`(); `code_dt_min`(); `code_dt_max`(); `code_time`(); `user_time`(); `istep`() i4; `reset`() i1; `vnmult`(2); `nl_istep_last`() i4; `nl_initializing`() i1; `exb_first`() i1 |
| `/dt_ladder/k1/collisions` | `collision_model_switch`() i4; `colls`() i1; `vnmult`(2); `ec1le`(S,Y,X,T,E,XI); `ebetaale`(S,Y,X,T,E,XI); `eqle`(S,Y,X,T,E,XI); `bz0le`(S,Y,X,T,E+1,XI); `bs0le`(S,Y,X,T,E+1,XI); `bw0le`(S,Y,X,T,E+1,XI); `c1le`(S,Y,X,T,E,XI); `betaale`(S,Y,X,T,E,XI); `qle`(S,Y,X,T,E,XI); `z0le`(S,Y,X,T,E+1,XI); `s0le`(S,Y,X,T,E+1,XI); `w0le`(S,Y,X,T,E+1,XI) |
| `/dt_ladder/k1/dist_fn` | `vpar`(S,E,L,Y,X,2,T); `wdrift`(S,E,L,Y,X,2,T); `wdriftttp`(S,E,L,Y,X,2,T); `a`(S,E,L,Y,X,2,T,ri) c (ri); `b`(S,E,L,Y,X,2,T,ri) c (ri); `r`(S,E,L,Y,X,2,T,ri) c (ri); `ainv`(S,E,L,Y,X,2,T,ri) c (ri); `g_h`(S,E,L,Y,X,2,T,ri) c (ri) |
| `/dt_ladder/k1/dist_fn_arrays` | `wstar`(S,E,Y) |
| `/dt_ladder/k1/fields_implicit` | `nfield`() i4; `nidx`() i4; `cls_ik`(C,Nmax,Mmax) i4; `cls_it`(C,Nmax,Mmax) i4; `sc_class`(Y,X) i4; `sc_member`(Y,X) i4; `sc_cell`(Y,X) i4; `aminv_c<c>`(M_c,R_c,R_c,ri) c (ri) |
| `/dt_ladder/k1/fields_local` | `resp_ik<ik>_is<is>`(n_sc,n_sc,ri) c (ri); `it_leftmost_ik<ik>_is<is>`() i4 |
| `/dt_ladder/k1/gs2_time` | `tunits`(Y); `wunits`(Y); `woutunits`(Y); `dt_not_set`(); `code_dt`(); `user_dt`(); `code_dt_prev1`(); `code_dt_prev2`(); `code_dt_old`(); `code_dt_cfl`(); `user_dt_cfl`(); `code_dt_min`(); `code_dt_max`(); `code_time`(); `user_time`(); `istep`() i4; `reset`() i1; `vnmult`(2); `nl_istep_last`() i4; `nl_initializing`() i1; `exb_first`() i1 |
| `/dt_ladder/k1/hyper` | `D_res`(Y,X); `hypervisc_filter`(Y,X,T) |
| `/dt_ladder/k1/nonlinear_terms` | `cflx`(); `cfly`() |
| `/fields` | `fieldopt_switch`() i4; `remove_zonal_flows_switch`() i1; `force_maxwell_reinit`() i1; `dump_response`() i1; `read_response`() i1; `field_subgath`() i1; `do_smart_update`() i1; `minnrow`() i4 |
| `/fields_implicit` | `nfield`() i4; `nidx`() i4; `cls_ik`(C,Nmax,Mmax) i4; `cls_it`(C,Nmax,Mmax) i4; `sc_class`(Y,X) i4; `sc_member`(Y,X) i4; `sc_cell`(Y,X) i4; `aminv_c<c>`(M_c,R_c,R_c,ri) c (ri) |
| `/fields_local` | `resp_ik<ik>_is<is>`(n_sc,n_sc,ri) c (ri); `it_leftmost_ik<ik>_is<is>`() i4 |
| `/gates` | `hard_single_precision`() i1; `hard_lowflow`() i1; `hard_multi_rank`() i1; `hard_do_eigsolve`() i1; `antenna_on`() i1; `eqzip`() i1; `afilter`() i1; `split_nonlinear`() i1; `zip`() i1; `field_option_test_or_gf_local`() i1; `gf_lo_integrate`() i1; `gf_local_fields`() i1; `fft_measure_plan`() i1; `new_field_init_false`() i1; `wstar_units`() i1; `use_old_diagnostics`() i1; `phiext_full_with_phi_ext`() i1; `opt_source`() i1; `linked_non_box`() i1; `experimental_boundary`() i1; `nonlinear_non_box`() i1; `exb_single_kx`() i1; `hyper_zero_normaliser`() i1; `scale_not_one`() i1; `hybrid_with_collisions`() i1; `nonmaxwellian_with_collisions`() i1; `use_le_layout_false`() i1; `no_trapped_lorentz_conserve`() i1; `conservative_false`() i1; `old_collision_scheme`() i1; `Q06_bmax_not_at_ends`() i1; `Q07_zero_shear_with_zonal`() i1; `Q09_bouncefuzz_edge`() i1; `diag_use_nonlin_convergence`() i1; `diag_save_distfn`() i1; `diag_append_old`() i1; `diag_make_movie`() i1; `diag_dump_fields_periodically`() i1; `diag_write_heating`() i1; `diag_write_collisional`() i1; `diag_write_verr`() i1; `diag_write_max_verr`() i1; `diag_write_cerr`() i1; `diag_write_symmetry`() i1; `diag_write_nl_flux_dist`() i1; `diag_write_parity`() i1; `diag_write_cross_phase`() i1; `diag_write_correlation_extend`() i1; `diag_write_jext`() i1; `diag_write_lorentzian`() i1; `diag_write_zonal_transfer`() i1; `diag_write_kinetic_energy_transfer`() i1; `diag_write_full_moments_notgc`() i1; `diag_write_g`() i1; `diag_write_gyx`() i1; `diag_write_gs`() i1; `diag_write_kpar`() i1; `diag_write_final_db`() i1 |
| `/gs2_layouts` | `g_lo_ik_ord`() i4; `g_lo_it_ord`() i4; `g_lo_il_ord`() i4; `g_lo_ie_ord`() i4; `g_lo_is_ord`() i4; `g_lo_ulim_world`() i4; `g_lo_blocksize`() i4; `g_lo_nproc`() i4; `iglo_of`(S,E,L,Y,X) i4; `local_field_solve`() i1; `intmom_sub`() i1; `intspec_sub`() i1; `gf_local_fields`() i1; `opt_local_copy`() i1; `fft_use_wisdom`() i1; `fft_measure_plan`() i1; `simple_gf_decomposition`() i1; `le_lo_ig_ord`() i4; `le_lo_it_ord`() i4; `le_lo_ik_ord`() i4; `le_lo_is_ord`() i4; `le_lo_ulim_world`() i4; `ile_of`(S,Y,X,T) i4 |
| `/gs2_reinit` | `delt_adj`(); `in_memory`() i1; `delt_cushion`(); `delt_minimum`(); `dt0`(); `abort_rapid_time_step_change`() i1 |
| `/gs2_time` | `tunits`(Y); `wunits`(Y); `woutunits`(Y); `dt_not_set`(); `code_dt`(); `user_dt`(); `code_dt_prev1`(); `code_dt_prev2`(); `code_dt_old`(); `code_dt_cfl`(); `user_dt_cfl`(); `code_dt_min`(); `code_dt_max`(); `code_time`(); `user_time`(); `istep`() i4; `reset`() i1; `vnmult`(2); `nl_istep_last`() i4; `nl_initializing`() i1; `exb_first`() i1 |
| `/hyper` | `D_res`(Y,X); `D_v`(); `D_eta`(); `nexp`() i4; `D_hypervisc`(); `D_hyperres`(); `omega_osc`(); `D_hyper`(); `akx4_max`(); `aky4_max`(); `aky_max`(); `akperp4_max`(); `kperp2_max`(); `kperp2_max_zonal`(); `hyper_option_switch`() i4; `const_amp`() i1; `include_kpar`() i1; `isotropic_shear`() i1; `damp_zonal_only`() i1; `hyper_on`() i1; `gridnorm`() i1; `D_hyper3D`(); `P_hyper3D`(); `ky_cut`(); `kx_cut`(); `isotropic_model`() i1; `hypervisc_filter`(Y,X,T) |
| `/init_g` | `tstart`(); `new_field_init`() i1; `initial_condition_is_nonadiabatic_dfn`() i1; `ginitopt_switch`() i4; `scale`(); `read_many`() i1; `include_explicit_source_in_restart`() i1 |
| `/kt_grids` | `aky`(Y); `akx`(X); `ikx`(X) i4; `theta0`(Y,X); `kperp2`(Y,X,T); `inv_kperp2`(Y,X,T); `kwork_filter`(Y,X) i1; `naky`() i4; `ntheta0`() i4; `nx`() i4; `ny`() i4; `jtwist_out`() i4; `gridopt_switch`() i4; `reality`() i1; `box`() i1 |
| `/le_grids` | `energy`(S,E); `w`(S,E); `speed`(S,E); `energy_maxw`(E); `w_maxw`(E); `speed_maxw`(E); `al`(L); `wl`(L,T); `jend`(T) i4; `forbid`(L,T) i1; `is_ttp`(L,T) i1; `is_bounce_point`(L,T) i1; `can_be_ttp`(L) i1; `xi`(T,2L); `wxi`(T,2L); `ixi_to_il`(T,2L) i4; `ixi_to_isgn`(T,2L) i4; `sgn`(2) i4; `is_lower_bounce_point`(L,T) i1; `is_upper_bounce_point`(L,T) i1; `il_is_wfb`(L) i1; `il_is_passing`(L) i1; `il_is_trapped`(L) i1; `is_passing_hybrid_electron`(S,Y,L) i1; `grid_has_trapped_particles`() i1; `negrid`() i4; `nlambda`() i4; `ng2`() i4; `lmax`() i4; `nxi`() i4; `nesub`() i4; `nesuper`() i4; `ngauss`() i4; `npassing`() i4; `wdim`() i4; `vcut`(); `bouncefuzz`(); `wfbbc_option_switch`() i4; `new_trap_int`() i1; `new_trap_int_split`() i1; `passing_wfb`() i1; `trapped_wfb`() i1; `mixed_wfb`() i1; `trapped_particles`() i1; `genquad`() i1; `radau_gauss_grid`() i1; `split_passing_region`() i1 |
| `/nonlinear_terms` | `cflx`(); `cfly`(); `nonlin`() i1; `split_nonlinear`() i1; `accelerated`() i1; `cfl`(); `max_vel_components`(3); `nonlinear_mode_switch`() i4; `use_cfl_limit`() i1; `use_2d_cfl`() i1; `error_target`(); `istep_error_start`() i4; `use_order_based_error`() i1; `include_phi`() i1; `include_apar`() i1; `include_bpar`() i1; `zip`() i1; `nl_forbid_force_zero`() i1; `nl_istep_last`() i4; `nl_initializing`() i1 |
| `/run_parameters` | `beta`(); `zeff`(); `tite`(); `fphi`(); `fapar`(); `fbpar`(); `k0`(); `rhostar`(); `delt`(); `max_sim_time`(); `avail_cpu_time`(); `margin_cpu_time`(); `has_phi`() i1; `has_apar`() i1; `has_bpar`() i1; `nstep`() i4; `ncheck_stop`() i4; `progress_frequency`() i4; `delt_option_switch`() i4; `immediate_reset`() i1; `wstar_units`() i1; `neo_test`() i1; `do_eigsolve`() i1; `use_old_diagnostics`() i1; `reset`() i1; `ieqzip`(Y,X) i1; `eqzip_option_switch`() i4; `seed`() i4 |
| `/species_` | `nspec`() i4; `spec_z`(S); `spec_mass`(S); `spec_dens`(S); `spec_temp`(S); `spec_tprim`(S); `spec_fprim`(S); `spec_vnewk`(S); `spec_stm`(S); `spec_zstm`(S); `spec_tz`(S); `spec_smz`(S); `spec_zt`(S); `spec_uprim`(S); `spec_uprim2`(S); `spec_dens0`(S); `spec_u0`(S); `spec_tpar0`(S); `spec_tperp0`(S); `spec_nu_h`(S); `spec_bess_fac`(S); `spec_vcrit`(S); `spec_vcprim`(S); `spec_type`(S) i4; `spec_f0type`(S) i4; `nonmaxw_corr`(S,E); `f0_values`(S,E); `dlnf0drho`(S,E); `has_electron_species`() i1; `has_ion_species`() i1; `has_hybrid_electron_species`() i1 |
| `/state_t0` | `code_dt`(); `user_dt`(); `code_dt_prev1`(); `code_dt_prev2`(); `code_dt_old`(); `code_dt_cfl`(); `user_dt_cfl`(); `code_dt_min`(); `code_dt_max`(); `code_time`(); `user_time`(); `istep`() i4; `reset`() i1; `vnmult`(2); `nl_istep_last`() i4; `nl_initializing`() i1; `exb_first`() i1; `g`(S,E,L,Y,X,2,T,ri) c (ri); `gnew`(S,E,L,Y,X,2,T,ri) c (ri); `phi`(Y,X,T,ri) c (ri); `apar`(Y,X,T,ri) c (ri); `bpar`(Y,X,T,ri) c (ri); `phinew`(Y,X,T,ri) c (ri); `aparnew`(Y,X,T,ri) c (ri); `bparnew`(Y,X,T,ri) c (ri); `theta0_shift`(Y); `gexp_1`(S,E,L,Y,X,2,T,ri) c (ri); `gexp_2`(S,E,L,Y,X,2,T,ri) c (ri); `gexp_3`(S,E,L,Y,X,2,T,ri) c (ri); `kx_shift`(Y) |
| `/state_t0_prediag` | `code_dt`(); `user_dt`(); `code_dt_prev1`(); `code_dt_prev2`(); `code_dt_old`(); `code_dt_cfl`(); `user_dt_cfl`(); `code_dt_min`(); `code_dt_max`(); `code_time`(); `user_time`(); `istep`() i4; `reset`() i1; `vnmult`(2); `nl_istep_last`() i4; `nl_initializing`() i1; `exb_first`() i1; `g`(S,E,L,Y,X,2,T,ri) c (ri); `gnew`(S,E,L,Y,X,2,T,ri) c (ri); `phi`(Y,X,T,ri) c (ri); `apar`(Y,X,T,ri) c (ri); `bpar`(Y,X,T,ri) c (ri); `phinew`(Y,X,T,ri) c (ri); `aparnew`(Y,X,T,ri) c (ri); `bparnew`(Y,X,T,ri) c (ri); `theta0_shift`(Y); `gexp_1`(S,E,L,Y,X,2,T,ri) c (ri); `gexp_2`(S,E,L,Y,X,2,T,ri) c (ri); `gexp_3`(S,E,L,Y,X,2,T,ri) c (ri); `kx_shift`(Y) |
| `/theta_grid` | `theta`(T); `theta2`(T); `delthet`(T); `delthet2`(T); `bset`(nbset); `bmag`(T); `gradpar`(T); `itor_over_B`(T); `IoB`(T); `gbdrift`(T); `gbdrift0`(T); `cvdrift`(T); `cvdrift0`(T); `cdrift`(T); `cdrift0`(T); `gds2`(T); `gds21`(T); `gds22`(T); `gds23`(T); `gds24`(T); `gds24_noq`(T); `grho`(T); `jacob`(T); `Rplot`(T); `Zplot`(T); `aplot`(T); `Rprime`(T); `Zprime`(T); `aprime`(T); `Bpol`(T); `bmin`(); `bmax`(); `eps_trapped`(); `shat`(); `drhodpsi`(); `kxfac`(); `qval`(); `surfarea`(); `dvdrhon`(); `rhoc`(); `field_line_average_weight`(); `field_line_average_weight_check`(); `cvdriftknob`(); `gbdriftknob`(); `ntheta`() i4; `ntgrid`() i4; `nperiod`() i4; `nbset`() i4; `eqopt_switch`() i4; `gb_to_cv`() i1; `is_effectively_zero_shear`() i1 |

Trace records (union over all 28 refdata traces; every record also holds the 0-D time block, hence
`vnmult` everywhere):

| tag | arrays |
|---|---|
| `advance_begin` | `aparnew`, `bparnew`, `gexp_1`, `gexp_2`, `gexp_3`, `gnew`, `kx_shift`, `phinew`, `theta0_shift`, `vnmult` |
| `after_add_explicit_terms` | `gexp_1`, `gexp_2`, `gexp_3`, `vnmult` |
| `source` | `g_work`, `vnmult` |
| `after_invert_rhs_1` | `gnew`, `vnmult` |
| `after_invert_rhs` | `gnew`, `vnmult` |
| `after_timeadv` | `apar`, `aparnew`, `bpar`, `bparnew`, `gnew`, `phi`, `phinew`, `vnmult` |
| `after_getfield` | `aparnew`, `bparnew`, `phinew`, `vnmult` |
| `after_field_update` | `aparnew`, `bparnew`, `phinew`, `vnmult` |
| `advance_end` | `apar`, `aparnew`, `bpar`, `bparnew`, `gexp_1`, `gexp_2`, `gexp_3`, `gnew`, `kx_shift`, `phi`, `phinew`, `theta0_shift`, `vnmult` |
| `step_end_pre_diag` | `aparnew`, `bparnew`, `gexp_1`, `gexp_2`, `gexp_3`, `gnew`, `kx_shift`, `phinew`, `theta0_shift`, `vnmult` |
| `step_end` | `aparnew`, `bparnew`, `g`, `gexp_1`, `gexp_2`, `gexp_3`, `gnew`, `kx_shift`, `phinew`, `theta0_shift`, `vnmult` |
| `after_parity` | `gnew`, `vnmult` |
| `after_exb_shear` | `aparnew`, `bparnew`, `gnew`, `ikx_indexed`, `jump`, `kx_shift`, `phinew`, `theta0_shift`, `vnmult` |
| `coll_h_in` | `h`, `vnmult` |
| `coll_gle_in` | `gle`, `vnmult` |
| `coll_after_ediffuse` | `gle`, `vnmult` |
| `coll_after_ediffuse_conserve` | `gle`, `vnmult` |
| `coll_gle_out` | `gle`, `vnmult` |
| `coll_h_out` | `h`, `vnmult` |
| `after_vspace_derivatives` | `gnew`, `vnmult` |
| `coll_after_drag` | `gle`, `vnmult` |
| `coll_after_lorentz` | `gle`, `vnmult` |
| `coll_after_lorentz_conserve` | `gle`, `vnmult` |
| `after_collisions_advance` | `gnew`, `vnmult` |
| `nl_before_center` | `g1`, `vnmult` |
| `nl_after_center` | `g1`, `vnmult` |
| `nl_cfl` | `max_vel_components`, `vnmult` |
| `after_reset_time_step` | `vnmult` |
| `after_hyper` | `gnew`, `hypervisc_filter`, `vnmult` |
| `after_remove_zonal_flows` | `aparnew`, `bparnew`, `phinew`, `vnmult` |

## 10. Verification actually run (2026-10-21)

All runs single rank (`mpirun -np 1`), no OpenMP, GCC 15.2.0, Open MPI 5.0.9, macOS arm64.

1. **Public-only edits, object code** (stock flags, `otool -tv` vs. the objects of the original build;
   §6 item 1): identical instruction streams apart from line-number immediates for 8 modules; extra
   out-of-line `get_rightmost_it` in `dist_fn.o`; I/O argument scheduling in `theta_grid.o`.
2. **Probe P-R1 rejected**: with P-R1 (either form) the Cyclone linear test output differed from pristine
   at ~1e-16 relative (129 of 227 variables); without it: identical (§3.2).
3. **Numerics unchanged (09 §4.4)**: pristine worktree of `HEAD` vs. a worktree with `HEAD` + all
   changes of this task, both built with the stock command; every variable of `<run>.out.nc` compared
   bitwise for **35 inputs**: `tests/linear_tests/cyclone_itg/cyclone_itg_low_res_base.in` (the required
   case), `fields_local_box.in`, `tests/nonlinear_tests/cyclone_itg/cyclone_itg_low_res.in`,
   `cyclone_itg_collisions_regular_low.in`, `cyclone_itg_collisions_split_low.in`,
   `cyclone_itg_timestep_change.in` (local fields + dt changes), a nonlinear case with hyperviscosity and
   flow shear, and all 28 `gpu/cases`: **all IDENTICAL** (166–306 variables each); the 4 restart files
   written by the dt-change cases are identical too. The original `bin/gs2` (built before any edit, in
   the main tree) and the pristine worktree build also agree bitwise (Cyclone linear test). Both
   worktrees were removed afterwards.
4. **Trajectory with probes active**: in every refdata case the traced `gs2_export` run reproduces the
   time-dependent variables of `bin/gs2`'s `.out.nc` bitwise at all common time records (2 records for
   linear, 3 for nonlinear cases with N = 10), including the dt-reset cases.
5. **Layout independence** (`compare_layouts.py`): exports of `lin_box_ae_linked`,
   `lin_single_coll_ke_full` (2 species, collisions, LE arrays) and `nl_box_ae` with `layout = lxyes`
   vs. `xyles`, `yxles` and `yxels`: all 517 / 553 static per-point arrays and scalars bitwise identical.
   For `xyles`/`yxles` (same relative order of `il, ie, is`) even the velocity-sum-derived arrays are
   identical; for `yxels` (energy before lambda) `gamtot*`, `aminv_c*`, `inv_phi_denominator_g`, the
   conservation vectors and the t = 0 fields/g differ at 2e-16 to 2e-15 relative, as expected from the
   changed summation order (`iglo_of` records it).
6. **Canonical map vs GS2's restart writer** (09 §12 T2, layout `yxels` and every refdata case):
   `iglo_of == iglo_kykxsel`; `gnew`, `vpa`, `vperp2` of the restart file re-indexed with `iglo_of`
   equal the export bitwise.
7. **LE6 map**: `aj0le`, `vperp_aj1le`, `vpa_aj0_le` rebuilt from `aj0`, `aj1`, `vpa` with
   `ixi_to_il/ixi_to_isgn` (forbidden points skipped) are bitwise equal in all collisional cases, and in
   the traces `coll_gle_in` = gather(`coll_h_in`), `coll_h_out` = scatter(`coll_gle_out`) bitwise.
8. **Lossless compression**: an NL trace and export written with `compression = 1` are bitwise identical
   to the uncompressed ones (trace 164 → 84 MB, export 27 → 5 MB).
9. **make_refdata.sh**: `gpu/scripts/make_refdata.sh` (reference build, N = 10, ladder 2, compression 1) on the 28
   `gpu/cases`: **28 of 28 PASS** (60–73 `check_export.py` checks per case, 0 failures; gate_status `ok`
   for all), 334 s, 4.0 GB in `gpu/refdata` (≈300–340 MB per nonlinear 16x16 case, mostly the trace).
   Afterwards `split_collisions` was made unconditional in `/collisions`; the 23 cases with inactive
   collisions were updated one at a time: `lin_box_ae_nonad_zero_f` fully regenerated, the other 22
   (`lin_box_ae_{fieldlineavg_lexys,linked,parity,wfb_passing,wfb_trapped,zero_forbid_yavg}`,
   `lin_range_ae_{exb_theta0,selfperiodic_np2,zero_bc}`, `lin_single_ae`, `lin_single_ae_local`,
   `lin_single_em_ke_miller`, `lin_single_hybrid`, `nl_box_ae`, `nl_box_ae_dt_{decrease_deferred,
   decrease_immediate,increase}`, `nl_box_ae_exb_hyper_iso`, `nl_box_ae_hyper_{aniso_constamp,simple3d}`,
   `nl_box_ae_local`, `nl_box_ae_remove_zonal`) with `GS2_REFDATA_ONLY=export` (their `gs2/` and `trace/`
   data unchanged; all checks against them pass). Final state: 28 of 28 `STATUS` PASS.

## 11. Known gaps

* P-R1 (pre-inversion `am`) not exported (§3.2); multi-rank export not implemented (§3.3).
* LZ6/E7 standard-layout collision matrices not exported (`use_le_layout = .false.` is gated).
* `fields_local`: only the dumped (pulled) response matrices are exported, not the pre-inversion local
  matrices (09 Q5); the response files of an empty supercell may be missing (then listed in `absent`).
* The dt ladder has downward rungs only (`code_dt/delt_adj^k`, 10 G-19); the case `nl_box_ae_dt_increase`
  increases dt above the initial value, which no rung covers (JAX must RECOMPUTE, 00 D3).
* With `immediate_reset = .false.` the CFL limit of step k (`code_dt_cfl`) is completed after the
  `step_end` probe and appears first in `advance_begin(k+1)` (no probe inside `check_time_step`).
* `/final` (trace) is taken after the forced final diagnostics call, so it is not equal to the last
  `step_end` when that call wrote a record.
* Trace size: ~30 MB per nonlinear 16x16 step after compression (≈36 full G7 dumps per step, several of
  them exact duplicates by construction, e.g. `gexp_*` at `advance_end`/`step_end_pre_diag`/`step_end`);
  use `trace_stages`, `trace_first_step`, `trace_every` or `GS2_REFDATA_NSTEP_TRACE` to reduce it.
* `check_export.py` checks shapes for every master-list name but recomputes only the subset listed in §8;
  coefficients such as `wdrift`, `a/b/r/ainv`, `g_h`, `aj0` and the collision matrices are only
  validated through determinism, layout independence and the trajectory comparison.
