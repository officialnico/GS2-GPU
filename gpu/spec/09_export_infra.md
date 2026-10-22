# Spec 09 — Export infrastructure: the `gs2_export` program

Status: contract for the implementers of the Fortran export program and of the Python/JAX loader.
Scope: how GS2 (Fortran, unchanged physics) is driven as the front end; how every array that the other
specs (01–08) mark EXPORT/BOTH/STATE is read out of GS2 memory, re-indexed into the canonical layout and
written to netCDF; how the Fortran time step is instrumented to dump intermediate stage results; which
private module variables need an access path and the exact source edits; the Makefile/dependency
mechanics for adding a program and modules; the GS2 restart and response-matrix file formats (for
round-trips JAX ↔ GS2).

Not in scope: the physics formulas (owned by 01–08). This spec cites them only to say *what* must be
exported and *when* it must be read.

All file:line citations are against commit `4d8c94bcf` (branch `gpu-port`). Line numbers in "insert after
line N" instructions refer to that commit; implementers must re-anchor on the quoted code text.

---------------------------------------------------------------------------------------------------

## 0. Build facts, conventions and terminology

### 0.1 Build facts that matter for the exporter

* `bin/gs2 --build-config`: GCC 15.2.0 (`/opt/homebrew/bin/gfortran`, wrapper `mpif90`, Open MPI 5.0.9),
  `-O3 -mcpu=apple-m1 -fPIC -fdefault-real-8 -fdefault-double-8 -ffree-line-length-none -cpp`,
  macros `-DMPI -DNEW_DIAG -DFFT=_FFTW3_ -DNETCDF -DGK_NETCDF_DEFAULT_COMPRESSION_ON -DF200X_INTRINSICS
  -DSPFUNC=_SPF200X_`; **no** OpenMP, **no** SHMEM, **no** LOWFLOW, **no** parallel netCDF, double precision.
  Hence default `real` = `real(real64)`, default `complex` = `complex(real64)` (16 bytes), default `integer` =
  `int32`, and `cmplx(x, y)` without KIND returns default complex = complex(8).
* `GK_SYSTEM=macosx-homebrew` (from the build config string). There is no `Makefile.local` and `GK_SYSTEM`
  is not in the environment, and the include paths in the build string (`-I/opt/homebrew/include`, Open MPI
  Cellar paths) do not come from `Makefiles/Makefile.macosx-homebrew` (which uses `/usr/local`), so `bin/gs2`
  was built with extra variables on the make command line. The exporter must be built with **exactly the
  same make command line** (Open question Q13).
* All `.f90` and `.fpp` sources are compiled with `-cpp` (F90FLAGS contains `$(FC_PREPROCESS_FLAG)`,
  `Makefile:548`; compile rule `Makefile.utils:64-67`), so `#ifdef` works in `.f90` files too.
* netCDF: neasy-f (`externals/neasyf/src/neasyf.f90`) over netCDF-Fortran 4.6.2. `neasyf_open(...,"w")`
  creates a **NETCDF4/HDF5** file with `NF90_CLOBBER` (`neasyf.in.f90`, function `neasyf_open`).
* Python side: `gpu/.venv` contains `netCDF4 1.7.5`, `numpy 2.5.4`, `f90nml` (checked).

### 0.2 Canonical index conventions (project-wide; restated so this spec is self-contained)

Every array is written to netCDF **in Fortran dimension order** (first index fastest). netCDF/NumPy
(C order) therefore shows the dimensions **reversed**. All index values (`it`, `ik`, `il`, `ie`, `is`,
`isgn`) are expanded from the layout compound indices with the layout's own index functions
(`it_idx/ik_idx/il_idx/ie_idx/is_idx`), independent of the `layout` string.

| tag | Fortran array in GS2 | canonical netCDF dims (Fortran order) | NumPy shape |
|---|---|---|---|
| **G7** | `x(-ntgrid:ntgrid, 2, g_lo)` | `(theta, sign, kx, ky, lambda, egrid, species)` | `(nspec, negrid, nlambda, naky, ntheta0, 2, nth)` |
| **G6** | `x(-ntgrid:ntgrid, g_lo)` | `(theta, kx, ky, lambda, egrid, species)` | `(nspec, negrid, nlambda, naky, ntheta0, nth)` |
| **G7m1** | `x(-ntgrid:ntgrid-1, 2, g_lo)` (`source_coeffs_*`) | `(theta_m1, sign, kx, ky, lambda, egrid, species)` | `(…, 2, 2*ntgrid)` |
| **GL** | `x(n_links_max, 2, g_lo)` (`g_adj`) | `(nlinks_max, sign, kx, ky, lambda, egrid, species)` | `(…, 2, n_links_max)` |
| **G5** | `x(g_lo)` | `(kx, ky, lambda, egrid, species)` | `(nspec, negrid, nlambda, naky, ntheta0)` |
| **G5s** | `x(2, g_lo)` (`save_h`) | `(sign, kx, ky, lambda, egrid, species)` | `(…, ntheta0, 2)` |
| **F3** | `x(-ntgrid:ntgrid, ntheta0, naky)` | `(theta, kx, ky)` | `(naky, ntheta0, nth)` |
| **K2** | `x(ntheta0, naky)` | `(kx, ky)` | `(naky, ntheta0)` |
| **LE6** | `x(nxi_lim or nxi+1, negrid or negrid+1, le_lo)` | `(xi_lim or xi_p1, egrid or egrid_p1, theta, kx, ky, species)` | `(nspec, naky, ntheta0, nth, NE, NX)` (spec 06 §2.4) |
| **LZ6** | `x(nxi_lim, lz_lo)` (`use_le_layout=F` only) | `(xi_lim, theta, kx, ky, egrid, species)` | `(nspec, negrid, naky, ntheta0, nth, nxi_lim)` |
| **E7** | `x(negrid, e_lo)` (`use_le_layout=F` only) | `(egrid, theta, sign, kx, ky, lambda, species)` | `(nspec, nlambda, naky, ntheta0, 2, nth, negrid)` |
| **V** | any velocity-grid / θ / per-species array | Fortran declared order, e.g. `wl(theta, lambda)`, `w(egrid, species)`, `wstar(ky, egrid, species)`, `xi(xi2, theta)` | reversed |

`nth = 2*ntgrid+1`; the netCDF `theta` index 0 corresponds to Fortran `ig = -ntgrid`. Lower bounds other
than 1 are recorded in the variable attribute `fortran_lbounds` (e.g. `"-ntgrid,1,0"`).

**Complex values** use a leading Fortran dimension `ri` (size 2: 1 = real, 2 = imaginary), the convention
of `gs2_io` (`complex_dim = "ri"`, `src/gs2_io.fpp:55`; `final_field_dims = [ri, theta, kx, ky]`,
`gs2_io.fpp:95`) and of `gs2_save_response` (`gs2_save.fpp:938-948`). In NumPy the `ri` axis is the
**last** axis: `z = a[..., 0] + 1j*a[..., 1]`. A G7 complex variable is therefore **8-dimensional** in
netCDF: `(ri, theta, sign, kx, ky, lambda, egrid, species)`.

**Logical values** are written as `int8` (0/1), attribute `fortran_type = "logical"` (neasy-f has no logical
writer; supported types are int8/int16/int32/real32/real64/character, ranks 0–7, `neasyf.in.f90`
`TYPE_NAMES`). **Integers** are `int32` and are stored with their Fortran values (1-based indices stay
1-based, `-1` sentinels stay `-1`); attribute `index_base = 1|0` documents index-valued arrays.
**Strings** are group or global attributes.

### 0.3 Terminology

* *proc0*: rank 0 of `mp_comm` (`externals/utils/mp.fpp:80-90`).
* *glo order*: GS2's native distributed order of a `g_lo` array (`iglo = 0 … ulim_world`).
* *export file*: `<run_name>.gs2export.nc` (static coefficients, grids, configuration, t = 0 state).
* *trace file*: `<run_name>.gs2trace.nc` (Fortran-stepped intermediate stage dumps).
* *probe*: a one-line call inserted in the Fortran time step that forwards a tag (and, where necessary,
  a local array) to a hook installed by `gs2_export` (§5). Probes are inert unless `probe_active`.

---------------------------------------------------------------------------------------------------

## 1. Call graph / control flow

### 1.1 The standard `gs2` program and the `gs2_main` state API

```
src/programs/gs2.fpp:46-58   program gs2
   52  call parse_command_line()          gs2_main.fpp:245-301 (--help/--version/--build-config/--check-input)
   53  call init_mp                        mp.fpp:640
   54  state%mp_comm = mp_comm
   55  if (.not. run_optimise_gs2(state)) stop      gs2_optimisation (optimisation mode)
   56  call run_gs2(state)                 gs2_main.fpp:1614-1661
   57  call finish_mp                      mp.fpp:831

gs2_main.fpp:1614 run_gs2(state, finalize, quiet)
  1640  initialize_gs2(state, header_out=local_header, quiet)      318-463
  1641  initialize_equations(state)                                  467-535
  1642  initialize_diagnostics(state, header=local_header)          542-614
  1644  write_used_inputs_file(header=local_header)                  1703-1737  (<run>.used_inputs.in, proc0)
  1646  state%print_times = .false.
  1647  if (state%do_eigsolve) run_eigensolver(state) else evolve_equations(state, state%nstep)   649-894
  1653  if (proc0) call state%exit_reason%write_exit_file()
  1655  if (state%set_outputs) calculate_outputs(state)
  1657  if (local_finalize) finish_gs2(state) -> finalize_diagnostics 977-1015,
                                              finalize_equations 1018-1043, finalize_gs2 1046-1099
```

`initialize_gs2` (318-463): `init_mp(state%mp_comm)` (358); optional `use_nproc` when an optimisation
override requests fewer ranks (374-380); on proc0 `init_file_utils(state%list, trin_run=…, name=…)`
(425-426) — this takes the input file from **command argument 1** (`file_utils.fpp:130-146`), strips the
`.in` suffix to form `run_name` (`file_utils.fpp:163-174`) and opens the comment-stripped input unit;
broadcast of `list` (430), `job_fork` for list/ensemble runs (433-441), broadcast of `run_name` (446-448);
`init_gs2_init()` (455); `state%init%level = basic` (457).

`initialize_equations` (467-535): `init_parameter_scan` (486); `init(state%init, init_level_list%full)`
(492, walks the 35-level ladder, `src/gs2_init_level_list.inc`, table in spec 01 §1.2);
`init_tstart(tstart)` (496); `state%nstep = nstep`, `state%do_eigsolve = do_eigsolve` (520-521);
`init_gs2_reinit` (531: `save_dt_min(delt_minimum)`, `save_dt_max(dt0)`, `gs2_reinit.f90:257-302`).

`initialize_diagnostics` (542-614): new diagnostics: `init_gs2_diagnostics_new(header)` (595) and
`run_diagnostics(0, state%exit)` (600) — **this performs in-place `g_adjust` round trips on `gnew` and `g`**
(spec 08 §0 item 3, §3.1.8); old diagnostics: `init_gs2_diagnostics`, `loop_diagnostics(0, …)` (604-607).

`evolve_equations(state, nstep_run, header)` (649-894), per step (`istep` from `state%istep_end+1`):

```
730  call advance(istep)                                   fields.f90:586-605
734  if (state%dont_change_timestep) reset = .false.
739  if (reset) then                                       immediate_reset path
741     prepare_initial_values_overrides(state)            1145-1163
743     set_initval_overrides_to_current_vals(...)         1166-1218 (in_memory copy or gs2_save_for_restart)
744     reset_time_step(state%init, istep, state%exit, …)  gs2_reinit.f90:104-225
751     if (.not. state%exit) call advance(istep)          second attempt of the same istep
764  should_save_restart = gnostics%nsave > 0 .and. mod(istep, nsave) == 0  -> gs2_save_for_restart (766-768)
771  call update_time                                      gs2_time.f90:161-174
777/780  loop_diagnostics(istep, exit) / run_diagnostics(istep, exit)   (more g_adjust round trips)
785  if (state%exit) state%converged = .true.
800  check_time_step(reset, exit)                          gs2_reinit.f90:228-252
811  reset_time_step(...)  (no re-advance)                 if reset
829  state%istep_end = istep
834  max_sim_time test
877/880  after the loop: loop/run_diagnostics(istep - offset, exit, force=.true.)
```

`gs2_program_state_type` (gs2_main.fpp:125-236) fields used by the exporter: `init` (`init_type`,
`gs2_init.fpp:50-96`: `level`, `diagnostics_initialized`, overrides), `mp_comm` (must be set before
`initialize_gs2`), `nstep` (upper bound for the cumulative number of steps), `istep_end`,
`skip_diagnostics`, `dont_change_timestep` (keep `.false.` to mirror gs2), `exit`, `exit_reason`,
`print_times`, `included`, `verb`. All data arrays are module globals (comment at gs2_main.fpp:46-53).

### 1.2 One time step (implicit fields, the default) with the probe points of §5

```
gs2_main::evolve_equations
 └ fields::advance(istep)                                fields.f90:586-605           [P advance_begin / advance_end]
    └ fields_implicit::advance_implicit                  fields_implicit.f90:315-371
       335 exb_shear(gnew, phinew, aparnew, bparnew, istep)   (iff kx_shift or theta0_shift allocated) [P after_exb_shear]
       337 copy(gnew, g); 338 copy(phinew,phi) …                                               [probe_pass = 1]
       341 timeadv(phi,…, istep)                          dist_fn.fpp:4014-4116
           4060 add_explicit_terms(gexp_1, gexp_2, gexp_3, phi, apar, bpar, istep, bkdiff(1))
                └ nonlinear_terms::add_explicit           nonlinear_terms.fpp:441-612
                    506 copy(g2,g3); copy(g1,g2)  520 add_nl(g, g1, …) [P nl_before_center]
                    523 nl_center(g1, bd)         [P nl_after_center]
                    567/569 max_allreduce + set_cfl_time_step_and_check_if_too_large  [P nl_cfl]
                                                  [P after_add_explicit_terms]  (before the reset return at 4066)
           4066 if (reset .and. immediate_reset) return
           4068 invert_rhs                                dist_fn.fpp:6152-6193
                6181 get_source_term_on_local_domain(g_work, …)   [P source]
                6184 invert_rhs_1(phinew, bparnew, source=g_work) [P after_invert_rhs_1]
                6190 apply_linked_boundary_conditions(…)          (linked only)
                                                  [P after_invert_rhs]
           4077 hyper_diff(gnew, phinew)          [P after_hyper]
           4084 vspace_derivatives(gnew, g, g_work, phi, bpar, phinew, bparnew, mode)   le_derivatives.fpp:15-…
                 71 g_adjust(g → h)   [P coll_h_in]   83 gather(g2le, g, gle) [P coll_gle_in]
                 90 solfp1(gle)  → collisions::solfp1_le_layout  collisions.fpp:3035-3090
                       3051 solfp_ediffuse_le_layout [P coll_after_ediffuse]  3052 conserve_diffuse [P coll_after_ediffuse_conserve]
                       3056-3085 drag [P coll_after_drag]  3087 solfp_lorentz [P coll_after_lorentz]
                       3088 conserve_lorentz [P coll_after_lorentz_conserve]
                 [P coll_gle_out]  94 scatter(g2le, gle, g) [P coll_h_out]  109 g_adjust(h → g)
                                                  [P after_vspace_derivatives]
           4091 enforce_parity                    [P after_parity]
                                                  [P after_timeadv]
       342 if (reset) return
       352 getfield(phinew, aparnew, bparnew)     [P after_getfield]  (phinew holds the increment u)
       354-356 phinew = phinew + phi …            [P after_field_update]
       362 remove_zonal_flows                     [P after_remove_zonal_flows]                  [probe_pass = 2]
       365 timeadv(…, istep, diagnostics)         (same inner probes; add_explicit returns early because istep == istep_last)
                                                  [P after_timeadv]
       369 collisions_advance(…)                  [P after_collisions_advance]
 ├ (reset path) reset_time_step                    [P after_reset_time_step]   then advance(istep) again (attempt 2)
 ├ update_time                                     [P step_end_pre_diag]
 └ run_diagnostics(istep)                          [P step_end]
```

`fields_local::advance_local` (fields_local.f90:4244-4342) has the same structure (exb 4280, copies
4283-4290, timeadv 4296-4297, getfield_local 4313, adds 4316-4320, remove_zonal_flows 4331, timeadv
4336-4337, collisions_advance 4340) and receives the same probes.

During response-matrix construction (`init_response_matrix`, fields_implicit.f90:459-647) `timeadv` is
called with `istep = 0` outside `fields::advance` (fields_implicit.f90:669); probes there are filtered out
by `probe_in_advance = .false.` (§5.3), except the optional `response_am` probe (§5.2, P-R1).

### 1.3 The existing dump programs (models)

`src/programs/dump_response.f90:15-192`: `init_mp` (29); proc0 `init_file_utils(list, name="gs")` (45);
broadcast `list`, `job_fork` if list (48-51); broadcast `run_name` (53-55); proc0 reads its own namelist
from the GS2 input file with `input_unit_exist("dump_response_knobs", exist)` (57-61) — the pattern we
reuse for `&gs2_export_knobs`; then it drives the **old module-level init chain** directly
(`fields_pre_init`, `set_dump_and_read_response(.false.,.false.)`, `fields_init_response`,
`dump_response_to_file`, 104-159) and steps through a dt sequence with `reduce_time_step`, `save_dt`,
`finish_dist_fn_level_3`, collisions `reset_init`, fields `reset_init` (83-101). It finishes the modules
one by one (162-191). It is listed in `parallel_targets` and `needs_git_hash` (Makefile.target_gs2:29,37).

`src/programs/dump_grids.fpp:3-309`: `init_mp`, `init_file_utils`, `init_theta_grid`, `init_kt_grids`,
`init_species`, `init_le_grids`, then proc0 writes `<run>.grids.nc` with neasy-f: `neasyf_dim` with
coordinate values for `kx, ky, theta, lambda`, `neasyf_dim(…, dim_size=…)` for `species, egrid`,
`neasyf_write(…, dim_names=[kx_dim, ky_dim])` for `theta0`, `gs2_io::nc_geo` and `nc_species`
(110-161). It sets `neasyf_default_compression` from `GK_NETCDF_DEFAULT_COMPRESSION_ON` (128-132).
It is a *serial* target (Makefile.target_gs2:4).

**Decision:** `gs2_export` does **not** copy dump_response's module-level init chain. It drives GS2 through
the `gs2_main` state API exactly like `run_gs2` (§1.4), so the initialisation path (init-level ladder,
overrides, `init_tstart`, `init_gs2_reinit`, diagnostics at istep 0) is byte-for-byte the one used by
`bin/gs2`. It reuses dump_response's namelist-reading pattern and dump_grids' neasy-f idioms.

### 1.4 `gs2_export` control flow (contract)

```
program gs2_export                                       src/programs/gs2_export.fpp
  call init_mp ; state%mp_comm = mp_comm                 (as gs2.fpp:53-54; no parse_command_line, no optimisation)
  call initialize_gs2(state, header_out=header)          proc0: init_file_utils reads argv(1)
  call read_export_knobs(knobs)                          §8.1: argv(2) file if present, else &gs2_export_knobs in the input file
  call check_hard_gates(knobs)                           §8.3 (abort before any heavy work)
  if (knobs%trace_response_build) then                   capture pre-inversion `am` during init
      call install_hooks(knobs, response_only=.true.); probe_active = .true.
  call initialize_equations(state)                       init(full) + init_tstart + init_gs2_reinit
  probe_active = .false.
  call collect_soft_gates(gates)                         §8.3 (recorded, never abort)
  call ex_open(export_file, header)                      proc0: neasyf_open("w"), create_metadata, root dims §6.1.1
  call export_static(knobs)                              §6.1.3 (collective; every rank calls)
  call export_state('state_t0_prediag')                  §7.6
  if (knobs%run_diagnostics) then
      call initialize_diagnostics(state, header=header)  istep-0 diagnostics (g_adjust round trips on gnew, g)
      call write_used_inputs_file(header=header)
  else
      state%skip_diagnostics = .true. ; gnostics%nsave = -1     (NEW_DIAG; see Pitfall 11.9)
  end if
  call export_state('state_t0')                          THE initial state of JAX step 1
  call export_config_text()                              used_inputs text, input file text
  if (knobs%write_restart_t0) call gs2_save_for_restart(gnew, user_time, vnmult, has_phi, has_apar,
        has_bpar, code_dt, code_dt_prev1, code_dt_prev2, code_dt_max, fileopt='.gs2export_t0',
        save_glo_info_and_grids=.true., save_velocities=.true., header=header)
  if (knobs%dt_ladder_n > 0) call export_dt_ladder(state, knobs)     §7.8 (re-initialises; no stepping afterwards)
  call ex_close()
  if (knobs%nstep_trace > 0) then                        mutually exclusive with dt_ladder_n > 0
      call tr_open(trace_file); call install_hooks(knobs); probe_active = .true.
      call evolve_equations(state, min(knobs%nstep_trace, state%nstep), header=header)
      probe_active = .false.
      call tr_snapshot('final')                          state after the forced final diagnostics call
      call tr_close()
  end if
  if (state%init%diagnostics_initialized) call finalize_diagnostics(state)
  call finalize_equations(state); call finalize_gs2(state)
  call finish_mp
end program gs2_export
```

`gs2_export input.in` with `nstep_trace = 0` and `run_diagnostics = .true.` also leaves `<run>.out.nc`
containing GS2's own static part and the istep-0 record (spec 08 §0 item 7 recommends the JAX side reuse
this file for the static part of its output).

### 1.5 Restart and response I/O (gs2_save.fpp) call graph

```
gs2_save_for_restart(g, t0, vnm, has_phi, has_apar, has_bpar, code_dt, code_dt_prev1, code_dt_prev2,
                     code_dt_max, fileopt, save_glo_info_and_grids, save_velocities, header)   86-469
  is_one_file_per_processor = save_many .or. .not. has_netcdf_parallel        (175; true in this build)
  file = get_file_proc(...)  = trim(restart_file)//fileopt//'.'//iproc          1248-1257
  neasyf_open("w"); create_metadata("GS2 restart file")                          185-195
  dims: glo (local count, unlimited if 0) 210, theta (values theta) 212, sign 213, aky (values) 214, akx (values) 215
        [lambda (values al), species, egrid, energy(egrid,species)] 220-230 if save_glo_info_and_grids
        nk_stir if ant_on 234-236
  scalars delt0, delt1, delt2, delt_max, t0, layout(char), ntheta0, naky, nlambda, negrid, nspec,
          current_scan_parameter_value, vnm1, vnm2                                246-265
  gr, gi       real (theta, sign, glo)                                           274-299
  gexp_1_r/_i, gexp_2_r/_i, gexp_3_r/_i  (theta, sign, glo), if allocated and include_explicit_source_in_restart  301-341
  iglo_llim_proc, iglo_kykxsel(aky,akx,species,egrid,lambda), procid_kykxsel(...), jtwist  343-381
  vpa(theta,sign,glo), vperp2(theta,glo) if save_velocities                     383-388
  phi_r/_i, apar_r/_i, bpar_r/_i (theta, akx, aky) from phinew/aparnew/bparnew   390-428
  kx_shift, theta0_shift (aky) (zeros if unallocated)                            430-453
gs2_restore_many(g, scale, istatus, has_phi, has_apar, has_bpar, fileopt)       604-888
  layout check: abort if restart layout /= current layout                        767-771
  g = cmplx(gr, gi); gexp_* likewise; kx_shift/theta0_shift if allocated; phinew… ; phi,apar,bpar zeroed
  scale: g = g*scale … if scale > 0, else fac = -scale/maxval(abs(phinew))       863-875
init_dt(delt0, delt1, delt2, delt_max, istatus)    1296-1386 (proc0 reads, broadcast)
init_vnm(vnm, istatus)                              1389-1442
read_t0_from_restart_file(tstart, istatus)          1517-1560
gs2_save_response(resp, fname, code_dt, condition, header)   892-969
  dims ri(2), ax1, ax2; vars dt, condition (default -1.0), response(ri, ax1, ax2) via convert::c2r
gs2_restore_response(resp, fname, code_dt, condition)        973-1035
```

---------------------------------------------------------------------------------------------------

## 2. Layouts: exact index arithmetic

### 2.1 `g_lo` (`gs2_layouts.fpp:826-1528`, type `layouts_type.f90:33-69`)

Construction (non-SHMEM branch, 861-865), with `N_g = naky*ntheta0*negrid*nlambda*nspec`:

```
llim_world = 0 ; ulim_world = N_g - 1
blocksize  = ulim_world/nproc + 1               (integer division; nproc = nproc_g_lo, see below)
llim_proc  = blocksize*iproc
ulim_proc  = min(ulim_world, llim_proc + blocksize - 1)
ulim_alloc = max(llim_proc, ulim_proc)
proc_id(i) = i / blocksize                       (1583-1603)
```

`nproc_g_lo = handle_processor_request(nproc_g_lo, nproc)` (1538, 479-487): input `layouts_knobs:nproc_g_lo`
(default 0 → all ranks). Ranks with `llim_proc > ulim_world` own nothing (`ulim_proc < llim_proc`) but
still allocate **one** slot (`ulim_alloc = llim_proc`); loops `do iglo = llim_proc, ulim_proc` are empty
there. Gathers must use `ulim_proc`, never `ulim_alloc` (Pitfall 11.6).

Dimension order from the layout string (`select case (layout)`, 917-930; `set_dimension_order_g`
4765-4781). Position 1 is the fastest-varying index of `iglo`:

| `layout` | it_ord | ik_ord | il_ord | ie_ord | is_ord | order (fast → slow) |
|---|---|---|---|---|---|---|
| `'lxyes'` (default, gs2_layouts.fpp:204) | 2 | 3 | 1 | 4 | 5 | il, it, ik, ie, is |
| `'xyles'` | 1 | 2 | 3 | 4 | 5 | it, ik, il, ie, is |
| `'yxles'` | 2 | 1 | 3 | 4 | 5 | ik, it, il, ie, is |
| `'yxels'` | 2 | 1 | 4 | 3 | 5 | ik, it, ie, il, is |
| `'lexys'` | 3 | 4 | 1 | 2 | 5 | il, ie, it, ik, is |
| `'lyxes'` | 3 | 2 | 1 | 4 | 5 | il, ik, it, ie, is |

Any other string aborts in `read_parameters` (566-571). `dim_size(*_ord) = (ntheta0, naky, nlambda, negrid,
nspec)` (932-936); `compound_count(1)=1, compound_count(k)=compound_count(k-1)*dim_size(k-1)` (938-942);
`x_comp = compound_count(x_ord)` (944-948). Index functions (1543-1580, `llim_world = 0`):

```
is_idx(i) = 1 + mod(i/is_comp, nspec)      il_idx(i) = 1 + mod(i/il_comp, nlambda)
ie_idx(i) = 1 + mod(i/ie_comp, negrid)     it_idx(i) = 1 + mod(i/it_comp, ntheta0)
ik_idx(i) = 1 + mod(i/ik_comp, naky)
idx(ik,it,il,ie,is) = i(1) + dim_size(1)*(i(2) + dim_size(2)*(i(3) + dim_size(3)*(i(4) + dim_size(4)*i(5))))
     with i(x_ord) = x - 1                                                                     (1606-1617)
```

The exporter writes `/gs2_layouts/iglo_of` (G5 int32, `index_base = 0`) = `idx(g_lo, ik, it, il, ie, is)`
and the five `*_ord` integers, so the JAX side can reproduce **velocity-space summation order** (sums run
over `iglo` ascending on one rank; spec 04 §3.1, spec 08 §0 item 5).

### 2.2 `le_lo`, `lz_lo`, `e_lo` (collision layouts)

`le_lo` (2003-2114): `ulim_world = nth*naky*ntheta0*nspec - 1`; same block formula; orders (2030-2043):
`ig_ord = 1` always, `is_ord = 4` always; `(it_ord, ik_ord) = (2,3)` for `'lexys','lxyes','xyles'` and
`(3,2)` for `'yxels','yxles','lyxes'`. `ig_idx_le(i) = -ntgrid + mod(i/ig_comp, nth)` (2152-2157), the other
functions as for g_lo; `idx_le(ig,ik,it,is)` (2160-2170) with `i(ig_ord) = ig + ntgrid`.
Canonical LE6 expansion: for each `ile`, `(ig,it,ik,is) = (ig_idx, it_idx, ik_idx, is_idx)(le_lo, ile)`.

`lz_lo` (2602-2691): orders (2655-2668) `ig_ord=1, is_ord=5`; `'lexys'`: it 3, ik 4, ie 2; `'lxyes','xyles'`:
it 2, ik 3, ie 4; `'yxels','yxles','lyxes'`: it 3, ik 2, ie 4. Its proc block is **not** the simple formula
when `lambda_local` (2629-2652, `check_llocal` 705-808); the exporter therefore never recomputes block
bounds for other ranks (§2.5). `e_lo` (2438-2497): orders (2462-2475), 6 dims (ig, isgn, it, ik, il, is).
Both are used only with `collisions_knobs:use_le_layout = .false.` (spec 06 §4.9).

### 2.3 `f_lo`, `jf_lo` and the storage of `aminv` (fields_implicit)

`nfield = count(has_phi, has_apar, has_bpar)`, `nidx = nth*nfield` (fields_implicit.f90:477-481).
`f_lo(c)`, `c = 1..i_class` (gs2_layouts.fpp:1650-1691): `nidx`, `nindex = nidx*N_class(c)`,
`ulim_world = nindex*M_class(c) - 1`; index functions 1755-1848: `if_idx = 1+mod(i, nindex)`,
`im_idx = 1+mod(i/nindex, M)`, `in_idx = 1+mod(i/nidx, N)`, `ig_idx = -ntgrid+mod(i/nfield, nth)`,
`ifield_idx = 1+mod(i, nfield)` (field fastest: "K" order of spec 04), `ij_idx_f(ig,if,n) = ig+ntgrid +
nth*(if-1 + nfield*(n-1))` ("B" order). Member tables `f_lo(c)%ik(m,n)`, `f_lo(c)%it(m,n)` are filled by
`fields_implicit::finish_fields_layouts` (fields_implicit.f90:1074-1133; ordering rules in spec 04 §3.11).

`jf_lo` (1855-1898): `ulim_world = nidx*ntheta0*naky - 1`, `ij_idx_jf(ig,if,ik,it) = ig+ntgrid +
nth*(if-1 + nfield*(it-1 + ntheta0*(ik-1)))` (1967-1972); `jf_lo%ij(jflo)` = class, `jf_lo%mj(jflo)` =
member, `jf_lo%dj(c, jflo)` = local dcell index (set in `init_inverse_matrix`, fields_implicit.f90:920-958).

`aminv(c)%dcell(dc)%supercell(1:nidx*N_class(c))` (types fields_implicit.f90:24-47, **private**) holds
one inverse **row** for output `jflo`; the row is in B order of the inputs, block `n` (columns
`1+nidx*(n-1) : nidx*n`) multiplying `fl(:, f_lo(c)%it(m,n), ik)` (getfield, 222-312;
fill 1036-1063: `jsc = ij_idx(f_lo(c), jg, jf, nn) + 1`).

### 2.4 Canonical expansion pseudo-code (exact; used by every writer)

All sizes and flat offsets are `integer(int64)` (Pitfall 11.12). `nth = 2*ntgrid+1`.

```
# G7, complex, written per species slab (bounded memory); full(:,:,0:ulim_world) is the glo-ordered array
# (the local array on 1 rank, or the proc0 gather of §2.5)
SLAB = 2*nth*2*ntheta0*naky*nlambda*negrid
for is = 1..nspec:
    slab(1:SLAB) = 0
    for ie = 1..negrid: for il = 1..nlambda: for ik = 1..naky: for it = 1..ntheta0:
        iglo = idx(g_lo, ik, it, il, ie, is)
        base = 2*nth*2*((it-1) + ntheta0*((ik-1) + naky*((il-1) + nlambda*(ie-1))))
        for isgn = 1..2: for ig = -ntgrid..ntgrid:
            k = 1 + base + 2*((ig+ntgrid) + nth*(isgn-1))
            slab(k)   = real (full(ig, isgn, iglo))
            slab(k+1) = aimag(full(ig, isgn, iglo))
    neasyf_write(gid, name, slab, dim_names=[ri,theta,sign,kx,ky,lambda,egrid,species],
                 start=[1,1,1,1,1,1,1,is], count=[2,nth,2,ntheta0,naky,nlambda,negrid,1])
```

The first call defines the 8-D variable from `dim_names`; later calls only write (neasy-f looks the
variable up first). Writing a rank-1 buffer into an 8-D variable with an explicit `count` was verified on
this machine with `objs/neasyf.o` + netCDF-Fortran 4.6.2 and read back with netCDF4-python: element
order and the `[..., 0] + 1j*[..., 1]` reconstruction are correct (Appendix A). **`count` must always be
passed** (without it netCDF-Fortran's 1-D wrapper uses `count = [size(values), 1, 1, …]`).

The other shapes are the obvious specialisations:

* **G7 real** (`vpa, vpac, vpar, wdrift, wdriftttp`): as G7 without the `ri` factor and `real()`.
* **G6** (`aj0, aj1, vperp2`): `k = 1 + (ig+ntgrid) + nth*((it-1) + ntheta0*((ik-1) + naky*((il-1) + nlambda*(ie-1))))` per species slab.
* **G7m1**: `theta_m1` axis `ig = -ntgrid..ntgrid-1`, length `2*ntgrid`.
* **GL** (`g_adj`): first axis `i = 1..n_links_max = size(g_adj,1)`.
* **G5 / G5s**: integer/int8 without the θ axis.
* **F3** (fields, `kperp2`, `gamtot`, …): already canonical in memory; complex ones get `ri` via
  `convert::c2r` (`externals/utils/convert.f90:16-30`, ranks 1–5): `a_ri(1,…) = real(a)`, `a_ri(2,…) = aimag(a)`.
* **LE6**: for `ile = 0..le_lo%ulim_world`: `(ig,it,ik,is)` from the le_lo functions;
  `out(ri, ixi, ie, ig, it, ik, is) = x(ixi, ie, ile)`, slab over `is`; 6 canonical dims, complex arrays (`gle`) add `ri` → 7-D.
* **LZ6 / E7**: same pattern with lz_lo / e_lo index functions (§2.2).

### 2.5 Gathering for nproc > 1 (collective; exact copies, no arithmetic)

Every writer is called by **all** ranks in the same order; only proc0 touches netCDF.

```
subroutine gather_glo(arr, full)            # arr(-ntgrid:ntgrid, 2, g_lo%llim_proc:g_lo%ulim_alloc)
  if nproc == 1:  full => arr(:,:,0:g_lo%ulim_world)          (no copy)
  elif proc0:
     full(:,:,llim_proc:ulim_proc) = arr(:,:,llim_proc:ulim_proc)
     for ip = 1..nproc-1:                                      (fixed rank order)
        receive(bounds(1:2), ip, tag=T1)                       # mp::receive_integer_array
        if bounds(2) >= bounds(1):
           receive(buf(1:nth*2*(bounds(2)-bounds(1)+1)), ip, tag=T2)    # receive_complex_array (1-D)
           full(:,:,bounds(1):bounds(2)) = reshape(buf, [nth, 2, bounds(2)-bounds(1)+1])
  else:
     send([g_lo%llim_proc, g_lo%ulim_proc], 0, tag=T1)          # mp::send_integer_array, mp.fpp:4132
     if ulim_proc >= llim_proc: send(reshape(arr(:,:,llim_proc:ulim_proc), [n]), 0, tag=T2)   # 4270
```

Each rank sends its *own* bounds, so the same routine works for `le_lo`, `lz_lo` (incl. `lambda_local`),
`e_lo` and `jf_lo`. Do **not** gather with `sum_reduce` of zero-padded buffers: it works numerically
(`x + 0 = x`) but costs nproc× memory and turns `-0.0` into `+0.0`. Field arrays (`phinew` etc.) and
other replicated arrays are taken from proc0 directly (they are identical on all ranks after the
`sum_allreduce`/`allgatherv` of the field solve).

**Rank-dependence warning.** Values that GS2 computes through MPI reductions (velocity-space integrals:
`gamtot*`, `inv_*_denominator_g`, `awgt`, the response matrix `aminv`, collision conservation vectors, every
field solve) depend on `nproc` at round-off level (order of `MPI_Allreduce`, `mp.fpp:3002-3018`). Reference
data for round-off comparison must come from **nproc = 1** (§8.3 hard gate unless `allow_multi_rank`).

---------------------------------------------------------------------------------------------------

## 3. Inventory: everything the exporter reads

Legend. **Access**: `pub` = public module variable today; `PRIV→§4` = private, made readable by the
edit of §4.2; `derived` = not stored in GS2, computed by the exporter with a *verbatim copy* of the GS2
expression or a public pure function (§7.7); `local` = local/SAVE variable, not exportable without
moving code (RECOMPUTE only). **When**: `I` = after `initialize_equations`, `D` = after
`initialize_diagnostics` (or the skip branch), `S` = per probe during stepping. **Decision** repeats the
owning spec's EXPORT/RECOMPUTE/BOTH choice (the exporter writes everything marked EXPORT or BOTH, and
the STATE arrays). **dt**: depends on `code_dt` (exported at the initial dt; also per ladder rung §7.8).
Shape tags as in §0.2. Group = netCDF group = Fortran module name.

### 3.1 `/theta_grid` (theta_grid.f90; geometry, finish_init 401-521)

| name | decl | computed at | dt | access | shape tag | decision (owner) |
|---|---|---|---|---|---|---|
| `theta, theta2, delthet, delthet2` | `(-ntgrid:ntgrid)` | finish_init; `delthet(ntgrid)=0` (510-511) | no | pub | V `(theta)` | EXPORT (01,02,08) |
| `bset` | `(nbset)` | geometry | no | pub | V `(bset)` | EXPORT (02/03 checks) |
| `bmag, gradpar, itor_over_B, IoB` | `(-ntgrid:ntgrid)` | geometry; `itor_over_B` overridden by `btor_slab` (dist_fn.fpp:1369) | no | pub | V | EXPORT effective (02) |
| `gbdrift, gbdrift0, cvdrift, cvdrift0, cdrift, cdrift0` | same | geometry, knobs applied (623-626) | no | pub | V | EXPORT effective (02,08) |
| `gds2, gds21, gds22, gds23, gds24, gds24_noq, grho, jacob` | same | `jacob` 516 | no | pub | V | EXPORT (04,07,08) |
| `Rplot, Zplot, aplot, Rprime, Zprime, aprime, Bpol` | same | geometry | no | pub | V | EXPORT (08 static output) |
| `bmin, bmax, eps_trapped, shat, drhodpsi, kxfac, qval, surfarea, dvdrhon, rhoc` | scalars | `kxfac` 592-619 | no | pub | scalar | EXPORT |
| `ntheta, ntgrid, nperiod, nbset, gb_to_cv` | int/logical | `ntgrid` 567 | no | pub | scalar | EXPORT |
| `shape` | char(8) | | no | pub | attr | EXPORT |
| `field_line_average_weight` | real | 520 `1.0/trapezoidal_integration(theta, jacob)` | no | PRIV→§4 | scalar | BOTH (04, 07, 08) |
| `cvdriftknob, gbdriftknob, eqopt_switch` | real/int | read_parameters | no | PRIV→§4 | scalar | EXPORT (metadata) |
| `is_effectively_zero_shear()` | function | | no | pub fn | scalar int8 | EXPORT (07: exb branch) |

### 3.2 `/kt_grids` (kt_grids.f90, kt_grids_box.f90)

| name | decl | computed at | dt | access | tag | decision |
|---|---|---|---|---|---|---|
| `aky` | `(naky)` | kt_grids_box.f90:378-380 / range 299-311 | no | pub | V `(ky)` | EXPORT (all) |
| `akx`, `ikx` | `(ntheta0)` | kt_grids_box.f90:382-389 | no | pub | V `(kx)` | EXPORT |
| `theta0` | `(ntheta0, naky)` | kt_grids_box.f90:391-401 | no | pub | K2 | EXPORT |
| `kperp2, inv_kperp2` | `(-ntgrid:ntgrid, ntheta0, naky)` | init_kperp2 313-349 | no | pub | F3 real | EXPORT |
| `kwork_filter` | logical `(ntheta0, naky)` | alloc `.false.` 188-189; `(1,1)=.true.` for box at fields.f90:376 (after the response build) | changes during init | pub | K2 int8 | EXPORT value after `init_fields` (01,02,03,04) |
| `naky, ntheta0, nx, ny, jtwist_out, gridopt_switch` | int | get_sizes 269-291 | no | pub | scalar | EXPORT |
| `reality, box` | logical | | no | pub | scalar int8 | EXPORT |
| `grid_option` | char(20) | | no | pub | attr | EXPORT |

### 3.3 `/le_grids` (le_grids.f90)

| name | decl | computed at | dt | access | tag | decision |
|---|---|---|---|---|---|---|
| `energy, w, speed` | `(negrid, nspec)` | set_vgrid 285-296 | no | pub | V `(egrid, species)` | EXPORT |
| `energy_maxw, w_maxw, speed_maxw` | `(negrid)` | 311-313, 1876 | no | pub | V `(egrid)` | EXPORT (06) |
| `al` | `(nlambda)` | lgridset 2120-2156 | no | pub | V `(lambda)` | EXPORT |
| `wl` | `(-ntgrid:ntgrid, nlambda)` | lgridset (+2160-2641) | no | pub | V `(theta, lambda)` | EXPORT |
| `jend` | int `(-ntgrid:ntgrid)` | calculate_jend 2871-2903 | no | pub | V int32 | EXPORT |
| `forbid, is_ttp, is_bounce_point` | logical `(-ntgrid:ntgrid, nlambda)` | 2715-2738, 2907-2946, 2746-2771 | no | pub | V int8 | EXPORT |
| `can_be_ttp` | logical `(nlambda)` | 2907-2946 | no | pub | V int8 | EXPORT |
| `xi, wxi` | `(2*nlambda, -ntgrid:ntgrid)` | 2955-3020 | no | pub | V `(xi2, theta)` | EXPORT (06) |
| `ixi_to_il, ixi_to_isgn` | int same | 2972-3011 | no | pub | V int32 | EXPORT (06) |
| `sgn` | int(2) | | no | pub | V | EXPORT |
| `negrid, nlambda, ng2, lmax, nxi, nesub, wdim` | int | 402-404, 1896-1904 | no | pub | scalar | EXPORT |
| `vcut` | real | | no | pub | scalar | EXPORT |
| `new_trap_int, passing_wfb, trapped_wfb, mixed_wfb` | logical | 413-427 | no | pub | scalar int8 | EXPORT (03) |
| `bouncefuzz, trapped_particles, wfbbc_option_switch, ngauss, npassing, nesuper, genquad, radau_gauss_grid, split_passing_region, new_trap_int_split` | | read_parameters | no | PRIV→§4 | scalar | EXPORT (metadata/gates) |
| `is_lower_bounce_point(ig,il)`, `is_upper_bounce_point(ig,il)` | elemental fns 2775-2819, 2823-2867 | | no | derived (pub fn) | V `(theta, lambda)` int8 | EXPORT-for-check (03) |
| `il_is_wfb, il_is_passing, il_is_trapped` | elemental fns 450-491 | | no | derived | V `(lambda)` int8 | EXPORT-for-check |
| `is_passing_hybrid_electron(is,ik,il)` | fn 4490-4499 | | no | derived | V `(lambda, ky, species)` int8 | EXPORT-for-check (03) |
| `grid_has_trapped_particles()` | fn 443-446 | | no | derived | scalar int8 | EXPORT |

### 3.4 `/species` (species.f90)

| name | decl | computed at | dt | access | tag | decision |
|---|---|---|---|---|---|---|
| `spec(:)%{z, mass, dens, temp, tprim, fprim, vnewk, stm, zstm, tz, smz, zt, uprim, uprim2, dens0, u0, tpar0, tperp0, nu_h, bess_fac, vcrit, vcprim}` | real per species (type 18-35) | 448-480, overrides 991-1013 | no | pub | V `(species)` one variable per component, named `spec_<comp>` | EXPORT (all) |
| `spec(:)%{type, f0type}` | int | | no | pub | V int32 | EXPORT |
| `nspec` | int | | no | pub | scalar | EXPORT |
| `nonmaxw_corr, f0_values, dlnf0drho` | `(negrid, nspec)` | calculate_f0_arrays 563-597 | no | pub | V `(egrid, species)` | EXPORT (02,03,06,08) |
| `has_electron_species, has_ion_species, has_hybrid_electron_species` | fns | 914-944 | no | derived | scalar int8 | EXPORT (04) |

### 3.5 `/run_parameters`, `/gs2_time`, `/gs2_reinit`, `/init_g`, `/fields`

| name | module | decl | computed at | dt | access | decision |
|---|---|---|---|---|---|---|
| `beta, zeff, tite, fphi, fapar, fbpar, k0, rhostar, delt, max_sim_time, avail_cpu_time, margin_cpu_time` | run_parameters | real | read_parameters 300-366 (`beta==0 ⇒ fapar=fbpar=0` 340-349) | no | pub | EXPORT effective |
| `has_phi, has_apar, has_bpar` | run_parameters | logical | 351-353 `f* > epsilon(0.0)` | no | pub | EXPORT (B) |
| `nstep, ncheck_stop, progress_frequency, delt_option_switch` | run_parameters | int | | no | pub | EXPORT |
| `immediate_reset, wstar_units, neo_test, do_eigsolve, use_old_diagnostics, reset` | run_parameters | logical | | no / runtime | pub | EXPORT (reset: state) |
| `ieqzip` | run_parameters | logical `(ntheta0, naky)` | init_run_parameters 259-272 | no | pub | EXPORT K2 int8 (gate) |
| `eqzip_option_switch, seed` | run_parameters | int | | no | PRIV→§4 | EXPORT |
| `code_dt, user_dt, code_dt_prev1, code_dt_prev2, code_dt_old, code_dt_cfl, user_dt_cfl, code_dt_min, code_dt_max, code_time, user_time` | gs2_time | real (19-35) | init_delt 140-151, update_time 161-174, init_tstart 129-136, save_dt_* | state | pub | STATE "time block" (§7.6) |
| `dt_not_set` | gs2_time | parameter −1.0 (19) | | | pub | EXPORT attr |
| `tunits, wunits, woutunits` | gs2_time | `(naky)` | adjust_time_norm 83-126 | no | pub | EXPORT (B) |
| `delt_adj, in_memory` | gs2_reinit | real/logical | init_reinit 257-302 | no | pub | EXPORT |
| `delt_cushion, delt_minimum, dt0, abort_rapid_time_step_change` | gs2_reinit | real/logical (16-20) | init_reinit | no | via `get_reinit_config()` (public, `gs2_reinit.f90:14`; `dt0` smart default resolved at 284-286 before `init`) | EXPORT |
| `tstart, new_field_init, initial_condition_is_nonadiabatic_dfn, restart_file` | init_g | | init_init_g | no | pub | EXPORT |
| `ginit_option, scale, restart_dir, read_many …` | init_g config | | | | via `get_init_g_config()` | EXPORT (attrs) |
| `fieldopt_switch, remove_zonal_flows_switch, force_maxwell_reinit` | fields | int/logical | read_parameters | no | pub | EXPORT |
| `field_option, response_dir, response_file, dump_response, read_response, field_subgath, …` | fields config | | | | via `get_fields_config()`; `fields_arrays::response_file` pub | EXPORT (attrs) |

### 3.6 `/gs2_layouts` and `/layout_maps` (derived)

| name | decl | access | decision |
|---|---|---|---|
| `layout` | char(5) | pub | EXPORT (global attr + group attr) |
| `g_lo%{ik_ord,it_ord,il_ord,ie_ord,is_ord, ulim_world, blocksize, nproc}` | int | pub (type components) | EXPORT |
| `le_lo%{ig_ord,it_ord,ik_ord,is_ord, ulim_world}` | int | pub | EXPORT |
| `iglo_of` (G5, `index_base=0`) | derived `idx(g_lo, ik,it,il,ie,is)` | derived | EXPORT (summation order) |
| `ile_of` `(theta, kx, ky, species)` (`index_base=0`) | derived `idx(le_lo, ig, ik, it, is)` | derived | EXPORT (06 summation order) |
| `local_field_solve, intmom_sub, intspec_sub, gf_local_fields, opt_local_copy, fft_use_wisdom, fft_measure_plan, fft_wisdom_file, simple_gf_decomposition` | | pub | EXPORT (metadata; see Pitfall 11.14 on wisdom) |
| `nproc_*_lo, unbalanced_*, max_unbalanced_*, opt_redist_*` | | via `get_layouts_config()` | EXPORT (attrs) |

### 3.7 `/dist_fn_arrays` (dist_fn_arrays.fpp:26-65; all public)

| name | decl | computed at | dt | tag | decision |
|---|---|---|---|---|---|
| `vpa, vpac` | real `(-ntgrid:ntgrid, 2, g_lo)` | init_vpa_vpac dist_fn.fpp:1592-1678; `vpac(ntgrid,:,:)=0` | no | G7 real | EXPORT |
| `vperp2, aj0, aj1` | real `(-ntgrid:ntgrid, g_lo)` | init_vperp2 1560-1589, init_bessel 1756-1823 | no | G6 | EXPORT |
| `wstar` | real `(naky, negrid, nspec)` | init_wstar 1732-1750 | **dt** | V `(ky, egrid, species)` | BOTH (02) |
| `kx_shift, theta0_shift` | real `(naky)`, allocated iff flow shear (3983-3990) | allocate_arrays, exb_shear 4301-4309, restart | state | V `(ky)` | STATE |
| `g, gnew` | complex `(-ntgrid:ntgrid, 2, g_lo)` | ginit / set_init_fields / stepping | state | G7 | STATE |
| `gexp_1, gexp_2, gexp_3` | same, allocated iff `nonlin .and. .not. split_nonlinear` (3968-3980) | add_explicit / restart | state | G7 | STATE |
| `g_work` | same | scratch | per call | G7 | trace only (`source` probe) |
| `antot, antota, antotp, fieldeq, fieldeqa, fieldeqp` | complex F3 | per call | | F3 | not exported (RECOMPUTE, 04) |
| `vpa_gf, vperp2_gf, aj0_gf, aj1_gf, c_rate, g_restart_tmp` | | gf_lo / diagnostics | | | not exported (GATE / diagnostics only) |

### 3.8 `/dist_fn` (dist_fn.fpp; most entries private)

| name | decl (line) | computed at | dt | access | tag | decision (owner) |
|---|---|---|---|---|---|---|
| `vpar` | real G7 (153) | init_vpar 1681-1729 | **dt** | PRIV→§4 | G7 real | BOTH (02,03) |
| `wdrift, wdriftttp` | real G7 (136,140; complex under LOWFLOW) | init_wdrift 1405-1498 | **dt** | PRIV→§4 | G7 real | BOTH (02,03) |
| `a, b, r, ainv` | complex G7 (147) | init_invert_rhs 1826-1929 | **dt** | PRIV→§4 | G7 | BOTH (02,03) |
| `g_h` | complex G7 (156/158) | init_homogeneous_g 1937-2125 | **dt** | PRIV→§4 | G7 | BOTH (03) |
| `g_adj` | complex `(n_links_max, 2, g_lo)` (162), linked only (3599) | zeroed at init; h-slots filled at init, p-slots every call | per call | PRIV→§4 | GL | not exported by default (03: RECOMPUTE as gathers); `write_g_adj` knob for debugging |
| `ufac` | real `(negrid, nspec)` (169) | init_source_term 2180-2187 | no | PRIV→§4 | V `(egrid, species)` | EXPORT (02) |
| `source_coeffs_phim, _phip, _aparm, _aparp` | complex G7m1 (172-173), iff `opt_source` (3951-3957) | init_source_term 2159-2270 | **dt** | PRIV→§4 | G7m1 | BOTH (02) |
| `gamtot, gamtot1, gamtot2` | real F3 (144) | init_fieldeq 6470-6572 | no | pub | F3 | EXPORT (04,07) |
| `gamtot3` | real F3 (144), iff fieldlineavg | 6511-6513, 6564-6565 | no | PRIV→§4 | F3 | EXPORT if allocated (04) |
| `gridfac1` | real F3 (150) | 6486-6506 | no | pub | F3 | EXPORT (04) |
| `inv_phi_denominator_g, inv_bpar_denominator_g` | real F3 (199), lazily at first get_fields_direct_from_dfn (7457-7585) | | no | PRIV→§4 | F3 | BOTH if allocated (04) |
| `awgt` | real `(ntheta0, naky)` (176), lazily at first calculate_flux_surface_average (8006-8020) | | frozen | PRIV→§4 | K2 | BOTH if allocated (04) |
| `fl_avg` | complex `(ntheta0, naky)` (177) | per call | | pub | — | not exported (RECOMPUTE, 04) |
| `itleft, itright` | int K2 (180), linked only (3534) | compute_itleft_and_itright 2441-2513 | no | pub | K2 int32 | EXPORT (03,04,07) |
| `l_links, r_links` | int K2 (223) | count_links 2645-2672 | no | pub | K2 int32 | EXPORT |
| `n_links` | int `(2, ntheta0, naky)` (224) | count_links | no | PRIV→§4 | V `(two, kx, ky)` | EXPORT-for-check (03 RECOMPUTE) |
| `M_class, N_class` | int `(i_class)` (233), `i_class` (234) | setup_fields_implicit_class_arrays 2676-2755 / 2144-2150 | no | pub | V `(n_class)` | EXPORT (04) |
| `supercell_labels, n_cells` | int K2 (188,196); `n_supercells` `(naky)` (191) | calculate_supercell_labels 2355-2397 | no | pub | K2 / V | EXPORT |
| `connections(iglo)%{iproc_left, iglo_left, iproc_right, iglo_right, neighbor}` | type private (205-211) | compute_connections 2518-2565 | no | PRIV→§4 (variable; components readable) | G5 int32 / int8 | EXPORT `connections_neighbor`, `connections_it_left/right` (= `it_idx` of `iglo_left/right`, or −1) for checks (03) |
| `save_h` | logical `(2, g_lo)` (225) | set_save_h 2570-2604 | no | PRIV→§4 | G5s int8 | EXPORT-for-check (03) |
| `no_connections` | logical (231) | setup_connected_bc_redistribute 2870-2876 | no | PRIV→§4 | scalar | EXPORT (03) |
| `jump, ikx_indexed` | int `(naky)`, `(ntheta0)` (166), lazily in exb_shear (4237-4249) | first shear call | state | PRIV→§4 | V | trace only (07: RECOMPUTE) |
| `exb_first` | logical (238) | `.true.` until first shear call | state | PRIV→§4 | scalar | STATE |
| knobs (effective values): `fexp, bkdiff` `(nspec)` (99); `gridfac, apfac, driftknob, tpdriftknob, poisfac, vparknob` (100); `t0, omega0, gamma0, source0` (101); `phi_ext, afilter, exponential_boundary_factor` (102); `wfb, g_exb, g_exbfac, omprimfac, btor_slab, mach` (103); `g_exb_start_time, g_exb_error_limit` (104); `g_exb_start_timestep` (105); `hyper_in_initialisation` (106); `start_from_previous_solution, boundary_off_grid, exponential_boundary` (107); `nonad_zero, lf_default, lf_decompose, esv, opt_source` (108); `adiabatic_option_switch` (110); `source_option_switch` (116); `boundary_option_switch` (120); `def_parity, even` (127); `zero_forbid` (128); `gf_lo_integrate` (129); `mult_imp` (130) | | read_parameters 1174-1400 (forced values: `bkdiff` 1319-1344, `esv` 1327, `fexp` 1397-1398) | no | public: `apfac, omega0, gamma0, g_exb, g_exbfac, adiabatic_option_switch, boundary_option_switch, def_parity, even, gf_lo_integrate`; the rest PRIV→§4 | scalars / V `(species)` | EXPORT effective (02,03,04,07) |
| `get_leftmost_it(it,ik)`, `get_rightmost_it(it,ik)` | elemental fns 3889-3920 | | no | leftmost pub, rightmost PRIV→§4 | K2 int32 `itleftmost`, `itrightmost` | EXPORT-for-check (03) |
| mixed-wfb end factors `aj0l, aj1l, aj0r, aj1r` per (ie,it,ik,is); `vperp2left, vperp2right` per (ie,is) | locals of apply_linked_boundary_conditions 6016-6028 | | no | derived (verbatim, §7.7) | V `(egrid, kx, ky, species)`, `(egrid, species)` | EXPORT (03 §3.11.4) |
| redistribute maps `links_p, links_h, wfb_p, wfb_h, pass_right/left, incoming_links, parity_redist, gc_from_left/right` | `redist_type` (215-221) | init_connected_bc | no | — | — | not exported (03, 07: RECOMPUTE as gather tables) |

### 3.9 `/fields_implicit` (fields_implicit.f90) and `/fields_local`

| name | decl | computed at | dt | access | tag | decision |
|---|---|---|---|---|---|---|
| `aminv(c)%dcell(dc)%supercell(:)` | private derived types (24-47) | init_inverse_matrix 722-1071 | **dt** (also vnmult, hyper, collisions at build time) | PRIV→§4 (variable; components readable) | per class `aminv_c<c>` `(ri, sc_R_c<c>, sc_R_c<c>, sc_M_c<c>)` = `aminv(col, row, m)` | BOTH (04 §8) |
| `nfield, nidx` | int (57) | init_response_matrix 477-481 | no | PRIV→§4 (or `jf_lo%nfield`, `f_lo(1)%nidx`) | scalar | EXPORT |
| `f_lo(c)%ik(m,n), f_lo(c)%it(m,n)` | gs2_layouts pub | finish_fields_layouts 1074-1133 | no | pub | `cls_ik`, `cls_it` `(sc_Mmax, sc_Nmax, n_class)` padded −1 | EXPORT (04 §8) |
| derived `sc_class, sc_member, sc_cell` | K2 | derived from `f_lo` | no | derived | K2 int32 | EXPORT (convenience) |
| `am` (pre-inversion response) | local of init_response_matrix (473, 548) | 548-635 | **dt** | probe P-R1 | `am_c<c>` `(ri, sc_R_c<c>, sc_RM_c<c>)` | optional (04 Q10) |
| `fieldmat` (fields_local) | public variable of private types; type-bound procedures `contains private` (fields_local.f90:126-127) | init_fields_matrixlocal | **dt** | pub variable, methods inaccessible | per supercell `resp_ik<ik>_is<is>` `(ri, ax1, ax2)` via `fields::dump_response_to_file` round trip (§7.5) | EXPORT iff `field_option='local'` (04 §7.3) |

### 3.10 `/collisions` (collisions.fpp; all arrays private, public list 17-23)

All arrays and scalars of spec 06 §10 "Export checklist", with the LE6 canonical layout of spec 06 §2.4:

| name | decl (line) | computed at | dt | access | tag | decision |
|---|---|---|---|---|---|---|
| `vnew, vnew_s, vnew_D, vnew_E, delvnew` | real `(naky, negrid, nspec)` (87) | init_vnew 1707-1876 | no | PRIV→§4 | V `(ky, egrid, species)` | EXPORT |
| `vnewh` | real `(-ntgrid:ntgrid, ntheta0, naky, nspec)` (95) | 1853-1863 iff hyper_colls | no | PRIV→§4 | V `(theta, kx, ky, species)` | EXPORT |
| `pitch_weights` | real `(nlambda, -ntgrid:ntgrid)` (133) | 2199-2210 | no | PRIV→§4 | V `(lambda, theta)` | EXPORT |
| `vpdiff` | real `(-ntgrid:ntgrid, 2, nlambda)` (91) | init_vpdiff 4088-4140 | no | PRIV→§4 | V `(theta, sign, lambda)` | EXPORT iff conservative |
| `vpdiffle` | real `(nxi_lim, -ntgrid:ntgrid)` (90) | 4142-4152 | no | PRIV→§4 | V `(xi_lim, theta)` | EXPORT iff le + conservative |
| `c1le, betaale, qle` | real `(nxi_lim, negrid, le_lo)` (121) | init_lorentz 2216, 2240-2285 | **dt**, vnmult | PRIV→§4 | LE6 `(xi_lim, egrid, …)` | BOTH |
| `ec1le, ebetaale, eqle` | same (123) | init_ediffuse 1968, 1981-1999 | **dt**, vnmult | PRIV→§4 | LE6 | BOTH |
| `z0le, s0le, w0le, bz0le, bs0le, bw0le` | real `(nxi_lim, negrid+1, le_lo)` (102,113) | 1158-1168, 1647-1654 | **dt**, vnmult | PRIV→§4 | LE6 `(xi_lim, egrid_p1, …)` | BOTH |
| `aj0le, vperp_aj1le, vpa_aj0_le` | same (105) | init_le_bessel 634-687 | no | PRIV→§4 | LE6 `(xi_lim, egrid_p1, …)` | EXPORT |
| `c1, betaa, ql, ec1, ebetaa, eql` | `(nxi_lim, lz_lo)` / `(negrid, e_lo)` (128,130) | 2289-2353, 2016-2033 | **dt** | PRIV→§4 | LZ6 / E7 | BOTH iff `use_le_layout=F` |
| `s0, w0, z0, bs0, bw0, bz0` | real G7 (99,109) | 1172-1192, 1657-1677 | **dt** | PRIV→§4 | G7 real | BOTH iff `use_le_layout=F` |
| `vnmult` | real(2) (72) | init_collisions 1731-1733 (≥ 1 after clamps), restart | state | pub | V `(two)` | EXPORT post-init |
| scalars: `collision_model_switch, colls, use_le_layout, adjust, heating, hyper_colls, split_collisions, vary_vnew, ncheck` | | read_parameters | no | pub | scalar | EXPORT |
| scalars: `has_lorentz, has_diffuse, drag, const_v, conserve_moments, conservative, resistivity, lorentz_switch, ediff_switch, ei_coll_only, special_wfb_lorentz, vpar_zero_mean, conserve_forbid_zero, force_collisions, hypermult, cfac, nxi_lim, vnfac, vnslow, vnm_init, timesteps_between_collisions` | (46-139) | read_parameters, 557-561 | no | PRIV→§4 | scalar | EXPORT |
| helper `e15(ie) = energy_maxw(ie)**1.5` | | | no | derived (§7.7) | V `(egrid)` | EXPORT (06 §3.4) |
| helper `capg_half(j)` | internal procedure `capg` | | no | derived (verbatim copy, §7.7) | V `(egrid_m1)` | EXPORT (06 §3.4) |

### 3.11 `/hyper` (hyper.f90)

| name | decl | computed at | dt | access | decision |
|---|---|---|---|---|---|
| `D_res` | real `(ntheta0, naky)` (38) | init_hyper 400-411 (0 at 615) | **frozen at initial dt** | pub | BOTH (02, 07) — K2 |
| `hypervisc_filter` | real F3 (41) | simple3D 419-445 (frozen); iso/aniso overwritten every call (738-743, 777-778) | yes | pub | EXPORT at I only for simple3D; per probe `after_hyper` (state for diagnostics, 08) |
| `D_v, D_eta, nexp` | scalars | 378-379 | no | pub | EXPORT |
| `D_hypervisc, D_hyperres, omega_osc, D_hyper, akx4_max, aky4_max, aky_max, akperp4_max, hyper_option_switch, const_amp, include_kpar, isotropic_shear, damp_zonal_only, hyper_on, gridnorm, kperp2_max, kperp2_max_zonal, D_hyper3D, P_hyper3D, ky_cut, kx_cut, isotropic_model` | scalars (22-49) | init_hyper / read_parameters 452-601 | no | PRIV→§4 | EXPORT effective (07 recomputes from config; export lets the loader assert) |

### 3.12 `/nonlinear_terms` (nonlinear_terms.fpp)

| name | decl (line) | computed at | dt | access | decision |
|---|---|---|---|---|---|
| `nonlin, split_nonlinear, accelerated, cfl` | (78, 72, 84, 56) | read_parameters, gs2_transforms | no | pub | EXPORT |
| `cflx, cfly` | real (56) | 357-358 | no | PRIV→§4 | BOTH (05) |
| `max_vel_components` | real(3) (62) | add_explicit | per step | PRIV→§4 (also passed by probe `nl_cfl`) | trace |
| `istep_last` | int (32), init 0 | add_explicit 497/611 | state | PRIV→§4 | STATE (05: R) |
| `initializing` | logical (80) | `.false.` after `nl_finish_init` (fields.f90:373) | state | PRIV→§4 | STATE |
| `nonlinear_mode_switch, use_cfl_limit, use_2d_cfl, error_target, istep_error_start, use_order_based_error, include_apar, include_bpar, include_phi, zip, nl_forbid_force_zero` | (35, 65, 68-70, 73, 82-83) | read_parameters | no | PRIV→§4 | EXPORT |
| transform scales `xb_fft%scale, yb_fft%scale`, dealias mask | gs2_transforms | | no | — | RECOMPUTE (05) |

### 3.13 Local / SAVE variables that cannot be exported (RECOMPUTE in JAX)

| variable | location | value at the start of JAX step 1 | owner rule |
|---|---|---|---|
| `timeadv::istep_last` (SAVE, init −1) | dist_fn.fpp:4031 | 0 if the response matrix was built (timeadv called with istep 0), −1 if it was read from file and never built | only used by `split_nonlinear` (GATED) |
| `exb_shear::istep_last` (SAVE, init 0), `itmin(1)`, `dkx`, `dtheta0` | dist_fn.fpp:4212-4219 | 0, unset, unset (set in the `exb_first` block) | 07 §3.3 |
| `reset_time_step::istep_last` (SAVE −1), `nconsec` (SAVE 0) | gs2_reinit.f90:120-121 | −1, 0 | 01 §2.2 |
| `getfield` `recvcnts, displs` | fields_implicit.f90:62 | MPI bookkeeping | not needed |
| fields_local `les_comm`, `first` | fields_local.f90:2610-2612 | MPI bookkeeping | not needed |

The trace (§6.2) lets the JAX side validate its recomputation of these through their observable effects.
If a later phase needs to *restart* JAX from a trace snapshot taken mid-run, these SAVE variables must be
moved to module scope (behaviour-preserving edit, out of scope here; Open question Q11).

---------------------------------------------------------------------------------------------------

## 4. Access to private module data — options, decision and exact edits

### 4.1 Options considered

| option | how | pros | cons | verdict |
|---|---|---|---|---|
| A. `public ::` statements | one block per module listing the needed private names | trivial, no logic in physics modules, the exporter (one place) owns all re-indexing; verified: changing PRIVATE→PUBLIC only changes the symbol visibility (`.private_extern` → external), the generated instructions of a test module are identical with `-O3 -fPIC -mcpu=apple-m1` (Appendix A.2) | widens the namespace (harmless: no module in `src/` or `tests/` does a bare `use dist_fn`/`collisions`/… without `only:`, checked by grep) | **chosen** |
| B. accessor routines inside each module (`get_dist_fn_coeff(name, buf)`) | copy-out with `select case` | keeps names private | ~15 accessors, duplicated re-indexing, more code in physics modules | rejected |
| C. export routine inside each module behind `#ifdef GS2_EXPORT` | per-module writer | — | requires a second object tree (the flag changes the objects), netCDF code in physics modules, circular-dependency risk (le_grids ↔ writer) | rejected |
| D. Fortran 2008 submodules | `submodule (dist_fn) dist_fn_export` | private access without changing declarations | parent module must still declare the separate module procedure interface; fortdep supports submodules (`scripts/fortdep` "submodule" handling, `get_dependency`) but none exist in GS2 yet; links the export code into every binary that uses the parent (fortdep makes the parent require its submodules) | rejected |
| E. pointers to module arrays | `type(view)` with `pointer` components | — | requires `TARGET` on the module arrays, which changes alias analysis and may change code generation | **forbidden** (Pitfall 11.13) |

**Rule for the implementation phase:** add only `public ::` statements (block of §4.2) and the probe call
lines of §5.2. Do not add `target`, `save`, `volatile` or `contiguous` attributes, do not reorder
declarations, do not change any executable statement.

### 4.2 Exact `public ::` blocks

Each block is inserted directly after the last existing `public` statement of the module's
specification part (anchor text quoted). Header comment for every block:
`! GS2 GPU export (gs2_export): read-only access, see gpu/spec/09_export_infra.md §4`.

**`src/dist_fn.fpp`** — after line 96 `public :: check_getan`:

```fortran
  public :: vpar, wdrift, wdriftttp, a, b, r, ainv, g_h, g_adj, ufac
  public :: source_coeffs_phim, source_coeffs_phip, source_coeffs_aparm, source_coeffs_aparp
  public :: gamtot3, awgt, inv_phi_denominator_g, inv_bpar_denominator_g
  public :: n_links, save_h, connections, no_connections, jump, ikx_indexed, exb_first
  public :: get_rightmost_it
  public :: fexp, bkdiff, gridfac, driftknob, tpdriftknob, poisfac, vparknob, t0, source0
  public :: phi_ext, afilter, exponential_boundary_factor, wfb, omprimfac, btor_slab, mach
  public :: g_exb_start_time, g_exb_error_limit, g_exb_start_timestep, hyper_in_initialisation
  public :: start_from_previous_solution, boundary_off_grid, exponential_boundary
  public :: nonad_zero, lf_default, lf_decompose, esv, opt_source, source_option_switch
  public :: zero_forbid, mult_imp
```

(`connections_type` itself stays private; the variable `connections` of that type may be public in
Fortran ≥ 2003 and its default-public components are readable from outside — verified with
`gfortran -std=f2008`, Appendix A.3.)

**`src/collisions.fpp`** — after line 23 `public :: collisions_config_type, set_collisions_config, get_collisions_config`:

```fortran
  public :: vnew, vnew_s, vnew_D, vnew_E, delvnew, vnewh, pitch_weights, vpdiff, vpdiffle
  public :: c1le, betaale, qle, ec1le, ebetaale, eqle
  public :: z0le, s0le, w0le, bz0le, bs0le, bw0le, aj0le, vperp_aj1le, vpa_aj0_le
  public :: c1, betaa, ql, ec1, ebetaa, eql, s0, w0, z0, bs0, bw0, bz0
  public :: nxi_lim, cfac, drag, has_lorentz, has_diffuse, const_v, conserve_moments
  public :: conservative, resistivity, lorentz_switch, ediff_switch, ei_coll_only
  public :: special_wfb_lorentz, vpar_zero_mean, conserve_forbid_zero, force_collisions
  public :: hypermult, vnfac, vnslow, vnm_init, timesteps_between_collisions
```

**`src/hyper.f90`** — after line 17 `public :: get_hyper_config`:

```fortran
  public :: D_hypervisc, D_hyperres, omega_osc, D_hyper
  public :: akx4_max, aky4_max, aky_max, akperp4_max, kperp2_max, kperp2_max_zonal
  public :: hyper_option_switch, const_amp, include_kpar, isotropic_shear, damp_zonal_only
  public :: hyper_on, gridnorm, D_hyper3D, P_hyper3D, ky_cut, kx_cut, isotropic_model
```

**`src/nonlinear_terms.fpp`** — after line 24 `public :: get_nonlinear_terms_config`:

```fortran
  public :: cflx, cfly, max_vel_components, istep_last, initializing, nonlinear_mode_switch
  public :: use_cfl_limit, use_2d_cfl, error_target, istep_error_start, use_order_based_error
  public :: include_apar, include_bpar, include_phi, zip, nl_forbid_force_zero
```

(The exporter imports `istep_last` renamed, `use nonlinear_terms, only: nl_istep_last => istep_last`.)

**`src/fields_implicit.f90`** — after line 13 `public :: dump_response_to_file_imp`:

```fortran
  public :: aminv, nfield, nidx
```

**`src/theta_grid.f90`** — after line 21 `public :: get_theta_grid_config`:

```fortran
  public :: field_line_average_weight, cvdriftknob, gbdriftknob, eqopt_switch
```

**`src/le_grids.f90`** — after line 22 `public :: le_grids_config_type, set_le_grids_config, get_le_grids_config`:

```fortran
  public :: bouncefuzz, trapped_particles, wfbbc_option_switch, ngauss, npassing, nesuper
  public :: genquad, radau_gauss_grid, split_passing_region, new_trap_int_split
```

**`src/run_parameters.f90`** — after line 23 `public :: set_overrides, get_knobs_config, knobs_config_type, set_run_parameters_config`:

```fortran
  public :: eqzip_option_switch, seed
```

No edits are needed in `dist_fn_arrays.fpp`, `fields_arrays.f90`, `kt_grids.f90`, `species.f90`,
`gs2_time.f90`, `gs2_layouts.fpp`, `gs2_reinit.f90` (use `get_reinit_config()`), `init_g.f90` (use
`get_init_g_config()`), `fields.f90` (use `get_fields_config()`), `fields_local.f90` (`fieldmat` is public
but its type-bound procedures are private, `fields_local.f90:126-127`; use the public
`fields::dump_response_to_file`, §7.5), `gs2_diagnostics_new` (`gnostics` public).

### 4.3 What `public` cannot reach, and the method used

| item | why | method |
|---|---|---|
| `aminv` element data | private *types* `dcell_type`, `field_matrix_type` | component access through the now-public variable (`aminv(c)%dcell(dc)%supercell(:)`), allowed because components are default-public |
| `fieldmat` matrices (fields_local) | type-bound procedures `pull_rows_to_arr` are `private` (`contains / private`) | `fields::dump_response_to_file(suffix)` (fields.f90:380-397) → per-supercell files → read back (§7.5) |
| `am` (pre-inversion response) | local to `init_response_matrix` | optional probe P-R1 (§5.2) |
| `gle` in collisions | local of `vspace_derivatives` | explicit-array probes (§5.2) |
| `g1` in `add_explicit` | dummy argument associated with `gexp_1` and modified through the dummy | explicit-array probe (§5.2), never read `gexp_1` through the module inside `add_explicit` (aliasing rule F2018 15.5.2.13 / F2008 12.5.2.13) |
| mixed-wfb end factors | locals of `apply_linked_boundary_conditions` | recomputed in the exporter with the verbatim expressions (§7.7) |
| `capg` | internal procedure of collisions | verbatim copy in the exporter (spec 06 §3.4) |
| SAVE locals (§3.13) | function-local | RECOMPUTE in JAX |

### 4.4 Proving the edits do not change the numerics

1. Before editing, keep copies of `bin/gs2` (`bin/gs2.orig`) and of `objs/{dist_fn,collisions,hyper,
   nonlinear_terms,fields_implicit,theta_grid,le_grids,run_parameters}.o`.
2. After adding only the §4.2 blocks and rebuilding, compare `otool -tv old.o` with `otool -tv new.o` for
   each of these objects: the text sections must be identical (only `nm -m` visibility changes).
3. After adding the probe lines (§5.2) — which do change the generated code — run `bin/gs2.orig` and the new
   `bin/gs2` (probes inactive) on the verification inputs of §12 and require **bitwise** equality of every
   numeric variable of `<run>.out.nc` and of the restart files (excluding `date_created`, `id`/UUID,
   timing attributes). If any difference appears, fall back to compiling the probe lines under
   `#ifdef GS2_EXPORT_PROBES` with a separate object tree (§9.3) so that `bin/gs2` is unaffected.

---------------------------------------------------------------------------------------------------

## 5. Stage probes (kernel-level verification hooks)

### 5.1 Module `gs2_export_probe` (new file `src/gs2_export_probe.f90`, no GS2 dependencies)

The physics modules call the probe; the probe forwards to procedure pointers that only `gs2_export`
installs. The module uses nothing from GS2, so it can be used from every level (gs2_main down to
collisions) without dependency cycles, and it costs one logical test per call site when inactive.

```fortran
module gs2_export_probe
  implicit none
  private
  public :: probe_active, probe_in_advance, probe_istep, probe_pass
  public :: probe_point, probe_g, probe_le, probe_real, probe_c2
  public :: probe_point_hook, probe_g_hook, probe_le_hook, probe_real_hook, probe_c2_hook

  !> Master switch, only ever set by gs2_export.
  logical :: probe_active = .false.
  !> True between the start and the end of fields::advance (set by the probes in fields.f90).
  logical :: probe_in_advance = .false.
  !> istep of the current advance, and which timeadv call (1 or 2; 0 outside timeadv).
  integer :: probe_istep = -1, probe_pass = 0

  abstract interface
     subroutine point_iface(tag, istep)
       character(len=*), intent(in) :: tag
       integer, intent(in) :: istep
     end subroutine point_iface
     subroutine g_iface(tag, istep, arr)           ! arr(1:nth, 1:2, 1:nloc), nloc = g_lo%ulim_alloc-g_lo%llim_proc+1
       character(len=*), intent(in) :: tag
       integer, intent(in) :: istep
       complex, dimension(:,:,:), intent(in) :: arr
     end subroutine g_iface
     subroutine le_iface(tag, istep, arr)          ! arr(1:nxi+1, 1:negrid+1, 1:nloc), le_lo
       character(len=*), intent(in) :: tag
       integer, intent(in) :: istep
       complex, dimension(:,:,:), intent(in) :: arr
     end subroutine le_iface
     subroutine real_iface(tag, istep, vec)
       character(len=*), intent(in) :: tag
       integer, intent(in) :: istep
       real, dimension(:), intent(in) :: vec
     end subroutine real_iface
     subroutine c2_iface(tag, ic, arr)             ! response build: am(1:R_c, 1:local columns)
       character(len=*), intent(in) :: tag
       integer, intent(in) :: ic
       complex, dimension(:,:), intent(in) :: arr
     end subroutine c2_iface
  end interface

  procedure(point_iface), pointer :: probe_point_hook => null()
  procedure(g_iface),     pointer :: probe_g_hook     => null()
  procedure(le_iface),    pointer :: probe_le_hook    => null()
  procedure(real_iface),  pointer :: probe_real_hook  => null()
  procedure(c2_iface),    pointer :: probe_c2_hook    => null()

contains
  subroutine probe_point(tag, istep)
    character(len=*), intent(in) :: tag
    integer, intent(in), optional :: istep
    if (.not. probe_active) return
    if (.not. associated(probe_point_hook)) return
    if (present(istep)) then
       call probe_point_hook(tag, istep)
    else
       call probe_point_hook(tag, probe_istep)
    end if
  end subroutine probe_point
  ! probe_g, probe_le, probe_real, probe_c2: identical pattern (tag, [istep defaulting to probe_istep], array)
end module gs2_export_probe
```

`complex`/`real` here are default kinds = 8 bytes under `-fdefault-real-8`, the same as the callers.
Assumed-shape dummies lose the lower bounds: the hook maps `arr(i, s, j)` to `ig = i - 1 - ntgrid`,
`isgn = s`, `iglo = g_lo%llim_proc + j - 1` (and `ile = le_lo%llim_proc + j - 1` for `probe_le`).

### 5.2 Insertion points (exact)

Call-site idiom (one line, plus a local `use gs2_export_probe, only: …` in the subroutine's use list):
`if (probe_active) call probe_point('tag')`. Explicit-array probes pass the dummy/local that holds the
data. "after L" = insert a new line after source line L (re-anchor on the quoted statement).

| tag | file : after line (statement) | probe | condition / notes |
|---|---|---|---|
| `advance_begin` | fields.f90 : after 593 (`integer, intent (in) :: istep`), before the `select case` at 595 | `if (probe_active) then; probe_in_advance=.true.; probe_istep=istep; probe_pass=0; call probe_point('advance_begin', istep); end if` | every attempt |
| `advance_end` | fields.f90 : after 604 (`end select`) | `if (probe_active) then; call probe_point('advance_end', istep); probe_in_advance=.false.; probe_pass=0; end if` | also reached when advance_* returned early on `reset` |
| `after_exb_shear` | fields_implicit.f90 : after 335 (`… call exb_shear (gnew, …)`); fields_local.f90 : after 4280 | `probe_point` | only meaningful if shear active |
| (pass 1) | fields_implicit.f90 : before 341 (`call timeadv (phi, … istep)`); fields_local.f90 : before 4296 | `if (probe_active) probe_pass = 1` | |
| `after_timeadv` | fields_implicit.f90 : after 341; fields_local.f90 : after 4297 | `probe_point` | pass 1 |
| `after_getfield` | fields_implicit.f90 : after 352 (`call getfield (phinew, aparnew, bparnew)`); fields_local.f90 : after 4313 | `probe_point` | implicit: `phinew…` hold the increment `u`; local: depends on `do_update` |
| `after_field_update` | fields_implicit.f90 : after 356 (`bparnew  = bparnew + bpar`); fields_local.f90 : after 4320 (`endif`) | `probe_point` | |
| `after_remove_zonal_flows` | fields_implicit.f90 : after 362; fields_local.f90 : after 4331 | `probe_point` | only if `remove_zonal_flows_switch` |
| (pass 2) | fields_implicit.f90 : before 365; fields_local.f90 : before 4336 | `if (probe_active) probe_pass = 2` | |
| `after_timeadv` | fields_implicit.f90 : after 365; fields_local.f90 : after 4337 | `probe_point` | pass 2 |
| `after_collisions_advance` | fields_implicit.f90 : after 369; fields_local.f90 : after 4340 | `probe_point` | |
| `after_add_explicit_terms` | dist_fn.fpp : after 4062 (`end if` closing the `split_nonlinear` branch), **before** the reset return at 4066 | `probe_point` | `gexp_*` are module arrays and no longer associated with a dummy here |
| `source` | dist_fn.fpp : after 6182 (end of `call get_source_term_on_local_domain(g_work, …)`) | `probe_point` | dumps `g_work`; must precede `invert_rhs_1`, which modifies `source` in place (`intent(in out)`, dist_fn.fpp:5452; 03 §2.1) |
| `after_invert_rhs_1` | dist_fn.fpp : after 6184 | `probe_point` | |
| `after_invert_rhs` | dist_fn.fpp : after 4068 (`call invert_rhs (…)`) | `probe_point` | includes `apply_linked_boundary_conditions` |
| `after_hyper` | dist_fn.fpp : after 4077 (`call hyper_diff (gnew, phinew)`) | `probe_point` | dumps `gnew`, `hyper::hypervisc_filter` |
| `after_vspace_derivatives` | dist_fn.fpp : after 4084 | `probe_point` | inside `if (.not. split_collisions)` |
| `after_parity` | dist_fn.fpp : after 4091 | `probe_point` | `def_parity` only |
| `coll_h_in` | le_derivatives.fpp : after 71 (le path) and after 133 (standard path) | `probe_g('coll_h_in', arr=g)` | `g` is the dummy for `gnew`, modified through the dummy |
| `coll_gle_in` | le_derivatives.fpp : after 83 (`if (g_to_c) call gather (g2le, g, gle)`) | `probe_le('coll_gle_in', arr=gle)` | local `gle(nxi+1, negrid+1, le_lo)` |
| `coll_after_ediffuse` | collisions.fpp : after 3051 (`call solfp_ediffuse_le_layout (gle)`) | `probe_le(…, gle)` | |
| `coll_after_ediffuse_conserve` | collisions.fpp : after 3052 (`if (conserve_moments) call conserve_diffuse (gle)`) | `probe_le` | |
| `coll_after_drag` | collisions.fpp : after 3085 (`end if` of the drag block) | `probe_le` | |
| `coll_after_lorentz` | collisions.fpp : after 3087 (`call solfp_lorentz (gle, diagnostics)`) | `probe_le` | |
| `coll_after_lorentz_conserve` | collisions.fpp : after 3088 (`if (conserve_moments) call conserve_lorentz (gle)`) | `probe_le` | |
| `coll_gle_out` | le_derivatives.fpp : after 90 (`call solfp1 (gle, diagnostics)`) | `probe_le` | |
| `coll_h_out` | le_derivatives.fpp : after 94 (`scatter`) and after 140 (standard path, `solfp1 (g, …)`) | `probe_g(…, g)` | |
| `nl_before_center` | nonlinear_terms.fpp : after 520 (`call add_nl (g, g1, phi, apar, bpar)`) | `probe_g('nl_before_center', arr=g1)` | only on the first `add_explicit` of an istep |
| `nl_after_center` | nonlinear_terms.fpp : after 523 (`call nl_center (g1, bd)`) | `probe_g(…, g1)` | |
| `nl_cfl` | nonlinear_terms.fpp : after 569 and after 582 (`call set_cfl_time_step_and_check_if_too_large`) | `probe_real('nl_cfl', vec=max_vel_components)` | `code_dt_cfl`, `reset` read by the hook |
| `after_reset_time_step` | gs2_main.fpp : after 744 and after 811 (`call reset_time_step (…)`) | `probe_point('after_reset_time_step', istep)` | outside advance |
| `step_end_pre_diag` | gs2_main.fpp : after 771 (`call update_time`) | `probe_point(…, istep)` | |
| `step_end` | gs2_main.fpp : after 785 (`if (state%exit) state%converged = .true.`) | `probe_point(…, istep)` | state seen by the next step |
| P-R1 `response_am` | fields_implicit.f90 : before 638 (`call init_inverse_matrix (am, i)`) | `probe_c2('response_am', i, am)` | only if `trace_response_build`; runs during init and on every dt reset |

All probe sites are at subroutine level (never inside loops over `iglo`/`ile`), and every rank executes
them in the same order (the time step is SPMD), so hooks may perform collective gathers (§2.5).

### 5.3 Hook semantics (implemented in the exporter, module `gs2_export_trace`)

* **Filtering.** Tags inside `timeadv` and below (`after_add_explicit_terms … coll_* … nl_*`) are recorded
  only when `probe_in_advance`; this drops the `istep = 0` calls of the response build (init and every
  reset) except P-R1. Tags `after_reset_time_step`, `step_end_pre_diag`, `step_end` are recorded
  unconditionally while `probe_active`. Knobs `trace_stages`, `trace_first_step`, `trace_every` (§8.1)
  filter further.
* **Attempts.** The hook keeps `attempt(istep)`; `advance_begin` for an `istep` already seen increments it
  (immediate-reset path, gs2_main.fpp:739-752). In attempt 2 the `nl_*` tags are absent because
  `add_explicit` skips the NL term when `istep == istep_last` (nonlinear_terms.fpp:495-497).
* **Sequence number.** Each recorded probe gets `seq` = running count within `(istep, attempt)`, so the two
  passes of the same tag are distinct and the order of execution is preserved.
* **Contents.** `probe_point` dumps module state chosen by tag (table below); explicit-array probes dump
  the array passed. Every record also writes the *time block* (§7.6) and `reset`.

| tag family | arrays dumped (canonical) |
|---|---|
| `advance_begin`, `advance_end`, `step_end*` | `gnew` (G7), `g` (G7, only `step_end`), `phinew, aparnew, bparnew` (F3), `phi, apar, bpar` (F3, `advance_end` only), `gexp_1..3` (G7, if allocated), `kx_shift`/`theta0_shift`, `exb_first`, `nl_istep_last`, `vnmult` |
| `after_exb_shear` | `gnew`, `phinew…`, `kx_shift`/`theta0_shift`, `jump`, `ikx_indexed` (if allocated) |
| `after_add_explicit_terms` | `gexp_1..3`, `code_dt_cfl`, `reset` |
| `source` | `g_work` |
| `after_invert_rhs_1`, `after_invert_rhs`, `after_vspace_derivatives`, `after_parity`, `after_timeadv`, `after_collisions_advance` | `gnew` |
| `after_hyper` | `gnew`, `hypervisc_filter` (F3 real) |
| `after_getfield`, `after_field_update`, `after_remove_zonal_flows` | `phinew, aparnew, bparnew` |
| `coll_h_*`, `nl_before/after_center` | the passed G7 array |
| `coll_gle_*`, `coll_after_*` | the passed LE array: LE6 complex `(ri, xi_p1, egrid_p1, theta, kx, ky, species)` |
| `nl_cfl` | `max_vel_components(3)`, `code_dt_cfl`, `reset` |
| `after_reset_time_step` | time block only (and, with P-R1, the rebuilt `am_c*`) |
| `response_am` | buffered in memory per class; written as `/fields_implicit/am_c<c>` in the export file (init) or into the trace record of the reset |

---------------------------------------------------------------------------------------------------

## 6. Output formats

### 6.1 Export file `<run_name>.gs2export.nc` (NETCDF4, written by proc0)

#### 6.1.1 Global attributes and root dimensions

Global attributes: everything from `gs2_metadata::create_metadata(ncid, "GS2 GPU export", header)`
(`src/gs2_metadata.fpp:20-57`: title, software_name, software_version, netcdf_version, date_created,
id (run UUID), Conventions, references) plus

| attribute | value |
|---|---|
| `export_format_version` | `"1.0"` (bump on any incompatible change) |
| `git_hash`, `git_state` | `get_git_hash()`, `get_git_state()` (`git_version_mod`) |
| `build_config` | `build_config::formatted_build_config()` (the `--build-config` text) |
| `run_name`, `input_file` | `file_utils::run_name`, `trim(run_name)//'.in'` |
| `nproc`, `nproc_g_lo`, `nproc_le_lo` | `mp::nproc`, `g_lo%nproc`, `le_lo%nproc` |
| `layout` | `gs2_layouts::layout` |
| `index_convention` | `"Fortran dimension order; NumPy shows reversed dims. Canonical g: (theta,sign,kx,ky,lambda,egrid,species). Complex: leading ri dim (1=re,2=im). Logical: int8."` |
| `real_kind_bytes` | `8` (abort if `kind(1.0) /= real64`) |
| `gate_status` | `"ok"` or a `;`-separated list of soft-gate violations (§8.3) |
| `diagnostics_run_at_t0` | `0/1` (whether `state_t0` is post-`initialize_diagnostics`) |
| `used_inputs` | text of `<run>.used_inputs.in` (written by `write_used_inputs_file`; read back by proc0), or the config dump of `gs2_config_type%write_to_unit` if diagnostics are skipped |
| `input_file_text` | raw text of the input file (proc0 reads it) |

Root dimensions (`neasyf_dim`; coordinate variables in `/` where values are given):

| dim | size | coordinate values |
|---|---|---|
| `ri` | 2 | 1..2 |
| `theta` | `2*ntgrid+1` | `theta_grid::theta` |
| `theta_m1` | `2*ntgrid` | 1..n |
| `sign` | 2 | 1..2 |
| `kx` | `ntheta0` | `akx` |
| `ky` | `naky` | `aky` |
| `lambda` | `nlambda` | `al` |
| `egrid`, `egrid_p1`, `egrid_m1` | `negrid`, `negrid+1`, `max(negrid-1,1)` | 1..n |
| `species` | `nspec` | 1..n |
| `xi2` | `2*nlambda` | 1..n |
| `xi_lim`, `xi_p1` | `nxi_lim`, `nxi+1` | 1..n (only if collisions are initialised) |
| `bset` | `nbset` (if > 0) | `bset` |
| `two`, `three` | 2, 3 | |
| `glo` | `g_lo%ulim_world+1` | 0-based iglo (only with `write_raw_glo`) |
| `nlinks_max` | `size(g_adj,1)` | (only if `g_adj` exported) |
| `n_class`, `sc_Mmax`, `sc_Nmax` | `i_class`, `maxval(M_class)`, `maxval(N_class)` | |
| `sc_R_c<c>`, `sc_M_c<c>`, `sc_RM_c<c>` | `nidx*N_class(c)`, `M_class(c)`, `nidx*N_class(c)*M_class(c)` | per class `c` |
| `lr_ax_ik<ik>_is<is>` | local-fields matrix size | per supercell (field_option local only) |

Groups inherit root dimensions (netCDF-4 searches ancestor groups in `nc_inq_dimid`; confirmed in
Appendix A.1), so `neasyf_write(gid, …, dim_names=[…])` works with a group id `gid` from
`nf90_def_grp(ncid, name, gid)`.

#### 6.1.2 Variable attributes (every variable)

`fortran_name` (e.g. `"dist_fn::ainv"`), `fortran_decl` (e.g. `"complex(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc)"`),
`canonical` (shape tag of §0.2), `fortran_type` (`real|complex|integer|logical`), `fortran_lbounds`,
`source` (file:line of the computing routine, as in §3), `dt_dependent` (`0/1`), `code_dt` (the `code_dt`
at export time, for `dt_dependent=1`), `index_base` (index-valued integers), `owner_spec` (e.g. `"03"`).
Arrays that are **not allocated** in GS2 are **not written**; the group attribute `absent` lists their
names (space separated) so the loader can distinguish "absent" from "forgotten".

#### 6.1.3 Groups and contents

| group | contents (definition and canonical tag in §3) |
|---|---|
| `/theta_grid` | §3.1 |
| `/kt_grids` | §3.2 (`kwork_filter` after `init_fields`) |
| `/le_grids` | §3.3 incl. derived bounce/wfb/hybrid masks |
| `/species` | §3.4 (`spec_<component>` vectors) |
| `/run_parameters`, `/gs2_time`, `/gs2_reinit`, `/init_g`, `/fields` | §3.5 (scalars as 0-D variables; strings as attributes) |
| `/gs2_layouts` | §3.6 incl. `iglo_of`, `ile_of` |
| `/dist_fn_arrays` | §3.7 static arrays (`vpa, vpac, vperp2, aj0, aj1, wstar`) |
| `/dist_fn` | §3.8 coefficients, topology, effective knobs, derived wfb factors, `itleftmost/itrightmost`, `connections_*`, `save_h` |
| `/fields_implicit` | `aminv_c<c>`, `cls_ik`, `cls_it`, `sc_class`, `sc_member`, `sc_cell`, `nfield`, `nidx`, optional `am_c<c>` |
| `/fields_local` | `resp_ik<ik>_is<is>`, with attributes `ik`, `is`, `it_leftmost`, `condition`, `dt` (field_option local only) |
| `/collisions` | §3.10 (LE6/LZ6/E7 arrays, vectors, scalars, helpers `e15`, `capg_half`) — only if `collision_model_switch /= none` and `colls` |
| `/hyper` | §3.11 |
| `/nonlinear_terms` | §3.12 |
| `/diagnostics` | parsed `gnostics` scalars/flags (spec 08 §2.3; NEW_DIAG; only if diagnostics initialised), `igomega`, `nwrite`, `navg`, `omegatol`, `omegatinst`, `exit_when_converged`, all `write_*` flags |
| `/state_t0_prediag` | time block + `g, gnew, phi, apar, bpar, phinew, aparnew, bparnew, gexp_1..3, kx_shift, theta0_shift, vnmult`, `exb_first`, `nl_istep_last`, `nl_initializing`, `reset` — immediately after `initialize_equations` |
| `/state_t0` | the same set after `initialize_diagnostics` (or identical to the previous group when diagnostics are skipped). **This is the JAX initial state.** |
| `/gates` | one int8 0-D variable per gate flag of §8.3 and of the other specs' gate lists (07 App. B, 01 §4, 05 §6, 06 §6, 08 §5) |
| `/dt_ladder/k<k>/<module>` | §6.3 |
| `/raw` | optional glo-ordered copies (`gnew_glo (ri, theta, sign, glo)`, `iglo_kykxsel` style maps) for debugging |

Naming rule: inside a group the variable name is the Fortran name spelled exactly as in its declaration
(e.g. `D_res`, `Rplot`, `vnew_D`; Fortran is case-insensitive, netCDF is not). Names that
the other specs suggested in flat form map to group paths as follows: `hyper_D_res` → `/hyper/D_res`,
`hyper_filter_s3d` → `/hyper/hypervisc_filter`, `spec_z` → `/species/spec_z`, `aminv_c<c>` →
`/fields_implicit/aminv_c<c>`, `cls_ik/cls_it` → `/fields_implicit/cls_ik|cls_it`, `wfb_aj0l…` →
`/dist_fn/wfb_aj0l…`, `vperp2left/right` → `/dist_fn/wfb_vperp2left|right`, `neighbor` →
`/dist_fn/connections_neighbor`, collisions names → `/collisions/<name>`.

### 6.2 Trace file `<run_name>.gs2trace.nc`

Same global attributes and root dimensions as the export file (plus `nstep_trace`, `trace_stages`).
One group per recorded probe:

```
/s<istep:06d>/a<attempt>/p<seq:03d>_<tag>        attributes: tag, istep, attempt, seq, pass (probe_pass),
                                                  in_advance (0/1)
   time block (0-D variables), reset (int8), plus the arrays of the §5.3 table for that tag
/final                                            snapshot after evolve_equations returned (after the forced
                                                  final diagnostics call, gs2_main.fpp:869-882)
```

Group names start with a letter (portable netCDF names). Consistency invariants the loader can assert:
`p…_step_end` of step k equals `p…_advance_begin` of step k+1 attempt 1 (same arrays, bitwise);
`after_timeadv` pass 2 of the last attempt equals `advance_end` when collisions are not split.

### 6.3 dt ladder (`dt_ladder_n > 0`)

For `k = 1..dt_ladder_n` the exporter re-initialises the dt-dependent levels exactly as `reset_time_step`
does (§7.8) and writes, under `/dt_ladder/k<k>/` (attribute `code_dt`), the dt-dependent arrays:
`/dist_fn/{vpar, wdrift, wdriftttp, a, b, r, ainv, g_h, source_coeffs_*}`, `/dist_fn_arrays/wstar`,
`/collisions/{c1le, betaale, qle, ec1le, ebetaale, eqle, z0le, s0le, w0le, bz0le, bs0le, bw0le}` (+ the
standard-layout ones), `/fields_implicit/aminv_c<c>` (or `/fields_local/resp_*`), `/hyper/D_res` and
`/hyper/hypervisc_filter` (simple3D; these stay at their initial-dt values, which the ladder must reproduce:
`init_hyper` is guarded by `initialized` and level 19 is not re-entered, spec 01 §1.2/§2.7), `/nonlinear_terms/{cflx, cfly}`.
This implements spec 01 §5.5 ("exported dt ladder"): JAX looks up the coefficients by the bit pattern of
`code_dt` and falls back to RECOMPUTE on a miss. The ladder run does not step afterwards.

### 6.4 GS2 restart file format (for JAX → GS2 and GS2 → JAX round trips)

With this build (`save_many`/`read_many` irrelevant because `has_netcdf_parallel = .false.`,
gs2_save.fpp:73-79) there is **one file per rank**: `trim(restart_file)//fileopt//'.'//iproc`
(get_file_proc 1248-1257; `restart_file` = `restart_dir//restart_file` from `init_g_knobs`,
init_g.f90:977-1036; default `"restart_file_not_set.nc"`, init_g.f90:314). On nproc = 1: `<restart>.0`.
Contents (gs2_save_for_restart 86-469):

| name | netCDF dims (Fortran order) | type | notes |
|---|---|---|---|
| `glo` | dim, size = local `ulim_proc-llim_proc+1` (unlimited if 0) | | **glo order of the writer's layout** |
| `theta`, `sign`, `aky`, `akx` | dims with coordinate values (`theta`, —, `aky`, `akx`) | f64 | note `aky/akx` names (not `ky/kx`) |
| `delt0, delt1, delt2, delt_max, t0` | scalar | f64 | `code_dt, code_dt_prev1, code_dt_prev2, code_dt_max, user_time` |
| `layout` | char, dim `layout_dim` (len 5) | char | restore aborts on mismatch (767-771) |
| `ntheta0, naky, nlambda, negrid, nspec` | scalar | i32 | checked loosely (warnings only) |
| `current_scan_parameter_value`, `vnm1`, `vnm2` | scalar | f64 | `vnmult(1:2)` |
| `gr`, `gi` | `(theta, sign, glo)` | f64 | real/imag parts of `g` (the caller passes `gnew`) |
| `gexp_1_r … gexp_3_i` | `(theta, sign, glo)` | f64 | iff allocated and `include_explicit_source_in_restart` (default `.true.`) |
| `phi_r, phi_i, apar_r, …, bpar_i` | `(theta, akx, aky)` | f64 | from `phinew…`, only present fields |
| `kx_shift, theta0_shift` | `(aky)` | f64 | zeros if unallocated |
| `a_ant_r/_i, b_ant_r/_i` | `(nk_stir)` | f64 | antenna only (GATED) |
| `lambda, species, egrid, energy, iglo_llim_proc, iglo_kykxsel(aky,akx,species,egrid,lambda), procid_kykxsel, jtwist` | | | iff `save_glo_info_and_grids` |
| `vpa (theta,sign,glo)`, `vperp2 (theta,glo)` | | f64 | iff `save_velocities` |

JAX → GS2: write `<restart>.0` with the glo axis in the **iglo order of the target layout**
(`full(:, :, idx(g_lo, ik,it,il,ie,is)) = g[is,ie,il,ik,it,:,:]`), set `ginit_option = 'many'`,
`scale = 1.0` (`g*scale` is exact for 1.0, gs2_save.fpp:865-870), `delt_option = 'check_restart'` (reads
`delt0..delt_max` via init_dt) and **`force_maxwell_reinit = .false.`** if the fields in the file must be
used as-is (otherwise `set_init_fields` recomputes them from `g`, gs2_init.fpp:318, which is not bitwise the
original fields after diagnostics round trips). GS2 → JAX: use the export file instead; the restart file is
only a cross-check of the glo↔canonical map (`iglo_kykxsel` gives `iglo` per `(ik,it,is,ie,il)`).

### 6.5 Response matrix files (`gs2_save_response`, gs2_save.fpp:892-969)

`<response_file>_ik_<ik>_is_<is>_dt_<code_dt in E14.6E3, trimmed><suffix>` (fields_arrays.f90:38-57;
default suffix `.response`). Dims `ri(2)`, `ax1`, `ax2` (square); variables `dt` (f64 scalar),
`condition` (f64, −1 if not given), `response(ri, ax1, ax2)` = `convert::c2r(resp)`; NumPy
`response[ax2, ax1, ri]`. Implicit (`dump_response_to_file_imp`, fields_implicit.f90:1172-1404) drops the
duplicate boundary points (lossless, spec 04 §3.15); local (`sc_dump_to_file`, fields_local.f90:987-1017)
writes `self%nrow × self%ncol` from `pull_rows_to_arr`.

### 6.6 Python loader contract (`gpu/gs2gpu/io/export.py`, to be written)

* Open with `netCDF4.Dataset(path)`; call `ds.set_auto_mask(False)` (no masked arrays, no fill-value
  surprises); access by path `ds['/group/var']`.
* `fortran_type == "complex"` → `a[..., 0] + 1j*a[..., 1]` as `complex128`; `"logical"` → `a != 0`; reals
  stay `float64` (never cast to float32); integers stay `int32` with their Fortran values.
* Expose `dims = {ntgrid, nth, ntheta0, naky, nlambda, negrid, nspec, ng2, nxi, nxi_lim, nfield, …}` from
  the root dimensions; verify `export_format_version` major == 1 and refuse if `gate_status != "ok"` unless
  the caller explicitly overrides.
* Missing variable ⇒ check the group's `absent` attribute; anything not listed there is an error.

---------------------------------------------------------------------------------------------------

## 7. Exact pseudo-code of the exporter

### 7.1 Files and modules (all new)

| file | module / program | uses | role |
|---|---|---|---|
| `src/gs2_export_probe.f90` | `gs2_export_probe` | nothing | §5.1 (used by physics modules) |
| `src/gs2_export_writer.fpp` | `gs2_export_writer` | `neasyf`, `netcdf`, `mp`, `gs2_layouts`, `convert`, `iso_fortran_env` | file/group handling, attributes, canonicalisers (§2.4), gathers (§2.5), time block |
| `src/gs2_export_static.fpp` | `gs2_export_static` | writer + every physics module of §3 | `export_static`, `export_state`, `export_dt_dependent`, derived arrays (§7.7) |
| `src/gs2_export_trace.fpp` | `gs2_export_trace` | writer, static (for `export_state`), probe | hook implementations (§5.3), trace file |
| `src/programs/gs2_export.fpp` | program `gs2_export` | `gs2_main`, the three modules above | control flow §1.4, knobs §8.1, gates §8.3 |

The last line of every file must be `end module <name>` / `end program gs2_export` with a single space
after `end`: `scripts/fortdep` discovers modules and programs only through these lines (regexes
`^[\t ]*end module[\t\n ]`, `^[\t ]*end program[\t\n ]`, scripts/fortdep `get_fileinfo`), and
dependencies through `use <name>` lines (`use, intrinsic :: …` is ignored, which is correct).

### 7.2 Writer API (all routines collective unless stated; only proc0 calls netCDF)

```fortran
integer function ex_open(fname, title, header)           ! proc0: neasyf_open(fname,"w"); create_metadata; returns ncid (-1 elsewhere)
subroutine ex_close(ncid)
integer function ex_group(parent, name)                  ! proc0: nf90_inq_ncid or nf90_def_grp; neasyf_error on failure
subroutine ex_root_dims(ncid)                            ! §6.1.1
subroutine ex_attr_text(gid, varname_or_global, name, value)
subroutine ex_attr_int / ex_attr_real
subroutine ex_scalar(gid, name, value, ...)              ! real/int/logical(int8) 0-D, plus §6.1.2 attributes
subroutine ex_vec(gid, name, values, dim_names, ...)     ! any V-type array (ranks 1-4, real/int/logical/complex)
subroutine ex_g7c(gid, name, arr, meta)                  ! complex(-ntgrid:,:,g_lo%llim_proc:)   → G7 + ri
subroutine ex_g7r(gid, name, arr, meta)                  ! real G7
subroutine ex_g6r(gid, name, arr, meta)                  ! real(-ntgrid:, g_lo%llim_proc:)
subroutine ex_g7m1c(gid, name, arr, meta)                ! complex(-ntgrid:ntgrid-1, 2, g_lo)
subroutine ex_glc(gid, name, arr, meta)                  ! complex(n_links_max, 2, g_lo)
subroutine ex_g5i / ex_g5l / ex_g5sl                     ! integer/logical per iglo (and (2,iglo))
subroutine ex_f3c(gid, name, arr, meta) / ex_f3r         ! (-ntgrid:ntgrid, ntheta0, naky)
subroutine ex_le6c / ex_le6r (gid, name, arr, nxdim, nedim, meta)   ! (:, :, le_lo%llim_proc:)
subroutine ex_lz6r / ex_e7r                              ! use_le_layout = .false. only
subroutine gather_glo_c / gather_glo_r / gather_le_c / gather_le_r / gather_jf_rows   ! §2.5
```

`meta` is a small derived type `(fortran_name, fortran_decl, source, owner_spec, dt_dependent)`; the
writer adds `canonical`, `fortran_type`, `fortran_lbounds`, `code_dt`. Compression:
`neasyf_default_compression = knobs%compression` (default 0; shuffle+deflate are lossless when enabled).

### 7.3 `export_static` — order of groups (deterministic; every rank executes the same sequence)

```
ex_root_dims
/theta_grid      §3.1     /kt_grids  §3.2     /le_grids §3.3 (+ derived masks §7.7)
/species         §3.4     /run_parameters, /gs2_time, /gs2_reinit, /init_g, /fields  §3.5
/gs2_layouts     §3.6 (iglo_of, ile_of: loops over (it,ik,il,ie,is) resp. (ig,it,ik,is) calling idx)
/dist_fn_arrays  §3.7 static part
/dist_fn         §3.8 (coefficients via ex_g7*, topology, knobs, derived §7.7)
/fields_implicit §7.4    or /fields_local §7.5   (by fieldopt_switch; neither for 'test'/'gf_local' — gated)
/collisions      §3.10   (skip if collision_model_switch == collision_model_none or .not. colls; record `absent`)
/hyper           §3.11   /nonlinear_terms §3.12   /diagnostics (if initialised later: written in export_config_text)
/gates           §8.3
```

For every allocatable array: `if (allocated(x)) then write else append name to the group's absent list`.
`g_h`, `a,b,r,ainv`, `vpar`, `wdrift(ttp)` are always allocated after init; `source_coeffs_*` iff
`opt_source` (`aparm/aparp` iff also `has_apar`), `gamtot3` iff fieldlineavg, `itleft/itright/connections/
save_h/g_adj` iff linked (3528-3534), `kx_shift`/`theta0_shift` iff flow shear, `gexp_*` iff nonlinear.

### 7.4 `aminv` export (field_option implicit)

```
nf = jf_lo%nfield ; nidx = (2*ntgrid+1)*nf            # = fields_implicit::nidx
for c = 1..i_class:
    R = nidx*N_class(c); M = M_class(c)
    buf(1:2, 1:R, 1:R, 1:M) = 0                         # Fortran (ri, col, row, m); proc0 only
    for m = 1..M:
        ik = f_lo(c)%ik(m,1)
        for n = 1..N_class(c):                           # output cell
            it = f_lo(c)%it(m,n)
            for f = 1..nf: for ig = -ntgrid..ntgrid:
                jflo = ij_idx(jf_lo, ig, f, ik, it)      # gs2_layouts.fpp:1967-1972
                row  = 1 + (ig+ntgrid) + nth*(f-1) + nidx*(n-1)     # B order (1-based)
                if idx_local(jf_lo, jflo):
                    vec = aminv(c)%dcell(jf_lo%dj(c, jflo))%supercell(1:R)
                    if not proc0: send(vec, 0)            # (nproc > 1 only, fixed loop order on all ranks)
                elif proc0: receive(vec, proc_id(jf_lo, jflo))
                if proc0: buf(1,:,row,m) = real(vec); buf(2,:,row,m) = aimag(vec)
    write /fields_implicit/aminv_c<c> (ri, sc_R_c<c>, sc_R_c<c>, sc_M_c<c>)   # NumPy [m, row, col, ri]
cls_ik(m,n,c) = f_lo(c)%ik(m,n), cls_it(m,n,c) = f_lo(c)%it(m,n) for m<=M_class(c), n<=N_class(c); −1 elsewhere
sc_class(it,ik), sc_member(it,ik), sc_cell(it,ik): inverse of cls_* (every (it,ik) appears exactly once)
```

This is identical to the loop of `dump_response_to_file_imp` (fields_implicit.f90:1251-1399) without its
duplicate-point removal, i.e. the raw runtime data consumed by `getfield` (222-312). The rows of the
`aky = akx = 0` mode are zero (`a_inv(:,jlo) = 0.0`, fields_implicit.f90:855-857) and the internal
`ig = ntgrid` rows of linked cells are copies (`ilo = jlo + nfield`, fields_implicit.f90:882-907); export them as stored.

### 7.5 Local-fields response export (field_option local)

```
saved = fields_arrays::response_file
fields_arrays::response_file = trim(scratch_prefix)            # e.g. trim(run_name)//'.gs2export_tmp'
call fields::dump_response_to_file(suffix='.resp')             # fields.f90:380-397 → dump_response_to_file_local
for ik = 1..naky: for is = 1..n_supercells(ik):                # dist_fn::n_supercells (public)
    fname = get_specific_response_file_name(ik, is, code_dt, '.resp')   # fields_arrays.f90:38-57
    proc0: nf90_open; inquire ax1, ax2; neasyf_read 'response' (ri,ax1,ax2), 'dt', 'condition'
           write /fields_local/resp_ik<ik>_is<is> (ri, lr_ax_ik<ik>_is<is>, same) + attributes ik, is, dt, condition
           delete the temporary file
fields_arrays::response_file = saved
```

The ky = kx = 0 supercell may have no file (fields_local.f90:1047-1057 comment); record it in `absent`.
The row/column order is fields_local's extended-domain order (spec 04 §7 owns its interpretation); the
exporter additionally writes `/fields_local/it_leftmost_ik<ik>_is<is>` from `supercell_labels`/
`get_leftmost_it` so the loader can identify the supercell.

### 7.6 Time block and state export

Time block (0-D f64 variables unless noted), read from `gs2_time` (public) and friends:
`code_dt, user_dt, code_dt_prev1, code_dt_prev2, code_dt_old, code_dt_cfl, user_dt_cfl, code_dt_min,
code_dt_max, code_time, user_time`, `istep` (i32; −1 outside stepping), `reset` (int8,
`run_parameters::reset`), `vnmult(two)`, `nl_istep_last` (i32), `nl_initializing` (int8), `exb_first`
(int8). These are bit-exact copies; the AB coefficients are **not** exported (RECOMPUTE, spec 01 §3.10).

```
subroutine export_state(gname)
    gid = ex_group(ncid, gname); write time block
    ex_g7c(gid,'g',g); ex_g7c(gid,'gnew',gnew)
    ex_f3c(gid,'phi',phi) … 'bparnew'          # all six, present or not; attribute present=has_*
    if allocated(gexp_1): ex_g7c 'gexp_1', 'gexp_2', 'gexp_3'
    if allocated(kx_shift): ex_vec 'kx_shift' (ky); if allocated(theta0_shift): ex_vec 'theta0_shift'
    ex_vec 'vnmult' (two)
```

### 7.7 Derived arrays (exporter-computed; verbatim expressions)

Compute on proc0 (they need only replicated data) with the **same expression text** as GS2, compiled with
the same flags; where the GS2 expression contains no `a*b+c` pattern, FMA contraction cannot differ.

```
# le_grids masks (public elemental functions)
is_lower_bounce_point(ig,il), is_upper_bounce_point(ig,il)   for ig, il   → int8 (theta, lambda)
il_is_wfb(il), il_is_passing(il), il_is_trapped(il)          for il       → int8 (lambda)
is_passing_hybrid_electron(is,ik,il)                         → int8 (lambda, ky, species)

# dist_fn topology helpers
itleftmost(it,ik)  = get_leftmost_it(it,ik) ; itrightmost(it,ik) = get_rightmost_it(it,ik)    (3889-3920)
connections_neighbor(iglo)  = connections(iglo)%neighbor                         → G5 int8
connections_it_left(iglo)   = it_idx(g_lo, connections(iglo)%iglo_left)  if iglo_left >= 0 else −1  → G5 int32
connections_it_right(iglo)  likewise                                              (2518-2565)
save_h(isgn, iglo)                                                                → G5s int8

# mixed-wfb end factors (only if boundary linked .and. mixed_wfb .and. nonad_zero), il = ng2+1,
# dist_fn.fpp:6016-6028, `use spfunc, only: j0, j1`:
for is, ie:
   wfb_vperp2left(ie,is)  = bmag(-ntgrid)*al(il)*energy(ie,is)
   wfb_vperp2right(ie,is) = bmag(ntgrid)*al(il)*energy(ie,is)
   for ik, it:
      itl = get_leftmost_it(it,ik); itr = get_rightmost_it(it,ik)
      argl = spec(is)%bess_fac*spec(is)%smz*sqrt(energy(ie,is)*al(il)/bmag(-ntgrid)*kperp2(-ntgrid,itl,ik))
      argr = spec(is)%bess_fac*spec(is)%smz*sqrt(energy(ie,is)*al(il)/bmag(ntgrid)*kperp2(ntgrid,itr,ik))
      wfb_aj0l(ie,it,ik,is) = j0(argl); wfb_aj0r = j0(argr); wfb_aj1l = j1(argl); wfb_aj1r = j1(argr)
   → (egrid, kx, ky, species) and (egrid, species)

# theta_grid
field_line_average_weight (public after §4.2); cross-check value = 1.0/trapezoidal_integration(theta, jacob)
   (integration.f90; theta_grid.f90:520) — must be bitwise equal

# collisions helpers (spec 06 §3.4)
e15(ie) = energy_maxw(ie)**1.5
capg_half(j) = capg((speed_maxw(j)+speed_maxw(j+1))*0.5),  j = 1..negrid-1, with capg/capg_kernel copied
   verbatim from collisions.fpp:2165-2176:
      capg_kernel(xe) = erf(xe) - 2 * xe * exp(-xe**2) / sqrt_pi
      capg(xe)        = 2 * exp(-xe**2) * capg_kernel(xe) / (xe * sqrt_pi)
   (and capg_old, 2178-2181, if ediff_switch selects the old scheme)
```

### 7.8 dt ladder (mirrors `reset_time_step`, gs2_reinit.f90:104-225, without state preservation)

```
use gs2_init, only: init, init_level_list ; use gs2_reinit, only: reduce_time_step
for k = 1..dt_ladder_n:
    call init(state%init, init_level_list%override_timestep)   # down handlers of levels 35..30 (finish fields_level_2,
                                                                 # dist_fn_level_3, collisions, …), gs2_init.fpp:237-269
    call reduce_time_step                                       # code_dt = code_dt/delt_adj  (gs2_reinit.f90:86-91)
    if (code_dt <= code_dt_min) then record and exit
    call save_dt(code_dt)                                       # gs2_time.f90:206-213
    call init(state%init, init_level_list%full)                # up handlers 30..35 (response matrix rebuilt)
    export_dt_dependent('/dt_ladder/k<k>')                      # §6.3 list; attribute code_dt
```

Level 19 (`hyper`) is not re-entered, so `D_res` / simple3D `hypervisc_filter` stay at the initial-dt values,
exactly as in a GS2 run that resets dt (spec 01 §2.7). `set_initial_values` (level 34) re-runs `ginit` with no
override, which changes `gnew` — irrelevant because the ladder run never steps or exports state afterwards.

### 7.9 Trace hooks (module `gs2_export_trace`)

```
subroutine tr_point(tag, istep)
    if (.not. selected_tag(tag)) return                                  # knobs, identical on all ranks
    if (inner_tag(tag) .and. .not. probe_in_advance) return
    if (tag == 'advance_begin') then
        attempt = merge(attempt+1, 1, istep == last_istep); last_istep = istep; seq = 0
    if (.not. selected_step(istep)) return                               # trace_first_step / trace_every
    seq = seq + 1
    gid = tr_group(istep, attempt, seq, tag)                             # proc0: /s<istep>/a<attempt>/p<seq>_<tag>
    write time block; dump arrays of §5.3 for `tag` (collective writers)
subroutine tr_g(tag, istep, arr) / tr_le(tag, istep, arr) / tr_real / tr_c2: same bookkeeping, then
    write the passed array with the canonicaliser that matches its layout (index remap of §5.1)
```

`inner_tag` = every tag of §5.2 between `after_exb_shear` and `after_collisions_advance` plus the `nl_*`
and `coll_*` tags. `tr_c2` (P-R1) stores `am` per class in memory (allocatable list); `export_static` (or the
reset record) writes it as `am_c<c>` `(ri, sc_R_c<c>, sc_RM_c<c>)`: Fortran `am(k, ilo)` with output `k` in B
order and probe column `ilo = idx(f_lo(c), if, im)` in K order (spec 04 §3.12).

---------------------------------------------------------------------------------------------------

## 8. Configuration options

### 8.1 New namelist `&gs2_export_knobs` (read by `gs2_export` only)

Source: the file named by command argument 2 if present (`gs2_export run.in export.nml`), otherwise the
`&gs2_export_knobs` group of the GS2 input file, read on proc0 with `file_utils::input_unit_exist`
exactly like `dump_response_knobs` (dump_response.f90:57-61), then broadcast. GS2 itself ignores unknown
namelist groups (no validation of foreign groups exists in `src/` or `externals/utils/`; checked).
`write_used_inputs_file` does not echo this group (it writes only GS2 config types).

| option | default | effect | decision |
|---|---|---|---|
| `export_file` | `''` → `trim(run_name)//'.gs2export.nc'` | export file name | SUPPORT |
| `write_response` | `.true.` | write `aminv_c*` / `resp_*` | SUPPORT |
| `write_collisions` | `.true.` | write `/collisions` | SUPPORT |
| `write_g_adj` | `.false.` | also write `g_adj` (GL) | SUPPORT (debug) |
| `write_raw_glo` | `.false.` | write glo-ordered copies of the state arrays in `/raw` | SUPPORT (debug) |
| `write_restart_t0` | `.false.` | `gs2_save_for_restart(…, fileopt='.gs2export_t0', save_glo_info_and_grids=.true., save_velocities=.true.)` after `state_t0` | SUPPORT (cross-check §12 T3) |
| `run_diagnostics` | `.true.` | call `initialize_diagnostics` (istep-0 round trips) and run per-step diagnostics while tracing, i.e. mirror `bin/gs2` exactly | SUPPORT |
| `nstep_trace` | `0` | number of Fortran steps to run with probes (`≤ nstep` of `&knobs`) | SUPPORT |
| `trace_file` | `''` → `trim(run_name)//'.gs2trace.nc'` | trace file name | SUPPORT |
| `trace_stages` | `'all'` | space-separated tag list of §5.2, or `'all'`, or `'step'` (= `advance_begin advance_end step_end`) | SUPPORT |
| `trace_first_step`, `trace_every` | `1`, `1` | record steps `istep >= first` with `mod(istep-first, every) == 0` | SUPPORT |
| `trace_response_build` | `.false.` | install hooks before `initialize_equations` and capture P-R1 `am` (init and resets) | SUPPORT (04 Q10) |
| `dt_ladder_n` | `0` | §7.8; incompatible with `nstep_trace > 0` (abort) | SUPPORT |
| `allow_multi_rank` | `.false.` | permit `nproc > 1` (results flagged rank-dependent, §2.5) | SUPPORT |
| `strict` | `.true.` | abort on hard gates (§8.3); soft gates are always only recorded | SUPPORT |
| `compression` | `0` | `neasyf_default_compression` for export/trace files (lossless) | SUPPORT |

### 8.2 Existing GS2 options that change the exporter's behaviour or output

| option | namelist | default | effect on export | decision | reason |
|---|---|---|---|---|---|
| `layout` | layouts_knobs | `'lxyes'` | glo order, summation order; canonical output is layout-independent, `iglo_of` is not | SUPPORT all six | exported explicitly |
| `local_field_solve` | layouts_knobs | `.false.` | `f_lo` blocksize (gs2_layouts.fpp:1678-1683) and broadcast path in `init_inverse_matrix` | SUPPORT | identical numerics on 1 rank |
| `intmom_sub`, `intspec_sub` | layouts_knobs | `.false.` | sub-communicator integrals | SUPPORT | trivial on 1 rank |
| `gf_local_fields` | layouts_knobs | `.false.` | forced `.false.` if `nproc < naky*ntheta0` (1007-1015) | SUPPORT (inert on 1 rank) | |
| `nproc_g_lo`, `nproc_le_lo`, … | layouts_knobs | 0 | fewer active ranks per layout | SUPPORT | gathers use each rank's own bounds |
| `opt_redist_*`, `opt_local_copy`, `unbalanced_*`, `max_unbalanced_*` | layouts_knobs | | communication only | SUPPORT | no numerics |
| `fft_measure_plan` | layouts_knobs | `.false.` | FFTW_MEASURE plans are not reproducible run-to-run | GATE (`.true.` → soft gate) | Pitfall 11.14 |
| `fft_use_wisdom`, `fft_wisdom_file` | layouts_knobs | `.true.`, `'default'` → `<run>.fftw_wisdom` or `$GK_FFTW3_WISDOM` (630-645) | wisdom from MEASURE runs is reused by ESTIMATE planning | SUPPORT, record in `/gs2_layouts` | Pitfall 11.14 |
| `field_option` | fields_knobs | `'default'` (= implicit) | which response is exported | SUPPORT implicit, local; GATE `gf_local`, `test` | 04 §4 |
| `dump_response`, `read_response`, `response_dir`, `response_file` | fields_knobs | `.false.` | `aminv` read from files instead of built (fields_implicit.f90:488-497) | SUPPORT | the exported `aminv` is what is used |
| `ginit_option`, `restart_file`, `restart_dir`, `read_many`, `scale` | init_g_knobs | | initial state source | SUPPORT | `state_t0` captures the result |
| `include_explicit_source_in_restart` | init_g / gs2_save | `.true.` | `gexp_*` from restart | SUPPORT | captured |
| `delt_option` | knobs | `'default'` | `'check_restart'` loads `delt0..delt_max`, vnmult | SUPPORT | time block captures |
| `nsave`, `save_many` | gs2_diagnostics_knobs | | restart files written during the trace | SUPPORT | no effect on numerics |
| `nwrite`, `write_moments`, `write_fluxes`, … | gs2_diagnostics_knobs | | g_adjust round trips during `run_diagnostics` (08) | SUPPORT | mirrored when `run_diagnostics=.true.` |
| `use_old_diagnostics` | knobs | `.false.` | old diagnostics instead of new | SUPPORT (recorded) | 08 owns |
| `do_eigsolve` | knobs | `.false.` | no time loop | GATE (hard) | not a time stepper |
| optimisation mode (`optimisation_knobs`) / `state%init%opt_ov%override_nproc` | | off | `run_optimise_gs2`, `use_nproc` | GATE (hard) | exporter never calls `run_optimise_gs2` |
| list mode (`*.list`), `nensembles > 1`, Trinity | command line / state | off | job forking | GATE (hard) | `initialize_gs2` would fork |

### 8.3 Gates

**Hard (abort before export unless noted; `strict=.false.` turns the starred ones into soft gates):**
`gs2_is_single_precision` or `kind(1.0) /= real64`; `-DSHMEM` build (g_lo blocks, `gexp_1` pointer,
proc_id via `proc_map`, gs2_layouts.fpp:866-905, 1583-1603); `gs2_has_lowflow` (`wdrift` complex);
`gs2_has_openmp .and. omp_get_max_threads() > 1`*; `nproc > 1 .and. .not. allow_multi_rank`*;
`do_eigsolve`; list/ensemble/Trinity run; `dt_ladder_n > 0 .and. nstep_trace > 0`; `nstep_trace > nstep`
(clamped with a warning instead).

**Soft (recorded in `/gates` and `gate_status`; the JAX runner refuses unless overridden):** the union of
the other specs' gate lists, evaluated after `initialize_equations` from exported effective values:
antenna on (`.not. no_driver`), `any(ieqzip)`, `afilter > epsilon(0.)`, `split_nonlinear`/`strang_split`,
`field_option` in {`gf_local`, `test`}, `gf_lo_integrate`, `fft_measure_plan`, `new_field_init = .false.`,
`wstar_units .and. (nonlin .or. any(aky == 0))`, flow shear with the kx branch and `ntheta0 == 1`,
hyper normaliser zero cases (07 App. B), `scale /= 1` with file-based reset (01 §6.12), hybrid or
non-Maxwellian species with collisions (07 Q6/Q7), `use_le_layout = .false.` if spec 06's standard-layout
path is not ported, plus any further gate a spec adds later (the `/gates` group is open-ended).

---------------------------------------------------------------------------------------------------

## 9. Build integration (exact edits)

### 9.1 How the build finds programs and modules

* Sources are discovered by `scripts/fortdep` over `VPATH` (`Makefile:551-554`: `SUBDIRS` = `src/geo`,
  `src/diagnostics`, `externals/utils`, plus `src`, `src/programs`, `generated_src`, `externals/neasyf/src`).
  New files in `src/` or `src/programs/` need **no** `Makefile` change.
* `make depend` (`Makefile:653-660`) runs `perl scripts/fortdep -m "$(MAKE) --no-print-directory" -1 -o -v=0
  $(VPATH) --generated_src_dir=$(GENERATED_SRC_DIR)`: it preprocesses `.fpp` files into `generated_src/`
  (by calling make), parses `use`/`end module`/`end program`/`submodule` lines, and writes
  `Makefile.depend` with one `<program>_mod = $(OBJDIR)/….o …` list per program (keyed by the program
  name) and one `$(OBJDIR)/x.o: deps` rule per file. `-o` overwrites, `-1` takes the first file for a
  module defined twice.
* `Makefile.depend` is git-ignored (`.gitignore:5`) and is regenerated automatically **only if it does
  not exist** (`Makefile:645-651`). After adding files or `use` statements, run `make depend` explicitly
  (the current file was generated on Oct 21 19:32 and knows nothing about the new modules).
* `Makefile.target_gs2` turns `<prog>_mod` into a link rule with `build_program` (`Makefile.utils:69-76`):
  `$(BINDIR)/<prog>: $(<prog>_mod)` linked with `$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)`.
* Programs that (transitively) use `git_version_mod` need the non-module object `git_hash.o` prepended via
  `needs_git_hash` (Makefile.target_gs2:21-31; fortdep cannot see it). `gs2_export` uses `gs2_main` and
  `gs2_metadata`, hence needs it.

### 9.2 Exact edits

`Makefile.target_gs2`:

```make
# line 29: append gs2_export
needs_git_hash := gs2 multigs2 ingen write_default_input_file dump_grids dump_response generate_fftw_wisdom ffttester box_to_ballooning_space gs2_export
```

and after line 39 (the last `parallel_targets +=` line):

```make
# GPU-port tooling: built only on request (`make gs2_export`), deliberately not part of gs2_all so that
# default builds, `make tests` and CI keep their current target set.
gpu_port_targets := gs2_export #Dumps GS2 grids/coefficients/state (+ optional stage traces) to netCDF for the JAX port
```

and after line 45 (`$(foreach prog,$(serial_targets),…)`):

```make
$(foreach prog,$(gpu_port_targets),$(eval $(call build_program,$(prog),$(LDFLAGS))))
```

No change to `Makefile`, `Makefile.utils`, `Makefiles/*`. Build sequence (use the same variables as for
the existing `bin/gs2`, Q13):

```
make GK_SYSTEM=macosx-homebrew <same overrides as bin/gs2> depend
make GK_SYSTEM=macosx-homebrew <same overrides as bin/gs2> gs2 gs2_export     # → bin/gs2, bin/gs2_export
mpirun -np 1 bin/gs2_export case.in [export.nml]                              # or plain bin/gs2_export case.in
```

Because the probe edits add `use gs2_export_probe` to `dist_fn`, `fields`, `fields_implicit`,
`fields_local`, `nonlinear_terms`, `le_derivatives`, `collisions`, `gs2_main`, every program's object list
gains `gs2_export_probe.o` after `make depend` — this is why `make depend` is mandatory before rebuilding
even `bin/gs2`. Without it, compilation of those modules fails with a missing `gs2_export_probe.mod`.

### 9.3 Fallback: isolated object tree (only if §4.4 step 3 fails)

Wrap every probe line in `#ifdef GS2_EXPORT_PROBES … #endif`, keep the `public` blocks unconditional, and
build the exporter into a separate tree so that `bin/gs2` objects are untouched:

```
MACRO_DEFS=-DGS2_EXPORT_PROBES make GK_SYSTEM=… OBJDIR=$PWD/objs_export MODDIR=$PWD/include_export gs2_export
```

(`OBJDIR`, `MODDIR` are plain `:=` variables, `Makefile:190-191`, so command-line values override them;
`Makefile.depend` refers to `$(OBJDIR)` symbolically. `MACRO_DEFS` must come from the **environment**:
`Makefile:172` is `MACRO_DEFS ?=` and all later `MACRO_DEFS += …` lines append to the environment value,
whereas a command-line `MACRO_DEFS…` would override and silently drop `-DMPI`, `-DNETCDF`, ….) The reference trajectory then comes from a binary
whose physics objects differ from `bin/gs2`'s only by the probe calls, and §4.4 step 3 must be repeated
between `bin/gs2` and `bin/gs2_export` with `nstep_trace = N`, `trace_stages = 'step'` vs. the `.out.nc`.

---------------------------------------------------------------------------------------------------

## 10. GPU mapping notes (loader → JAX)

* **Host-side, once.** The loader runs on the host (netCDF4 + NumPy), converts `ri` → `complex128`,
  `int8` → `bool`, and `jax.device_put`s each array once. Nothing in the export is float32; enable
  `jax_enable_x64` before loading.
* **Axis order is already batch-friendly.** In NumPy every g-like array is `(S, E, L, Y, X, σ, θ)`: the five
  leading axes are the independent batch axes of the per-`iglo` kernels (source, invert_rhs, g_adjust,
  coefficients) → flatten them to one batch axis (or `vmap` over them); θ is last (contiguous), which suits
  the θ recurrences of invert_rhs (`lax.scan` over a moved θ axis) and the per-θ vector work. Fields
  `(Y, X, θ)` broadcast against the `(Y, X, σ?, θ)` tail.
* **Summation order.** Velocity-space reductions must follow `iglo` order on one rank: use
  `/gs2_layouts/iglo_of` to build, once, the permutation from canonical `(S,E,L,Y,X)` to glo order, and the
  per-(it,ik) ordered list of contributions; the reduction itself is then a `lax.scan`/`fori_loop` over that
  static order (owner specs 04, 08). LE arrays follow `ile_of` (06).
* **Irregular structures.** Supercell classes: `i_class` is small and static → keep one dense array per class
  (`aminv_c<c>` `(M_c, R_c, R_c)`), loop over classes in Python (static), batch members with `vmap` or a
  batched matmul whose inner sums are ordered as spec 04 §5.3 requires. `cls_it/cls_ik` are padded with −1 to
  `(i_class, Nmax, Mmax)` for gather tables; mask with `cls_it >= 1`. Linked chains: `itleft/itright` with
  `-1` sentinels, `itleftmost/itrightmost`, `l_links/r_links` → static gather indices (03 §5.3).
* **Masks** come straight from the export: `forbid`, `is_ttp`, `is_bounce_point`, derived lower/upper bounce
  masks, wfb masks, `kwork_filter` (after init), `ieqzip`, `connections_neighbor`, `save_h`.
* **dt changes.** Either rebuild from the exported dt-independent inputs (RECOMPUTE paths of 01–07) or look
  up `/dt_ladder/k<k>` by the exact `code_dt` bit pattern (6.3).
* **Verification on device.** The trace groups are small enough to load lazily per probe; compare each JAX
  stage against the matching `p<seq>_<tag>` group (same istep, attempt, pass).

---------------------------------------------------------------------------------------------------

## 11. Pitfalls (things that would break round-off agreement or silently corrupt the export)

11.1 **Reference runs must be nproc = 1 and single-threaded.** All integrals and the field solve use
`MPI_Allreduce` (mp.fpp:3002-3018); the reduction order across ranks is implementation-defined. With one
rank the reduction is an identity copy. Build has no OpenMP; with an OpenMP build require
`OMP_NUM_THREADS=1` (OMP loops in get_source_term_on_local_domain etc. are per-iglo, but reductions in
integrals would reorder).

11.2 **Export the t = 0 state after `initialize_diagnostics`.** `run_diagnostics(0)` changes `gnew` (and `g`)
by `(g + adj) − adj` round trips (spec 08 §0 item 3). `state_t0_prediag` vs `state_t0` makes this visible;
JAX must start from `state_t0` *and* replicate the round trips on its own write steps. If
`run_diagnostics=.false.`, the Fortran trace is a different (diagnostics-free) trajectory than `bin/gs2`'s.

11.3 **`kwork_filter` changes during init.** It is `.false.` at allocation, toggled by fields_local during
the response build (fields_local.f90:2006-2076), and `(1,1)` is set `.true.` for box grids only after the
response build (fields.f90:376). Export it after `initialize_equations`, never earlier.

11.4 **Lazily allocated arrays.** `awgt` (first `calculate_flux_surface_average`, normally during the response
build), `inv_phi/bpar_denominator_g` (first `get_fields_direct_from_dfn`, normally in `set_init_fields`),
`jump`, `ikx_indexed` (first flow-shear step). The exporter must not call physics routines to force their
creation (that would change module state, e.g. `awgt` is frozen at first use, spec 04 P8); write if allocated,
else list in `absent`.

11.5 **Do not perturb state while exporting.** No call into dist_fn/fields/collisions routines from the
exporter or the hooks (several keep persistent work arrays — `antot*`, `fieldeq*`, `g_work`, `g_adj` —
and lazily initialised data). Only reads, pure/elemental functions (`il_is_*`, `is_*_bounce_point`,
`get_leftmost_it`, `get_rightmost_it`, `spfunc::j0/j1`, `trapezoidal_integration`) and the public dump
routine of §7.5 (which only reads `fieldmat`).

11.6 **`ulim_alloc` vs `ulim_proc`.** Empty ranks allocate one slot (`ulim_alloc = llim_proc` when
`ulim_proc < llim_proc`, gs2_layouts.fpp:864); never send/receive or canonicalise that slot. The same holds for
le_lo, jf_lo, f_lo, lz_lo, e_lo.

11.7 **Assumed-shape bounds in hooks.** Arrays passed to `probe_g`/`probe_le` arrive with lower bounds 1
(§5.1); forgetting the `-ntgrid` / `llim_proc` offsets shifts θ by `ntgrid+1` points.

11.8 **8-D complex variables.** neasy-f stops at rank 7 (`RANKS = range(0, 8)`); G7 complex with `ri` is
rank 8. Use the rank-1 buffer + explicit `count` write of §2.4 (verified, Appendix A.1). Never call the
rank-1 writer without `count`.

11.9 **`gnostics%nsave` when diagnostics are skipped.** `evolve_equations` reads `gnostics%nsave`
(gs2_main.fpp:764) even if diagnostics were never initialised; the component has no default
initialiser (diagnostics_config.f90:150). Set `gnostics%nsave = -1` in the skip branch.

11.10 **`finalize_diagnostics` aborts if diagnostics were not initialised** (gs2_main.fpp:996-998); call it
only when `state%init%diagnostics_initialized`.

11.11 **Logical and integer conventions.** Logical → int8 0/1 (never rely on Fortran's internal logical bit
pattern). Index arrays keep 1-based Fortran values and −1 sentinels (`itleft`, `itright`, `ixi_to_il`,
`cls_*`); `iglo_of`/`ile_of` are 0-based by construction. The loader must not "helpfully" shift them.

11.12 **64-bit sizes.** `2*nth*2*ntheta0*naky*nlambda*negrid*nspec` exceeds 2³¹ for large grids
(e.g. nth=129, ntheta0=127, naky=64, nlambda=40, negrid=16, nspec=3 → 8.0e9). Use `integer(int64)` for buffer
sizes and offsets, and the per-species slab writes of §2.4 to bound memory.

11.13 **No `TARGET`/pointer views.** Adding `target` to module arrays changes alias analysis and can change
the generated code; giving the exporter pointer views is forbidden (§4.1 E). Also do not add
`contiguous`/`volatile`, and do not move declarations.

11.14 **FFTW plans and wisdom.** With `fft_use_wisdom = .true.` (default) GS2 imports
`<run>.fftw_wisdom` (or `$GK_FFTW3_WISDOM`); FFTW applies wisdom created with MEASURE to ESTIMATE planning,
so the presence of a wisdom file from an earlier run can change the FFT algorithm and its round-off.
For reproducible Fortran references delete stale wisdom files (or set `fft_use_wisdom = .false.`) and keep
`fft_measure_plan = .false.`; the exporter records both settings. (JAX's FFT matches only at FFT round-off
level anyway, spec 05 §5.)

11.15 **Restart-file round trips are not bit-transparent by default.** On restart, `set_init_fields`
recomputes the fields from `g` when `force_maxwell_reinit = .true.` (default, gs2_init.fpp:318), and
`scale ≠ 1` rescales (gs2_save.fpp:863-875). Use `force_maxwell_reinit = .false.`, `scale = 1.0` for
bitwise round trips (§6.4).

11.16 **Layout check on restart.** `gs2_restore_many` aborts if the file's `layout` differs (767-771); a
JAX-written restart must use the target layout's glo order.

11.17 **Response-matrix rebuild calls `timeadv(…, 0)`.** Inner probes must be filtered by
`probe_in_advance` or the trace fills with `istep = 0` junk (and costs a factor ~`nidx*N` in time).

11.18 **Reset attempts reuse the NL term.** After an immediate reset, attempt 2 uses the `gexp_*` already
shuffled in attempt 1 (`istep_last` guard, nonlinear_terms.fpp:495-497); a JAX implementation that
recomputes or re-shuffles on retry diverges. The trace's attempt structure exposes this.

11.19 **Dummy/module aliasing in hooks.** Inside `add_explicit` and `vspace_derivatives` the data live in
dummy arguments that alias module arrays and are modified through the dummy; read them only through the
array passed to the probe (§4.3), otherwise the standard does not guarantee the values.

11.20 **Strings.** `layout` is `character(len=5)` and is compared verbatim on restart; write it untrimmed in
the restart writer, trimmed in attributes. neasy-f writes character scalars with an extra `<name>_dim`
dimension (`neasyf.in.f90`, character rank-0 branch); in the export use attributes instead.

11.21 **Fill values / masked arrays.** netCDF4-python masks values equal to `_FillValue`; call
`set_auto_mask(False)` (§6.6). Unwritten slabs of a partially written variable would read as fill values —
the slab loop must cover every species.

11.22 **`-0.0`.** Copies preserve signed zeros; reductions (`sum_reduce` gathers) would not (`-0.0 + 0.0 =
+0.0`). Signed zeros can matter for `atan2`-like diagnostics and for bitwise comparisons; use the
send/receive gathers of §2.5.

11.23 **Do not trust `dump_response_to_file_imp` for the runtime matrix layout.** It removes duplicate
boundary rows/columns and reorders (fields_implicit.f90:1251-1399); the runtime `getfield` uses the raw
`aminv` rows exported in §7.4.

11.24 **The two `timeadv` calls and two collision applications per step** (spec 06 §8.9) produce the same
tags twice; only `probe_pass` and `seq` distinguish them.

---------------------------------------------------------------------------------------------------

## 12. Verification plan (to run once the exporter exists)

Test inputs (single rank, no wisdom file): (a) linear electrostatic, adiabatic electrons, `grid_option =
'single'` or `'range'`, s-alpha; (b) linear electromagnetic (`beta > 0`, `fbpar = 1`), kinetic electrons,
collisions (`collision_model = 'default'`), Miller; (c) nonlinear `box`, `nx = ny = 16`, `nstep = 10`,
hyperviscosity on, flow shear `g_exb ≠ 0`; (d) case (c) with `field_option = 'local'`; (e) a case with
`delt` large enough to trigger a CFL reset within the traced steps.

* **T0 build identity** (§4.4): `otool -tv` identical for the `public`-only edits; `.out.nc` and restart
  files bitwise identical between `bin/gs2.orig` and the rebuilt `bin/gs2` for (a)–(e).
* **T1 export self-consistency**: loader recomputes `vperp2 = (bmag*al)*energy` (0 where forbid),
  `tunits/wunits`, `has_*`, `field_line_average_weight` (Klein trapezoid) and requires bitwise equality.
* **T2 canonical map**: `iglo_of` vs the restart file's `iglo_kykxsel` (`write_restart_t0 = .true.`);
  `gnew` from `/state_t0` vs `gr/gi` of the same restart file re-indexed with `iglo_kykxsel` — bitwise.
* **T3 restart round trip**: JAX loader writes `/state_t0` back as a GS2 restart (§6.4); `bin/gs2` with
  `ginit_option='many'`, `scale=1`, `delt_option='check_restart'`, `force_maxwell_reinit=.false.` must
  reproduce the trace of step 1 bitwise.
* **T4 trace invariants** (§6.2): `step_end(k)` = `advance_begin(k+1)`; for (e) attempt 2 exists and has no
  `nl_*` records.
* **T5 ladder**: `/dt_ladder/k1` coefficients equal those of a GS2 run that resets once (trace of (e),
  `after_reset_time_step` followed by a fresh export) — bitwise.

---------------------------------------------------------------------------------------------------

## 13. Open questions

* **Q1** Source edits: this phase is read-only; the next phase must touch 11 existing files (§4.2, §5.2) and add
  5. Is editing upstream GS2 sources on branch `gpu-port` acceptable, or should the edits live in a patch
  applied to a copy of `src/` (then `VPATH`/`OBJDIR` overrides as in §9.3)?
* **Q2** Is bitwise identity between `bin/gs2.orig` and the instrumented `bin/gs2` a hard requirement (§4.4
  step 3)? If the probe calls perturb code generation, do we accept the §9.3 isolated build (reference =
  `bin/gs2_export`, not `bin/gs2`)?
* **Q3** Default `run_diagnostics`: `.true.` reproduces `bin/gs2` (diagnostics round trips included) but forces
  the JAX side to replicate diagnostics-induced perturbations from step 0 (spec 08). Confirm.
* **Q4** Should the JAX side consume `used_inputs` (via f90nml) at all, or only exported effective values?
  (This spec exports both; the effective values are authoritative.)
* **Q5** `field_option='local'`: is the file round trip of §7.5 sufficient, or does spec 04 §7.3 need the
  pre-inversion local matrices too (not reachable without making fields_local type-bound procedures public)?
* **Q6** Multi-rank references (spec 04 Q8): out of scope? The exporter supports `allow_multi_rank` for static
  data, but the values are rank-dependent at round-off.
* **Q7** Trace size: full G7 dumps at ~30 probes × 2 passes per step are large for production grids; is a
  ky/kx subset option needed, or are traces restricted to small verification grids?
* **Q8** Naming: groups-per-module (this spec) vs. the flat names suggested in specs 03/04/06/07 (mapping in
  §6.1.3). Confirm groups.
* **Q9** Contract for dt changes (spec 01 Q2): dt ladder lookup, JAX RECOMPUTE, or both?
* **Q10** P-R1 (`am` before inversion) is included as optional; spec 04 Q10 asked whether it is wanted.
* **Q11** If JAX must restart from a mid-run trace snapshot, the SAVE locals of §3.13 must be moved to module
  scope (behaviour-preserving). Needed?
* **Q12** Should `gs2_export` also produce `<run>.out.nc` "static part only" without running diagnostics (spec 08
  §0 item 7)? Currently `.out.nc` exists only with `run_diagnostics = .true.`.
* **Q13** The exact make command line used to build the existing `bin/gs2` is not recorded in the repository
  (no `Makefile.local`, `GK_SYSTEM` not in the environment, include paths differ from
  `Makefiles/Makefile.macosx-homebrew`). It must be recovered (or a `Makefile.local` written) so that
  `bin/gs2_export` is built with identical flags.

---------------------------------------------------------------------------------------------------

## Appendix A — Checks performed for this spec (scratch programs outside the repository)

A.1 **8-D write through neasy-f.** A program linked against `objs/neasyf.o` and netCDF-Fortran 4.6.2 wrote a
`real64` rank-1 buffer of 2880 values into a variable with dims `(ri, theta, sign, kx, ky, lambda, egrid,
species) = (2,3,2,4,5,2,3,2)` inside a group (`nf90_def_grp`), with dims defined only in the root group, using
`neasyf_write(gid, "g", flat, dim_names=[…8 names…], count=d)`. netCDF4-python read it as shape
`(2,3,2,5,4,2,3,2)` with dims reversed; Fortran element `(2,3,1,4,5,2,3,2)` (flat index 2874) appeared at
NumPy `[1,2,1,4,3,0,2,1]`; `a[...,0] + 1j*a[...,1]` gave a `complex128` array of shape `(2,3,2,5,4,2,3)`. An
`int8` logical mask round-tripped as int8.

A.2 **PRIVATE vs PUBLIC code generation.** A module with two module variables compiled with
`gfortran -O3 -fPIC -mcpu=apple-m1 -fdefault-real-8 -S` once with the variables private and once public: the
assembly differs only in the `.private_extern` directives (symbol visibility); the instruction stream is
identical. `nm -m` shows `private external` vs `external`.

A.3 **Public variable of a private type.** `type, private :: ct` with a `public` allocatable array of type `ct`:
components readable from a program (`conn(2)%iglo_left`, `conn%iglo_left`) under `gfortran -std=f2008`.

A.4 **No bare `use`.** `grep` over `src/` and `tests/` found no `use dist_fn|collisions|hyper|nonlinear_terms|
fields_implicit|theta_grid|le_grids|run_parameters` without `only:`, so the new public names cannot clash.

## Appendix B — Change list for the implementation phase (summary)

| file | change | lines touched |
|---|---|---|
| `src/gs2_export_probe.f90` | new (§5.1) | ~90 |
| `src/gs2_export_writer.fpp` | new (§7.2, §2.4, §2.5) | ~900 |
| `src/gs2_export_static.fpp` | new (§7.3–7.8) | ~1200 |
| `src/gs2_export_trace.fpp` | new (§7.9) | ~400 |
| `src/programs/gs2_export.fpp` | new (§1.4, §8.1, §8.3) | ~350 |
| `src/dist_fn.fpp` | `public` block after 96; probes after 4062, 4068, 4077, 4084, 4091, 6182, 6184; `use gs2_export_probe` in `timeadv`, `invert_rhs` | +20 |
| `src/fields.f90` | probes in `advance` after 593 and 604 | +12 |
| `src/fields_implicit.f90` | `public :: aminv, nfield, nidx` after 13; probes after 335, 341, 352, 356, 362, 365, 369 (+ pass markers before 341, 365); optional P-R1 before 638 | +16 |
| `src/fields_local.f90` | probes after 4280, 4297, 4313, 4320, 4331, 4337, 4340 (+ pass markers) | +14 |
| `src/nonlinear_terms.fpp` | `public` block after 24; probes after 520, 523, 569, 582 | +8 |
| `src/le_derivatives.fpp` | probes after 71, 83, 90, 94, 133, 140 | +8 |
| `src/collisions.fpp` | `public` block after 23; probes after 3051, 3052, 3085, 3087, 3088 | +14 |
| `src/gs2_main.fpp` | probes after 744, 771, 785, 811 | +6 |
| `src/hyper.f90`, `src/theta_grid.f90`, `src/le_grids.f90`, `src/run_parameters.f90` | `public` blocks (§4.2) | +2..4 each |
| `Makefile.target_gs2` | §9.2 (3 lines) | +4 |
| `Makefile.depend` | regenerate with `make depend` (git-ignored) | — |
