# Spec 01 — The time loop and time-step control

Status: contract for implementers of the JAX (float64/complex128) time stepper.
Scope: `gs2_main::evolve_equations`, the `fields::advance` dispatch, `advance_implicit` /
`advance_local`, the step-level structure of `dist_fn::timeadv` and `collisions_advance`,
`exb_shear`, `remove_zonal_flows`, the antenna, `gs2_time`, `gs2_reinit` (CFL checks,
`reset_time_step`, what is rebuilt on a dt change), stopping conditions, the meaning of `istep`
(including `istep = 0` during response-matrix construction), and the t = 0 state that must be
exported to start the JAX run identically.

Not in scope (owned by other specs, only referenced here): the body of `invert_rhs` /
`get_source_term` (source + parallel solve), the response-matrix algebra (`init_inverse_matrix`,
`fields_local` blocks), `getfieldeq` / `getan`, the collision operator (`solfp1`), the nonlinear
term (`add_nl`, transforms), `hyper_diff` internals, `get_fields_direct_from_dfn`, diagnostics.

All file:line citations are against the tree at commit `4d8c94bcf` (branch `gpu-port`).

---------------------------------------------------------------------------------------------------

## 0. Build facts and conventions that matter for round-off

* `bin/gs2 --build-config` reports: GCC 15.2.0, `-O3 -mcpu=apple-m1`, `-fdefault-real-8
  -fdefault-double-8`, `MPI`, **no OpenMP**, **no LOWFLOW**, `NEW_DIAG` on, FFTW3, netCDF,
  `F200X_INTRINSICS`, **no `-DRANDOM=_RANMT_`** (so `ranf` = Fortran intrinsic `random_number`,
  externals/utils/ran.fpp:21-45).
* Because of `-fdefault-real-8`, every literal (`0.5`, `1.e8`, `1.e-5`, `3.`) is a double, and
  `epsilon(0.)` = `epsilon(0.0)` = 2.220446049250313e-16 (used as a threshold in
  run_parameters.f90:351-353, fields_local.f90:4280, dist_fn.fpp:3983, antenna.f90:213,397).
* `-O3` on aarch64 with GCC's default `-ffp-contract=fast` will contract `a*b+c` into FMA
  instructions wherever the optimiser chooses. Bitwise identity with the Fortran binary is
  therefore not achievable; the target is ~1e-12 relative. See Pitfalls §6.1.
* `real :: x = <literal>` local declarations (`fac = 1.0` in gs2_reinit.f90:124,236) carry an
  implicit SAVE; they are never modified so this is harmless.
* `-O3` without `-ffast-math`: gfortran honours source evaluation order for `a*b*c*d`
  (left-to-right: `((a*b)*c)*d`) and does not reassociate `sum` intrinsics (sequential loop).

### Index conventions used in this document

* g-like arrays: Fortran `g(ig=-ntgrid:ntgrid, isgn=1:2, iglo)` with `iglo` expanded to
  `(it, ik, il, ie, is)` via `it_idx/ik_idx/il_idx/ie_idx/is_idx` →
  NumPy shape `(nspec, negrid, nlambda, naky, ntheta0, 2, 2*ntgrid+1)`.
* g-like arrays without sigma (e.g. `aj0(ig, iglo)`) → `(nspec, negrid, nlambda, naky, ntheta0, 2*ntgrid+1)`.
* Fields `phi(ig, it, ik)` → NumPy `(naky, ntheta0, 2*ntgrid+1)`.
* `(it, ik)` 2-D arrays (`kwork_filter`, `ieqzip`, `theta0`) Fortran `(ntheta0, naky)` → NumPy `(naky, ntheta0)`.
* 1-D arrays keep their index: `aky(ik)`, `akx(it)`, `ikx(it)`, `kx_shift(ik)`, `tunits(ik)`,
  `theta(ig)`, `jacob(ig)`.
* `wstar(ik, ie, is)` Fortran → NumPy `(nspec, negrid, naky)`.

In pseudo-code below, Fortran 1-based / `-ntgrid`-based indices are kept; "`gnew[.., it, ik, ..]`"
means the canonical array.

---------------------------------------------------------------------------------------------------

## 1. Call graph / control flow

### 1.1 Program entry and initialisation

```
src/programs/gs2.fpp:52-57      parse_command_line; init_mp; run_gs2(state); finish_mp
gs2_main.fpp:1614  run_gs2
  1640 initialize_gs2           (mp, file utils, init_gs2_init; level := basic)
  1641 initialize_equations     gs2_main.fpp:467-535
         486  init_parameter_scan
         492  init(state%init, init_level_list%full)      -> walks levels 2..35 upwards (table 1.2)
         496  init_tstart(tstart)                           gs2_time.f90:129-136  user_time = code_time = tstart
         523  state%nstep = nstep ; do_eigsolve
         531  init_gs2_reinit -> init_reinit               gs2_reinit.f90:257-302
                 smart default dt0 = code_dt_max (282-286); read reinit_knobs;
                 save_dt_min(delt_minimum) (300); save_dt_max(dt0) (301)
  1642 initialize_diagnostics   gs2_main.fpp:542-614
         600  run_diagnostics(0, state%exit)    (new diags)  or  607 loop_diagnostics(0, ...)
  1647 if do_eigsolve -> run_eigensolver (1648)  else  evolve_equations(state, state%nstep)  (1650)
```

### 1.2 Init-level ladder (gs2_init_level_list.inc:1-35, handlers gs2_init_subroutines.inc)

| lvl | name | up handler (file:line) | relevant to time loop |
|---|---|---|---|
| 7 | init_g | init_init_g (subroutines.inc:75-85) | reads tstart, scale, new_field_init, ginit_option |
| 13 | kt_grids | init_kt_grids | allocates `kwork_filter=.false.` (kt_grids.f90:188-189) |
| 14 | run_parameters | init_run_parameters (run_parameters.f90:247-273) | has_*, nstep, delt, ieqzip |
| 16 | gs2_time | init_gs2_time (gs2_time.f90:45-69) | code_dt, prev1/2, code_dt_max, tunits/wunits |
| 19 | hyper | init_hyper (hyper.f90:360-449) | **guarded by `initialized`; never rebuilt on dt change** |
| 20 | antenna | init_antenna (antenna.f90:196-224) | no_driver, beta_s |
| 24 | nonlinear_terms | init_nonlinear_terms | |
| 26 | dist_fn_arrays | init_dist_fn_arrays (dist_fn.fpp:885-925) | allocates g, gnew, g_work, g_h, gexp_1..3, kx_shift/theta0_shift (3943-3991) |
| 27 | dist_fn_level_1 | vperp2, vpa/vpac, init_bc (928-936) | `finish` resets `exb_first=.true.` (1080) |
| 28 | dist_fn_level_2 | init_bessel, init_fieldeq (939-952) | |
| 29 | override_timestep | run_parameters::set_overrides (immediate_reset override) | **reset target level** |
| 30 | collisions | init_collisions / finish_collisions (collisions.fpp:436-463, 4165-) | dt-dependent |
| 31 | dist_fn_level_3 | init_dist_fn_level_3 (dist_fn.fpp:960-1013) | vpar, wdrift, wstar, invert_rhs coeffs, g_h, source coeffs — dt-dependent |
| 32 | fields_level_2 | init_fields (fields.f90:353-377) | response matrix, nl_finish_init, kwork_filter(1,1) |
| 33 | override_initial_values | no-op | |
| 34 | set_initial_values | set_initial_field_and_dist_fn_values (gs2_init.fpp:292-375) | ginit + set_init_fields |
| 35 | full | no-op | |

`init_pass_int` (gs2_init.fpp:237-269) walks one level at a time; going down from level L calls
level L's `down` handler and lands on L-1.

### 1.3 Initialisation of the equations (what produces the t = 0 state)

```
init_dist_fn_level_3 (dist_fn.fpp:960)     init_vpar, init_wdrift, init_wstar, init_collisions,
                                           init_invert_rhs (-> init_homogeneous_g), init_hyper,
                                           init_source_term, [init_g2gf]
init_fields (fields.f90:353)
  fields_pre_init (241-279)               init_fields_parameters, theta_grid, run_parameters,
                                           init_dist_fn, gf_local->local fallback (261-269),
                                           init_antenna, allocate_arrays (fields; zeros phi..bparnew)
  fields_init_response (323-350)          implicit: init_fields_implicit -> init_response_matrix
                                                     (fields_implicit.f90:459-647)
                                           local:    init_fields_local -> init_fields_matrixlocal
  nl_finish_init (373)                     nonlinear_terms.fpp:1262-1265  initializing := .false.
  if (gridopt_switch == gridopt_box) kwork_filter(1,1) = .true.     (fields.f90:376)
set_initial_field_and_dist_fn_values (gs2_init.fpp:292)   [override = .false. at first init]
  ginit(restarted)                         init_g.f90:1075-1243
  if (.not.restarted .and. initial_condition_is_nonadiabatic_dfn) return   (311)
  if (.not.restarted .or. force_maxwell_reinit) call set_init_fields       (318)
     -> init_allfields_implicit (fields_implicit.f90:102-136), also used by 'local'
        (fields_local.f90:4234-4241):
        new_field_init (default .true.):  get_init_field(phinew,aparnew,bparnew)  (dist_fn.fpp:7792-7798
                                           = get_fields_direct_from_dfn(gnew, ...))
                                          phi=phinew; apar=aparnew; bpar=bparnew; g=gnew
        else:                             getfield(phinew,...) (implicit solve); phi=phinew ...
init_tstart(tstart)                        gs2_main.fpp:496
init_gs2_reinit                            gs2_main.fpp:531
initialize_diagnostics -> run_diagnostics(0)   (may do in-place g_adjust round trips on g/gnew, §6.14)
```

### 1.4 Response-matrix construction: the meaning of `istep = 0`

`init_response_matrix` (fields_implicit.f90:459-647), per class `i = i_class..1` (541):
zero `am, g, phi, phinew, apar, aparnew, bpar, bparnew` (553-557) — note **`gnew` is not zeroed
here** (it is zero at first init from `allocate_arrays`, and is zeroed by `reset_time_step`
gs2_reinit.f90:163 before a rebuild). For each cell `n`, theta point `ig`, field: place unit
impulses in `phinew` (all members `m` of the class at once, 584-588; plus the connected
end point `phinew(ntgrid, it(m,n-1), ik(m,n-1)) = 1` when `n>1 .and. ig==-ntgrid`, 575-582),
then `init_response_row` (650-719) calls

```
timeadv(phi, apar, bpar, phinew, aparnew, bparnew, 0)      fields_implicit.f90:669
getfieldeq(phinew, aparnew, bparnew, fieldeq, fieldeqa, fieldeqp)   670
```

`fields_local` does the same via `fm_init_next_field_points` → `timeadv(..., 0)`
(fields_local.f90:2184), manipulating `kwork_filter` itself (2006, 2032, 2058) and finally
resetting `kwork_filter = .false.` (2076).

What `istep = 0` changes inside one `timeadv` call (dist_fn.fpp:4014-4116):

| place | effect of istep == 0 |
|---|---|
| `add_explicit_terms` (nonlinear_terms.fpp:424-435) | nonlinear: `add_explicit` skipped (`istep /= 0` guard, 425); linear: `save_dt_cfl(1e8)` executed as for every call |
| `get_source_term_on_local_domain` (dist_fn.fpp:5382-5389) | `include_explicit = istep > 0` → no Adams–Bashforth explicit terms |
| `hyper_diff` (4073) | only if `hyper_in_initialisation` (default `.true.`) |
| `vspace_derivatives` (4083-4087) | collisions applied (if `.not. split_collisions`), `mode` absent → no heating |
| split nonlinear (4049, 4097) | skipped (`istep > 0` guards) |
| `exb_shear` | not called (only called from `advance_*`) |
| `timeadv` SAVE `istep_last` (4031, 4114) | set to 0 |

The time loop itself uses `istep = 1..nstep` (gs2_main.fpp:715-717). `istep = 0` additionally
appears as the argument of `run_diagnostics(0)` / `loop_diagnostics(0)` at diagnostics init.

### 1.5 The main loop (`evolve_equations`, gs2_main.fpp:649-894)

```
715 istep_loop_start = state%istep_end + 1 ; istep_loop_max = state%istep_end + nstep_run
717 do istep = istep_loop_start, istep_loop_max
718    if (istep > nstep) exit                          ! only for repeated evolve calls
730    call advance(istep)                              ! fields.f90:586-605
734    if (state%dont_change_timestep) reset = .false.  ! library flag, default .false.
739    if (reset) then                                  ! only possible with immediate_reset (CFL too large)
740       prepare_initial_values_overrides(state)       ! init level -> 32 (no-op teardown), alloc override store
742       set_initval_overrides_to_current_vals(...)    ! save gnew, vnmult, gexp_*, [phinew..] or write restart files
743       state%init%initval_ov%override = .true.
744       reset_time_step(state%init, istep, state%exit, job_id)
747       reset = .false.
751       if (.not. state%exit) call advance(istep)     ! redo the step with the new dt
       end if
754    if (state%exit) exit
759-769  if (nsave>0 .and. mod(istep,nsave)==0) gs2_save_for_restart(gnew, user_time, ...)  ! BEFORE update_time
771    call update_time                                 ! gs2_time.f90:161-174
775-783  if (.not. skip_diagnostics) run_diagnostics(istep, state%exit) [or loop_diagnostics]
785    if (state%exit) state%converged = .true.
792    if (.not. state%exit) then
800       call check_time_step(reset, state%exit)        ! gs2_reinit.f90:228-254
806       if (dont_change_timestep) reset = .false.
807       if (reset) then
808          prepare_initial_values_overrides; 809 set_initval_overrides_to_current_vals
810          override = .true. ; 811 reset_time_step(...) ; 815 reset = .false.
          end if                                        ! NOTE: no re-advance here; next step uses new dt
818       if (mod(istep,ncheck_stop)==0 .and. .not.exit) checkstop(exit, exit_reason)   ! <run_name>.stop
822       if (.not.exit) checktime(avail_cpu_time, exit, margin_cpu_time)               ! wall clock
       end if
829    state%istep_end = istep
833    if (.not.exit .and. user_time + user_dt > max_sim_time) -> exit_reason = EXIT_MAX_SIM_TIME; exit = .true.
839-858 progress messages (cosmetic)
860    if (state%exit) exit
    end do
870-881 final forced diagnostics at istep - istep_offset (offset 1 if loop completed, 0 if exited)
```

`run_gs2` sets `exit_reason = EXIT_NSTEP` if nothing else was requested (gs2_main.fpp:1651).

### 1.6 `fields::advance` dispatch (fields.f90:586-605) and `field_option` mapping

`fieldopts` (fields.f90:42-48): `'default'`→`fieldopt_implicit`, `'implicit'`→implicit,
`'test'`→test, `'local'`→local, `'gf_local'`→gf_local, `'implicit_local'`→local.
`gf_local` silently falls back to `local` if `nproc < ntheta0*naky` (fields.f90:261-269) or if
`gf_lo_integrate`/`gf_local_fields` are not both set (set_init_fields, 425-435).

```
select case (fieldopt_switch)
  implicit : advance_implicit(istep, remove_zonal_flows_switch)   fields_implicit.f90:315-371
  test     : advance_test(istep)
  local    : advance_local(istep, remove_zonal_flows_switch)      fields_local.f90:4244-4342
  gf_local : advance_gf_local(istep, remove_zonal_flows_switch)
```

### 1.7 `advance_implicit` (fields_implicit.f90:315-371)

```
333  if (.not. no_driver) call antenna_amplitudes(apar_ext)
335  if (allocated(kx_shift) .or. allocated(theta0_shift)) call exb_shear(gnew, phinew, aparnew, bparnew, istep)
337  g = gnew
338  phi = phinew ; apar = aparnew ; bpar = bparnew          ! all three, unconditionally
341  call timeadv(phi, apar, bpar, phinew, aparnew, bparnew, istep)            ! predictor, mode absent
343  if (reset) return                                         ! immediate CFL reset: gnew still == g
345  if (.not. no_driver) aparnew = aparnew + apar_ext
352  call getfield(phinew, aparnew, bparnew)   ! fl = fieldeq(gnew_pred, phinew=phi_old); phinew := -aminv*fl (delta)
354  phinew  = phinew  + phi
355  aparnew = aparnew + apar
356  bparnew = bparnew + bpar
362  if (remove_zonal_flows_switch) call remove_zonal_flows    ! modifies phinew only
365  call timeadv(phi, apar, bpar, phinew, aparnew, bparnew, istep, diagnostics=1)  ! corrector
369  call collisions_advance(phi, bpar, phinew, aparnew, bparnew, diagnostics=1)    ! no-op unless split_collisions
```

### 1.8 `advance_local` (fields_local.f90:4244-4342) — differences from implicit

```
4270 do_update = do_smart_update ; 4272 if (g_lo%x_local .and. g_lo%y_local) do_update = .false.
4275 antenna_amplitudes (same)
4280 if (abs(g_exb*g_exbfac) > epsilon(0.)) call exb_shear(gnew, phinew, aparnew, bparnew, istep, field_local_allreduce_sub)
4283 g = gnew
4284-4290 if (do_update) fieldmat%update_fields_newstep  else  {if has_phi: phi=phinew; if has_apar: apar=aparnew; if has_bpar: bpar=bparnew}
4296 timeadv(..., istep)          ; 4298 if (reset) return
4302 aparnew = aparnew + apar_ext (if driver)
4313 getfield_local(phinew, aparnew, bparnew, do_gather=.true., do_update)
4316-4320 if (.not.do_update) { if has_phi: phinew = phinew + phi ; has_apar ; has_bpar }
4331 remove_zonal_flows (same routine, from fields_implicit)
4336 timeadv(..., istep, diagnostics=1)
4340 collisions_advance(...)
```

On a single process (the JAX reference configuration) `do_update` is always `.false.`.
The exb condition `abs(g_exb*g_exbfac) > epsilon(0.)` is the same condition under which
`kx_shift`/`theta0_shift` were allocated (dist_fn.fpp:3983), so both paths call `exb_shear`
in the same circumstances. Copying only `has_*` fields is equivalent because the arrays of
absent fields are identically zero (they are never written), except for `aparnew + apar_ext`
when the antenna is on with `has_apar = .false.` (antenna is GATEd).

### 1.9 `timeadv` structure (dist_fn.fpp:4014-4116), as seen from the loop

```
4036 is_first_half_step  = istep /= istep_last       ! SAVE istep_last, initial -1
4037 is_second_half_step = istep == istep_last
4039 nl_dt = strang_split ? code_dt*0.5 : code_dt
4046 if (split_nonlinear .and. nonlin) then  [GATED]
     else
4060    call add_explicit_terms(gexp_1, gexp_2, gexp_3, phi, apar, bpar, istep, bkdiff(1))
     end if
4066 if (reset .and. immediate_reset) return          ! before invert_rhs, before istep_last update
4068 call invert_rhs(phi, apar, bpar, phinew, aparnew, bparnew, istep)   ! source + parallel solve (spec 02/03)
4073 if ((istep == 0 .and. hyper_in_initialisation) .or. istep > 0) call hyper_diff(gnew, phinew)
4083 if (.not. split_collisions) call vspace_derivatives(gnew, g, g_work, phi, bpar, phinew, bparnew, mode)
4091 if (def_parity) call enforce_parity(parity_redist)
4094 [split nonlinear Strang second half — GATED]
4114 istep_last = istep
```

Consequences for the loop:
* **Collisions are applied in both `timeadv` calls of a step** (predictor and corrector) when
  `split_collisions = .false.` (default), and in every response-matrix `timeadv(...,0)` call.
* `add_explicit_terms` (nonlinear_terms.fpp:407-439) → `add_explicit` (442-612): the explicit
  history is shifted and the nonlinear term computed **only when `istep /= istep_last`**
  (module variable, nonlinear_terms.fpp:32, initial 0), i.e. once per step, in the first
  `timeadv` call; `istep_last` is set at 497 *before* the CFL check, so after an immediate
  reset the redone step does not recompute or re-shift. The `zip` block (599-609) runs on every
  call (idempotent zeroing).
* Linear runs: `add_explicit_terms` calls `save_dt_cfl(dt_cfl_default_large = 1.e8)`
  (nonlinear_terms.fpp:429-430, 90) on every call.

### 1.10 `collisions_advance` (dist_fn.fpp:4130-4167)

```
4142 if (.not. split_collisions) return
4145 call vspace_derivatives(gnew, g, g_work, phi, bpar, phinew, bparnew, mode)
4164 call g_adjust(gnew, phinew, bparnew, from_g_gs2)
4165 call calculate_potentials_from_nonadiabatic_dfn(gnew, phinew, aparnew, bparnew)
4166 call g_adjust(gnew, phinew, bparnew, to_g_gs2)
```

Note the fields are recomputed after collisions; `remove_zonal_flows` is not re-applied.

### 1.11 Reset of the time step (`reset_time_step`, gs2_reinit.f90:104-225)

```
127  if (first) call init_reinit ; first = .false.
133  if (istep_last + 1 == istep) then nconsec = nconsec + 1 else nconsec = 0     ! SAVE istep_last=-1, nconsec=0
139  if (nconsec > 4 .and. abort_rapid_time_step_change) -> my_exit = .true.; return
145  if (code_dt/delt_adj <= code_dt_min) -> code_dt = code_dt_min; my_exit = .true.; return
         ! NOTE: tested even when the reset would INCREASE dt
157  reset_in = reset ; reset = .false.          ! so that timeadv(...,0) does not early-return
163  gnew = 0
166  init(current_init, override_timestep=29)   ! tears down 32 (fields), 31 (dist_fn_level_3), 30 (collisions)
169  original_time_step = code_dt
174  do while (attempts < 10)
176     if (code_dt*fac > code_dt_cfl) then           ; reduce_time_step:  code_dt = code_dt/delt_adj        (86-91)
179     else if (code_dt*fac < min(code_dt_max, code_dt_cfl/delt_adj/delt_cushion)) then
                                                     ; increase_time_step: code_dt = min(code_dt*delt_adj, code_dt_max) (96-101)
182     else exit
184     attempts = attempts + 1
189  if (attempts == 10) -> reset = reset_in; my_exit = .true.; return     ! init level left < full
198  call save_dt(code_dt)                              ! user_dt = code_dt
201  if (code_dt <= code_dt_min) -> code_dt = code_dt_min; my_exit = .true.; return
216  init(current_init, full)                           ! rebuild 30, 31, 32, then 34 restores state (§1.12)
220  istep_last = istep
223  reset = reset_in
```

Evaluation order of the thresholds: `code_dt_cfl/delt_adj/delt_cushion` = `(code_dt_cfl/delt_adj)/delt_cushion`.
`fac = 1.0` (multiplication by 1.0 is exact).

### 1.12 State restoration after a dt change (`set_initial_field_and_dist_fn_values`, gs2_init.fpp:292-375, override = .true.)

```
in_memory (reinit_knobs, default .false.):
  T: g = initval_ov%g (= gnew saved at reset time) ; gnew = g ; new_vnmult = saved vnmult    (321-325)
  F: ginit(restarted, ginitopt_restart_many) -> ginit_restart_many (init_g.f90:3984-4004)
        -> gs2_restore(g, scale, ...) (gs2_save.fpp:604-888): g, phinew.. from the restart files just
           written; phi/apar/bpar zeroed (861-863); g and phinew.. multiplied by `scale` (865-876)!
        -> gnew = g ; phi = phinew ...    ; back in ginit: scale = 1. (init_g.f90:1177)
     init_vnmult(new_vnmult) from file                                                   (326-327)
  set_vnmult(new_vnmult)                                                                  (329)
  in_memory: gexp_1..3 = saved copies (332-336)          [file path: restored in gs2_restore if include_explicit_source_in_restart]
if (force_maxwell_reinit)  (fields_knobs, default .true.)                                 (347)
  set_init_fields -> init_allfields_implicit: phinew.. = get_fields_direct_from_dfn(gnew); phi=phinew..; g=gnew
else
  in_memory: phinew/aparnew/bparnew = saved (or zero if field absent); phi=phinew; apar=aparnew; bpar=bparnew  (351-370)
  file:      fields already read in gs2_restore
```

Module state that is **not** touched by a reset (levels below 29 stay up): `kx_shift`,
`theta0_shift`, `jump`, `ikx_indexed`, `exb_first`, exb `istep_last`, nonlinear
`istep_last`, `gexp_*` (restored with identical values anyway), `kwork_filter`, `ieqzip`,
antenna amplitudes, hyper arrays (`D_res`, `hypervisc_filter`), code_dt_prev1/2, code_dt_old,
code_time, user_time, code_dt_cfl, code_dt_max, code_dt_min.

### 1.13 What is recomputed on a dt change (exhaustive)

Torn down and rebuilt by `init(…, override_timestep)` then `init(…, full)`:

| level | routine (file:line) | arrays recomputed | dt enters via |
|---|---|---|---|
| 30 collisions | finish_collisions (collisions.fpp:4165-) / init_collisions (436-463 → init_arrays) | `vnew*`, `delvnew`, `vnewh`, `z0,w0,s0`, `bz0,bw0,bs0`, `c1,betaa,ql,d1,h1`, `ec1,ebetaa,eql` (or `*le` versions), `vpdiff`, `dtot,fdf,fdb`, `pitch_weights`, `c_rate` | `code_dt` in z0tmp/s0tmp/w0tmp (758-761,877-880,999), bz0tmp etc. (1281,1380,1495), energy diffusion tridiagonal (2078-2161), Lorentz tridiagonal (2469-2620); `vnmult` preserved through `vnm_init` (4174-4175, 1731-1733) |
| 31 dist_fn_level_3 | init_dist_fn_level_3 (dist_fn.fpp:960-1013) | `vpar` (init_vpar 1681-1729), `wdrift`, `wdriftttp` (init_wdrift 1405-1495 via wdrift_func 1500-1530 and wcoriolis_func 1534-1557), `wstar` (init_wstar 1732-1750), `a, b, r, ainv` (init_invert_rhs 1826-1929), `g_h`, `g_adj` (init_homogeneous_g 1937-2125), `source_coeffs_*` if `opt_source` (init_source_term 2159-2270) | `code_dt*tunits`, `code_dt*wunits`, `code_dt/2` |
| 32 fields_level_2 | reset_init (fields.f90:608-637) / init_fields (353-377) | implicit: `aminv` (init_response_matrix); local: `fieldmat` blocks; fields zeroed then restored (§1.12); `kwork_filter(1,1)=.true.` reasserted | through `timeadv(…,0)` with the new coefficients |
| 34 set_initial_values | §1.12 | `g, gnew, phi.., phinew..` | fields recomputed from g |

Recomputed every step anyway (runtime, no init): AB coefficients
(`get_adams_bashforth_coefficients`, gs2_time.f90:250-312, called from dist_fn.fpp:5388-5389),
`gdt` in exb_shear (dist_fn.fpp:4283), antenna `dt` (antenna.f90:365), iso/aniso hyperviscosity
filter (hyper.f90:738-743, 777-778), heating diagnostics `1/code_dt`.

**Not recomputed although dt-dependent (stale after a dt change — GS2 behaviour to replicate):**
* `D_res(it,ik) = D_hyperres*code_dt*(aky^2+akx^2)^nexp/akperp4_max` (hyper.f90:400-406), used in the
  source term (dist_fn.fpp:2227,2252,5087,5097). `init_hyper` returns early (`initialized`, 370-371)
  because level 19 is not torn down.
* `hypervisc_filter` for `hyper_option = 'simple3D'` (hyper.f90:419-445), used as-is by
  `simple3Dfilter` (677-683).
* `kwork_filter` is *not* reset to `.false.` before an `implicit` rebuild (see §6.11).

### 1.14 CFL detection (where `reset` is set)

* `run_parameters::reset` (run_parameters.f90:29, initial `.false.`) is the global flag.
* immediate_reset = .true. (default): in `add_explicit` after computing the NL term, if
  `use_cfl_limit .or. use_order_based_error`: `max_allreduce(max_vel_components)`
  (nonlinear_terms.fpp:567) then `set_cfl_time_step_and_check_if_too_large` (1229-1252):
  `max_vel = get_max_vel(...)` (1214-1224: `use_2d_cfl` ? c1+c2 : max(c1,c2); then max with c3);
  `dt_cfl = max_vel <= 0 ? 1e8 : 1./max_vel`; `save_dt_cfl(dt_cfl)`;
  `reset = code_dt > code_dt_cfl` (gs2_time.f90:154-158). If neither limit is active the same
  routine is called with `max_vel = epsilon` (582). Then `if (reset) return` (585-587).
* immediate_reset = .false.: `nb_max_allreduce` started (573); finished in
  `check_time_step` → `nb_check_time_step_too_large` (1191-1197). `check_time_step`
  (gs2_reinit.f90:228-254) then: `reset = .false.`; `if (exit) return`;
  `nb_check…` (writes `run_parameters::reset`, which is *the same variable* passed as the
  dummy argument — aliasing, harmless because the next two lines recompute it);
  `if (code_dt*fac > code_dt_cfl) reset = .true.`;
  `if (code_dt*fac < min(code_dt_max, code_dt_cfl/delt_adj/delt_cushion)) reset = .true.`.
* `parameter_scan.f90:317` can also set `reset` but `update_scan_parameter_value` is commented
  out in the loop (gs2_main.fpp:801) — GATE parameter scans.

Linear runs: `code_dt_cfl = 1e8` (re-saved every timeadv) so `check_time_step` resets only if
`code_dt < min(code_dt_max, 1e8/delt_adj/delt_cushion)`, i.e. if `code_dt < code_dt_max`. With
default `dt0` (= initial `delt`) this never fires; it does fire if the user sets `dt0 > delt`
or restarts with `delt_option='check_restart'` and a restart `delt_max > delt0`.

### 1.15 Stopping conditions (in loop order)

1. `istep > nstep` (gs2_main.fpp:718) — only on repeated `evolve_equations` calls.
2. Exit set by `reset_time_step` (rapid change, dt < delt_minimum, >10 attempts), checked at 754.
3. Diagnostics exit (run_diagnostics): omega converged with `exit_when_converged`
   (diagnostics_omega.f90:117-124; default `.true.`, `omegatol = 1e-3`, `navg = 10`) — linear
   runs only (`calculate_omega` returns early for nonlinear runs without ieqzip, 80-82); numerical
   instability `abs(omega_average)*code_dt > omegatinst` (126-132, default 1e6); nonlinear
   convergence check. Owned by the diagnostics spec; the time loop must expose a hook.
4. `.stop` file every `ncheck_stop` steps (job_manage.fpp checkstop), wall-clock `checktime`.
5. `user_time + user_dt > max_sim_time` after `update_time` and after any end-of-step reset (833).
6. Loop completes `nstep` steps.

### 1.16 Antenna (`antenna_amplitudes`, antenna.f90:341-443) — GATED

Called at the start of every `advance_*` (also again on the redone step after an immediate
reset, so the amplitudes are advanced twice in that step). Uses `ranf()` (intrinsic
`random_number`) on proc0, `code_dt`, `code_time`, `w_antenna + w_dot*(time-t0)`,
`gradpar(0)`, `beta_s`, and writes `apar_ext`, `apar_new/apar_old`, `a_ant/b_ant`. The random
stream cannot be reproduced in JAX → GATE (`driver` namelist absent or `ant_off = .true.`).

---------------------------------------------------------------------------------------------------

## 2. Inventory of arrays and scalars read by this subsystem

Legend: **E** = EXPORT from Fortran, **R** = RECOMPUTE in JAX (formula given), **B** = BOTH
(export at the initial dt for verification, recompute in JAX after a dt change), **S** = runtime
state carried by the JAX loop (export the t = 0 value).

### 2.1 Time-control scalars (module `gs2_time`, gs2_time.f90:19-37)

| name | type | set by | dt-dep / runtime | rec. | notes |
|---|---|---|---|---|---|
| `code_dt` | real | init_delt (140-151) from `delt` or restart; reduce/increase (gs2_reinit.f90:90,100); save_dt (206-213) | runtime | S | the step size |
| `user_dt` | real | init_delt, save_dt | runtime | S | equals `code_dt` except transiently inside reset_time_step (user2code is identity, 223-228) |
| `code_dt_prev1`, `code_dt_prev2` | real | init_delt (`dt_not_set = -1.0` unless check_restart); update_time (164-165) | runtime | S | sentinel −1.0 drives AB order |
| `code_dt_old` | real | module init `0.` (31); update_time (169) | runtime | S | used by exb `gdt`; 0 at first step even after restart |
| `code_time`, `user_time` | real | init_tstart (129-136) = `tstart` (init_g, or restart t0); update_time (170-171) | runtime | S | accumulated by repeated addition |
| `code_dt_cfl`, `user_dt_cfl` | real | module init `1e8` (23-24); save_dt_cfl (177-183) | runtime | S | |
| `code_dt_min` | real | save_dt_min(delt_minimum) (gs2_reinit.f90:300) | const | E | default 1e-5 |
| `code_dt_max` | real | init_gs2_time (65) = delt or restart delt_max; save_dt_max(dt0) (gs2_reinit.f90:301) | const | E | |
| `tunits(ik)`, `wunits(ik)`, `woutunits(ik)` | real(naky) | adjust_time_norm (83-126) | const | E (or R) | `wstar_units`: tunits=2/aky (0 where aky=0), wunits=1; else tunits=1, wunits=aky/2.0 |
| `dt_not_set` | parameter −1.0 | (19) | const | R | |

### 2.2 Reset-control scalars (module `gs2_reinit`, gs2_reinit.f90:16-23)

| name | default | rec. | notes |
|---|---|---|---|
| `delt_adj` | 2.0 (40) | E | |
| `delt_cushion` | 1.5 (47) | E | |
| `delt_minimum` | 1e-5 (50) | E | → code_dt_min |
| `dt0` | smart default = code_dt_max (282-286) | E | → code_dt_max |
| `abort_rapid_time_step_change` | .true. (31) | E | |
| `in_memory` | .false. (61) | E | only affects the scale quirk (§6.12) |
| `reset_time_step::istep_last` | −1 (SAVE, 120) | S | |
| `reset_time_step::nconsec` | 0 (SAVE, 121) | S | |
| `timestep_change_attempts_limit` | 10 (126) | R | |

### 2.3 Run parameters (module `run_parameters`)

| name | shape | set at | rec. | notes |
|---|---|---|---|---|
| `nstep` | int | read_parameters (322) | E | default 100 |
| `max_sim_time` | real | (319) | E | default 1e6 |
| `immediate_reset` | logical | (316), override (389) | E | default .true. |
| `reset` | logical | module (29) | S | initial .false. |
| `has_phi/has_apar/has_bpar` | logical | (351-353) `f* > epsilon(0.0)`, after beta logic (342-366) | E | |
| `fphi, fapar, fbpar, beta` | real | (308-315, 342-366) | E | used by g_adjust etc. (other specs) |
| `ieqzip(it,ik)` | logical(ntheta0,naky) → NumPy (naky,ntheta0) | init_run_parameters (259-272) | E | masks in invert_rhs_1 (dist_fn.fpp:5477), collisions |
| `ncheck_stop`, `avail_cpu_time`, `margin_cpu_time` | | (320,307,318) | E | host-side only |
| `wstar_units` | logical | (331) | E | via tunits/wunits |

### 2.4 Distribution-function and field state (S = must be exported at t = 0)

| name | module | Fortran shape | NumPy shape | set at init | rec. | notes |
|---|---|---|---|---|---|---|
| `g` | dist_fn_arrays | (-ntgrid:ntgrid,2,iglo) complex | (nspec,negrid,nlambda,naky,ntheta0,2,2nt+1) | ginit / set_init_fields (`g = gnew`) | S | == gnew at t=0 |
| `gnew` | dist_fn_arrays | same | same | ginit | S | **the** state |
| `gexp_1, gexp_2, gexp_3` | dist_fn_arrays | same (allocated iff `nonlin .and. .not. split_nonlinear`, dist_fn.fpp:3970-3980) | same | zero, or restart (gs2_save include_explicit_source_in_restart) | S | NL history, shifted once per step |
| `g_work` | dist_fn_arrays | same | same | zero | scratch | zeroed in invert_rhs (6179) before use |
| `phi, apar, bpar` | fields_arrays | (-ntgrid:ntgrid,ntheta0,naky) complex | (naky,ntheta0,2nt+1) | set_init_fields copy | S | == phinew at t=0 |
| `phinew, aparnew, bparnew` | fields_arrays | same | same | get_init_field / restart | S | |
| `apar_ext` | fields_arrays | same | same | antenna | GATE | |
| `kx_shift(ik)` | dist_fn_arrays | (naky) real, allocated iff `abs(g_exb*g_exbfac)>eps` and (`box` or `shat==0`) | (naky) | 0 or restart (gs2_save.fpp:820-824) | S | |
| `theta0_shift(ik)` | dist_fn_arrays | (naky), allocated iff exb and not box and shat≠0 | (naky) | 0 or restart | S | |
| `jump(ik)` | dist_fn | int(naky) | (naky) | exb_first block | S (scratch) | recomputed each shear step |
| `vnmult(2)` | collisions | real(2) | (2) | init_collisions (1731-1733), restart | S | constant at runtime (adjust_vnmult never called) |
| `kwork_filter(it,ik)` | kt_grids | logical(ntheta0,naky) | (naky,ntheta0) | `.false.` (kt_grids.f90:189); fields_local manipulates during build (fields_local.f90:2006-2076); `(1,1)=.true.` for box (fields.f90:376) | S/E | export **after** init_fields |

### 2.5 Hidden SAVE / module counters (S)

| name | file:line | initial | role |
|---|---|---|---|
| `timeadv::istep_last` | dist_fn.fpp:4031 | −1 (becomes 0 after response build) | split-NL half-step logic |
| `nonlinear_terms::istep_last` | nonlinear_terms.fpp:32 | 0 | NL term once per step |
| `nonlinear_terms::initializing` | nonlinear_terms.fpp:80 | .true. → .false. at nl_finish_init | |
| `exb_shear::istep_last` | dist_fn.fpp:4217 | 0 | shear applied once per step |
| `exb_first` | dist_fn.fpp:238 | .true. | one-time setup of ikx_indexed, dkx/dtheta0 |
| `exb_shear::itmin, dkx, dtheta0` | dist_fn.fpp:4212,4219 | set in exb_first block | |
| `ikx_indexed(ntheta0)` | dist_fn.fpp:166 | set in exb_first block | |
| `cfl_req_hand` | nonlinear_terms.fpp:64 | null | MPI only |

### 2.6 Geometry/grid inputs read here (E)

| name | module | shape | used in |
|---|---|---|---|
| `aky(ik)`, `akx(it)`, `ikx(it)` | kt_grids | (naky), (ntheta0), int(ntheta0) | exb_shear |
| `theta0(it,ik)` | kt_grids | (ntheta0,naky) → (naky,ntheta0) | exb (`dtheta0 = theta0(2,1)-theta0(1,1)`, 4262) |
| `box` (`gridopt_switch == gridopt_box`), `reality` | kt_grids | logical | exb branch, kwork_filter, antenna |
| `shat` | theta_grid | real | exb branch (`should_use_kx = box .or. shat == 0.0`, 4228) |
| `ntheta`, `ntgrid` | theta_grid | int | exb theta0 branch (`ileft=-ntgrid+ntheta`, `iright=ntgrid-ntheta`, 4354) |
| `theta(ig)`, `jacob(ig)`, `field_line_average_weight` | theta_grid | (2nt+1), scalar (theta_grid.f90:520) | remove_zonal_flows |
| `g_exb, g_exbfac, g_exb_start_time, g_exb_start_timestep` | dist_fn | real, int | exb_shear |
| `hyper_in_initialisation`, `split_collisions` | dist_fn / collisions | logical | timeadv |
| `nonlin`, `split_nonlinear`, `strang_split` | nonlinear_terms | logical | timeadv |
| `remove_zonal_flows_switch`, `force_maxwell_reinit`, `fieldopt_switch` | fields | | advance, reset |
| `new_field_init`, `scale`, `tstart`, `ginitopt_switch` | init_g | | reset path |

### 2.7 dt-dependent coefficient arrays (owned by other specs; listed for the reset contract)

All are **B** (export at the initial dt; recompute in JAX on a dt change, or look up an exported
dt ladder, §5.5). Formulas quoted where short.

| name | module | shape (canonical) | init routine | formula (dt part) |
|---|---|---|---|---|
| `vpar` | dist_fn | g-like real | init_vpar dist_fn.fpp:1681-1729 | `vpar(ig,isgn) = spec(is)%zstm*tunits(ik)*code_dt*0.5/delthet(ig)*(abs(gradpar(ig))+abs(gradpar(ig+1)))*vpac(ig,isgn)` for ig<ntgrid; `vpar(ntgrid,:)=0` (evaluated left to right) |
| `wdrift`, `wdriftttp` | dist_fn | g-like real | init_wdrift 1405-1495 | `wdrift_func` (1500-1530): aky==0: `akx/shat*(cvdrift0*vpa^2 + gbdrift0*0.5*vperp2)*code_dt/2.0`; else `((cvdrift+theta0*cvdrift0)*vpa^2 + (gbdrift+theta0*gbdrift0)*0.5*vperp2)*code_dt*wunits`; Coriolis (1534-1557); then `wdrift(ig)=0.5*(wdrift(ig)+wdrift(ig+1))` for ig<ntgrid (1486-1490) |
| `wstar` | dist_fn_arrays | Fortran (naky,negrid,nspec) | init_wstar 1732-1750 | `wstar(ik,ie,is) = - code_dt*wunits(ik) * dlnf0drho(ie,is)` |
| `a, b, r, ainv` | dist_fn | g-like complex | init_invert_rhs 1826-1929 | from `wdrift, wdriftttp, vpar, bkdiff, fexp, spec%tz` (1885-1921) |
| `g_h`, `g_adj` | dist_fn | g-like complex; `g_adj` (n_links,2,iglo) | init_homogeneous_g 1937-2125 | recurrences in `r` (2071-2073, 2108-2114) |
| `source_coeffs_phim/phip/aparm/aparp` | dist_fn | (-nt:nt-1,2,iglo) | init_source_term 2159-2270 (only if `opt_source`) | `code_dt*wunits`, `D_res` |
| collision coefficients | collisions | various | init_collisions | see collisions spec |
| `aminv` / `fieldmat` | fields_implicit / fields_local | per supercell | init_response_matrix / init_fields_matrixlocal | via `timeadv(…,0)` |
| `D_res(it,ik)` | hyper | (ntheta0,naky) | init_hyper 400-406 | **stale**: keep initial-dt values |
| `hypervisc_filter` (simple3D) | hyper | (-nt:nt,ntheta0,naky) | init_hyper 419-445 | **stale**: keep initial-dt values |

---------------------------------------------------------------------------------------------------

## 3. Exact pseudo-code (order of operations and loop order preserved)

### 3.1 JAX-side main loop (mirrors gs2_main.fpp:717-861 without file I/O)

```
state  = load_export_t0()            # §3.10
coeffs = load_export_coeffs(code_dt) # §2.7
for istep in 1 .. nstep:
    # ---- advance(istep) ----
    state, pending_reset = advance(state, coeffs, istep)            # §3.2
    if pending_reset:                                               # immediate_reset path
        state, exit = reset_time_step(state, istep)                 # §3.6 (rebuilds coeffs)
        if not exit:
            state, _ = advance(state, coeffs, istep)                # cannot raise reset again (NL skipped)
    if exit: break
    # (restart writing: optional host feature, uses user_time BEFORE update)
    update_time(state)                                              # §3.4
    exit = diagnostics_hook(state, istep)                           # owned by diagnostics spec
    if not exit:
        reset = check_time_step(state, exit)                        # §3.5
        if reset:
            state, exit = reset_time_step(state, istep)             # no re-advance
        if istep % ncheck_stop == 0 and not exit: exit = stop_file_exists()
        if not exit: exit = wall_clock_exceeded()
    if not exit and state.user_time + state.user_dt > max_sim_time: exit = True   # uses post-reset user_dt
    if exit: break
```

### 3.2 `advance` (implicit / local, single process)

```
def advance(S, C, istep):
    # 1. antenna: GATED (no_driver must be .true.)
    # 2. flow shear (if kx_shift or theta0_shift allocated)
    if exb_active: S = exb_shear(S, istep)                          # mutates gnew, phinew, aparnew, bparnew (§3.7)
    # 3. copies
    S.g = S.gnew
    implicit: S.phi, S.apar, S.bpar = S.phinew, S.aparnew, S.bparnew
    local:    copy only fields with has_* (identical result)
    # 4. predictor
    S, reset = timeadv(S, C, phi=S.phi.., new=S.phinew.., istep, mode=None)   # §3.3
    if reset: return S, True                                         # S.gnew == S.g, gexp shifted, NL computed
    # 5. field solve (field spec): returns delta in phinew..
    S.phinew, S.aparnew, S.bparnew = getfield(S.phinew, S.aparnew, S.bparnew)
    # 6. add old fields (implicit: all three; local: has_* only)
    S.phinew  = S.phinew  + S.phi
    S.aparnew = S.aparnew + S.apar
    S.bparnew = S.bparnew + S.bpar
    # 7. zonal flow removal
    if remove_zonal_flows_switch: S.phinew = remove_zonal_flows(S.phinew)   # §3.8
    # 8. corrector
    S, _ = timeadv(S, C, phi=S.phi.., new=S.phinew.., istep, mode=1)
    # 9. split collisions
    if split_collisions: S = collisions_advance(S, C, mode=1)               # §1.10
    return S, False
```

### 3.3 `timeadv` step-level skeleton

```
def timeadv(S, C, phi, new, istep, mode):
    first_half = (istep != S.timeadv_istep_last)
    # explicit terms
    if nonlin:                                   # (split_nonlinear GATED)
        if istep != 0 and istep != S.nl_istep_last:
            S.nl_istep_last = istep
            S.gexp_3 = S.gexp_2 ; S.gexp_2 = S.gexp_1            # copy(g2,g3); copy(g1,g2) — order matters
            S.gexp_1 = NL(S.g, phi..)            # add_nl + nl_center (nonlinear spec)
            if not use_cfl_limit: max_vel[0:2] = eps
            if use_order_based_error and istep > istep_error_start: max_vel[2] = 1/(code_dt*sqrt(error_target/error))
            if use_cfl_limit or use_order_based_error:
                if immediate_reset:
                    S.code_dt_cfl = cfl_from(max_vel)  ;  S.reset = S.code_dt > S.code_dt_cfl
                else:
                    S.pending_maxvel = max_vel          # resolved in check_time_step
            else:
                S.code_dt_cfl = cfl_from(max_vel) ; S.reset = S.code_dt > S.code_dt_cfl
            if S.reset: (skip zip) -> fallthrough to return below
        if istep != 0 and not S.reset and zip:                              # every call with istep>0
            zero S.g and S.gexp_1 for all ik != 1   (nonlinear_terms.fpp:599-609)
    else:
        S.code_dt_cfl = 1e8 ; S.user_dt_cfl = 1e8                           # linear (every call)
    if S.reset and immediate_reset: return S, True                          # timeadv_istep_last NOT updated
    S.gnew = invert_rhs(S, C, phi, new, istep, sourcefac(S.code_time))      # spec 02/03
    if (istep == 0 and hyper_in_initialisation) or istep > 0: S.gnew = hyper_diff(S.gnew, new.phi)
    if not split_collisions: S.gnew = vspace_derivatives(S.gnew, S.g, phi, bpar, new.phi, new.bpar, mode)
    if def_parity: enforce_parity(S.gnew)
    S.timeadv_istep_last = istep
    return S, False
```

Note `add_explicit`'s exact control flow: when `reset` becomes true at 585-587 it returns
**before** the `zip` block and before line 611 (`istep_last = istep`, already set at 497).

### 3.4 `update_time` (gs2_time.f90:161-174)

```
code_dt_prev2 = code_dt_prev1
code_dt_prev1 = code_dt
code_dt_old   = code_dt
code_time     = code_time + code_dt
user_time     = user_time + user_dt
```

### 3.5 `check_time_step` (gs2_reinit.f90:228-254)

```
reset = False
if exit: return reset
if not immediate_reset and pending_maxvel is not None:
    code_dt_cfl = cfl_from(pending_maxvel)          # nb_check_time_step_too_large
    (reset = code_dt > code_dt_cfl  — overwritten below, harmless)
if code_dt*1.0 > code_dt_cfl: reset = True
if code_dt*1.0 < min(code_dt_max, (code_dt_cfl/delt_adj)/delt_cushion): reset = True
return reset
```

`cfl_from(mv)`: `max_vel = (use_2d_cfl ? mv[0]+mv[1] : max(mv[0],mv[1]))`; `max_vel = max(max_vel, mv[2])`;
`dt_cfl = 1e8 if max_vel <= 0 else 1./max_vel` (nonlinear_terms.fpp:1214-1246).

### 3.6 `reset_time_step` (scalar part, exact) + rebuild

```
def reset_time_step(S, istep):
    if R.istep_last + 1 == istep: R.nconsec += 1 else: R.nconsec = 0
    if R.nconsec > 4 and abort_rapid_time_step_change: return S, exit=True
    if S.code_dt/delt_adj <= code_dt_min: S.code_dt = code_dt_min; return S, exit=True
    S.gnew_saved = S.gnew                       # (Fortran zeroes gnew then restores it)
    attempts = 0
    while attempts < 10:
        if S.code_dt > S.code_dt_cfl:                                   S.code_dt = S.code_dt/delt_adj
        elif S.code_dt < min(code_dt_max, (S.code_dt_cfl/delt_adj)/delt_cushion):
                                                                        S.code_dt = min(S.code_dt*delt_adj, code_dt_max)
        else: break
        attempts += 1
    if attempts == 10: return S, exit=True      # user_dt NOT updated
    S.user_dt = S.code_dt                       # save_dt
    if S.code_dt <= code_dt_min: S.code_dt = code_dt_min; return S, exit=True
    C = rebuild_coeffs(S.code_dt)               # §1.13 table; response matrix with current kwork_filter (§6.11)
    S = restore_state(S)                        # §3.9
    R.istep_last = istep
    return S, exit=False
```

### 3.7 `exb_shear` (dist_fn.fpp:4183-4600), single-process form

```
# one-time setup (exb_first)
if S.exb_first:
    S.exb_first = False ; jump[:] = 0
    if should_use_kx:                                       # box .or. shat == 0.0
        itmin = argmin_first(ikx)                           # minloc: first occurrence, 1-based
        for it in itmin..ntheta0:      ikx_indexed[it+1-itmin] = it
        for it in 1..itmin-1:          ikx_indexed[ntheta0-itmin+1+it] = it
        dkx = akx[2]-akx[1]  if ntheta0 > 1 else <unset>   # GATE ntheta0 == 1
    else:
        dtheta0 = theta0[2,1]-theta0[1,1] if ntheta0 > 1 else 2*pi
# guards
if istep == S.exb_istep_last: return
if g_exb_start_timestep > istep: return
if g_exb_start_time >= 0 and S.code_time < g_exb_start_time: return
gdt = 0.5*(S.code_dt + S.code_dt_old)
S.exb_istep_last = istep
if should_use_kx:
    for ik in 1..naky:
        kx_shift[ik] = kx_shift[ik] - aky[ik]*g_exb*g_exbfac*gdt*tunits[ik]     # (((aky*g_exb)*g_exbfac)*gdt)*tunits
        jump[ik]     = NINT(kx_shift[ik]/dkx)                                   # round half AWAY from zero
        kx_shift[ik] = kx_shift[ik] - real(jump[ik])*dkx
else:
    for ik in 1..naky:
        theta0_shift[ik] = theta0_shift[ik] - g_exb*g_exbfac*gdt/shat*tunits[ik]   # (((g_exb*g_exbfac)*gdt)/shat)*tunits
        jump[ik]         = NINT(theta0_shift[ik]/dtheta0)
        theta0_shift[ik] = theta0_shift[ik] - dtheta0*real(jump[ik])
if all(jump == 0): return
```

kx branch (`should_use_kx`), `for ik = naky, 2, -1` (**ik = 1 never shifted**):

```
j = jump[ik]
pos(it) := ikx_indexed[it]      # it-index of the it-th kx in ascending-kx order
if j < 0:  for it in 1..ntheta0+j:            X[pos(it), ik] = X[pos(it-j), ik]
           for it in ntheta0+j+1..ntheta0:    X[pos(it), ik] = 0
if j > 0:  for it in ntheta0 down to 1+j:     X[pos(it), ik] = X[pos(it-j), ik]
           for it in j down to 1:             X[pos(it), ik] = 0
```
applied, in this order, to `phinew` (if has_phi), `aparnew` (if has_apar), `bparnew` (if has_bpar),
then to `gnew` for every `(is, ie, il)` (both sigma, all theta). Because the source index is
always ahead of the destination in the loop direction, this equals the out-of-place gather
`Xnew[pos(it)] = (1 <= it-j <= ntheta0) ? Xold[pos(it-j)] : 0`.

theta0 branch (`.not. should_use_kx`), `for ik = 1..naky` (ascending, **order matters**):

```
j = jump[ik]; if j == 0: continue
if abs(j) >= ntheta0: ABORT ("jump(ik) > ntheta0 ... reduce timestep or increase ntheta0")
iib = |j| ; iit = ntheta0+1-|j| ; ileft = -ntgrid+ntheta ; iright = ntgrid-ntheta
for F in [phinew if has_phi, aparnew if has_apar, bparnew if has_bpar]:
    if j < 0:
        temp2 = F[:, 1:iib, ik]                              # copy
        for i in 1..iit-1:   F[:, i, ik] = F[:, i-j, ik]
        F[ileft:ntgrid, iit:ntheta0, ik]  = temp2[-ntgrid:iright, :]
        F[-ntgrid:ileft-1, iit:ntheta0, ik] = 0
    else:
        temp2 = F[:, iit:ntheta0, ik]
        for i in ntheta0 down to iib+1:   F[:, i, ik] = F[:, i-j, ik]
        F[-ntgrid:iright, 1:iib, ik] = temp2[ileft:ntgrid, :]
        F[iright+1:ntgrid, 1:iib, :] = 0            # <-- ALL ky, not just ik (dist_fn.fpp:4370,4387,4404)
for is, ie, il:                                      # gnew, both sigma
    if j < 0:
        temp[:, :, it] = gnew[it] for it in 1..iib
        for it in 1..iit-1:          gnew[it] = gnew[it-j]
        for it in iit..ntheta0:      gnew[ileft:, :, it] = temp[:iright, :, it-iit+1]; gnew[:ileft-1, :, it] = 0
    else:
        temp[:, :, it] = gnew[iit+it-1] for it in 1..j
        for it in ntheta0 down to iib+1: gnew[it] = gnew[it-j]
        for it in 1..iib:                gnew[:iright, :, it] = temp[ileft:, :, it]; gnew[iright+1:, :, it] = 0
```
(“`gnew[it]`” = all theta, both sigma, at that `(ik, il, ie, is)`; the zeroing of g in the j>0
branch is per-ik — only the field statements carry the all-ky bug.)

### 3.8 `remove_zonal_flows` (fields_implicit.f90:374-425, theta_grid.f90:665-670, externals/utils/integration.f90)

```
for it in 1..ntheta0:                        # ik = 1 only; all other ky untouched (avg array zero)
    q   = phinew[:, it, 1] * jacob[:]        # complex*real
    d   = theta[2:] - theta[:-1]             # real, length 2*ntgrid
    y   = d * (q[:-1] + q[1:])               # real*complex
    I   = complex(klein_sum(real(y)), klein_sum(imag(y))) * 0.5
    avg = I * field_line_average_weight      # weight = 1.0/trapezoidal_integration(theta, jacob)  (theta_grid.f90:520)
    phinew[:, it, 1] = phinew[:, it, 1] - avg
```
`klein_sum` (externals/utils/summation.f90) is a sequential second-order compensated
(Kahan–Babuška–Klein) sum; implement as a `lax.scan` over the 2*ntgrid terms, final
`(result + c1) + c2`.

### 3.9 State restoration after a dt change

```
def restore_state(S):
    S.g = S.gnew_saved ; S.gnew = S.g                 # gexp_*, vnmult, kx_shift unchanged
    if (not in_memory) and scale_pending:             # §6.12 quirk; scale_pending is True only before the
        if scale > 0:  g*=scale; phinew*=scale; aparnew*=scale; bparnew*=scale    # first file-based reset
        else: fac = -scale/max(abs(phinew)); g*=fac; phinew*=fac; ...             # and only if the initial
        S.gnew = S.g ; scale_pending = False                                       # ginit was not restart-type
    if force_maxwell_reinit:                         # default
        (requires new_field_init == .true.)
        S.phinew, S.aparnew, S.bparnew = get_fields_direct_from_dfn(S.gnew)   # field spec
    else:
        absent fields set to 0; present fields keep saved values
    S.phi, S.apar, S.bpar = S.phinew, S.aparnew, S.bparnew
    S.g = S.gnew
    return S
```

### 3.10 Adams–Bashforth coefficients (gs2_time.f90:234-312), exact expressions

```
order = 3 if prev2 != -1.0 else (2 if prev1 != -1.0 else 1)      # time_history_available
this_dt = code_dt
order 3:
  c1 = (((this_dt + prev1)*(this_dt + prev1)/prev1)*((this_dt+prev1)/3. + prev2/2.)
         - prev1*(prev1/3. + prev2/2.)) / (this_dt*(prev1 + prev2))
  c2 = - this_dt*(this_dt/3. + (prev1+prev2)/2.) / (prev1*prev2)
  c3 = this_dt*(this_dt/3. + prev1/2.) / (prev2*(prev1+prev2))
order 2: c1 = (0.5*this_dt + prev1)/prev1 ; c2 = -0.5*this_dt/prev1
order 1: c1 = 1.0
```
`(x)**2` is computed by gfortran as `x*x`. `/3.` and `/2.` are true divisions. Compute these
scalars on the host (Python float) or in float64 on device without algebraic rewriting; they
are evaluated in every `timeadv` call with `istep > 0` (twice per step, same values).

---------------------------------------------------------------------------------------------------

## 4. Namelist / configuration options that change behaviour here

| option | namelist | default | effect | decision | reason |
|---|---|---|---|---|---|
| `nstep` | knobs | 100 | number of steps | SUPPORT | |
| `delt` | knobs | 0.1 | initial dt (and default dt0) | SUPPORT | |
| `delt_option` | knobs | 'default' | 'check_restart' reads dt history, delt_max, vnmult from restart | SUPPORT | init only; export captures result |
| `max_sim_time` | knobs | 1e6 | stop when next step would exceed | SUPPORT | |
| `immediate_reset` | knobs | .true. | CFL-too-large handling (redo step vs. finish step) | SUPPORT | both paths |
| `avail_cpu_time`, `margin_cpu_time` | knobs | 1e10, 300 | wall-clock stop | SUPPORT | host-side, non-deterministic |
| `ncheck_stop` | knobs | 5 | `.stop` file polling | SUPPORT | host-side |
| `progress_frequency` | knobs | 0 | cosmetic | SUPPORT | no numerics |
| `wstar_units` | knobs | .false. | tunits/wunits | SUPPORT | via exported tunits/wunits |
| `fphi, fapar, fbpar, beta` | knobs | 1, 0, −1, 0 | has_* fields; beta=0 disables apar/bpar | SUPPORT | |
| `eqzip_option` | knobs | 'none' | ieqzip masks (secondary/tertiary/equilibrium) | SUPPORT | simple mask |
| `do_eigsolve` | knobs | .false. | eigensolver instead of time loop | GATE | not a time stepper |
| `use_old_diagnostics` | knobs | .false. | which diagnostics run (round trips, exits) | SUPPORT | diagnostics spec |
| `seed` | knobs | 0 | RNG seed | SUPPORT | only affects init (exported) / antenna (gated) |
| `delt_adj` | reinit_knobs | 2.0 | dt change factor | SUPPORT | |
| `delt_cushion` | reinit_knobs | 1.5 | increase hysteresis | SUPPORT | |
| `delt_minimum` | reinit_knobs | 1e-5 | abort threshold | SUPPORT | |
| `dt0` | reinit_knobs | = code_dt_max (= delt) | max dt | SUPPORT | |
| `abort_rapid_time_step_change` | reinit_knobs | .true. | abort after >4 consecutive-step resets | SUPPORT | |
| `in_memory` | reinit_knobs | .false. | reset via memory vs restart files | SUPPORT | identical numerics except scale quirk §6.12 |
| `field_option` | fields_knobs | 'default' (=implicit) | field algorithm | SUPPORT 'default','implicit','local','implicit_local'; GATE 'gf_local','test' | gf_lo integration order differs; test is debug |
| `remove_zonal_flows_switch` | fields_knobs | .false. | subtract ky=0 field-line average of phinew | SUPPORT | |
| `force_maxwell_reinit` | fields_knobs | .true. | recompute fields from g after a reset/restart | SUPPORT | |
| `dump_response`, `read_response`, `response_dir`, `response_file` | fields_knobs | .false./'' | response I/O | SUPPORT | JAX uses exported matrix |
| `field_subgath`, `do_smart_update`, `minnrow`, `field_local_allreduce(_sub)`, `field_local_nonblocking_collectives`, `field_local_tuneminnrow` | fields_knobs | | parallel layout of field solve | SUPPORT | no single-process numerics effect |
| `g_exb` | dist_fn_knobs | 0.0 | flow shear | SUPPORT | requires ntheta0 > 1 (GATE otherwise) |
| `g_exbfac` | dist_fn_knobs | 1.0 | multiplies g_exb in the perpendicular shift | SUPPORT | |
| `g_exb_start_time`, `g_exb_start_timestep` | dist_fn_knobs | −1, −1 | delay shear | SUPPORT | |
| `g_exb_error_limit` | dist_fn_knobs | 0.0 | unused in code | SUPPORT | no effect |
| `hyper_in_initialisation` | dist_fn_knobs | .true. | hyperviscosity in response build | SUPPORT | |
| `gf_lo_integrate` | dist_fn_knobs | .false. | gf_lo integrals | GATE | |
| `nonlinear_mode` | nonlinear_terms_knobs | 'default'(off) | NL term, CFL | SUPPORT | |
| `split_nonlinear` | nonlinear_terms_knobs | .false. | operator-split NL; rewrites code_dt_prev1/2 (split_nonlinear_terms.fpp:313-391) | GATE | experimental, different control flow |
| `strang_split`, `advance_nonadiabatic_dfn` | split_nonlinear_terms_knobs | .false. | only with split_nonlinear | GATE | |
| `cfl`, `use_cfl_limit`, `use_2d_cfl` | nonlinear_terms_knobs | 0.95, .true., .true. | CFL limit | SUPPORT | nonlinear spec |
| `use_order_based_error`, `error_target`, `istep_error_start` | nonlinear_terms_knobs | .false., 0.1, 30 | error-based dt limit | SUPPORT | uses AB coeffs §3.10 |
| `zip` | nonlinear_terms_knobs | .false. | zero all ky≠1 in g, gexp_1 each call | SUPPORT | mask |
| `split_collisions` | collisions_knobs | .false. | collisions_advance path | SUPPORT | components ported anyway |
| `timesteps_between_collisions` | collisions_knobs | 1 | no effect (commented out, collisions.fpp:1777) | SUPPORT | no-op |
| `vary_vnew` | collisions_knobs | .false. | no numerical effect (adjust_vnmult never called) | SUPPORT | no-op |
| `ant_off` / `driver` namelist | driver | ant_off=.false., absent ⇒ off | Langevin antenna with random forcing | GATE (require off) | RNG not reproducible |
| `ginit_option` | init_g_knobs | 'default' | initial condition | SUPPORT | export post-init state |
| `tstart` | init_g_knobs | 0.0 | initial time | SUPPORT | exported user_time/code_time |
| `scale` | init_g_knobs | 1.0 | restart scaling; quirk on first file-based reset | SUPPORT if 1.0 or ginit restart-type or in_memory; GATE otherwise | §6.12 |
| `new_field_init` | init_g_knobs | .true. | field init method on reset | GATE if .false. | `local` cannot run it; non-standard |
| `initial_condition_is_nonadiabatic_dfn` | init_g_knobs | .false. | init only | SUPPORT | export |
| `restart_file`, `include_explicit_source_in_restart` | init_g / gs2_save | | restart I/O | SUPPORT | export captures |
| `nsave` | gs2_diagnostics_knobs | | write restart files | SUPPORT | host-side optional |
| `exit_when_converged`, `omegatol`, `navg`, `omegatinst`, `nwrite` | gs2_diagnostics_knobs | .true., 1e-3, 10, 1e6, 10 | diagnostics-driven exit; nwrite controls in-place round trips | SUPPORT | via diagnostics hook |
| `grid_option` | kt_grids_knobs | | box sets kwork_filter(1,1); box/shat=0 → kx exb branch | SUPPORT box/range/single/specified | |
| `hyper_option` | hyper_knobs | 'default' | iso/aniso recomputed each call; simple3D/D_res stale after dt change | SUPPORT | replicate staleness |
| `layout` | layouts_knobs | 'lxyes' | iglo ordering | SUPPORT | canonical reindex at export |
| `nensembles`, list mode, Trinity | (command/state) | | multiple jobs | GATE | host-level |
| parameter scans | parameter_scan_knobs | | (update call commented out) | GATE | |
| `opt_source` | dist_fn_knobs/optimisations | .false. | source coefficient precompute | SUPPORT | dt-dependent coeffs |

---------------------------------------------------------------------------------------------------

## 5. GPU mapping notes

### 5.1 Overall structure

* **Sequential**: the time loop (`istep`), the scalar time-control recurrences
  (`code_time`, `user_time`, `code_dt_prev*`, `kx_shift`, reset counters). Use a jitted
  `lax.while_loop` whose carry holds all S-state of §2; the predicate stops on
  `istep > nstep`, `exit`, or `pending_reset`. Return to Python for a reset (coefficient
  rebuild is a host-driven, shape-stable re-entry; the coefficient pytree is a jit argument so
  no retrace).
* **Data-dependent branches inside a step** (AB order, `istep == istep_last` guards, exb
  start conditions, `all(jump == 0)`, reset flags) → `lax.cond` / `jnp.where` on carry scalars.
  Static configuration (field_option, has_*, nonlin, split_collisions, should_use_kx, box,
  remove_zonal_flows_switch) → Python-level branches at trace time.
* **Independent (vmap/broadcast)**: everything indexed by `(is, ie, il, ik, it, isgn, ig)` in
  copies, `phinew + phi`, zonal removal (per it), kx-branch shifts (per ik, per (is,ie,il)).

### 5.2 Masks

* `kwork_filter[ik, it]` and `ieqzip[ik, it]`: modes where `invert_rhs_1` cycles
  (dist_fn.fpp:5475-5477) keep `gnew == g` (copied at the start of the step, after exb shift);
  apply as `gnew = where(mask, g, gnew_computed)`. Note collisions and other routines also
  consult these masks (collisions.fpp:3015,3068,…).
* exb out-of-range targets → zero fill masks.

### 5.3 `exb_shear` on GPU

* kx branch: per-ik integer shift `j[ik]` (carry), apply
  `take_along_axis` on the `ntheta0` axis with index `pos_inv[pos(it)-j]` and a validity mask;
  skip ik = 1 (`mask ik >= 2`). Fully vectorised over ky (rows are independent).
* theta0 branch: must be a `lax.fori_loop` over `ik` in ascending order (or an unrolled
  Python loop, naky is small) because of the all-ky zeroing bug for j > 0; within each ik the
  operation is a gather in (it, theta) with masks. Error flag for `|j| >= ntheta0` carried and
  checked on the host after each chunk (cannot raise inside jit).
* `NINT`: implement round-half-away-from-zero exactly:
  `r = jnp.round(x); r = jnp.where(jnp.abs(x - jnp.trunc(x)) == 0.5, jnp.trunc(x) + jnp.sign(x), r)`.

### 5.4 Response-matrix rebuild in JAX (if RECOMPUTE is chosen)

* Each `timeadv(…,0)` call in GS2 processes **all members `m` of one class simultaneously**
  plus the duplicated endpoint for `n > 1` (fields_implicit.f90:575-588). Operations that are
  nonlinear in the impulse (hyperviscosity `iso_shear`/`aniso_shear` shear rates sum
  `|phi|^2` over all modes, hyper.f90:695-778) couple the members, so to replicate exactly the
  JAX batch must reproduce the same impulse grouping: one batch element per
  `(class i, cell n, ig, ifield)` containing all M_class(i) impulses, not one per supercell.
  Without hyperviscosity (or with `const_amp`) every impulse is independent and can be vmapped
  freely.
* The `timeadv` used must be the `istep = 0` variant (§1.4), with `g = 0`, `phi = 0`, and
  `gnew` zero, `kwork_filter` as it is at that moment (§6.11).

### 5.5 Alternative to RECOMPUTE: exported dt ladder

dt only moves by `÷delt_adj`, `×delt_adj` capped at `code_dt_max`, or to `code_dt_min` (exit).
The set of reachable values is `{dt_init·adj^-k}` ∪ `{code_dt_max·adj^-k}` (exactly, when
`delt_adj` is a power of two; otherwise track the exact op sequence). GS2 already has a tool for
this (`src/programs/dump_response.f90:62-99`, `n_time_steps`, uses `reduce_time_step` +
`finish_dist_fn_level_3` + collisions/fields `reset_init`). The export program can dump all §2.7
arrays for `k = 0..K` (`K = ceil(log(dt_init/dt_min)/log(adj))`), keyed by the exact float64
bit pattern of `code_dt`; JAX looks up, falls back to RECOMPUTE on a miss. This gives
bit-faithful coefficients and matrices.

### 5.6 Irregular structures

The time loop itself has no linked chains; the field solve (supercells) and linked boundary
conditions are handled in their specs. exb kx shifts move data between `it` indices of fixed
linkage topology (no re-chaining).

---------------------------------------------------------------------------------------------------

## 6. Pitfalls (things that break round-off agreement)

6.1 **FMA / compiler contraction.** gfortran `-O3` on aarch64 contracts `a*b±c` (e.g. the
    `kx_shift` update, AB coefficients, `phinew + phi` is safe). XLA GPU also contracts. For a
    validation build use `-ffp-contract=off` on the Fortran side; on the JAX side compute the
    time-control scalars on the host in Python floats.

6.2 **Decision flips.** Integer/boolean decisions computed from floats — CFL reset
    (`code_dt > code_dt_cfl`), `jump = NINT(kx_shift/dkx)`, `user_time + user_dt > max_sim_time`,
    omega convergence — can flip at a different step under 1-ulp differences, after which
    trajectories diverge. Provide a **decision replay** mode: the Fortran reference run (or an
    instrumented export program) writes per-step `(code_dt, code_dt_cfl, reset events, jump[ik])`;
    JAX can assert or force identical decisions.

6.3 **NINT semantics.** Fortran `NINT` rounds half away from zero; `jnp.round` rounds half to
    even. Use the exact replacement in §5.3.

6.4 **First-step factor 2 in exb.** `code_dt_old = 0` initially (gs2_time.f90:31) so the first
    applied shear uses `gdt = 0.5*code_dt` — also after a restart (code_dt_old is not saved).

6.5 **exb applied once per step.** `exb_shear::istep_last` (dist_fn.fpp:4217,4270,4286)
    prevents a second shift when the step is redone after an immediate reset; it is only updated
    when the start conditions pass. The restored `g` after a reset already contains the shift.

6.6 **exb all-ky zeroing bug** in the theta0 branch, j > 0: `phi(iright+1:ntgrid,:iib,:) = 0.0`
    (dist_fn.fpp:4370, also 4387 apar, 4404 bpar) zeroes those entries for every ky; process
    ik in ascending order. The kx branch skips ik = 1 (`do ik = naky, 2, -1`, 4492).

6.7 **NL term computed once per step and never recomputed on redo.** After an immediate
    reset, `gexp_1` contains the NL term evaluated with the pre-reset fields; the history was
    shifted once. The redone step uses the new dt and new AB coefficients
    (`this_dt` = new dt, `prev1/prev2` unchanged).

6.8 **Fields after a dt change are recomputed from g** (force_maxwell_reinit = .true., default,
    via `get_fields_direct_from_dfn`), not taken from the field solve — they differ from the
    pre-reset phinew at round-off level. JAX must use the same routine (field spec).

6.9 **AB bootstrap.** Order is chosen by comparing `code_dt_prev2`, then `code_dt_prev1`, with
    the sentinel −1.0: step 1 Euler, step 2 AB2, step 3+ AB3 (unless check_restart restored the
    history). On a dt change prev1/prev2 keep the old dts (variable-step coefficients).

6.10 **Stale hyper arrays.** `D_res` and simple3D `hypervisc_filter` are not rebuilt on dt
    change (hyper.f90:370-371). Do not "fix" this.

6.11 **kwork_filter during rebuilds.** `implicit`: the initial response matrix is built with
    `kwork_filter` all `.false.` (set to (1,1)=.true. only afterwards, fields.f90:376), but a
    rebuild after a dt change runs with `kwork_filter(1,1) = .true.` already set, so the
    (it=1, ik=1) mode is skipped in source/solve during that build and its matrix block differs.
    `local`: `fm_populate` sets its own filter and resets it to `.false.` at the end
    (fields_local.f90:2006-2076), then (1,1) is set again — consistent between builds.
    Replicate per field_option.

6.12 **`scale` quirk on file-based reset.** With `in_memory = .false.` (default) a reset re-reads
    the just-written restart file through `gs2_restore(g, scale, …)` which multiplies g and the
    fields by the init_g `scale` (gs2_save.fpp:865-876). `scale` is set to 1 only after a
    restart-type ginit (init_g.f90:1172-1209), so for e.g. `ginit_option='noise'` with
    `scale /= 1` the first reset rescales the solution. GATE that combination (or replicate with a
    one-shot flag).

6.13 **delt_minimum test before the direction is known.** `code_dt/delt_adj <= code_dt_min`
    aborts even when the reset would increase dt (gs2_reinit.f90:145).

6.14 **Diagnostics modify the state.** Several diagnostics do in-place
    `g_adjust(from_g) … g_adjust(to_g)` round trips on `gnew` and `g` (e.g.
    diagnostics_fluxes.fpp:262-280, diagnostics_moments.fpp:223/329, gs2_heating) at `nwrite`
    steps, including `run_diagnostics(0)` during `initialize_diagnostics`. `(g+a)-a ≠ g` in
    floating point. (a) Export the t = 0 state **after** `initialize_diagnostics`; (b) either
    run the reference with those diagnostics disabled, or have the JAX diagnostics hook replay
    the same round trips at the same steps.

6.15 **Restart write timing.** `gs2_save_for_restart` is called before `update_time`
    (gs2_main.fpp:759-771), so a restart's `t0` is the start-of-step time while g is the
    end-of-step state. Irrelevant for JAX (export is taken after init), but restart-based
    comparisons are off by one step in time.

6.16 **MPI reduction order.** Field solve `sum_allreduce` (fields_implicit.f90:300) and velocity
    integrals depend on the decomposition. Produce reference data with 1 MPI rank (no OpenMP in
    this build).

6.17 **Sequential sums.** Fortran `sum` intrinsics are sequential; `jnp.sum` is a tree. Expect
    ~1e-16 relative differences; for `klein_sum` (compensated) agreement is ~1 ulp regardless.

6.18 **Collisions twice per step.** With `split_collisions = .false.`, collisions act on the
    predictor and on the corrector, and inside every response-matrix impulse.

6.19 **`reset` flag disabled during rebuild.** `reset_time_step` sets `reset = .false.` while
    rebuilding (gs2_reinit.f90:157-158) so `timeadv(…,0)` does not early-return; restores it at
    the end (223).

6.20 **`gnew` after an early return.** When the predictor returns because of `reset`, `gnew` is
    still equal to `g` (post-exb). That is the state saved and restored by the reset.

6.21 **`user_dt` vs `code_dt`.** Only `save_dt` updates `user_dt`; the abort branches of
    `reset_time_step` can leave them different (only matters for the exit path).

6.22 **Linear runs and dt0.** A linear run with `dt0 > delt` increases dt at the first step
    (`code_dt_cfl = 1e8`); support it through the general reset path.

6.23 **`check_time_step` aliasing** (dummy `reset` == module `reset`) is harmless (§1.14).

6.24 **`theta0` single-mode abort.** With `ntheta0 = 1` and `g_exb /= 0`, `dtheta0 = 2π` and
    the run aborts when `|theta0_shift|` reaches ~π (dist_fn.fpp:4347-4351); with box and
    `ntheta0 = 1`, `dkx` is never set (4251-4257). GATE `g_exb /= 0` with `ntheta0 = 1`.

6.25 **`advance_local` vs `advance_implicit`** differ only in which fields are copied/added
    (identical results for zero arrays) and in the exb guard (equivalent).

---------------------------------------------------------------------------------------------------

## 7. The t = 0 export manifest (state needed to start the JAX run identically)

Take the snapshot after `initialize_gs2`, `initialize_equations` and `initialize_diagnostics`
(i.e. exactly where `run_gs2` calls `evolve_equations`, gs2_main.fpp:1640-1650).

State (S):
`gnew` (and `g`, must be equal), `phinew, aparnew, bparnew` (and `phi, apar, bpar`, equal),
`gexp_1, gexp_2, gexp_3` (if allocated), `kx_shift` or `theta0_shift` (if allocated),
`vnmult(2)`, `kwork_filter`, `ieqzip`,
`code_dt, user_dt, code_dt_prev1, code_dt_prev2, code_dt_old, code_time, user_time,
code_dt_cfl, user_dt_cfl, code_dt_min, code_dt_max`,
`reset` (.false.), `exb_first` (.true. unless an exb call happened), the counters of §2.5
(`timeadv istep_last` = 0, `nonlinear istep_last` = 0, `exb istep_last` = 0,
`reset_time_step istep_last` = −1, `nconsec` = 0), `initializing` (.false.), and the init_g
`scale`/`ginitopt_switch` (for §6.12).

Constants (E): §2.1-2.3, §2.6, plus all dt-dependent coefficients of §2.7 at the initial dt
(and optionally the ladder of §5.5), the response matrix, and `istep` semantics (JAX starts at
`istep = 1`).

---------------------------------------------------------------------------------------------------

## 8. Open questions

1. Validation configuration: will reference runs keep the default diagnostics (which perturb
   g/gnew via g_adjust round trips every `nwrite` steps), or switch them off? (§6.14)
2. RECOMPUTE in JAX vs exported dt ladder for coefficients and response matrices after a dt
   change — which is the contract? (§5.5) RECOMPUTE of the inverse matrix will differ from
   Fortran by ~cond·eps.
3. Does `phi(:,1,1)` (box, kx=ky=0) end up non-zero, and does the different response block after
   an `implicit` rebuild (§6.11) visibly change it? Needs a numerical check.
4. XLA GPU FMA contraction policy for float64 and whether it can be disabled for validation.
5. Should the `scale` quirk (§6.12) be replicated or gated? (proposed: gate)
6. Antenna: permanently GATE, or export a pre-drawn random stream from Fortran? (§1.16)
7. Is the all-ky zeroing in the theta0 exb branch (§6.6) known upstream? (replicate regardless)
8. Diagnostics-driven exits (omega convergence for linear runs) — owned by the diagnostics spec;
   the time loop only provides a hook. Confirm.
9. For `field_option='local'`, which inversion algorithm will JAX use when recomputing (the
   local code inverts per supercell block; implicit uses `init_inverse_matrix`)?
