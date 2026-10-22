# 05 — The nonlinear term (E×B Poisson bracket), GS2 FFTs, CFL control, split-NL

Status: contract for the JAX port. READ-ONLY analysis of GS2 8.2.1 (branch `gpu-port`, commit 4d8c94bcf).
All `file:line` citations are relative to `/Users/nico/Workspace/gs2`. `src/` unless otherwise noted.

Scope of this document:

* `src/nonlinear_terms.fpp` (all of it): `add_explicit_terms`, `add_explicit`, `add_nl`, `load_k{x,y}_{phi,apar,bpar}`,
  `calculate_bracket`, `nl_center`, CFL limit (`get_max_vel`, `set_cfl_time_step_and_check_if_too_large`,
  non-blocking variant), `estimate_error`, `zip`, `nonlinear_mode`.
* `src/gs2_transforms.fpp`: `transform2_5d`, `inverse2_5d`, `transform2_5d_accel`, `inverse2_5d_accel`, the
  g→xxf→yxf redistributions, dealiasing, FFTW plan creation (`externals/utils/fft_work.fpp`).
* `src/kt_grids*.f90`: `nx, ny, naky, ntheta0, akx, aky, ikx, theta0`.
* `src/split_nonlinear_terms.fpp`, `src/rk_schemes.f90`: `split_nonlinear`, `strang_split`, `advance_nonadiabatic_dfn`.
* The pieces of `gs2_time.f90`, `gs2_reinit.f90`, `gs2_main.fpp`, `fields_implicit.f90`, `dist_fn.fpp`,
  `dist_fn_arrays.fpp`, `le_grids.f90` that this subsystem reads or that consume its output.

Build facts that matter for round-off (checked on this machine):

* `bin/gs2` links Homebrew FFTW 3.3.x (`/opt/homebrew/opt/fftw/lib/libfftw3.3.dylib`), OpenMPI, no OpenMP
  (`Makefile:66 USE_OPENMP ?=`), compiled with `gfortran 15.2 -O3 -fdefault-real-8 -fdefault-double-8`
  (`Makefiles/Makefile.gnu-gfortran:18,24`) on **arm64**.
* Consequently every literal (`0.5`, `2.0`, `1.e8`, `3.`) is float64 and `epsilon(0.0)` = `2.220446049250313e-16`
  (= `np.finfo(np.float64).eps`).
* gfortran `-O3` on arm64 **does contract `a - b*c` into FMA** (verified: `bracket = bracket - ba*gb` compiles to
  `fmls`/`fmsub`; see §8 P1).

---------------------------------------------------------------------------------------------------------------------

## 0. Executive summary (what the JAX code must do each step)

For a nonlinear (`nonlinear_mode='on'`, `split_nonlinear=.false.`) run, once per time step (first `timeadv` call of
step `istep`, only if `istep /= istep_last`):

1. Shuffle history: `gexp_3 ← gexp_2`, `gexp_2 ← gexp_1`.
2. `gexp_1 ← NL(g, φ, A∥, B∥)` where `g, φ, …` are the values **at the start of the step** (after `exb_shear`):
   * χ-derivatives `i kx χ`, `i ky χ` with `χ = J0 φ fφ + 2 vperp2 J1 (T/Z) B∥ fB − J0 (√(T/m)) v∥ A∥ fA` (GS2 sign/ordering
     exactly as §4.3),
   * `h_w = g + (Z/T) nonmaxw_corr fφ φ J0 + 2 fB vperp2 J1 B∥` (g_wesson, `g_adjust(from_g_gs2)`),
   * real-space bracket `B = (∂xχ)(∂y h_w) − (∂yχ)(∂x h_w)` on an `nx × ny` grid (2/3-rule zero padding),
   * back to spectral, scale `(1/nx)(1/ny)` (×2 for ky>0), times `kxfac`,
   * theta cell-centring with bakdif upwinding (`nl_center`), giving **2×NL at cell centres** (×1 at ttp points).
3. CFL: `max_vel_x = max|∂yχ|·(max akx·0.5/cfl)`, `max_vel_y = max|∂xχ|·(max aky·0.5/cfl)`;
   `dt_cfl = 1/(max_vel_x + max_vel_y)` (if `use_2d_cfl`); `reset = code_dt > dt_cfl`.
4. If `reset .and. immediate_reset`: abort this step, change dt (`gs2_reinit::reset_time_step`), redo the step
   **without recomputing the NL term** (gexp arrays retained).

`gexp_1..3` are then consumed by `dist_fn::get_source_term` as `+0.5*code_dt*(c1 gexp_1 + c2 gexp_2 + c3 gexp_3)` with
variable-step AB3 coefficients (§4.11). The NL term itself has **no dt dependence**; nothing in this subsystem needs
recomputation when dt changes except the AB coefficients (computed every step anyway) and the CFL comparisons.

---------------------------------------------------------------------------------------------------------------------

## 1. Index conventions and notation

Canonical export convention (fixed by the project, restated):

| Object | Fortran order | NumPy (C order) shape |
|---|---|---|
| dist-fn-like `g(ig,isgn,iglo)` | `g(ig=-ntgrid:ntgrid, isgn=1:2, it=1:ntheta0, ik=1:naky, il=1:nlambda, ie=1:negrid, is=1:nspec)` | `(nspec, negrid, nlambda, naky, ntheta0, 2, nth)` with `nth = 2*ntgrid+1` |
| `aj0, aj1, vperp2 (ig,iglo)` | `(ig, it, ik, il, ie, is)` | `(nspec, negrid, nlambda, naky, ntheta0, nth)` |
| fields `phi(ig,it,ik)` | `(ig, it, ik)` | `(naky, ntheta0, nth)` |
| `forbid(ig,il)`, `is_ttp(ig,il)` | `(-ntgrid:ntgrid, 1:nlambda)` | `(nlambda, nth)` |
| `jend(ig)` | `(-ntgrid:ntgrid)` | `(nth,)` |
| `akx(it)`, `ikx(it)` | `(1:ntheta0)` | `(ntheta0,)` |
| `aky(ik)` | `(1:naky)` | `(naky,)` |
| `nonmaxw_corr(ie,is)` | `(1:negrid, 1:nspec)` | `(nspec, negrid)` |

In this document, 0-based NumPy indices are written `ig0 = ig+ntgrid`, `it0 = it-1`, `ik0 = ik-1`, `isgn0 = isgn-1`.
"Line" = a fixed `(is, ie, il, isgn, ig)`; the FFTs act on the `(it, ik)` plane of one line.

`m ≡ (ntheta0-1)/2` (ntheta0 is forced odd in box mode, `kt_grids_box.f90:203-205`).

**kx index → FFT position** (identical in all transform paths):

```
J(it) = it-1                 if it <= (ntheta0+1)/2        (kx = 0, dkx, …, m dkx)
J(it) = it-1 + nx - ntheta0  if it >  (ntheta0+1)/2        (kx = -m dkx, …, -dkx)
```
(`gs2_transforms.fpp:485-489` non-accel; `gs2_layouts.fpp:4700-4710` accel dealias mask; `gs2_transforms.fpp:1593-1603`
for the 3d variant). Equivalently `J(it) = ikx(it) mod nx` with `ikx` from `kt_grids_box.f90:382-390`.

**ky index → FFT position**: `jy = ik-1`, `ik = 1..naky`; positions `naky..ny/2` are zero padding. Only `ky ≥ 0` is
stored; negative ky follows from reality (c2r/r2c transforms).

---------------------------------------------------------------------------------------------------------------------

## 2. Call graph / control flow

### 2.1 Initialisation (once; also re-entered partially on dt reset)

```
gs2_init level 'nonlinear_terms' (24)  src/gs2_init_subroutines.inc:290-300
 └ init_nonlinear_terms                 nonlinear_terms.fpp:287-361
    ├ init_gs2_layouts, init_theta_grid, init_kt_grids, init_le_grids, init_species, init_dist_fn_layouts (308-318)
    ├ read_parameters                   nonlinear_terms.fpp:364-404   (nonlin = nonlinear_mode=='on', 403)
    ├ max_vel_components = 0            :323
    └ if nonlin:
       ├ init_transforms(…, nx, ny, accelerated)  gs2_transforms.fpp:82-151
       │   ├ abort if nx==0 or ny==0                 :106-108   (→ only grid_option='box' can run NL)
       │   ├ measure_plan = fft_measure_plan; load_wisdom if fft_use_wisdom  :113-114
       │   ├ accel = (pe_layout char=='v') .and. mod(negrid*nlambda*nspec,nproc)==0   :119
       │   │      ('v' ⇔ layout ∈ {'yxels','yxles'}, gs2_layouts.fpp:4713-4735)
       │   ├ accel:  init_accel_transform_layouts, init_y/x_transform_layouts, init_accel_fft (:124-134, :156-283)
       │   └ else:   init_y_redist_local → init_x_redist_local (g2x, x2y redistribution maps) (:141, :553-569, :807-827)
       │             init_xy_fft (:286-314): xf_fft c2c BACKWARD n=nx, xb_fft c2c FORWARD n=nx,
       │                                     yf_fft c2r n=ny,           yb_fft r2c n=ny
       ├ allocate ba,gb,bracket (ny, yxf_lo)  or aba,agb,abracket (nth,2,accelx_lo), zeroed  :330-355
       ├ cfly = maxval(aky) * 0.5 / cfl     :357
       └ cflx = maxval(akx) * 0.5 / cfl     :358
fields::init_fields                      fields.f90:353-377
 ├ fields_init_response  (response matrices; calls timeadv with istep=0 → add_explicit_terms is a no-op for nonlin)
 ├ nl_finish_init → nonlinear_terms::finish_init: initializing = .false.   fields.f90:373, nonlinear_terms.fpp:1262-1265
 └ if box: kwork_filter(1,1) = .true.    fields.f90:376
```

`initializing` is never set back to `.true.` (`nonlinear_terms::reset_init`, :1255, has no caller); the
`initializing` branch of `add_explicit` (:466-476) is therefore dead during time stepping.

### 2.2 One time step (implicit fields, the default)

```
gs2_main::evolve_equations             gs2_main.fpp:649-894
 loop istep = istep_end+1, …          :717
  ├ fields::advance(istep)             fields.f90:586-605 → fields_implicit::advance_implicit  fields_implicit.f90:315-371
  │   ├ antenna_amplitudes                              :333
  │   ├ exb_shear(gnew, phinew, …, istep)  (if flow shear; guarded to once per istep, dist_fn.fpp:4270)  :335
  │   ├ g = gnew ; phi = phinew ; apar = aparnew ; bpar = bparnew          :337-338
  │   ├ timeadv(phi,apar,bpar, phinew,…, istep)  [1st call]  dist_fn.fpp:4014
  │   │   ├ (split_nonlinear .and. nonlin) ? advance_nonlinear_term(...)  (GATED, §4.13)   :4046-4058
  │   │   │   : add_explicit_terms(gexp_1,gexp_2,gexp_3, phi,apar,bpar, istep, bkdiff(1))   :4060-4061
  │   │   │        nonlinear_terms.fpp:407-439
  │   │   │        └ if nonlin .and. istep/=0: add_explicit(…, nl=.true.)                    :424-427
  │   │   │             nonlinear_terms.fpp:442-612
  │   │   │             ├ if istep /= istep_last:                                            :495
  │   │   │             │   ├ copy(gexp_2→gexp_3); copy(gexp_1→gexp_2)                      :505-506
  │   │   │             │   ├ max_vel_components = epsilon(0.0)                              :518
  │   │   │             │   ├ add_nl(g, gexp_1, phi, apar, bpar)                             :520, :788-929
  │   │   │             │   ├ nl_center(gexp_1, bd)                                          :523, :677-760
  │   │   │             │   ├ CFL / order-based error / reductions                           :532-583
  │   │   │             │   │   └ set_cfl_time_step_and_check_if_too_large                   :1229-1252
  │   │   │             │   │        ├ save_dt_cfl(dt_cfl)            gs2_time.f90:177-183
  │   │   │             │   │        └ check_time_step_too_large(reset) gs2_time.f90:154-158
  │   │   │             │   └ if reset: return                                               :585-587
  │   │   │             ├ if zip: zero g and gexp_1 for ik/=1                                 :599-609
  │   │   │             └ istep_last = istep                                                 :611
  │   │   │        (if .not. nonlin: save_dt_cfl(1e8))                                       :428-430
  │   │   ├ if (reset .and. immediate_reset) return                                            :4066
  │   │   ├ invert_rhs → get_source_term[_opt] consumes gexp_1..3 with AB coeffs  (§4.11)      :4068
  │   │   ├ hyper_diff, vspace_derivatives (collisions)
  │   ├ if (reset) return                                                     fields_implicit.f90:343
  │   ├ getfield; phinew = phinew + phi …                                     :352-356
  │   └ timeadv(…, istep, diagnostics) [2nd call] → add_explicit: istep == istep_last ⇒ NL NOT recomputed;
  │                                              only the zip block runs.   :365
  ├ if reset (only possible with immediate_reset):                         gs2_main.fpp:739-752
  │    prepare_initial_values_overrides; set_initval_overrides_to_current_vals (saves gnew, phinew, gexp_*)
  │    reset_time_step(state%init, istep, exit)          gs2_reinit.f90:104-225
  │    reset = .false.; if .not.exit: advance(istep)     (redo the step; add_explicit skips NL: istep==istep_last)
  ├ gs2_save_for_restart (every nsave)                                       :767
  ├ update_time: code_dt_prev2=code_dt_prev1; code_dt_prev1=code_dt; code_dt_old=code_dt; time += dt  gs2_time.f90:161-174
  ├ diagnostics
  └ check_time_step(reset, exit)                                             gs2_main.fpp:800, gs2_reinit.f90:228-254
       ├ nb_check_time_step_too_large (only if .not.immediate_reset: wait on MPI_Iallreduce, set dt_cfl, reset) nonlinear_terms.fpp:1191-1197
       ├ reset = code_dt > code_dt_cfl            (too large)
       ├ reset = reset .or. code_dt < min(code_dt_max, code_dt_cfl/delt_adj/delt_cushion)   (too small → increase)
       └ if reset: same override dance + reset_time_step (dt changed for the NEXT step; current step is kept)  :807-816
```

`fields_local::advance_local` (`fields_local.f90:4244`) and `fields_gf_local::advance_gf_local`
(`fields_gf_local.f90:4239`) have the same structure (exb_shear → copy gnew→g, fields → timeadv ×2); the NL term is
identical for all `field_option`s.

### 2.3 dt reset path (what is and is not re-initialised)

* `prepare_initial_values_overrides` drops to level `override_initial_values-1` = 32 (`gs2_main.fpp:1145-1163`),
  `reset_time_step` drops to `override_timestep` = 29 (`gs2_reinit.f90:166`) and comes back to `full` (35)
  (`gs2_reinit.f90:216`). Level list: `src/gs2_init_level_list.inc`.
* Levels re-run: collisions (30), dist_fn_level_3 (31), fields_level_2 (32, response matrices, `init_fields`),
  set_initial_values (34).
* **Not** re-run: nonlinear_terms (24) — `istep_last`, `cflx/cfly`, FFT plans, `accelerated` survive;
  dist_fn_arrays (26) — `gexp_1..3` survive in memory (and are also copied/restored through `initval_ov` when
  `in_memory=.true.`, `gs2_init.fpp:332-336`; with the default `in_memory=.false.` g is round-tripped through the
  restart file, exact for float64, `gs2_main.fpp:1209-1215`).
* `set_initial_field_and_dist_fn_values` (`gs2_init.fpp:292-375`): `g = gnew = saved gnew`; with the default
  `force_maxwell_reinit=.true.` (`fields.f90:118`) the fields are **recomputed** from g (`set_init_fields`,
  `gs2_init.fpp:347-348`) rather than restored. (Fields spec; see §8 P10.)
* gexp_1 computed in the aborted first call is the one used in the redone step.

---------------------------------------------------------------------------------------------------------------------

## 3. Data inventory (everything the NL subsystem reads)

Legend: **E** = EXPORT from Fortran, **R** = RECOMPUTE in JAX (formula given), **B** = export AND recompute/assert
equal. "dt?" = depends on `code_dt` or other runtime-changing quantity.

### 3.1 State (time-evolving) inputs

| Fortran name | Module | Declared shape | Set by | dt? | Decision |
|---|---|---|---|---|---|
| `g` | dist_fn_arrays | `(-ntgrid:ntgrid,2,g_lo%llim_proc:ulim_alloc)` complex | `advance_implicit` copy from `gnew` (`fields_implicit.f90:337`) | evolves | E (initial state), then JAX-owned |
| `phi, apar, bpar` | fields_arrays | `(-ntgrid:ntgrid,ntheta0,naky)` complex | copy from `phinew…` (`fields_implicit.f90:338`) | evolves | E (initial), JAX-owned |
| `gexp_1, gexp_2, gexp_3` | dist_fn_arrays | `(-ntgrid:ntgrid,2,g_lo)` complex; allocated only if `nonlin .and. .not.split_nonlinear` (`dist_fn.fpp:3970-3980`) | `add_explicit` (`nonlinear_terms.fpp:505-523`); zero-initialised; restored from restart if `include_explicit_source_in_restart` (`init_g.f90:240`, `gs2_save.fpp:301-340,728-800`) | no (values); consumed with dt-dependent coefficients | E (initial: zeros or restart values), JAX-owned thereafter |
| `istep_last` | nonlinear_terms (module, `:32`) | integer, init 0 | `add_explicit` :497, :611 | — | R (JAX state, init 0; never reset on dt change) |
| `code_dt` | gs2_time | real | `init_delt`, `reset_time_step` | yes | E (initial) |
| `code_dt_prev1, code_dt_prev2` | gs2_time | real, init `dt_not_set=-1.0` | `init_delt` (from restart if `delt_option='check_restart'`), `update_time` | yes | E (initial), JAX-owned |
| `code_dt_cfl` | gs2_time | real, init `1e8` (`gs2_time.f90:24`) | `save_dt_cfl` | derived | R (per step) |
| `reset` | run_parameters | logical | `check_time_step_too_large`, `check_time_step` | — | R |
| `max_vel_components(3)` | nonlinear_terms `:62` | real(3) | `add_explicit`/`add_nl` | — | R (transient) |

### 3.2 Static coefficient arrays (no dt dependence)

| Fortran name | Module | Shape (Fortran) | Computed in | Notes | Decision |
|---|---|---|---|---|---|
| `aj0` | dist_fn_arrays | `(-ntgrid:ntgrid, g_lo)` real | `init_bessel` `dist_fn.fpp:1756-1823` (`aj0 = j0(arg)`, `arg = bess_fac*smz*sqrt(energy*al*kperp2/bmag)`; **0 where forbid**, :1784-1791) | uses `spfunc::j0` (local implementation, `USE_LOCAL_SPFUNC=on`, `externals/utils/spfunc.fpp`) | E |
| `aj1` | dist_fn_arrays | same | same (`aj1 = j1(arg)`; GS2's local `j1` is **J1(z)/z**, `externals/utils/spfunc.fpp:369-371`) | 0 where forbid | E |
| `vperp2` | dist_fn_arrays | `(-ntgrid:ntgrid, g_lo)` real | `init_vperp2` `dist_fn.fpp:1560-1589`: `bmag*al1*e1` where not forbid, else 0 | | E (could R: `bmag(ig)*al(il)*energy(ie,is)`, same rounding order `(bmag*al1)*e1`) |
| `vpa` | dist_fn_arrays | `(-ntgrid:ntgrid,2,g_lo)` real | `init_vpar` `dist_fn.fpp:1646-1652`: `vpa(:,1)=sqrt(e1*max(0,1-al1*bmag))*vparknob; vpa(:,2)=-vpa(:,1)`; 0 where forbid | | E |
| `spec(is)%stm` | species | real | `species.f90:1009` `sqrt(temp/mass)` | | E |
| `spec(is)%tz` | species | real | `species.f90:1011` `temp/z` | | E |
| `spec(is)%zt` | species | real | `species.f90:1012` `z/temp` | used in g_adjust | E |
| `nonmaxw_corr(ie,is)` | species | `(negrid,nspec)` | `species.f90:571-871` (1.0 for Maxwellian) | used in g_adjust | E |
| `fphi, fapar, fbpar` | run_parameters | real | `run_parameters.f90:313-315`; `beta==0 ⇒ fapar=fbpar=0` (:340-349) | defaults 1, 0, −1 | E |
| `has_phi, has_apar, has_bpar` | run_parameters | logical | `run_parameters.f90:351-353`: `f* > epsilon(0.0)` | | B |
| `kxfac` | theta_grid | real scalar | geometry (`theta_grid.f90:592-619`, `kxfac = abs(kxfac)` :619) | | E |
| `akx` | kt_grids | `(ntheta0)` | `box_get_grids` `kt_grids_box.f90:326-404` | `akx(i)=real(i-1)*dkx` / `real(i-ntheta0-1)*dkx` | B |
| `aky` | kt_grids | `(naky)` | same, `aky(i)=real(i-1)*dky`, `dky=2π/ly` | | B |
| `ikx` | kt_grids | `(ntheta0)` int | same | | B |
| `nx, ny, naky, ntheta0` | kt_grids | int | `init_kt_grids_box` `kt_grids_box.f90:175-274` | defaults: `naky=(ny-1)/3+1`, `ntheta0=2*((nx-1)/3)+1`, `ny=(naky-1)*3+1`, `nx=((ntheta0-1)/2)*3+1` | E |
| `ntgrid, nlambda, negrid, nspec, ng2` | theta_grid/le_grids/species | int | | `grid_has_trapped_particles() = nlambda > ng2` (`le_grids.f90:443-446`) | E |
| `forbid` | le_grids | `(-ntgrid:ntgrid,nlambda)` logical | `calculate_forbidden_region` `le_grids.f90:2716-2739`: `1-al(il)*bmag(ig) < -bouncefuzz` | | E |
| `jend` | le_grids | `(-ntgrid:ntgrid)` int | `calculate_jend` `le_grids.f90:2871-2904`: count of `il` with `1-al*bmag > -bouncefuzz`; **0 for all ig if no trapped particles** | | E |
| `is_ttp` | le_grids | `(-ntgrid:ntgrid,nlambda)` logical | `calculate_ittp` `le_grids.f90:2907-2945` | all false if untrapped | E |
| `bkdiff(1)` (`bd`) | dist_fn | real (nspec) | `dist_fn.fpp:1319-1340` (forced equal across species when `.not.mult_imp`; `mult_imp` forced off for nonlin :1330) | default `bakdif=0.0` | E |
| `cfl` | nonlinear_terms | real | read :382 | default 0.95 | E |
| `cflx, cfly` | nonlinear_terms | real | `:357-358` `maxval(akx)*0.5/cfl`, `maxval(aky)*0.5/cfl` (left-to-right) | | B |
| `accelerated` | nonlinear_terms / gs2_transforms `accel` | logical | `gs2_transforms.fpp:98,119,148` | selects scale constant (§4.6) | E |
| `xb_fft%scale, yb_fft%scale` | gs2_transforms | real | `fft_work.fpp:219-221` `1/real(n)` | non-accel: `1/real(nx)`, `1/real(ny)`; accel: `yb_fft%scale = 1/real(ny*nx)` | R |
| `aidx` (dealias mask, accel) | gs2_transforms | logical(accel_lo) | `gs2_layouts.fpp:4700-4710` | equivalent to padding | R |
| `code_dt_min, code_dt_max` | gs2_time | real | `init_reinit` `gs2_reinit.f90` (`save_dt_min(delt_minimum)`, `save_dt_max(dt0)`, dt0 default = `code_dt_max` = `delt` or restart value) | | E |
| `delt_adj, delt_cushion` | gs2_reinit | real | reinit_knobs (2.0, 1.5) | | E |
| `immediate_reset` | run_parameters | logical | knobs, default `.true.` (`run_parameters.f90:86`) | | E |

Not read by the NL computation itself (but relevant downstream): `kwork_filter` (mode (1,1) skipped in
`get_source_term_on_local_domain`/`invert_rhs`), `ieqzip` (frozen modes in `invert_rhs`), `theta0`/`kperp2`
(only through `aj0/aj1`), `wunits/tunits` (not used here).

### 3.3 Workspaces (never export; optional debug dumps)

`ba, gb, bracket` `(yxf_lo%ny, yxf_lo)` real (`nonlinear_terms.fpp:40`), `aba, agb, abracket` `(nth,2,accelx_lo)`
(`:43`), `xxf (nx, xxf_lo)` complex, `fft (ny/2+1, yxf_lo)` complex (`gs2_transforms.fpp:72,295-296`), `ag`
`(-ntgrid:ntgrid,2,accel_lo)` complex, local `g2` in `add_nl` (`:804`).

Debug-dump layout (for verification only): non-accel `ba(iy, iyxf)` with
`iyxf = ix + nx*(ig0 + nth*(isgn0 + 2*(il-1 + nlambda*(ie-1 + negrid*(is-1)))))` for all layouts except `'yxels'`
(where `ie` and `il` swap places), `gs2_layouts.fpp:3901-3914`; i.e. NumPy C-order
`(nspec, negrid, nlambda, 2, nth, nx, ny)` with `x` index `ix-1`, `y` index `iy-1`.

---------------------------------------------------------------------------------------------------------------------

## 4. Exact algorithms (order of operations preserved)

Conventions in pseudo-code: complex numbers written `(re, im)`. `×` is float64 multiply, rounding after every
operation unless stated. A Fortran product `(complex) * (real)` is evaluated by gfortran either componentwise or as
a full complex multiply with a zero imaginary part — both give bit-identical finite results (the cross terms are
exact zeros), so it is specified componentwise here. Where an FMA may be emitted by gfortran it is flagged **[FMA]**.

### 4.1 `add_explicit_terms` / `add_explicit` (`nonlinear_terms.fpp:407-612`)

```
def add_explicit_terms(istep, bd):                           # :407
    if nonlin:
        if istep != 0: add_explicit(istep, bd)               # :425-427
    else:
        code_dt_cfl = 1e8                                    # :429-430 (dt_cfl_default_large, :90)

def add_explicit(istep, bd):                                 # :442
    if initializing:                                         # :466-476 (dead during stepping)
        code_dt_cfl = 1e8; max_vel_components[:] = 1.0/1e8; return
    if istep != istep_last:                                  # :495
        istep_last = istep                                   # :497
        gexp_3 = gexp_2.copy()                               # :505  (order matters: 3<-2 first)
        gexp_2 = gexp_1.copy()                               # :506
        max_vel_components[0:3] = EPS                        # :518  EPS = 2.220446049250313e-16
        gexp_1, (mvx, mvy) = add_nl(g, phi, apar, bpar)      # :520  (module g = start-of-step g)
        max_vel_components[0] = mvx; max_vel_components[1] = mvy   # set inside add_nl :913-914
        gexp_1 = nl_center(gexp_1, bd)                       # :523
        if not use_cfl_limit:                                # :532-534
            max_vel_components[0:2] = EPS
        if use_order_based_error and istep > istep_error_start:      # :540-556
            error = estimate_error(gexp_1, gexp_2, gexp_3)
            target_dt = code_dt * sqrt(error_target / error)
            max_vel_components[2] = 1.0 / target_dt
        if use_cfl_limit or use_order_based_error:           # :561-583
            if immediate_reset:
                max_vel_components = MPI_max_allreduce(max_vel_components)   # elementwise max, exact
                set_cfl_time_step_and_check_if_too_large()   # may set reset
            else:
                start MPI_Iallreduce(max) ; finished in gs2_reinit::check_time_step (§4.9.3)
        else:
            set_cfl_time_step_and_check_if_too_large()       # with EPS components
        if reset: return                                      # :585-587
    if zip:                                                   # :599-609  (GATED)
        for every iglo with ik != 1:  g[:, :, iglo] = 0 ; gexp_1[:, :, iglo] = 0
    istep_last = istep                                        # :611
```

Notes:
* The block is executed exactly once per `istep` (first `timeadv` call). It is skipped in the 2nd `timeadv` call and
  in both calls of a redone step after an immediate reset (`:484-494` comment).
* `reset` is `run_parameters::reset`; `check_time_step_too_large` assigns it (it can be set to `.false.` as well).
* The `else zero_array(g1)` branch (:589-591) is unreachable for nonlinear runs (`nl` always present).

### 4.2 `add_nl` (`nonlinear_terms.fpp:788-929`) — step by step

Inputs: `g_in` (= module `g`), `phi, apar, bpar`. Output `g1` (= `gexp_1`), max velocities.
`should_adjust = .true.` (default of optional `adjust`, :813).

```
# (A) g1 = i kx chi                                                   :816-818
g1 = load_kx_phi(phi)        # overwrites (or zeros) g1
g1 = load_kx_bpar(g1, bpar)  # adds   (if include_bpar and has_bpar)
g1 = load_kx_apar(g1, apar)  # subtracts (if include_apar and has_apar)
# (B) to real space                                                    :823-829
ba = transform2(g1)          # NB transform2 HALVES g1 in place for ik>1 (side effect, harmless)
max_vel_y = max(|ba|) * cfly                                          # :825/828, :830
# (C) g_wesson, i ky g_wesson                                          :835-847
g1 = g_in.copy()
g1 = g_adjust(g1, phi, bpar, facphi=+fphi, facbpar=+fbpar)            # from_g_gs2, §4.4
g2 = g1.copy()                                                        # local g2 holds g_wesson
g1 = (g1 * zi) * aky[ik]   per iglo   →  (re,im) = (-im_g*aky, re_g*aky)  (exact transform of zi*)
# (D)                                                                  :850-854
gb = transform2(g1)          # g1 halved in place again (ik>1)
# (E) bracket = ba*gb                                                  :858, :1141-1165
bracket = ba * gb            # elementwise real
# (F) g1 = i ky chi                                                    :863-865
g1 = load_ky_phi(phi); g1 = load_ky_bpar(g1,bpar); g1 = load_ky_apar(g1,apar)
# (G)                                                                  :870-877
ba = transform2(g1)          # g1 halved in place for ik>1  ← g1 now holds 0.5*(i ky chi) for ik>1
max_vel_x = max(|ba|) * cflx
# (H) g2 = i kx g_wesson                                               :885-889
g2 = (g2 * zi) * akx[it]   →  (-im*akx, re*akx)
# (I)                                                                  :892-896
gb = transform2(g2)
# (J) bracket -= ba*gb                                                 :900, :1166-1185
bracket = bracket - ba * gb  # [FMA] in gfortran (fmls/fmsub)
# (K) back to spectral space                                           :906-910
g1 = inverse2(bracket, g1)   # non-accel: forbidden lines are NOT written (keep content from (G)); see 4.5.3
max_vel_components[0] = max_vel_x ; max_vel_components[1] = max_vel_y   # :913-914
# (L)                                                                  :921-923
g1 = g1 * kxfac              # complex*real, componentwise
```

Physical meaning: `NL = kxfac * [ (∂xχ)(∂y h) − (∂yχ)(∂x h) ]` in GS2 Fourier-coefficient convention, with
`∂x ↔ i kx`, `∂y ↔ i ky`.

Variable naming trap: `max_vel_y` is built from `∂xχ` (the E×B velocity in y) and multiplied by `cfly` (uses `aky`);
`max_vel_x` from `∂yχ` and `cflx` (uses `akx`). `max_vel_components(1)=max_vel_x`, `(2)=max_vel_y`.

### 4.3 `load_k{x,y}_{phi,bpar,apar}` (`nonlinear_terms.fpp:932-1134`)

For every `(is, ie, il, ik, it)` (any order — independent) and every `ig`; `k = akx(it)` for `load_kx_*`,
`k = aky(ik)` for `load_ky_*`.

**phi** (`:932-962`, `:965-996`): Fortran `fac = zi*k*aj0(:,iglo)*phi(:,it,ik)*fphi` evaluated left-to-right:
```
if has_phi and include_phi:
    c      = k * aj0[ig]                       # ((zi*k)*aj0) = (0, c)
    t_re   = -(c * phi_im[ig])                 # (0,c)*(pr,pi) = (0*pr - c*pi, 0*pi + c*pr)  (exact zeros; FMA-safe)
    t_im   =   c * phi_re[ig]
    fac    = (t_re * fphi, t_im * fphi)
    g1[ig, isgn=1] = fac ; g1[ig, isgn=2] = fac
else:
    g1[...] = 0                                 # zero_array, :960/:993
```
**bpar** (`:1067-1099`, `:1102-1134`), only `if include_bpar and has_bpar` (else no-op, :1081/:1116):
Fortran `fac = zi*k*aj1*2.0*vperp2*spec(is)%tz*bpar*fbpar`
```
    c   = (((k * aj1[ig]) * 2.0) * vperp2[ig]) * tz[is]
    t   = (-(c * bpar_im[ig]),  c * bpar_re[ig])
    fac = (t_re * fbpar, t_im * fbpar)
    g1[ig, 1] = g1[ig, 1] + fac ; g1[ig, 2] = g1[ig, 2] + fac        # complex add
```
**apar** (`:1001-1031`, `:1034-1064`), only `if include_apar and has_apar` (else no-op, :1015/:1048):
Fortran `fac = zi*k*aj0*spec(is)%stm*apar*fapar`
```
    c   = (k * aj0[ig]) * stm[is]
    t   = (-(c * apar_im[ig]),  c * apar_re[ig])
    fac = (t_re * fapar, t_im * fapar)
    g1[ig, isgn] = g1[ig, isgn] - fac * vpa[ig, isgn]               # componentwise: re - fac_re*vpa  [FMA]
```
Order across loads is phi → bpar → apar (`:816-818`, `:863-865`), i.e. `((phi_term + bpar_term) − apar_term·vpa)`.
At forbidden points `aj0 = aj1 = vperp2 = vpa = 0` exactly ⇒ every load is exactly zero there.

### 4.4 `g_adjust(…, from_g_gs2)` (`dist_fn_arrays.fpp:101-175`)

`g_adjust_direction(from_g_gs2)` → `g_adjust_floats(g, phi, bpar, facphi=fphi, facbpar=fbpar)` (:111).
```
trapped     = nlambda > ng2
with_bpar   = (facbpar != 0) ; with_phi = (facphi != 0)
bpar_factor = 2.0*facbpar if with_bpar else 0.0                       # :139-143
for iglo (ik,it,il,ie,is):
    phi_factor = (zt[is] * nonmaxw_corr[ie,is]) * facphi  if with_phi else 0.0     # :155-159
    for ig in -ntgrid..ntgrid:
        if trapped and il > jend[ig]: continue                        # :165  (jend-based, NOT forbid-based)
        r1    = (bpar_factor * vperp2[ig]) * aj1[ig]
        left  = (r1 * bpar_re, r1 * bpar_im)
        right = ((phi_factor * phi_re) * aj0[ig], (phi_factor * phi_im) * aj0[ig])
        adj   = left + right                                          # [FMA] possible when bpar term present
        g[ig,1] = g[ig,1] + adj ; g[ig,2] = g[ig,2] + adj
```
Note: `phi_factor` is a loop-private variable initialised to 0 only once (:137); with `with_phi` true it is set
every iteration. If `fbpar = 0` the bpar term is an exact `+0` (no rounding change).
Note: if `beta /= 0` and `fbpar=-1` (default) were left, `with_bpar` is true even though `has_bpar` is false; the term
is then `-2*vperp2*aj1*bpar` with `bpar ≡ 0` → exact zero. Harmless.

### 4.5 Non-accelerated transforms (default path)

#### 4.5.1 `transform2_5d(g, yxf)` (`gs2_transforms.fpp:919-955`)

```
for iglo: if ik_idx(iglo) != 1: g[:, :, iglo] = g[:, :, iglo] * 0.5        # IN PLACE, exact (power of 2)  :947-951
transform_x5d(g, xxf)                                                       # :830-855
    xxf[:, :] = 0                                                           # zero_array (padding + skipped points)
    gather g2x:  xxf[J(it)+1, ixxf(ik,ig,isgn,il,ie,is)] = g[ig, isgn, iglo(ik,it,il,ie,is)]   # :318-550
                 for it = 1..ntheta0, ik = 1..naky  — SKIPPED (left 0) if forbid(ig,il) and
                 (destination proc /= iproc or .not. opt_local_copy)       # :786-802
    xxf[:, ixxf] = FFTW c2c BACKWARD length nx, unnormalised, in place      # xf_fft, :853, :308
transform_y5d(xxf, yxf)                                                     # :875-899
    fft[:, :] = 0                                                           # :890
    gather x2y:  fft[ik, iyxf(ix,ig,isgn,il,ie,is)] = xxf[ix, ixxf(ik,ig,isgn,il,ie,is)]
                 for ix = 1..nx, ik = 1..naky; rows naky+1..ny/2+1 stay 0; forbidden lines skipped   # :573-782
    yxf[:, iyxf] = FFTW c2r length ny, unnormalised                         # yf_fft, :897, :312
```
Mathematically, per line, with `Ĝ[ik,it] = G[ik,it]·(1 if ik==1 else 0.5)`:
```
X[jy, ix] = Σ_{it} Ĝ[jy+1, it] · exp(+2πi · J(it) · ix / nx)                    (jy < naky; else 0)
f[iy, ix] = Re X[0, ix] + Σ_{jy=1}^{⌊(ny-1)/2⌋} 2 Re( X[jy, ix] · exp(+2πi jy iy / ny) )  (+ Nyquist term, always 0)
```
The `Re X[0,ix]` (imaginary part of the ky=0 / DC element discarded) is FFTW's c2r behaviour; **verified** on the
installed FFTW (scratch test: 1D and 2D c2r with random non-Hermitian DC input agree with the "Re-part / Hermitian
projection" formula to 5e-15 and 9e-16 relative). Because GS2's ky=0 data need not be exactly Hermitian in kx, the
JAX recipe must reproduce this explicitly (§5).

Coefficient convention (`:929-942`): `F^gs2_k = f_k` for ky=0 and `2 f_k` for ky>0, so that
`f(x,y) = Σ_{kx} F_{kx,0} e^{i kx x} + Σ_{kx, ky>0} 2 Re( (F/2) e^{i(kx x + ky y)} )`.

Real-space grid: `x_ix = (ix-1)·Lx/nx`, `y_iy = (iy-1)·Ly/ny` (only pointwise products and maxima are taken, so the
physical coordinates never enter).

#### 4.5.2 `inverse2_5d(yxf, g)` (`gs2_transforms.fpp:958-985`)

```
inverse_y5d(yxf, xxf)                                                       # :902-916
    fft[:, iyxf] = FFTW r2c length ny, unnormalised (exp(-2πi jy iy/ny)), jy = 0..ny/2    # yb_fft :912
    scatter x2y: xxf[ix, ixxf(ik,…)] = fft[ik, iyxf(ix,…)]  for ik = 1..naky only (truncation), forbidden lines skipped
inverse_x5d(xxf, g)                                                         # :858-872
    xxf[:, ixxf] = FFTW c2c FORWARD length nx, unnormalised, in place      # xb_fft :868
    scatter g2x: g[ig, isgn, iglo(ik,it,…)] = xxf[J(it)+1, ixxf]  for it = 1..ntheta0 (truncation), forbidden skipped
scale = xb_fft%scale * yb_fft%scale        = (1.0/real(nx)) * (1.0/real(ny))   # :968 (two roundings + one)
for iglo: if ik != 1: g = g * (scale*2.0)   else: g = g * scale           # :977-983  (scale*2.0 computed first)
```
Mathematically, per line: `Y[jy, ix] = Σ_iy f[iy,ix] e^{-2πi jy iy/ny}`, `Z[ik, jx] = Σ_ix Y[ik-1, ix] e^{-2πi jx ix/nx}`,
`G_out[ik, it] = Z[ik, J(it)] · s_ik`, `s_1 = scale`, `s_{ik>1} = scale*2.0`.

#### 4.5.3 Forbidden region in the non-accel path

Forbidden lines (`forbid(ig,il)`) are neither gathered nor scattered (unless local and `opt_local_copy`), so:
* real-space `ba/gb` are 0 on those lines;
* in `inverse2`, `g1` at forbidden `(ig, il)` keeps its content from step (G) of §4.2, i.e. `0.5·(i ky χ)` (ik>1) or
  `i ky χ` (ik=1), then is multiplied by `s_ik` and `kxfac`.

Because `aj0 = aj1 = vperp2 = vpa = 0` at every forbidden point (`dist_fn.fpp:1584-1586, 1649-1652, 1784-1788`),
`χ` and its loads are **exactly zero** there, so: (i) the stale values are zero, (ii) the CFL maxima do not depend
on whether forbidden lines are transformed, (iii) the accel path, `opt_local_copy=.true.` and
`nl_forbid_force_zero=.false.` all produce identical numbers at forbidden points (±0). The JAX port therefore needs no
forbidden mask in the FFT path (it may still apply one; it is exact either way).

### 4.6 Accelerated transforms (`layout ∈ {'yxels','yxles'}` and `mod(negrid*nlambda*nspec, nproc)==0`)

`transform2_5d_accel` (`gs2_transforms.fpp:990-1278`, non-SHMEM parts :1033-1065, :1105-1111, :1120-1128):
```
idx = g_lo%llim_proc
for k in accel_lo (ik_pad fastest (ndky=ny/2+1), then it_pad (nx), then il/ie/is per layout):
    if aidx(k):  ag[:, :, k] = 0                                   # dealiased slot: it_pad in (ntheta0/2+1, nx-ntheta0/2] or ik_pad > naky
    else:
        if ik_idx(g_lo, idx) != 1: g[:,:,idx] = 0.5*g[:,:,idx]     # in place, exact
        ag[:, :, k] = g[:, :, idx];  idx += 1
for each (il,ie,is) block: axf = FFTW 2D c2r, logical size (nx, ny) (C order; y is the r2c dim),
                            howmany = 2*nth interleaved (stride = 2*nth, dist = 1)   # init_crfftw_2d fft_work.fpp:512-573
```
FFTW rank-2 c2r = c2c along x then c2r along y (same semantics as §4.5.1, verified numerically).
`inverse2_5d_accel` (`:1283-1352`):
```
ag = FFTW 2D r2c of axf (unnormalised)                             # yb_fft :1310
for non-dealiased k (in order), idx over g_lo:
    g[:,:,idx] = (2.0 * yb_fft%scale) * ag[:,:,k]   if ik != 1      # :1326-1328
    g[:,:,idx] =         yb_fft%scale  * ag[:,:,k]  if ik == 1      # :1333-1335
yb_fft%scale = 1/real(ny*nx)   (fft%init(m*n), fft_work.fpp:481, :220)    ← single rounding, differs from (1/nx)*(1/ny)
```
Differences vs non-accel (all at round-off level only): 2D FFTW plan instead of two 1D plans, scale constant
`1/(nx*ny)` instead of `(1/nx)*(1/ny)`, all lines (incl. forbidden, which are zero) transformed. CFL max:
`gs2_max_abs(aba)` over all lines — identical values by §4.5.3. SHMEM variants (`#ifdef SHMEM`) are GATED.

### 4.7 `calculate_bracket` (`nonlinear_terms.fpp:1141-1187`)
```
first : bracket[iy, iyxf] = ba[iy, iyxf] * gb[iy, iyxf]
second: bracket[iy, iyxf] = bracket[iy, iyxf] - ba[iy, iyxf] * gb[iy, iyxf]     # gfortran -O3 arm64: FMA (fmls)
```
Accel: same with `abracket/aba/agb`.

### 4.8 `nl_center(gtmp, bd)` (`nonlinear_terms.fpp:677-760`)

```
trapped = nlambda > ng2                                          # :699
for iglo (il = il_idx(iglo)):
    if nl_forbid_force_zero:                                     # default .true. (:83, :116)
        gtmp[ig, 1:2] = 0  where forbid[ig, il]                  # :718-724   (all ig incl. ±ntgrid)
    for ig in -ntgrid .. ntgrid-1 (ascending, in place):         # :726
        if is_ttp[ig, il]: continue                              # :739   keep raw grid value (×1, not ×2!)
        if trapped and (il > jend[ig] or il > jend[ig+1]):       # :741
            gtmp[ig, 1:2] = 0
        else:
            gtmp[ig,1] = (1.+bd)*gtmp[ig+1,1] + (1.-bd)*gtmp[ig,1]      # :753  [FMA possible; exact if bd==0]
            gtmp[ig,2] = (1.-bd)*gtmp[ig+1,2] + (1.+bd)*gtmp[ig,2]      # :754
    # gtmp[ntgrid, :] is never touched by the loop (keeps raw grid value, or 0 if forbidden)
```
In-place semantics: when `ig` is processed, `gtmp[ig+1]` has not been modified yet (ascending loop) and `gtmp[ig]` is
read before being written, so every read is of the **pre-loop** array. Fully vectorisable:
```
A   = where(forbid, 0, G)                         (if nl_forbid_force_zero)
C1  = (1+bd)*A[ig+1, s=1] + (1-bd)*A[ig, s=1] ;  C2 = (1-bd)*A[ig+1, s=2] + (1+bd)*A[ig, s=2]     for ig < ntgrid
out[ig] = A[ig]           if is_ttp[ig,il]
        = 0               elif trapped and (il > jend[ig] or il > jend[ig+1])
        = (C1, C2)        else
out[ntgrid] = A[ntgrid]
```
`(1.+bd)` and `(1.-bd)` are float64 scalars computed once; `bd = bkdiff(1)` (passed from `dist_fn.fpp:4061`).
The wfb (il = ng2+1) needs no special treatment (never forbidden). Downstream (`dist_fn.fpp:4987-5010, 5106-5126`)
uses only `ig ≤ ntgrid-1`; `gexp(ntgrid)` matters only for `estimate_error` and restart files.

### 4.9 CFL limit and time-step reset

#### 4.9.1 `get_max_vel` / `set_cfl_time_step_and_check_if_too_large` (`nonlinear_terms.fpp:1214-1252`)
```
def get_max_vel(c):                       # c = max_vel_components after reduction
    max_vel_cfl   = c[0] + c[1]  if use_2d_cfl else max(c[0], c[1])
    max_vel_error = c[2]
    return max(max_vel_cfl, max_vel_error)
def set_cfl_time_step_and_check_if_too_large():
    max_vel = get_max_vel(max_vel_components)
    dt_cfl  = 1e8 if max_vel <= 0.0 else 1.0/max_vel
    code_dt_cfl = dt_cfl ; user_dt_cfl = dt_cfl           # save_dt_cfl, gs2_time.f90:177-183
    reset = (code_dt > code_dt_cfl)                       # gs2_time.f90:157
```
Since `c[2] = EPS` unless the order-based error is active, `max_vel ≥ EPS`; with all fields zero
`dt_cfl = 1/EPS ≈ 4.5e15`; with `use_cfl_limit=.false.` and no error control `dt_cfl = 1/(2 EPS)` (2D) or `1/EPS`.
The `max` reductions (`gs2_max_abs`, `array_utils.fpp:262-316`; MPI max) are exact and order-independent.

#### 4.9.2 Immediate reset (default `immediate_reset=.true.`)
`add_explicit` returns early with `reset=.true.` → `timeadv` returns (`dist_fn.fpp:4066`) → `advance_implicit`
returns (`fields_implicit.f90:343`) → `gs2_main.fpp:739-752`: `reset_time_step`, then `advance(istep)` again.

#### 4.9.3 Deferred reset (`immediate_reset=.false.`)
`nb_max_allreduce` (`:573`); the step completes with the old dt; at the end of the step
`gs2_reinit::check_time_step` (`gs2_reinit.f90:228-254`) → `nb_check_time_step_too_large` (`nonlinear_terms.fpp:1191-1197`:
wait, then `set_cfl_time_step_and_check_if_too_large`) → then:
```
reset = False                                         # NB: dummy 'reset' aliases run_parameters::reset (same storage)
if code_dt*1.0 > code_dt_cfl: reset = True
if code_dt*1.0 < min(code_dt_max, code_dt_cfl/delt_adj/delt_cushion): reset = True     # ((cfl/adj)/cushion)
```
This end-of-step check runs for every non-exiting step in both modes (it is the only place where dt is *increased*).

#### 4.9.4 `reset_time_step` (`gs2_reinit.f90:104-225`) — needed to reproduce the dt sequence
```
nconsec = nconsec+1 if istep_last_reset + 1 == istep else 0          # save vars, istep_last init -1
if nconsec > 4 and abort_rapid_time_step_change: exit run
if code_dt/delt_adj <= code_dt_min: code_dt = code_dt_min; exit run
attempts = 0
while attempts < 10:
    if   code_dt > code_dt_cfl:                                   code_dt = code_dt/delt_adj          # reduce
    elif code_dt < min(code_dt_max, code_dt_cfl/delt_adj/delt_cushion): code_dt = min(code_dt*delt_adj, code_dt_max)
    else: break
    attempts += 1
if attempts == 10: exit run
user_dt = code_dt (save_dt)
if code_dt <= code_dt_min: code_dt = code_dt_min; exit run
re-initialise dt-dependent levels 29..35 (response matrices etc. — other specs); istep_last_reset = istep
```
`code_dt_prev1/2` are NOT touched by a reset; the next AB3 step uses the variable-step formula with the new `code_dt`.

### 4.10 `estimate_error` (`nonlinear_terms.fpp:620-672`) — only if `use_order_based_error`
```
n_avail = time_history_available()          # 3 if prev2 set, 2 if prev1 set, else 1   (gs2_time.f90:234-244)
if min(n_avail,3) != 3 or min(n_avail,2) != 2:
    return error_target * cfl * cfl          # ((error_target*cfl)*cfl)  → target_dt = code_dt*sqrt(et/(et*cfl*cfl))
c3 = get_adams_bashforth_coefficients(maximum_order=3)    # this_dt = code_dt
c2 = get_adams_bashforth_coefficients(maximum_order=2)
s3 = (c3[0]*g1 + c3[1]*g2) + c3[2]*g3       # elementwise complex, left-to-right [FMA]
s2 =  c2[0]*g1 + c2[1]*g2
e1 = max |s3 - s2|   (complex abs = hypot)  over the WHOLE local array (incl. ig=ntgrid, forbidden, ttp)
e2 = max |s2|
(e1, e2) = MPI_max_allreduce
return e1 / e2
```
Called with `g1,g2,g3 = gexp_1,gexp_2,gexp_3` after shuffle and `nl_center`.

### 4.11 Consumer interface: AB coefficients and source (`gs2_time.f90:250-312`, `dist_fn.fpp:5363-5420`)

`include_explicit = istep > 0 .and. nonlin .and. .not. split_nonlinear` (`dist_fn.fpp:5381-5386`);
coefficients recomputed on every `invert_rhs` call with `this_dt = code_dt`, `p1 = code_dt_prev1`, `p2 = code_dt_prev2`:
```
order 3 (p2 set):
  c1 = ( (((dt + p1)**2)/p1) * ((dt + p1)/3. + p2/2.)  -  p1*(p1/3. + p2/2.) ) / (dt*(p1 + p2))
  c2 = -( (dt * (dt/3. + (p1 + p2)/2.)) / (p1*p2) )
  c3 =    (dt * (dt/3. + p1/2.)) / (p2*(p1 + p2))
order 2 (only p1 set):  c1 = (0.5*dt + p1)/p1 ;  c2 = -((0.5*dt)/p1)
order 1:                c1 = 1.0
```
(`x**2` with integer exponent is `x*x` in gfortran; divisions by `3.`/`2.` are true divisions.)
Use in `get_source_term` (`set_source`, `dist_fn.fpp:5106-5126`; ttp `:4987-5010`; opt variant `:5268-5290, 5332-5351`):
`source(ig) = source(ig) + (0.5*code_dt) * ((c1*gexp_1(ig) + c2*gexp_2(ig)) + c3*gexp_3(ig))` for `ig = -ntgrid..ntgrid-1`.
Fresh start: steps 1, 2, 3 use orders 1, 2, 3 (prev1/prev2 become set by `update_time`). With
`gexp_2 = gexp_3 = 0` initially the order-2/3 terms multiply zeros on the first step(s).

### 4.12 `zip` (GATED)
`nonlinear_terms.fpp:599-609`: on **every** `add_explicit` call that does not return early (both `timeadv` calls),
zero `g(:,:,iglo)` (the start-of-step g used by `invert_rhs`) and `gexp_1(:,:,iglo)` for all `ik ≠ 1`.

### 4.13 `split_nonlinear` path (GATED; documented for a later phase)

Selected by `split_nonlinear=.true.` in `nonlinear_terms_knobs` (requires `nonlin`). Then `gexp_*` are not
allocated, `include_explicit=.false.`, and in `timeadv` (`dist_fn.fpp:4038-4109`):
```
nl_dt = code_dt*0.5 if strang_split else code_dt
first timeadv call of a step (istep /= timeadv's own save istep_last, istep>0):
    advance_nonlinear_term(g, istep, phi, apar, bpar, nl_dt, fields_func, calculate_current_nl_source_and_cfl_limit)
    fields_func = calculate_potentials_from_nonadiabatic_dfn (dist_fn.fpp:7322) if advance_nonadiabatic_dfn
                  else get_fields_direct_from_dfn (dist_fn.fpp:7590)   — direct (non-response-matrix) field solves
... linear implicit step ...
if strang_split, second timeadv call (istep == istep_last, istep>0), after collisions/parity:
    advance_nonlinear_term(gnew, istep, phinew, aparnew, bparnew, nl_dt, …)
```
NB: dist_fn's `timeadv` has its own `save istep_last` (init −1, `dist_fn.fpp:4031`); response-matrix construction
calls `timeadv` with `istep=0` and resets it.

`calculate_current_nl_source_and_cfl_limit` (`nonlinear_terms.fpp:764-785`): **sets the module variable
`use_cfl_limit = calculate_cfl_limit` (side effect!)**, calls `add_nl(g_in, g1, phi, apar, bpar, need_to_adjust)`
(no `nl_center`), returns `max_vel_local = get_max_vel(max_vel_components)` (with whatever component 3 holds) or
`1e-8` if not computing the limit.

`advance_nonlinear_term` (`split_nonlinear_terms.fpp:205-280`): copies fields to locals; if
`advance_nonadiabatic_dfn`: `g_adjust(g_state, phinew, bparnew, from_g_gs2)` before and
`g_adjust(g_state, phi, bpar, to_g_gs2)` after (note: *original* fields used for the back-conversion); dispatches:

* `split_method='AB3'` (`:284-396`): `code_dt_prev1 = code_dt_prev2 = dt_not_set` (**clobbers gs2_time globals**);
  loop `while internal_time < dt`: source with CFL → `max_allreduce(max_vel)`; `time_step = (1/max_vel)*time_step_safety_factor`;
  abort if `< code_dt_min`; `time_step = min(time_step, dt - internal_time)`; `ab = get_adams_bashforth_coefficients(target_dt=time_step)`;
  `g += (0.5*time_step)*(ab1 gnl_1 [+ ab2 gnl_2 [+ ab3 gnl_3]])` by `counter`; fields from g (`local_only`);
  `prev2 = prev1; prev1 = time_step`; shift `gnl_3←gnl_2←gnl_1`.
* `'beuler'` (`:400-599`): fixed-point backward Euler with `iteration_limit=20`, `convergence_tolerance`,
  error estimate vs CN, `time_step_last` saved across calls, reject/accept with `immediate_reset` and CFL.
* `'RK'`/`'default'` (`:604-814`): embedded RK (`rk_schemes.f90`, default `heun` = 2-stage order 2(1)),
  stage update `g_state = g_cur + Σ_{j<i} (0.5*time_step*coeffs(j,i)) * stage_j`, error
  `errors(1) = max|Σ (high−low)_i stage_i| * half_time`, `errors(2) = max|g_cur|`, `errors(3) = max_vel(stage 1)`,
  `inverse_error = (relative_tolerance*errors(2) + absolute_tolerance)/(errors(1) + EPS)`,
  `new_time_step = min(time_step*inverse_error**(1.0/order), 1/errors(3)) * time_step_safety_factor`,
  reject if `inverse_error < 1 .or. (cfl_limited .and. immediate_reset)`; `time_step_last` saved.
* `'picard'` (`:818-982`): Picard iteration, `npic_limit=15`, step doubling/halving.

Every source evaluation is a full `add_nl` (5 FFT passes); field solves are direct velocity integrals. Data-dependent
trip counts ⇒ `lax.while_loop` with dynamic termination; reductions `max_allreduce` exact.

---------------------------------------------------------------------------------------------------------------------

## 5. JAX FFT recipe (reproduces GS2/FFTW semantics; round-off agreement ≈ FFT-library level)

All arrays float64 / complex128 (`jax.config.update("jax_enable_x64", True)`). Two-stage 1D transforms in the same
order as GS2 (x then y forward; y then x inverse), explicit DC-imaginary drop, explicit GS2 normalisation.

```python
import numpy as np
import jax, jax.numpy as jnp
jax.config.update("jax_enable_x64", True)

def kx_fft_index(ntheta0, nx):
    """J(it) for it0 = 0..ntheta0-1 (0-based), gs2_transforms.fpp:485-489."""
    it0 = np.arange(ntheta0)
    return np.where(it0 <= (ntheta0 - 1) // 2, it0, it0 + nx - ntheta0)

def gs2_transform2(G, nx, ny):
    """G: canonical complex128 (..., naky, ntheta0, 2, nth)  -> real (..., 2, nth, nx, ny).
    == transform2_5d / transform2_5d_accel (values), excluding the in-place halving side effect."""
    naky, ntheta0 = G.shape[-4], G.shape[-3]
    half = np.where(np.arange(naky) == 0, 1.0, 0.5)            # ik==1 unscaled, else *0.5  (exact)
    G = G * half[:, None, None, None]
    G = jnp.moveaxis(G, (-4, -3), (-1, -2))                     # (..., 2, nth, ntheta0, naky)
    jx = kx_fft_index(ntheta0, nx)
    P = jnp.zeros(G.shape[:-2] + (nx, ny // 2 + 1), jnp.complex128)
    P = P.at[..., jx, :naky].set(G)                             # 2/3-rule zero padding in kx and ky
    X = jnp.fft.ifft(P, axis=-2, norm="forward")                # Σ_j P e^{+2πi j x/nx}, NO 1/nx (FFTW_BACKWARD)
    X = X.at[..., 0].set(jax.lax.complex(X[..., 0].real, jnp.zeros_like(X[..., 0].real)))  # FFTW c2r drops Im(DC)
    # (Nyquist column ny//2 for even ny is always zero-padded since naky <= (ny-1)//3+1)
    return jnp.fft.irfft(X, n=ny, axis=-1, norm="forward")      # unnormalised c2r  (FFTW c2r)

def gs2_inverse2(f, naky, ntheta0, nx, ny, accelerated=False):
    """f: real (..., 2, nth, nx, ny) -> canonical complex (..., naky, ntheta0, 2, nth). == inverse2_5d[_accel]."""
    Y = jnp.fft.rfft(f, axis=-1)                                # unnormalised r2c (FFTW r2c, exp(-i))
    Y = Y[..., :naky]                                           # keep ky = 0..naky-1
    Z = jnp.fft.fft(Y, axis=-2)                                 # unnormalised c2c FORWARD along x
    Z = Z[..., kx_fft_index(ntheta0, nx), :]                    # (..., ntheta0, naky)
    if accelerated:
        scale = 1.0 / float(nx * ny)                            # fft_work.fpp:220 with n = ny*nx
    else:
        scale = (1.0 / float(nx)) * (1.0 / float(ny))           # gs2_transforms.fpp:968
    fac = np.where(np.arange(naky) == 0, scale, scale * 2.0)    # (scale*2.0) formed first, :979
    Z = jax.lax.complex(Z.real * fac, Z.imag * fac)             # complex * real, componentwise
    return jnp.moveaxis(Z, (-1, -2), (-4, -3))

def gs2_nl_bracket(chi_x, chi_y, h_y, h_x, nx, ny, kxfac, accelerated=False):
    """chi_x = i kx chi, chi_y = i ky chi, h_y = i ky g_w, h_x = i kx g_w (canonical spectral arrays)."""
    ba1 = gs2_transform2(chi_x, nx, ny); max_abs_dchidx = jnp.max(jnp.abs(ba1))
    gb1 = gs2_transform2(h_y,  nx, ny)
    br  = ba1 * gb1
    ba2 = gs2_transform2(chi_y, nx, ny); max_abs_dchidy = jnp.max(jnp.abs(ba2))
    gb2 = gs2_transform2(h_x,  nx, ny)
    br  = br - ba2 * gb2                                        # Fortran: FMA (see P1)
    naky, ntheta0 = chi_x.shape[-4], chi_x.shape[-3]
    nl = gs2_inverse2(br, naky, ntheta0, nx, ny, accelerated)
    nl = jax.lax.complex(nl.real * kxfac, nl.imag * kxfac)      # :921-923
    return nl, max_abs_dchidy, max_abs_dchidx                   # -> max_vel_x = [1]*cflx, max_vel_y = [2]*cfly
```

Why not `jnp.fft.irfft2/rfft2`: mathematically identical, but (a) the DC-imaginary drop must happen between the x
and y stages to be backend-independent (cuFFT C2R behaviour on non-Hermitian input is undefined; ducc/pocketfft drop
it), and (b) the 1D-per-axis form keeps the operation order of GS2's non-accel path. `irfftn` on JAX moves axes to
the end and calls a single multi-dim IRFFT — acceptable only if Hermitian projection is applied beforehand; the
two-stage form above is preferred.

Mapping table (one line, NumPy 0-based):

| GS2 | JAX |
|---|---|
| `g(ig, isgn, iglo(ik,it,…))`, ik>1 halved | `G[..., ik0, it0, isgn0, ig0] * 0.5` |
| `xxf(J(it)+1, ixxf(ik,…))` | `P[..., isgn0, ig0, J(it), ik0]` |
| `xf_fft` (c2c BACKWARD, nx) | `jnp.fft.ifft(axis=-2, norm="forward")` |
| `fft(ik, iyxf(ix,…))`, `ik ≤ ny/2+1` | `X[..., isgn0, ig0, ix0, ik0]` |
| `yf_fft` (c2r, ny) | `jnp.fft.irfft(n=ny, axis=-1, norm="forward")` after Im(DC)=0 |
| `ba(iy, iyxf(ix,…))` | `f[..., isgn0, ig0, ix0, iy0]` |
| `yb_fft` (r2c) | `jnp.fft.rfft(axis=-1)` (norm "backward" = unnormalised forward) |
| `xb_fft` (c2c FORWARD) | `jnp.fft.fft(axis=-2)` |
| `scale*2.0` / `scale` | `fac[ik0]` |

---------------------------------------------------------------------------------------------------------------------

## 6. Namelist / configuration options

### 6.1 `nonlinear_terms_knobs` (`nonlinear_terms.fpp:93-150`)

| Option | Default | Effect | Decision |
|---|---|---|---|
| `nonlinear_mode` | `'default'` (= none) | `'on'` ⇒ `nonlin=.true.`; `'none'/'off'/'default'` linear (`:369-373, :403`) | SUPPORT |
| `cfl` | 0.95 | `cflx/cfly = max(k)*0.5/cfl` | SUPPORT |
| `use_cfl_limit` | `.true.` | if false, components 1:2 = EPS | SUPPORT |
| `use_2d_cfl` | `.true.` | sum (true) vs max (false) of x/y velocities | SUPPORT |
| `use_order_based_error` | `.false.` | AB2 vs AB3 error control (§4.10) | SUPPORT (cheap) |
| `error_target` | 0.1 | target relative error | SUPPORT |
| `istep_error_start` | 30 | error control only for `istep > istep_error_start` | SUPPORT |
| `include_phi/apar/bpar` | `.true.` | drop field contributions to χ in the NL term only | SUPPORT |
| `nl_forbid_force_zero` | `.true.` | zero gexp at forbidden points before centring | SUPPORT (both values give identical numbers, §4.5.3) |
| `zip` | `.false.` | zero g and gexp_1 for ky≠0 every call (§4.12) | GATE (expert/secondary-mode option, modifies g in place) |
| `split_nonlinear` | `.false.` | IMEX split NL integrators (§4.13) | GATE (experimental; needs direct field solves and adaptive sub-stepping) |

### 6.2 `split_nonlinear_terms_knobs` (`split_nonlinear_terms.fpp:50-112`) — only relevant if `split_nonlinear`
`split_method` ('default'=RK), `rk_method` ('default'=heun), `strang_split` (F), `advance_nonadiabatic_dfn` (F),
`absolute_tolerance` (1e-3), `relative_tolerance` (1e-2), `convergence_tolerance` (0.1), `time_step_safety_factor` (0.9),
`show_statistics` (F). All GATE (with `split_nonlinear`).

### 6.3 Other namelists that change behaviour here

| Option | Namelist | Default | Effect | Decision |
|---|---|---|---|---|
| `immediate_reset` | knobs | `.true.` | reset inside step vs at end of step | SUPPORT |
| `fphi, fapar, fbpar` | knobs | 1.0, 0.0, −1.0 (beta=0 ⇒ fapar=fbpar=0) | field scalings, `has_*` | SUPPORT |
| `beta` | knobs | 0 | zero ⇒ electrostatic | SUPPORT |
| `wstar_units` | knobs | `.false.` | must be false for NL (`check_nonlinear_terms` :169-176) | GATE if `nonlin .and. wstar_units` |
| `delt` | knobs | 0.1 | initial dt, default `code_dt_max` | SUPPORT |
| `delt_adj` | reinit_knobs | 2.0 | dt change factor | SUPPORT (warn if not a power of 2, P9) |
| `delt_cushion` | reinit_knobs | 1.5 | hysteresis for dt increase | SUPPORT |
| `delt_minimum` | reinit_knobs | 1e-5 | `code_dt_min` | SUPPORT |
| `dt0` | reinit_knobs | `code_dt_max` (=delt) | `code_dt_max` | SUPPORT |
| `abort_rapid_time_step_change` | reinit_knobs | `.true.` | abort after >4 consecutive resets | SUPPORT |
| `in_memory` | reinit_knobs | `.false.` | reset state via memory vs restart file (both exact) | SUPPORT |
| `force_maxwell_reinit` | fields_knobs | `.true.` | after reset, fields recomputed from g | SUPPORT (fields spec) |
| `field_option` | fields_knobs | implicit | NL identical for implicit/local/gf_local | SUPPORT |
| `grid_option` | kt_grids_knobs | 'default'(single) | NL needs `'box'`/`'nonlinear'` (nx,ny>0, `gs2_transforms.fpp:106`) | SUPPORT box; GATE `nonlin` with non-box |
| `nx, ny, naky, ntheta0` | kt_grids_box_parameters | 0 (derived) | FFT sizes, padding | SUPPORT (ntheta0 odd) |
| `y0, ly, jtwist, x0, rtwist, n0, rhostar_box` | kt_grids_box_parameters | 2.0, 0, smart, 0, 0, 0, 0 | `akx, aky` (export) | SUPPORT |
| `bakdif` | dist_fn_species_knobs | 0.0 | `bd` in nl_center | SUPPORT |
| `vparknob` | dist_fn_knobs | 1.0 | scales exported `vpa` | SUPPORT (via export) |
| `layout` | layouts_knobs | 'lxyes' | 'yxels'/'yxles' may select accel path | SUPPORT (export `accelerated`) |
| `fft_measure_plan` | layouts_knobs | `.false.` | FFTW_PATIENT ⇒ run-to-run non-reproducible | GATE `.true.` for reference runs |
| `fft_use_wisdom` / `fft_wisdom_file` | layouts_knobs | `.true.` / 'default' | may import PATIENT wisdom (P3) | SUPPORT; recommend `.false.` for reference data |
| `opt_local_copy` | layouts_knobs | `.false.` | local copy also moves forbidden points (zeros) | SUPPORT (no numeric effect, §4.5.3) |
| `unbalanced_xxf/yxf`, `nproc_xxf_lo/yxf_lo`, `opt_redist_*` | layouts_knobs | off | data distribution only | SUPPORT (no numeric effect except P4) |
| `include_explicit_source_in_restart` | init_g_knobs | `.true.` | gexp saved/restored in restart | SUPPORT (export gexp) |
| `nperiod` | theta_grid_parameters | — | >1 only warns | SUPPORT |
| `mult_imp` | dist_fn_knobs | — | forced off for nonlin | n/a |
| `LOWFLOW` | compile flag | off | extra branch in add_explicit_terms (:431-434) | GATE |
| `SHMEM` | compile flag | off | SHMEM transform paths | GATE |
| `OPENMP` | compile flag | off | threads; FFTW threads | SUPPORT (no reductions affected; P4) |

---------------------------------------------------------------------------------------------------------------------

## 7. GPU mapping notes

* **Loads, g_adjust, ×(i k), ×kxfac**: pure elementwise broadcasts over `(is, ie, il, ik, it, isgn, ig)`.
  `aj0/aj1/vperp2` broadcast along `isgn`; `phi/apar/bpar` broadcast along `(is, ie, il, isgn)`; `akx` along `it`,
  `aky` along `ik`. No masks needed (forbidden points are exact zeros in the coefficient arrays). `g_adjust` keeps
  its `trapped and il > jend[ig]` mask (exactness at the `1-al*bmag == -bouncefuzz` edge, P8).
* **FFTs**: batched over the leading line dimensions `(is, ie, il, isgn, ig)`; transform axes must be trailing
  `(…, kx/x, ky/y)`. Either transpose once into an internal "NL layout" `(nspec, negrid, nlambda, 2, nth, ntheta0,
  naky)` or use `moveaxis` per call (XLA materialises the transpose). 4 forward (c2r) and 1 inverse (r2c) 2D
  transforms per NL evaluation, 2 stages each.
* **Memory**: one real-space array = `nspec·negrid·nlambda·2·nth·nx·ny` float64 (e.g. 2·16·24·2·33·96·96 ≈ 4.7e8 ⇒
  3.7 GB). Peak live set ≈ 3 real arrays (`ba, gb, bracket`). Chunk over `(is, ie)` (or `il`) with `lax.map`; the
  only coupling between chunks is the CFL `max`, reduced at the end (exact).
* **Electrostatic optimisation (bit-exact)**: when `apar` is absent χ is independent of `isgn`; the two `isgn` lines
  of `chi_x, chi_y` are identical inputs, so transforming once and broadcasting gives identical results (provided the
  FFT is deterministic per line — true for cuFFT/ducc with a fixed plan).
* **Bracket**: elementwise; **CFL**: `jnp.max(jnp.abs(.))` over full real-space arrays (no masks).
* **nl_center**: vectorised shift along `ig` (§4.8) with three precomputed masks of shape `(nlambda, nth)`:
  `forbid`, `is_ttp`, `cell_forbidden[il, ig] = trapped & (il > jend[ig] | il > jend[ig+1])` for `ig < ntgrid`.
* **No sequential recurrences** inside the NL term (no `lax.scan`); no linked-boundary chains (the term is local
  in θ; `nl_center` couples only `ig, ig+1` within one line, no wrap-around).
* **History shuffle**: functional `(gexp_1, gexp_2, gexp_3) ← (NL, gexp_1, gexp_2)`; keep `istep_last` in the carry
  and execute the NL evaluation under `lax.cond(istep != istep_last, …)` so that the second implicit half-step and a
  redone (reset) step reuse it.
* **Reset control flow**: the dt-reset decision is a host-side branch (changes array shapes/coefficients of other
  subsystems); structure the JAX driver as a Python loop over steps with jitted step bodies, returning
  `(code_dt_cfl, reset)` to the host.
* **Split NL** (if ever un-gated): `lax.while_loop` with dynamic termination; stage arrays `(nstages, …)` stacked.

---------------------------------------------------------------------------------------------------------------------

## 8. Pitfalls (round-off faithfulness)

* **P1 — FMA contraction in the Fortran reference.** gfortran 15 `-O3` on arm64 contracts `bracket - ba*gb`
  (verified `fmls/fmsub` in a test compile of the same expression), and probably also `g1 - fac*vpa`,
  `(1+bd)*a + (1-bd)*b`, `adj = left + right` (bpar case), the AB sums. XLA may or may not contract. Expect 1-ulp
  differences in real space ⇒ ~1e-16 relative in the NL term. Options: rebuild the export GS2 with
  `-ffp-contract=off`, or accept. Products with an exactly-zero partner (`zi*k`, `(0,c)*φ`) are FMA-safe.
* **P2 — FFT algorithms differ.** FFTW (ESTIMATE plan, codelets) vs cuFFT/ducc give different round-off,
  ~`eps·log2(n)·‖input‖`. Bitwise agreement is impossible; compare the NL term normwise
  (`max|Δ|/max|ref|`), expect ≲1e-14. Element-wise relative error on near-zero modes will be larger.
* **P3 — FFTW plan reproducibility.** `fft_measure_plan=.true.` (FFTW_PATIENT) is non-deterministic run-to-run.
  With `fft_use_wisdom=.true.` (default) wisdom from an earlier PATIENT run can be imported and used by ESTIMATE
  requests. For reference data use `fft_measure_plan=.false., fft_use_wisdom=.false.` (or a fresh wisdom file).
* **P4 — nproc/layout dependence of the reference.** NL has no sum reductions across processes, but the FFTW plan's
  `howmany` (`nb_ffts`, `gs2_transforms.fpp:307-313`) and the accel switch depend on nproc/layout, which can change
  FFTW's algorithm (and FFTW threads with OpenMP). Generate reference dumps with a fixed configuration
  (recommended nproc=1, `layout='lxyes'`).
* **P5 — Chaotic divergence.** Nonlinear trajectories amplify 1e-16 differences exponentially; 1e-12 agreement of
  `g` is only achievable for a limited number of steps (O(10–100) depending on growth rate). Verify single-step NL
  evaluations against dumps rather than long trajectories.
* **P6 — ky=0 Hermiticity.** FFTW c2r discards `Im` of the (x-transformed) ky=0 column (equivalently Hermitian
  projection in kx). Must be explicit in JAX (§5). The inverse r2c returns exactly real DC; the subsequent x c2c then
  produces ky=0 coefficients that are conjugate-symmetric only to round-off — do not "fix" them.
* **P7 — in-place halving.** `transform2_5d` halves its input for `ik≠1` in place; `add_nl` relies on overwriting
  afterwards except for the stale forbidden entries of `g1` (which are zero). The JAX implementation must not halve
  `g`/state arrays.
* **P8 — two forbidden criteria.** `forbid` uses `1-al*bmag < -bouncefuzz`, `jend` counts `> -bouncefuzz`
  (`bouncefuzz = 10*epsilon(0.0)`, `le_grids.f90:131`). At exact equality `forbid=.false.` but `il > jend`. Use
  `forbid` where Fortran uses `forbid` (nl_center `where`) and `jend` where Fortran uses `jend` (g_adjust,
  nl_center cell test). Note `jend ≡ 0` when the grid has no trapped particles and is only consulted under `trapped`.
* **P9 — dt sequence.** dt changes by repeated `/delt_adj` and `min(*delt_adj, code_dt_max)`; with non-power-of-2
  `delt_adj` the values are path-dependent (`(x/1.5)*1.5 ≠ x`). Replicate the exact arithmetic (§4.9.4) and key any
  dt-dependent cache by the exact float64 value. Knife-edge CFL decisions (`code_dt ≈ code_dt_cfl` within FFT
  round-off) can flip between codes: log `|code_dt-code_dt_cfl|/code_dt` and flag < 1e-10.
* **P10 — after a reset the NL term is not recomputed** but the fields are (default `force_maxwell_reinit`); the
  redone step uses the pre-reset `gexp_1` together with recomputed `phi`. `istep_last` must survive the reset.
* **P11 — ttp and ig=ntgrid entries are not cell-centred.** ttp points keep the raw grid value (×1, not ×2);
  `gexp(ntgrid)` keeps the raw grid value. Both must be reproduced (ttp is consumed at `dist_fn.fpp:4987-5010`).
* **P12 — `max_vel_components(3)=EPS` floor** keeps `max_vel>0`; reproduce it, including `1/(2·EPS)` when limits are
  off, and `dt_cfl=1e8` for linear runs (every `add_explicit_terms` call).
* **P13 — scale constant.** Non-accel `(1/nx)*(1/ny)` vs accel `1/(nx*ny)`; and the `×2` is applied to the scale
  before multiplying (`g*(scale*2.0)`), then `×kxfac` separately.
* **P14 — filtered modes participate.** `kwork_filter(1,1)` and `ieqzip` modes are not evolved, but their (frozen)
  `g` and fields enter the NL products; the JAX state must carry them unchanged. The NL output for these modes is
  computed and stored but not used.
* **P15 — `use_cfl_limit` side effect** in `calculate_current_nl_source_and_cfl_limit` (split path only).
* **P16 — AB3 history at start/restart.** Orders 1,2,3 on steps 1,2,3 of a fresh run; on restart with
  `delt_option='check_restart'` and gexp in the restart file, order 3 immediately. `estimate_error` returns
  `error_target*cfl*cfl` until 3 history levels exist.
* **P17 — exb_shear before the NL term.** Flow-shear remapping (`fields_implicit.f90:335`) acts on `gnew/phinew`
  before they are copied into `g/phi`; the NL term uses unshifted `akx(it)` (no `kx_shift` in the NL term).
* **P18 — MPI aliasing/non-blocking.** With `immediate_reset=.false.` the max-reduction completes at end of step;
  results identical to the blocking path (max is exact).

---------------------------------------------------------------------------------------------------------------------

## 9. Open questions

1. Does XLA (CPU and CUDA backends, the versions to be used) contract `a - b*c` into FMA, and can it be disabled?
   Decide whether to rebuild the reference GS2 with `-ffp-contract=off` (affects all subsystems, not just NL).
2. cuFFT C2R with non-zero `Im(DC)`: the recipe zeroes it explicitly; confirm on the target GPU that the explicit
   zeroing yields identical output to the CPU path (unit test).
3. Is the ky=0 sector of GS2's `g`/fields kept exactly Hermitian in kx by any other routine? (Not needed for the
   recipe, but relevant for diagnostics comparisons.)
4. Whether to un-gate `split_nonlinear` later: it needs `get_fields_direct_from_dfn` /
   `calculate_potentials_from_nonadiabatic_dfn` (fields spec) and it clobbers `code_dt_prev1/2`.
5. FFTW plan dependence on `howmany`/strides for the 1D x/y transforms with FFTW_ESTIMATE on this build: not
   verified empirically; affects only cross-nproc reproducibility of the reference.
6. (Resolved) `aj1` convention: the local `spfunc::j1` returns J1(z)/z (`externals/utils/spfunc.fpp:369-371`); the NAG
   build variant (`:695`) should be checked if ever used. Irrelevant for the port since `aj1` is exported.
7. Equality edge `1-al*bmag == -bouncefuzz` (P8): confirm no exported test case hits it (assert
   `forbid == (il > jend)` for trapped grids in the exporter).

---------------------------------------------------------------------------------------------------------------------

## 10. Suggested verification hooks for the export program

1. Dump after `add_nl` (before `nl_center`) and after `nl_center` for `istep = 1, 2, 3` (canonical layout).
2. Dump `ba` after step (B) and `bracket` after (J) once (layout §3.3) to validate `gs2_transform2` separately.
3. Dump `max_vel_components` before/after reduction, `code_dt_cfl`, `reset`, `code_dt`, `code_dt_prev1/2` per step.
4. A synthetic unit test like `tests/unit_tests/nonlinear_terms/test_nonlinear_terms.fpp` (cos(jx x)cos(jy y)
   inputs) to check the it/ik ↔ FFT-position mapping and the 0.5/2 convention.
5. Test a non-Hermitian ky=0 input to confirm the Im(DC) drop matches FFTW (§4.5.1, P6).
