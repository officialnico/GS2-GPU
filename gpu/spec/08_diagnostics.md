# 08 — Diagnostics and `<run>.out.nc` output (contract for the JAX port)

Status: read-only analysis of GS2 8.2.1 (`gpu-port` branch, commit 4d8c94bcf). Nothing in `src/` was modified.

Scope: everything the **new diagnostics** (`src/diagnostics/gs2_diagnostics_new.f90`, the default since
`knobs:use_old_diagnostics = .false.`, `src/run_parameters.f90:122`; built with `-DNEW_DIAG`, `Makefile:388-389`)
computes during a run and writes to `<run>.out.nc`, plus any diagnostic **side effects on the evolved state**.
The JAX code must write a netCDF file with the **same variable names, dimension names, dimension order and
normalisation**, so that existing GS2 analysis scripts work.

All `file:line` citations are relative to `/Users/nico/Workspace/gs2/`.

Index conventions used here (same as the other specs):

* Distribution-function-like arrays use the canonical Fortran order
  `g(ig=-ntgrid:ntgrid, isgn=1:2, it=1:ntheta0, ik=1:naky, il=1:nlambda, ie=1:negrid, is=1:nspec)`.
  In NumPy (C order) that is `g[is, ie, il, ik, it, isgn, ig+ntgrid]`, with shape `(nspec, negrid, nlambda, naky, ntheta0, 2, nth)`
  where `nth = 2*ntgrid+1`. Arrays without `isgn` (for example `aj0(ig,iglo)`) drop that axis:
  `aj0[is, ie, il, ik, it, ig]`.
* Fields `phi/apar/bpar(ig,it,ik)` → NumPy `phi[ik, it, ig]`, shape `(naky, ntheta0, nth)`.
* Per-mode flux arrays (`real flux(it,ik,is)` in Fortran) → NumPy `flux[is, ik, it]`.
* Moments `complex m(ig,it,ik,is)` → NumPy `m[is, ik, it, ig]`.
* Velocity and theta arrays keep their Fortran index order and say so explicitly: `w(ie,is)` → `w[is, ie]`,
  `wl(ig,il)` → `wl[il, ig]`, `energy(ie,is)` → `energy[is, ie]`, `forbid(ig,il)` → `forbid[il, ig]`,
  `jend(ig)` → `jend[ig]`, theta-grid arrays `x(ig)` → `x[ig+ntgrid]`.
* netCDF variables are written by Fortran (neasyf) with the dimension list in Fortran order, so `ncdump`
  and netCDF4-python show them **reversed**. Example: Fortran `dim_names=[ri, kx, ky, t]` → netCDF
  `omega(t, ky, kx, ri)`. That reversed order equals the NumPy canonical order above with `ri` as a
  trailing real/imaginary axis. The JAX writer must create variables with the reversed (C-order) tuple
  `("t","ky","kx","ri")` and write NumPy arrays of shape `(nt, naky, ntheta0, 2)`.

---

## 0. Executive summary for implementers

1. Diagnostics are run once at `istep = 0` (in `initialize_diagnostics`), then **after every time step**
   (`istep = 1..nstep`), then once more after the loop with `force=.true.` (`src/gs2_main.fpp:600, 780, 879`).
   Output records (`nout`) are written when `mod(istep,nwrite)==0`, when `exit` is set, or on the forced final
   call, but never twice for the same `istep` (`src/diagnostics/gs2_diagnostics_new.f90:426`).
2. **In linear runs the instantaneous complex frequency is computed on every step, not only on write steps**
   (`gs2_diagnostics_new.f90:417-420`), filling a ring buffer of length `navg`. Its convergence and instability tests
   can end the run on any step. The JAX loop must therefore run the omega update and the exit test every step.
3. **Diagnostics change the evolved state at round-off level.** `write_moments` (default `.true.`) and the flux
   calculation (default `.true.` for nonlinear runs) convert `gnew` from g to h and back with `g_adjust`, which
   leaves `gnew ← (gnew + adj) − adj`. That is not bit-identical to `gnew`, and it happens on every write step,
   **including istep 0**. The JAX code must apply the same round trips in the same order (Section 3.1.8, Pitfall P1).
4. Every θ integral is a trapezoid rule with a **second-order compensated (Klein) sum**
   (`externals/utils/integration.f90:55-71`, `externals/utils/summation.f90:71-103`). It is a strictly sequential
   recurrence, so implement it with `lax.scan` and never with `jnp.sum` or `jnp.trapz`.
5. Velocity-space integrals (`integrate_moment`) accumulate in **`iglo` order**, which depends on the `layout`
   string (`src/le_grids.f90:1321-1422`, `src/gs2_layouts.fpp:917-930`). Reference runs for round-off comparison
   must use **nproc = 1** and **no OpenMP**.
6. The binary is built for aarch64 with `-O3` and no `-ffp-contract=off`, so gfortran **fuses a*b+c into FMA**
   (confirmed with `otool -tv objs/diagnostics_fluxes.o`: 269 `fmadd`/`fmsub`). XLA on GPU also contracts, but in
   different places. Bitwise agreement is therefore not achievable. The target is about 1e-12 relative, with
   absolute floors where values are dominated by cancellation or are near zero (Pitfall P10).
7. Recommendation for static content: have the Fortran export program create `<run>.out.nc` with the
   **static part** (metadata, dimensions, grids, species, geometry, input file), using the same `gs2_io` routines
   (Section 4.1). The JAX code then opens it in append mode and writes only the time-dependent and "final" variables.

---

## 1. Control flow and call graph

### 1.1 Program level (`run_gs2`, `src/gs2_main.fpp:1640-1670`)

```
initialize_gs2 → initialize_equations
initialize_diagnostics(state)                         gs2_main.fpp:542-614
   init_gs2_diagnostics_new(header)                   gs2_diagnostics_new.f90:30-185
   run_diagnostics(0, state%exit)                     gs2_main.fpp:600   ← writes record nout=1
evolve_equations(state, nstep)                        gs2_main.fpp:652-…
   do istep = istep_end+1, istep_end+nstep_run       gs2_main.fpp:717
      advance(istep)                                  :730  (fields_implicit.f90:315-371 etc.)
      [reset path: reset_time_step + advance again]   :739-752
      if (state%exit) exit                            :754
      [restart checkpoint gs2_save_for_restart]       :759-769
      update_time                                     :771  (gs2_time.f90:161-174: user_time += user_dt)
      run_diagnostics(istep, state%exit)              :780
      if (.not. exit) check_time_step / reset / checkstop / checktime   :791-825
      exit if user_time + user_dt > max_sim_time      :833-836
      if (state%exit) exit                            :856-859
   end do
   istep_offset = 1; if (state%exit) istep_offset = 0 :873-874
   run_diagnostics(istep - istep_offset, exit, force=.true.)   :879
finalize_diagnostics → finish_gs2_diagnostics_new     gs2_main.fpp:977-1015, gs2_diagnostics_new.f90:235-296
```

Notes:

* After a normal loop completion, Fortran's DO variable is `istep_loop_max+1`, so the forced call passes
  `istep = nstep` (offset 1). After an explicit `exit`, it passes the step on which the run stopped (offset 0).
* The forced call recomputes omega for that same `istep`. The result is idempotent: same `phi`, `phinew`, `code_dt`,
  same history slot. It writes a record only if that `istep` was not already written. This produces the familiar
  extra last record when `nstep` is not a multiple of `nwrite`. In the verification run below, `nstep=45`,
  `nwrite=10` gave `t = 0, 0.5, 1, 1.5, 2, 2.25`.
* `state%skip_diagnostics` (`gs2_main.fpp:213`) bypasses all calls. It is not settable from an input file.

### 1.2 State at diagnostic time

`advance_implicit` (`src/fields_implicit.f90:315-371`) does:

1. `exb_shear(gnew, phinew, aparnew, bparnew)` when flow shear is on (:335).
2. `g ← gnew`, `phi ← phinew`, `apar ← aparnew`, `bpar ← bparnew` (:337-338).
3. Advance to the new `gnew` and `phinew/aparnew/bparnew`.

So when `run_diagnostics(istep)` runs:

| Fortran name | meaning at diagnostic time |
|---|---|
| `gnew`, `phinew`, `aparnew`, `bparnew` | state at the new time level n+1 |
| `g`, `phi`, `apar`, `bpar` | state at time level n (after the ExB remap of the previous step, if any) |
| `code_dt` | the dt actually used for the step just taken (after any in-advance reset) |
| `user_time` | t(n+1) = accumulated sum of `user_dt` (`gs2_time.f90:171`), **not** `istep*dt` |

At `istep = 0`, init sets `phi = phinew` and `g = gnew` (`src/init_g.f90:3979, 4003, 4029`,
`src/gs2_init.fpp:366-368`).

### 1.3 `init_gs2_diagnostics_new` (`gs2_diagnostics_new.f90:30-185`)

1. `init_diagnostics_config` → `read_parameters` (`src/diagnostics/diagnostics_config.f90:267-395`):
   * "Smart defaults" (:279-285) apply only if the namelist object has not already been initialised. When
     `nonlin`, it sets `write_fluxes = write_fluxes_by_mode = .true.`. `apply_new_defaults` (:252-265) then sets
     `navg=10, nmovie=-1, nwrite=10, ob_midplane=F, omegatinst=1e6, print_line=F, write_correlation=T, write_moments=T`.
   * Namelist `gs2_diagnostics_knobs` is read (:287).
   * Overrides: `use_le_layout ⇒ write_cerr=F` (:385). `.not.(save_for_restart .or. save_distfn) ⇒ nsave=-1` (:387-389).
2. `check_parameters` (:188-208):
   * `print_line .or. write_line` forces `write_fields=T` and `write_omega=T`. **`write_line` defaults to `.true.`**,
     so by default `write_fields` and `write_omega` cannot be switched off.
   * `print_flux_line .or. write_flux_line` forces `write_fields=T`.
   * `write_jext` without `fapar > epsilon` stops the run.
3. `if (.not. write_any) return` (:87). With `write_any=F` there is no file, no omega, no convergence test and no g_adjust side effects.
4. `fluxfac(naky)` is allocated and never used (:109-113).
5. Heating interplay (:115-121): `write_heating .and. .not. collisions:heating` disables `write_heating`.
   `heating .and. .not.(write_heating .or. write_collisional)` calls `set_heating(.false.)`, so the collision
   operator stops accumulating heating. **The collisions spec must take `heating` as effectively false unless one
   of these two diagnostics is on.**
6. ASCII files are opened if `write_ascii` (default T). `.out` is always opened when `write_ascii`
   (`gs2_diagnostics_new.f90:216-232`, `src/diagnostics/diagnostics_ascii.f90:43-60`).
7. Sub-module init:
   * `init_diagnostics_fluxes` zeroes all flux arrays (`src/diagnostics/diagnostics_fluxes.fpp:36-75`).
   * `init_diagnostics_omega` allocates `omegahist(0:navg-1,ntheta0,naky)=0`, `domega(navg,ntheta0,naky)` and
     `omega_average=0` (`src/diagnostics/diagnostics_omega.f90:26-36`).
   * velocity-space, antenna, transfer and collisional modules; heating and nonlinear-convergence modules when enabled.
8. File creation (:144-174). With `append_old .and. exists` the file is reopened and `nout = len(t)`. Otherwise:
   `create_metadata`, `define_dims`, `nc_norms`, `nc_species`, `nc_geo`, `save_input`, `nout = 1`, and with
   `make_movie` also `init_transforms` and `nc_grids_mymovie`.

### 1.4 `run_diagnostics(istep_in, exit, force)` (`gs2_diagnostics_new.f90:339-519`)

```
if (.not. write_any) return
gnostics%exit = exit
writing = proc0 ; create = (istep_in==-1 .and. .not.appending)   ! istep_in is never -1 in gs2_main
istep = max(istep_in, 0)
calculate_fluxes = write_fluxes .or. print_flux_line .or. write_flux_line      (:410-412)
user_time recorded; if istep==0: start_time = user_time                        (:414-415)
if (istep > 0):
     calculate_omega(gnostics)              ! every step (linear, or nonlinear with any(ieqzip))
     if write_heating: calculate_heating    ! GATED
if (istep /= istep_last) .and. (mod(istep,nwrite)==0 .or. gnostics%exit .or. force):   (:426)
     if istep_in > -1: istep_last = istep
     write "t"(nout) = user_time ;  if wstar_units: "t_wstar"(ky,nout) = user_time*tunits   (:434-440)
     vary_vnew_only = F
     -- ORDER MATTERS (side effects on gnew) --
     write_omega        → write_omega                    (:445)
     write_fields       → write_fields                   (:447)
     dump_fields_periodically → GATED                    (:449)
     write_ql_metric    → write_ql_metric                (:450)
     calculate_fluxes   → calculate_fluxes  [g_adjust round trip on gnew AND g]   (:451)
     write_symmetry     → GATED                          (:453-456)
     write_nl_flux_dist → GATED  [g_adjust round trip on gnew]                     (:457)
     write_parity       → GATED  [g_adjust round trip on gnew]                     (:458)
     write_verr         → GATED  [round trips + may change vnmult]                 (:459)
     write_cerr         → GATED  [round trip]                                      (:460)
     write_moments      → write_moments  [g_adjust round trip on gnew]             (:462)
     write_full_moments_notgc → GATED                    (:463)
     make_movie         → GATED                          (:464)
     write_heating      → GATED                          (:465)
     nonlin .and. use_nonlin_convergence → check_nonlin_convergence  (:466)
     write_cross_phase .and. has_electron_species → GATED [round trip via getemoms]  (:467)
     write_jext         → GATED                          (:468)
     write_correlation  → do_write_correlation           (:469)
     write_correlation_extend .and. istep>nstep/4 .and. mod(istep,nwrite_mult*nwrite)==0 → GATED (:470-473)
     write_lorentzian   → GATED                          (:474)
     write_eigenfunc    → write_eigenfunc                (:475)
     print_line / write_line (ascii)                     (:477-478)
     proc0: print_flux_line / write_flux_line (ascii)    (:479-482)
     write_zonal_transfer / write_kinetic_energy_transfer / write_collisional → GATED (:484-497)
     run_diagnostics_to_be_updated  (write_gyx/write_g ascii, every nwrite*nwrite_mult) → GATED (:499, 299-316)
     nc_sync every nc_sync_freq records; flush ascii
     user_time_old = user_time                          (:506)
     nout = nout + 1                                    (:508)
else if mod(istep, ncheck)==0:                          (:509-514)
     vary_vnew_only = T; if write_verr .and. vary_vnew → write_velocity_space_checks   ! GATED
exit = gnostics%exit
```

`istep_last` is a `save` variable initialised to −1 (:377).

### 1.5 Finalisation (`finish_gs2_diagnostics_new`, `gs2_diagnostics_new.f90:235-296`)

1. `save_restart_dist_fn(...)`: restart and dist-fn files (not in `.out.nc`; see the restart spec).
2. `run_old_final_routines` (:544-602):
   * `write_final_fields` → `do_write_final_fields` (writes `phi/apar/bpar` = new fields; `.fields` ascii uses the OLD fields).
   * `write_kpar` (ascii).
   * `write_final_epar` → `epar`.
   * `write_final_db` (ascii).
   * `phi0 = get_phi0()` (always).
   * `write_final_moments` → `do_write_final_moments` (getmoms round trip, netCDF `ntot`, `density`, …).
   * `write_final_antot` → `getan` then `antot`, `antota`, `antotp`.
   * `dump_ant_amp` (antenna).
   * `nonlin .and. write_gs` (ascii).
3. `finish_*` of the sub-modules; `neasyf_close`.
4. ASCII: `write_gyx`, `write_g` dumps; `do_write_geom` (writes `<run>.g`) when `write_ascii`.

---

## 2. Inventory of everything read (arrays and scalars)

Legend: **EXPORT** = dump from Fortran. **RECOMPUTE** = port the formula (given). **STATE** = part of the evolving
JAX state (produced by the time stepper). "dt?" = depends on `code_dt` or changes at runtime.

### 2.1 Static grids and geometry

| Fortran name | module | decl. shape / ranges | computed in | dt? | decision / formula |
|---|---|---|---|---|---|
| `theta` | theta_grid | `(-ntgrid:ntgrid)` | geometry init (`theta_grid.f90`) | no | EXPORT |
| `jacob` | theta_grid | `(-ntgrid:ntgrid)` | `theta_grid.f90:516`: `jacob = 1.0/(drhodpsi*gradpar*bmag)` | no | EXPORT (recompute only for checking) |
| `field_line_average_weight` | theta_grid (PRIVATE) | scalar | `theta_grid.f90:520`: `1.0/trapezoidal_integration(theta, jacob)` | no | BOTH: the export program evaluates the same expression through the public `integration` module and dumps it; JAX recomputes it with the Klein trapezoid and asserts bitwise equality |
| `grho` | theta_grid | `(-ntgrid:ntgrid)` | geometry | no | EXPORT |
| `bmag`, `gradpar`, `Bpol`, `Rplot` | theta_grid | `(-ntgrid:ntgrid)` | geometry | no | EXPORT |
| `gds21`, `gds22` | theta_grid | `(-ntgrid:ntgrid)` | geometry | no | EXPORT |
| `gbdrift`, `gbdrift0`, `cvdrift`, `cvdrift0` | theta_grid | `(-ntgrid:ntgrid)` | geometry | no | EXPORT (used by `wdrift_func` in the energy exchange) |
| `delthet` | theta_grid | `(-ntgrid:ntgrid)` | `theta_grid.f90:510-511`: `delthet(:ntgrid-1)=theta(-ntgrid+1:)-theta(:ntgrid-1)`; `delthet(ntgrid)=0` | no | EXPORT |
| `rhoc`, `qval`, `shat`, `drhodpsi` | theta_grid | scalars | geometry | no | EXPORT |
| `aky(naky)`, `akx(ntheta0)` | kt_grids | 1-D | kt_grids init | no | EXPORT |
| `theta0(ntheta0,naky)` | kt_grids | 2-D | kt_grids init | no | EXPORT |
| `kperp2(-ntgrid:ntgrid,ntheta0,naky)` | kt_grids | 3-D | `kt_grids.f90:313-349` | no (flow shear remaps indices, not kperp2) | EXPORT |
| `w(negrid,nspec)` | le_grids | `w(ie,is)` | le_grids init | no | EXPORT |
| `wl(-ntgrid:ntgrid,nlambda)` | le_grids | `wl(ig,il)` | le_grids init | no | EXPORT |
| `energy(negrid,nspec)` | le_grids | `energy(ie,is)` | le_grids init | no | EXPORT |
| `al(nlambda)` | le_grids | 1-D | le_grids init | no | EXPORT (`lambda` coordinate) |
| `forbid(-ntgrid:ntgrid,nlambda)` | le_grids | logical `(ig,il)` | le_grids init | no | EXPORT |
| `jend(-ntgrid:ntgrid)` | le_grids | integer | le_grids init | no | EXPORT |
| `ng2`, `nlambda` | le_grids | ints | le_grids | no | EXPORT; `trapped = nlambda > ng2` (`le_grids.f90:443-446`) |
| `aj0`, `aj1`, `vperp2` | dist_fn_arrays | `(-ntgrid:ntgrid, iglo)` → `(ig,it,ik,il,ie,is)` | dist_fn init (Bessel) | no | EXPORT |
| `vpa`, `vpac` | dist_fn_arrays | `(-ntgrid:ntgrid,2,iglo)` → canonical | `dist_fn.fpp:1592-1678`; `vpac(ntgrid,:,:)=0` | no | EXPORT |
| `spec(is)%{z,mass,dens,temp,tprim,fprim,stm,zt,tz,smz}` | species | per species | `species.f90:476-480, 1009-1013` | no | EXPORT the derived values (`stm=sqrt(temp/mass)`, `zt=z/temp`, `tz=temp/z`, `smz=abs(sqrt(temp*mass)/z)`) |
| `spec(is)%type` | species | int | species | no | EXPORT (`has_electron_species`; static netCDF) |
| `nonmaxw_corr(negrid,nspec)` | species | `(ie,is)` | `species.f90:571-591` | no | EXPORT |
| `fphi, fapar, fbpar` | run_parameters | scalars | input | no | EXPORT |
| `has_phi/has_apar/has_bpar` | run_parameters | logical | `run_parameters.f90:351-353`: `f* > epsilon(0.0)` | no | RECOMPUTE (as written) |
| `nonlin` | nonlinear_terms | logical | input | no | EXPORT |
| `ieqzip(ntheta0,naky)` | run_parameters | logical | `run_parameters.f90:259-272` | no | EXPORT (`any(ieqzip)` re-enables omega in nonlinear runs) |
| `tunits(naky)`, `woutunits(naky)`, `wunits(naky)` | gs2_time | 1-D | `gs2_time.f90:83-126`: `tunits=1`, `wunits=aky/2.0`, `woutunits=1.0/tunits` | no (yes only if `wstar_units`; gated) | RECOMPUTE (trivial) or EXPORT |
| `def_parity`, `even` | dist_fn | logical | input | no | EXPORT (`get_phi0`) |
| `layout` | gs2_layouts | char(5) | input | no | EXPORT (sets summation order) |
| `igomega` | diagnostics | int (default 0) | input | no | EXPORT; must satisfy `-ntgrid ≤ igomega ≤ ntgrid` |

### 2.2 Dynamic quantities (state)

| Fortran name | module | shape | produced by | dt? | decision |
|---|---|---|---|---|---|
| `gnew`, `g` | dist_fn_arrays | `(-ntgrid:ntgrid,2,iglo)` | time stepper | yes | STATE (diagnostics **mutate** `gnew` and `g` at round-off level; Section 3.1.8) |
| `phinew, aparnew, bparnew, phi, apar, bpar` | fields_arrays | `(-ntgrid:ntgrid,ntheta0,naky)` | time stepper | yes | STATE |
| `code_dt`, `user_dt`, `user_time` | gs2_time | scalars | stepper / reset | yes | STATE (`user_time` accumulated exactly as `user_time = user_time + user_dt`) |
| `hypervisc_filter(-ntgrid:ntgrid,ntheta0,naky)` | hyper | 3-D | `hyper.f90:612-613` (=1.0 when hyper is off); overwritten each step by `hyper_diff` (`hyper.f90:435-441, 738-741, 777-778`) | **yes** (code_dt, and \|φ\| for `iso_shear`/`aniso_shear`) | STATE: the value left by the hyper step of the current advance. Needed for `es_energy_exchange` |
| `omegahist(0:navg-1,ntheta0,naky)` | diagnostics_omega | complex | `calculate_omega` | yes | STATE (ring buffer, initialised to 0; **not** saved in restart files) |
| `omega_average(ntheta0,naky)` | diagnostics_omega | complex | `calculate_omega` | yes | STATE |
| `istep_last`, `nout`, `user_time_old`, `start_time` | diagnostics | scalars | run_diagnostics | — | STATE |
| `current_results%*` | diagnostics_config | scalars and `(nspec)` | fluxes and fields | — | STATE (used for ascii flux line, Trinity outputs) |
| `vnmult(2)` | collisions | — | collisions | yes | only `vspace_vnewk` (gated) |

### 2.3 Scalars and parameters of the diagnostics namelist

`nwrite, navg, omegatol, omegatinst, exit_when_converged, igomega, write_*` flags, `nwrite_mult, ncheck, nc_sync_freq`
and the convergence-test parameters: EXPORT the parsed `gnostics` values *after* `read_parameters` and
`check_parameters`, so that forced flags are captured. Defaults are in Section 5.

---

## 3. Exact computations

The Fortran expression trees are quoted, then given as left-to-right pseudo-code. Fortran evaluates `a*b*c` as
`(a*b)*c` and `a*b/c` as `(a*b)/c`. Unary minus binds weaker than `*` and `/`: `-a*b` = `-(a*b)`, but negation
is exact so its position does not matter. `x**2` with a literal integer 2 is `x*x`. `x**2.0` with a real exponent
2.0 is folded by gfortran to `x*x`; use `x*x` in JAX, not `jnp.power`.

### 3.1 Primitives (shared by everything below)

#### 3.1.1 Klein compensated sum — `summation.f90:71-103`

```
res = 0; c1 = 0; c2 = 0
for i in 0..n-1:                 # strictly sequential
    t = res + a[i]
    if abs(res) >= abs(a[i]): lc1 = (res - t) + a[i]
    else:                     lc1 = (a[i] - t) + res
    res = t
    t = c1 + lc1
    if abs(c1) >= abs(lc1):   lc2 = (c1 - t) + lc1
    else:                     lc2 = (lc1 - t) + c1
    c1 = t
    c2 = c2 + lc2
return res + c1 + c2             # = (res + c1) + c2
```

`summation.o` contains no FMA (checked with `otool`), so this sequence is exact IEEE and can be reproduced bitwise.
In JAX: `lax.scan` over the summation axis carrying `(res, c1, c2)`, vmapped or broadcast over the other axes,
with `jnp.where` for the two branches. Do **not** enable XLA fast-math, which would allow reassociation and
destroy the compensation.

#### 3.1.2 Trapezoid — `integration.f90:55-71`

```
real:    trapz(x, y)  = klein_sum( (x[1:]-x[:-1]) * (y[:-1]+y[1:]) ) * 0.5
complex: trapz(x, yc) = complex( klein_sum(Re(d)), klein_sum(Im(d)) ) * 0.5,   d = (x[1:]-x[:-1]) * (yc[:-1]+yc[1:])
```

The real×complex product is evaluated by gfortran as a full complex multiply with a promoted `(dx,0)`, using FMA
(`otool` on `integration.o`). Numerically this equals `(dx*Re, dx*Im)` except for the sign of zeros.
`integrate_moment` results and fields are never ±0 sensitive here.

#### 3.1.3 Field-line (θ) average — `theta_grid.f90:657-670`

```
W  = 1.0 / trapz(theta, jacob)                      # theta_grid.f90:520
flavg(q) = trapz(theta, q * jacob) * W              # q real or complex; q*jacob elementwise first
```

The θ grid is the full `-ntgrid:ntgrid` extended grid, including both ±π end points of each segment, with no
special treatment of the duplicated boundary point.

#### 3.1.4 Flux-surface weight used by `get_flux` — `diagnostics_fluxes.fpp:825`

```
wgt = 1.0 / trapz(theta, grho * jacob)     # "essentially 2 pi / surfarea"
```

#### 3.1.5 ky/kx "volume averages" — `src/diagnostics/volume_averages.f90`

GS2 stores only ky ≥ 0 with non-standard Fourier coefficients. Every mode with `aky(ik) /= 0` gets weight 0.5 and
the `aky(ik) == 0.` mode gets weight 1.0. The test is an **exact compare on `aky`**, not on `ik==1`.

```
fac(ik) = 1.0 if aky[ik] == 0.0 else 0.5

average_ky_xy(f[it,ik]) -> favg[it]:      (volume_averages.f90:373-404)
    favg = 0
    for ik in 1..naky:  for it in 1..ntheta0:  favg[it] = favg[it] + f[it,ik]*fac(ik)
average_kx_xy(f[it,ik]) -> favg[ik]:      (:457-489)
    favg = 0
    for ik: for it: favg[ik] = favg[ik] + f[it,ik]*fac(ik)
average_ky_xys(f[it,ik,is]) -> favg[it,is] (:407-441)   for is: for ik: for it: favg[it,is] += f*fac(ik)
average_kx_xys -> for is: average_kx_xy(f[:,:,is])     (:444-454)
average_all_real_xys(a[it,ik,is]) -> avg[is]: (:112-126)
    for is: tmp[it] = average_ky_xy(a[:,:,is]); avg[is] = sum_seq_it(tmp)
average_theta_txy(a,b)[it,ik] = flavg( Re(conj(a[:,it,ik]) * b[:,it,ik]) )        (:258-278)
average_theta_txys(a,b)[it,ik,is] = flavg( Re(conj(a[:,it,ik,is]) * b[:,it,ik,is]) ) (:170-193)
average_theta_complex_complex_txys(a)[it,ik,is] = flavg( a[:,it,ik,is] )  (complex)  (:145-167, 237-243)
```

Here `sum_seq` is Fortran intrinsic `SUM`. gfortran evaluates it as a plain left-to-right loop with no
reassociation (no `-ffast-math`), so reproduce it as a sequential accumulation starting from 0.
`Re(conj(a)*a)` = `ar*ar + ai*ai`, computed as the real part of a complex product, possibly with FMA.
Multiplication by 0.5 or 1.0 is exact.

`distributed` is always `.false.` (`gnostics%parallel` is hard-wired false, `gs2_diagnostics_new.f90:20, 89-101, 408`),
so `field_k_local` and `sum_allreduce` branches never run.

#### 3.1.6 `integrate_moment` (complex, rank-3 in, rank-4 out) — `le_grids.f90:1321-1422`

```
total[ig,it,ik,is] = 0
for iglo in llim..ulim (layout order):        # (it,ik,il,ie,is) = decode(iglo)
    total[:,it,ik,is] = total[:,it,ik,is] + (w[ie,is]*wl[:,il]) * (g[:,1,iglo] + g[:,2,iglo])
```

* Expression: `w(ie,is)*wl(:,il)*(g(:,1,iglo)+g(:,2,iglo))` = `((w*wl) * (g1+g2))`, then added to the accumulator.
* For fixed `(it,ik,is)` the contributions arrive in increasing `iglo`. Using the dimension orders in
  `gs2_layouts.fpp:917-930` (first letter = fastest):
  * layouts `lxyes`, `xyles`, `yxles`, `lexys`, `lyxes`: `for ie: for il:` (il fastest);
  * layout `yxels`: `for il: for ie:`.
* All calls made by the supported diagnostics use `all=.false.` or the default. With `nproc=1` nothing is reduced.
  With `nproc>1` the result is a sum of per-rank partials in an MPI-implementation-defined order (not reproducible).
* JAX: `lax.scan` over the flattened `(ie,il)` sequence in the layout order, carrying
  `acc[is, ik, it, ig]` (complex) and adding `(w[is,ie]*wl[il,:]) * (g[is,ie,il,:,:,0,:] + g[is,ie,il,:,:,1,:])`.

#### 3.1.7 `get_flux(g_in, fld) -> flx[it,ik,is]` — `diagnostics_fluxes.fpp:801-838`

```
total = integrate_moment(g_in)                 # (ig,it,ik,is)
wgt   = 1.0/trapz(theta, grho*jacob)
for is: for ik: for it:
    y   = Im( total[:,it,ik,is] * conj(fld[:,it,ik]) ) * jacob
    flx[it,ik,is] = trapz(theta, y) * wgt * aky[ik]      # ((trapz*wgt)*aky)
flx = flx * 0.5
```

`Im(t*conj(f))` = `tr*(-fi) + ti*fr` (complex multiply; FMA possible).

#### 3.1.8 `g_adjust` (g ↔ h) — `src/dist_fn_arrays.fpp:101-175`

```
direction from_g_gs2:  facphi =  fphi, facbpar =  fbpar
direction to_g_gs2:    facphi = -fphi, facbpar = -fbpar
trapped     = nlambda > ng2
bpar_factor = 2.0*facbpar      if facbpar /= 0 else 0.0
phi_factor  = 0.0 ; if facphi /= 0: phi_factor = spec(is)%zt * nonmaxw_corr(ie,is) * facphi   # (zt*nm)*facphi
for each iglo (it,ik,il,ie,is), ig:
    if trapped and il > jend(ig): skip                  # forbidden region left untouched
    adj = ((bpar_factor*vperp2(ig,iglo))*aj1(ig,iglo))*bpar(ig,it,ik) + (phi_factor*phi(ig,it,ik))*aj0(ig,iglo)
    g(ig,1,iglo) += adj ; g(ig,2,iglo) += adj
```

Because negation is exact, `adj_to == -adj_from` bit for bit, so a round trip leaves **`g ← (g + adj) − adj`**,
which is not the identity in floating point.

**Round trips that run in the supported configuration, per write step, in this order:**

1. `calculate_fluxes` (only if `calculate_fluxes`), `diagnostics_fluxes.fpp:262-263, 279-280`:
   * `gnew` from (with `phinew, bparnew`), `g` from (with `phi, bpar`), compute;
   * then `gnew` to, `g` to.
2. `write_moments` → `getmoms` (`diagnostics_moments.fpp:223, 329`): `gnew` from then to (with `phinew, bparnew`).
3. `write_final_moments` at finalisation: once more (affects only final outputs).

The perturbation of `g` is harmless: `g` is overwritten by `copy(gnew, g)` at the start of the next advance,
`fields_implicit.f90:337`. The perturbation of **`gnew` propagates into the evolution**.
Gated diagnostics that would add further round trips: `write_nl_flux_dist`, `write_parity`, `write_verr`,
`write_cerr`, `write_cross_phase`, `write_heating`, `write_zonal_transfer`, `write_collisional`.

### 3.2 Time `t`

* `t(nout) = user_time`, written first in each write block (`gs2_diagnostics_new.f90:435`). The netCDF variable
  `t` has **no attributes**: it is created by `neasyf_write` on an unlimited dimension that `define_dims` created
  without a coordinate variable (`src/gs2_io.fpp:311`).
* `user_time` starts at `tstart` (from input or restart, `gs2_time.f90:129-136`) and evolves as
  `user_time = user_time + user_dt` after each accepted step (`gs2_time.f90:171`).
* `t_wstar(t, ky) = user_time * tunits` only if `wstar_units` (gated).

### 3.3 Complex frequency `omega`, `omega_average`, convergence — `src/diagnostics/diagnostics_omega.f90:64-136`

Called each step with `istep > 0`, before the write test:

```
if nonlin and not any(ieqzip): return            # omegahist stays 0; zeros are written
j   = igomega
fac = 1.0 if omegatol == 0 else 1000/abs(omegatol)          # real; 1000 integer → real
den = phi(j,:,:)    + apar(j,:,:)    + bpar(j,:,:)          # ((phi+apar)+bpar), OLD fields
num = phinew(j,:,:) + aparnew(j,:,:) + bparnew(j,:,:)       # NEW fields
k   = mod(istep, navg)
where (abs(num) < tiny(0.0)*fac .or. abs(den) < tiny(0.0)*fac):
    omegahist(k,:,:) = 0
elsewhere:
    omegahist(k,:,:) = log(num/den) * zi / code_dt          # ((log(q)*zi)/code_dt)
omega_average = sum( omegahist / real(navg), dim=1 )        # divide each slot first, then sequential sum over k=0..navg-1
current_results%omega_average = omega_average
if istep > navg:
    domega(m,:,:) = omega_average - omegahist(m-1,:,:)      # m = 1..navg
    if all( sqrt( sum_{m=1..navg} abs(domega)**2 / real(navg) ) <= min(abs(omega_average), 1.0) * omegatol ):
        write "*** omega converged" to .out (if write_ascii)
        gnostics%exit = exit_when_converged                  # NB assignment, can also set .false.
    if any( abs(omega_average) * code_dt > omegatinst ):
        print "*** numerical instability detected"; gnostics%exit = .true.
```

Details:

* `tiny(0.0)` = 2.2250738585072014e-308 (double, because of `-fdefault-real-8`). Unused fields are zero arrays,
  so the unconditional sum equals the `has_*`-guarded sum of `calculate_instantaneous_omega`.
* `log` is the principal complex log: Im in (−π, π]. `log(q)*zi = (−Im log q, Re log q)`, so
  Re(ω) = −arg(q)/dt (real frequency) and Im(ω) = ln|q|/dt (growth rate). The GS2 convention is
  φ ∝ exp(−iωt), so a positive Im(ω) means growth.
* `abs(domega)**2` is the complex modulus (hypot) squared, not `re²+im²`.
* **A negative `omegatol` disables the convergence exit**: the right-hand side is negative and `sqrt ≥ 0`.
  `fac` still uses `abs(omegatol)`. With `omegatol > 0`, a mode whose fields stay below the tiny threshold has
  `omegahist ≡ 0`, so `domega ≡ 0 ≤ 0` holds and the run is declared "converged" once `istep > navg`.
  Replicate this behaviour.
* For `istep ≤ navg` the average includes the zero-initialised slots. That is intentional bias; replicate it.
* The written values (`write_omega`, `src/gs2_io.fpp:1259-1294`):
  * `omega(t, ky, kx, ri)` = `omegahist(mod(istep,navg),it,ik) * woutunits(ik)`;
  * `omegaavg(t, ky, kx, ri)` and `omega_average(t, ky, kx, ri)` (identical duplicates) = `omega_average*woutunits`.
  * Units attribute `"aref/vtref"`. With `wstar_units=F`, `woutunits=1.0` and the multiply is exact.
* At `istep=0`, slot 0 is 0, so the first record of `omega` and `omega_average` is exactly 0.

### 3.4 Fields: `phi2`, `phi2_by_mode`, `phi0`, … — `src/gs2_io.fpp:1318-1390` via `diagnostics_fields.f90:13-40`

For each field `X ∈ {phi, apar, bpar}` with `has_X`, the input is `field = Xnew` (new time level):

```
X2_by_mode[it,ik] = flavg( Re(conj(field[:,it,ik]) * field[:,it,ik]) )    # average_theta_txy
if write_X_over_time: write X_t(t, ky, kx, theta, ri) = field
X2_by_kx[it] = average_ky_xy(X2_by_mode)          → var "X2_by_kx"(t, kx)
X2_by_ky[ik] = average_kx_xy(X2_by_mode)          → var "X2_by_ky"(t, ky)
X2           = sum_seq_it( X2_by_kx )             → var "X2"(t)   units "(T_ref/q_ref rho_ref/a_ref)**2" (same string for all 3 fields)
var "X0"(t, ky, kx, ri)                 = field(igomega,:,:)    # "phi(t) at outboard midplane" when igomega=0
var "X_igomega_by_mode"(t, ky, kx, ri)  = field(igomega,:,:)    # duplicate
var "X2_by_mode"(t, ky, kx)             = X2_by_mode
current_results%X2 = X2
nc_final_fields(file): overwrite "phi","apar","bpar"(ky, kx, theta, ri) with phinew/aparnew/bparnew (all present fields)
```

* long_names: `phi` "Electrostatic potential", `apar` "Parallel magnetic potential", `bpar` "Parallel magnetic field".
  Derived names are suffixed: " squared and averaged over theta and ky", " squared and averaged over theta and kx",
  " at ig=igomega", " per Fourier mode", " volume average", " over time".
* The "final" `phi/apar/bpar` variables use long_names "Electrostatic Potential" (units "T_ref/q_ref rho_ref/a_ref"),
  "Parallel Magnetic Potential", and "delta B Parallel" (`gs2_io.fpp:680-706`).
* **Final eigenfunctions in `.out.nc` are therefore `phinew` at the last written record, unnormalised.**
  Normalised versions exist only with `write_eigenfunc` (Section 3.9).
* `phi2` is the "volume average" ⟨|φ|²⟩ = Σ_kx Σ_ky fac(ky) · ∫|φ|² J dθ / ∫J dθ.

### 3.5 Quasilinear metric `ql_metric_by_mode` — `gs2_diagnostics_new.f90:319-333`, `src/gs2_diagnostics.fpp:2807-2957`

Computed on every write step when `write_ql_metric` (default `.true.`, linear and nonlinear alike):

```
# calculate_instantaneous_omega(ig=igomega, tolerance=omegatol)   (gs2_diagnostics.fpp:2897-2957)
fs  = (has_phi? phi(ig): 0)    [+ apar(ig) if has_apar] [+ bpar(ig) if has_bpar]     # OLD
fsn = (has_phi? phinew(ig): 0) [+ aparnew(ig) ...] [+ bparnew(ig) ...]               # NEW
too_small = tiny(0.0) if tol == 0 else tiny(0.0)*1000./abs(tol)    # NB (tiny*1000)/|tol| differs in rounding from calculate_omega
omega_inst = 0 where (abs(fs) < too_small or abs(fsn) < too_small) else log(fsn/fs)*zi/code_dt
growth = Im(omega_inst)

# calculate_simple_quasilinear_flux_metric_by_k(growth)          (:2827-2893)
ql = 0
norm = maxval_theta(abs(phinew)) if has_phi elif has_apar maxval(abs(aparnew)) elif has_bpar maxval(abs(bparnew)) else return 0
where norm > 0: norm = 1.0/norm
for X in (phi, apar, bpar) if has_X:                                     # in this order
    avg[it,ik] = get_field_weighted_average_kperp2(Xnew)
    where avg > 0: ql = (ql +) maxval_theta(abs(Xnew)) * norm / avg       # first term assigns, later ones add
lim = growth where growth > 0 else 0
ql[:, ik] = 0 for every ik with aky(ik) == 0
ql = ql * lim

# get_field_weighted_average_kperp2(field)                       (:2807-2824)
fsq = abs(field)*abs(field)                   # NB abs*abs, not re^2+im^2
avg = 0
for ik: for it:
    if maxval(abs(fsq[:,it,ik])) <= epsilon(0.0): cycle         # epsilon = 2.22e-16
    avg[it,ik] = flavg(fsq[:,it,ik]*kperp2[:,it,ik]) / flavg(fsq[:,it,ik])
```

At `istep=0`, `phi == phinew`, so `fsn/fs ≈ 1` and growth is round-off-sized (exactly 0 only if the complex
division happens to return exactly (1,0)). Compare record 1 with an **absolute** tolerance.

### 3.6 Fluxes — `src/diagnostics/diagnostics_fluxes.fpp:101-242`

Runs only when `calculate_fluxes = write_fluxes .or. print_flux_line .or. write_flux_line`.

#### 3.6.1 Sequence

```
zero current_results totals and species_{heat,momentum,particle}_flux
if not nonlin: write_diffusive_estimates          # uses omega_average (Section 3.6.6)
common_calculate_fluxes()                          # Section 3.6.2
if has_phi:  standard-properties for: es_heat_flux(qheat1), es_heat_flux_par(qheat2), es_heat_flux_perp(qheat3),
             es_part_flux(pflux), es_part_tormom_flux(pflux_tormom), es_mom_flux(vflux), es_mom_flux_par(vflux_par),
             es_energy_exchange(exchange), es_mom_flux_perp(vflux_perp)        -- in exactly this order
             write "es_flux_vs_e"(t, species, egrid)    ← written whenever calculate_fluxes (NOT gated on write_fluxes)
             [LOWFLOW build only: es_mom0, es_mom1 (bug: both use vflux0)] -- not in this build
if has_apar: apar_heat_flux(qmheat1), apar_heat_flux_par(qmheat2), apar_heat_flux_perp(qmheat3),
             apar_part_flux(pmflux), apar_mom_flux(vmflux); write "apar_flux_vs_e"
if has_bpar: bpar_heat_flux(qbheat1), bpar_heat_flux_par, bpar_heat_flux_perp, bpar_part_flux(pbflux), bpar_mom_flux(vbflux);
             write "bpar_flux_vs_e"
species_*_flux_avg += species_*_flux * (user_time - user_time_old)     # Trinity only, not in .out.nc
if write_fluxes: heat_flux_tot, hflux_tot (dup), mom_flux_tot, vflux_tot (dup), part_flux_tot, zflux_tot (dup)   [all (t)]
```

`gs2_has_lowflow` is false (no `-DLOWFLOW` in the build flags), so `pflux_tormom ≡ 0` and **`es_part_tormom_flux*`
variables are written as exact zeros**. This was verified in a test run.

#### 3.6.2 `common_calculate_fluxes` (:247-330)

```
g_adjust(gnew, phinew, bparnew, from) ; g_adjust(g, phi, bpar, from)      # gnew, g now hold h
flux(phinew, aparnew, bparnew, ...)                                       # 3.6.3
flux_vs_e(phinew, aparnew, bparnew, esflux_vs_e, apflux_vs_e, bpflux_vs_e) # 3.6.4
[lowflow: lf_flux]                                                         # not in this build
if has_phi: eexchange(phinew, phi, exchange_dummy, exchange)              # 3.6.5
g_adjust(gnew, phinew, bparnew, to) ; g_adjust(g, phi, bpar, to)
# species normalisation (proc0); associate-expressions are evaluated ONCE:
for is:
   DT  = dens*temp ; MOM = dens*sqrt(mass*temp) ; MOMP = (dens*mass)*stm
   if has_phi:
      qheat[:,:,is,k] *= DT  (k=1,2,3)
      pflux *= dens ; pflux_tormom *= dens ; esflux_vs_e[:,is] *= dens
      vflux *= MOM ; vflux_par *= MOM ; vflux_perp *= MOMP          # NB vflux was formed as par+perp BEFORE scaling
      exchange = (exchange*dens)*z
   if has_apar: qmheat *= DT ; pmflux *= dens ; vmflux *= MOMP ; apflux_vs_e *= dens
   if has_bpar: qbheat *= DT ; pbflux *= dens ; vbflux *= MOMP ; bpflux_vs_e *= dens
```

Because `vflux = vflux_par + vflux_perp` is formed before scaling and is then scaled with `MOM` while `vflux_perp`
uses `MOMP`, `es_mom_flux` ≠ `es_mom_flux_par + es_mom_flux_perp` at round-off level, since
dens·sqrt(m·T) and dens·m·sqrt(T/m) round differently.

#### 3.6.3 `flux(phi, apar, bpar)` (:539-739), with `gnew` = h

All arrays are zeroed first. Each block forms `g_work` and calls `get_flux` with the indicated field:

```
has_phi (field = phinew):
  pflux      : g_work = gnew*aj0
  qflux(.,1) : g_work = (gnew*aj0)*energy(ie,is)                    # in place on the previous g_work
  qflux(.,2) : g_work = ((gnew*2.)*(vpa*vpa))*aj0
  qflux(.,3) : g_work = (gnew*vperp2)*aj0
  vflux_par  : g_work = (((gnew*aj0)*vpac)*Rplot)*sqrt(1.0 - (Bpol*Bpol)/(bmag*bmag))        [per ig]
  vflux_perp : g_work = -( (((((((zi*aky)*gnew)*aj1)*rhoc)*(gds21+theta0*gds22))*vperp2)*smz) / ((qval*shat)*(bmag*bmag)) )
  vflux      = vflux_par + vflux_perp
has_apar (field = aparnew):
  pmflux     : g_work = -(((gnew*aj0)*stm)*vpa)
  qmflux(.,1): g_work = g_work*energy
  qmflux(.,2): g_work = -(((((gnew*aj0)*stm)*vpa)*2.)*(vpa*vpa))
  qmflux(.,3): g_work = -((((gnew*aj0)*stm)*vpa)*vperp2)
  vmflux     : g_work = -((((gnew*aj0)*stm)*vpa)*vpac)
has_bpar (field = bparnew):
  pbflux     : g_work = (((gnew*aj1)*2.0)*vperp2)*tz
  qbflux(.,1): g_work = g_work*energy
  qbflux(.,2): g_work = ((((gnew*aj1)*2.0)*vperp2)*tz)*2.*(vpa*vpa)      # ((…*2.)*vpa**2)
  qbflux(.,3): g_work = ((((gnew*aj1)*2.0)*vperp2)*tz)*vperp2
  vbflux     : g_work = ((((gnew*aj1)*2.0)*vperp2)*tz)*vpac
```

Quoted Fortran for the subtle term (:613-614):

```fortran
g_work(:,isgn,iglo) = -zi*aky(ik)*gnew(:,isgn,iglo)*aj1(:,iglo) &
     *rhoc*(gds21+theta0(it,ik)*gds22)*vperp2(:,iglo)*spec(is)%smz/(qval*shat*bmag**2)
```

Here `zi*aky` = (0, aky). `(0,a)*(gr,gi)` = (−a·gi, a·gr) up to the sign of zero.
`vpac(ntgrid,:,:)=0`, so the `vflux_par` and `vmflux`/`vbflux` integrands vanish at the last θ point.

#### 3.6.4 `flux_vs_e` / `get_flux_vs_e` (:742-798, 841-889) → `es_flux_vs_e` etc.

The `g_work` definitions are as for `pflux`, `pmflux`, `pbflux`. Then:

```
wgt = 1.0/trapz(theta, grho*jacob); total[ie,is] = 0
for iglo in layout order:
    fac = 1.0 if aky(ik)==0 else 0.5
    total[ie,is] += fac*trapz(theta, (Im(g_work[:,1,iglo]*conj(fld[:,it,ik]))*wl[:,il])*jacob)*aky(ik)*wgt
    total[ie,is] += fac*trapz(theta, (Im(g_work[:,2,iglo]*conj(fld[:,it,ik]))*wl[:,il])*jacob)*aky(ik)*wgt
flx = 0.5*total            # (double 1/2 as in the source comment); then *= dens in common_calculate_fluxes
```

Each term is `((fac*T)*aky)*wgt`. For fixed `(ie,is)` the order over `(il,it,ik)` follows the layout string, with the
dimension of highest `*_ord` outermost:

| layout | order |
|---|---|
| `lxyes`, `lexys` | `ik → it → il` |
| `xyles` | `il → ik → it` |
| `yxles`, `yxels` | `il → it → ik` |
| `lyxes` | `it → ik → il` |

Variable `es_flux_vs_e(t, species, egrid)` (Fortran dims `[egrid, species, t]`), long_name "Electrostatic flux".
The apar and bpar variants have long_names "Perpendicular magnetic flux, as function of energy" and
"Parallel magnetic flux, as function of energy".

#### 3.6.5 Energy exchange `eexchange` (:1063-1211), with `gnew` = h_new, `g` = h_old

```
total = total2 = 0
if nonlin: exchange_dummy(1,1,:) = 0; exchange(1,1,:) = 0
for iglo (it,ik,il,ie,is):
   if nonlin and it==1 and ik==1: g_work(:,:,iglo)=0; cycle
   for isgn in 1,2:
      for ig: g_work(ig) = 0 if forbid(ig,il) else aj0(ig)*( (((zi*wdrift_func(ig,iglo))/code_dt)*gnew(ig))*tz )
      for ig: if abs(hypervisc_filter(ig,it,ik)-1.0) > epsilon(0.0):
                  g_work(ig) = g_work(ig) - ((1.0-1./hypervisc_filter(ig,it,ik))*gnew(ig))/code_dt
      for ig in -ntgrid..ntgrid-1:            # reads g_work(ig+1) BEFORE it is updated → no recurrence
          g_work(ig) = 0.5*(g_work(ig)+g_work(ig+1))
                     + (((0.5*vpac(ig))*(gradpar(ig)+gradpar(ig+1)))/delthet(ig))*(gnew(ig+1)-gnew(ig))*stm
      # g_work(ntgrid) keeps its grid-point value
total = integrate_moment(g_work)
exchange_dummy[it,ik,is] = flavg(Re(total*conj(phinew)))           # computed, never used
g_work = (aj0*0.25)*(gnew + g) ; total  = integrate_moment(g_work)
g_work = (aj0*0.25)*(gnew - g) ; total2 = integrate_moment(g_work)
for is: for ik: for it: (skip (1,1) if nonlin)
   exchange[it,ik,is] = flavg( Re( total*conj(phinew-phi) - (phinew+phi)*conj(total2) ) ) / code_dt
```

`wdrift_func` (`src/dist_fn.fpp:1500-1530`) is evaluated with the current `code_dt` and then divided by `code_dt`.
Do **not** cancel the dt algebraically:

```
vp2 = vperp2(ig,iglo); vpa2 = vpa(ig,1,iglo)**2           # isgn=1 always
if aky(ik)==0: w = akx(it)/shat*(cvdrift0*vpa2 + gbdrift0*0.5*vp2)*code_dt/2.0
else:          w = ((cvdrift+theta0*cvdrift0)*vpa2 + (gbdrift+theta0*gbdrift0)*0.5*vp2)*code_dt*wunits(ik)
```

Exact left-to-right reading of the `aky==0` branch: `((akx/shat)*(cvdrift0*vpa2 + (gbdrift0*0.5)*vp2))*code_dt/2.0`.

#### 3.6.6 Linear diffusive estimates `write_diffusive_estimates` (:334-430), only when `.not. nonlin`

```
j = igomega ; D = 0; Q = 0; P = 0
for ik: for it:
   if kperp2(j,it,ik) == 0.0: cycle
   D[it,ik] = (max(Im(omega_average[it,ik]), 0.0)/kperp2(j,it,ik))*2.0
   for is: Q[it,ik,is] = (((D*dens)*temp)*tprim)*grho(j)
           P[it,ik,is] = ((D*(dens*dens))*fprim)*grho(j)          # dens**2.0
diffusivity = maxval(D)*grho(j)
standard-properties('heat_flux_diff', Q) ; standard-properties('part_flux_diff', P)
heat_flux_max_by_spec[is] = maxval_ik( maxval_it(Q) ) ; heat_flux_max = sum_seq_is(...)
(same for P)
if write_fluxes: es_heat_flux_diff_max(t,species), heat_flux_diff_max(t), es_particle_flux_diff_max(t,species),
                 particle_flux_diff_max(t), diffusivity(t)
```

`omega_average` is the time-averaged complex frequency from Section 3.3. In record 1 it is 0, so all these are 0.

#### 3.6.7 `calculate_standard_flux_properties(name, F[it,ik,is])` (:434-531)

```
flux_by_species[is] = average_all_real_xys(F)                     # sum_it( sum_ik F*fac )
total_by_mode[it,ik] = 0 + F[it,ik,1] + F[it,ik,2] + ...           # sequential over is
total_by_ky[ik] = average_kx_xy(total_by_mode)
total_by_kx[it] = average_ky_xy(total_by_mode)
if write_fluxes:
   name              (t, species)       = flux_by_species        long_name: desc//" averaged over kx and ky"
   "total_"//name//"_by_ky" (t, ky)     = total_by_ky            desc//" summed over species and averaged over kx"
   "total_"//name//"_by_kx" (t, kx)     = total_by_kx            desc//" summed over species and averaged over ky"
   "total_"//name     (t)               = sum_seq_it(total_by_kx) desc//" summed over species and averaged over kx and ky"
   if write_fluxes_by_mode:
      name//"_by_mode"       (t, species, ky, kx) = F            desc//" as a function of species, kx and ky"
      "total_"//name//"_by_mode" (t, ky, kx)      = total_by_mode desc//" summed over species"
current_results bookkeeping:
   name == es_heat_flux       → species_es_heat_flux = flux_by_species
   name == es_energy_exchange → species_energy_exchange = ...
   name == apar_heat_flux / bpar_heat_flux → species_apar/bpar_heat_flux
   name ∈ {es,apar,bpar}_heat_flux → total_heat_flux += sum_seq_is(flux_by_species); species_heat_flux += flux_by_species
   name ∈ {es,apar,bpar}_mom_flux  → total_momentum_flux += ... ; species_momentum_flux += ...
   name ∈ {es,apar,bpar}_part_flux → total_particle_flux += ... ; species_particle_flux += ...
```

So `heat_flux_tot = ((0 + Σ_s es) + Σ_s apar) + Σ_s bpar`, accumulated in the order of the calls in 3.6.1.
It includes only the `*_heat_flux` names, not `_par` or `_perp`. The same holds for the momentum and particle totals.

The (description, units) strings are verbatim from `diagnostics_fluxes.fpp:127-208, 394-400`, and some have a
trailing "= ". For example, `es_heat_flux` is ('Turbulent flux of heat', 'Q_gB = ').
Linear-only names: `heat_flux_diff`, `part_flux_diff`.

### 3.7 Moments (`write_moments`, default ON) — `src/diagnostics/diagnostics_moments.fpp:19-330`

`getmoms(phinew, bparnew)` (:183-330):

```
g_adjust(gnew, phinew, bparnew, from)
g_work = aj0*gnew                          → density = integrate_moment
g_work = vpa*g_work                        → upar
g_work = (2.*vpa)*g_work                   → tpar ; tpar = tpar - density
g_work = (vperp2*gnew)*aj0                 → tperp ; tperp = tperp - density
g_work = ((vpa*gnew)*aj0)*energy           → qparflux
g_work = (((gnew*aj1)*2.0)*vperp2)*tz      → pperpj1
g_work = g_work*energy                     → qpperpj1
g_work = (aj0*gnew) - ((phinew*zt)*nonmaxw_corr)   → ntot
for is: ntot*=dens; density*=dens; upar*=stm; tpar*=temp; tperp*=temp;
        qparflux=((q*dens)*temp)*stm; pperpj1=(p*dens)*tz; qpperpj1=((q*dens)*temp)*tz
g_adjust(gnew, phinew, bparnew, to)
```

Then, for each moment `M` in the order `ntot, density, upar, tpar, tperp, qparflux, pperpj1, qpperpj1`,
`write_standard_moment_properties` (:83-180) writes:

```
m2[it,ik,is] = average_theta_txys(M, M)        = flavg(Re(conj(M)*M))
M2_by_ky (t, species, ky)  = average_kx_xys(m2)
M2_by_kx (t, species, kx)  = average_ky_xys(m2)
M2       (t, species)      = sum_seq_it(M2_by_kx)
mbm[it,ik,is] = flavg(M[:,it,ik,is])  (complex)
M_flxsurf_avg (t, species, kx, ri) = mbm[:, ik=1, :]           # "at ky=0 (actually ik==1)"
M         (species, ky, kx, theta, ri) = M                      # overwritten each record (no t)
M_igomega_by_mode (t, species, ky, kx, ri) = M(igomega,:,:,:)
M_by_mode (t, species, ky, kx, ri) = mbm
M2_by_mode (t, species, ky, kx)    = m2
M_t (t, species, ky, kx, theta, ri) = M   if write_{ntot,density,upar,tperp}_over_time (only these four)
```

(description, units): ntot ('The total perturbed species density', 'n_ref'); density ('The non-adiabatic part of the
perturbed species density', 'n_ref'); upar ('The perturbed parallel flow', 'v_ths'); tpar ('…parallel temperature',
'T_ref'); tperp ('…perpendicular temperature', 'T_ref'); qparflux ('The parallel heat flux', 'n_ref T_ref v_thr');
pperpj1 (…, 'n_ref T_ref / q_ref'); qpperpj1 (…, 'n_ref T_ref^2 / q_ref'). Units on squared variables are
`"("//units//")^2"`, except `M2_by_mode`, which carries the unsquared units string as in the source.

### 3.8 Parallel correlation `phi_corr_2pi` (`write_correlation`, default ON) — `src/gs2_diagnostics.fpp:1961-1976, 3389-3412`

```
cf[ig,ik] = 0
for ik: fac = 1.0 if ik==1 else 0.5          # NB index test, not aky==0
  for it: for ig: cf[ig,ik] = cf[ig,ik] + (phinew(0,it,ik)*conj(phinew(ig,it,ik)))*fac   # theta index 0 hard-coded
if has_phi: write "phi_corr_2pi"(t, ky, theta, ri)
```

### 3.9 Eigenfunction normalisation (`write_eigenfunc`, default OFF) — `src/diagnostics/diagnostics_fields.f90:57-123`, `gs2_io.fpp:560-603`

`get_phi0()` uses the **old-time-level** `phi` and `apar`:

```
phi0 = phi(0,:,:)                                            # index 0 hard-coded, not igomega
if def_parity and fapar > 0 and not even: phi0 = apar(0,:,:)
where abs(phi0) < 10*epsilon(0.0): phi0 = phi(1,:,:)/(theta(1)-theta(0))
where abs(phi0) < 10*epsilon(0.0): phi0 = 1.0
```

`nc_eigenfunc` writes `phase(ky, kx, ri) = phi0`. For each present field it writes
`X_norm(ky, kx, theta, ri) = X(ig,:,:)/phi0` (old fields, complex division). These variables are overwritten on every
write step, and on the final write as well. The ascii `.eigenfunc` file holds θ, θ0, ky and the normalised fields
for `|ig| ≤ ntheta/2 + (nperiod-1)*ntheta`.

### 3.10 ASCII `.out` and screen lines (defaults: `write_ascii=T`, `write_line=T`)

* `write_line`/`print_line` (`src/diagnostics/diagnostics_printout.f90:43-99`) produce one line per mode, with formats
  given verbatim in the source. `phtot(it,ik) = 0.5/pi*trapz(theta, abs(phinew)**2+abs(aparnew)**2+abs(bparnew)**2)`
  (no jacobian). The lines show `omegahist` and `omega_average` times `woutunits`.
* `print_flux_line`/`write_flux_line` (:123-225) use `flux_format`
  `"('t= ', e17.10, ' <', a, '**2>= ', e13.6, ' ', a, ':', *(1x, e13.6))"`.
* These files have no feedback on the run. Porting them is optional (Section 5).

### 3.11 Final-time routines (default OFF) — `src/diagnostics/diagnostics_final_routines.f90`

* `write_final_fields`: same `phi/apar/bpar` variables as 3.4 (new fields). `.fields` ascii uses OLD fields.
* `write_final_epar` (:157-230) writes `epar(ky, kx, theta, ri)` for `ig < ntgrid`:

  ```
  phi_m  = 0.5*(phi(ig+1)-phi(ig)+phinew(ig+1)-phinew(ig))*fphi
  apar_m = 0.5*(aparnew(ig+1)+aparnew(ig)-apar(ig+1)-apar(ig))*fapar
  epar   = -phi_m/delthet(ig)*0.5*(abs(gradpar(ig))+abs(gradpar(ig+1))) - apar_m/tunits(ik)/code_dt
  ```

  `epar(ntgrid,:,:)` is left uninitialised in Fortran; mask it out in comparisons.
* `write_final_moments` (:292-415): `getmoms` (one more round trip) → `nc_final_moments` writes `ntot, density, upar,
  tpar, tperp, qparflux, pperpj1, qpperpj1 (species, ky, kx, theta, ri)`. These are the same names as 3.7, overwritten
  with identical values. Ascii `.moments`, `.mom2` and `.amoments` are normalised by `phi0` only if `write_eigenfunc`.
* `write_final_antot`: `getan(antot, antota, antotp)` from the fields module → `antot/antota/antotp(ky,kx,theta,ri)`.

### 3.12 Velocity-space error diagnostics (brief; GATED)

* `write_verr` → `write_velocity_space_checks` (`src/diagnostics/diagnostics_velocity_space.fpp:50-128`) calls:
  * `get_verr` (:355-511). It round-trips `gnew`, integrates `J0 h` (and `J0 v∥ h` for apar) with the standard
    weights and with less accurate weights, dropping one energy or pitch-angle point at a time (`eint_error`,
    `lint_error`, `trap_error`, weights from `le_grids:init_weights`), and computes `errest(5,2)` and `erridx(5,3)`
    with `estimate_error`. It then calls **`adjust_vnmult(errest, …)`, which changes the collision frequency
    multiplier when `collisions:vary_vnew`**.
  * `get_gtran` (:586-707): fraction of free energy in high-order Legendre coefficients.
* netCDF (write steps only): `vspace_lpcfrac(t, 3)`, `vspace_err(t, 2, 5)`, `vspace_vnewk(t, 2)`,
  `vspace_err_maxindex(t, 3, 5)` (if `write_max_verr`). ascii: `.lpc`, `.vres`, `.vres2`.
* On non-write steps with `mod(istep,ncheck)==0` it runs again when `vary_vnew` (`gs2_diagnostics_new.f90:509-514`).
* `write_cerr` → `.cres` ascii only. It needs the Lorentz operator and `use_le_layout=F`.
* These feed back into the dynamics (`vnmult`) and are expensive. GATE them for phase 1.

### 3.13 Nonlinear convergence (`use_nonlin_convergence`, default OFF) — `src/diagnostics/diagnostics_nonlinear_convergence.f90:49-117`

This is pure scalar logic on `total_heat_flux`, run on write steps, and can set `exit`. Port it verbatim if supported.
Note the integer divisions `conv_nstep_av/nwrite` and `conv_nsteps_converged/nwrite`, and the ring buffer
`conv_heat(0:conv_nstep_av/nwrite-1)`.

---

## 4. netCDF file contract (`<run>.out.nc`)

The format is netCDF-4 (HDF5), chunked, with no deflate applied (verified with `ncdump -hs`). Every variable is
`double` except integer coordinates and counts, and the `char` variables `layout` and `input_file`.
Complex data use a trailing `ri` dimension (size 2: real, imaginary).

### 4.1 Static part (written at creation): recommended to be produced by the Fortran export program

| routine | content |
|---|---|
| `create_metadata` (`src/gs2_metadata.fpp:20-57`) | global attrs `title="GS2 Simulation Data"`, `software_name`, `software_version`, `netcdf_version`, `date_created`, `id` (uuid), `Conventions`, `references`, `comment` |
| `define_dims` (`src/gs2_io.fpp:276-364`) | dims+coords `kx(akx)`, `ky(aky)`, `theta`, `lambda(al)`, `species(1..n)`, `egrid(1..n)`, `energy(species, egrid)`, `sign`, `t` (unlimited), `ri`, `nheat(7)`, `vpa(negrid*nlambda)`, `char5`, `char10`, `char200`, `nlines`, `2`,`3`,`4`,`5`, `bset`; scalars `layout`, `nproc`, `ntheta_tot`, `nkx`, `nky`, `nspecies`, `nenergy`, `nlambda`, `nperiod`, `ntheta`, `nbset`, `theta0(ky, kx)`, `nmesh` |
| `nc_norms` (:368-396) | `mref, zref, nref, tref, aref, vref, bref, rhoref` |
| `nc_species` (:1601-1653) | `charge, mass, dens, temp, tprim, fprim, uprim, uprim2, vnewk, type_of_species, type, f0type, nu_h, bess_fac, vcrit, vcprim, dens0, u0, tpar0, tperp0 (species)` |
| `nc_geo` (:1656-1722) | `bmag, bpol, kxfac, gradpar, gbdrift, gbdrift0, cvdrift, cvdrift0, gds2, gds21, gds22, gds23, gds24, gds24_noq, grho, jacob, shat, drhodpsi, eps_trapped, cdrift, cdrift0, qval, beta, rplot, zplot, aplot, rprime, zprime, aprime, kperp2(ky, kx, theta)` |
| `save_input` (:453-480) | `input_file(nlines, char200)` |

`kx` is stored in GS2's FFT order (0, +, …, −), not sorted. `nproc` is the Fortran run's value. Recommendation: the
export program calls these same routines on a fresh `<run>.out.nc`. JAX opens it in `"a"` mode and appends the
time-dependent and "final" variables with attributes as listed. This guarantees byte-identical static content.

### 4.2 Time-dependent and final variables written by the supported configuration

`nout` is 1-based in Fortran; in JAX the record index is `nout-1`. The netCDF dims below are in C order.

| variable | netCDF dims | condition | source (Section) |
|---|---|---|---|
| `t` | (t) | always | 3.2 |
| `omega`, `omegaavg`, `omega_average` | (t, ky, kx, ri) | write_omega | 3.3 |
| `X2_by_kx` | (t, kx) | write_fields & has_X | 3.4 |
| `X2_by_ky` | (t, ky) | ″ | 3.4 |
| `X0`, `X_igomega_by_mode` | (t, ky, kx, ri) | ″ | 3.4 |
| `X2_by_mode` | (t, ky, kx) | ″ | 3.4 |
| `X2` | (t) | ″ | 3.4 |
| `X` (final, overwritten) | (ky, kx, theta, ri) | write_fields (all present X) or write_final_fields | 3.4 |
| `X_t` | (t, ky, kx, theta, ri) | write_fields & write_X_over_time | 3.4 |
| `ql_metric_by_mode` | (t, ky, kx) | write_ql_metric | 3.5 |
| `F`, `total_F_by_ky`, `total_F_by_kx`, `total_F` | (t,species), (t,ky), (t,kx), (t) | calculate_fluxes & write_fluxes | 3.6.7 |
| `F_by_mode`, `total_F_by_mode` | (t,species,ky,kx), (t,ky,kx) | … & write_fluxes_by_mode | 3.6.7 |
| `es_flux_vs_e`, `apar_flux_vs_e`, `bpar_flux_vs_e` | (t, species, egrid) | calculate_fluxes & has_* (not gated on write_fluxes) | 3.6.4 |
| `heat_flux_tot`, `hflux_tot`, `mom_flux_tot`, `vflux_tot`, `part_flux_tot`, `zflux_tot` | (t) | calculate_fluxes & write_fluxes | 3.6.1 |
| `es_heat_flux_diff_max`, `es_particle_flux_diff_max` | (t, species) | linear & calculate_fluxes & write_fluxes | 3.6.6 |
| `heat_flux_diff_max`, `particle_flux_diff_max`, `diffusivity` | (t) | ″ | 3.6.6 |
| moment family for `M` | see 3.7 | write_moments | 3.7 |
| `phi_corr_2pi` | (t, ky, theta, ri) | write_correlation & has_phi | 3.8 |
| `phase`; `phi_norm`, `apar_norm`, `bpar_norm` | (ky,kx,ri); (ky,kx,theta,ri) | write_eigenfunc | 3.9 |
| `epar` | (ky, kx, theta, ri) | write_final_epar (end) | 3.11 |
| `antot`, `antota`, `antotp` | (ky, kx, theta, ri) | write_final_antot (end) | 3.11 |

In the flux rows, `F` ranges over:

* has_phi: `es_heat_flux, es_heat_flux_par, es_heat_flux_perp, es_part_flux, es_part_tormom_flux, es_mom_flux, es_mom_flux_par, es_energy_exchange, es_mom_flux_perp`;
* has_apar: `apar_heat_flux, apar_heat_flux_par, apar_heat_flux_perp, apar_part_flux, apar_mom_flux`;
* has_bpar: the five `bpar_*` names;
* linear only: `heat_flux_diff, part_flux_diff`.

Verified by running `bin/gs2` on a linear electrostatic single-mode case and a nonlinear EM two-species box case.
The variable sets in `ncdump -h` match this table exactly.

**Old-diagnostics names** (`use_old_diagnostics=T`, gated) that some older scripts read, such as `es_heat_by_k`,
`es_part_by_k`, `es_mom_by_k`, `es_heat_par`, `apar_heat_by_x`, `es_exchange`, `phtot`, `ntot2`, `phi00`, `ntot00`,
are **not** produced by the new diagnostics. GS2 8.2.1 does not produce them by default either, so the JAX code
should not produce them. The new diagnostics keep `hflux_tot`, `vflux_tot`, `zflux_tot`, `omegaavg` and `phi0`
for backwards compatibility.

---

## 5. Options that change behaviour here

Namelist `gs2_diagnostics_knobs` unless stated otherwise. Defaults are from
`src/diagnostics_base_configuration.f90:13-536` as modified by `apply_new_defaults` and the nonlinear smart defaults.

| option | default | effect | decision | reason |
|---|---|---|---|---|
| write_any | T | master switch (no file, no omega or exit tests, no round trips) | SUPPORT | trivial |
| nwrite | 10 | write period | SUPPORT | core |
| navg | 10 | omega ring-buffer length | SUPPORT | core |
| omegatol | 1e-3 | convergence tolerance; ≤0 disables; tiny-threshold factor | SUPPORT | core |
| omegatinst | 1e6 | instability abort threshold | SUPPORT | core |
| exit_when_converged | T | omega/nonlinear convergence sets exit | SUPPORT | core |
| igomega | 0 | θ index for omega, X0, M_igomega, diffusive estimates | SUPPORT | core |
| write_omega | T (forced T by write_line) | omega output | SUPPORT | core |
| write_fields | T (forced) | phi2 family + final fields | SUPPORT | core |
| write_phi/apar/bpar_over_time | F | X_t | SUPPORT | cheap |
| write_ql_metric | T | ql_metric_by_mode | SUPPORT | default on |
| write_fluxes | F (T if nonlinear) | flux family; triggers calculate_fluxes | SUPPORT | core |
| write_fluxes_by_mode | F (T if nonlinear) | *_by_mode | SUPPORT | core |
| print_flux_line / write_flux_line | F / F | triggers calculate_fluxes (side effects!) + text | SUPPORT | standard |
| print_line / write_line | F / T | text; force write_fields/write_omega | SUPPORT | defaults; text optional |
| write_ascii | T | text files | SUPPORT | text optional, never refuse |
| write_moments | T | moment family (+round trip) | SUPPORT | default on |
| write_ntot/density/upar/tperp_over_time | F | M_t | SUPPORT | cheap |
| write_correlation | T | phi_corr_2pi | SUPPORT | default on |
| write_eigenfunc | F | phase, *_norm, .eigenfunc | SUPPORT | standard linear |
| write_final_fields | F | final phi/apar/bpar (+.fields) | SUPPORT | trivial |
| write_final_moments | F | final moments | SUPPORT | reuses getmoms |
| write_final_epar | F | epar | SUPPORT | trivial |
| write_final_antot | F | antot/antota/antotp | SUPPORT | needs getan (fields spec) |
| use_nonlin_convergence (+conv_*) | F | can end nonlinear runs | SUPPORT | small scalar logic |
| nc_sync_freq | 10 | I/O flush only | SUPPORT | no-op semantics |
| nwrite_mult | 10 | only gated diagnostics | SUPPORT | no effect when they are gated |
| enable_parallel, serial_netcdf4 | F | ignored (`parallel_io_capable=.false.`) | SUPPORT | no-op |
| write_avg_moments, write_omavg, ob_midplane, dump_check1/2, nmovie, write_pflux_sym, write_pflux_tormom | — | non-functional in new diagnostics (warnings only) | SUPPORT | no-op |
| save_for_restart, nsave, save_many, file_safety_check, save_glo_info_and_grids, save_velocities | F, -1, … | restart files (not .out.nc) | SUPPORT | per restart spec; warn and continue acceptable |
| save_distfn | F | dist-fn dump at end (with g_adjust) | GATE | not needed |
| append_old | F | append to existing file | GATE | record-index bookkeeping; rarely used |
| make_movie | F | real-space fields | GATE | needs transforms; heavy |
| dump_fields_periodically | F | dump.fields.t= files | GATE | ascii, rarely used |
| write_heating | F | heating diagnostics | GATE | needs collision heating accumulators |
| write_collisional | F | collisional heating | GATE | same |
| write_verr / write_max_verr / ncheck | F / F / 10 | v-space errors; may change vnmult | GATE | feedback, expensive |
| write_cerr | F | .cres | GATE | Lorentz-only, ascii |
| write_symmetry | F | momentum flux vs θ,v∥ | GATE | non-standard |
| write_nl_flux_dist | F | poloidal flux distributions (+round trip) | GATE | non-standard |
| write_parity | F | .parity (+round trip) | GATE | non-standard |
| write_cross_phase | F | electron n–T⊥ cross phase (+round trip) | GATE | non-standard |
| write_correlation_extend | F | extended-domain correlation | GATE | large memory, index juggling |
| write_jext / write_lorentzian | F | antenna | GATE | antenna gated |
| write_zonal_transfer / write_kinetic_energy_transfer | F | transfer diagnostics | GATE | non-standard |
| write_full_moments_notgc | F | non-guiding-centre moments | GATE | non-standard |
| write_g / write_gyx / write_gs / write_kpar / write_final_db | F | ascii-only | GATE | ascii-only |
| knobs:use_old_diagnostics | F | old module | GATE | different outputs |
| knobs:wstar_units | F | tunits/woutunits per ky, t_wstar | GATE | changes time normalisation |
| knobs:eqzip_option (ieqzip) | none | re-enables omega in nonlinear runs | GATE | rarely used |
| knobs:fphi/fapar/fbpar | 1/0/0 | has_*; g_adjust factors | SUPPORT | core |
| layouts:layout | 'lxyes' | summation order | SUPPORT | order parameter |
| collisions:vary_vnew | F | runtime vnmult changes | GATE | feedback |
| collisions:heating | F | forced off unless write_heating/write_collisional | SUPPORT | as rule |
| dist_fn:def_parity / even | F / T | get_phi0 | SUPPORT | trivial |
| nonlinear_terms:nonlinear_mode | off | defaults, omega skip, exchange (1,1) skip, diffusive estimates | SUPPORT | core |
| knobs:nstep, max_sim_time, avail_cpu_time | — | exit triggers → final forced write | SUPPORT | core |

---

## 6. GPU mapping notes

* **Per-step work** (linear runs, or nonlinear with any(ieqzip)): `omega_inst[ik,it]` and the ring-buffer update are
  elementwise over `(ik,it)`. The average over `navg` and the convergence sum over `navg` are tiny sequential
  sums (unrolled Python loop or `lax.fori_loop`). Run them inside the step function. The exit decision is a boolean
  scalar: use `lax.while_loop`, or run fixed-size chunks with a host-side check of a per-step `exit_step` record,
  then truncate.
* **Write-step work**: run a diagnostics function only when `write_step(istep)`. With `lax.cond` inside a scan, or by
  splitting the scan into chunks of `nwrite` steps on the host, make sure the g_adjust round trips modify the carried
  `gnew` only on write steps. Recommended host loop: `for chunk: state = scan(step, state, nwrite); state, rec = diag(state); write rec`.
  The final forced record and exit-triggered records need host logic.
* **Batching**: everything below is independent over `(is, ik, it)` and can be broadcast or `vmap`-ed:
  * `g_adjust` (elementwise with mask `trapped & il > jend(ig)`);
  * all `g_work` constructions (elementwise);
  * the `trapz` integrands.
* **Sequential recurrences** (`lax.scan`, carrying state vectorised over the batch axes):
  1. Klein sums along θ (length `2*ntgrid`), used by every `trapz`/`flavg`;
  2. `integrate_moment` accumulation over `(ie, il)` in layout order (length `negrid*nlambda`);
  3. `get_flux_vs_e` accumulation over `(il,it,ik)` in layout order (length `nlambda*ntheta0*naky`, 2 adds per element);
  4. the small `sum_seq` reductions over `it`, `ik`, `is` and the accumulators in `average_ky/kx`;
  5. `correlation` accumulation over `it`.
  Implement one generic `seq_sum(x, axis)` and one `klein_sum(x, axis)` and reuse them. Never use `jnp.sum`,
  `jnp.cumsum` (associative scan) or `jnp.trapz` where round-off faithfulness matters.
* **Masks**: forbidden region (`forbid[il,ig]` for the exchange; `trapped & il>jend[ig]` for g_adjust); `aky==0`
  (fac, ql zeroing, wdrift branch, exact float compare); `nonlin & (it,ik)==(1,1)` (exchange); `kperp2(igomega)==0`
  (diffusive); `hypervisc_filter ≠ 1` (exchange); tiny and epsilon thresholds (omega, ql). No linked chains or
  ragged structures appear in the supported diagnostics. The only ragged one, the extended-domain correlation, is gated.
* `eexchange` cell-centring is a shifted elementwise op (`0.5*(a[:-1]+a[1:]) + …`, last element kept). It is not a
  recurrence.
* Memory: `integrate_moment` outputs are `(nspec, naky, ntheta0, nth)` complex. A write step needs about 20 of them
  (9 flux moments + 3 exchange + 8 moments). Compute them in sequence and avoid materialising extra copies of `g`.

---

## 7. Pitfalls for round-off agreement

* **P1 — Diagnostics mutate `gnew`.** On every write step, including istep 0, `gnew ← (gnew + adj) − adj` once for
  `calculate_fluxes` (if on) and once for `write_moments` (if on), in that order, with `adj` from
  `phinew, bparnew`. Skipping this gives about 1-ulp differences in `gnew` that grow in nonlinear runs. If the
  Fortran export program dumps state *after* `initialize_diagnostics`, the istep-0 round trips are already applied
  and must not be applied again. If it dumps *before*, JAX must perform the istep-0 diagnostics itself. Choose one
  and document it in the export spec. Recommendation: dump before, and let JAX produce record 1.
* **P2 — Omega every step.** The ring-buffer slot is `mod(istep,navg)`. The convergence test (`istep > navg`) and the
  instability test run every step, and convergence forces a write on that step. The final forced call recomputes
  omega for the same `istep`, an idempotent no-op in value.
* **P3 — Old and new time levels.** omega uses `phi` (old) and `phinew`. The eigenfunction `phase` and `*_norm` and
  `get_phi0` use OLD `phi`. `phi`, `phi2` and `phi0` in the file are NEW. The exchange uses both g levels and both φ
  levels. `epar` uses both.
* **P4 — Time.** `t` is the accumulated `user_time`, not `istep*delt`. Accumulate in the same way (sequential float64
  adds of `user_dt`, which may change after resets).
* **P5 — Summation orders.** Klein sums for every θ integral. Sequential `SUM`. `iglo`-ordered accumulation in
  `integrate_moment` and `get_flux_vs_e`, which depends on `layout`. MPI reductions are rank-order dependent, so
  reference runs must use nproc = 1. OpenMP `REDUCTION` clauses (`le_grids.f90:1363-1379`) reorder sums if OpenMP is
  enabled; this build has no OpenMP.
* **P6 — Mode weights.** `fac = 0.5/1.0` from `aky == 0.` in the volume averages and `flux_vs_e`, but from `ik == 1`
  in `correlation`.
* **P7 — Thresholds.** `tiny(0.0)*(1000/|omegatol|)` in `calculate_omega` versus `(tiny(0.0)*1000)/|omegatol|` in
  `calculate_instantaneous_omega`. `epsilon(0.0)` on `max|field|²` in the ql metric. `10*epsilon` in `get_phi0`.
  `epsilon` on `|filter−1|` in the exchange. Branch decisions near a threshold can flip with 1-ulp input differences,
  for example the omega convergence step. When they flip, report the step instead of failing the numeric comparison.
* **P8 — Complex intrinsics.** Complex division (gfortran Smith-style with `-fcx-fortran-rules` versus XLA), `log`
  (macOS libm `clog` versus XLA `log|z| + i·atan2`) and `abs` (hypot) differ at ulp level. The principal branch and
  signed zeros agree with C99 semantics.
* **P9 — Duplicate and quirky outputs.** `omegaavg ≡ omega_average`; `X0 ≡ X_igomega_by_mode`; `hflux_tot ≡ heat_flux_tot`;
  `vflux_tot ≡ mom_flux_tot`; `zflux_tot ≡ part_flux_tot`. `es_part_tormom_flux*` ≡ 0. Fields `phi/apar/bpar` are
  overwritten each record. Moment `M` variables have no `t` dimension. `es_flux_vs_e` is written even without
  `write_fluxes`. `es_mom_flux ≠ es_mom_flux_par + es_mom_flux_perp` at round-off level.
* **P10 — FMA.** gfortran -O3 on aarch64 contracts multiply-adds; FMA instructions are present in `diagnostics_fluxes.o`,
  `diagnostics_moments.o`, `le_grids.o`, `dist_fn_arrays.o`, `volume_averages.o`, `integration.o`
  (complex trapz) and `diagnostics_omega.o`. XLA GPU (NVPTX) also contracts, at different sites. Exact emulation would
  need Dekker two-product at each fused site, which is not recommended. Use relative tolerance around 1e-12, plus an
  absolute floor of about 1e-14 × the typical magnitude of the quantity. That floor matters for fluxes, which are
  differences of large cancelling θ integrals, and for record-1 quantities that are round-off-sized: `heat_flux_tot`
  at t=0 was 3.7e-19 in the test run.
* **P11 — Exchange arithmetic.** `wdrift_func` includes `code_dt` and is then divided by `code_dt`; do not simplify.
  The `hypervisc_filter` used is the one left by the current step's hyper routine (state, dt-dependent). The
  `(1,1)` mode is skipped only for nonlinear runs. `delthet(ntgrid)=0` is never used as a divisor here
  (loop stops at `ntgrid-1`).
* **P12 — `vpac(ntgrid)=0`** and `epar(ntgrid)` undefined.
* **P13 — Forced flags.** `write_line=T` (default) forces `write_fields=T` and `write_omega=T`. `print_flux_line` and
  `write_flux_line` trigger flux computation and therefore the round trip. Export the post-`check_parameters` flags.
* **P14 — Restart.** `omegahist` is not in restart files. After a restart the omega history restarts at zero and
  `istep` restarts at 0, so record 1 of a restarted run is written with zero omega. `user_time_old` starts at 0.0, so
  `species_*_flux_avg` (Trinity only) receives `flux*tstart` from the first record.
* **P15 — XLA simplifications.** Make sure XLA does not fold `(g + adj) − adj` to `g` or reassociate the Klein
  compensation. Keep fast-math off (`--xla_cpu_enable_fast_math=false` on CPU; GPU does not reassociate by default),
  and wrap the round trip in `jax.lax.optimization_barrier` if a simplification is ever observed.

---

## 8. Open questions

1. Should the reference `.out.nc` used for validation come from a run with `nproc=1`? The integration order is only
   reproducible serially. If production references must be multi-rank, agree on a looser tolerance (about 1e-14
   relative per reduction).
2. The export point relative to `initialize_diagnostics` (P1) must be fixed jointly with the export-program spec.
3. Exact FMA fusion sites in `g_adjust` and `flux`: which operand gets fused could be read from the disassembly
   (`otool -tv objs/dist_fn_arrays.o`) if bitwise agreement of the round-trip perturbation ever becomes a goal. This
   spec assumes it does not.
4. `hypervisc_filter` must be in the JAX state at diagnostic time. For the `iso_shear` and `aniso_shear` options it
   depends on `phi` at the time `hyper_diff` ran within the step. The hyper spec must say which time level that is.
5. `write_final_antot` needs `getan`. Confirm that the fields spec exposes it with identical summation order
   (`integrate_species`).
6. Whether existing user analysis scripts (outside this repo) rely on old-diagnostics names. Only the in-repo tests
   were checked; they read `t`, `theta`, `ky`, `kx`, `jacob`, `bmag`, `gbdrift`, `bpol`.
7. The `units` attribute strings contain literal oddities (`"Q_gB = "`, `"n_r? "`, `"Pi_gB =  "`). Should the JAX writer
   reproduce them verbatim? Recommended yes, copying them from the source.

---

## 9. Suggested verification

1. Single-mode linear electrostatic case (for example `tests/linear_tests/cyclone_itg/cyclone_itg_low_res_base.in`
   with `grid_option="single"`, `nstep=45`, `nwrite=10`, `write_fluxes=T`, `write_eigenfunc=T`). Compare every
   variable in `ncdump -h` at all records. This case exercises omega, the diffusive estimates, fluxes, moments,
   correlation and eigenfunctions.
2. A `range` or `box` linear case with several ky: checks `fac` handling and the all-mode convergence test.
3. A nonlinear EM two-species box case (as in `tests/nonlinear_tests/cyclone_itg_all_diagnostics/cyclone_itg_low_res_base.in`
   with default diagnostics and `nwrite=2`): checks the apar and bpar families, the exchange `(1,1)` skip and the
   omega zeros.
4. A layout sweep (`lxyes`, `yxels`, `xyles`) on case 1: checks the integrate_moment and flux_vs_e orders.
5. An ablation on case 3 with the g_adjust round trips disabled in JAX: should show measurable divergence of `gnew`
   after N steps, confirming P1 is needed.
