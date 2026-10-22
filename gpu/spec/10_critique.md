# 10 — Completeness and consistency critique (merged contract)

Status: merged contract over specs 01–09. Read-only analysis of GS2 8.2.1 (branch `gpu-port`,
commit `4d8c94bcf`); no source file was modified. All `file:line` citations are relative to
`/Users/nico/Workspace/gs2/src` unless another directory is given.

Inputs reviewed in full: `01_time_loop.md`, `02_source_term.md`, `03_invert_rhs.md`, `04_fields.md`,
`05_nonlinear.md`, `06_collisions.md`, `07_hyper_misc.md`, `08_diagnostics.md`, `09_export_infra.md`.
Wherever two specs disagreed, or a claim decided control flow or round-off, the Fortran was re-read; the
evidence is cited in §D.

Where this document and an earlier spec disagree, **this document wins** (it records the resolution
determined from the Fortran). It does not repeat formulas that are correct and unambiguous in the
owner specs; it points to them.

---------------------------------------------------------------------------------------------------

## 0. Conventions and summary of findings

### 0.1 Shape symbols (NumPy C order; netCDF files store the reversed Fortran order, spec 09 §0.2)

| symbol | meaning |
|---|---|
| `S, E, L, Y, X` | `nspec, negrid, nlambda, naky, ntheta0` |
| `T` | `nth = 2*ntgrid+1` theta points (index 0 ↔ `ig = -ntgrid`) |
| `N` | `T-1 = 2*ntgrid` theta cells |
| `G7` | `(S,E,L,Y,X,2,T)` canonical g-like array: Fortran `g(ig, isgn, it, ik, il, ie, is)` |
| `G6` | `(S,E,L,Y,X,T)`; `G5` = `(S,E,L,Y,X)` |
| `F3` | `(Y,X,T)` field-like; `K2` = `(Y,X)` |
| `LE6(NE,NX)` | `(S,Y,X,T,NE,NX)`, Fortran `(ixi, ie, ig, it, ik, is)` (06 §2.4) |
| `XI` | `nxi_lim`; `XP` = `nxi+1`; `E1` = `negrid+1`; `L2` = `2*nlambda` |
| `R_c, M_c` | `nidx*N_class(c)`, `M_class(c)` (04 §0.1) |

Complex arrays: a trailing `ri` axis of length 2 in NumPy (leading Fortran `ri` dimension), as in spec 09
§0.2 and in `<run>.out.nc` (see G-03). Logical arrays are int8. Index-valued integers keep their Fortran
1-based values and the `-1` sentinel.

### 0.2 Headline findings (details in §D)

1. **Wrong claim, behaviour-relevant (G-01).** Spec 06 says `kwork_filter` is all `.false.` while stepping
   and asks for an assertion. In box runs, `kwork_filter(1,1) = .true.` from the end of `init_fields` on
   (fields.f90:376), and every collision stage skips that mode (collisions.fpp:3014, 3067, 3280, 3309, 3333,
   3361, 3562, 3588, …). The assertion would fail in every nonlinear box run. The masks already written in 06 §4
   (`active(it,ik)`) are correct; delete the "assert all false" instruction.
2. **Wrong claim, provenance only (G-02).** Specs 02 and 05 say `aj0/aj1` come from GS2's local `spfunc`. The
   build uses `-DSPFUNC=_SPF200X_` (`bin/gs2 --build-config`; Makefile:455), so `j0 = bessel_j0` and
   `j1 = bessel_j1(x)/x`, with `0.5` only at `x` exactly 0 (externals/utils/spfunc.fpp:737-790,
   externals/utils/warning_helpers.f90:87-90). The decision is unchanged: EXPORT all Bessel values.
3. **Unneeded state removed (G-06).** Spec 08 keeps `hypervisc_filter` as diagnostic state for
   `es_energy_exchange`. The filter and `wdrift_func` enter only `exchange_dummy`, which is never written
   (diagnostics/diagnostics_fluxes.fpp:273, 1127-1128, 1156). The written `es_energy_exchange` comes from
   `exchange_sym`, which uses only `g, gnew, phi, phinew` (1169, 1181, 1199). The first half of `eexchange` can
   be skipped.
4. **Rebuild concern resolved (G-17).** Specs 01, 03, 04 and 07 say the dt-reset rebuild of the implicit
   response matrix runs with `kwork_filter(1,1) = .true.`, so its block "differs". This has no numerical
   effect. The only filtered mode is box `(1,1)`, which has `aky = akx = 0`. Its inverse rows are set to 0 at
   every pivot regardless of `am` (fields_implicit.f90:834-857), and no other member couples to it: its
   shear-rate contribution is an exact 0 because `k = 0` (hyper.f90:702, 766-768). JAX may rebuild with either
   filter state.
5. **Option conflicts resolved** (G-07 … G-16, G-24, G-25): `eqzip_option` GATE, `zip` GATE,
   `use_old_diagnostics` GATE, `split_collisions` SUPPORT, `use_order_based_error` SUPPORT, `def_parity`
   SUPPORT, `boundary_option='alternate-zero'` SUPPORT as an exact alias of `zero`,
   `source_option='phiext_full'` GATE when `phi_ext /= 0`, `vary_vnew` SUPPORT as a no-op (the gate is on
   `write_verr`), `wstar_units` GATE, `use_le_layout=F` GATE, hybrid electrons SUPPORT except with
   collisions.
6. **Export contract conflicts resolved** (G-03, G-04, G-05): complex values use the `ri` axis; the mixed-wfb
   end Bessel factors are pre-gathered at the chain ends with NumPy shape `(S,Y,X,E)`; the JAX initial state
   is `/state_t0` (after `initialize_diagnostics`).
7. **Exporter gaps** (G-19, G-20): the dt ladder covers only downward rungs, and some `get_*_config()`
   objects are reset by the level-30/31 teardown.

---------------------------------------------------------------------------------------------------

## A. Master export list

Legend. **Decision**: `EXPORT` = dump from Fortran (JAX never recomputes); `RECOMPUTE` = JAX computes it
with the formula in the owner spec; `BOTH` = exported at the initial dt for validation, and JAX ports the
formula to rebuild it when dt changes. "STATE" marks evolving state; its t = 0 value comes from
`/state_t0` (spec 09 §6.1.3). **dt**: depends on `code_dt`. "frozen" means the formula contains `code_dt`
but GS2 never recomputes the array after the first init (keep the initial-dt values forever).
Group = netCDF group of the export file (spec 09 §6.1.3).

### A.1 Evolving state (`/state_t0`; time block of spec 09 §7.6)

| name | group | NumPy shape | set at (file:line) | dt | decision | notes / owner |
|---|---|---|---|---|---|---|
| `gnew` | dist_fn_arrays | G7 c128 | ginit / set_init_fields; diagnostics round trips at istep 0 (08 §3.1.8) | state | EXPORT (STATE) | the state; JAX sets `g := gnew` at the start of each step (fields_implicit.f90:337) |
| `g` | dist_fn_arrays | G7 c128 | copy of gnew; round-tripped by istep-0 flux diagnostics | state | EXPORT (check only) | not bitwise `gnew` after istep-0 `write_moments`; never read before step 1 overwrites it |
| `phinew, aparnew, bparnew` | fields_arrays | F3 c128 | get_init_field (dist_fn.fpp:7792-7799) | state | EXPORT (STATE) | absent fields are exactly 0 |
| `phi, apar, bpar` | fields_arrays | F3 c128 | copies of `*new` | state | EXPORT (check only) | equal to `*new` at t0 |
| `gexp_1, gexp_2, gexp_3` | dist_fn_arrays | G7 c128 | allocated iff `nonlin .and. .not. split_nonlinear` (dist_fn.fpp:3970-3980); zero or restart | state | EXPORT (STATE) | not remapped by exb_shear (02 §3.8) |
| `kx_shift` / `theta0_shift` | dist_fn_arrays | (Y,) f64 | allocated iff `abs(g_exb*g_exbfac) > epsilon(0.)` (3983-3991); restart (gs2_save.fpp:820-830) | state | EXPORT (STATE) | host-side scalar recurrence (07 §3.3.2) |
| `vnmult` | collisions | (2,) | init_vnew (collisions.fpp:1723-1733), clamps (2225, 2298, 1971, 2006), set_vnmult (gs2_init.fpp:329) | state | EXPORT | constant at runtime (see G-15) |
| `code_dt, user_dt` | gs2_time | () | init_delt, save_dt | state | EXPORT (STATE) | |
| `code_dt_prev1, code_dt_prev2` | gs2_time | () | init_delt (sentinel −1.0), update_time | state | EXPORT (STATE) | set the AB order |
| `code_dt_old` | gs2_time | () | static 0 (gs2_time.f90:31), update_time | state | EXPORT (STATE) | 0 at the first step, including after a restart |
| `code_time, user_time` | gs2_time | () | init_tstart, update_time | state | EXPORT (STATE) | accumulated by repeated addition |
| `code_dt_cfl, user_dt_cfl` | gs2_time | () | 1e8; save_dt_cfl | state | EXPORT (STATE) | |
| `code_dt_min, code_dt_max` | gs2_time | () | save_dt_min(delt_minimum), save_dt_max(dt0) (gs2_reinit.f90:300-301) | no | EXPORT | |
| `reset` | run_parameters | () bool | `.false.` | state | EXPORT | |
| `exb_first` | dist_fn | () bool | `.true.` until the first shear call (dist_fn.fpp:238, 4235) | state | EXPORT | |
| `nl_istep_last`, `nl_initializing` | nonlinear_terms | () | 0 / `.false.` after nl_finish_init | state | EXPORT (09) or RECOMPUTE (05) | must survive resets |
| `kwork_filter` | kt_grids | K2 bool | `.false.` (kt_grids.f90:188-189); `(1,1)=.true.` for box (fields.f90:376) | no | EXPORT (value after init_fields) | read by source, invert_rhs_1, linked BCs, integrals, collisions (G-01) |
| SAVE counters: `timeadv::istep_last`, `exb_shear::istep_last`, `reset_time_step::{istep_last,nconsec}`, `run_diagnostics::istep_last`, `nout`, `user_time_old`, `start_time` | various | () | dist_fn.fpp:4031, 4217; gs2_reinit.f90:120-121; gs2_diagnostics_new.f90:377 | state | RECOMPUTE | initial values: 0 (−1 if the response was read from file), 0, (−1, 0), 0 (istep-0 record written), nout = 2, 0 |
| `omegahist, omega_average` | diagnostics_omega | (navg,Y,X) c, (Y,X) c | zero (diagnostics_omega.f90:26-36) | state | RECOMPUTE | still zero after istep 0 (omega runs only for `istep > 0`) |
| `trin_istep, conv_*` | diagnostics_nonlinear_convergence | | | state | — | GATED (`use_nonlin_convergence`) |

### A.2 Theta grid and geometry (`/theta_grid`)

| name | NumPy shape | source | dt | decision | read by |
|---|---|---|---|---|---|
| `theta` | (T,) | theta_grid finish_init | no | EXPORT | remove_zonal_flows, Klein trapezoids, NL none |
| `delthet` | (T,) | theta_grid.f90:510-511 (`delthet(ntgrid)=0`) | no | EXPORT | vpar rebuild (02 §3.1), epar |
| `bmag, gradpar` | (T,) | geometry | no | EXPORT | vperp2 check, vpar rebuild, collisions rebuild, field eqs, diagnostics |
| `gbdrift, gbdrift0, cvdrift, cvdrift0` | (T,) | geometry with drift knobs applied (theta_grid.f90:623-626) | no | EXPORT | wdrift rebuild |
| `cdrift, cdrift0` | (T,) | geometry | no | EXPORT | Coriolis rebuild |
| `itor_over_B` | (T,) | geometry, `btor_slab` override (dist_fn.fpp:1369) | no | EXPORT (effective) | source (g_exb drive) |
| `gds2, gds21, gds22` | (T,) | geometry | no | EXPORT | simple3D init only (frozen filter is exported); vflux_perp diag |
| `grho, jacob, Rplot, Bpol` | (T,) | `jacob = 1/(drhodpsi*gradpar*bmag)` (theta_grid.f90:516) | no | EXPORT | field-line averages, fluxes |
| `field_line_average_weight` | () | `1.0/trapezoidal_integration(theta, jacob)` (theta_grid.f90:520) | no | BOTH | remove_zonal_flows, diagnostics; recompute must match bitwise (Klein sum) |
| `shat, kxfac, qval, rhoc, drhodpsi` | () | geometry | no | EXPORT | exb branch, NL scale, fluxes |
| `ntheta, ntgrid, nperiod` | () | `ntgrid = ntheta/2+(nperiod-1)*ntheta` | no | EXPORT | exb theta0 branch (`ileft/iright`), eigenfunc |
| `is_effectively_zero_shear()` | () bool | theta_grid.f90 | no | EXPORT | boundary default |

### A.3 Wavenumber grids (`/kt_grids`)

| name | NumPy shape | source | dt | decision | notes |
|---|---|---|---|---|---|
| `aky` | (Y,) | kt_grids_box.f90:378-380 / range | no | EXPORT | exact `== 0.0` tests everywhere (07 P14) |
| `akx`, `ikx` | (X,), (X,) int | kt_grids_box.f90:382-389 | no | EXPORT | kx order 0,+,−; FFT position `J(it)` (05 §1) |
| `theta0` | K2 | kt_grids_box.f90:391-401 | no | EXPORT | wdrift rebuild, exb `dtheta0`, vflux_perp |
| `kperp2, inv_kperp2` | F3 | init_kperp2 (kt_grids.f90:313-349) | no | EXPORT | collisions rebuild, fields direct, ql metric |
| `nx, ny, box, reality, gridopt_switch` | () | get_sizes (kt_grids.f90:269-291) | no | EXPORT | |
| `kwork_filter` | K2 bool | see A.1 | no | EXPORT | |

### A.4 Velocity grids (`/le_grids`, `/species`)

| name | NumPy shape | source | dt | decision | notes |
|---|---|---|---|---|---|
| `energy, w` | (S,E) | set_vgrid (le_grids.f90:285-296) | no | EXPORT | |
| `energy_maxw, speed_maxw, w_maxw` | (E,) | le_grids.f90:311-313 | no | EXPORT | collisions only |
| `al` | (L,) | lgridset | no | EXPORT | |
| `wl` | (L,T) | lgridset; 0 where forbidden (le_grids.f90:2129) | no | EXPORT | |
| `forbid, is_ttp, is_bounce_point` | (L,T) bool | le_grids.f90:2715-2771, 2907-2946 | no | EXPORT | |
| `can_be_ttp` | (L,) bool | le_grids.f90:2907-2946 | no | EXPORT | |
| `is_lower_bounce_point, is_upper_bounce_point` | (L,T) bool | elemental (le_grids.f90:2775-2867) | no | BOTH (exported for checks) | 03 §3.1 |
| `jend` | (T,) int | le_grids.f90:2871-2903 | no | EXPORT | 0 on untrapped grids; g_adjust uses `jend`, nl_center uses both (05 P8) |
| `xi, wxi` | (T, L2) | le_grids.f90:2955-3020 | no | EXPORT | Fortran `(2*nlambda, theta)` |
| `ixi_to_il, ixi_to_isgn` | (T, L2) int | le_grids.f90:2972-3011 | no | EXPORT | 1-based values |
| `ng2, nlambda, lmax, nxi, negrid` | () | le_grids.f90:402-404, 1896-1904 | no | EXPORT | |
| `passing_wfb, trapped_wfb, mixed_wfb, new_trap_int, bouncefuzz` | () | le_grids.f90:413-427 | no | EXPORT | |
| `il_is_wfb/passing/trapped` | (L,) bool | le_grids.f90:450-491 | no | RECOMPUTE (export for checks) | |
| `is_passing_hybrid_electron` | (S,Y,L) bool | le_grids.f90:4490-4499 | no | RECOMPUTE (export for checks) | |
| `spec_{z,mass,dens,temp,tprim,fprim,uprim,uprim2,vnewk,nu_h,bess_fac,stm,zstm,tz,zt,smz,type,f0type}` | (S,) | species.f90:448-480 | no | EXPORT | derived values exported, not recomputed |
| `nonmaxw_corr, dlnf0drho` | (S,E) | calculate_f0_arrays (species.f90:563-597) | no | EXPORT | `wstar` rebuild needs `dlnf0drho` |
| `has_electron_species, has_ion_species, has_hybrid_electron_species` | () bool | species.f90:914-944 | no | EXPORT | adiabatic logic (07 §3.9.1) |

### A.5 Run control (`/run_parameters`, `/gs2_time`, `/gs2_reinit`, `/init_g`, `/fields`, `/nonlinear_terms`)

| name | NumPy shape | source | dt | decision | notes |
|---|---|---|---|---|---|
| `fphi, fapar, fbpar, beta, tite, zeff` | () | run_parameters.f90:300-366 (`fbpar` default −1.0; `beta == 0` forces `fapar = fbpar = 0`) | no | EXPORT (effective) | |
| `has_phi, has_apar, has_bpar` | () bool | run_parameters.f90:351-353 | no | EXPORT | |
| `ieqzip` | K2 bool | run_parameters.f90:259-272 | no | EXPORT (gate check) | must be all false (G-07) |
| `nstep, max_sim_time, immediate_reset, ncheck_stop, delt_option_switch, wstar_units` | () | read_parameters | no | EXPORT | |
| `tunits, wunits, woutunits` | (Y,) | adjust_time_norm (gs2_time.f90:83-126) | no | EXPORT | |
| `delt_adj, delt_cushion, delt_minimum, dt0, abort_rapid_time_step_change, in_memory` | () | init_reinit (gs2_reinit.f90:257-302) | no | EXPORT | |
| `dt_not_set = -1.0` | () | gs2_time.f90:19 | no | RECOMPUTE (constant) | |
| `tstart, new_field_init, scale, ginitopt_switch` | () | init_g | no | EXPORT | scale quirk (01 §6.12) |
| `fieldopt_switch, remove_zonal_flows_switch, force_maxwell_reinit` | () | fields read_parameters | no | EXPORT | |
| `nonlin, split_nonlinear, accelerated, cfl, use_cfl_limit, use_2d_cfl, use_order_based_error, error_target, istep_error_start, include_phi/apar/bpar, zip, nl_forbid_force_zero` | () | nonlinear_terms read_parameters | no | EXPORT | |
| `cflx, cfly` | () | `maxval(akx)*0.5/cfl`, `maxval(aky)*0.5/cfl` (nonlinear_terms.fpp:357-358) | no | BOTH | not dt-dependent |
| FFT scale `(1/nx)*(1/ny)` or `1/(nx*ny)` | () | gs2_transforms.fpp:968, fft_work.fpp:220 | no | RECOMPUTE | depends on `accelerated` (05 §4.6) |
| AB coefficients | (≤3,) | gs2_time.f90:250-312, every `invert_rhs` call with `istep > 0` | **dt** | RECOMPUTE | host float64, literal formula (01 §3.10) |
| `sourcefac` | () c | dist_fn.fpp:6166-6176 | `code_time` | RECOMPUTE | only used by `phiext_full` (gated when `phi_ext /= 0`) |

### A.6 Layout / summation-order tables (`/gs2_layouts`)

| name | NumPy shape | source | decision | notes |
|---|---|---|---|---|
| `layout` + `g_lo%*_ord` | attr / () | gs2_layouts.fpp:917-948 | EXPORT | fixes velocity-sum order (04 §3.1, 06 §5.5, 08 §3.1.6) |
| `iglo_of` | G5 int (0-based) | `idx(g_lo, ik,it,il,ie,is)` | EXPORT (derived) | builds the sequential reduction order |
| `ile_of`, `le_lo%*_ord` | (S,Y,X,T) int | `idx_le` (gs2_layouts.fpp:2160-2170) | EXPORT (derived) | collision rebuild integrals |

### A.7 Static per-point coefficients (`/dist_fn_arrays`, `/dist_fn`)

| name | NumPy shape | source | dt | decision | notes |
|---|---|---|---|---|---|
| `vpa, vpac` | G7 f64 | init_vpa_vpac (dist_fn.fpp:1592-1678) | no | EXPORT | `vparknob` already applied; `vpac(ntgrid)=0`; 0 at forbidden points |
| `vperp2` | G6 f64 | init_vperp2 (dist_fn.fpp:1560-1589) | no | EXPORT | |
| `aj0, aj1` | G6 f64 | init_bessel (dist_fn.fpp:1756-1823) via libm `bessel_j0/j1` | no | EXPORT | never recompute (G-02) |
| `ufac` | (S,E) | init_source_term (dist_fn.fpp:2180-2187), allocated once | no | EXPORT | `energy**1.5` → `pow` |
| `fexp, bkdiff` | (S,) | dist_fn.fpp:1319-1344, 1397-1398 | no | EXPORT (effective) | |
| scalar knobs `driftknob, tpdriftknob(eff), mach, g_exb, g_exbfac, g_exb_start_time, g_exb_start_timestep, omprimfac, phi_ext, t0, omega0, gamma0, source0, source_option_switch, wfb, gridfac, apfac, poisfac, afilter, esv(eff), nonad_zero, zero_forbid, boundary_option_switch, adiabatic_option_switch, def_parity, even, hyper_in_initialisation, start_from_previous_solution, boundary_off_grid, exponential_boundary, opt_source, mult_imp, gf_lo_integrate` | () | read_parameters (dist_fn.fpp:1174-1400) | no | EXPORT (effective values; must be read before any ladder rung, see G-20) | |
| `gamtot, gamtot1, gamtot2` | F3 | init_fieldeq (dist_fn.fpp:6470-6572) | no | EXPORT | final values, incl. `+tite`, poisfac and esv |
| `gamtot3` | F3 | dist_fn.fpp:6564-6565 (fieldlineavg only) | no | EXPORT if allocated | 07 RECOMPUTE is only a check |
| `gridfac1` | F3 | dist_fn.fpp:6486-6506 | no | EXPORT | |
| `awgt` | K2 | first `calculate_flux_surface_average` (dist_fn.fpp:8006-8020) | frozen | EXPORT if allocated | computed with the filter state of its first call |
| `inv_phi_denominator_g, inv_bpar_denominator_g` | F3 | dist_fn.fpp:7457-7585, lazily at the first get_fields_direct_from_dfn | no | BOTH if allocated | persist through dt resets (deallocated only at level 28, dist_fn.fpp:1131-1132) |

### A.8 dt-dependent `dist_fn` coefficients (rebuilt at level 31 on every dt change)

| name | NumPy shape | source | dt | decision | formula owner |
|---|---|---|---|---|---|
| `vpar` | G7 f64 | init_vpar (dist_fn.fpp:1681-1729) | **dt** | BOTH | 02 §3.1 (checked: `((((((zstm*tunits)*code_dt)*0.5)/delthet)*(|gp|+|gp+1|))*vpac)`) |
| `wdrift` | G7 f64 | init_wdrift (dist_fn.fpp:1405-1498), wdrift_func 1500-1530, wcoriolis_func 1534-1557 | **dt** | BOTH | 02 §3.1 (centring reads the uncentred `ig+1`, checked) |
| `wdriftttp` | G7 f64 | dist_fn.fpp:1451-1480 | **dt** | BOTH | 02 §3.1 (not centred) |
| `wstar` | (S,E,Y) | init_wstar (dist_fn.fpp:1746) `- code_dt*wunits(ik)*dlnf0drho(ie,is)` | **dt** | BOTH | 02 §3.1 |
| `a, b` | G7 c128 | init_invert_rhs (dist_fn.fpp:1892-1897, ttp 1919) | **dt** | BOTH | 02 §3.1 / 03 §3.3 (component forms agree) |
| `r, ainv` | G7 c128 | dist_fn.fpp:1885-1891, 1899-1921 | **dt** | BOTH | 03 §3.3 (Smith division) |
| `g_h` | G7 c128 | init_homogeneous_g (dist_fn.fpp:1937-2125) | **dt** | BOTH | 03 §3.4 |
| `g_adj` | (S,E,L,Y,X,2,n_links_max) | p-slots every call, h-slots from `g_h` (dist_fn.fpp:2120-2123) | dt (h-slots) | RECOMPUTE (gathers) | 03 §3.11; never export except for debugging |
| `source_coeffs_{phim,phip,aparm,aparp}` | (S,E,L,Y,X,2,N) c | init_source_term (2194-2264) | **dt** | GATED (`opt_source`) | 02 §3.6 |

### A.9 Linked-boundary and supercell topology (static)

| name | NumPy shape | source | decision | notes |
|---|---|---|---|---|
| `itleft, itright` | K2 int (1-based, −1) | compute_itleft_and_itright (dist_fn.fpp:2441-2513) | EXPORT | ky=0 box cells are self-linked |
| `l_links, r_links` | K2 int | count_links (2608-2672) | EXPORT | |
| `n_links` | (Y,X,2) int | count_links | RECOMPUTE (export for checks) | `2*links-1` or 0 |
| `itleftmost, itrightmost` | K2 int | get_leftmost/rightmost_it (3889-3920) | EXPORT (derived) | |
| `supercell_labels, n_cells` / `n_supercells` | K2 / (Y,) int | calculate_supercell_labels (2355-2397) | EXPORT | |
| `connections_neighbor` | G5 bool | compute_connections (2518-2565) | RECOMPUTE (export for checks) | |
| `save_h` | (S,E,L,Y,X,2) bool | set_save_h (2570-2604) | RECOMPUTE (export for checks) | |
| `no_connections` | () bool | dist_fn.fpp:2870-2876 | EXPORT | |
| `wfb_aj0l, wfb_aj1l, wfb_aj0r, wfb_aj1r` | **(S,Y,X,E)** | verbatim recompute of dist_fn.fpp:6016-6028 at `itl/itr` (09 §7.7) | EXPORT | already evaluated at the chain ends; JAX does **not** gather again (G-04) |
| `wfb_vperp2left, wfb_vperp2right` | (S,E) | dist_fn.fpp:6019-6020 | EXPORT (or exact recompute) | |
| parity partner `par_it` | K2 int | get_parity_conn (dist_fn.fpp:4742-4793) | RECOMPUTE | 07 §3.6 (verified: walks `itleft/itright` for all `il`) |
| `jump, ikx_indexed, dkx, dtheta0, itmin` | (Y,), (X,), () | exb_shear first call (4235-4267) | RECOMPUTE | host side |
| `i_class, N_class, M_class` | (), (C,), (C,) | setup_fields_implicit_class_arrays (2676-2755) | EXPORT | |
| `cls_ik, cls_it` | (C, Nmax, Mmax) int, −1 pad | finish_fields_layouts (fields_implicit.f90:1074-1133) | EXPORT | values 1-based |
| `sc_class, sc_member, sc_cell` | K2 int | derived | EXPORT | |

### A.10 Field-solve matrices

| name | NumPy shape | source | dt | decision | notes |
|---|---|---|---|---|---|
| `aminv_c<c>` | (M_c, R_c, R_c) c128 = `[m, B_out, B_in]` | init_inverse_matrix (fields_implicit.f90:722-1071) | **dt** (also collisions, hyper, vnmult at build time) | BOTH | raw runtime rows (09 §7.4); zonal-member rows are 0 |
| `nfield, nidx` | () | fields_implicit.f90:477-481 | no | EXPORT | |
| `am_c<c>` | (R_c, R_c·M_c) c | probe P-R1 (09 §5.2) | **dt** | optional EXPORT | validates the JAX probe builder |
| `resp_ik<ik>_is<is>` (`fieldmat`) | per supercell `(ax2, ax1)` c | `fields::dump_response_to_file` round trip (09 §7.5) | **dt** | BOTH iff `field_option='local'` | extended-domain order (04 §7) |
| `antot*, fieldeq*, fl_avg, g_work` | F3 / K2 / G7 | per call | — | RECOMPUTE | scratch |

### A.11 Collisions (`/collisions`; only if `colls`)

| name | NumPy shape | source (collisions.fpp) | dt | decision | notes |
|---|---|---|---|---|---|
| `vnew, vnew_s, vnew_D, vnew_E, delvnew` | (S,E,Y) | init_vnew 1707-1876 | no | EXPORT | erf/exp/pow |
| `vnewh` | (S,Y,X,T) | 1853-1863 (iff hyper_colls) | no | EXPORT | |
| `pitch_weights` | **(T,L)** (Fortran `(nlambda, theta)`, transposed relative to `wl`) | 2199-2210 | no | EXPORT | |
| `vpdiff` | (L,2,T) | init_vpdiff 4088-4140 | no | EXPORT (conservative) | |
| `vpdiffle` | (T,XI) | 4142-4152 | no | EXPORT | |
| `c1le, betaale, qle` | LE6(E, XI) | init_lorentz 2240-2285 | **dt**, vnmult(1) | BOTH | 06 §5.1 |
| `ec1le, ebetaale, eqle` | LE6(E, XI) | init_ediffuse 1981-1999 | **dt**, vnmult(2) | BOTH | 06 §5.2 |
| `z0le, s0le, w0le` | LE6(E1, XI) | init_lorentz_conserve 691-1194 | **dt**, vnmult(1) | BOTH | 06 §5.3 |
| `bz0le, bs0le, bw0le` | LE6(E1, XI) | init_diffuse_conserve 1198-1679 | **dt**, vnmult(2) | BOTH | 06 §5.4 |
| `aj0le, vperp_aj1le, vpa_aj0_le` | LE6(E1, XI) | init_le_bessel 634-687 | no | EXPORT | |
| `e15` | (E,) | `energy_maxw**1.5` (exporter, 09 §7.7) | no | EXPORT | |
| `capg_half` | (E−1,) | verbatim `capg` (2165-2176) | no | EXPORT | |
| scalars `collision_model_switch, colls, has_lorentz, has_diffuse, drag, use_le_layout, conserve_moments, conservative, resistivity, const_v, ei_coll_only, cfac, hypermult, hyper_colls, special_wfb_lorentz, vpar_zero_mean, conserve_forbid_zero, force_collisions, adjust, heating, split_collisions, lorentz_switch, ediff_switch, nxi_lim, vary_vnew` | () | read_parameters 466-562 | no | EXPORT | |
| standard-layout `c1/betaa/ql, ec1/ebetaa/eql, z0..bw0` | LZ6 / E7 / G7 | 2289-2353, 2016-2033, 1172-1192 | **dt** | GATED (`use_le_layout=F`, G-24) | |

### A.12 Hyperviscosity (`/hyper`)

| name | NumPy shape | source | dt | decision | notes |
|---|---|---|---|---|---|
| `D_res` | K2 | init_hyper (hyper.f90:400-411) | **frozen** at the first-init dt | EXPORT | never rebuilt on dt change (level 19 not torn down; hyper.f90:370) |
| `hypervisc_filter` (simple3D) | F3 | hyper.f90:419-445 | **frozen** | EXPORT (simple3D only) | |
| `hypervisc_filter` (iso/aniso) | F3 | hyper.f90:738-743, 777-778, every call | dt (current) | RECOMPUTE, transient | **not** diagnostic state (G-06) |
| `akx4_max, aky_max, aky4_max, akperp4_max, kperp2_max, kperp2_max_zonal` | () | hyper.f90:386-431 | no | BOTH | |
| effective `D_hypervisc, D_hyperres, nexp, hyper_on, hyper_option_switch, const_amp, isotropic_shear, damp_zonal_only, omega_osc, gridnorm, D_hyper3D, P_hyper3D, ky_cut, kx_cut, isotropic_model` | () | read_parameters (hyper.f90:452-601) | no | EXPORT | |

### A.13 Diagnostics (`/diagnostics` and the Fortran-written static part of `<run>.out.nc`)

| name | shape | source | decision | notes |
|---|---|---|---|---|
| parsed `gnostics` flags and scalars (`nwrite, navg, omegatol, omegatinst, exit_when_converged, igomega, write_*`, `calculate_fluxes`) | () | diagnostics_config read_parameters + check_parameters (08 §1.3) | EXPORT (post-`check_parameters`) | forced flags (`write_line` forces `write_fields` and `write_omega`) |
| static `.out.nc` content (metadata, dims, `nc_norms`, `nc_species`, `nc_geo`, `input_file`) and record 1 (istep 0) | | written by `gs2_export` with `run_diagnostics=.true.` | EXPORT (file) | JAX appends from record index 1 (`nout = 2`) (G-05) |

---------------------------------------------------------------------------------------------------

## B. Merged option gate table

`SUPPORT` = port; `GATE` = the JAX loader refuses a non-default value (soft gate recorded by the exporter,
09 §8.3, unless marked hard). "Owner" is the spec holding the exact algorithm. Where specs disagreed, the
resolution and its evidence are in §D.

### B.1 Time loop, run control, initial state

| option | namelist | default | decision | reason / owner |
|---|---|---|---|---|
| `nstep` | knobs | 100 | SUPPORT | 01 |
| `delt` | knobs | 0.1 | SUPPORT | 01 |
| `delt_option` | knobs | `'default'` | SUPPORT | `'check_restart'` state captured by the export (01, 02) |
| `max_sim_time` | knobs | 1e6 | SUPPORT | 01 §1.15 |
| `immediate_reset` | knobs | `.true.` | SUPPORT | both paths (01 §3.1, §3.5) |
| `avail_cpu_time, margin_cpu_time, ncheck_stop, progress_frequency` | knobs | 1e10, 300, 5, 0 | SUPPORT | host-side only |
| `seed` | knobs | 0 | SUPPORT | inert (init state exported; antenna gated) |
| `do_eigsolve` | knobs | `.false.` | GATE (hard) | not a time stepper |
| `use_old_diagnostics` | knobs | `.false.` | GATE | G-09 |
| `wstar_units` | knobs | `.false.` | GATE | G-16 |
| `eqzip_option` | knobs | `'none'` | GATE | G-07 |
| `fphi, fapar, fbpar, beta, tite, zeff` | knobs | 1, 0, −1, 0, 1, 1 | SUPPORT | export effective values (02, 04) |
| `delt_adj, delt_cushion, delt_minimum, dt0, abort_rapid_time_step_change` | reinit_knobs | 2, 1.5, 1e-5, delt, T | SUPPORT | 01 §3.6; a non-power-of-2 `delt_adj` disables ladder hits (G-19) |
| `in_memory` | reinit_knobs | `.false.` | SUPPORT | identical numerics (01 §1.12), apart from the scale quirk |
| `ginit_option` | init_g_knobs | `'default'` | SUPPORT | post-init state exported |
| `tstart, restart_file, restart_dir, read_many, include_explicit_source_in_restart` | init_g_knobs | | SUPPORT | captured by the export |
| `scale` | init_g_knobs | 1.0 | SUPPORT if 1.0, or a restart-type ginit, or `in_memory`; GATE otherwise | 01 §6.12 (init_g.f90:1168-1209 sets `scale = 1` only after restart-type options; gs2_save.fpp:865-876) |
| `new_field_init` | init_g_knobs | `.true.` | GATE if `.false.` | 01, 04 |
| `initial_condition_is_nonadiabatic_dfn` | init_g_knobs | `.false.` | SUPPORT | init only |
| `nensembles`, list mode, Trinity, parameter scans, optimisation mode | | | GATE (hard) | host-level job control (01, 09) |
| MPI ranks > 1, OpenMP threads > 1 | run config | | GATE (hard for reference data) | reduction order (04 P3, 09 §11.1) |

### B.2 Fields

| option | namelist | default | decision | reason / owner |
|---|---|---|---|---|
| `field_option` = `default/implicit` | fields_knobs | `'default'` | SUPPORT | 04 |
| `field_option` = `local/implicit_local` | fields_knobs | | SUPPORT (separate exact path; strict mode refuses a dt-rebuild miss until 04 §7.3 is ported) | G-23 |
| `field_option` = `gf_local`, `test` | fields_knobs | | GATE | 04 |
| `remove_zonal_flows_switch` | fields_knobs | `.false.` | SUPPORT | 01 §3.8 / 04 §3.7 |
| `force_maxwell_reinit` | fields_knobs | `.true.` | SUPPORT | 01 §3.9 |
| `dump_response, read_response, response_dir, response_file` | fields_knobs | | SUPPORT | exported `aminv` is what Fortran used (04 P23) |
| `field_subgath, do_smart_update, minnrow, field_local_allreduce(_sub), field_local_nonblocking_collectives, field_local_tuneminnrow` | fields_knobs | | SUPPORT | no numeric effect on 1 rank |
| `adiabatic_option` | dist_fn_knobs | `'default'` | SUPPORT (all: default, fieldlineavg, yavg) | 04 §3.4, 07 §3.9 |
| `esv` | dist_fn_knobs | `.false.` | SUPPORT | 04 §3.3 |
| `apfac, poisfac, gridfac` | dist_fn_knobs | 1, 0, 1 | SUPPORT | 04 |
| `afilter` | dist_fn_knobs | 0 | GATE | debug knob (04, 07) |
| antenna (`driver` namelist present and `ant_off = .false.`) | driver | off | GATE (hard in export) | RNG stream (01 §1.16, 07 §3.13) |

### B.3 Source, drifts, parallel solve, boundaries

| option | namelist | default | decision | reason / owner |
|---|---|---|---|---|
| `source_option` = `default/full`, `homogeneous` | source_knobs | `'default'` | SUPPORT | 02 §3.4 |
| `source_option` = `phiext_full` | source_knobs | | SUPPORT iff `phi_ext == 0`; GATE otherwise | G-14 |
| `phi_ext, t0, omega0, gamma0, source0` | source_knobs | 0, 100, 0, 0, 1 | SUPPORT (inert unless `phiext_full`) | 02 §3.2 |
| `fexpr, bakdif` | dist_fn_species_knobs | 0.4, 0 | SUPPORT | 02, 03 |
| `mult_imp` | dist_fn_knobs | `.false.` | SUPPORT | export effective `bkdiff` |
| `driftknob, tpdriftknob, vparknob, mach, omprimfac, btor_slab` | dist_fn_knobs | 1, sentinel, 1, 0, 1, 0 | SUPPORT | 02 §3.1 |
| `opt_source` | dist_fn_knobs (+ optimisation override) | `.false.` | GATE | different results (02 §3.6) |
| `boundary_option` = `default, zero, unconnected, self-periodic, periodic, kperiod=1, linked` | dist_fn_knobs | `'default'` | SUPPORT | 03 §3.7, 04 §3.11 |
| `boundary_option` = `alternate-zero` | dist_fn_knobs | | SUPPORT (exact alias of `zero`) | G-13 |
| `boundary_option='linked'` with `grid_option /= 'box'` | | | GATE | 03 §4 |
| `nonad_zero` | dist_fn_knobs | `.true.` | SUPPORT (both values) | 03 §3.8, §3.10 |
| `zero_forbid` | dist_fn_knobs | `.false.` | SUPPORT | 03 §3.7 |
| `wfb` | dist_fn_knobs | 1.0 | SUPPORT | 03 §3.4 |
| `wfbbc_option` | le_grids_knobs | `'default'` (mixed) | SUPPORT (`mixed, passing, trapped`) | 03 |
| `boundary_off_grid, exponential_boundary(_factor), start_from_previous_solution` | dist_fn_knobs | `.false.` | GATE | experimental (03 §3.12) |
| `def_parity, even` | dist_fn_knobs | `.false.`, `.true.` | SUPPORT | G-12 (03 §3.12.2, 07 §3.6) |
| `gf_lo_integrate` | dist_fn_knobs | `.false.` | GATE | different integration grouping (04) |
| `lf_default, lf_decompose, neo_test` | | | SUPPORT (inert) | LOWFLOW not built |
| LOWFLOW, SHMEM build macros | build | off | GATE (hard) | 07 §3.14, 09 §8.3 |
| species `type = 'hybrid_electron'` | species_parameters_N | | SUPPORT without collisions; GATE with collisions | G-25 |
| species `type = 'trace'` | | | SUPPORT | zero integral weights (07 §3.15) |
| species `f0type /= maxwellian` | species_parameters_N | | SUPPORT without collisions; GATE with collisions | collision operator uses the Maxwellian grid; not verified (07 Q7) |
| `trapped_particles, bouncefuzz, ngauss, npassing, radau_gauss_grid, new_trap_int, negrid, vcut, genquad` | le_grids_knobs | | SUPPORT | exported grids |
| `grid_option` = `single, range, specified, box/nonlinear` | kt_grids_knobs | `'default'` (single) | SUPPORT | box sets `kwork_filter(1,1)` |
| geometry (s-alpha, Miller, numerical/file) and drift knobs | theta_grid_* | | SUPPORT | exported effective geometry |

### B.4 Nonlinear term and dt control

| option | namelist | default | decision | reason / owner |
|---|---|---|---|---|
| `nonlinear_mode` | nonlinear_terms_knobs | `'default'` (off) | SUPPORT; GATE `on` with a non-box grid | 05 |
| `cfl, use_cfl_limit, use_2d_cfl` | nonlinear_terms_knobs | 0.95, T, T | SUPPORT | 05 §4.9 |
| `use_order_based_error, error_target, istep_error_start` | nonlinear_terms_knobs | F, 0.1, 30 | SUPPORT | G-11 (05 §4.10, verified) |
| `include_phi, include_apar, include_bpar` | nonlinear_terms_knobs | T | SUPPORT | 05 §4.3 |
| `nl_forbid_force_zero` | nonlinear_terms_knobs | T | SUPPORT | 05 §4.8 |
| `zip` | nonlinear_terms_knobs | F | GATE | G-08 |
| `split_nonlinear` + all `split_nonlinear_terms_knobs` (`strang_split, advance_nonadiabatic_dfn, split_method, rk_method, tolerances`) | | F | GATE | different algorithm; clobbers `code_dt_prev1/2` (05 §4.13) |
| `fft_measure_plan` | layouts_knobs | F | GATE | non-reproducible plans (05 P3) |
| `fft_use_wisdom, fft_wisdom_file` | layouts_knobs | T | SUPPORT (recommend F or no stale wisdom) | 09 §11.14 |
| `layout` | layouts_knobs | `'lxyes'` | SUPPORT (all six) | summation order exported |
| `local_field_solve, intmom_sub, intspec_sub, opt_local_copy, opt_redist_*, unbalanced_*, nproc_*_lo` | layouts_knobs | | SUPPORT | no numeric effect on 1 rank |
| `gf_local_fields` | layouts_knobs | F | GATE (with gf_local) | |

### B.5 Collisions

| option | namelist | default | decision | reason / owner |
|---|---|---|---|---|
| `collision_model` (`default/full, lorentz, ediffuse, lorentz-test, none/collisionless`) | collisions_knobs | `'default'` | SUPPORT | 06 §4-5 |
| `use_le_layout` | collisions_knobs | `.true.` | SUPPORT T; GATE F | G-24 |
| `conserve_moments` | collisions_knobs | `.true.` | SUPPORT (except the no-trapped-grid gate) | 06 §4.5, §4.8 |
| no trapped particles ∧ Lorentz active ∧ (`conserve_moments` ∨ ¬`use_le_layout`) | | | GATE | out-of-bounds bug, confirmed (06 §8.1) |
| `conservative` | collisions_knobs | `.true.` | SUPPORT T; GATE F | uninitialised endpoints (06 §0) |
| `lorentz_scheme, ediff_scheme` | collisions_knobs | `'default'` | SUPPORT default; GATE `'old'` | 06 |
| `resistivity, const_v, ei_coll_only, cfac, hypermult, special_wfb_lorentz, vpar_zero_mean, conserve_forbid_zero, force_collisions, adjust` | collisions_knobs | | SUPPORT | 06 |
| `nu_h, vnewk` | species | 0 | SUPPORT | 06 |
| `split_collisions` | collisions_knobs | `.false.` | SUPPORT | G-10 |
| `timesteps_between_collisions` | collisions_knobs | 1 | SUPPORT (ignored by GS2) | 06 §4.10 |
| `heating` | collisions_knobs | `.false.` | SUPPORT (forced off unless gated diagnostics; inert otherwise) | G-26 |
| `vary_vnew, ncheck, vnfac, vnslow, etol, etola, ewindow, ewindowa` | collisions_knobs | F | SUPPORT (no-op while `write_verr` is gated) | G-15 |
| `test` | collisions_knobs | F | SUPPORT (ignored) | |

### B.6 Hyperviscosity, flow shear, miscellaneous

| option | namelist | default | decision | reason / owner |
|---|---|---|---|---|
| `hyper_option` (`default/none, visc_only, res_only, both, simple3D`) | hyper_knobs | `'default'` | SUPPORT | 07 §3.2 |
| `d_hypervisc, d_hyperres, d_hyper, const_amp, isotropic_shear, omega_osc, nexp, gridnorm, d_hyper3d, p_hyper3d, ky_cut, kx_cut, isotropic_model` | hyper_knobs | | SUPPORT | 07; GATE the zero-normaliser cases (07 App. B) |
| `damp_zonal_only` | hyper_knobs | F | SUPPORT (replicate the bug: no mode filtered; `hypervisc_filter` not updated) | 07 P12 |
| `include_kpar` | hyper_knobs | F | SUPPORT (unused) | |
| `hyper_in_initialisation` | dist_fn_knobs | T | SUPPORT | 07 §1.5, 04 §5.4 |
| `g_exb, g_exbfac, g_exb_start_time, g_exb_start_timestep` | dist_fn_knobs | 0, 1, −1, −1 | SUPPORT; GATE `g_exb /= 0` with `ntheta0 == 1` | 01 §3.7, 07 §3.3 |
| `g_exb_error_limit` | dist_fn_knobs | 0 | SUPPORT (unused) | |

### B.7 Diagnostics (`gs2_diagnostics_knobs`)

| option | default | decision | reason / owner |
|---|---|---|---|
| `write_any, nwrite, navg, omegatol, omegatinst, exit_when_converged, igomega` | T, 10, 10, 1e-3, 1e6, T, 0 | SUPPORT | 08 §3.3 |
| `write_omega, write_fields, write_phi/apar/bpar_over_time, write_ql_metric` | T (forced), T (forced), F, T | SUPPORT | 08 §3.3-3.5 |
| `write_fluxes, write_fluxes_by_mode, print_flux_line, write_flux_line` | F (T if NL), F (T if NL), F, F | SUPPORT (triggers g_adjust round trips) | 08 §3.6 |
| `write_moments, write_ntot/density/upar/tperp_over_time` | T, F | SUPPORT (round trip) | 08 §3.7 |
| `write_correlation, write_eigenfunc` | T, F | SUPPORT | 08 §3.8-3.9 |
| `write_final_fields, write_final_moments, write_final_epar, write_final_antot` | F | SUPPORT | 08 §3.11 |
| `print_line, write_line, write_ascii` | F, T, T | SUPPORT (text optional, never refuse) | 08 §3.10 |
| `nc_sync_freq, nwrite_mult, enable_parallel, serial_netcdf4`, legacy no-op flags | | SUPPORT (no-op) | 08 |
| `save_for_restart, nsave, save_many, file_safety_check, save_glo_info_and_grids, save_velocities` | | SUPPORT (host I/O; inert) | 09 §6.4 |
| `use_nonlin_convergence` (+ `conv_*`) | F | GATE | G-22 (Trinity-only experimental logic) |
| `save_distfn, append_old, make_movie, dump_fields_periodically, write_heating, write_collisional, write_verr, write_max_verr, write_cerr, write_symmetry, write_nl_flux_dist, write_parity, write_cross_phase, write_correlation_extend, write_jext, write_lorentzian, write_zonal_transfer, write_kinetic_energy_transfer, write_full_moments_notgc, write_g, write_gyx, write_gs, write_kpar, write_final_db` | F | GATE | 08 §5 (extra round trips, feedback into `vnmult`, or non-standard) |

---------------------------------------------------------------------------------------------------

## C. Full ordered call sequence of one time step (with the spec that owns each call)

Single rank, `field_option='implicit'` (local differences in brackets). Numbers in `[..]` are owner-spec
sections. Every call listed was traced in the Fortran; nothing in the step is unowned.

```
C.0  gs2_main.fpp:717-718  istep loop; exit if istep > nstep (repeated evolve calls only)             [01 §1.5]
C.1  gs2_main.fpp:730  fields::advance(istep) (fields.f90:586-605): select fieldopt_switch            [01 §1.6]
C.1.1  fields_implicit.f90:333  antenna_amplitudes(apar_ext) iff .not.no_driver                      GATED [01 §1.16, 07 §3.13]
C.1.2  :335  exb_shear(gnew, phinew, aparnew, bparnew, istep) iff kx_shift/theta0_shift allocated    [01 §3.7, 07 §3.3]
       [local: fields_local.f90:4280, test abs(g_exb*g_exbfac) > eps -- equivalent]
         exb_first setup (dist_fn.fpp:4235-4267) -> guards istep_last / start_timestep / start_time (4270-4272)
         -> gdt = 0.5*(code_dt+code_dt_old) (4283) -> istep_last = istep (4286)
         -> kx_shift/theta0_shift and jump update (4300-4310) -> return if all(jump==0) (4314)
         -> theta0 branch (4341-4489, ik ascending, all-ky field-zeroing bug) or kx branch (4490-4599, ik = naky..2)
C.1.3  :337-338  g = gnew ; phi = phinew ; apar = aparnew ; bpar = bparnew   [local: only has_* fields] [01 §3.2]
C.1.4  :341  timeadv(phi..., phinew..., istep)  -- pass 1, phinew == phi                              [01 §3.3]
   a  dist_fn.fpp:4046-4058  split_nonlinear branch                                                   GATED [05 §4.13]
   b  :4060  add_explicit_terms (nonlinear_terms.fpp:407-439)                                        [05 §4.1, 02 §3.7]
        linear: save_dt_cfl(1e8)
        nonlinear and istep /= istep_last (:495):
          istep_last = istep; gexp_3 = gexp_2; gexp_2 = gexp_1 (:497, 505-506)
          max_vel_components = eps (:518)
          add_nl (:520, 788-929): loads i*kx*chi (+bpar, -apar) -> transform2 -> max|.|*cfly
                                   -> g_adjust(from) of g -> *i*ky -> transform2 -> bracket = ba*gb
                                   -> i*ky*chi -> transform2 -> max|.|*cflx -> i*kx*h_w -> transform2
                                   -> bracket -= ba*gb -> inverse2 -> *kxfac                          [05 §4.2-4.7]
          nl_center (:523, 677-760)                                                                   [05 §4.8 = 02 §3.7]
          CFL / order-based error -> set_cfl_time_step_and_check_if_too_large (:529-583)             [05 §4.9-4.10, 01 §1.14]
          if reset: return (:585-587)
        zip block (:599-609)                                                                         GATED
        istep_last = istep (:611)
   c  :4066  if (reset .and. immediate_reset) return                                                 [01 §3.3]
   d  :4068  invert_rhs (6152-6193)                                                                  [03 §3.6]
        sourcefac (6166-6176)                                                                         [02 §3.2]
        g_work = 0 (6179)
        get_source_term_on_local_domain (5363-5420): AB coefficients (5388-5389); kwork_filter skip (5401);
          get_source_term per isgn (4856-5130): phigavg/apargavg, set_source, phiext, b*g + a*g, ttp, explicit
                                                                                                      [02 §3.3-3.5]
        invert_rhs_1 (5440-5707): kwork_filter / ieqzip skip; hybrid set_h_zero; flags; BCs;
          sigma=2 sweep; lower-bounce source; trapped wfb; sigma=1 sweep; self_periodic; trapped beta;
          trapped-wfb beta; zero_forbid; enforce_parity (5706)                                       [03 §3.7-3.10, 07 §3.6]
        apply_linked_boundary_conditions (5919-6149) iff linked                                      [03 §3.11]
   e  :4073-4080  hyper_diff(gnew, phinew)  (istep > 0, or istep == 0 with hyper_in_initialisation)  [07 §3.2]
   f  :4083-4087  vspace_derivatives iff .not.split_collisions (le_derivatives.fpp:15-174):
        g_adjust(from, phinew, bparnew) -> gather g2le -> solfp1_le_layout (ediffuse -> conserve_diffuse
        -> drag -> lorentz -> conserve_lorentz) -> scatter -> g_adjust(to)                           [06 §4.1-4.8]
   g  :4091  enforce_parity iff def_parity                                                           [07 §3.6]
   h  :4094-4110  Strang second half                                                                  GATED [05 §4.13]
   i  :4114  timeadv::istep_last = istep                                                             [01 §2.5]
C.1.5  fields_implicit.f90:343  if (reset) return                                                     [01 §3.2]
C.1.6  :345  aparnew += apar_ext                                                                      GATED
C.1.7  :352  getfield (222-312): get_field_vector -> getfieldeq -> getan_from_dfn(gnew)
         (integrate_species, esv) -> getfieldeq1 (fl_avg) -> u = -aminv.fl (per-cell sums)
         -> get_field_solution                                                                       [04 §3.1-3.6]
       [local: getfield_local (fields_local.f90:4193-4232), matvec order of G-23]                    [04 §7]
C.1.8  :354-356  phinew += phi ; aparnew += apar ; bparnew += bpar  [local: has_* only]              [04 §3.7]
C.1.9  :362  remove_zonal_flows iff remove_zonal_flows_switch (374-425; ik = 1 only)                 [01 §3.8, 04 §3.7]
C.1.10 :365  timeadv(..., istep, diagnostics=1) -- pass 2: as C.1.4, but add_explicit skips the
             NL block (istep == istep_last); collisions run with heating_flag possible               [01 §3.3]
C.1.11 :369  collisions_advance iff split_collisions (dist_fn.fpp:4130-4167):
             vspace_derivatives -> g_adjust(from) -> calculate_potentials_from_nonadiabatic_dfn
             -> g_adjust(to) with the new fields                                                     [06 §4.10, 04 §3.9]
C.2  gs2_main.fpp:739-752  immediate-reset path (only if reset is set by C.1.4.b):
       prepare_initial_values_overrides; set_initval_overrides_to_current_vals (gs2_main.fpp:1166-1218)
       reset_time_step (gs2_reinit.f90:104-225): rebuild of levels 30-32 (§E)                        [01 §3.6]
       restore state (gs2_init.fpp:292-375), fields from get_init_field                              [01 §3.9, 04 §3.8]
       reset = .false.; advance(istep) again: exb skipped (istep_last); NL skipped; gexp reused      [01 §1.11]
C.3  :754  exit if state%exit
C.4  :759-769  gs2_save_for_restart iff nsave > 0 .and. mod(istep,nsave) == 0 (pure write)         host, inert
C.5  :771  update_time (gs2_time.f90:161-174)                                                        [01 §3.4]
C.6  :775-783  run_diagnostics(istep) (gs2_diagnostics_new.f90:339-519):                             [08 §1.4]
       calculate_omega every step (istep > 0) -> may set exit                                        [08 §3.3]
       write block iff mod(istep,nwrite)==0 .or. exit: t; write_omega; write_fields; write_ql_metric;
         calculate_fluxes (g_adjust round trips on gnew AND g); write_moments (round trip on gnew);
         correlation; eigenfunc; text lines                                                          [08 §3.2-3.10]
C.7  :785  if exit: converged
C.8  :800  check_time_step (gs2_reinit.f90:228-254) -> :811 reset_time_step (no re-advance)         [01 §3.5-3.6]
C.9  :818-822  checkstop (every ncheck_stop), checktime                                             host
C.10 :829-836  istep_end = istep; exit if user_time + user_dt > max_sim_time                         [01 §1.15]
C.11 :860  exit if exit
After the loop: gs2_main.fpp:870-881 forced run_diagnostics(istep - offset, force=.true.);
finalize_diagnostics final routines (08 §1.5).
```

Response-matrix construction (init and every dt rebuild) uses the same C.1.4 chain with `istep = 0`
(fields_implicit.f90:669; fields_local.f90:2184): no explicit/NL terms (nonlinear_terms.fpp:424-426,
dist_fn.fpp:5382), hyper iff `hyper_in_initialisation`, collisions iff `.not.split_collisions`, parity, and
the linked BCs. Owners: 04 §3.12-3.13 and §5.4-5.5. Inside a rebuild, the drag term reads the **module**
`aparnew`, which holds the probe there (G-27).

---------------------------------------------------------------------------------------------------

## D. Gaps and contradictions, with resolutions from the Fortran

Each entry gives the specs involved, the problem, the evidence, and the binding resolution.

**G-01 `kwork_filter` while stepping (06 vs 01/02/03/04/07/09).** Spec 06 §3.3 says the filter is
"always `.false.` while stepping; only set during fields_local response build. Port: assert all false". The
Fortran: `if(gridopt_switch.eq.gridopt_box) kwork_filter(1,1)=.true.` at the end of `init_fields`
(fields.f90:376), and level 13 (`kt_grids`) is never torn down, so the flag persists through every dt reset.
Collisions skip filtered modes in every stage (collisions.fpp:3014, 3067, 3130, 3148, 3280, 3309, 3333,
3361, 3419, 3442, 3562, 3588, …). **Resolution:** JAX reads the exported `kwork_filter` and applies
06's `active(it,ik) = .not.kwork_filter .and. .not.ieqzip` mask in all collision stages; the assertion is
removed. With the mask, the `(1,1)` mode stays bit-exactly frozen. Every operation that touches it is
exact: the source is zero there; `invert_rhs_1` passes it through (dist_fn.fpp:5475); hyper filter = 1;
exb never shifts `ik = 1`; and the g_adjust round trips add `adj = 0` because `phi(1,1) = 0` (its integral
is skipped, so `get_init_field` returns 0).

**G-02 Bessel implementation (02 §2.2, 05 §3.2 vs 03 §0.2, 07 §0.2, 09 §0.1).** Build macro
`-DSPFUNC=_SPF200X_` (`bin/gs2 --build-config`; Makefile:455; `USE_LOCAL_SPFUNC` empty, Makefile:64);
spfunc.fpp:737-790: `dj0 = bessel_j0(x)`, `dj1 = bessel_j1(x)/x`, `0.5` iff `is_zero(x)` (exact compare,
warning_helpers.f90:87-90). **Resolution:** `aj0`, `aj1` and the mixed-wfb factors are always EXPORTed;
the exporter's verbatim recompute (09 §7.7) calls `spfunc::j0/j1`. Specs 02 and 05 are corrected on
provenance only.

**G-03 Complex storage in the export file (02 §8 vs 09 §0.2; 03 Q8).** Spec 02 proposes `_r/_i` variable
pairs; spec 09 uses a leading Fortran `ri` dimension (trailing in NumPy), as `gs2_io` and
`gs2_save_response` do (gs2_io.fpp:55, 95; gs2_save.fpp:938-948). **Resolution:** the 09 convention for the
export and trace files. Only GS2 restart files keep `gr/gi`, `phi_r/phi_i`, `gexp_*_r/_i` (gs2_save.fpp:274-428).

**G-04 Mixed-wfb end Bessel factors (03 §3.11.4 vs 09 §7.7).** Spec 03 defines `wfb_aj0l…` per cell
`it`, evaluated with the cell's own `kperp2`, shape NumPy `(S,E,Y,X)`, "JAX gathers at itl/itr". Spec 09
evaluates them already at `itl = get_leftmost_it(it,ik)` / `itr`, Fortran `(egrid,kx,ky,species)`.
Fortran (dist_fn.fpp:6016-6028) uses `kperp2(-ntgrid,itl,ik)` and `kperp2(ntgrid,itr,ik)`, i.e. a function
of `(ie,is,ik,itl)` only. **Resolution:** the 09 definition is binding. The arrays are pre-gathered and have
NumPy shape `(S,Y,X,E)`; JAX uses `wfb_aj0l[is,ik,it,ie]` directly for cell `(it,ik)` with **no second
gather**. Following 03's text literally would gather twice. Note the association
`argl = (bess_fac*smz)*sqrt(((energy*al)/bmag(-nt))*kperp2(-nt,itl,ik))`; it differs from `init_bessel`
(03 §3.11.4).

**G-05 Export point of the t = 0 state (08 P1 vs 01 §7, 09 §11.2).** Spec 08 recommends dumping before
`initialize_diagnostics` and letting JAX write record 1; specs 01 and 09 dump after it. Fortran:
`run_diagnostics(0)` (gs2_main.fpp:600) mutates `gnew` and `g` through `g_adjust` round trips (08 §3.1.8).
**Resolution:** `/state_t0` (after `initialize_diagnostics`) is the JAX initial state. `gs2_export` with
`run_diagnostics=.true.` writes the static part and record 1 of `<run>.out.nc`; JAX opens it in append mode
and starts at `nout = 2`, with `run_diagnostics::istep_last = 0`. `/state_t0_prediag` exists only to unit-test
a JAX implementation of istep-0 diagnostics. JAX must replicate the round trips on its own write steps
(08 P1).

**G-06 `hypervisc_filter` is not diagnostic state (08 §2.2, P11, Q4; 09 §3.11).** Fortran `eexchange`
(diagnostics_fluxes.fpp:1063-1211) uses `hypervisc_filter` (1127-1128) and `wdrift_func/code_dt` only to
build the integrand of `exchange_dummy` (1156). Its caller says "we don't use / write exchange_dummy
anywhere" (273). The written `es_energy_exchange` is `exchange_sym` (1199), built from
`aj0*0.25*(gnew ± g)` (1169, 1181), `phinew` and `phi`. **Resolution:** JAX carries no hyper filter for
diagnostics and may skip the first half of `eexchange`. The rest of 08 §3.6.5 stands; the `(1,1)` skip
for nonlinear runs remains (1092-1101).

**G-07 `eqzip_option` (SUPPORT in 01 §4, 06 §6 vs GATE in 03 §4, 07 §3.5, 08 §5).** Fortran honours
`ieqzip` inconsistently. `invert_rhs_1` skips such modes (dist_fn.fpp:5477), and collisions skip them. The
source still computes them (5401 tests only `kwork_filter`), the linked BCs still modify them (5996-6010, no
`ieqzip` test), `hyper_diff` filters them twice per step, the field solve solves them, and omega is
re-enabled in nonlinear runs (08 §3.3). **Resolution: GATE** (`any(ieqzip)` refused). No spec gives the
complete interplay.

**G-08 `zip` (SUPPORT in 01 vs GATE in 02, 05).** Fortran nonlinear_terms.fpp:599-609 zeroes `g` and
`gexp_1` for `ik /= 1` on every `add_explicit` call that does not return early, i.e. in both `timeadv`
calls, before `invert_rhs` reads `g`. **Resolution: GATE** (expert option; the algorithm is in 05 §4.12 if it
is ever un-gated).

**G-09 `use_old_diagnostics` (SUPPORT in 01 vs GATE in 08).** The old module produces different
outputs and its own round trips (`loop_diagnostics`, gs2_main.fpp:604-607, 777). No spec covers it.
**Resolution: GATE.**

**G-10 `split_collisions` (SUPPORT in 01, 06 vs GATE in 04 §4).** Fortran dist_fn.fpp:4130-4167 (read in
full): `vspace_derivatives`; `g_adjust(gnew, phinew, bparnew, from_g_gs2)`;
`calculate_potentials_from_nonadiabatic_dfn` (7322-7455, checked against 04 §3.9: `phi_weight =
1/(Σ_s ((dens*z)*z)/temp [+ tite])` summed sequentially, `apar = antota*inv_kperp2`,
`bpar = antotp*(-beta/bmag**2)`, no `fl_avg`); `g_adjust(..., to_g_gs2)` with the **updated** fields.
`timeadv` then skips collisions (4083), including in the response build. **Resolution: SUPPORT**, with the
algorithm in 06 §4.10 plus 04 §3.9. Not re-applied after it: `remove_zonal_flows`.

**G-11 `use_order_based_error` (SUPPORT in 01, 05 vs GATE in 02 §4).** Fortran
nonlinear_terms.fpp:540-556 and 620-672 match 05 §4.10 exactly, including the
`error_target*cfl*cfl` early return and `target_dt = code_dt*sqrt(error_target/error)`. The max-abs
reductions are exact. **Resolution: SUPPORT** (05 §4.10).

**G-12 `def_parity` (GATE in 03 vs SUPPORT in 04, 07).** Fortran enforce_parity (dist_fn.fpp:4796-4853),
parity redistribute (4686-4737) and get_parity_conn (4742-4793) were read. `scatter` moves data from the
`to` list `(ig, 2, iglo_conn)` to the `from` list `(-ig, 1, iglo)`, so `gnew(-ig,1,iglo) = gnew(ig,2,conn)`
for every `ig`, including 0, in linked mode. The connection walks `itleft/itright` directly
(get_left_connection 3837-3860) for every `il`, including trapped ones. The non-linked branch leaves
`gnew(0,1)` untouched. Both 03 §3.12.2 and 07 §3.6 describe this correctly. **Resolution: SUPPORT**, applied
at both call sites (5706 inside `invert_rhs_1` before the linked BCs, and 4091 after collisions), also during
response builds.

**G-13 `boundary_option = 'alternate-zero'` (SUPPORT alias in 03 vs GATE in 04).** The switch value
`boundary_option_alternate_zero` (dist_fn.fpp:124) appears only in a `write` (770). Every behavioural
`select`/`if` treats it as the default case: init_fieldeq `case default` (6500-6506), invert_rhs_1
(5504-5506), init_homogeneous_g (1965-1967), linked tests (3529, 3641, 4621, 4819). **Resolution: SUPPORT
as an exact alias of `zero`.**

**G-14 `source_option = 'phiext_full'` (SUPPORT in 02 vs GATE in 04).** Fortran dist_fn.fpp:4910-4920
adds `-zi*wdrift*nonmaxw*2*phi_ext*sourcefac*aj0` for `aky < eps` with no `istep` guard. The term is
therefore present in every `timeadv(…,0)` probe, and the response columns become affine. **Resolution:
GATE when `phi_ext /= 0`; with `phi_ext == 0` the option is numerically identical to `full` (SUPPORT).**

**G-15 `vary_vnew` (no-op in 01 vs GATE in 04, 06, 08).** `adjust_vnmult` is reachable only through
`get_verr` ← `write_velocity_space_checks`. In the new diagnostics that runs only if `write_verr`, on write
steps (gs2_diagnostics_new.f90:459) or every `ncheck` steps when `write_verr .and. vary_vnew` (509-514).
**Resolution:** `write_verr` is GATED, so `vary_vnew` has no effect: **SUPPORT (inert)**. With old
diagnostics (gated) it would be live.

**G-16 `wstar_units` (SUPPORT linear-only in 07; units flow through 01/02/03/06 vs GATE in 08).** Every
stepper formula consumes the exported `tunits/wunits`, but diagnostics need `t_wstar` and `woutunits` scaling
(gs2_diagnostics_new.f90:436-440; 08 §3.3), `tunits = 0` for `aky = 0` ("garbage", gs2_time.f90:115-120),
and 05 requires it to be false for nonlinear runs. **Resolution: GATE in v1.** Un-gating needs only the two
diagnostic changes plus a linear-only and `aky /= 0` check.

**G-17 `kwork_filter` during the dt-reset response rebuild (01 §6.11, 03 §6.6, 04 P6/P11, 07 P11).**
The specs claim that the `(1,1)` block of the rebuilt implicit matrix differs from the initial build. The
Fortran shows otherwise. (i) In init_inverse_matrix every pivot sets `a_inv(:,jlo) = 0` for members with
`aky == 0 .and. akx == 0` (fields_implicit.f90:834, 856), so that member's `aminv` is identically zero
whatever `am` holds. (ii) No other member couples to it. Its `iso/aniso_shear` contribution is
`|phi|²·(akx²+aky²)²·fac = 1·0·1 = 0` exactly (hyper.f90:766-768, 702-713). Source, invert, collisions and
fl_avg are per mode; the `(1,1)` cell is an `N = 1` chain skipped by the linked fix (dist_fn.fpp:6000-6002).
**Resolution:** the asymmetry is real but numerically irrelevant in box runs. Non-box runs never set the
filter. JAX may rebuild with either filter state; no special handling is required.

**G-18 Hyperviscosity grouping in rebuilds (04 P9, 07 P11; consistent).** When the shear rate depends
on amplitude (`const_amp = F`, iso or aniso), the probe pattern matters: implicit perturbs all members of
one class at one `(n, ig, f)`. This is stated consistently; no change. Linear hyper (`const_amp = T`) and
simple3D allow free batching.

**G-19 dt-ladder coverage (01 §5.5 vs 09 §6.3/§7.8).** Spec 09 generates only `dt_init/delt_adj^k`. The
reachable set also contains `min(dt*delt_adj, code_dt_max)` moves (gs2_reinit.f90:86-101, 176-183). These
go **above** `dt_init` when `dt0 > delt`; a linear run with `dt0 > delt` does so at step 1
(01 §6.22). For non-power-of-2 `delt_adj` the values are path-dependent (05 P9). **Resolution:** the ladder
is a cache only, and the RECOMPUTE path (02/03/06/04 formulas) is mandatory. Optionally extend the ladder
with upward rungs `min(dt_init*adj^j, dt_max)` and with `dt_max/adj^k`. Key all lookups by the exact float64
bit pattern. Ladder rungs must contain every array marked **dt** in §A.8, §A.10 and §A.11, and must reproduce
the frozen `D_res`/simple3D filter (09 already states the latter).

**G-20 Config objects reset by teardown (gap in 09).** `finish_dist_fn_level_3` calls
`dist_fn_config%reset()`, `source_config%reset()` and deallocates `dist_fn_species_config`
(dist_fn.fpp:1151-1153). `finish_collisions` calls `collisions_config%reset()` (collisions.fpp:4211). After
any dt reset or ladder rung, `get_dist_fn_config()`, `get_source_config()` and the collisions config may
return defaults or re-read values, while the module variables keep the effective ones (collisions re-read the
namelist on re-init, collisions.fpp:466-562). **Resolution:** the exporter reads every knob from the module
variables made public in 09 §4.2, or from `get_*_config()` **before** the first ladder rung. It must never
read a config object after a teardown.

**G-21 Spec 02 export timing (02 §8 "after init(full) and before the first advance").** This conflicts
with G-05. **Resolution:** use `/state_t0` (post-diagnostics). Coefficients are unaffected by the
diagnostics.

**G-22 Nonlinear convergence (SUPPORT in 08 §5, but no algorithm).** Fortran
diagnostics_nonlinear_convergence.f90:49-117 uses Trinity-specific counters (`trin_istep += nwrite` on every
write call even when the write was forced by exit), integer divisions `conv_nstep_av/nwrite`, a ring buffer
`conv_heat` (allocated, initialised only through `reset_averages_and_counters`), and `exit = .true. .and.
exit_when_converged`. **Resolution: GATE** (`use_nonlin_convergence = .true.` refused). The above is the
algorithm summary if it is ever ported.

**G-23 `field_option = 'local'` exactness (01, 04 §7, 09 §7.5).** SUPPORT is consistent across specs, but
the exact algorithm is only partly written. Fortran: per cell and per present field,
`rb%tmp_sum = matmul(vect, rb%data)` (fields_local.f90:370-384); cell `tmp_sum = 0;
tmp_sum(rows_f) = tmp_sum(rows_f) - rb(f)%tmp_sum` for `f` in phi, apar, bpar order (512-546); supercell
`local = 0; local = local + cell%tmp_sum` over cells in order (735-750); then the gather and duplicate fill
(04 §7.2). So `u = ((0 + T_1) + T_2) + …` with `T_c = ((0 - S_phi,c) - S_apar,c) - S_bpar,c`. **Resolution:**
SUPPORT for the initial dt with the exported `fieldmat` (09 §7.5) and this matvec order. On a dt change, use a
ladder hit, or the literal `invert_block` port (04 §7.3, pending Q-04). Strict mode refuses a ladder miss
until that port is validated. The order of the matmul inner sum is still open (Q-04).

**G-24 `use_le_layout = .false.` (06 "run le algorithm" vs 09 soft gate).** The standard-layout path
sums in `iglo` order through `integrate_moment_c34` with different associations (06 §4.9), and on untrapped
grids it carries the out-of-bounds bug (06 §8.1). **Resolution: GATE** until the standard-layout path is
ported. Running the le algorithm instead would be a ≈1e-15 per-step deviation, not a replication.

**G-25 Hybrid electrons (SUPPORT in 03, 07 vs GATE in 06).** The dist_fn parts are fully specified
(`set_h_zero` 03 §3.9 / 07 §3.8; linked-fix skip). The collision hooks (`set_hzero_lorentz_collisions_matrix`
with asymmetric ranges, `zero_out_passing_hybrid_electrons`) are written in 06 §5.1/§5.3 but flagged unverified
(06 Q8). **Resolution: SUPPORT without collisions; GATE hybrid together with collisions.**

**G-26 `heating` (06 inert vs 08 forced-off rule).** With `heating = .true.`, the pass-2 collision call
round-trips `gold = g` (le_derivatives.fpp:72, 108). `g` is not read again in that step and is overwritten at
the next step (fields_implicit.f90:337). Diagnostics turn heating off unless `write_heating` or
`write_collisional` is set (gs2_diagnostics_new.f90:115-121), and both are gated. **Resolution: SUPPORT
(inert); `c_rate` is not ported.**

**G-27 Drag term input during rebuilds (gap, 04/06).** The drag term reads the module `aparnew`
(06 §3.3, collisions.fpp:3041, 3074). `init_response_row` calls
`timeadv(phi,apar,bpar,phinew,aparnew,bparnew,0)` on the module arrays (fields_implicit.f90:651, 669), so
during a rebuild the drag sees the **probe** `aparnew`. **Resolution:** the JAX collision kernel takes
`aparnew` as an explicit argument equal to the `aparnew` passed to `timeadv`. In the time step that value is
`apar` in pass 1 and the solved `aparnew` in pass 2; in a rebuild it is the probe.

**G-28 Response-build state that is not zeroed (01 §1.4, 04 P20; consistent, made explicit).**
`init_response_matrix` zeroes `g`, the fields and `am`, but not `gnew` (fields_implicit.f90:553-557).
`reset_time_step` zeroes `gnew` first (gs2_reinit.f90:163); the 09 ladder does not. **Resolution:** this
is irrelevant to the matrices. `invert_rhs_1` starts each solved mode from 0 (dist_fn.fpp:5563-5575, with
sfps gated). The only pass-through mode is the filtered zonal `(1,1)` (G-17), and integrals skip it.
After any rebuild, JAX must restore `gnew` and the fields from saved state (01 §3.9); a rebuild must never
touch the evolving buffers (04 P20).

**G-29 Index-order traps across specs (consistency check; no contradiction).** `pitch_weights` is
Fortran `(nlambda, theta)` → NumPy `(T, L)`, whereas `wl` and `forbid` are NumPy `(L, T)`. `xi, wxi,
ixi_to_*` are NumPy `(T, L2)`. `vpdiff` is NumPy `(L, 2, T)`. `vnewh` is NumPy `(S, Y, X, T)`. `wstar` is
NumPy `(S, E, Y)`. LE6 arrays are NumPy `(S, Y, X, T, NE, NX)`. `aminv_c` is NumPy `[m, row(B, output),
col(B, input)]` (04 §0.1 and 09 §7.4 agree). `cls_it/cls_ik` are NumPy `(C, Nmax, Mmax)` with 1-based values.
The loader must not "fix" these.

**G-30 Sign, order and dt-dependence cross-checks (no contradiction found).** The following were checked
against the Fortran and agree in every spec that states them:
`a/b/r/ainv` component forms (02 §3.1 vs 03 §3.3; dist_fn.fpp:1885-1897); `vpar`, `wstar`, `wdrift_func`,
`wcoriolis_func`, and the wdrift centring reading the uncentred `ig+1` (1405-1557, 1715-1746); AB
coefficients (gs2_time.f90:282-310); the explicit term
`source + (0.5*code_dt)*((ab1*g1 + ab2*g2) + ab3*g3)` (dist_fn.fpp:5106-5127, ttp 4987-5010); g_adjust
association `((bf*vperp2)*aj1)*bpar + ((pf*phi)*aj0)` with `pf = (zt*nm)*facphi` (dist_fn_arrays.fpp:155-168);
getfieldeq1 associations (dist_fn.fpp:6620-6659); field update `phinew = u + phi`; the NL bracket sign and
`max_vel_x/y` pairing with `cflx/cfly` (nonlinear_terms.fpp:820-914); exb `gdt` and `NINT`
(dist_fn.fpp:4283-4310); frozen `D_res`/simple3D filter (hyper.f90:370); `ufac` allocated once
(dist_fn.fpp:2180); `inv_*_denominator_g` persisting through resets (dist_fn.fpp:1131-1132); fields zeroed at
the level-32 teardown (fields.f90:618-620).

**G-31 Owner of `nl_center` (02 Q6).** It is specified identically in 02 §3.7 and 05 §4.8. **Resolution:**
05 owns it; 02 is a consistent copy.

**G-32 Private-data access (02 Q9, 03 Q3, 06 Q2, 07 Q10, 04 Q10).** Resolved by 09 §4 (`public ::` blocks,
code generation proven unchanged) and probe P-R1 for `am`. The JAX side treats SAVE locals as RECOMPUTE
(09 §3.13).

**G-33 Two collision applications and the drag in pass 1 (06 §8.9).** In pass 1 the module
`aparnew == apar` (copied at fields_implicit.f90:338), so the drag uses the old `apar`. This is consistent
with passing `aparnew` explicitly (G-27).

**G-34 `kx_shift`/`jump` host-side (01 §3.7, 07 §3.3.2).** The specs are consistent: these depend only on dt
history, `aky`, `tunits` and inputs. The bookkeeping is updated in the failed attempt of an immediate reset
with the pre-reset `gdt` and is not repeated (dist_fn.fpp:4270, 4286). JAX state must carry the exb
`istep_last` across the reset.

---------------------------------------------------------------------------------------------------

## E. Total recompute set on a dt change (all subsystems)

Trigger: `reset_time_step` (gs2_reinit.f90:104-225) on either path. The rebuild goes down to level 29
(`override_timestep`, :166) and back up to `full` (:216). Levels ≤ 28 survive (gs2_init_level_list.inc).

| step | what | file:line | JAX action |
|---|---|---|---|
| E.1 | scalar update: `code_dt` (reduce `/delt_adj` or increase `min(*delt_adj, dt_max)`, ≤ 10 attempts), `user_dt = code_dt` (save_dt); exits on `nconsec > 4`, `code_dt/delt_adj <= code_dt_min` (tested even when increasing), 10 attempts, `code_dt <= code_dt_min` | gs2_reinit.f90:133-206 | host, literal (01 §3.6) |
| E.2 | **unchanged**: `code_dt_prev1/2`, `code_dt_old`, `code_time`, `user_time`, `code_dt_cfl`, `kx_shift/theta0_shift`, `jump`, exb `istep_last`, `exb_first`, `nl istep_last`, `gexp_*` (restored with identical values), `kwork_filter`, `vnmult` (via `vnm_init`, collisions.fpp:4174-4175, 1731-1733), `omegahist` | | keep |
| E.3 | level 30 collisions: `vnew*`, `vnewh`, `pitch_weights`, `vpdiff(le)`, `aj0le…` rebuilt with **identical** values; **dt-dependent**: `c1le, betaale, qle, ec1le, ebetaale, eqle, z0le, s0le, w0le, bz0le, bs0le, bw0le` | collisions.fpp:436-631, 4165-4212 | RECOMPUTE 06 §5.1-5.4 (host NumPy, layout summation order §5.5), or ladder |
| E.4 | level 31 dist_fn: **dt-dependent** `vpar, wdrift, wdriftttp, wstar, a, b, r, ainv, g_h` (+ `g_adj` h-slots); `source_coeffs_*` (gated); **not** rebuilt: `ufac`, `D_res`, simple3D `hypervisc_filter` (frozen at the first-init dt), `aj0/aj1/vperp2/vpa/vpac`, `gamtot*`, `gridfac1`, `awgt`, `inv_*_denominator_g` | dist_fn.fpp:960-1013, 1140-1154; hyper.f90:370 | RECOMPUTE 02 §3.1, 03 §3.3-3.4, or ladder |
| E.5 | level 32 fields: fields zeroed (fields.f90:618-620); response matrix rebuilt by `timeadv(…,0)` probes using E.3/E.4 coefficients, current-dt iso/aniso hyper (if `hyper_in_initialisation`), frozen `D_res`/simple3D, collisions (unless split), the current `kwork_filter` (irrelevant, G-17) and the probe `aparnew` in the drag (G-27); `kwork_filter(1,1) = .true.` re-asserted; the iso/aniso `hypervisc_filter` is left holding the last probe's filter (harmless, G-06) | fields.f90:353-377, 608-637; fields_implicit.f90:459-1071 | RECOMPUTE 04 §5.4-5.5 (implicit) / 04 §7.3 (local), or ladder |
| E.6 | level 34 state restore: `g = gnew = saved gnew` (or re-read from the restart file just written, with the `scale` quirk), `vnmult` set, `gexp_*` restored; fields: `force_maxwell_reinit` → `get_init_field(gnew)` (dist_fn.fpp:7590-7799) else the saved/file fields; `phi = phinew` … | gs2_init.fpp:292-375 | 01 §3.9, 04 §3.8 |
| E.7 | counters: `reset_time_step::istep_last = istep`; `timeadv::istep_last` becomes 0 through the probes (used only by gated split NL) | gs2_reinit.f90:220; dist_fn.fpp:4114 | |
| E.8 | immediate path only: `advance(istep)` again; exb is skipped and the NL is not recomputed or reshuffled; AB coefficients recomputed with the new `code_dt` and old `prev1/prev2` | gs2_main.fpp:751 | 01 §3.2 |
| E.9 | every step regardless: AB coefficients, exb `gdt`, iso/aniso hyper filter, omega `/code_dt`, drag `code_dt` factor, `wdrift_func/code_dt` (dead, G-06) | | already in the step |

---------------------------------------------------------------------------------------------------

## F. Remaining open questions

1. **Q-01 FMA / reference build policy** (01 Q4, 02 Q1, 03 Q1-Q2, 04 Q1, 05 Q1, 06 Q3, 07 Q1, 08 P10).
   Validate against the shipped `-O3` FMA-contracted binary at ~1e-12, or rebuild the reference
   (`bin/gs2` and `bin/gs2_export`) with `-ffp-contract=off` and run JAX with the strict
   `optimization_barrier` helpers (03 §3.0, verified bit-exact on XLA:CPU)? XLA:GPU behaviour is unverified.
2. **Q-02 dt-change contract** (01 Q2, 09 Q9): RECOMPUTE is mandatory (G-19). Is the ladder cache also
   wanted, and with upward rungs?
3. **Q-03 Bug replication policy**: theta0-branch all-ky zeroing (dist_fn.fpp:4370, 4387, 4404),
   `damp_zonal_only` (hyper.f90:776), frozen `D_res`/simple3D filter, stale `aminv` under `vary_vnew`, the
   no-trapped-grid collision OOB bug (gated). This document assumes "replicate everything not gated"; confirm.
4. **Q-04 `field_option='local'`**: summation order of gfortran's inlined `matmul` / libgfortran
   `matmul_c8` (04 Q3); is exact replication of `invert_block` required for dt rebuilds, or is the ladder
   plus strict refusal enough (04 Q2)?
5. **Q-05 Diagnostics as perturbation** (01 Q1, 09 Q3): reference runs keep the default diagnostics, so
   JAX replicates the `g_adjust` round trips on write steps (G-05). Confirm that no diagnostics-free
   reference is wanted.
6. **Q-06 Geometries with `bmax` not at `±ntgrid`** (03 Q6, §6.9): the trapped bounce matching silently
   skips `ig = ntgrid`. Does any supported input hit this?
7. **Q-07 `shat == 0` exactly with zonal modes** (02 Q2): `wdrift_func` divides by `shat`
   (dist_fn.fpp:1519, 1550). The proposal is to GATE `shat == 0` with any `aky == 0` mode unless
   `boundary_option` makes the term irrelevant; needs a decision.
8. **Q-08 Non-box `aky = akx = 0` modes** (04 Q7): GS2 freezes their fields (zero `aminv` rows) while
   evolving `g`; replicate (assumed yes).
9. **Q-09 Edge case `1 - al*bmag == -bouncefuzz`** (05 P8, 06 Q10): `forbid` and `jend` disagree. Should
   the exporter assert that no exported grid hits it?
10. **Q-10 Restart from a mid-run trace snapshot** (09 Q11): would need the SAVE counters of §A.1 moved to
    module scope (behaviour-preserving source edit, out of scope for this phase).
11. **Q-11 Exact build command for `bin/gs2`** (09 Q13): must be recovered so that `bin/gs2_export` is
    bitwise-compatible.
12. **Q-12 Non-Maxwellian `f0type` and hybrid electrons with collisions** (06 Q8, 07 Q6/Q7): currently
    GATED (B.3); un-gate only after an empirical check against `bin/gs2`.
