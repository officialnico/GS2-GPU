# 07 — Hyperviscosity and miscellaneous time-step pieces

GS2 8.2.1 (`gpu-port` branch, HEAD `4d8c94bcf`) → JAX float64/complex128 port.
Status: **contract for implementers** (read-only analysis; no source was modified).

Subsystems covered:

| § | Item | Decision |
|---|------|----------|
| 3.2 | `hyper.f90` — hyperviscosity / hyperresistivity (`hyper_diff`, `D_res`, all `hyper_option`s) | SUPPORT (all options) |
| 3.3 | ExB flow shear (`exb_shear`, `g_exb`, `g_exbfac`, `kx_shift`, `theta0_shift`) | SUPPORT (both kx and theta0 branches) |
| 3.4 | `kwork_filter` | SUPPORT |
| 3.5 | `ieqzip` / `eqzip_option` | GATE |
| 3.6 | `def_parity` / `even` / `enforce_parity` | SUPPORT |
| 3.7 | `wstar_units` (`tunits`, `wunits`) | SUPPORT (linear, all aky≠0) |
| 3.8 | Hybrid electrons (`is_passing_hybrid_electron`, `set_h_zero`) | SUPPORT (dist_fn part; collisions caveat) |
| 3.9 | Adiabatic response (`adiabatic_option`, `tite`, `fl_avg`) | SUPPORT |
| 3.10 | `poisfac` | SUPPORT (via exported `gamtot`) |
| 3.11 | `esv` | SUPPORT |
| 3.12 | `remove_zonal_flows_switch` | SUPPORT |
| 3.13 | Antenna (`antenna.f90`, `apar_ext`) | GATE |
| 3.14 | `lowflow.f90` (LOWFLOW build) | GATE |
| 3.15 | `species.f90` (types, adiabatic detection, f0 arrays) | SUPPORT |
| 3.16 | `afilter` (debug knob) | GATE |

---

## 0. Scope, conventions, build facts

### 0.1 Sources read (all of them, in full where relevant)

`src/hyper.f90` (all 816 lines), `src/antenna.f90` (all), `src/antenna_data.f90`, `src/species.f90` (all),
`src/lowflow.f90` (header + `init_lowflow`), and from `src/dist_fn.fpp`: module header/config (1–527),
`read_parameters` (1174–1375), init levels (849–1154), `init_wdrift`/`init_vpar`/`init_wstar` (1405–1750),
`init_invert_rhs` (1826–1929), `init_bc` (2128–2156), `init_source_term` (2159–2265), `compute_jshift0`,
`compute_itleft_and_itright` (2401–2512), connections (3837–3920), `allocate_arrays` (3923–4003),
`timeadv` (4014–4116), `exb_shear` (4183–4600), parity (4603–4853), source terms (4856–5420),
`invert_rhs_1` (5440–5707), `apply_linked_boundary_conditions` (5919–6018), `invert_rhs` (6152–6193),
`ensure_single_val_fields_pass[_r]` (6201–6261), `getan_from_dfn` (6273–6360), `init_fieldeq` (6470–6572),
`getfieldeq1[_nogath]` (6575–6832), `calculate_potentials_from_nonadiabatic_dfn` (7322–7448),
`get_fields_direct_from_dfn` (7457–7789), `calculate_flux_surface_average` (7994–8035).
Also: `fields.f90`, `fields_implicit.f90` (advance/remove_zonal/response row), `fields_local.f90`
(advance, `fm_populate`, `fm_getfieldeq1_nogath`), `fields_gf_local.f90` (advance), `kt_grids.f90`,
`kt_grids_box.f90`, `kt_grids_range.f90`, `run_parameters.f90`, `gs2_time.f90`, `gs2_main.fpp` (time loop),
`gs2_reinit.f90` (`reset_time_step`), `gs2_init.fpp` + `gs2_init_level_list.inc` + `gs2_init_subroutines.inc`,
`gs2_save.fpp` (kx_shift restart), `dist_fn_arrays.fpp` (`set_h_zero`), `le_grids.f90`
(`integrate_species_original`, `il_is_*`, `is_passing_hybrid_electron`), `collisions.fpp` (hybrid hooks),
`theta_grid.f90` (`field_line_average`), `externals/utils/integration.f90`, `externals/utils/summation.f90`.

### 0.2 Build facts that matter for round-off (from `bin/gs2 --build-config`)

* GCC/gfortran 15.2.0, `-O3 -mcpu=apple-m1` (aarch64), `-fdefault-real-8 -fdefault-double-8`.
  No `-ffast-math`, but **GCC's default `-ffp-contract=fast` is in force → FMA contraction of `a*b±c` is
  likely in the reference binary** (see §6 P1).
* MPI: T, **OpenMP: F**, Lowflow: **F**, F200X intrinsics (Bessel via intrinsic), FFTW3, NetCDF.
* `epsilon(0.0)` with `-fdefault-real-8` is the **double** epsilon `2.220446049250313e-16` everywhere below.
* Literals like `0.5`, `2.0`, `epsilon(0.)` are all real(8).

### 0.3 Index conventions used in this document

* Fortran: `ig ∈ [-ntgrid, ntgrid]`, `isgn ∈ {1,2}`, `it ∈ [1, ntheta0]`, `ik ∈ [1, naky]`,
  `il ∈ [1,nlambda]`, `ie ∈ [1,negrid]`, `is ∈ [1,nspec]`. `nth := 2*ntgrid+1`.
* Canonical distribution array in NumPy/JAX: `G[is, ie, il, ik, it, isgn, ig']` with
  `ig' = ig + ntgrid` (0-based), shape `(nspec, negrid, nlambda, naky, ntheta0, 2, nth)`.
* Fields: Fortran `phi(ig,it,ik)` → NumPy `PHI[ik, it, ig']`, shape `(naky, ntheta0, nth)`.
* `(ntheta0,naky)` Fortran arrays (`D_res`, `kwork_filter`, `ieqzip`, `itleft`, `awgt`, `fl_avg`, ...)
  → NumPy `(naky, ntheta0)`.
* `(naky)` arrays (`kx_shift`, `theta0_shift`, `tunits`, `wunits`, `aky`, `jump`) → NumPy `(naky,)`.
* Velocity/species arrays keep Fortran order in the netCDF metadata and are reversed in NumPy:
  Fortran `nonmaxw_corr(negrid,nspec)` → NumPy `(nspec, negrid)`; Fortran `wstar(naky,negrid,nspec)` → NumPy
  `(nspec, negrid, naky)`; Fortran theta arrays `x(-ntgrid:ntgrid)` → NumPy `(nth,)`.
* In pseudo-code, Fortran 1-based indices are kept unless written `p0`/`[...]` (0-based).
  Parenthesisation `((a*b)*c)` is written explicitly wherever Fortran's left-to-right rule fixes the order.

### 0.4 Mixed real/complex arithmetic (applies to every formula below)

gfortran lowers `complex*real`, `real*complex` and `complex/real` **component-wise** (GCC's complex-lowering
pass treats the promoted real as "ONLY_REAL"): `z*r = (z.re*r, z.im*r)`, `z/r = (z.re/r, z.im/r)`.
**In JAX always do these component-wise** (`lax.complex(z.real*r, z.imag*r)`); never let XLA promote `r` to
complex and run a complex multiply/divide (XLA complex divide is not guaranteed to reduce to `re/r, im/r`).
`real(conjg(z)*z)` must be computed as `z.re*z.re + z.im*z.im` (the Fortran complex product's real part is
`zr*zr - (-zi)*zi`, numerically identical up to FMA, see P1).

---

## 1. Control flow

### 1.1 Initialisation (levels, `gs2_init_level_list.inc`)

`gs2_init.fpp:237-268` (`init_pass_int`) walks levels one by one; going *down* to level L finishes every level
> L (`offset = 0`, loop at 263-266).

| Level | Routine | What this subsystem gets there |
|---|---|---|
| 13 `kt_grids` | `kt_grids.f90:158-191` | `akx, aky, theta0, ikx` (`get_grids`); **`kwork_filter = .false.`** (188-189); `kperp2`, `inv_kperp2` (`init_kperp2` 313-349) |
| 14 `run_parameters` | `run_parameters.f90:247-273` | `ieqzip` from `eqzip_option` (259-272); `wstar_units`, `tite`, `fphi/fapar/fbpar` (`read_parameters` 276-383) |
| 15 `species` | `species.f90:381-481` | `spec(:)` incl. derived `stm,zstm,tz,zt,smz` (476-480) |
| 16 `gs2_time` | `gs2_time.f90:45-69` → `adjust_time_norm` 83-126 | `code_dt`, `tunits(naky)`, `wunits(naky)`; `code_dt_old = 0.` (static init, line 31) |
| 18 `le_grids` | (velocity-grid spec) | `energy`, `al`, `ng2`, `forbid`; `calculate_f0_arrays` via `egrid.f90:66/106` → `nonmaxw_corr`, `dlnf0drho`, `f0_values` |
| **19 `hyper`** | `hyper.f90:360-449` `init_hyper` | `D_res(ntheta0,naky)` (uses **current `code_dt`**), `D_v`, `D_eta`, scalars `akx4_max...`; simple3D `hypervisc_filter` (uses current `code_dt`) |
| 20 `antenna` | `antenna.f90:196-224` | `no_driver` (259); antenna data (gated) |
| 21 `dist_fn_parameters` | `dist_fn.fpp:1174-1375` | `g_exb, g_exbfac, g_exb_start_*`, `def_parity, even, esv, poisfac, adiabatic_option_switch, hyper_in_initialisation, boundary_option_switch`; `esv = esv .and. linked` (1327) |
| 26 `dist_fn_arrays` | `dist_fn.fpp:885-925` → `allocate_arrays` 3923-4003 | `kx_shift(naky)=0` **or** `theta0_shift(naky)=0` iff `abs(g_exb*g_exbfac) > epsilon(0.)` (3983-3991) |
| 27 `dist_fn_level_1` | `init_bc` 2128-2156 | `l_links, r_links, itleft, itright` (`init_connected_bc`), `parity_redist` (2142) |
| 28 `dist_fn_level_2` | `init_fieldeq` 6470-6572 | `gamtot` (+`kperp2*poisfac`, esv, +`tite`), `gamtot1/2`, `gamtot3`, `gridfac1` |
| 29 `override_timestep` | (no teardown work) | — |
| 31 `dist_fn_level_3` | `dist_fn.fpp:960-1013` | `init_wstar` (988), `init_hyper` again at 1006 (**no-op**, `initialized` guard hyper.f90:370), `init_source_term` (1009) |
| 32 `fields_level_2` | `fields.f90:353-377` `init_fields` | response matrices (calls `timeadv(...,istep=0)`, see §1.6); **then `kwork_filter(1,1) = .true.` if box (376)** |
| 34 `set_initial_values` | init_g / restart | `g`, `gnew`, fields; restart restores `kx_shift`/`theta0_shift` (`gs2_save.fpp:820-830`) |

### 1.2 One time step (`gs2_main.fpp:717-829`)

```
do istep = istep_loop_start, istep_loop_max          ! gs2_main.fpp:717 (starts at 1)
   call advance(istep)                                ! :730 -> fields.f90:586-605 dispatch on field_option
   if (reset) then                                    ! :739-751 immediate_reset path (CFL)
      set_initval_overrides_to_current_vals; reset_time_step(...)   ! gs2_reinit.f90:104-225
      call advance(istep)                             ! retry SAME istep
   end if
   ...
   call update_time                                   ! :771 -> gs2_time.f90:161-174
                                                      !   code_dt_old = code_dt ; code_time += code_dt
   check_time_step / reset_time_step (non-immediate path)   ! :800-816
end do
```

`advance_implicit` (`fields_implicit.f90:315-371`; `advance_local` `fields_local.f90:4244-4342` and
`advance_gf_local` `fields_gf_local.f90:4239-4350` have the same order):

```
1  if (.not.no_driver) call antenna_amplitudes(apar_ext)            ! :333   [GATED]
2  if (allocated(kx_shift).or.allocated(theta0_shift)) &
       call exb_shear(gnew, phinew, aparnew, bparnew, istep)        ! :335   (§3.3)
       ! local/gf_local test abs(g_exb*g_exbfac) > epsilon(0.) instead (fields_local.f90:4280) -- equivalent
3  g = gnew ; phi = phinew ; apar = aparnew ; bpar = bparnew          ! :337-338
4  call timeadv(phi,apar,bpar,phinew,aparnew,bparnew,istep)          ! :341   first half (predictor)
5  if (reset) return                                                 ! :343
6  if (.not.no_driver) aparnew = aparnew + apar_ext                  ! :345   [GATED]
7  call getfield(phinew,aparnew,bparnew)                             ! :352   field solve (fl_avg, esv used inside)
8  phinew = phinew + phi ; aparnew = aparnew + apar ; bparnew = bparnew + bpar   ! :354-356
9  if (remove_zonal_flows_switch) call remove_zonal_flows            ! :362   (§3.12)
10 call timeadv(phi,apar,bpar,phinew,aparnew,bparnew,istep,diagnostics)  ! :365  second half
11 call collisions_advance(...)                                      ! :369   (only if split_collisions)
```

### 1.3 Inside `timeadv` (`dist_fn.fpp:4014-4116`)

```
NL terms (4046-4062); if (reset .and. immediate_reset) return (4066)
call invert_rhs(...)                                   ! 4068 -> 6152-6193
   sourcefac; g_work = 0; get_source_term_on_local_domain   ! 6181 -> 5363-5420
        per iglo: if (kwork_filter(it,ik)) cycle             ! 5401  (§3.4)
        set_source uses D_res(it,ik)                         ! 5087/5097, opt: 2227/2252 (§3.2.6)
   invert_rhs_1(phinew, bparnew, g_work)                     ! 6184 -> 5440-5707
        per iglo: if (kwork_filter(it,ik)) cycle             ! 5475  (§3.4)
                  if (ieqzip(it,ik)) cycle                   ! 5477  (§3.5)
                  if (is_passing_hybrid_electron(is,ik,il)) then set_h_zero; cycle   ! 5487-5490 (§3.8)
                  ... solve ...
        if (def_parity) call enforce_parity(parity_redist)   ! 5706  (§3.6)
   if linked: apply_linked_boundary_conditions               ! 6190 -> 5919
        skips kwork_filter (6000), passing hybrid (6005)
if ((istep == 0 .and. hyper_in_initialisation) .or. istep > 0) &
     call hyper_diff(gnew, phinew)                          ! 4073-4080 (§3.2)
if (.not. split_collisions) call vspace_derivatives(...)     ! 4083-4087 (collisions spec)
if (def_parity) call enforce_parity(parity_redist)          ! 4091  (§3.6)
strang-split NL second half (4094-4112)
istep_last = istep (4114)
```

So **hyper damping is applied in both `timeadv` calls of every step**, each time to the freshly solved `gnew`
with the `phinew` passed in (first call: `phinew == phi` = old fields; second call: new fields).

### 1.4 Reset / retry path — what is *not* re-initialised

`reset_time_step` (`gs2_reinit.f90:104-225`): `zero_array(gnew)` (163), `init(..., override_timestep)` (166,
tears down levels 35…30), change `code_dt` (`save_dt` 198), `init(..., full)` (216, rebuilds 30…35).

Consequences for this subsystem:

* **Level 19 `hyper` is never torn down** ⇒ `D_res` and the simple3D `hypervisc_filter`, which contain
  `code_dt`, keep the value of `code_dt` **at the first `init_hyper` call for the whole run** (the
  `init_hyper` call at `dist_fn.fpp:1006` is a no-op because of `hyper.f90:370`). iso/aniso filters are
  recomputed every call with the current `code_dt`, so they are not affected.
* `kx_shift`, `theta0_shift` (level 26), `jump`, `ikx_indexed`, `dkx`, `dtheta0`, `exb_first` (level 27
  resets `exb_first` only, `dist_fn.fpp:1080`), `awgt` (level 27), `kwork_filter` (level 13), `tunits/wunits`
  (level 16) all survive a timestep reset.
* `exb_shear`'s saved `istep_last` (`dist_fn.fpp:4217`) survives, so the retry `advance(istep)` does
  **not** shift again (4270). The `kx_shift` increment for that step used the **pre-reset** `gdt`.
* `code_dt_old` is only written in `update_time` (`gs2_time.f90:169`), so the step after a reset uses
  `gdt = 0.5*(new_dt + previous_dt)`.

### 1.5 Response-matrix construction hooks (istep = 0)

`init_response_row` (`fields_implicit.f90:650-719`) calls `timeadv(...,0)` (669) with unit-delta `phinew`
etc.; `fm_init_next_field_points` (fields_local) does the same. With the default
`hyper_in_initialisation = .true.` (`dist_fn.fpp:368`), **`hyper_diff` runs during response construction**,
and for iso/aniso with `const_amp = .false.` the shear rate is computed from the delta-function `phinew`
of that row. Response matrices are normally EXPORTed (fields spec); if the JAX side ever recomputes them
(e.g. after a dt change) it must replicate this exactly (§6 P11).

---

## 2. Inventory of everything this subsystem reads

Legend: **E** = EXPORT from Fortran, **R** = RECOMPUTE in JAX (formula given in §3), **B** = both (export, and
formula given for recompute/validation). "dt" = depends on `code_dt`. "priv" = module-private in Fortran, so
the export program must recompute it or a tiny accessor must be added (this phase is read-only).

### 2.1 Grids / geometry (owned by other specs, listed because read here)

| Fortran name | Module | Fortran shape | NumPy shape | Computed at | dt | Rec. |
|---|---|---|---|---|---|---|
| `akx` | kt_grids | `(ntheta0)` | `(ntheta0,)` | `box_get_grids` kt_grids_box.f90:382-389; range 330-347 | no | E |
| `aky` | kt_grids | `(naky)` | `(naky,)` | kt_grids_box.f90:378-380; range 299-311 | no | E |
| `theta0` | kt_grids | `(ntheta0,naky)` | `(naky,ntheta0)` | kt_grids_box.f90:391-401; range 313-326 | no | E |
| `ikx` | kt_grids | `(ntheta0)` int | `(ntheta0,)` | kt_grids_box.f90:382-389; range 350-352; single 151 | no | E |
| `kperp2` | kt_grids | `(-ntgrid:ntgrid,ntheta0,naky)` | `(naky,ntheta0,nth)` | `init_kperp2` kt_grids.f90:313-329 | no | E |
| `box`, `reality`, `gridopt_switch` | kt_grids | scalars | — | `get_sizes` kt_grids.f90:269-291 | no | E |
| `gds2, gds21, gds22` | theta_grid | `(-ntgrid:ntgrid)` | `(nth,)` | geometry | no | E (simple3D only) |
| `shat` | theta_grid | scalar | — | geometry | no | E |
| `theta`, `jacob` | theta_grid | `(-ntgrid:ntgrid)` | `(nth,)` | `finish_init` theta_grid.f90:401-521 (`jacob` 516) | no | E |
| `field_line_average_weight` | theta_grid (priv) | scalar | — | theta_grid.f90:520 | no | R: `1.0/trapz_klein(theta, jacob)` |
| `ntheta, ntgrid, nperiod` | theta_grid | int | — | `ntgrid = ntheta/2 + (nperiod-1)*ntheta` (theta_grid.f90:567) | no | E |
| `ng2` | le_grids | int | — | velocity-grid spec | no | E |
| `vperp2`, `aj0`, `aj1` | dist_fn_arrays | `(-ntgrid:ntgrid, iglo)` | `(nspec,negrid,nlambda,naky,ntheta0,nth)` | 1560-1589, 1756-1823 | no | E (other spec) |

### 2.2 Time normalisation (`gs2_time`)

| Fortran name | Shape | Computed at | dt | Rec. | Formula |
|---|---|---|---|---|---|
| `tunits` | `(naky)` | `adjust_time_norm` gs2_time.f90:108-124 | no | B | `wstar_units ? (aky/=0 ? 2.0/aky : 0.0) : 1.0` |
| `wunits` | `(naky)` | same | no | B | `wstar_units ? 1.0 : aky/2.0` |
| `code_dt` | scalar | `init_delt` 140-151, `save_dt` 206-213 | — | E (initial) | state |
| `code_dt_old` | scalar | static `0.`(31), `update_time` 169 | yes | E (initial = 0 for fresh start **and** restart) | state |
| `code_time` | scalar | `init_tstart`, `update_time` 170 | — | E | state |

### 2.3 Hyper (`hyper.f90`)

| Fortran name | Shape | Computed at | dt | Rec. | Note |
|---|---|---|---|---|---|
| `D_res` (public) | `(ntheta0,naky)` → NumPy `(naky,ntheta0)` | hyper.f90:400-411 (zeroed 615) | **yes, frozen at initial dt** | B | consumed by source term (§3.2.6) |
| `hypervisc_filter` (public), simple3D | `(-ntgrid:ntgrid,ntheta0,naky)` → `(naky,ntheta0,nth)` | hyper.f90:432-444 (init 1.0 at 613) | **yes, frozen at initial dt** | E (B) | never modified by `simple3Dfilter` |
| `hypervisc_filter`, iso/aniso | same | hyper.f90:738-743 / 777-778, every call | yes (current dt) | R | do **not** export: after response construction it holds a garbage row |
| `akx4_max, aky_max, aky4_max, akperp4_max` (priv) | scalars | hyper.f90:386-396 | no | R | §3.2.2 |
| `kperp2_max, kperp2_max_zonal, ky_cut, kx_cut` (priv, effective) | scalars | hyper.f90:421-431 | no | R | §3.2.2 |
| effective `D_hypervisc, D_hyperres, nexp, hyper_on, hyper_option_switch` (priv) | scalars | `read_parameters` hyper.f90:452-601 | no | R | resolve from `get_hyper_config()` with the switch logic of §3.2.1 |
| `D_v, D_eta` (public) | scalars | hyper.f90:378-379, 407, 414 | no | — | diagnostics only (`gs2_heating.fpp:972`); not read by the step |

### 2.4 ExB flow shear

| Fortran name | Module | Shape | Computed at | dt | Rec. |
|---|---|---|---|---|---|
| `g_exb, g_exbfac` (public) | dist_fn | scalars | dist_fn.fpp:1256, 1260 (overridable 1162) | no | E |
| `g_exb_start_time, g_exb_start_timestep` (priv) | dist_fn | scalars | 1258-1259 | no | E via `get_dist_fn_config()` |
| `kx_shift` | dist_fn_arrays | `(naky)` | allocated/zeroed 3985-3986; updated 4301-4303; restart 820-824 | yes (state) | E (initial) |
| `theta0_shift` | dist_fn_arrays | `(naky)` | 3988-3989; 4307-4309; restart 826-830 | yes (state) | E (initial) |
| `jump` (priv) | dist_fn | `(naky)` int | 4237-4238, 4302/4308 | yes | R each step |
| `ikx_indexed` (priv) | dist_fn | `(ntheta0)` int | 4240-4249 | no | R |
| `dkx`, `dtheta0` (priv, saved) | dist_fn | scalars | 4251-4265 | no | R |
| `should_use_kx` | — | logical | 4228 `box .or. (shat .eq. 0.0)` | no | R |
| `istep_last` (saved) | exb_shear local | int, init 0 | 4217, 4286 | — | R (JAX state) |
| `has_phi, has_apar, has_bpar` | run_parameters | logical | run_parameters.f90:351-353 | no | E |

### 2.5 Filters, parity, links

| Fortran name | Module | Shape | Computed at | dt | Rec. |
|---|---|---|---|---|---|
| `kwork_filter` | kt_grids | `(ntheta0,naky)` logical | false (kt_grids.f90:188-189); `(1,1)=.true.` if box (fields.f90:376); toggled inside `fm_populate` (fields_local.f90:2006-2076) | no | E (value after `init_fields`), formula B |
| `ieqzip` | run_parameters | `(ntheta0,naky)` logical | run_parameters.f90:259-272 | no | E (gate check) |
| `itleft, itright` | dist_fn | `(ntheta0,naky)` int | `compute_itleft_and_itright` 2441-2512 | no | E |
| `l_links, r_links` | dist_fn | `(ntheta0,naky)` int | `init_connected_bc` | no | E |
| parity mapping `par_it(it,ik)` | — (derived from `parity_redist`, priv) | `(ntheta0,naky)` int | `get_parity_conn` 4742-4793 | no | R (§3.6) |
| `def_parity, even, boundary_option_switch` | dist_fn | scalars | 1250, 1253, 1290-1310 | no | E |

### 2.6 Field-equation coefficients touched here

| Fortran name | Module | Shape | Computed at | dt | Rec. |
|---|---|---|---|---|---|
| `gamtot` | dist_fn (public) | `(-ntgrid:ntgrid,ntheta0,naky)` | 6523-6524 (+poisfac), 6550-6554 (esv), 6557-6569 (+tite) | no | **E** (final value) |
| `gamtot1, gamtot2` | dist_fn | same | 6536, 6547, esv 6552-6553 | no | E |
| `gamtot3` (priv) | dist_fn | same (fieldlineavg only) | 6564-6565 | no | R: `(gamtot - tite)/gamtot`, `<2eps → 1` from **exported final `gamtot`** |
| `awgt` (priv, lazily on first `calculate_flux_surface_average`) | dist_fn | `(ntheta0,naky)` | 8006-8020 | no | R (§3.9.3) |
| `fl_avg` (public) | dist_fn | `(ntheta0,naky)` complex | 8022-8034, per field solve | no | R |
| `tite` | run_parameters | scalar | run_parameters.f90:120 | no | E |
| `poisfac, esv, adiabatic_option_switch` | dist_fn | scalars | 1271, 1252+1327, 1315-1317 | no | E |
| `gridfac1` | dist_fn | same as gamtot | 6487-6506 | no | E (field spec) |

### 2.7 Species (`species.f90`)

| Fortran name | Shape | Computed at | dt | Rec. |
|---|---|---|---|---|
| `spec(is)%{z,mass,dens,temp,tprim,fprim,uprim,uprim2,vnewk,nu_h,bess_fac,vcrit,vcprim,type,f0type}` | `(nspec)` | species.f90:448-474 (overrides 991-995) | no | E |
| `spec(is)%{stm,zstm,tz,zt,smz}` | `(nspec)` | species.f90:476-480 / 1009-1013 | no | E (B) |
| `nonmaxw_corr, f0_values, dlnf0drho` | `(negrid,nspec)` → NumPy `(nspec,negrid)` | `calculate_f0_arrays` 563-597 | no | E |
| `wstar` | `(naky,negrid,nspec)` → NumPy `(nspec,negrid,naky)` | `init_wstar` dist_fn.fpp:1732-1750 | **yes** | R: `wstar(ik,ie,is) = -((code_dt*wunits(ik))*dlnf0drho(ie,is))` (source spec owns it; see P20) |
| `fphi, fbpar` | scalars | run_parameters.f90 | no | E |

### 2.8 Gate flags to export

`no_driver` (antenna, public), `gs2_has_lowflow` (build_config, public), `eqzip_option_switch`/`any(ieqzip)`,
`afilter` (via config), `nonlin`, `wstar_units`, `hyper` effective switch, `spec%type`, `spec%f0type`,
`field_option`, `split_collisions`, `collisions on/off`.

---

## 3. Algorithms

### 3.1 Shared numerical primitives (must be ported bit-exactly)

#### 3.1.1 `klein_sum` (`externals/utils/summation.f90:71-103`) — sequential, compensated

```python
def klein_sum(a):                 # a: 1-D, summed strictly i = 1..n in order
    s = 0.0; c1 = 0.0; c2 = 0.0
    for x in a:
        t = s + x
        if abs(s) >= abs(x): lc1 = (s - t) + x
        else:                lc1 = (x - t) + s
        s = t
        t2 = c1 + lc1
        if abs(c1) >= abs(lc1): lc2 = (c1 - t2) + lc1
        else:                   lc2 = (lc1 - t2) + c1
        c1 = t2
        c2 = c2 + lc2
    return (s + c1) + c2           # Fortran: the_result + compensation_1 + compensation_2
```
Complex version (`integration.f90:42-52`): `cmplx(klein_sum(real(a)), klein_sum(aimag(a)))`.
No multiplications → FMA-proof. JAX: `lax.scan` over the summation axis, vectorised over all other axes,
branches via `jnp.where`. XLA must not reassociate (default; never enable fast-math).

#### 3.1.2 `trapezoidal_integration(x, y)` (`integration.f90:55-71`)

```python
n = len(x)                         # n = nth = 2*ntgrid+1
integral = local_sum((x[1:] - x[:-1]) * (y[:-1] + y[1:])) * 0.5    # local_sum = klein_sum (use_compensated_sum = .true., line 24)
```
For complex `y`: `dx*(y1+y2)` is real×complex component-wise; final `*0.5` component-wise.

#### 3.1.3 `field_line_average(q)` (`theta_grid.f90:657-670`)

`trapezoidal_integration(theta, q*jacob) * field_line_average_weight`, with
`field_line_average_weight = 1.0/trapezoidal_integration(theta, jacob)` (theta_grid.f90:520).
`q*jacob` is complex×real (component-wise) for complex `q`.

#### 3.1.4 Fortran `NINT`

Round half **away from zero**: `nint(x) = int(sign(x) * floor(abs(x) + 0.5))`.
**Do not use `jnp.round`** (half-to-even). Used in `exb_shear` (4302, 4308).

#### 3.1.5 `real ** integer-variable` (e.g. `**nexp`, `**(2*nexp)`)

gfortran → `__builtin_powi` → libgcc `__powidf2`: `y = (n odd ? x : 1); while (n >>= 1) { x = x*x; if (n odd) y = y*x; }`.
`jax.lax.integer_pow(x, n)` (static `n`) performs the identical sequence of multiplications (`acc*x` with
`acc` first) → bit-identical. Constant integer powers (`**2`, `**4`) are expanded the same way
(`x*x`, `(x*x)*(x*x)`). `real ** real` (`**0.5`, `**(P_hyper3D/2.0)`) is a libm `pow` call (not folded to
`sqrt` without `-funsafe-math-optimizations`); XLA may lower `pow(x,0.5)` to `sqrt` → ≤1 ulp difference (P3).

---

### 3.2 Hyperviscosity / hyperresistivity (`hyper.f90`) — SUPPORT

There is **no "gyro" or separate "nonlinear" hyper option in 8.2.1**: the options are exactly
`default|none|visc_only|res_only|both|simple3D` (hyper.f90:458-464). The amplitude-dependent
(Belli–Hammett "nonlinear") damping is `const_amp = .false.`; `include_kpar` is read but unused (106, 481).

#### 3.2.1 Parameter resolution (replicate `read_parameters`, hyper.f90:452-601)

Input defaults (hyper.f90:63-126): `const_amp=F, d_hyper=-10, d_hyper3d=-10, d_hyperres=-10, d_hypervisc=-10,
damp_zonal_only=F, gridnorm=T, hyper_option='default', include_kpar=F, isotropic_model=T,
isotropic_shear=T, kx_cut=-10, ky_cut=-10, nexp=2, omega_osc=0.4, p_hyper3d=4`.

```
switch = {'default':NONE,'none':NONE,'visc_only':VISC,'res_only':RES,'both':BOTH,'simple3D':S3D}[hyper_option]
hyper_on = False                                   # module default (35); never reset except finish_hyper
if (not isotropic_shear) and nexp != 2: nexp = 2   # 498-501 (forced for aniso)
NONE: if D_hyperres > 0: D_hyperres = -10 ; if D_hypervisc > 0: D_hypervisc = -10        # 506-516
VISC: hyper_on = True ; if D_hyperres > 0: D_hyperres = -10 ; if D_hypervisc < 0: hyper_on = False   # 518-534
RES:  hyper_on = True ; if D_hyperres < 0: hyper_on = False ; if D_hypervisc > 0: D_hypervisc = -10  # 536-552
BOTH: hyper_on = True
      if D_hyper >= 0: D_hyperres = D_hyper ; D_hypervisc = D_hyper                 # 557-566 (`if (D_hyper < 0.)` else-branch)
      if D_hypervisc < 0 and D_hyperres < 0: hyper_on = False                          # 582
S3D:  hyper_on = True ; if P_hyper3D < 4: hyper_on = False ; if D_hyper3D < 0: hyper_on = False   # 584-599
```
(`check_hyper`, 141-302, is only called by `ingen`, never by `gs2`.) Note `D_hyper == 0` counts as "set"
(`D_hyper < 0.` test), giving `D_hyperres = D_hypervisc = 0` and `hyper_on = True`.

#### 3.2.2 `init_hyper` (hyper.f90:360-449), executed **once** with `dt0 := code_dt` at that moment

```
D_v = 0; D_eta = 0                                                       # 378-379
D_res(:,:) = 0 ; hypervisc_filter(:,:,:) = 1.0                            # allocate_arrays 604-617
if switch in (BOTH, RES, VISC):                                          # 382
    if gridnorm:                                                         # 386-391
        kxm = akx(ntheta0/2 + 1)            # Fortran integer division; 1-based index
        akx4_max    = powi(kxm, 2*nexp)
        aky_max     = aky(naky)
        aky4_max    = powi(aky(naky), 2*nexp)                            # unused by the step
        akperp4_max = powi(kxm*kxm + aky(naky)*aky(naky), nexp)
    else:  akx4_max = aky_max = aky4_max = akperp4_max = 1.0              # 392-395
    if D_hyperres > 0:                                                   # 400-411
        for ik in 1..naky: for it in 1..ntheta0:
            D_res(it,ik) = ((D_hyperres*dt0) * powi(aky(ik)*aky(ik) + akx(it)*akx(it), nexp)) / akperp4_max
        D_eta = D_hyperres/akperp4_max
    else: D_res = 0 ; D_eta = 0
    D_v = (D_hypervisc > 0) ? D_hypervisc/akperp4_max : 0.0               # 413-417
elif switch == S3D:                                                      # 419-445
    if not (ky_cut > 0): ky_cut = aky(naky)                              # 421
    if not (kx_cut > 0): kx_cut = akx((ntheta0+1)/2)                     # 422  NB different index from kxm above
    kperp2_max = max_ig( ((ky_cut*ky_cut)*gds2(ig)
                          + ((2.0*ky_cut)*kx_cut)*gds21(ig)/shat)
                          + (kx_cut*kx_cut)*gds22(ig)/(shat*shat) )      # 424-426 (left-to-right sum)
    kperp2_max_zonal = isotropic_model ? kperp2_max
                       : max_ig( (kx_cut*kx_cut)*gds22(ig)/(shat*shat) ) # 427-431
    ik = 1   (ZONAL BY INDEX, regardless of aky(1)):                     # 433-437
        for it: filter(:,it,1) = exp(-( (D_hyper3D*dt0) * pow(kperp2(:,it,1)/kperp2_max_zonal, P_hyper3D/2.0) ))
    for ik in 2..naky: for it:                                           # 439-444
        filter(:,it,ik) = exp(-( (D_hyper3D*dt0) * pow(kperp2(:,it,ik)/kperp2_max, P_hyper3D/2.0) ))
```
`aky(ik)**2 + akx(it)**2` → `aky*aky + akx*akx` (constant power 2). Note `init_hyper` does **not** look at
`hyper_on`, so `D_res`/`filter` are computed even if `hyper_on` is false; only `hyper_diff` checks it.

**Recommendation:** EXPORT `D_res` and (for S3D) `hypervisc_filter` exactly as left by `init_hyper`
(both carry `dt0`); recompute the scalars (`akx4_max`, `aky_max`, `akperp4_max`) with the formulas above.
**Never recompute `D_res`/S3D filter when dt changes** (§1.4) — use `dt0`.

#### 3.2.3 `hyper_diff(g0 = gnew, phi = phinew)` dispatcher (hyper.f90:620-669)

```
if not hyper_on: return                                            # 652
if switch in (BOTH, RES, VISC):
    if D_hypervisc < 0: return                                     # 657  (RES always returns here)
    iso_shear() if isotropic_shear else aniso_shear()              # 659-663
elif switch == S3D: simple3Dfilter()                               # 665-667
```

#### 3.2.4 `iso_shear` (hyper.f90:753-784) — uses current `code_dt`

```
if const_amp:  S(ig) = 1.0  for all ig                              # 757-758
else:                                                               # 759-771
    fac(1) = 1.0 ; fac(2:naky) = 0.5
    S(:) = 0.0
    for ik in 1..naky:                    # ORDER MATTERS: ik outer, it inner, sequential accumulation
        for it in 1..ntheta0:
            k2 = akx(it)*akx(it) + aky(ik)*aky(ik)
            S(:) = S(:) + ((real(conjg(phi(:,it,ik))*phi(:,it,ik)) * (k2*k2)) * fac(ik))
    S = pow(S, 0.5)
for every (ik,it) [Fortran: per local iglo, same value for all il,ie,is]:   # 773-782
    if damp_zonal_only and not (aky(ik) == epsilon(0.0)): skip       # 776  (see P12: skips EVERYTHING)
    q = aky(ik)*aky(ik) + akx(it)*akx(it)
    F(:,it,ik) = exp(-( (D_hypervisc*code_dt) * ((S(:) * powi(q, nexp)) / akperp4_max) ))
    g0(:,1,iglo) = g0(:,1,iglo) * F(:,it,ik)      # complex * real, component-wise
    g0(:,2,iglo) = g0(:,2,iglo) * F(:,it,ik)
```
`phi` is the **full replicated field** (no MPI reduction; deterministic, nproc-independent).

#### 3.2.5 `aniso_shear` (hyper.f90:688-750), `nexp == 2` forced

```
fac(1) = 1.0 ; fac(2:) = 0.5
Snz(:) = 0 ; for ik in 2..naky: for it in 1..ntheta0:                       # 699-705
    k2 = akx(it)*akx(it) + aky(ik)*aky(ik)
    Snz = Snz + ((|phi(:,it,ik)|^2 * (k2*k2)) * fac(ik))
Snz = 0.5 * ( -omega_osc + pow(omega_osc*omega_osc + 2*Snz, 0.5) )          # 706
Sz(:) = 0 ; ik = 1: for it in 1..ntheta0:                                   # 709-715
    k2 = akx(it)*akx(it) + aky(1)*aky(1)
    Sz = Sz + ((|phi(:,it,1)|^2 * (k2*k2)) * fac(1))
Sz = 0.5 * ( -omega_osc + pow(omega_osc*omega_osc + 2*Sz, 0.5) )            # 717
Sznz(:) = 0 ; for ik in 2..naky: for it in 1..ntheta0:                      # 720-726
    a2 = aky(ik)*aky(ik)
    Sznz = Sznz + ((|phi(:,it,ik)|^2 * (a2*a2)) * fac(ik))                   # aky**4 = (a*a)*(a*a)
Sznz = 0.5 * ( -omega_osc + pow(omega_osc*omega_osc + 2*Sznz, 0.5) )        # 728
for each (ik,it):                                                           # 733-748
    x2 = akx(it)*akx(it) ; x4 = x2*x2
    if aky(ik) == 0.0:                                                      # EXACT compare (737)
        arg = (D_hypervisc*code_dt) * ((Sznz * x4) / akx4_max)
    else:
        q  = aky(ik)*aky(ik) + akx(it)*akx(it)
        arg = (D_hypervisc*code_dt) * ( ((Snz * (q*q)) / akperp4_max)
                                        + ((((Sz * x4) / akx4_max) * aky(ik)) / aky_max) )
    F(:,it,ik) = exp(-arg) ; g0(:,1:2,iglo) *= F(:,it,ik)
```
`|phi|^2 := phi.re*phi.re + phi.im*phi.im`; `2*Snz` is `2.0*Snz` (exact). `(q*q)` is `powi(q,2)`.

#### 3.2.6 `simple3Dfilter` (hyper.f90:673-685) and `D_res` consumption

```
g0(:,isgn,iglo) = g0(:,isgn,iglo) * hypervisc_filter(:,it,ik)   for isgn = 1,2   (filter from init, dt0)
```
`D_res(it,ik)` enters the source term only (owned by the source spec; reproduced here for completeness):
non-opt `set_source` (dist_fn.fpp:5094-5102): `... - spec(is)%zstm*vpac*nonmaxw_corr*((aj0(ig+1)+aj0(ig))*0.5*apar_m + D_res(it,ik)*apar_p) ...`;
opt path `source_coeffs_aparp = -fapar*D_res(it,ik)*spec(is)%zstm*vpac*nonmaxw_corr - ...` (2252-2258).
`D_res = 0` unless `switch ∈ {RES, BOTH}` and effective `D_hyperres > 0`.

#### 3.2.7 Masks / special cases

* The filter multiplies **all** theta points, including forbidden-region points (g is 0 there unless
  `zero_forbid` is false and junk accumulates; multiply anyway).
* `kwork_filter`, `ieqzip`, passing-hybrid cells are **not** excluded by `hyper_diff` (§6 P10).
* For `(it,ik) = (1,1)` in box: iso `q = 0 → powi(0,nexp)=0 → F = exp(-0) = 1`; aniso `aky==0` branch with
  `akx(1)=0 → F=1`; S3D `kperp2=0 → pow(0,P/2)=0 → F=1`. Exact identity.
* Division hazards → GATE: `gridnorm` and (`akperp4_max == 0` or, for aniso, `akx4_max == 0` or
  `aky_max == 0`); S3D with `shat == 0` (`gds21/shat`), or `kperp2_max(_zonal) == 0`.

---

### 3.3 ExB flow shear (`exb_shear`, dist_fn.fpp:4183-4600) — SUPPORT

Called once per `advance`, **before** `g = gnew` and before the first `timeadv`, on `gnew, phinew, aparnew,
bparnew` (fields_implicit.f90:335). Active iff `abs(g_exb*g_exbfac) > epsilon(0.)` (3983 / fields_local 4280).
It is a pure **permutation + zero-fill** (no arithmetic on g or fields) except for the scalar bookkeeping
of `kx_shift`/`theta0_shift`. It does **not** shift the Adams–Bashforth history `gexp_1/2/3`, `g_h`, the
antenna arrays, or anything else (P6).

#### 3.3.1 One-time setup (first call, `exb_first`, 4235-4267)

```
should_use_kx = box .or. (shat .eq. 0.0)                   # 4228, exact equality
jump(:) = 0
if should_use_kx:
    itmin = minloc(ikx)                                     # first index of the minimum
    ikx_indexed(p) for p = 1..ntheta0:  [itmin, itmin+1, ..., ntheta0, 1, 2, ..., itmin-1]   # 4243-4249
    # box: ikx = [0,1,..,(nt-1)/2, -(nt-1)/2,..,-1] -> ikx_indexed lists it in increasing kx.
    # range/specified: ikx = it-1 -> identity. single: ikx=0 -> identity.
    if ntheta0 > 1: dkx = akx(2) - akx(1)                   # 4252
    else: error message only, dkx UNINITIALISED              # 4254 -> GATE (g_exb!=0, kx branch, ntheta0==1)
else:
    dtheta0 = (ntheta0 > 1) ? theta0(2,1) - theta0(1,1) : twopi    # 4261-4265
```

#### 3.3.2 Per-call gating and shift bookkeeping (4270-4314)

```
if istep == istep_last: return                              # 4270 (istep_last saved, initial 0)
if g_exb_start_timestep > istep: return                     # 4271 (default -1: never)
if g_exb_start_time >= 0 and code_time < g_exb_start_time: return   # 4272
gdt = 0.5*(code_dt + code_dt_old)                           # 4283  (first step: code_dt_old = 0 -> dt/2)
istep_last = istep                                          # 4286  (only updated when active!)
if should_use_kx:
    for ik in 1..naky:                                      # 4300-4304
        kx_shift(ik) = kx_shift(ik) - ((((aky(ik)*g_exb)*g_exbfac)*gdt)*tunits(ik))
        jump(ik)     = NINT(kx_shift(ik)/dkx)
        kx_shift(ik) = kx_shift(ik) - real(jump(ik))*dkx
else:
    for ik in 1..naky:                                      # 4306-4310
        theta0_shift(ik) = theta0_shift(ik) - ((((g_exb*g_exbfac)*gdt)/shat)*tunits(ik))
        jump(ik)         = NINT(theta0_shift(ik)/dtheta0)
        theta0_shift(ik) = theta0_shift(ik) - dtheta0*real(jump(ik))
if all(jump == 0): return                                   # 4314
```
`code_time` here is the time at the *start* of the step (`update_time` runs after `advance`).
The bookkeeping depends only on dt history, `aky`, `tunits` and inputs — never on g — so it can be done on the
host in NumPy float64 with the same operation order (recommended), carrying `kx_shift/theta0_shift/
istep_last/code_dt_old` as Python state and passing `jump` (int32) into the jitted step.

#### 3.3.3 kx branch (`should_use_kx`, box or shat==0) — 4490-4599

For each `ik = naky, naky-1, ..., 2` (**ik = 1 is never shifted**, even if `aky(1) ≠ 0`), with `J = jump(ik)`
and sorted positions `p` (`it = ikx_indexed(p)`):

```
J < 0: for p = 1 .. ntheta0+J        : X(:, ikx_indexed(p), ik) = X(:, ikx_indexed(p-J), ik)   (ascending)
       for p = ntheta0+J+1 .. ntheta0 : X(:, ikx_indexed(p), ik) = 0
J > 0: for p = ntheta0 .. 1+J        : X(:, ikx_indexed(p), ik) = X(:, ikx_indexed(p-J), ik)   (descending)
       for p = J .. 1                 : X(:, ikx_indexed(p), ik) = 0
```
applied to `phi` (if `has_phi`), `apar` (if `has_apar`), `bpar` (if `has_bpar`) (4494-4517, 4547-4570), and to
`g0(:,:,iglo)` for every `(il, ie, is)`, both `isgn` (4518-4543, 4571-4596). The loop directions make each
in-place loop equivalent to a **simultaneous gather**:

```
q = p - J ;  valid = (1 <= q <= ntheta0)
X_new[ik, ikx_indexed(p)] = valid ? X_old[ik, ikx_indexed(q)] : 0        (for ik >= 2; ik = 1 unchanged)
```
If `|J| ≥ ntheta0` everything for that ik is zeroed (no abort in this branch).

#### 3.3.4 theta0 branch (ballooning space, finite shat, not box) — 4341-4489

Per `ik = 1..naky` **in increasing order**, `j = jump(ik)`; skip if `j == 0`;
`if abs(j) >= ntheta0: mp_abort` (4347-4351, host check).
`m = |j|, iit = ntheta0+1-m, iib = m, ileft = -ntgrid+ntheta, iright = ntgrid-ntheta` (4353-4354).
Unified simultaneous-gather form (Fortran 1-based `i`, `c := i - j`):

```
if 1 <= c <= ntheta0:      X_new(ig, i) = X_old(ig, c)
elif c < 1   (j > 0):      X_new(ig, i) = (ig <= iright) ? X_old(ig + ntheta, c + ntheta0) : 0
elif c > nt  (j < 0):      X_new(ig, i) = (ig >= ileft)  ? X_old(ig - ntheta, c - ntheta0) : 0
```
(derivation: j<0: 4357-4363 / 4414-4447 copy `temp2(:iright,:)` into `(ileft:, iit:)`; j>0: 4364-4370 /
4449-4484 copy `temp2(ileft:,:)` into `(:iright, :iib)`; i.e. a wrap in theta0 by one 2π period is a shift of
`ntheta` theta points.) For `g` this is applied per `(il,ie,is)` and both `isgn` with no cross-ik effects.

**Field-only bug that must be replicated (4370, 4387, 4404):** in the `j > 0` branch the zeroing statement
is `phi(iright+1:ntgrid, :iib, :) = 0.0` — the last index is `:` (**all ky**), not `ik`. So after processing
ky index `ik` with `j>0`, rows `ig > iright` of columns `it = 1..m` are zeroed for **every** ky (including
ones processed later, whose shift then reads those zeros into columns `m+1..2m`, and ones with `j = 0`).
Faithful algorithm for each field `X ∈ {phi (has_phi), apar (has_apar), bpar (has_bpar)}`:

```
for ik in 1..naky (sequential):
    j = jump(ik); if j == 0: continue
    X[:, :, ik] = theta0_gather(X[:, :, ik], j)     # unified formula above
    if j > 0: X[iright+1:ntgrid, 1:j, 1:naky] = 0    # ALL ik
```
With `wstar_units = F`, `tunits = 1` so all `jump(ik)` are equal; with `wstar_units = T` they differ.
Note the physical meaning requires the theta0 grid to span exactly 2π (`ntheta0*dtheta0 = 2π`); GS2 does not
check this (P8).

#### 3.3.5 Restart

`kx_shift`, `theta0_shift` are written/read in the restart file (`gs2_save.fpp:428-450`, `820-830`);
`code_dt_old` is **not** saved → first step of a restarted run again uses `gdt = code_dt/2`.

---

### 3.4 `kwork_filter` (kt_grids) — SUPPORT

* Value during time stepping: all `.false.` except `kwork_filter(1,1) = .true.` when
  `gridopt_switch == gridopt_box` (fields.f90:376, set at the end of `init_fields`, after response matrices).
* Read (skip `(it,ik)`) in: source term (`dist_fn.fpp:5401` → `g_work` stays 0 there), `invert_rhs_1`
  (5475 → **`gnew` for that cell is not touched at all: neither zeroed nor solved**, so it keeps its value from
  the start of the step = `g`), linked BCs (6000), `getfieldeq1_nogath` (6750…6823), `getan_nogath`,
  `integrate_species_original` (le_grids.f90:744-761: filtered cells contribute 0 / stay 0), collisions
  (collisions.fpp many sites), `calculate_flux_surface_average` (8016, 8030), fields_local / gf_local loops.
* Not read by: `hyper_diff` (filter is exactly 1 there anyway), `exb_shear` (ik=1 never shifted),
  `enforce_parity`, `remove_zonal_flows` (it does subtract the field-line average of `phinew(:,1,1)`).
* JAX: boolean mask `KW[ik,it]`; after `invert_rhs_1`: `gnew = where(KW, gnew_in, gnew_solved)`; source:
  `where(KW, 0, source)`; same mask in the other consumers listed.
* Transient states: during `fm_populate` (fields_local.f90:2006, 2032, 2058, 2076, 2114, 2125) the filter is
  reset/raised per response-matrix pass and reset to `.false.` at the end; for implicit fields after a dt
  reset, the rebuild runs **with `(1,1)` already `.true.`** (P11).

### 3.5 `ieqzip` / `eqzip_option` (run_parameters) — GATE

`ieqzip(ntheta0,naky)`: `'secondary' → (1,2)`, `'tertiary' → (2,1),(ntheta0,1)`, `'equilibrium' → (1:ntheta0,1)`
(run_parameters.f90:259-272; note the docstring wrongly says "harris"). Effect: `invert_rhs_1` skips the cell
(5477: `gnew` frozen), collisions skip it, **but** the source is still computed, linked BCs still modify it
(no ieqzip test in 5996-6010), `hyper_diff` is applied in **both** `timeadv` calls (so a frozen cell decays by
`F²` per step), and the field solve still solves it. Too many interacting special cases for a research option
(secondary-instability analysis) → **GATE: refuse if `eqzip_option /= 'none'` (any(ieqzip))**.

### 3.6 `def_parity` / `even` / `enforce_parity` (dist_fn.fpp:4603-4853) — SUPPORT

`mult = even ? +1 : -1` (integer; `mult*z` is exact). Called at the end of `invert_rhs_1` (5706, **before**
linked BCs) and in `timeadv` after collisions (4091). Operates on `gnew` only.

**Non-linked (ballooning / zero / periodic)** (4832-4851):
```
for every iglo, for ig in 1..ntgrid:
    gnew(-ig, 1) = mult * gnew( ig, 2)      # 4848
    gnew( ig, 1) = mult * gnew(-ig, 2)      # 4849
# ig = 0 of isgn=1 is NOT touched; isgn=2 never touched.
```
**Linked** (4819-4831): redistribute built in `init_enforce_parity` (4603-4738):
`gnew(-ig, 1, iglo) = gnew(ig, 2, iglo_conn)` for **all** `ig ∈ [-ntgrid, ntgrid]` (including 0), with
`iglo_conn = idx(ik, par_it(it,ik), il, ie, is)`; then, if `mult ≠ 1`, `gnew(:,1,:) = mult*gnew(:,1,:)`.
Source (isgn=2) and destination (isgn=1) never overlap → simultaneous gather.

Parity partner `par_it(it,ik)` (`get_parity_conn` 4742-4793), recompute from exported links:
```
L = l_links(it,ik); R = r_links(it,ik)
if L == R: par_it = it
elif L > R: walk itleft from it, k steps, stopping when r_links(it_k, ik) == L  (k = L-R)
else:       walk itright from it, k steps, stopping when l_links(it_k, ik) == R  (k = R-L)
```
i.e. the mirror cell of the supercell chain. JAX: static gather
`G1_new[..., ik, it, ig'] = mult * G2[..., ik, par_it[ik,it], nth-1-ig']` (linked: all `ig'`;
non-linked: `par_it = it`, all `ig' ≠ ntgrid`, keep `G1[..., ntgrid]`).
Note: with up–down asymmetric geometry, `forbid(θ) ≠ forbid(-θ)` and parity can write non-zero values into
forbidden points — replicate, do not mask.

### 3.7 `wstar_units` — SUPPORT (linear only)

`gs2_time.f90:108-124`: `wunits = 1; tunits = 2/aky (aky≠0) else 0` vs default `tunits = 1; wunits = aky/2`.
It enters every dt-dependent coefficient via `tunits`/`wunits` (vpar 1716, wdrift 1528, coriolis 1554,
wstar 1746, source 5100-5101, exb 4301/4307). EXPORT `tunits`, `wunits` and use them everywhere.
GATE: `wstar_units .and. nonlin` (nonlinear_terms.fpp:169-172 says it is required false) and
`wstar_units .and. any(aky == 0)` (tunits = 0 → "garbage results", gs2_time.f90:115-120).

### 3.8 Hybrid electrons — SUPPORT (dist_fn part)

Species `type = 'hybrid_electron'` (=4, species.f90:40). Electron-like (`has_electron_species` true, so no
adiabatic electron response). Mask (le_grids.f90:4490-4499):
```
HYB(is, ik, il) = (spec(is)%type == 4) .and. (aky(ik) /= 0.0) .and. (il <= ng2 + 1)   # passing OR wfb
```
In `invert_rhs_1` (5487-5490), for such cells (after the kwork/ieqzip skips), `gnew` is **overwritten**
by `set_h_zero(gnew, phinew, bparnew, iglo)` (dist_fn_arrays.fpp:246-294), i.e. h = 0:
```
phi_factor  = (fphi /= 0)  ? (spec(is)%zt*nonmaxw_corr(ie,is))*(-fphi) : 0.0
bpar_factor = (fbpar /= 0) ? 2.0*(-fbpar) : 0.0
for ig in -ntgrid..ntgrid:
    adj = ((bpar_factor*vperp2(ig))*aj1(ig))*bpar(ig,it,ik) + ((phi_factor*phi(ig,it,ik))*aj0(ig))
    gnew(ig,1) = adj ; gnew(ig,2) = adj
```
(`phi`, `bpar` = the `phinew`, `bparnew` passed to `timeadv`.) These cells are skipped by linked BCs (6005),
but **are** subsequently multiplied by the hyper filter, touched by parity enforcement, collisions etc.
JAX: `gnew = where(HYB[is,:,il,ik], adj, gnew_solved)` inside `invert_rhs_1`, before linked BCs.
Other hybrid hooks (owned by the collisions spec; must be ported there or GATE `hybrid + collisions`):
`zero_out_passing_hybrid_electrons` (collisions.fpp:1683-1703), `set_hzero_lorentz_collisions_matrix`
(2265, 2334, 2370-2391). Initial state: `init_g.f90:1236-1264` zeroes passing hybrid `gnew` (exported state).
GS2 warns that EM + hybrid is unphysical (run_parameters.f90:205-213) — allowed, not gated.

### 3.9 Adiabatic response — SUPPORT

#### 3.9.1 Detection (species.f90:914-959)

`has_electron_species = any(type ∈ {2 electron, 4 hybrid})`; `has_ion_species = any(type == 1)` (tracer=3 is
neither). Two **different** conditions are used:

* `init_fieldeq` adds the adiabatic term if `(.not. has_electron) .or. (.not. has_ion)` (6557) — i.e. adiabatic
  electrons (ITG) *or* adiabatic ions (ETG); `tite` is "n q²/T of the adiabatic species" (run_parameters.f90:118-120).
* `fl_avg` (field-line average) is used only if `.not. has_electron` **and** `adiabatic_option_switch ==
  fieldlineavg` (6605-6609, 6734-6738, 7637, fields_local.f90:3392-3396).
* `calculate_potentials_from_nonadiabatic_dfn` uses `phi_weight = 1/(Σ dens z²/temp + tite)` under the
  first condition and **ignores fl_avg, poisfac and yavg** (7360-7366).

#### 3.9.2 `init_fieldeq` adiabatic part (6556-6569) — EXPORT final `gamtot`

```
gamtot = real(tot) + kperp2*poisfac                            # 6524 (tot from integrate_species, field spec)
if esv: ensure_single_val_fields_pass_r(gamtot/gamtot1/gamtot2)  # 6550-6554 (local it/ik range! P15)
if (.not.has_e) .or. (.not.has_i):
    yavg ('iphi00=3'):   for ik: if aky(ik) > epsilon(0.0): gamtot(:,:,ik) += tite
    fieldlineavg:        gamtot = gamtot + tite ; gamtot3 = (gamtot - tite)/gamtot ; where gamtot3 < 2*eps: gamtot3 = 1.0
    default:             gamtot = gamtot + tite
```
Option strings (dist_fn.fpp:1205-1213): `default, no-field-line-average-term, iphi00=0, iphi00=1` → default;
`field-line-average-term, iphi00=2` → fieldlineavg; `iphi00=3` → yavg.
`gamtot3` must be recomputed from the **final** exported `gamtot` as `(gamtot - tite)/gamtot` (not from the
pre-tite value: `(a+t)-t ≠ a` in floating point).

#### 3.9.3 `calculate_flux_surface_average(fl_avg, antot)` (7994-8035)

```
fl_avg(:,:) = 0
if awgt not allocated (first call ever):                       # 8006-8020
    awgt = 0
    for ik: if aky(ik) > epsilon(0.0): continue
        for it: if kwork_filter(it,ik) [value at that first call]: continue
            awgt(it,ik) = 1.0 / trapz_klein(theta, jacob*gamtot3(:,it,ik))
for ik in 1..naky: for it in 1..ntheta0:                       # 8027-8033
    if kwork_filter(it,ik): continue
    y(ig) = (jacob(ig)*antot(ig,it,ik)) / gamtot(ig,it,ik)       # real*complex then complex/real, component-wise
    fl_avg(it,ik) = (tite * trapz_klein(theta, y)) * awgt(it,ik)
```
The first call happens during response-matrix construction (getfieldeq with all-false `kwork_filter`) or in
the initial `get_init_field`; either way only `awgt(1,1)` could depend on it and `fl_avg(1,1)` is skipped
(box) during stepping. Recompute `awgt` with all-false filter semantics (P16).
Consumers add `fl_avg(it,ik)` (constant along theta) to `fieldeq` (6620-6623; nogath 6751-6770; fields_local
3422-3426 with a **different operand order**) or to the `phi` numerator `antot + fl_avg`
(`get_fields_direct_from_dfn` 7655/7670). These sums are owned by the field-solve spec.

### 3.10 `poisfac` — SUPPORT

Only effect in the time step: `gamtot += kperp2*poisfac` (6524) before esv/tite. Covered by exporting final
`gamtot`. `calculate_potentials_from_nonadiabatic_dfn` ignores it (comment 7361) — replicate (field spec).

### 3.11 `esv` — SUPPORT

Effective `esv = esv_input .and. (boundary_option_switch == linked)` (1327). `ensure_single_val_fields_pass`
(6201-6227, complex) / `_r` (6235-6261, real):
```
for ik in ik_min..ik_max: for it in it_min..it_max:      # force_all -> 1..naky, 1..ntheta0; else g_lo local range
    lit = itright(it,ik); if lit < 0: continue
    arr(-ntgrid, lit, ik) = arr(ntgrid, it, ik)
```
Reads `ig = +ntgrid`, writes `ig = -ntgrid`, and `itright` is injective (ky=0 cells map to themselves) →
simultaneous scatter: `A[ik, itright[ik,it], 0] = A[ik, it, nth-1]` where `itright ≥ 0`.
Used on `antot` (6317), `antota` (6337), `antotp` (6357) with `force_all=.true.` (after `afilter`), in the
nogath variants (7087, 7164, 7245) **without** `force_all`, and on `gamtot/1/2` at init (6550-6554) **without**
`force_all` (P15).

### 3.12 `remove_zonal_flows` (fields_implicit.f90:374-425) — SUPPORT

After `phinew = phinew + phi` (step 9 in §1.2), only if `fields_knobs:remove_zonal_flows_switch` (default F):
```
phi_avg(:,:,:) = 0
for it in 1..ntheta0: phi_avg(:, it, 1) = field_line_average(phinew(:, it, 1))     # ik = 1 BY INDEX
phinew = phinew - phi_avg                                                         # whole array
```
Only `phinew`; `aparnew`, `bparnew`, `gnew` untouched (gnew is then recomputed by the 2nd `timeadv` with the
modified `phinew`). `field_line_average` = §3.1.3 (Klein-compensated trapezoid). gf_local prints a warning
but applies it.

### 3.13 Antenna (`antenna.f90`, `apar_ext`) — GATE

What it is: Langevin antenna (Hammett/Dorland/Quataert) driving `apar` at `nk_stir` modes for Alfvén-wave
turbulence studies. Active iff namelist `driver` exists and `ant_off = .false.`
(`no_driver = ant_off .or. .not.exist`, antenna.f90:259). Per step (`antenna_amplitudes` 341-443):
`w_stir = gradpar(0)*|kz|*sqrt(1/beta_s)*wtmp`; on proc0, random forcing
`force = cmplx(ranf()-0.5, ranf()-0.5)*sqrt(12*|Im w|/dt)*amplitude` (Mersenne-Twister `ranf`,
`externals/utils/ran.fpp:21-46`, state shared with `init_g`); `a_ant = a_ant*exp(-i w dt) + force*dt` etc.;
`apar(:,it,ik) += (a+b)/sqrt(2)*exp(i kz theta)`; reality fix for ky=0; added to `aparnew` before the field
solve (fields_implicit.f90:345). Reasons to gate: RNG stream replication, dt- and reset-dependent RNG
consumption (called at the top of every `advance`, including failed attempts), not supported for gf_local,
not a standard flux-tube feature. **Export program must abort if `no_driver == .false.`.**

### 3.14 `lowflow.f90` — GATE

What it is: neoclassical "low-flow" (O(ρ*)) corrections from NEO output files (H^neo, φ^neo and their
derivatives) entering as modified drifts (`wdrift` becomes **complex**, dist_fn.fpp:132-138), `wstar_neo`,
`vparterm`, `wdfac`, `wstarfac`, `hneoc`, `wdttpfac` (lowflow.f90:57-327) and altered source/invert_rhs code
paths (`#ifdef LOWFLOW` throughout dist_fn.fpp, nonlinear_terms, le_derivatives, collisions). It is a
**compile-time** option (`-DLOWFLOW`, Makefile:471-475), **deprecated in 8.2** ("will be removed"), and the
reference build has `Lowflow: F`. Related inputs `lf_default`, `lf_decompose`, `neo_test`, `rhostar` are inert
in a non-LOWFLOW build. **GATE: export program aborts if `build_config::gs2_has_lowflow`.**

### 3.15 Species (`species.f90`) — SUPPORT

* Types (37-43): `ion=1`, `electron=2`, `trace=3` (zero weight in all velocity integrals:
  `where (spec%type == tracer_species) weights = 0`, le_grids.f90:741/814/901/961/1061 — it is still evolved),
  `hybrid_electron=4`. Type also selects collision-operator treatment (collisions spec).
* Derived (476-480, identical in 1009-1013): `stm = sqrt(temp/mass)`, `zstm = z/sqrt(temp*mass)`,
  `tz = temp/z`, `zt = z/temp`, `smz = abs(sqrt(temp*mass)/z)`. EXPORT (all are IEEE-exact to recompute too).
* Background `f0type` (563-597): Maxwellian → `nonmaxw_corr = 1`, `dlnf0drho = -fprim - (E-1.5)*tprim`,
  `f0_values = exp(-v²)/π^1.5`; `sdanalytic`/`tabulated` fill the same arrays (and change the energy grid,
  egrid.f90:28-38). All downstream use is through these arrays → EXPORT them and **always multiply by
  `nonmaxw_corr` exactly where GS2 does even when it is 1.0** (×1.0 is exact, but omitting it changes nothing
  only if done identically). Collisions use `energy_maxw`/`w_maxw` (collisions spec decides whether to gate
  non-Maxwellian + collisions).
* Adiabatic species are *not* species entries: they are implied by the absence of electrons (or ions), with
  strength `tite` (§3.9). `me`, `zi_fac` only feed `sdanalytic` parameters (484-560).

### 3.16 `afilter` (dist_fn_knobs, "for debugging only") — GATE

`if (afilter > epsilon(0.0)) antot = antot * exp(-afilter**4*kperp2**2/4.)` (6316, 7126), before esv. Changes
every `getan` variant and the response matrices; debug-only → refuse if `afilter > epsilon(0.0)`.

---

## 4. Options (namelist knobs) that change behaviour here

| Name | Namelist | Default | Decision | Effect / reason |
|---|---|---|---|---|
| `hyper_option` | hyper_knobs | `'default'` | SUPPORT | none/visc_only/res_only/both/simple3D (§3.2) |
| `d_hypervisc` | hyper_knobs | -10 | SUPPORT | viscosity coefficient; <0 disables |
| `d_hyperres` | hyper_knobs | -10 | SUPPORT | hyperresistivity → `D_res` (frozen at dt0) |
| `d_hyper` | hyper_knobs | -10 | SUPPORT | `both`: if ≥0 sets both |
| `const_amp` | hyper_knobs | F | SUPPORT | F: Belli–Hammett amplitude-dependent shear rate |
| `isotropic_shear` | hyper_knobs | T | SUPPORT | F: anisotropic zonal/non-zonal model, forces nexp=2 |
| `omega_osc` | hyper_knobs | 0.4 | SUPPORT | aniso only |
| `nexp` | hyper_knobs | 2 | SUPPORT | power of k⊥²; use `lax.integer_pow` |
| `gridnorm` | hyper_knobs | T | SUPPORT | normalise by grid max k (GATE if a normaliser is 0) |
| `damp_zonal_only` | hyper_knobs | F | SUPPORT | replicate bug: T ⇒ iso filter applied to no mode |
| `include_kpar` | hyper_knobs | F | SUPPORT | unused (no-op) |
| `d_hyper3d` | hyper_knobs | -10 | SUPPORT | simple3D coefficient (<0 disables) |
| `p_hyper3d` | hyper_knobs | 4 | SUPPORT | simple3D power (<4 disables) |
| `ky_cut`, `kx_cut` | hyper_knobs | -10 | SUPPORT | simple3D cut-offs (≤0 → grid max) |
| `isotropic_model` | hyper_knobs | T | SUPPORT | simple3D zonal normalisation |
| `hyper_in_initialisation` | dist_fn_knobs | T | SUPPORT | hyper inside response-matrix build (istep=0) |
| `g_exb` | dist_fn_knobs | 0.0 | SUPPORT | flow shear rate (perp shift × g_exbfac; parallel drive in source) |
| `g_exbfac` | dist_fn_knobs | 1.0 | SUPPORT | multiplies perp shift only |
| `g_exb_start_time` | dist_fn_knobs | -1 | SUPPORT | delay activation (time) |
| `g_exb_start_timestep` | dist_fn_knobs | -1 | SUPPORT | delay activation (step) |
| `g_exb_error_limit` | dist_fn_knobs | 0.0 | SUPPORT | read but never used (no-op) |
| `omprimfac`, `btor_slab`, `mach` | dist_fn_knobs | 1, 0, 0 | SUPPORT | source / geometry / drift specs (cross-ref) |
| `def_parity` | dist_fn_knobs | F | SUPPORT | enforce parity (§3.6); GS2 warns for nonlinear |
| `even` | dist_fn_knobs | T | SUPPORT | parity sign |
| `esv` | dist_fn_knobs | F | SUPPORT | single-valued antot / gamtot (linked only) |
| `poisfac` | dist_fn_knobs | 0.0 | SUPPORT | via exported gamtot |
| `adiabatic_option` | dist_fn_knobs | `'default'` | SUPPORT | default / fieldlineavg / yavg |
| `afilter` | dist_fn_knobs | 0.0 | GATE | debug only |
| `lf_default`, `lf_decompose` | dist_fn_knobs | T, F | SUPPORT (inert) | only meaningful in LOWFLOW builds |
| `wstar_units` | knobs | F | SUPPORT | linear only; GATE with nonlin or any aky==0 |
| `eqzip_option` | knobs | `'none'` | GATE | frozen-mode research option (§3.5) |
| `tite` | knobs | 1.0 | SUPPORT | adiabatic species strength |
| `immediate_reset` | knobs | T | SUPPORT | retry semantics: exb not re-applied on retry (§1.4) |
| `seed` | knobs | 0 | SUPPORT (inert) | RNG only used by init_g (exported state) and antenna (gated) |
| `neo_test`, `rhostar` | knobs | F, 3e-3 | GATE / inert | lowflow only |
| `remove_zonal_flows_switch` | fields_knobs | F | SUPPORT | §3.12 |
| `driver` namelist / `ant_off` | driver | absent / F | GATE | antenna active iff present and ant_off=F |
| `amplitude, w_antenna, w_dot, t0, nk_stir, restarting, write_antenna` | driver | 0, (1,0), 0, -1, 1, F, F | GATE | antenna |
| `stir_N: kx, ky, kz, travel, a, b` | stir | 1,1,1,T,-1,-1 | GATE | antenna |
| `nspec` | species_knobs | 2 | SUPPORT | |
| `me`, `zi_fac` | species_knobs | -1, -1 | SUPPORT | sdanalytic f0 only (exported arrays) |
| `type` | species_parameters_N | `'default'`(ion) | SUPPORT | ion/electron/e/trace/hybrid_electron |
| `f0type` | species_parameters_N | `'default'` | SUPPORT | via exported nonmaxw_corr/dlnf0drho/f0_values/energy grid; collisions spec may gate |
| `z, mass, dens, temp, tprim, fprim, uprim, uprim2, vnewk, nu_h, bess_fac, vcrit, vcprim` | species_parameters_N | 1,1,1,1,6.9,2.2,0,0,0,0,1,-1,-999.9 | SUPPORT | exported |
| `dens0, u0, tpar0, tperp0` | species_parameters_N | 1,1,0,0 | SUPPORT (inert) | init_g only |
| `LOWFLOW` (build macro) | Makefile | off | GATE | §3.14 |

---

## 5. GPU mapping notes

* **Hyper filter (iso/aniso)**: shear rates are reductions over `(ik,it)` for each `ig`. Faithful mode:
  `lax.fori_loop`/`lax.scan` over the flattened `k = (ik-1)*ntheta0 + (it-1)` (= C-order flatten of
  `PHI[ik,it,:]`), carrying `S[nth]` — sequential depth `naky*ntheta0`. Fast mode: `jnp.sum` (tree) — differs by
  ≲ `naky*ntheta0*eps` relative in S. A Pallas kernel (one thread per `ig`, sequential loop over k) gives
  faithful order at full speed. The filter `F[ik,it,ig]` is then a vectorised elementwise expression;
  apply as `G *= F[None,None,None,:,:,None,:]` component-wise (both `isgn`).
* **simple3D / D_res**: precomputed constants (exported); elementwise multiply / pass to source.
* **ExB kx branch**: per-ky dynamic gather along the `it` axis. Static: `ikx_indexed[p]`, inverse
  `pos[it]`. Dynamic (traced int32): `jump[ik]`. `q = pos[it] - jump[ik]`; `valid = (0<=q<nt) | (ik==0)`;
  `src = where(ik==0, it, ikx_indexed[clip(q)])`; `X_new = where(valid, take_along_axis(X, src, axis=it), 0)`.
  Same index arrays for `G` (axis 4) and fields (axis 1). Cheap; no host sync needed if `jump` is passed in.
* **ExB theta0 branch**: per ky gather along `it` plus a static `±ntheta` shift along `ig` selected by
  `s = (c<1) ? -1 : (c>nt) ? +1 : 0`: compute `G1 = take(G, wrap(c))`, `Gp = shift_ig(G1, +ntheta)`,
  `Gm = shift_ig(G1, -ntheta)` (static slices + zero pad), select by `s`. Fields: Python loop over static
  `naky` (or `lax.fori_loop`) to replicate the all-ky zeroing bug sequentially.
* Precompute `kx_shift/theta0_shift/jump` on the host (NumPy float64, same op order, Fortran NINT) each step;
  host check `|jump| < ntheta0` (theta0 branch abort).
* **kwork_filter / hybrid**: boolean masks broadcast over `G`; `jnp.where` — no irregular structure.
* **Parity**: static index arrays (`par_it[ik,it]`, theta reversal `nth-1-ig'`, isgn swap); one gather.
* **esv**: static scatter along the boundary (`.at[ik, itright, 0].set(A[ik, it, -1])` with the `itright<0`
  entries dropped — use a validity mask and a dummy index).
* **fl_avg / remove_zonal / awgt / field_line_average_weight**: Klein sum = `lax.scan` over the `nth-1`
  trapezoid terms, vectorised over `(naky, ntheta0)` and real/imag. Depth ~ 2*ntgrid — cheap.
* Everything here is independent across `(is,ie,il)`; the only sequential pieces are the shear-rate sum
  (faithful mode), the Klein sums (theta), and the theta0-branch field loop over ky.
* No MPI reductions exist in this subsystem (phi replicated; exb comms are pure copies); results are
  nproc-independent **except** esv on `gamtot` (P15).

---

## 6. Pitfalls (round-off and semantics)

* **P1 FMA / libm.** Reference binary: gfortran `-O3` on aarch64 with GCC default `-ffp-contract=fast`, so
  `a*b+c`, `kx_shift - prod*tunits`, `zr*zr+zi*zi`, etc. may be fused; `exp`, `pow` from Apple libm vs CUDA
  libdevice differ by ≤1 ulp. Bitwise equality is impossible; 1e-12 is fine. For tighter validation rebuild
  the reference with `-ffp-contract=off` (open question).
* **P2 NINT** is half-away-from-zero; `jnp.round` is half-to-even (`exb_shear` 4302/4308).
* **P3 `**0.5`** is `pow(x,0.5)` in Fortran; XLA may rewrite to `sqrt` (≤1 ulp).
* **P4 `code_dt_old = 0`** at program start and after restart ⇒ first active ExB step uses `gdt = dt/2`;
  after a dt reset `gdt = (new+old)/2`. Not saved in restart files.
* **P5 ExB once per istep.** `istep_last` blocks re-application on the reset retry; the shift (and the
  `kx_shift` increment with the pre-reset `gdt`) happened in the failed attempt and is carried by the
  initial-value overrides. `istep_last` is updated only when the shear is active (start-time gating).
* **P6 ExB does not shift NL history** `gexp_1/2/3` (nor anything else besides `gnew` and the three fields).
* **P7 kx branch never shifts ik = 1**, even if `aky(1) ≠ 0` (range/specified grids with `shat == 0`), but
  `kx_shift(1)`/`jump(1)` are still updated and `jump(1)` participates in `all(jump==0)`.
* **P8 theta0 branch**: field zeroing bug across all ky (§3.3.4); abort if `|jump| ≥ ntheta0`; physics assumes
  `ntheta0*dtheta0 = 2π` (unchecked); `should_use_kx` uses exact `shat .eq. 0.0` while boundaries use
  `is_effectively_zero_shear()` — a tiny non-zero shat gives periodic BCs *and* the theta0 branch.
* **P9 D_res and simple3D filter are frozen at the initial dt** (hyper level 19 < override_timestep 29);
  iso/aniso use the current dt.
* **P10 Hyper in both timeadv calls.** Normal cells are re-solved so only the 2nd application survives;
  cells skipped by `kwork_filter`/`ieqzip` keep `gnew` and get the filter twice per step (exactly 1.0 for
  box (1,1); ≠1 for ieqzip cells → another reason ieqzip is gated); passing-hybrid cells are overwritten by
  `set_h_zero` then filtered.
* **P11 Response-matrix rebuild state.** During construction `timeadv(istep=0)` applies hyper
  (`hyper_in_initialisation`) with an amplitude-dependent filter from the unit perturbation; for implicit
  fields a rebuild after a dt reset sees `kwork_filter(1,1) = .true.`, while fields_local resets it to false
  during population. Relevant only if JAX rebuilds response matrices itself.
* **P12 `damp_zonal_only` bug** (hyper.f90:776): test is `aky(ik) == epsilon(0.0)` (exact 2.22e-16), so with
  `damp_zonal_only = T` the iso filter is applied to **no** mode (S still computed).
* **P13 Different max-k indices**: `akx(ntheta0/2+1)` for visc/res normalisation vs `akx((ntheta0+1)/2)` for
  simple3D `kx_cut`; simple3D "zonal" means `ik == 1` by index.
* **P14 Exact-comparison zoo**: `aky(ik) == 0.0` (aniso branch, kperp2, tunits), `aky(ik) /= 0.0` (hybrid),
  `aky > epsilon(0.0)` (awgt, yavg), `aky < epsilon(0.0)` (phi_ext), `aky == epsilon(0.0)` (damp_zonal_only),
  `shat .eq. 0.0` (exb). Replicate each literally.
* **P15 esv on gamtot uses the local g_lo `(it,ik)` range** (no `force_all`, 6550-6554): under MPI, each rank's
  replicated `gamtot` is single-valued only for its own `(it,ik)` range → rank-dependent. Produce reference data
  with **1 MPI rank** when `esv = T` (or confirm the layout makes every rank own all `(it,ik)`).
* **P16 `awgt` is lazily computed on the first `calculate_flux_surface_average` call** and never refreshed;
  compute it with the filter state of that first call (all false in practice). `gamtot3` must come from the
  final `gamtot`.
* **P17 Adiabatic conditions differ**: `tite` added for no-electrons *or* no-ions, `fl_avg` only for
  no-electrons + fieldlineavg; `calculate_potentials_from_nonadiabatic_dfn` ignores fl_avg/poisfac/yavg.
  Field-equation operand order for `+fl_avg` differs between implicit (dist_fn) and local (fields_local) —
  field-solve spec must follow the active `field_option`; gf_local also writes `*beta*apfac/bmag(:)**2`
  (fields_gf_local.f90:3721) vs `*(beta*apfac/bmag**2)` elsewhere.
* **P18 Mixed complex/real ops** must be component-wise in JAX (§0.4), incl. `phinew - phi_avg`,
  `jacob*antot/gamtot`, `g*F`.
* **P19 Klein summation** in every trapezoid (fl_avg, awgt, field_line_average(_weight)); a plain or
  pairwise sum is off by O(eps·nth).
* **P20 Operator precedence.** `-code_dt*wunits(ik)*dlnf0drho` (1746) parses as `-((code_dt*wunits)*dlnf0drho)`
  (unary minus has the precedence of binary minus; negation is exact so the value is identical either way).
  `exp(- ( ... ))` negates the full product. `a*b/c*d` is `((a*b)/c)*d`.
* **P21 Parity**: non-linked leaves `gnew(0,1)` untouched; linked overwrites `ig = 0` too; parity runs before
  linked BCs inside `invert_rhs_1` and again after collisions.
* **P22 kwork_filter skip means "do not touch"** (not "zero"): `gnew(1,1)` stays equal to `g(1,1)` forever;
  `g_work` (source) is zero there because it is zeroed before the loop (6179).
* **P23 wfb/ttp**: hybrid mask treats the wfb (`il = ng2+1`) as passing; hyper/exb/parity have no wfb/ttp/forbid
  special cases (they act on every theta point).
* **P24 Tracer weights**: `integrate_species_original` zeroes the *caller's* weights array in place
  (`intent(in out)`, le_grids.f90:741) — harmless for GS2's local `wgt` arrays; just zero tracer weights.

---

## 7. Open questions

1. Can the 1e-12 target be met against the FMA-contracted `-O3` aarch64 reference over long runs, or should
   the validation reference be rebuilt with `-ffp-contract=off` (and XLA GPU contraction controlled)?
2. Replicate the theta0-branch field-zeroing bug (`phi(iright+1:,:iib,:)`) — recommended yes; confirm policy
   for replicating known GS2 bugs (also `damp_zonal_only`, frozen `D_res`/simple3D filter).
3. If JAX changes dt itself (adaptive CFL), must it keep `D_res`/simple3D filter at the export-time dt0
   (faithful) — recommended yes.
4. Response-matrix recomputation on dt change (fields spec): `kwork_filter(1,1)` already true for implicit
   fields; amplitude-dependent hyper filter inside the rebuild — confirm the fields spec replicates.
5. esv + MPI: confirm reference data is produced on one rank, or define which rank's `gamtot` is exported.
6. Hybrid electrons with collisions: does the collisions spec port `set_hzero_lorentz_collisions_matrix` and
   `zero_out_passing_hybrid_electrons`? Otherwise GATE `hybrid_electron + collisions`.
7. Non-Maxwellian `f0type` with collisions (operator uses `energy_maxw`): gate in the collisions spec?
8. Is the theta0 (ballooning) ExB branch actually exercised by target studies? It requires a 2π-periodic theta0
   grid and `nperiod ≥ 2` to be meaningful; SUPPORTed for faithfulness.
9. `wstar_units` with `aky(1) = 0` in a range grid produces `tunits = 0` (GS2 only warns) — gated; confirm.
10. The export program cannot read private variables (`gamtot3`, `awgt`, `akx4_max`, effective hyper
    parameters, `field_line_average_weight`, `ikx_indexed`, `dkx`, `parity_redist`): this spec gives exact
    recompute formulas; confirm no accessor patches are wanted.

---

## Appendix A — Export manifest for this subsystem

| netCDF variable | Fortran source | Fortran dims (order) | NumPy shape | dtype |
|---|---|---|---|---|
| `hyper_D_res` | `hyper::D_res` | `(ntheta0, naky)` | `(naky, ntheta0)` | f64 |
| `hyper_filter_s3d` | `hyper::hypervisc_filter` (only if simple3D) | `(theta, ntheta0, naky)` | `(naky, ntheta0, nth)` | f64 |
| `hyper_*` config | `get_hyper_config()` raw inputs + `code_dt` at init | scalars | — | f64/i32/bool |
| `tunits`, `wunits` | `gs2_time` | `(naky)` | `(naky,)` | f64 |
| `code_dt`, `code_dt_old`, `code_time` | `gs2_time` | scalars | — | f64 |
| `kx_shift` / `theta0_shift` | `dist_fn_arrays` (if allocated) | `(naky)` | `(naky,)` | f64 |
| `g_exb`, `g_exbfac`, `g_exb_start_time`, `g_exb_start_timestep` | dist_fn / config | scalars | — | f64/i32 |
| `kwork_filter` | `kt_grids` (after init_fields) | `(ntheta0, naky)` | `(naky, ntheta0)` | i8 |
| `ieqzip` | `run_parameters` (gate check) | `(ntheta0, naky)` | `(naky, ntheta0)` | i8 |
| `itleft`, `itright`, `l_links`, `r_links` | dist_fn | `(ntheta0, naky)` | `(naky, ntheta0)` | i32 |
| `def_parity`, `even`, `boundary_option_switch`, `esv`(effective) | dist_fn | scalars | — | i32 |
| `gamtot`, `gamtot1`, `gamtot2` | dist_fn | `(theta, ntheta0, naky)` | `(naky, ntheta0, nth)` | f64 |
| `tite`, `poisfac`, `adiabatic_option_switch` | run_parameters / dist_fn | scalars | — | f64/i32 |
| `theta`, `jacob`, `gds2`, `gds21`, `gds22`, `shat`, `ntheta`, `ntgrid`, `nperiod` | theta_grid | `(theta)` / scalars | `(nth,)` | f64/i32 |
| `akx`, `aky`, `ikx`, `theta0`, `kperp2`, `box`, `gridopt_switch` | kt_grids | as §2.1 | as §2.1 | f64/i32 |
| `spec_*` (z, mass, dens, temp, tprim, fprim, uprim, uprim2, vnewk, nu_h, bess_fac, stm, zstm, tz, zt, smz, type, f0type) | species | `(nspec)` | `(nspec,)` | f64/i32 |
| `nonmaxw_corr`, `dlnf0drho`, `f0_values` | species | `(negrid, nspec)` | `(nspec, negrid)` | f64 |
| `ng2`, `fphi`, `fapar`, `fbpar`, `has_phi/apar/bpar` | le_grids / run_parameters | scalars | — | |
| `remove_zonal_flows_switch` | fields | scalar | — | i32 |
| gate flags: `no_driver`, `gs2_has_lowflow`, `afilter`, `wstar_units`, `nonlin`, `eqzip_option_switch` | various | scalars | — | |

## Appendix B — Gate checks (export program and JAX loader)

Refuse to run if any of: `.not. no_driver` (antenna); `gs2_has_lowflow`; `any(ieqzip)`; `afilter > epsilon(0.)`;
`wstar_units .and. (nonlin .or. any(aky == 0.))`; ExB active with kx branch and `ntheta0 == 1`;
hyper active with a zero normaliser (`gridnorm` and `akperp4_max == 0`; aniso and `akx4_max == 0` or
`aky_max == 0`; simple3D with `shat == 0` or `kperp2_max(_zonal) == 0`); (pending Q6/Q7) hybrid or
non-Maxwellian species together with collisions if the collisions spec does not port their hooks.
