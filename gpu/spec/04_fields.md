# 04 — Field equations and the implicit field solve

Status: contract for implementers of the JAX (float64/complex128) port. Read-only analysis of GS2 8.2.1
(branch `gpu-port`, commit `4d8c94bcf`). Nothing in `src/` was modified. All `file:line` citations are relative to
`/Users/nico/Workspace/gs2/` (`src/` unless another directory is given).

Scope (this document):

* `src/dist_fn.fpp`: `getan`, `getan_from_dfn`, `getan_nogath`, `getan_nogath_from_dfn`, `getfieldeq`,
  `getfieldeq1`, `getfieldeq_nogath`, `getfieldeq1_nogath`, `init_fieldeq` (`gamtot`, `gamtot1`, `gamtot2`,
  `gamtot3`, `gridfac1`, adiabatic terms), `calculate_flux_surface_average` (`fl_avg`, `awgt`),
  `ensure_single_val_fields_pass(_r)`, `get_init_field` / `get_fields_direct_from_dfn` (+ denominators),
  `calculate_potentials_from_nonadiabatic_dfn`, the link / class set-up that defines the supercells
  (`count_links`, `setup_fields_implicit_class_arrays`, `get_leftmost_it`, `init_bc`).
* `src/le_grids.f90`: `integrate_species` (all variants) and the weights `w`, `wl` it uses.
* `src/fields_implicit.f90`: entirely (`get_field_vector`, `get_field_solution`, `getfield`, `advance_implicit`,
  `remove_zonal_flows`, `init_response_matrix`, `init_response_row`, `init_inverse_matrix`,
  `finish_fields_layouts`, `kt2ki`, dump/read response, `alloc_response_objects`).
* `src/gs2_layouts.fpp`: `f_lo` / `jf_lo` (field layouts) and their index functions.
* `src/fields_local.f90`: enough to state precisely where `field_option='local'` is and is not identical to
  `'implicit'`, and what an exact replication of the local path needs.
* `externals/utils/{integration,summation,matrix_inverse}.f*`: trapezoid rule, Klein compensated sum, block inverse.

Out of scope (owned by other specs, only referenced): the body of `timeadv` / `invert_rhs` / source terms
(spec 02 and the parallel-solve spec), collisions (`vspace_derivatives`), `hyper_diff` internals (spec 07),
the time loop / reset mechanics (spec 01), `exb_shear` remap (spec 07), diagnostics (spec 08).

---------------------------------------------------------------------------------------------------------------------

## 0. Conventions, notation, build facts

### 0.1 Index conventions (identical to the other specs)

* Distribution-function-like arrays are exported in canonical Fortran order
  `g(ig=-ntgrid:ntgrid, isgn=1:2, it=1:ntheta0, ik=1:naky, il=1:nlambda, ie=1:negrid, is=1:nspec)`
  (iglo expanded with `it_idx/ik_idx/il_idx/ie_idx/is_idx`, independent of `layout`).
  NumPy (C order): `g[is, ie, il, ik, it, isgn, ig+ntgrid]`, shape `(nspec, negrid, nlambda, naky, ntheta0, 2, nth)`.
  Arrays without `isgn` (`aj0(ig,iglo)`, `aj1`, `vperp2`) drop that axis: `aj0[is, ie, il, ik, it, ig+ntgrid]`.
* Field-like arrays `X(ig,it,ik)` (`phi`, `apar`, `bpar`, `*new`, `antot*`, `fieldeq*`, `gamtot*`, `gridfac1`,
  `kperp2`, `inv_kperp2`, `apar_ext`) → NumPy `X[ik, it, ig+ntgrid]`, shape `(naky, ntheta0, nth)`.
* `(it,ik)` arrays (`kwork_filter`, `awgt`, `fl_avg`, `itleft`, `itright`, `l_links`, `r_links`) → NumPy `[ik, it]`.
* Velocity / theta arrays keep the Fortran index order: `w(ie,is)` → NumPy `w[is, ie]`; `wl(ig,il)` → NumPy
  `wl[il, ig+ntgrid]`; `forbid(ig,il)` → `forbid[il, ig+ntgrid]`; `bmag(ig)`, `theta(ig)`, `jacob(ig)` → `[ig+ntgrid]`.
* Field-solve objects (new, defined in §3.11/§3.15): per class `c`, NumPy `aminv_c[m, r, k]` (member, output row,
  input column), both axes in the *blocked* order `B` defined below.

Notation used throughout:

| symbol | meaning |
|---|---|
| `nth` | `2*ntgrid+1` |
| `igp` | `ig+ntgrid` (0-based theta index) |
| `nf` (=`nfield`) | number of present fields, `has_phi+has_apar+has_bpar` (fields_implicit.f90:477-481); present fields are ordered phi, apar, bpar |
| `f` | field slot `1..nf` among present fields (Fortran) / `0..nf-1` (NumPy) |
| `nidx` | `nth*nf` (cell size, fields_implicit.f90:481) |
| class `c` | set of supercells with the same number of cells; `c = 1..i_class`, ascending `N_class` |
| `N_c`, `M_c` | `N_class(c)` cells per supercell, `M_class(c)` supercells (members) in class `c` |
| `ik_c[m]`, `it_c[m,n]` | `f_lo(c)%ik(m,n)` (constant in `n`), `f_lo(c)%it(m,n)`; cell `n=1` is leftmost |
| `nindex_c` | `nidx*N_c` (supercell size incl. duplicated boundary points) |
| **K** (interleaved, GS2 `am`/`f_lo` order) | `K(ig,f,n) = (f-1) + nf*(igp + nth*(n-1))`, 0-based, field fastest |
| **B** (blocked, `aminv`/`fl` order) | `B(ig,f,n) = igp + nth*((f-1) + nf*(n-1))`, 0-based, theta fastest |
| skip set `S_c` | `{(ig=ntgrid, n<N_c), all f}` — the duplicated right-boundary points (empty if `N_c=1`) |

### 0.2 Build facts that matter here (checked on this machine)

* `bin/gs2`: gfortran 15.2.0 (Homebrew), `-O3 -mcpu=apple-m1 -fdefault-real-8 -fdefault-double-8`, MPI build,
  **no OpenMP** (`Makefile:66`), **no `-DLAPACK`** (so `mat_inv::invert_serial` uses `invert_block`,
  externals/utils/matrix_inverse.fpp:194-198), arm64. `epsilon(0.0)` = `2.220446049250313e-16`.
* Disassembly of the shipped objects (`otool -tV objs/*.o`) shows, for the routines of this subsystem:
  1. **`integrate_species_original`** (inlined into `___le_grids_MOD_integrate_species_master`): the term
     `(weights(is)*w(ie,is))*wl(:,il)*(g(:,1,iglo)+g(:,2,iglo))` is computed as `fac=weights*w` (`fmul`),
     `fac2=wl*fac` (`fmul.2d`), `s=g1+g2` (`fadd.2d`), then a *full* complex multiply by `(fac2, 0.0)`
     (`fmul`, `fmls`/`fmla` against a zero register `movi v27,#0`), and finally `total = term + total` as a
     **separate `fadd`** (not fused). For finite operands this equals componentwise `fac2*s.re`, `fac2*s.im`
     followed by a rounded add. ⇒ this kernel is bitwise reproducible in JAX if XLA does not contract the add.
  2. **`getfield`** (`___fields_implicit_MOD_getfield.constprop.0`): the `sum(aminv*fl)` is a sequential
     accumulation starting from `+0.0` (`movi d18,#0`), one element after the other, real and imaginary
     accumulators separate (`fadd d17,d17,d4`…). The complex products use FMA: in the 2-wide vector body
     `p.re = RN(RN(x.re*y.re) − x.im*y.im)` (`fmls.2d`), `p.im = RN(RN(x.im*y.re) + x.re*y.im)` (`fmla.2d`); in the
     scalar remainder (odd element count) the *other* product is the fused one (`fmadd`/`fnmsub`). Then
     `u = u − S` (`fsub`).
  3. **`init_inverse_matrix`** (inlined into `init_fields_implicit`): complex division is GCC's "wide"
     Smith algorithm (`-fcx-fortran-rules`, see §6 P2) with FMA in the denominator; `am − lhscol*fac` uses
     ARMv8.3 `fcmla` (complex FMA) followed by `fsub`.
* **Round-off policy adopted by this spec (Tier A, mandatory):** same operations, same association
  (Fortran left-to-right), same loop and summation order, same algorithm for complex division (Smith, same
  branch rule), no re-association, no tree reductions in "exact mode". FMA placement is *not* reproduced
  (Tier B, optional; see open question Q1). Expected discrepancy: ≤ a few ulp per operation.

### 0.3 Executive summary — what the JAX code must do

Per time step (field_option='implicit', fields_implicit.f90:315-371):

1. After the first `timeadv` (predictor, with `phinew=phi`), compute the field-equation residuals
   `fieldeq, fieldeqa, fieldeqp` from `gnew` and the *old* fields (§3.1-§3.5).
2. For every supercell: `u = −A⁻¹·fl` with the exported/rebuilt inverse `aminv`, exact per-cell summation (§3.6).
3. `phinew = u + phi`, `aparnew = u_A + apar`, `bparnew = u_B + bpar` (§3.7); optional `remove_zonal_flows`.

At t=0 and after every dt change: fields are obtained **directly** from g (`get_init_field`, §3.8), and after
a dt change `aminv` must be rebuilt: unit field perturbations → `timeadv(istep=0)` → residual → GS2's
pivot-free Gauss–Jordan (§3.12-§3.13, GPU algorithm §5.4-§5.5).

---------------------------------------------------------------------------------------------------------------------

## 1. Call graph / control flow

### 1.1 Initialisation (init levels, src/gs2_init_level_list.inc:1-35)

```
level 27 dist_fn_level_1 -> init_bc (dist_fn.fpp:928-937, 2128-2157)
    linked:  init_connected_bc (3513-3610): jshift0, itleft/itright (2441-2515), count_links (2645-2672),
             setup_fields_implicit_class_arrays (2676-2755) -> i_class, N_class(:), M_class(:)
    other :  i_class=1, M_class=naky*ntheta0, N_class=1 (2144-2150)
    calculate_supercell_labels (2355-2397)
level 28 dist_fn_level_2  -> dist_fn::init_dist_fn_level_2 (dist_fn.fpp:939-953)
    init_bessel                                   (aj0, aj1)               [spec 02]
    init_fieldeq          dist_fn.fpp:6470-6572   gridfac1, gamtot, gamtot1, gamtot2, (gamtot3), adiabatic terms
level 32 fields_level_2 -> fields::init_fields (fields.f90:353-377)
    fields_pre_init (241-279): init_dist_fn (full), init_antenna, allocate_arrays (zeroes phi.. phinew..)
    fields_init_response (323-350)
      implicit: init_fields_implicit (fields_implicit.f90:65-99)
        init_response_matrix (459-647)
          nfield, nidx                                   (477-481)
          init_fields_layouts  (gs2_layouts.fpp:1650-1705)  f_lo(1:i_class) sizes / blocks
          init_jfields_layouts (gs2_layouts.fpp:1855-1908)  jf_lo
          finish_fields_layouts (fields_implicit.f90:1074-1133) f_lo(c)%ik(m,n), %it(m,n), jf_lo%ij
          [read_response -> read_response_from_file_imp (1408-1713)]
          do c = i_class, 1, -1                          (541)
             zero am, g, phi, phinew, apar, aparnew, bpar, bparnew   (553-557)
             do n = 1,N_c ; do ig = -ntgrid,ntgrid ; per present field:
                set unit perturbation (566-633)
                init_response_row (650-719): timeadv(...,istep=0) (669); getfieldeq (670); store columns
                zero the perturbed field
             init_inverse_matrix(am, c) (722-1071): Gauss-Jordan, skip fill, reorder into aminv(c)
          [dump_response -> dump_response_to_file_imp (1172-1404)]
        init_antenna
      local: init_fields_local (fields_local.f90:3937-3977)  (see §7)
    nl_finish_init
    if (gridopt_switch == gridopt_box) kwork_filter(1,1) = .true.     (fields.f90:376)   <-- AFTER the response
level 34 set_initial_values -> gs2_init.fpp:292-375 -> fields::set_init_fields (fields.f90:401-440)
    implicit/local: init_allfields_implicit (fields_implicit.f90:102-136)
       new_field_init (init_g default .true., init_g.f90:262):
          get_init_field(phinew,aparnew,bparnew) (dist_fn.fpp:7792-7799) -> get_fields_direct_from_dfn (7590-7790)
          phi=phinew, apar=aparnew, bpar=bparnew ; g = gnew
       else: getfield(phinew,aparnew,bparnew) ; copies           (fields_implicit.f90:131-132)
```

### 1.2 One time step (`fields::advance`, fields.f90:586-605 → `advance_implicit`, fields_implicit.f90:315-371)

```
advance_implicit(istep):
  if (.not.no_driver) antenna_amplitudes(apar_ext)                            (333)  [GATE]
  if (allocated(kx_shift) .or. allocated(theta0_shift)) exb_shear(gnew, phinew, aparnew, bparnew, istep) (335) [spec 07]
  g = gnew ; phi = phinew ; apar = aparnew ; bpar = bparnew                     (337-338)  (all three, regardless of has_*)
  timeadv(phi,apar,bpar, phinew,aparnew,bparnew, istep)                         (341)  -> gnew = g*  (predictor, phinew==phi)
  if (reset) return                                                             (343)  (CFL reset inside timeadv, spec 01)
  if (.not.no_driver) aparnew = aparnew + apar_ext                              (345)  [GATE]
  getfield(phinew, aparnew, bparnew)                                            (352)
     get_field_vector(fl, phinew, aparnew, bparnew)                             (139-171)
        getfieldeq(phinew,aparnew,bparnew, fieldeq,fieldeqa,fieldeqp)          (dist_fn.fpp:6667-6680)
           getan(antot,antota,antotp) -> getan_from_dfn(gnew,...)              (6372-6380, 6273-6360)
              integrate_species -> integrate_species_master -> integrate_species_original (le_grids.f90:997-1139, 724-786)
              [ensure_single_val_fields_pass if esv]                           (6201-6233)
           getfieldeq1 (6575-6665)
              [calculate_flux_surface_average if no electrons & fieldlineavg]  (7994-8034)
        fl(:,it,ik) = [fieldeq ; fieldeqa ; fieldeqp] (present fields only)
     u(jflo) = - sum_n sum( aminv(c)%dcell(dc)%supercell(cell n) * fl(:,it_n,ik) )  (261-289)
     sum_allreduce(u) / allgatherv                                              (295-301)
     get_field_solution(u) -> phinew/aparnew/bparnew = u (present fields)        (174-219)
  phinew = phinew + phi ; aparnew = aparnew + apar ; bparnew = bparnew + bpar    (354-356)
  if (remove_zonal_flows_switch) remove_zonal_flows                              (362, 374-387)
  timeadv(phi,apar,bpar, phinew,aparnew,bparnew, istep, diagnostics)             (365)  corrector -> final gnew
  collisions_advance(phi,bpar, phinew,aparnew,bparnew, diagnostics)              (369)  (no-op unless split_collisions)
     -> g_adjust ; calculate_potentials_from_nonadiabatic_dfn (dist_fn.fpp:7322-7455) ; g_adjust   [GATE]
```

### 1.3 dt change (spec 01 §reset; gs2_reinit.f90:104-225)

`gs2_main.fpp:739-751` / `:807-811` → `reset_time_step` → `zero_array(gnew)` (gs2_reinit.f90:163) →
`init(current_init, override_timestep)` (166) finishes levels ≥ 29 → new `code_dt` → `init(..., full)` (216)
re-runs **collisions (30), dist_fn_level_3 (31), fields_level_2 (32)** — i.e. `fields_implicit::reset_init`
(fields_implicit.f90:428-456) then `init_response_matrix` with the new dt — and set_initial_values (34):
g restored from memory/restart, then (`force_maxwell_reinit`, default `.true.`, fields.f90:118; gs2_init.fpp:347)
`set_init_fields` → **fields recomputed with `get_init_field` from g** (not via `aminv`).
Level 28 (`init_fieldeq`) and level 13 (`kt_grids`, holder of `kwork_filter`) are **not** re-run: `gamtot*`
persist, `kwork_filter(1,1)` stays `.true.` during the rebuild, `awgt` persists (deallocated only at level 27,
dist_fn.fpp:1100).

### 1.4 field_option='local' (fields_local.f90)

```
advance_local (4244-4345): as advance_implicit but: exb_shear if |g_exb*g_exbfac|>eps (4280); copies and
  "phinew=phinew+phi" only for present fields (4284-4291, 4316-4320); getfield_local (4193-4232):
     fm_get_field_update (2298-2343): fm_getfieldeq_nogath (3084-3106) [has_to_gather=.false. for decomp 3,
        3788] -> getan_nogath (7294-7302) + fm_getfieldeq1_nogath (3374-3449)
     per supercell/cell/field: matmul(vect, rb%data) (370-384), c_get_field_update (512-546),
        sc_reduce_tmpsum (723-777); fm_gather_fields (2593-2804) incl. duplicate fill
init: init_fields_matrixlocal (3979-4174): fm_populate (1973-2078) / fm_init_next_field_points (2081-2273),
      fm_prepare -> sc_invert_local (1238-1287): invert_serial = invert_block (matrix_inverse.fpp:337-391)
```

### 1.5 Other callers of the field primitives

* `get_fields_direct_from_dfn` is also used by the split nonlinear scheme (dist_fn.fpp:4049-4057, 4097-4105)
  and `calculate_potentials_from_nonadiabatic_dfn` by `collisions_advance` (4130-4170) and split NL with
  `advance_nonadiabatic_dfn`. Both are specified in §3.8-§3.9; the options that reach them are GATEd here.
* Diagnostics call `getan`, `get_field_inconsistency` etc. (spec 08).

---------------------------------------------------------------------------------------------------------------------

## 2. Inventory of everything the field solve reads

Decision legend: **EXPORT** = dump from Fortran after full initialisation; **RECOMPUTE** = port the formula;
**BOTH** = export for t=0 / validation and port the formula for runtime rebuilds.

| Fortran name | module | declared shape / ranges | computed at (file:line) | dt / runtime dependence | decision | notes |
|---|---|---|---|---|---|---|
| `aminv(c)%dcell(dc)%supercell(:)` | fields_implicit (private) | per class: one complex vector of length `nindex_c` per output row `jflo` | init_inverse_matrix fields_implicit.f90:920-1066 | **dt (code_dt), vnmult at build time, hyper, collisions, kwork_filter** | **BOTH** | export as `aminv_c[m, B_out, B_in]` (§3.15, §8); rebuilt on dt change (§5.4-5.5). Not rebuilt when vnmult changes (P10) |
| `am` (response, pre-inversion) | local to init_response_matrix | `(nindex_c, M_c*nindex_c)` | init_response_row 650-719 | dt | RECOMPUTE (optional EXPORT for validation, needs a hook, Q10) | `am(k, idx(j,m))` = response k to perturbation j |
| `i_class`, `N_class(:)`, `M_class(:)` | dist_fn | scalar, `(i_class)` | dist_fn.fpp:2676-2755 / 2144-2150 | static | EXPORT | `N_class` strictly ascending |
| `f_lo(c)%ik(m,n)`, `f_lo(c)%it(m,n)` | gs2_layouts | `(M_c, N_c)` per class | fields_implicit.f90:1081-1131 | static | EXPORT (+RECOMPUTE check) | ordering rule §3.11 |
| `itleft`, `itright` | dist_fn | `(ntheta0, naky)`, `-1` = no link | dist_fn.fpp:2441-2515 | static | EXPORT | ky=0 in box: `itright(it,1)=it` (self) |
| `l_links`, `r_links` | dist_fn | `(ntheta0, naky)` | 2608-2672 | static | EXPORT | `N = 1+l+r`, cell index `n = 1+l` |
| `gamtot` | dist_fn | real `(-ntgrid:ntgrid, ntheta0, naky)` | init_fieldeq 6517-6568 | static (species/geometry) | EXPORT | includes `poisfac*kperp2`, `+tite` adiabatic terms, esv |
| `gamtot1`, `gamtot2` | dist_fn | same | 6527-6547 | static | EXPORT | |
| `gamtot3` | dist_fn | same, allocated only if `adiabatic_option_switch==fieldlineavg` | 6511-6513, 6564-6565 | static | EXPORT (if allocated) | |
| `gridfac1` | dist_fn | real `(-ntgrid:ntgrid, ntheta0, naky)` | 6486-6506 | static | EXPORT | ≡1 for default `gridfac=1` |
| `awgt` | dist_fn | real `(ntheta0, naky)` | first call of calculate_flux_surface_average 8006-8020 | set once, then frozen (P8) | BOTH | needed only if no electrons & fieldlineavg |
| `fl_avg` | dist_fn | complex `(ntheta0, naky)` | 8022-8033 | per call | RECOMPUTE | |
| `antot`, `antota`, `antotp` | dist_fn_arrays | complex `(-ntgrid:ntgrid, ntheta0, naky)` persistent | getan_from_dfn 6273-6360 | per call | RECOMPUTE | absent-field arrays keep old values (unused) |
| `fieldeq`, `fieldeqa`, `fieldeqp` | dist_fn_arrays | same, persistent | getfieldeq1 6575-6665 | per call | RECOMPUTE | |
| `g_work` | dist_fn_arrays | g-like | scratch | per call | RECOMPUTE | |
| `gnew` | dist_fn_arrays | g-like | timeadv | state | input from timeadv (spec 02/03) | forbidden-region values must be finite (P18) |
| `phi, apar, bpar, phinew, aparnew, bparnew` | fields_arrays | field-like | state | state | BOTH | EXPORT t=0; RECOMPUTE via get_init_field on reset |
| `aj0`, `aj1` | dist_fn_arrays | real `(-ntgrid:ntgrid, g_lo)` | init_bessel (spec 02) | static | EXPORT | canonical `[is,ie,il,ik,it,ig]` |
| `vpa` | dist_fn_arrays | real `(-ntgrid:ntgrid, 2, g_lo)` | init_vpa_vpac (spec 02) | static | EXPORT | g-like |
| `vperp2` | dist_fn_arrays | real `(-ntgrid:ntgrid, g_lo)` | init_vperp2 (spec 02) | static | EXPORT | |
| `w` | le_grids | real `w(negrid, nspec)` | set_vgrid le_grids.f90:285-296 | static | EXPORT | NumPy `w[is, ie]` |
| `wl` | le_grids | real `wl(-ntgrid:ntgrid, nlambda)` | lgridset 2120-2156 (+2160-2641) | static | EXPORT | `wl = 0` where forbidden (2129, 2630-2638) |
| `forbid` | le_grids | logical `(-ntgrid:ntgrid, nlambda)` | 2715-2738 | static | EXPORT (validation) | assert `wl[forbid]==0` |
| `kperp2`, `inv_kperp2` | kt_grids | real field-like | kt_grids.f90:312-347 | static | EXPORT | `inv_kperp2=0` where `|kperp2|≤eps` |
| `kwork_filter` | kt_grids | logical `(ntheta0, naky)` | kt_grids.f90:188-189, fields.f90:376 | **changes after init_fields** (P6) | EXPORT (state after init) | `(1,1)=.true.` in box after response build |
| `aky`, `akx` | kt_grids | `(naky)`, `(ntheta0)` | kt_grids | static | EXPORT | zonal tests §3.13 |
| `bmag`, `theta`, `jacob` | theta_grid | real `(-ntgrid:ntgrid)` | theta_grid.f90:516 (jacob) | static | EXPORT | |
| `field_line_average_weight` | theta_grid | real | theta_grid.f90:520 | static | EXPORT | used by remove_zonal_flows |
| `spec%z, dens, temp, mass, stm, type` | species | `(nspec)` | species.f90:476 (stm) | static | EXPORT | `stm = sqrt(temp/mass)` |
| `has_electron_species(spec)`, `has_ion_species(spec)` | species | logical | species.f90:914-936 | static | EXPORT | electron-like = electron + hybrid (species.f90:42,944) |
| `beta`, `tite` | run_parameters | real | | static | EXPORT | |
| `has_phi, has_apar, has_bpar` | run_parameters | logical | run_parameters.f90:351-353 (`f* > eps`) | static | EXPORT | |
| `apfac`, `poisfac`, `gridfac`, `afilter`, `esv` | dist_fn | real/logical | read_parameters dist_fn.fpp:1240-1330 | static | EXPORT | `esv = esv .and. linked` (1327) |
| `adiabatic_option_switch`, `boundary_option_switch` | dist_fn | integer | 1294-1316 | static | EXPORT | |
| `inv_phi_denominator_g`, `inv_bpar_denominator_g` | dist_fn | real field-like | 7457-7585 (lazily, first get_fields_direct call) | static | BOTH | |
| `nonmaxw_corr(ie,is)` | species | real | | static | EXPORT (only if recomputing gamtot) | |
| `layout` | gs2_layouts | `character(5)` | gs2_layouts.fpp:204,540 | static | EXPORT | fixes velocity summation order (§3.1) |
| `apar_ext` | fields_arrays | field-like | antenna | time | GATE | antenna |
| `code_dt` | gs2_time | real | | dt | input to rebuild | |
| `jf_lo%ij/mj/dj`, `recvcnts`, `displs` | gs2_layouts / fields_implicit | MPI bookkeeping | 1110-1118, 949-958 | — | not needed | |
| `supercell_labels`, `n_cells`, `n_supercells` | dist_fn | `(ntheta0,naky)`, `(naky)` | 2355-2397 | static | EXPORT (local path only) | |

---------------------------------------------------------------------------------------------------------------------

## 3. Exact algorithms (Tier A order of operations)

Fortran evaluates `a*b*c` as `(a*b)*c`, `a+b-c+d` as `((a+b)-c)+d`, `x**2` as `x*x` (gfortran expands small
integer powers). `real*complex` and `complex/real` are numerically componentwise for finite operands (GCC
emits a full complex op with a zero imaginary part, see 0.2.1; results differ at most in the sign of zero).
In JAX write them explicitly componentwise (`lax.complex(r*z.real, r*z.imag)`) rather than relying on dtype
promotion, and write complex×complex and complex÷complex in real arithmetic (P1, P2).

### 3.1 Velocity-space integral `integrate_species` (le_grids.f90:724-786, dispatch 997-1139)

Dispatch: `integrate_species_master(g, weights, total[, nogath, gf_lo])` → if `gf_lo` → `integrate_species_gf`
(GATE); if `.not.intspec_sub` → **`integrate_species_original`** (always the case with 1 MPI rank:
gs2_layouts.fpp:1351-1378 force `intspec_sub=.false.` when x,y are local); else sub-communicator paths
(multi-rank only, GATE for strict validation).

```
integrate_species_original(g, weights, total):            # le_grids.f90:724-786
  total = 0                                                # zero_array (+0.0)
  where (spec%type == tracer_species) weights = 0          # 741
  for iglo = llim_proc .. ulim_proc:                       # iglo order of the *layout*
      (it,ik,il,ie,is) = idx(iglo)
      if any(kwork_filter) and kwork_filter(it,ik): cycle   # 751-753 (branch only differs by the skip)
      total(:,it,ik) = total(:,it,ik) + (weights(is)*w(ie,is))*wl(:,il)*(g(:,1,iglo)+g(:,2,iglo))   # 759-760
  sum_allreduce(total)                                     # 785 (identity with 1 rank)
```

**Summation order.** For a fixed `(ig,it,ik)` the terms are added in iglo order restricted to that `(it,ik)`,
i.e. in the relative order of `l`, `e`, `s` in the layout string (gs2_layouts.fpp:917-930):

| layout | order of the velocity loop (outer → inner) |
|---|---|
| `lxyes` (default), `lyxes`, `lexys`, `xyles`, `yxles` | `is` outer, `ie`, `il` innermost |
| `yxels` | `is` outer, `il`, `ie` innermost |

Exact per-element arithmetic (verified by disassembly, 0.2.1):

```
fac  = weights(is) * w(ie,is)            # real
fac2 = fac * wl(ig,il)                   # real (commutative product; any operand order is identical)
s    = g(ig,1,iglo) + g(ig,2,iglo)       # complex add (componentwise)
term = (fac2*s.re , fac2*s.im)           # two rounded products
total(ig,it,ik) = total + term           # rounded add (NOT fused in the shipped binary)
```

JAX (exact mode):

```python
def integrate_species_exact(gw, wgt, w, wl, kwork_filter, spec_is_tracer, layout):
    # gw: complex (ns, ne, nl, nky, nkx, 2, nth); wgt: (ns,); w: (ns, ne) == Fortran w(ie,is); wl: (nl, nth)
    wgt = jnp.where(spec_is_tracer, 0.0, wgt)
    fac  = wgt[:, None] * w                                # (ns, ne)
    fac2 = fac[:, :, None, None] * wl[None, None, :, :]    # (ns, ne, nl, nth)   == fac*wl(ig,il)
    s    = gw[..., 0, :] + gw[..., 1, :]                   # (ns, ne, nl, nky, nkx, nth)
    term = fac2[:, :, :, None, None, :] * s                # real*complex componentwise
    if layout == 'yxels':
        term = term.transpose(0, 2, 1, 3, 4, 5)            # (ns, nl, ne, ...) -> ie innermost
    term = term.reshape((-1,) + term.shape[3:])            # flatten velocity: C order => last velocity axis fastest
    def body(tot, t): return tot + t, None                 # sequential, one rounded add per term
    tot, _ = jax.lax.scan(body, jnp.zeros(term.shape[1:], term.dtype), term)
    return jnp.where(kwork_filter[:, :, None], 0.0 + 0.0j, tot)   # skipped modes stay exactly +0
```

(The products `fac`, `fac2`, `term` may be computed vectorised — they involve no summation; only the
accumulation must be sequential. A Pallas kernel with one thread per output and an in-register sequential loop
is the fast exact implementation, §5.2.)

`total` for filtered modes is exactly `+0`. `wl(ig,il)=0` at forbidden points (lgridset sets `wl=0` first,
le_grids.f90:2129, and only allowed points receive weights), so forbidden-region `g` contributes `0*s`
— exact zero **provided `s` is finite** (P18).

`integrate_species_sub` (le_grids.f90:789-879, used by `getan_nogath` when `intspec_sub`): identical per-element
arithmetic into `total_small`, reduction over the xy-block sub-communicator — only multi-rank. `integrate_species_gf`
(882-936), `integrate_species_gf_nogather` (944-993): loop order `igf; il; ie; is` (species innermost) — GATE
(`gf_lo_integrate`).

### 3.2 `getan_from_dfn` (dist_fn.fpp:6273-6360) — `getan(…)` = `getan_from_dfn(gnew, …)` (6372-6380)

```
if has_phi:                                                     # 6290-6318
   g_work(ig,isgn,iglo) = aj0(ig,iglo)*g_in(ig,isgn,iglo)       # real*complex
   wgt = spec%z*spec%dens                                       # 6313
   integrate_species(g_work, wgt, antot)                        # 6314
   if (afilter > epsilon(0.0)) antot = antot*exp(-afilter**4*kperp2**2/4.)   # 6316  [GATE afilter]
   if (esv) ensure_single_val_fields_pass(antot, force_all=.true.)           # 6317
if has_apar:                                                    # 6320-6338
   g_work = aj0*vpa*g_in            == (aj0(ig,iglo)*vpa(ig,isgn,iglo))*g_in(ig,isgn,iglo)
   wgt = 2.0*beta*spec%z*spec%dens*sqrt(spec%temp/spec%mass)    == (((2.0*beta)*z)*dens)*sqrt(temp/mass)   # 6335
   integrate_species(g_work, wgt, antota) ; esv as above        # 6336-6337
if has_bpar:                                                    # 6340-6358
   g_work = aj1*vperp2*g_in         == (aj1(ig,iglo)*vperp2(ig,iglo))*g_in
   wgt = spec%temp*spec%dens                                    # 6355
   integrate_species(g_work, wgt, antotp) ; esv                 # 6356-6357
```

Absent fields: the corresponding `antot*` arrays are *not* touched (keep their previous value, zero from
allocation, dist_fn.fpp:3994-3997). They are never read in a way that matters (see 3.5).

`getan_nogath_from_dfn` (6969-7290) with `gf_lo_integrate=.false.` performs the same arithmetic: the `g_work`
loops skip `kwork_filter` modes (7059-7080, 7136-7160, 7217-7240; irrelevant, integrate skips them too), `integrate_species(…, nogath=.true.)`
(→ `integrate_species_original` when `.not.intspec_sub`), apar weight uses `spec%stm` (7133), bitwise equal to
`sqrt(spec%temp/spec%mass)` (species.f90:476), `esv` without `force_all` (local `g_lo` ranges = full ranges
with one rank), afilter identical (7126). With `gf_lo_integrate=.true.` (7046-7047, 7094-7120, 7171-7200, 7252-7282) the
integral is fused with a different grouping `wgt(is)*w(ie,is)*wl(:,il)*((aj0*gf1)+(aj0*gf2))` and loop order
`il, ie, is` — GATE.

### 3.3 `ensure_single_val_fields_pass(arr, force_all)` (dist_fn.fpp:6201-6233; `_r` 6235-6271)

```
if (boundary_option_switch /= boundary_option_linked) return
ranges = force_all ? full : g_lo%{ik,it}_{min,max}       # local ranges (identical with 1 rank)
do ik ; do it
   link_it = itright(it,ik); if (link_it < 0) cycle
   arr(-ntgrid, link_it, ik) = arr(ntgrid, it, ik)
```

Destinations are always `ig=-ntgrid`, sources always `ig=+ntgrid`, and `itright(:,ik)` is injective ⇒
order-independent ⇒ one vectorised scatter. For ky=0 in box `itright(it,1)=it` ⇒ the zonal modes get
`arr(-ntgrid)=arr(+ntgrid)` (periodic copy).

### 3.4 Flux-surface-average term `fl_avg` (dist_fn.fpp:7994-8034)

Called (only) when `.not. has_electron_species(spec) .and. adiabatic_option_switch == adiabatic_option_fieldlineavg`
(getfieldeq1 6605-6609, getfieldeq1_nogath 6734-6738, fm_getfieldeq1_nogath fields_local.f90:3390-3394,
get_fields_direct_from_dfn 7637-7639).

```
fl_avg = 0
if (.not. allocated(awgt)):                               # first call only (8006)
   awgt = 0
   do ik: if (aky(ik) > epsilon(0.0)) cycle
      do it: if (kwork_filter(it,ik)) cycle               # filter state AT FIRST CALL
         awgt(it,ik) = 1.0/trapezoidal_integration(theta, jacob*gamtot3(:,it,ik))     # 8017
do ik ; do it
   if (kwork_filter(it,ik)) cycle
   fl_avg(it,ik) = tite*trapezoidal_integration(theta, jacob*antot(:,it,ik)/gamtot(:,it,ik))*awgt(it,ik)   # 8031
```

Order: integrand `y(ig) = (jacob(ig)*antot(ig,it,ik))/gamtot(ig,it,ik)` (real×complex, then ÷real, componentwise);
result `(tite*T)*awgt(it,ik)`. For `aky>eps` modes `awgt=0` ⇒ `fl_avg = ±0` (the JAX code may skip the
integral there and use `0`).

Trapezoid (externals/utils/integration.f90:55-71) with `n = nth`:

```
v(i) = (x(i+1) - x(i)) * (y(i) + y(i+1)),  i = 1..n-1      # real*complex componentwise for the complex version
integral = local_sum(v) * 0.5
local_sum (complex) = cmplx(klein_sum(real(v)), klein_sum(aimag(v)))     # integration.f90:42-52 (use_compensated_sum=.true., :24)
```

Klein (second-order compensated) sum, externals/utils/summation.f90:71-103, sequential over `i`:

```
s = 0; c1 = 0; c2 = 0
for a in array:
    t = s + a
    lc1 = (s - t) + a      if |s| >= |a|   else   (a - t) + s
    s = t
    t = c1 + lc1
    lc2 = (c1 - t) + lc1   if |c1| >= |lc1| else  (lc1 - t) + c1
    c1 = t
    c2 = c2 + lc2
return (s + c1) + c2
```

(`cmplx` without KIND is complex(8) under `-fdefault-real-8`.) JAX: `lax.scan` over `i` (length `2*ntgrid`)
vectorised over all `(it,ik)`; no multiplications inside ⇒ FMA-proof.

### 3.5 Field-equation residuals `getfieldeq1` (dist_fn.fpp:6575-6665)

`getfieldeq(phi,apar,bpar,…)` = `getan` + `getfieldeq1(phi,apar,bpar, antot,antota,antotp, …)` (6667-6680).
All `(it,ik)` (no `local_only` in the implicit path), no `kwork_filter` skip.

```
fl_avg = 0 ; if (no electrons & fieldlineavg) fl_avg = calculate_flux_surface_average(antot)    # 6602-6609
if has_phi:   fieldeq  = ((antot + bpar*gamtot1) - (gamtot*gridfac1)*phi) + fl_avg(it,ik)      # 6620-6623
if has_apar:  fieldeqa = antota - (kperp2*gridfac1)*apar                                       # 6638-6639
if has_bpar:  fieldeqp = ((antotp + bpar*gamtot2) + (0.5*phi)*gamtot1) * ((beta*apfac)/(bmag*bmag)) + bpar*gridfac1   # 6656-6659
```

Notes: the phi residual always contains `bpar*gamtot1` (bpar ≡ 0 when absent ⇒ `antot + 0`, exact); the bpar
residual always contains `0.5*phi*gamtot1` (phi ≡ 0 when absent). `bmag(ig)` broadcasts over `(it,ik)`.
Sign convention: these are residuals `F(fields, g)`; the solve finds `δ` with `F(fields+δ, g+Rδ) = 0`.

`getfieldeq1_nogath` (6703-6830, used by diagnostics/`check_getan`): same formulas, skips `kwork_filter` modes,
and when the companion field is absent drops the term (`fieldeq = (antot − gamtot*gridfac1*phi) + fl_avg`,
6768-6770; `fieldeqp = (antotp + bpar*gamtot2)*(…) + bpar*gridfac1`, 6824-6826) — bitwise identical for finite
values. `fm_getfieldeq1_nogath` (fields_local, §7) has a **different association** for the phi residual when
`has_bpar`.

### 3.6 Applying the inverse: `get_field_vector`, `getfield`, `get_field_solution`

`get_field_vector` (fields_implicit.f90:139-171): `fl(nidx, ntheta0, naky)` with
`fl(1+(f-1)*nth : f*nth, it, ik) = fieldeq_f(:, it, ik)` (field-blocked; `f` over present fields).

`getfield` (222-312):

```
u_small(:) = 0                                              # 258
do jflo (output index, jf_lo order):                        # 261
   c  = jf_lo%ij(jflo) ; m = jf_lo%mj(jflo) ; ik = f_lo(c)%ik(m,1) ; dc = jf_lo%dj(c,jflo)
   do n = 1, N_c                                            # 276
      it = f_lo(c)%it(m,n)
      u_small(jflo) = u_small(jflo) - sum( aminv(c)%dcell(dc)%supercell(1+nidx*(n-1) : nidx*n) * fl(:, it, ik) )   # 286
sum_allreduce(u_small)   [or allgatherv if field_subgath]   # 295-301  (exact: each entry is nonzero on exactly one rank)
get_field_solution(u)                                       # 174-219
```

`sum(...)` is a sequential left-to-right accumulation from `+0.0` over the `nidx` products (verified, 0.2.2).
Hence for each output row `r` of member `m`:

```
u_r = ((( 0 − S_1 ) − S_2 ) − … ) − S_N ,     S_n = (((0 + p_1) + p_2) + …) + p_nidx ,
p_k = aminv_row_r[cell n, k] * fl_k(it_n, ik)     (k in blocked order: f outer, ig inner)
```

`get_field_solution` (174-219): with `jflo = ij_idx(jf_lo, ig, f, ik, it) = igp + nth*((f-1) + nf*((it-1) + ntheta0*(ik-1)))`
(gs2_layouts.fpp:1967-1972), `phinew(:,it,ik) = u(ll:ll+2*ntgrid)` for `f`=phi slot, likewise apar/bpar; absent
fields are not written (they keep the value copied at fields_implicit.f90:338, i.e. the old field, zero).

### 3.7 Field update and `remove_zonal_flows` (fields_implicit.f90:352-387)

```
phinew  = phinew  + phi        # u + phi^n   (complex add)   354
aparnew = aparnew + apar                                      355
bparnew = bparnew + bpar                                      356
```

All three unconditionally: for an absent field `X` this is `X^n + X^n` (stays 0 because absent fields are 0).

`remove_zonal_flows` (374-387; `fieldline_average_phi`, 398-425):

```
phi_avg = 0 ; for it = 1..ntheta0:  phi_avg(:, it, 1) = field_line_average(phinew(:, it, 1))     # ik = 1 ASSUMED zonal
phinew = phinew - phi_avg                                       # other ky: subtract exact 0
field_line_average_complex(q) = trapezoidal_integration(theta, q*jacob) * field_line_average_weight   # theta_grid.f90:665-670
field_line_average_weight     = 1.0/trapezoidal_integration(theta, jacob)                              # theta_grid.f90:520
```

### 3.8 Direct field solve `get_fields_direct_from_dfn` (dist_fn.fpp:7590-7790), `get_init_field` (7792-7799)

Used at t=0, after every reset (force_maxwell_reinit), and by split NL. Inverts the local (no parallel
coupling) field equations:

```
if (.not. allocated(inv_phi_denominator_g)) setup_get_fields_direct_from_dfn_denominators   # 7620-7623
getan_from_dfn(g_in, antot, antota, antotp)          # integrals of g_in (kwork_filter applies)
fl_avg = 0 ; if (no electrons & fieldlineavg) fl_avg = calculate_flux_surface_average(antot)  # 7632-7639
has_phi:
   has_bpar: num = (beta*gamtot2 + bmag**2)*(antot + fl_avg(it,ik)) - (beta*gamtot1)*antotp   # 7655
   else    : num = antot + fl_avg(it,ik)                                                       # 7670
   phi = num * inv_phi_denominator_g                                                           # 7684
else: phi = 0
has_apar: apar = antota * inv_kperp2      else 0                                               # 7711
has_bpar:
   has_phi: num = (-(beta*gamtot)*antotp) - ((beta/2.0)*gamtot1)*antot                         # 7746
   else   : num = -beta*antotp                                                                  # 7758
   bpar = num * inv_bpar_denominator_g                                                          # 7772
else: bpar = 0
```

Denominators (7457-7585), computed once, dt-independent:

```
phi_denom = has_bpar ? (beta*gamtot2 + bmag**2)*gamtot + ((beta/2.0)*gamtot1)*gamtot1   # 7481-7483
                     : gamtot                                                            # 7491
inv_phi_denominator_g = |phi_denom| > epsilon(0.0) ? 1.0/phi_denom : 0.0                # 7503-7508
bpar_denom = has_phi ? gamtot*(beta*gamtot2 + bmag**2) + ((beta/2.0)*gamtot1)*gamtot1   # 7537-7539
                     : beta*gamtot2 + bmag**2                                            # 7552
inv_bpar_denominator_g = |bpar_denom| > epsilon(0.0) ? 1.0/bpar_denom : 0.0             # 7567-7572
```

(`bmag**2 = bmag*bmag`; real multiplication is commutative so `gamtot*(x)` = `(x)*gamtot` bitwise.)
Note: `gridfac1` is not used here (consistent only for `gridfac=1`, P27). `init_allfields_implicit` then copies
`phinew→phi`, `aparnew→apar`, `bparnew→bpar` and `gnew→g` (fields_implicit.f90:126-129).

### 3.9 `calculate_potentials_from_nonadiabatic_dfn` (dist_fn.fpp:7322-7455)

Fields from `h` (used by `collisions_advance` with `split_collisions`, and split NL with `advance_nonadiabatic_dfn`):

```
getan_from_dfn(f_in, ...)
phi_weight = 1.0/(sum(spec%dens*spec%z*spec%z/spec%temp) + tite)   if (no electrons .or. no ions)   # 7363
           = 1.0/ sum(spec%dens*spec%z*spec%z/spec%temp)            otherwise                         # 7365
   (sum: sequential over is=1..nspec of ((dens*z)*z)/temp)
phi  = antot*phi_weight                      (has_phi, else 0)
apar = antota*inv_kperp2                     (has_apar, else 0)
bpar = antotp*bpar_weight, bpar_weight(ig) = -beta/bmag(ig)**2   (has_bpar, else 0)    # 7422
```

No `fl_avg` here. Status: specified; the options reaching it (`split_collisions`, `split_nonlinear`) are GATEd
here pending their own specs.

### 3.10 `init_fieldeq` (dist_fn.fpp:6470-6572) — for validation / optional recompute (EXPORT is the decision)

```
gridfac1 = 1                                                         # 6487
linked : do it; do ik; if (aky(ik)==0.0) cycle
             if (itleft(it,ik)<0)  gridfac1(-ntgrid,it,ik) = gridfac
             if (itright(it,ik)<0) gridfac1( ntgrid,it,ik) = gridfac          # 6491-6499
self_periodic: nothing ; default (zero, alternate-zero): both ends = gridfac for aky/=0   # 6500-6506
g_work(:,isgn,iglo) = (1.0 - aj0**2)*nonmaxw_corr(ie,is)           ; wgt = ((z*z)*dens)/temp
   gamtot  = real(integrate_species(...)) + kperp2*poisfac          # 6516-6524
g_work = ((aj0*aj1)*2.0)*vperp2*nonmaxw_corr                       ; wgt = z*dens   -> gamtot1   # 6526-6536
g_work = ((aj1**2)*2.0)*(vperp2**2)*nonmaxw_corr                   ; wgt = temp*dens -> gamtot2  # 6538-6547
if (esv) ensure_single_val_fields_pass_r(gamtot, gamtot1, gamtot2)   (local g_lo ranges!)        # 6550-6554
if (.not.has_electron_species .or. .not.has_ion_species):          # 6557-6570
   yavg (iphi00=3):  gamtot(:,:,ik) += tite  for aky(ik) > epsilon(0.0)
   fieldlineavg:     gamtot += tite ; gamtot3 = (gamtot - tite)/gamtot ; where (gamtot3 < 2*eps) gamtot3 = 1
   default:          gamtot += tite
```

(Both `isgn` copies of `g_work` are equal, so the integral sums `x+x`.) `kwork_filter` is all-false when this
runs (level 28 < level 32).

### 3.11 Supercell / class structure

* Links: `itleft/itright` from `compute_itleft_and_itright(jshift0)` (dist_fn.fpp:2441-2515); `l_links/r_links`
  count the chain to the left/right, stopping at `-1` or a self-link (2608-2643). Supercell size
  `N(it,ik) = 1 + l_links + r_links`; cell position `n = 1 + l_links` (`kt2ki`, fields_implicit.f90:1136-1169).
* Classes: `N_class` = sorted unique values of `N(it,ik)` (n_k filled ik outer/it inner 2693-2700, quicksort 2705, unique 2709-2718);
  `M_class(j) = count(N==N_class(j))/N_class(j)` (2730). Non-linked: `i_class=1, N=1, M=naky*ntheta0` (2144-2150).
* Member ordering (`finish_fields_layouts`, fields_implicit.f90:1074-1133):

```
linked:  do c=1..i_class: m=1
           do it=1..ntheta0 ; do ik=1..naky          # it OUTER, ik INNER
              (cls, n) = kt2ki(ik,it); if (cls==c .and. n==1):
                 ik(m,1)=ik; it(m,1)=it; itr=it
                 do n=2..N_c: itr=itright(itr,ik); ik(m,n)=ik; it(m,n)=itr
                 m=m+1
other:   m=0; do it=1..ntheta0; do ik=1..naky: m=m+1; ik(m,1)=ik; it(m,1)=it      # m = ik + naky*(it-1)
```

* `f_lo(c)` world index `i ∈ [0, nindex_c*M_c)` (gs2_layouts.fpp:1650-1691, functions 1735-1849):
  `if_idx = 1+mod(i,nindex)`, `im_idx = 1+mod(i/nindex, M)`, `in_idx = 1+mod(i/nidx, N)`,
  `ig_idx = -ntgrid+mod(i/nf, nth)`, `ifield_idx = 1+mod(i,nf)`, `idx(if,im) = if-1+nindex*(im-1)`,
  `ij_idx_f(ig,if,n) = igp + nth*((if-1)+nf*(n-1))`. So within a member the f_lo index is **K** (0.1) and
  `ij_idx_f` is **B**.
* `jf_lo` (1855-1972): `ij_idx_jf(ig,if,ik,it) = igp + nth*((if-1)+nf*((it-1)+ntheta0*(ik-1)))`; only MPI
  bookkeeping besides the output scatter of 3.6.

### 3.12 Response matrix construction (fields_implicit.f90:459-719)

For each class `c` (processed `i_class` down to 1; order irrelevant, classes are independent and every
array is re-zeroed, 553-557):

```
am = 0 ; g = 0 ; phi=apar=bpar=0 ; phinew=aparnew=bparnew=0
do n = 1, N_c ; do ig = -ntgrid, ntgrid
   endpoint = (n > 1) .and. (ig == -ntgrid)                                          # 566-567
   f = 0
   for X in (phi, apar, bpar) if present:   f = f+1                                  # 572-633
      if endpoint: Xnew(ntgrid, it(m,n-1), ik(m)) = 1.0   for all m=1..M_c           # duplicate copy
      Xnew(ig, it(m,n), ik(m)) = 1.0                       for all m=1..M_c
      init_response_row(ig, f, am, c, n):
         timeadv(phi,apar,bpar, phinew,aparnew,bparnew, 0)     # 669: istep=0 (no explicit/NL terms, 5382;
                                                               #  hyper_diff iff hyper_in_initialisation, 4073;
                                                               #  collisions iff .not.split_collisions, 4082)
         getfieldeq(phinew,aparnew,bparnew, fieldeq,fieldeqa,fieldeqp)   # 670 (kwork_filter state applies)
         irow = f + nf*(igp + nth*(n-1))                     # 685  (1-based) = K(ig,f,n)+1
         do nn=1..N_c ; do m=1..M_c:
            it'=it(m,nn); ik'=ik(m,nn); iflo = idx(irow, m)
            am(nidx*(nn-1) + f' + nf*igp', iflo) = fieldeq_{f'}(ig', it', ik')   for all present f', ig'   # 694-712
      Xnew = 0
```

So `am(k+1, idx(j+1, m))` (0-based `k = K(ig',f',nn)`, `j = K(ig,f,n)`) is the residual of equation `k` of
member `m` for a unit perturbation `j`; perturbation `j` with `ig=-ntgrid, n>1` sets both copies of the shared
point, perturbation `j ∈ S_c` (`ig=+ntgrid, n<N`) sets only one copy (non-physical; that column is never used).
All `M_c` members are perturbed **simultaneously** (relevant if `hyper_diff` is amplitude dependent, P9).

### 3.13 Inversion `init_inverse_matrix` (fields_implicit.f90:722-1071) — exact transcription

Per class, `n = nindex_c`. Using `AMT[m, j, k] := am(k+1, idx(j+1,m))` (row `j` = response vector of perturbation
`j`) and `AIT[m, j, k] := a_inv(k+1, idx(j+1,m))`, `iskip(j) ≡ jskip(j) ≡ (j ∈ S_c)`
(for 1-based `i`: `N>1 ∧ i ≤ n−nf ∧ mod((i+nf−1)/nf, nth)==0 ∧ i>nf`, 770-774, 816-821 — equivalent to
`ig=ntgrid ∧ cell<N`), `zonal[m] ≡ (aky(ik)==0.0 .and. akx(it)==0.0)` exactly (834; `N=1` for such members):

```
AIT[m] = I                                                 # 760-762
for i in 0..n-1 (ascending):                               # 767
    if iskip(i): continue                                  # 774
    L  = AMT[:, i, :].copy()      # lhscol(:,m) = am(:, idx(i,m))   (778-807; broadcast or local copy: exact)
    Rr = AIT[:, i, :].copy()      # rhsrow(:,m) = a_inv(:, idx(i,m))
    for j in 0..n-1 (all members), not jskip(j):           # 814-821
        if not zonal[m]:                                   # 834
            fac = AMT[m,j,i] / L[m,i]                      # 836  complex ÷ (Smith, P2)
            AMT[m,j,i] = fac                               # 839
            AMT[m,j,k] = AMT[m,j,k] - L[m,k]*fac   (k ≠ i) # 842-843
            if j == i:  AIT[m,j,:] = AIT[m,j,:] / L[m,i]                    # 849-850  (elementwise ÷)
            else:       AIT[m,j,:] = AIT[m,j,:] - (Rr[m,:]*L[m,j]) / L[m,i] # 852-853  ((rhsrow*lhscol(irow))/lhscol(i), per element)
        else:
            AIT[m,j,:] = 0                                 # 856  (every pivot)
AIT[:, :, k] = 0  for k in S_c                             # 868-880  (columns of the inverse rows = coefficients of fl at duplicates)
AIT[:, j, :] = AIT[:, j+nf, :]  for j in S_c (ascending)   # 883-908  (row of (ntgrid,n) := row of (-ntgrid,n+1); j+nf ∉ S_c)
am = a_inv                                                 # 914
```

Then the reorder into the runtime storage (920-1063): for each output `j` (f_lo index → `ig,f,n` → `jflo`),
`aminv(c)%dcell(dj(c,jflo))%supercell(B(jg,jf,nn)+1) = AIT[m, j, K(jg,jf,nn)]`. In matrix form:

```
aminv_c[m, B(o), B(i)] = AIT[m, K(o), K(i)]        for all output points o and input points i of the supercell
```

The `L`, `Rr` copies are taken before the `j` loop, so within one pivot step all `j` updates are independent
(vectorisable exactly). GCC hoists the Smith ratio/denominator of `L[m,i]` out of the element loop
(same values, P2).

### 3.14 What 3.13 computes (verified)

Let `A_c,m` be the square matrix `A[k, j] = am(k, idx(j,m))` restricted to `k, j ∉ S_c` (size `nf*(1+2*ntgrid*N_c)`).
3.13 is a pivot-free Gauss–Jordan on `Aᵀ` with in-place multiplier storage; its result is
`AIT[m, j, k] = (A⁻¹)[j, k]` for `j,k ∉ S_c`, zero for `k ∈ S_c`, duplicated rows for `j ∈ S_c`. Verified with a
pure-Python transcription of 3.13 on random complex matrices (nf ∈ {1,2,3}, N ∈ {1,2,3}): `max|X·A − I| ≤ 3.4e-15`.
Hence `u = −aminv·fl` solves `A·u = −fl` on the reduced (single-valued) system, and the field at each duplicate
point is a bitwise copy of the field at its partner (before the `+phi` of 3.7).

### 3.15 Storage and dump/read of the response

* Runtime storage: `aminv(1:i_class)%dcell(1:#local rows)%supercell(1:nindex_c)` = one inverse row per output
  `jflo` (alloc 962-1000, `alloc_response_objects` 1716-1774); row layout **B** of the inputs.
* `dump_response_to_file_imp` (1172-1404): one netCDF file per supercell, name
  `<response_file>_ik_<ik>_is_<is>_dt_<code_dt formatted E14.6E3><suffix>` (fields_arrays.f90:37-58), `is` = per-ky
  supercell counter in leftmost-it order (1193-1234). Content: square complex `response(ax1, ax2)`
  (gs2_save.fpp:892-969; netCDF dims Fortran `[ri, ax1, ax2]` → NumPy `[ax2, ax1, ri]`), size
  `nf*(1+2*ntgrid*N)` (duplicates dropped): column `ax2 = icount` enumerates outputs in order
  (field outer, cell `in` following `itright` from the leftmost `it`, `ig` inner, skipping `ig=ntgrid` for `in<N`);
  `ax1 = cur_idx = ig_tmp + 2*ntgrid*(in_tmp−1) + (ifld_tmp−1)*(1+2*ntgrid*N)` (1350), same skipping. Also `dt`,
  `condition=-1`. This is **lossless** w.r.t. `aminv` (dropped inputs are zero columns, dropped outputs are copies).
* `read_response_from_file_imp` (1408-1713) inverts the mapping, re-creates duplicate rows by copying
  (1666-1703), and leaves zero coefficients at duplicate inputs (tmp_arr zeroed, 1590-1591).

---------------------------------------------------------------------------------------------------------------------

## 4. Options that change behaviour here

| option | namelist | default | effect | decision | reason |
|---|---|---|---|---|---|
| `field_option` = `default`/`implicit` | fields_knobs | `'default'` (= implicit, fields.f90:42-48,106) | algorithm of this spec | **SUPPORT** | reference path |
| `field_option` = `local`/`implicit_local` | fields_knobs | – | same math, different inverse/matvec round-off and response grouping (§7) | **SUPPORT** (separate replication path, §7.3) | widely used for flux-tube runs; until §7.3 is implemented the runner must refuse in strict mode |
| `field_option` = `gf_local` | fields_knobs | – | gf_lo decomposition + gf integration order | GATE | different integration grouping; needs `gf_lo_integrate` |
| `field_option` = `test` | fields_knobs | – | debugging solver | GATE | not physical |
| `field_subgath` | fields_knobs | F | allgatherv vs allreduce of `u` | SUPPORT (no-op) | both exact |
| `remove_zonal_flows_switch` | fields_knobs | F | §3.7 | SUPPORT | cheap |
| `force_maxwell_reinit` | fields_knobs | T | fields from g after reset/restart (§1.3) | SUPPORT | default behaviour |
| `dump_response`, `response_dir`, `response_file` | fields_knobs | F, '', '' | write files only | SUPPORT (no-op) | no numerics |
| `read_response` | fields_knobs | F | initial `aminv` from files (dt-matched by name) | SUPPORT | export whatever Fortran used; JAX rebuilds on dt change (P23) |
| `do_smart_update`, `minnrow`, `field_local_allreduce(_sub)`, `field_local_nonblocking_collectives`, `field_local_tuneminnrow` | fields_knobs | F,16,F,F,F,F | parallel decomposition of the local path | SUPPORT (no-op at 1 rank) | no numeric effect with 1 rank |
| `adiabatic_option` | dist_fn_knobs | `'default'` | adiabatic species term: default (`+tite`), fieldlineavg (`iphi00=2`, `fl_avg`), yavg (`iphi00=3`) | SUPPORT (all) | enters via exported `gamtot`/`gamtot3` + §3.4 |
| `esv` | dist_fn_knobs | F | §3.3 on antot*, gamtot* (only if linked) | SUPPORT | cheap scatter |
| `apfac` | dist_fn_knobs | 1.0 | factor in bpar residual | SUPPORT | scalar |
| `poisfac` | dist_fn_knobs | 0.0 | `+kperp2*poisfac` in gamtot | SUPPORT | via exported gamtot |
| `gridfac` | dist_fn_knobs | 1.0 | `gridfac1` at unconnected ends | SUPPORT | via exported gridfac1 |
| `afilter` | dist_fn_knobs | 0.0 | `antot *= exp(-afilter⁴kperp2²/4)` | GATE | debug knob; libm `exp` mismatch |
| `boundary_option` (default/zero/unconnected/periodic/self-periodic/kperiod=1/linked) | dist_fn_knobs | `'default'` (linked for box, zero otherwise, periodic at zero shear; dist_fn.fpp:1294-1310) | supercell structure, gridfac1, esv | SUPPORT | |
| `boundary_option='alternate-zero'` | dist_fn_knobs | – | "Ignored"; treated as default case in init_fieldeq | GATE | undocumented semantics |
| `hyper_in_initialisation` | dist_fn_knobs | T | hyper_diff inside the response build | SUPPORT | must replicate grouping (P9) |
| `gf_lo_integrate` | dist_fn_knobs | F | gf_lo integration order | GATE | different summation grouping |
| `def_parity`, `even` | dist_fn_knobs | F, T | enforce_parity inside timeadv (also in response) | SUPPORT (as per spec 07) | timeadv-side only; response follows automatically |
| `source_option='phiext_full'` | source_knobs | `'default'` | adds a field-independent source → response becomes affine | GATE | would corrupt the probe construction |
| `fphi`, `fapar`, `fbpar` (→ `has_*`) | knobs | 1, 0, −1 | which fields are solved (`> eps`) | SUPPORT | non-unit scale factors act in source terms (spec 02) |
| `beta`, `tite` | knobs | 0, 1 | coefficients | SUPPORT | |
| species `type` (ion/electron/hybrid/trace) | species_knobs | – | adiabatic detection; tracer weights → 0 | SUPPORT | |
| `layout` | layouts_knobs | `'lxyes'` | velocity summation order (§3.1) | SUPPORT | order table |
| `local_field_solve` | layouts_knobs | F | f_lo block layout (broadcast vs local copy) | SUPPORT (no-op) | exact copies |
| `intspec_sub` | layouts_knobs | F | sub-communicator integrals | SUPPORT (no-op at 1 rank) | forced off with 1 rank |
| `gf_local_fields` | layouts_knobs | F | gf layouts | GATE | with gf_local |
| `new_field_init` | init_g_knobs | T | t=0 fields via get_init_field (T) or `getfield` without `+phi` (F) | SUPPORT T / GATE F | F path is non-standard |
| antenna (`no_driver=.false.`) | driver namelist | no driver | `apar_ext` added before the solve | GATE | own spec needed |
| `split_collisions` | collisions_knobs | F | §3.9 after collisions | GATE | pending collisions spec |
| `vary_vnew` | collisions_knobs | F | changes vnmult without rebuilding aminv (P10) | GATE | pending collisions spec |
| `split_nonlinear` / `advance_nonadiabatic_dfn` | nonlinear / split_nonlinear knobs | F | §3.8/§3.9 inside NL substeps | GATE | pending spec 05 |
| `grid_option` (single/range/specified/box) | kt_grids_knobs | – | classes; `kwork_filter(1,1)` only for box | SUPPORT | |
| `skip_initialisation` (internal, gs2_optimisation.f90:79) | – | F | garbage response | GATE | benchmark only |
| MPI ranks > 1 / OpenMP build | run configuration | – | integral reduction order | GATE for strict validation | §6 P3 |

---------------------------------------------------------------------------------------------------------------------

## 5. GPU mapping

### 5.1 Data structures

* Per class `c` (static, from export): `N_c`, `M_c`, `ik_c[m]` (M,), `it_c[m,n]` (M,N) (0-based in JAX),
  `zonal_c[m]` (bool, exact `aky==0 and akx==0`), `aminv_c` complex `(M_c, R_c, R_c)`, `R_c = N_c*nidx`, axes in
  order **B**. Gather/scatter maps between field arrays `(naky, ntheta0, nth)` and supercell arrays
  `(M_c, N_c, nth)`: `X_c = X[ik_c[:,None], it_c, :]`. Every `(it,ik)` belongs to exactly one `(c,m,n)` ⇒ the
  scatter of 3.6 is a permutation.
* **Bucketing**: classes may be padded to a common `N_b` (e.g. power-of-two buckets) by appending cells whose
  `aminv` rows/columns are zero and whose gathered `fl` is zero: `S_n = +0` for a padded cell and `u − (+0) = u`
  exactly, so padding is exact. Memory grows ≤ 2× with power-of-two buckets; full padding to `N_max` is
  infeasible for box runs (Σ M_c R_c² vs (#supercells)·R_max²).

### 5.2 Velocity integrals (§3.1) — exact and fast

* Exact: `lax.scan` over the flattened velocity index (length `nspec*negrid*nlambda`), vectorised over
  `(ik,it,ig)` (code in 3.1); or a Pallas kernel: one program per tile of outputs, each thread loops
  sequentially over velocity points in layout order (one kernel, memory-bound, exact). Three integrals per field
  solve (phi, apar, bpar): fuse them in one pass over `gnew` (each keeps its own sequential accumulator).
* Fast (opt-in): `jnp.einsum` / `jnp.sum` (pairwise/tree) — ~1e-16 relative differences per integral.

### 5.3 Applying `aminv` (§3.6) — exact with depth `nidx + N`

The `S_n` of different cells are independent sums; only their combination is sequential:

```python
def apply_aminv_exact(aminv_c, fl_c):
    # aminv_c: (M, R, N, nidx) [B order: column index = n*nidx + f*nth + igp]; fl_c: (M, N, nidx) = blocked residuals
    M, R, N, nidx = aminv_c.shape
    def body(k, S):                                    # sequential over k (nidx steps), parallel over (m, r, n)
        return S + aminv_c[:, :, :, k] * fl_c[:, None, :, k]
    S = jax.lax.fori_loop(0, nidx, body, jnp.zeros((M, R, N), complex))   # S_n = ((0 + p1) + p2) + ...
    u = jnp.zeros((M, R), complex)
    for n in range(N):                                 # static unroll / fori_loop: u = ((0 - S1) - S2) ...
        u = u - S[:, :, n]
    return u                                           # u[m, B(ig,f,n)] -> field f at (ig, it_c[m,n], ik_c[m])
```

Complex products inside `body` must be the plain `(a.re*b.re − a.im*b.im, a.re*b.im + a.im*b.re)` (Tier A); see
Q1 for FMA. `fl_c` must be built from **all** `nth` theta points of every cell (the duplicate inputs carry zero
coefficients; any finite value is fine). Fast mode: `u = -jnp.einsum('mrc,mc->mr', aminv_flat, fl_flat)`.
Then `phinew = scatter(u_phi) + phi` etc. (§3.7).

### 5.4 Rebuilding the response after a dt change (§3.12) — batched probes

A *probe* is `(c, n, ig, f)`; its perturbation sets `X_f = 1` at `(ig, it_c[m,n], ik_c[m])` for **all** members `m`
(plus `(ntgrid, it_c[m,n−1], ik_c[m])` when `n>1, ig=−ntgrid`), all other fields 0, `g = 0`, old fields `0`.

```
for each class c (any order):
    probes = [(n, ig, f) for n in 1..N_c for ig in -ntgrid..ntgrid for f in present fields
              if not (N_c > 1 and ig == ntgrid and n < N_c)]          # S_c probes are never used -> may be skipped
    for batch in chunks(probes):                                       # vmap over the batch
        Xnew  = unit_perturbation(batch)                               # (B, nf, M_c, N_c, nth) in supercell view
        gnew  = timeadv_linear(g=0, old_fields=0, new_fields=Xnew, istep=0)   # spec 02/03/06/07 kernels
        fe    = getfieldeq(Xnew, gnew)                                 # §3.1-3.5, current kwork_filter
        AMT[m, K(ig,f,n), K(ig',f',nn)] = fe_{f'}[m, nn, ig']          # scatter rows
    if zonal member: its block is never used (rows of the inverse are 0) -> may be skipped
```

* `timeadv_linear(istep=0)` must reproduce GS2's istep=0 semantics: no explicit/nonlinear terms
  (dist_fn.fpp:5382, nonlinear_terms.fpp:424-427), `hyper_diff(gnew, phinew)` iff `hyper_in_initialisation`
  (4073-4079), collisions iff `.not.split_collisions` (4082-4086), `enforce_parity` if `def_parity` (4091),
  linked BCs, and the **current** `kwork_filter` (dist_fn.fpp:5401, 5475, 6000 skip filtered modes).
* Mode independence: every istep=0 operation acts on one `(it,ik)` (or one supercell for the linked BC)
  **except** amplitude-dependent hyperviscosity (`iso_shear` with `const_amp=F`, `aniso_shear`; hyper.f90:688-784)
  whose shear rate sums `|phinew|²` over all modes (ik outer, it inner). Therefore:
  - hyper off / `simple3D` / `const_amp=T`: probes of *different classes* at the same `(n,ig,f)` may be superposed
    in one full-grid `timeadv` (exactly what fields_local does), giving `nf*nth*N_max` timeadv calls in total;
  - amplitude-dependent hyper on: each probe must contain exactly GS2's pattern (all members of one class, nothing
    else), and the shear-rate sum must run over the full `(ik,it)` grid in Fortran loop order (zeros elsewhere add
    exactly 0) — see spec 07.
* Memory per probe: `M_c*N_c` modes × velocity grid × `2*nth`. Batch size chosen to fit memory.

### 5.5 Gauss–Jordan on GPU (§3.13)

Per class (or bucket), sequential over pivots, fully parallel over members and rows `j` inside a pivot:

```python
def gs2_gauss_jordan(AMT, skip, zonal, nf):      # AMT: (M, n, n) complex; skip: (n,) numpy bool; zonal: (M,) bool
    M, n, _ = AMT.shape
    AIT = jnp.broadcast_to(jnp.eye(n, dtype=AMT.dtype), (M, n, n))
    jj = jnp.arange(n)
    def step(i, carry):
        AMT, AIT = carry
        L, Rr = AMT[:, i, :], AIT[:, i, :]               # copies of pivot row (lhscol / rhsrow)
        Lii = L[:, i]
        fac = cdiv(AMT[:, :, i], Lii[:, None])           # Smith (P2); (M, n)
        AMn = AMT - L[:, None, :] * fac[:, :, None]
        AMn = AMn.at[:, :, i].set(fac)
        AIp = cdiv(AIT, Lii[:, None, None])               # used for j == i
        AIo = AIT - cdiv(Rr[:, None, :] * L[:, :, None], Lii[:, None, None])   # (rhsrow*lhscol(j))/lhscol(i)
        AIn = jnp.where((jj == i)[None, :, None], AIp, AIo)
        act = (~skip)[None, :] & (~zonal)[:, None]        # (M, n)
        AMT2 = jnp.where(act[..., None], AMn, AMT)
        AIT2 = jnp.where(act[..., None], AIn,
                         jnp.where(((~skip)[None, :] & zonal[:, None])[..., None], 0.0, AIT))
        keep = skip[i]                                    # iskip pivots are no-ops
        return (jnp.where(keep, AMT, AMT2), jnp.where(keep, AIT, AIT2))
    AMT, AIT = jax.lax.fori_loop(0, n, step, (AMT, AIT))
    AIT = AIT.at[:, :, skip].set(0.0)
    sk = jnp.nonzero(skip)[0]                             # static
    AIT = AIT.at[:, sk, :].set(AIT[:, sk + nf, :])
    return AIT                                            # then aminv_c[m, B(o), B(i)] = AIT[m, K(o), K(i)]
```

* `cdiv` must implement GCC's Smith expansion (P2) in real arithmetic; do **not** use a different inverse
  (LAPACK/`jnp.linalg.inv`) — pivot-free elimination order determines the round-off.
* Use `jnp.where` (not multiplicative masks) so NaN/Inf from zonal/skip lanes (e.g. `gamtot(1,1)=0` for kinetic
  electrons) never propagate.
* Depth: `n − |S_c|` pivots per class; work `O(M n³)`. Classes run one after another (or bucketed); this only
  happens on dt changes.
* Alternative for the initial dt: use the exported `aminv` directly (no rebuild error at t=0).

### 5.6 Masks and irregular structures

| mask / structure | where | GPU handling |
|---|---|---|
| `kwork_filter(it,ik)` | integrals (3.1), fl_avg (3.4), timeadv skips | `where(filter, 0, total)`; exported state; flips `(1,1)` after init (P6) |
| forbidden region (`wl=0`) | 3.1 | multiply by exact 0; ensure `g` finite there |
| tracer species | 3.1 | `wgt=0` |
| `aky(ik) > eps` | awgt (3.4), yavg | static masks |
| `aky==0 .and. akx==0` (exact) | 3.13 | `zonal_c[m]` |
| `aky==0.0` | gridfac1 (3.10) | exported |
| skip set `S_c` | 3.13, 3.6 | static index sets; zero columns, copied rows |
| linked chains | classes | per-class gathers `ik_c`, `it_c`; padding exact (5.1) |
| `itright<0` sentinel | esv (3.3) | precomputed scatter list of valid `(it,ik)` |
| `|denom| > eps` | 3.8 | `where` |

---------------------------------------------------------------------------------------------------------------------

## 6. Pitfalls (round-off and semantics)

* **P1 FMA.** The shipped binary contracts complex products in `getfield` (with *different* fusion choice in the
  2-wide body and in the odd remainder) and in the Gauss–Jordan (`fcmla`); it does **not** contract the
  accumulation in `integrate_species`. XLA (GPU) may contract where Fortran did not. Tier A ignores FMA; for
  bitwise agreement either rebuild the reference with `-ffp-contract=off` or emulate per site (Q1).
* **P2 Complex division = GCC "wide" Smith, `-fcx-fortran-rules`.** For `(a+ib)/(c+id)`:
  if `|d| > |c|` (strict): `r = c/d; den = c*r + d; re = (a*r + b)/den; im = (b*r − a)/den`;
  else: `r = d/c; den = d*r + c; re = (b*r + a)/den; im = (b − a*r)/den`. Observed in the binary: `den` fused in both
  branches; numerator fused in the `|d|≤|c|` branch only. `cdiv` in JAX must follow this branch rule; XLA's native
  complex divide is not guaranteed to. `x/x` is not exactly `(1,0)` (imaginary residual); replicate, don't
  "clean". Ratio/denominator are computed once per divisor — identical values, so per-element recomputation is fine.
* **P3 Integral summation order and parallelism.** Sequential, layout-dependent (3.1). With >1 MPI rank each rank
  sums its contiguous iglo block and `MPI_Allreduce` combines partials in an implementation-defined order; with
  OpenMP the `REDUCTION(+:total)` reorders. Validation references must be **1 rank, no OpenMP**.
* **P4 Per-cell grouping in `getfield`.** `u = ((0 − S_1) − S_2) − …`, each `S_n` summing all fields of cell `n`
  (field-major, theta-minor). A single long dot product, or per-field sums (the local path), round differently.
* **P5 Duplicate points.** `aminv` has zero coefficients at duplicate inputs and copied rows at duplicate outputs, so
  `u` is single-valued; `phinew = u + phi` is single-valued only if `phi` is (not enforced; `esv` acts on antot* /
  gamtot*, not on fields).
* **P6 `kwork_filter` timeline.** All-false during the first response build and during `init_fieldeq`;
  `(1,1)=.true.` (box) from the end of `init_fields` on, including during `get_init_field` at t=0 and during every
  rebuild after a dt change. Integrals return exact 0 for filtered modes; the `(1,1)` aminv block is all zero, so
  the `(1,1)` fields never change from their initial value.
* **P7 Zonal test variants.** GS2 implicit uses exact `aky(ik)/=0.0 .or. akx(it)/=0.0` (3.13); fields_local uses
  `abs(.)<epsilon(0.0)` (fields_local.f90:580); fl_avg uses `aky > epsilon(0.0)`; gridfac1 uses `aky == 0.0`.
  Non-box grids that contain a mode with `aky=akx=0` get a frozen field for it (Q7).
* **P8 `awgt` is computed once**, at the first `calculate_flux_surface_average` call (normally inside the first
  response build with an all-false filter; with `read_response` the first call is in `get_init_field`, where
  `kwork_filter(1,1)` is already true ⇒ `awgt(1,1)=0`, harmless). Never recomputed across dt changes.
* **P9 Hyperviscosity makes the response grouping-dependent.** `hyper_diff(gnew, phinew)` with amplitude-dependent
  shear rates during `istep=0` uses the unit perturbation pattern: implicit perturbs all members of one class at
  one `(n,ig)`; local perturbs one point of every supercell of every class. The two response matrices differ by
  more than round-off when such hyperviscosity is on. apar/bpar probes have `phinew=0` ⇒ filter `exp(−0)=1`.
* **P10 `vary_vnew`** changes `vnmult` and the collision coefficients (collisions.fpp:1892-1930) **without**
  rebuilding `aminv` — the stale response must be kept to match GS2.
* **P11 After a reset** fields come from `get_init_field` (3.8), not from the last solve, and `gnew` is restored from
  the saved g (spec 01). The rebuild runs with `kwork_filter(1,1)=.true.`.
* **P12 Absent fields.** `antot*`/`fieldeq*` of absent fields keep stale values; `phinew = phinew + phi` doubles an
  absent field (0 stays 0). The residual formulas still contain the absent field (multiplied, exactly zero).
* **P13 Association.** Keep `((antot + bpar*gamtot1) − (gamtot*gridfac1)*phi) + fl_avg` exactly; fields_local's
  `((antot − gamtot*gridfac1*phi) + fl_avg) + bpar*gamtot1` rounds differently when `has_bpar`.
* **P14 `esv` changes which value enters the solved equation.** The equation kept at a junction is the `ig=−ntgrid`
  one of the right cell; with `esv` its `antot` is overwritten by the left cell's `ig=+ntgrid` value.
* **P15 `ensure_single_val_fields_pass_r` on gamtot*** uses local `g_lo` ranges (not `force_all`): with several
  ranks different ranks can hold different gamtot. Exported gamtot (1 rank) is the reference.
* **P16** `epsilon(0.0)` is the double epsilon (2.22e-16) everywhere here.
* **P17** Fortran `sum` in `getfield` is sequential (verified); `jnp.sum`, `jnp.dot`, `einsum`, cuBLAS are not.
* **P18 Forbidden-region `g`.** Contributes `0*s`; only exact if `s` is finite. `zero_forbid=.false.` (default) leaves
  whatever `invert_rhs_1` produced there; the JAX timeadv must keep it finite (ideally identical).
* **P19 MPI in the solve itself is exact** (each `u(jflo)` nonzero on one rank; allreduce adds zeros).
* **P20 Response build zeroes g and all fields** (553-557); after the build `gnew` holds the last probe's response.
  JAX must build with separate buffers and leave the evolving state untouched.
* **P21** The perturbation enters `getfieldeq` as "phi" (so `−gamtot*gridfac1*δ` is part of the column) and as
  `phinew` in `timeadv`, with `phi_old = 0`.
* **P22** `x/x` and similar near-identities appear in pivot columns of `am`; those residuals never reach `a_inv`
  (3.13), but an implementation that "simplifies" them changes nothing either — keep the transcription literal.
* **P23 `read_response`** matches files by `code_dt` printed with 7 significant digits; a file for a slightly
  different dt would be accepted by Fortran.
* **P24** `remove_zonal_flows` assumes `ik=1` is ky=0 (fields_implicit.f90:384).
* **P25** The order of field slots is phi, apar, bpar among present fields, both in `am` (interleaved K) and in
  `aminv`/`fl` (blocked B).
* **P26 Self-periodic zonal modes** (box, ky=0) are `N=1` supercells: the field solve treats `ig=±ntgrid` as
  independent unknowns; periodicity comes only from the g boundary condition and `esv`.
* **P27** `get_fields_direct_from_dfn` ignores `gridfac1` and `fl_avg` is the only adiabatic coupling there.
* **P28** `advance_implicit` copies all three fields and adds all three; `advance_local` only present ones —
  identical when absent fields are zero.
* **P29 Antenna**: `aparnew = apar + apar_ext` enters only the residual, then `aparnew = u + apar` (GATEd).

---------------------------------------------------------------------------------------------------------------------

## 7. `field_option='local'` versus `'implicit'`

### 7.1 Identical

* The linear system and its reduction: same unknowns (duplicate inputs dropped: `ignore_boundary` for cells
  `ic<ncell`, fields_local.f90:660, `ncol=2*ntgrid`, 474-475), same equations kept (`ig=−ntgrid` of the right
  cell), same perturbation of both copies at junctions (2161-2171), same `timeadv(…,0)`.
* Response entries: computed per mode with the same arithmetic as implicit (getan_nogath ≡ getan for 1 rank, 3.2),
  **except** the phi residual association when `has_bpar` (P13) and grouping effects of amplitude-dependent hyper (P9).
* Zonal `(kx=ky=0)` cell: treated as empty (c_set_locality, fields_local.f90:580) → its field update is 0 (gather fills zeros),
  as in implicit.

### 7.2 Different (all at round-off level unless noted)

| aspect | implicit | local |
|---|---|---|
| matrix order | interleaved K, per member block | field-major extended domain `iex = (ig+ntgrid+1)+(ic−1)*2*ntgrid+(ifl−1)*nextend` (1113) |
| response grouping | all members of one class at one `(n,ig,f)` | one point in every supercell of all ky at once, `kwork_filter` set for finished supercells (2081-2273) — **differs beyond round-off with amplitude-dependent hyper** |
| residual | getfieldeq (3.5) | fm_getfieldeq1_nogath: `((antot − gamtot*gridfac1*phi) + fl_avg) + bpar*gamtot1` (3423-3426); bpar eq `bpar*gridfac1 + (…)*(…)` (3435-3437; commutative ⇒ identical) |
| inversion | pivot-free GJ on Aᵀ (3.13) | `invert_serial` → `invert_block` (matrix_inverse.fpp:337-391, no LAPACK in this build), recursive 2×2 block inverse with `invert_gauss_jordan` (266-295) for n≤4, then `transpose` (1263-1264) |
| matvec | per output: `((0−S_1)−S_2)…`, `S_n` over all fields of cell n | per (cell, field): `matmul(vect, data)` (379), cell: `tmp=0; tmp−=S_phi; tmp−=S_apar; tmp−=S_bpar` (538-544), supercell: `local=0; local+=tmp_cell` over cells (747-750) |
| update | `phinew = u + phi` all fields | only present fields |
| exb trigger | `allocated(kx_shift/theta0_shift)` | `abs(g_exb*g_exbfac) > eps` (equivalent) |

### 7.3 What an exact replication of 'local' requires

1. Response columns: same probes as 5.4 but with fields_local's grouping (needed only with amplitude-dependent hyper),
   residual with the association of 7.2.
2. Inverse: port `invert_block` literally (static recursion on sizes; `m = max(n/2,1)`; `tmp = −matmul(bb,dd)`;
   `aa = aa + matmul(tmp,cc)`; invert `aa`; `bb = matmul(aa,tmp)`; `tmp = −matmul(dd,cc)`; `cc = matmul(tmp,aa)`;
   `dd = dd + matmul(tmp,bb)`), `invert_gauss_jordan` for `n≤4` (`fac=1/a(i,i); a(i,i)=1; a(:,i)=a(:,i)*fac;
   for k≠i: tmp=a(i,k); a(i,k)=0; a(:,k)=a(:,k)−a(:,i)*tmp`), each `matmul` as a sequential-in-k sum from 0
   (assumed for both gfortran's inlined matmul and libgfortran `matmul_c8`; Q3), negation applied after the sum.
3. Matvec and accumulation exactly as in 7.2; duplicate outputs filled by copy (2689-2694).
Alternative accepted by this spec: run the implicit algorithm for a 'local' input and compare with a looser,
documented tolerance (≈ cond(A)·1e-16), or export the local inverse via `dump_response` for the initial dt.

---------------------------------------------------------------------------------------------------------------------

## 8. Export checklist for this subsystem (names are suggestions for the export program)

Scalars: `ntgrid, ntheta0, naky, nlambda, negrid, nspec, nfield, has_phi, has_apar, has_bpar, beta, tite, apfac,
poisfac, gridfac, afilter, esv, adiabatic_option_switch, boundary_option_switch, has_electron_species,
has_ion_species, layout, code_dt, i_class, field_line_average_weight, gridopt_switch, hyper_in_initialisation`.

Arrays (Fortran shape → NumPy shape):

| name | Fortran | NumPy |
|---|---|---|
| `N_class`, `M_class` | `(i_class)` | `(i_class,)` |
| `cls_ik`, `cls_it` (padded −1) | `f_lo(c)%ik(m,n)` as `(Mmax, Nmax, i_class)` | `(i_class, Nmax, Mmax)`: `cls_it[c,n,m] = f_lo(c+1)%it(m+1,n+1)` |
| `aminv_c<c>` (one per class) | `(R_c, R_c, M_c)` complex, `aminv(col,row,m)` | `(M_c, R_c, R_c)`: `[m, B_out, B_in]` |
| `itleft`, `itright`, `l_links`, `r_links` | `(ntheta0, naky)` | `(naky, ntheta0)` |
| `gamtot`, `gamtot1`, `gamtot2`, `gamtot3`, `gridfac1`, `kperp2`, `inv_kperp2`, `inv_phi_denominator_g`, `inv_bpar_denominator_g` | `(-ntgrid:ntgrid, ntheta0, naky)` | `(naky, ntheta0, nth)` |
| `awgt`, `kwork_filter` (after init_fields) | `(ntheta0, naky)` | `(naky, ntheta0)` |
| `aj0`, `aj1`, `vperp2` | canonical `(ig,it,ik,il,ie,is)` | `(ns, ne, nl, nky, nkx, nth)` |
| `vpa` | canonical g-like | `(ns, ne, nl, nky, nkx, 2, nth)` |
| `w` | `(negrid, nspec)` | `(nspec, negrid)` |
| `wl`, `forbid` | `(-ntgrid:ntgrid, nlambda)` | `(nlambda, nth)` |
| `bmag`, `theta`, `jacob` | `(-ntgrid:ntgrid)` | `(nth,)` |
| `aky`, `akx` | `(naky)`, `(ntheta0)` | same |
| `spec_z, spec_dens, spec_temp, spec_mass, spec_stm, spec_type` | `(nspec)` | `(nspec,)` |
| t=0 `phi, apar, bpar` (= `*new`) | field-like | `(naky, ntheta0, nth)` complex |

How to obtain `aminv` (private in fields_implicit): either (a) add a public read-only accessor in the export build,
iterating `c, m, n, f, ig` → `jflo = ij_idx(jf_lo, ig, f, ik_c(m), it_c(m,n))`, `dc = jf_lo%dj(c,jflo)`, row =
`aminv(c)%dcell(dc)%supercell(:)`; or (b) call the public `fields::dump_response_to_file` and reassemble from the
per-supercell files (lossless, §3.15). Optional validation exports: `am` before inversion (needs a hook) and the
post-first-`timeadv` residual `fl` and `u` of step 1.

---------------------------------------------------------------------------------------------------------------------

## 9. Open questions

1. **Q1 FMA policy** (shared with specs 01/02/05/07): validate against the `-O3` FMA-contracted binary at Tier A
   tolerance, or rebuild the reference with `-ffp-contract=off` and control XLA contraction on the target GPU?
2. **Q2** Is exact replication of `field_option='local'` (§7.3) required, or is "implicit algorithm + looser
   tolerance" acceptable for 'local' inputs?
3. **Q3** libgfortran `matmul_c8` (blocked path) and the `-O3` inline matmul are assumed to sum sequentially in k from
   0 (with whatever FMA use); needs a micro-test before porting 7.3.
4. **Q4** XLA's complex multiply/divide lowering on GPU (Smith? FMA?) — this spec mandates explicit real-arithmetic
   implementations; confirm the performance cost is acceptable.
5. **Q5** Can the JAX `timeadv` operate on a subset of modes (supercell view of one class) for probe batching, or only
   on the full grid (then use cross-class superposition when hyper is linear, 5.4)?
6. **Q6** `vary_vnew` keeps a stale `aminv` (P10): confirm we replicate rather than "fix".
7. **Q7** Non-box grids containing an `aky=akx=0` mode (e.g. single/range with `aky=0, theta0=0`): GS2 freezes that
   mode's fields at their initial value; replicate (assumed yes).
8. **Q8** Multi-rank references: is matching the MPI_Allreduce order (OpenMPI algorithm selection) in scope, or are
   all references single-rank?
9. **Q9** Is any istep=0 operation other than amplitude-dependent hyperviscosity mode-coupling? From the code read
   here none is (collisions, invert_rhs, parity, fl_avg are per mode / per supercell); a probe-superposition test in
   the JAX port should confirm bitwise independence.
10. **Q10** The pre-inversion `am` cannot be exported without a code hook (it is local to `init_response_matrix` and
    overwritten in place); decide whether the export phase adds such a hook for validating the JAX probe builder.
