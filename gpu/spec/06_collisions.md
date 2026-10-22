# 06 — COLLISIONS subsystem: port contract (GS2 → JAX, float64/complex128)

Status: read-only analysis of GS2 8.2.1 (`gpu-port` branch, commit 4d8c94bcf).
Sources read in full: `src/collisions.fpp` (4230 lines), `src/le_derivatives.fpp` (419 lines).
Sources read in part: `src/gs2_layouts.fpp` (le_lo / e_lo / lz_lo / g_lo ordering), `src/le_grids.f90`
(maps, integrals, xi grid, jend/forbid, weights), `src/dist_fn.fpp` (timeadv, collisions_advance,
init order), `src/dist_fn_arrays.fpp` (g_adjust), `src/fields_implicit.f90`, `src/fields_local.f90`,
`src/gs2_reinit.f90`, `src/gs2_init*.inc`, `src/gs2_time.f90`, `src/run_parameters.f90`,
`externals/utils/redistribute.f90`, `src/diagnostics/diagnostics_velocity_space.fpp`.

Every `file:line` below refers to these files. "C:" means `src/collisions.fpp`, "LD:" means
`src/le_derivatives.fpp`, "LG:" means `src/le_grids.f90`, "GL:" means `src/gs2_layouts.fpp`,
"DF:" means `src/dist_fn.fpp`, "DFA:" means `src/dist_fn_arrays.fpp`.

---------------------------------------------------------------------------------------------------

## 0. Executive summary

* The collision operator is a **linear, implicit (backward-Euler, full `code_dt`) operator applied to
  h = g + (q φ/T) J0 F0 + … (after `g_adjust`)** at the end of every `timeadv` call
  (`DF:4083-4087`). Order inside one application (`C:3050-3089`):
  1. energy diffusion tridiagonal solve in energy (`solfp_ediffuse_le_layout`), then
  2. energy-diffusion momentum/energy restoring terms (`conserve_diffuse_le_layout`), then
  3. electron–ion drag (resistivity) explicit term (only if `drag`), then
  4. pitch-angle (Lorentz) tridiagonal solve in ξ (`solfp_lorentz_le_layout`), then
  5. Lorentz momentum/energy restoring terms (`conserve_lorentz_le_layout`).
* **Default data path is `use_le_layout = .true.`** (`C:349`). In that path all solves and conservation
  integrals for one spatial point `(ig, it, ik, is)` are local and independent ⇒ embarrassingly parallel
  over `(is, ik, it, ig)`; sequential only along `ie` (energy solve) and `ixi` (pitch-angle solve).
* **dt-dependent arrays** (must be rebuilt whenever `code_dt` changes, and when `vnmult` changes):
  `c1le, betaale, qle` (Lorentz LU), `ec1le, ebetaale, eqle` (energy LU), and the six conservation
  vectors `z0le, s0le, w0le, bz0le, bs0le, bw0le`. Everything else the step reads is dt-independent.
  Recommendation: **BOTH** — export them at init (for validation and for the initial dt) and port the
  exact rebuild formulas (§5) so JAX/NumPy can rebuild on a dt change. All transcendental
  sub-expressions are removable from the rebuild by exporting two helper tables
  (`e15 = energy_maxw**1.5`, `capg_half`), making the rebuild purely `+ - * / sqrt` (IEEE-exact
  given the same operation order).
* **Confirmed GS2 memory bug (no trapped particles)**: when the λ grid has no trapped particles
  (`nlambda == ng2`, e.g. `eps = 0` or `trapped_particles = F`) several collision arrays are indexed one
  past their first extent (`nxi_lim = 2*ng2` but the solve touches row `2*ng2+1`). Verified empirically
  with `bin/gs2` (single ky, ntheta=12, 40 steps, see §8.1):
  * trapped grid (eps=0.18): `use_le_layout=T` vs `F` agree to 3.2e-16 (conserve_moments=T) and
    bit-exactly (conserve_moments=F);
  * no-trapped grid (eps=0): `use_le_layout=T` vs `F` differ by **21 % in φ** (conserve_moments=F),
    and `use_le_layout=T, conserve_moments=T` **aborts with heap corruption** (SIGTRAP in malloc).
  ⇒ GATE: (no trapped particles) ∧ (Lorentz operator active) ∧ (`conserve_moments` ∨ ¬`use_le_layout`).
* `vary_vnew` (adaptive collisionality) is GATED: it rebuilds collision matrices from diagnostics
  without rebuilding the field response matrix.
* `conservative = .false.` is GATED: `init_vnew` leaves `vnew_E(:,1,:)`, `vnew_E(:,negrid,:)`,
  `delvnew(:,1,:)`, `delvnew(:,negrid,:)` **uninitialised** (`C:1837-1845`).

---------------------------------------------------------------------------------------------------

## 1. Call graph and control flow

### 1.1 Initialisation

```
gs2_init level "collisions" (=30, src/gs2_init_level_list.inc)          src/gs2_init_subroutines.inc:362-372
  └─ init_collisions                                                     C:436-463
       ├─ init_gs2_layouts, init_species
       ├─ hyper_colls = any(spec%nu_h > epsilon(0.0))                    C:451-452
       ├─ init_theta_grid, init_kt_grids, init_le_grids, init_run_parameters, init_dist_fn_layouts
       ├─ read_parameters                                                C:466-562
       │    ├─ copy namelist collisions_knobs into module vars           C:496-523
       │    ├─ collision_model → collision_model_switch                  C:529-531, options C:475-481
       │    ├─ lorentz_scheme → lorentz_switch; ediff_scheme → ediff_switch   C:533-539
       │    ├─ has_lorentz / has_diffuse                                 C:541-550
       │    ├─ drag = has_lorentz ∧ resistivity ∧ beta>eps ∧ nspec>1 ∧ fapar>0   C:552-553
       │    └─ nxi_lim = nxi+1 if nxi > 2*ng2 else nxi                   C:557-561
       └─ init_arrays                                                    C:576-631
            ├─ if model == none: colls=F, return                          C:599-602
            ├─ init_vnew                                                 C:1707-1876
            ├─ if all(|vnew(:,1,:)| <= 2 eps) ∧ ¬force_collisions: model=none, colls=F, return  C:605-609
            ├─ if heating: allocate/zero c_rate                          C:611-614
            ├─ init_map(use_lz_layout, use_e_layout, use_le_layout, test) LG:3323-3364
            ├─ if has_lorentz: init_lorentz                              C:2185-2358
            │      ├─ pitch_weights (once)                               C:2199-2210
            │      ├─ init_vpdiff                                        C:4079-4155
            │      └─ per (ile,ie): get_lorentz_matrix (C:2392-2639) [+ set_hzero_lorentz_collisions_matrix C:2366-2389]
            │                       forward elimination → c1le, betaale, qle
            │   if conserve_moments: init_lorentz_conserve               C:691-1194
            ├─ if has_diffuse: init_ediffuse                             C:1945-2038
            │      └─ per (ile,ixi): get_ediffuse_matrix (C:2041-2182), forward elimination → ec1le, ebetaale, eqle
            │   if conserve_moments: init_diffuse_conserve               C:1198-1679
            └─ if use_le_layout ∧ (conserve_moments ∨ drag): init_le_bessel   C:634-687
```

Prerequisites that are already initialised before level 30: `aj0, aj1` (`init_bessel`, DF:1756-1823,
dist_fn level 2), `vpa` (`init_vpa_vpac`, DF:1592-1678, dist_fn level 1), `vperp2`. `init_collisions`
is called again from `init_dist_fn_level_3` (DF:991) but returns immediately (`initialized` flag, C:446).

### 1.2 Re-initialisation when dt changes

`gs2_reinit.f90:104-225 reset_time_step` → `init(current_init, init_level_list%override_timestep)` (level
29) tears down level 30 (`collisions`) and above: `finish_collisions` (C:4165-4212) deallocates every
collision array and **saves `vnm_init = vnmult; vnmult = -1`** (C:4174-4175). Going back up,
`init_collisions` rebuilds everything with the new `code_dt`; `init_vnew` restores `vnmult = vnm_init`
(C:1731-1733). Then `dist_fn_level_3` and `fields_level_2` (response matrices) are rebuilt.
Note: with `delt_option = 'check_restart'` `init_vnew` first overwrites `vnm_init` from the restart file
(C:1723-1726) — on a mid-run dt reset this re-reads the file value (equals 1 unless `vary_vnew` was used).

### 1.3 Time step (one call of the collision operator)

```
advance_implicit (src/fields_implicit.f90:315-371)  [same structure in fields_local.f90:4296-4340]
  ├─ timeadv(phi,apar,bpar, phinew,aparnew,bparnew, istep)            # 1st call, phinew==phi
  │     invert_rhs → hyper_diff → [if ¬split_collisions] vspace_derivatives(gnew, g, g_work, phi, bpar, phinew, bparnew, mode)   DF:4083-4087
  ├─ getfield; phinew += phi …
  ├─ timeadv(..., istep, diagnostics=1)                               # 2nd call: heating diagnostics only here
  │     invert_rhs → hyper_diff → [if ¬split] vspace_derivatives(..., mode=1)
  └─ collisions_advance(phi,bpar,phinew,aparnew,bparnew,diagnostics)  DF:4130-4167 (no-op unless split_collisions)
         vspace_derivatives(...) ; g_adjust(from) ; calculate_potentials_from_nonadiabatic_dfn ; g_adjust(to)

vspace_derivatives (LD:15-174), use_le_layout branch (LD:67-112):
  if colls (or LOWFLOW build):
    if adjust: g_adjust(gnew, phinew, bparnew, from_g_gs2)     # g → h      DFA:101-175
               [if heating_flag: g_adjust(gold, phi, bpar, from_g_gs2)]
    [if heating_flag: gc3 = copy(gnew)]
    gle(1:nxi+1, 1:negrid+1, le_lo) = 0                          LD:79-80
    gather(g2le, gnew, gle)                                      LD:83   (pure permutation, allowed points only)
    solfp1_le_layout(gle, diagnostics)                           C:3035-3090
        has_diffuse: solfp_ediffuse_le_layout (C:4037-4076); if conserve_moments: conserve_diffuse_le_layout (C:3541-3628)
        has_lorentz: if drag: drag term (C:3056-3085)
                     solfp_lorentz_le_layout (C:3794-3981); if conserve_moments: conserve_lorentz_le_layout (C:3256-3376)
    scatter(g2le, gle, gnew)                                     LD:94   (allowed points only)
    [heating_flag: c_rate(:,:,:,:,3) from gc3 — diagnostic]
    if adjust: g_adjust(gnew, phinew, bparnew, to_g_gs2)        # h → g
               [if heating_flag: g_adjust(gold, phi, bpar, to_g_gs2)]
```

`heating_flag = heating .and. present(diagnostics)` (LD:65); `diagnostics` is present only in the 2nd
`timeadv` call (fields_implicit.f90:365) and in `collisions_advance`. Optional args `gtoc/ctog` are
never passed by any caller ⇒ always `.true.`.

The **field response matrix** is built by `init_response_row` calling `timeadv(..., istep=0)`
(fields_implicit.f90:669, fields_local.f90:2184): collisions **are included** in the response (only
`hyper_diff` is excluded at istep 0, DF:4073), unless `split_collisions`. Any JAX rebuild of response
matrices on a dt change must therefore include the (rebuilt) collision operator.

### 1.4 Adaptive collisionality (GATED)

`get_verr` (diagnostics_velocity_space.fpp, call at line 510) always calls `adjust_vnmult(errest, …)`
(C:1892-1942), which returns immediately unless `vary_vnew`. If active, it calls
`init_ediffuse(vnmult_target)` + `init_diffuse_conserve` and `init_lorentz(vnmult_target)` +
`init_lorentz_conserve` — **even when vnmult is unchanged and even for operators not in use**
(`has_diffuse` not checked). The response matrix is **not** rebuilt. Cadence: every `ncheck` steps
(gs2_diagnostics.fpp:1113-1114; new diagnostics use their own `ncheck`, gs2_diagnostics_new.f90:509-513).

---------------------------------------------------------------------------------------------------

## 2. Index sets, layouts and the (il, isgn) ↔ ixi map

### 2.1 Grid integers

| symbol | definition | source |
|---|---|---|
| `ng2` | `2*ngauss` or `npassing` | LG:402-404 |
| `nlambda` | `ng2 + nbset` if trapped particles, else `ng2` | LG:1896-1902 |
| trapped grid | `nlambda > ng2` (`grid_has_trapped_particles`) | LG:443-446 |
| `nxi` | `max(2*nlambda-1, 2*ng2)` | LG:1904 |
| `nxi_lim` | `nxi+1` if `nxi > 2*ng2` (trapped grid) else `nxi`; i.e. `2*nlambda` (trapped) or `2*ng2` (no trapped) | C:557-561 |
| `NX` (runtime `gle` rows) | `nxi+1` = `nxi_lim` (trapped) or `2*ng2+1` (no trapped) | LD:79 |
| `jend(ig)` | 0 if no trapped grid; else `#{il : 1 - al(il)*bmag(ig) > -bouncefuzz}` (≥ ng2+1) | LG:2871-2903 |
| `forbid(ig,il)` | `1 - al(il)*bmag(ig) < -bouncefuzz` | LG:2715-2738 |
| wfb | `il == ng2+1` (`il_is_wfb`), passing `il <= ng2`, trapped `il > ng2+1` | LG:450-491 |

`al` is strictly increasing in `il` (passing from Gauss–Legendre/Radau points 1→0 mapped by
`(1-x²)/bmax`, LG:2209; trapped `1/bset`); hence the allowed set at `ig` is `il = 1..jend(ig)` (trapped
grid) or all `il` (no trapped grid). At the `bmax` theta points `jend = ng2+1` (the wfb bounces there).

### 2.2 Layout semantics (only the index sets matter for the port)

* `g_lo`: index `iglo` ↔ `(it, ik, il, ie, is)`, order given by the `layout` string, **first letter
  fastest** (GL:917-930; default `'lxyes'`: il fastest, then it, ik, ie, is slowest). Arrays
  `(-ntgrid:ntgrid, 2, iglo)`.
* `le_lo`: index `ile` ↔ `(ig, it, ik, is)`; ig fastest, then it/ik per layout, is slowest
  (GL:2003-2114, `idx_le` GL:2160-2170). Arrays `(ixi, ie, ile)`; **all of (ξ, E) is local** to one
  `ile`. This is the default collision layout.
* `lz_lo` (standard layout, Lorentz only): `ilz` ↔ `(ig, it, ik, ie, is)` (GL:2602-2691); arrays `(ixi, ilz)`.
* `e_lo` (standard layout, energy only): `ielo` ↔ `(ig, isgn, it, ik, il, is)` (GL:2438-2497); arrays `(ie, ielo)`.

### 2.3 The ξ index `ixi` (identical for g2le and lambda_map)

`ixi_to_il(ixi,ig)`, `ixi_to_isgn(ixi,ig)` (LG:2972-3011); the gather direction used by `g2le`
(LG:3519-3541) and `lambda_map` (LG:4127-4147):

* No trapped grid (`je = jend(ig) = 0`): `isgn=1 → ixi = il` (1..ng2); `isgn=2 → ixi = 2*nlambda+1-il`
  (ng2+1..2ng2). All `2*ng2` entries filled.
* Trapped grid (`je = jend(ig) ≥ ng2+1`):

```
 ixi :  1 ..  je-1 | je   | je+1 .. 2je-1 | 2je  | 2je+1 .. nl+je | nl+je+1 .. 2nl
 il  :  1 ..  je-1 | je   | je-1 .. 1     | je   | je+1  .. nl    | nl .. je+1
 isgn:  1          | 2    | 2             | 1    | 1              | 2
                     ^ bounce point, isgn=2      ^ duplicate bounce point, isgn=1       (forbidden region)
```
  i.e. `il < je`: isgn1→`il`, isgn2→`2je-il`; `il == je`: **isgn1→`2je`, isgn2→`je`**; `il > je`:
  isgn1→`il+je`, isgn2→`2nl+1-il+je`.
* `xi(ixi,ig)`: `xi(1:je) = sqrt(max(1-al*bmag,0))`, `xi(je+1:2je-1) = -xi(je-1:1:-1)`, rest 0 (LG:2965-2969).
* `wxi(ixi,ig) = wl(ig, ixi_to_il(ixi,ig))` (LG:3014-3020); `wl(ig,il)` is 0 at forbidden points.

**g2le** (`setup_g2le_redistribute_local`, LG:3379-3599) **skips forbidden points**
(`skip_forbidden_region = .true.`, LG:3403, 3422, 3487, 3523): `gather` writes only allowed entries
(the destination is zeroed beforehand, LD:80), `scatter` writes back only allowed entries (forbidden
entries of `g` keep their values). The redistribute is a pure element copy (redistribute.f90:6425-6480);
there is no arithmetic. `lambda_map`/`energy_map` (standard layout) map **all** points, incl. forbidden.

### 2.4 Canonical export conventions for collision-specific arrays

* g-like arrays (`aj0, aj1, vpa, vperp2, s0, w0, z0, bs0, bw0, bz0`, the field-shaped
  `aparnew/phinew/bparnew`): use the project canonical convention (g: Fortran
  `(ig, isgn, it, ik, il, ie, is)` → NumPy `(nspec, negrid, nlambda, naky, ntheta0, 2, nθ)`; `aj0/aj1/vperp2`
  have no isgn: Fortran `(ig, it, ik, il, ie, is)` → NumPy `(nspec, negrid, nlambda, naky, ntheta0, nθ)`).
* **LE-canonical** (all `*le` arrays): expand `ile → (ig, it, ik, is)` via `ig_idx/it_idx/ik_idx/is_idx(le_lo, ile)`:
  Fortran `(ixi=1:nxi_lim, ie=1:NE, ig=-ntgrid:ntgrid, it=1:ntheta0, ik=1:naky, is=1:nspec)` → NumPy
  `(nspec, naky, ntheta0, nθ, NE, nxi_lim)`, with `NE = negrid` for LU arrays
  (`c1le, betaale, qle, ec1le, ebetaale, eqle`) and `NE = negrid+1` for vectors
  (`s0le, w0le, z0le, bs0le, bw0le, bz0le, aj0le, vperp_aj1le, vpa_aj0_le`).
* Velocity/θ arrays keep Fortran order: `pitch_weights(nlambda, -ntgrid:ntgrid)`,
  `vpdiff(-ntgrid:ntgrid, 2, nlambda)`, `vpdiffle(nxi_lim, -ntgrid:ntgrid)`,
  `vnew*(naky, negrid, nspec)`, `vnewh(-ntgrid:ntgrid, ntheta0, naky, nspec)`, `xi/wxi/ixi_to_*(2*nlambda, -ntgrid:ntgrid)`,
  `wl(-ntgrid:ntgrid, nlambda)`, `w(negrid, nspec)`, `kperp2(-ntgrid:ntgrid, ntheta0, naky)`.
  State in the netCDF attributes that these are Fortran order (NumPy reads them reversed).

---------------------------------------------------------------------------------------------------

## 3. Data inventory (everything the collision step and its dt-rebuild read)

Legend: **dt** = depends on `code_dt`; **vnm** = depends on `vnmult`; decision EXPORT / RECOMPUTE / BOTH.
All collision module arrays are **PRIVATE** (only the names in C:17-23 are public) — the exporter needs
an accessor added to `collisions.fpp` (see §10).

### 3.1 Collision-module arrays

| name | shape (Fortran) | computed in | dt | vnm | decision / notes |
|---|---|---|---|---|---|
| `vnew` | `(naky,negrid,nspec)` | `init_vnew` C:1761-1799 | no | no | EXPORT (erf/exp/pow). Also used as the per-(ik,is) "skip" test `abs(vnew(ik,1,is)) < 2 eps` |
| `vnew_s` | same | C:1780-1782 (0 if ei_coll_only) | no | no | EXPORT; used in ediffuse `ee` |
| `vnew_D` | same | C:1788-1796 | no | no | EXPORT; Lorentz conservation |
| `vnew_E` | same | C:1801-1845 | no | no | EXPORT; ediffuse conservation |
| `delvnew` | same | C:1801-1845 | no | no | EXPORT; ediffuse conservation |
| `vnewh` | `(-ntgrid:ntgrid,ntheta0,naky,nspec)` | C:1853-1863 (if hyper_colls) | no | no | EXPORT (or recompute `nu_h*kperp2**2/k4max`) |
| `pitch_weights` | `(nlambda,-ntgrid:ntgrid)` | C:2199-2210 | no | no | EXPORT. ≠ `transpose(wl)` in general (always old-FD trapped weights, even if `new_trap_int`) |
| `vpdiff` | `(-ntgrid:ntgrid,2,nlambda)` | `init_vpdiff` C:4088-4140 | no | no | EXPORT (only if `conservative`) |
| `vpdiffle` | `(nxi_lim,-ntgrid:ntgrid)` | C:4142-4152 | no | no | EXPORT (le + conservative) |
| `c1le, betaale, qle` | `(nxi_lim,negrid,ile)` | `init_lorentz` C:2240-2285 | **yes** | v1 | **BOTH** (rebuild §5.1) |
| `ec1le, ebetaale, eqle` | `(nxi_lim,negrid,ile)` | `init_ediffuse` C:1981-1999 | **yes** | v2 | **BOTH** (rebuild §5.2) |
| `z0le, s0le, w0le` | `(nxi_lim,negrid+1,ile)` | `init_lorentz_conserve` C:691-1194 | **yes** | v1 | **BOTH** (rebuild §5.3) |
| `bz0le, bs0le, bw0le` | `(nxi_lim,negrid+1,ile)` | `init_diffuse_conserve` C:1198-1679 | **yes** | v2 | **BOTH** (rebuild §5.4) |
| `aj0le` | `(nxi_lim,negrid+1,ile)` | `init_le_bessel` C:651-662 | no | no | EXPORT (or recompute = gather of `aj0`) |
| `vperp_aj1le` | same | C:663-678 | no | no | EXPORT (= `aj1*energy_maxw(ie)*al(il)`, order `(aj1*E)*al`) |
| `vpa_aj0_le` | same | C:680-686 | no | no | EXPORT (= `vpa*aj0`) |
| `c1, betaa, ql` | `(nxi_lim, ilz)` | C:2313-2353 | yes | v1 | only if `use_le_layout=F` (see §4.8) |
| `ec1, ebetaa, eql` | `(negrid, ielo)` | C:2016-2033 | yes | v2 | only if `use_le_layout=F` |
| `z0,s0,w0,bz0,bs0,bw0` | `(-ntgrid:ntgrid,2,iglo)` | C:1172-1192, 1657-1677 | yes | yes | only if `use_le_layout=F` (g-canonical) |
| `d1le,h1le,d1,h1` | as c1le / c1 | `get_lorentz_matrix` dd,hh | no | v1 | heating diagnostic only — do not port |
| `c_rate` (dist_fn_arrays) | `(-ntgrid:ntgrid,ntheta0,naky,nspec,3)` | LD:104,143-157; C:3842-3891 | – | – | diagnostic only — do not port |
| `dtot, fdf, fdb` | `(-ntgrid:ntgrid,nlambda,5)`,`(…,nlambda)` | `init_lorentz_error` C:2646-2965 | no | no | diagnostic only (collision_error) — do not port |
| `vnmult(2)` | real(2) | C:72, C:1731-1733, 2225/2298, 1971/2006 | – | – | EXPORT the value **after** init (≥1 after the `max(1,·)` clamps) |

### 3.2 Collision-module scalars / flags (EXPORT all; most come from `collisions_knobs`)

`collision_model_switch` (int: 1 lorentz, 3 none, 5 lorentz-test, 6 full, 7 ediffuse; C:60-64),
`colls`, `has_lorentz`, `has_diffuse`, `drag`, `use_le_layout`, `conserve_moments`, `conservative`,
`resistivity`, `const_v`, `ei_coll_only`, `cfac`, `hypermult`, `hyper_colls`, `special_wfb_lorentz`,
`vpar_zero_mean`, `conserve_forbid_zero`, `force_collisions`, `adjust`, `heating`, `split_collisions`,
`lorentz_switch`, `ediff_switch`, `nxi_lim`, `vary_vnew`.

### 3.3 Arrays from other modules read by collisions

| name | module | shape | dt | decision / notes |
|---|---|---|---|---|
| `code_dt` | gs2_time | scalar | – | runtime value; drives rebuild |
| `tunits` | gs2_time | `(naky)` | no | EXPORT (1, or `2/aky` with `wstar_units`; gs2_time.f90:107-125) |
| `aparnew` | fields_arrays | `(-ntgrid:ntgrid,ntheta0,naky)` | runtime | **module variable** read directly by the drag term (C:3041, 3074) |
| `phinew, bparnew` | fields_arrays | same | runtime | passed to g_adjust (LD:71,109) |
| `kperp2` | kt_grids | `(-ntgrid:ntgrid,ntheta0,naky)` | no | EXPORT (unchanged by flow shear: exb_shear remaps `it`, not kperp2) |
| `aky` | kt_grids | `(naky)` | no | EXPORT (hybrid electron tests only) |
| `kwork_filter` | kt_grids | `(ntheta0,naky)` logical | – | always `.false.` while stepping; only set during fields_local response build (fields_local.f90:2114-2125). Port: assert all false |
| `ieqzip` | run_parameters | `(ntheta0,naky)` logical | no | EXPORT (from `eqzip_option`, run_parameters.f90:259-272) |
| `beta, fapar, zeff, fphi, fbpar` | run_parameters | scalars | no | EXPORT |
| `bmag` | theta_grid | `(-ntgrid:ntgrid)` | no | EXPORT |
| `al` | le_grids | `(nlambda)` | no | EXPORT |
| `wl` | le_grids | `(-ntgrid:ntgrid,nlambda)` | no | EXPORT (rebuild integrals) |
| `wxi, xi` | le_grids | `(2*nlambda,-ntgrid:ntgrid)` | no | EXPORT |
| `ixi_to_il, ixi_to_isgn` | le_grids | `(2*nlambda,-ntgrid:ntgrid)` int | no | EXPORT |
| `jend` | le_grids | `(-ntgrid:ntgrid)` int | no | EXPORT |
| `forbid` | le_grids | `(-ntgrid:ntgrid,nlambda)` logical | no | EXPORT |
| `w` | le_grids | `(negrid,nspec)` | no | EXPORT (conservation integrals) |
| `energy_maxw, speed_maxw, w_maxw` | le_grids | `(negrid)` | no | EXPORT (Maxwellian species grid; collisions require ≥1 Maxwellian species, LG:1887-1888) |
| `ng2, nlambda, negrid, nxi, ntgrid` | le_grids / theta_grid | int | no | EXPORT |
| `aj0, aj1` | dist_fn_arrays | `(-ntgrid:ntgrid, iglo)` | no | EXPORT (aj1 = J1(x)/x; both 0 where forbid, DF:1786-1792) |
| `vpa` | dist_fn_arrays | `(-ntgrid:ntgrid,2,iglo)` | no | EXPORT |
| `vperp2` | dist_fn_arrays | `(-ntgrid:ntgrid,iglo)` | no | EXPORT (g_adjust) |
| `nonmaxw_corr` | species | `(negrid,nspec)` | no | EXPORT (g_adjust) |
| `spec(is)%{vnewk,z,dens,stm,zstm,smz,zt,nu_h,type}` | species | `(nspec)` | no | EXPORT |
| `layout` | gs2_layouts | char(5) | no | EXPORT (fixes the summation order of rebuild integrals, §5.5) |

### 3.4 Helper tables to add to the export (remove transcendental functions from the rebuild)

| name | definition (Fortran, compute in the exporter with the *same* expression) | used by |
|---|---|---|
| `e15(ie)` | `energy_maxw(ie)**1.5` | drag term (C:3075), `vns(:,:,:,3)` (C:740) |
| `capg_half(j)`, j=1..negrid-1 | `capg((speed_maxw(j)+speed_maxw(j+1))*0.5)` with `capg(x) = 2*exp(-x**2)*(erf(x) - 2*x*exp(-x**2)/sqrt_pi)/(x*sqrt_pi)` (C:2165-2176) | ediffuse matrix: `capgr(ie) = capg_half(ie)` (ie ≤ negrid-1), `capgl(ie) = capg_half(ie-1)` (ie ≥ 2) |

`xer` at ie and `xel` at ie+1 are bit-identical (`(xe1+xe2)*0.5` vs `(xe0+xe1)*0.5` with the same two
operands; IEEE `+` is commutative), so a single table suffices. Copy the `capg_kernel/capg` code
verbatim (it is an internal procedure, not callable from outside).

---------------------------------------------------------------------------------------------------

## 4. Exact runtime computation (default `use_le_layout = .true.`)

Notation: Fortran index ranges are inclusive; `for a in x..y` ascending, `for a in y..x step -1`
descending. `E(ie) = energy_maxw(ie)`, `V(ie) = speed_maxw(ie)`, `W(ie,is) = w(ie,is)`,
`pw = pitch_weights`. "real×complex" and "complex/real" are **component-wise** (see §8.6).
Left-to-right association is written explicitly with parentheses where it matters.

Per-`ile` constants (depend on `ig` only unless noted):
```
cur_jend = jend(ig)
JE       = max(2*cur_jend, 2*ng2+1)        # = 2*cur_jend on a trapped grid, 2*ng2+1 otherwise
is_wfb   = (cur_jend == ng2+1)             # il_is_wfb(cur_jend); false when cur_jend==0
skip(is,ik)   = (abs(vnew(ik,1,is)) < 2*epsilon(0.0)) .and. .not. force_collisions
active(it,ik) = .not. kwork_filter(it,ik) .and. .not. ieqzip(it,ik)
```
`epsilon(0.0)` is the **double** epsilon (2.220446049250313e-16) because GS2 is built with
`-fdefault-real-8` (Makefiles/Makefile.gnu-gfortran:18).

### 4.1 `g_adjust` (DFA:117-175) — h ↔ g

```
from_g_gs2: facphi = fphi,  facbpar = fbpar        to_g_gs2: facphi = -fphi, facbpar = -fbpar
bpar_factor = 2.0*facbpar if facbpar /= 0 else 0.0
for iglo (it,ik,il,ie,is):
   phi_factor = (spec(is)%zt*nonmaxw_corr(ie,is))*facphi  if facphi /= 0 else 0.0
   for ig in -ntgrid..ntgrid:
      if trapped_grid and il > jend(ig): cycle            # forbidden: untouched
      adj = ((bpar_factor*vperp2(ig,iglo))*aj1(ig,iglo))*bpar(ig,it,ik) + (phi_factor*phi(ig,it,ik))*aj0(ig,iglo)
      g(ig,1,iglo) += adj ;  g(ig,2,iglo) += adj
```
`adj(to)` is the exact IEEE negation of `adj(from)` (all factors flip sign exactly), so implement the
reverse as `h + (-adj)`; it equals `h - adj` bitwise. `(phi_factor*phi)*aj0` is (real×complex)×real.
`with_phi=F` ⇒ `phi_factor = 0.0` literally multiplies (still evaluated); `0*finite = ±0`, harmless.
**Does not apply when `adjust = .false.`** (LD:70, 108).

### 4.2 Gather / scatter

```
gle[ixi, ie, (ig,it,ik,is)] = 0  for all ixi in 1..NX, ie in 1..negrid+1
for (ig, il, isgn) allowed (not forbid(ig,il)), ie in 1..negrid:
    gle[ixi(ig,il,isgn), ie] = h[ig, isgn, it, ik, il, ie, is]
... solfp1 ...
for (ig, il, isgn) allowed: h[ig, isgn, ...] = gle[ixi(ig,il,isgn), ie]      # forbidden h untouched
```
Row `ie = negrid+1` and rows `ixi > 2*cur_jend` (trapped) are padding/forbidden: 0 on entry.

### 4.3 `solfp1_le_layout` stage order (C:3035-3090)

```
if has_diffuse:  solfp_ediffuse_le(gle);  if conserve_moments: conserve_diffuse_le(gle)
if has_lorentz:  if drag: drag_le(gle)
                 solfp_lorentz_le(gle);   if conserve_moments: conserve_lorentz_le(gle)
```

### 4.4 Energy diffusion solve `solfp_ediffuse_le_layout` (C:4037-4076)

```
for each ile:
   if skip(is,ik): cycle ; if not active(it,ik): cycle
   max_xi = max(2*jend(ig), 2*ng2)
   for ixi in 1..max_xi:
      delta(1) = gle(ixi,1)
      for ie in 1..negrid-1:
         delta(ie+1) = gle(ixi,ie+1) - eqle(ixi,ie+1)*delta(ie)
      gle(ixi,negrid+1) = 0.0
      for ie in negrid..1 step -1:
         gle(ixi,ie) = (delta(ie) - ec1le(ixi,ie)*gle(ixi,ie+1)) * ebetaale(ixi,ie)
```
Both bounce-point copies (`ixi = cur_jend` and `2*cur_jend`) are solved independently (same matrix,
since the matrix depends on `il = ixi_to_il` only). `ec1le(ixi,negrid) = 0`.

### 4.5 Energy-diffusion conservation `conserve_diffuse_le_layout` (C:3541-3628)

Three sequential passes over all `ile` with `active(it,ik)` (no `skip` test):
```
# pass 1: y1 = y0 - v0y0 * z0
v0y0 = (0,0)
for ie in 1..negrid:
   vnE = (vnmult(2)*vnew_E(ik,ie,is))*W(ie,is)
   for ixi in 1..nxi_lim:
      v0y0 = v0y0 + ((vnE*aj0le(ixi,ie))*gle(ixi,ie))*wxi(ixi,ig)
gle(:,:) = gle(:,:) - bz0le(:,:)*v0y0          # all ixi in 1..nxi_lim, all ie in 1..negrid+1
# pass 2: y2 = y1 - v1y1 * s1
v1y1 = 0
for ie: delnu = (vnmult(2)*delvnew(ik,ie,is))*W(ie,is)
        for ixi: v1y1 = v1y1 + ((vpa_aj0_le(ixi,ie)*delnu)*gle(ixi,ie))*wxi(ixi,ig)
gle = gle - bs0le*v1y1
# pass 3: x = y2 - v2y2 * w2
v2y2 = 0
for ie: delnu = (vnmult(2)*delvnew(ik,ie,is))*W(ie,is)
        for ixi: v2y2 = v2y2 + ((delnu*vperp_aj1le(ixi,ie))*gle(ixi,ie))*wxi(ixi,ig)
gle = gle - bw0le*v2y2
```
Note the different factor orders: pass 2 is `(vpa_aj0_le*delnu)`, pass 3 is `(delnu*vperp_aj1le)`.

### 4.6 Drag term (C:3056-3085), only if `drag`

```
for each ile with is_electron_species(spec(is)) (electron or hybrid) and active(it,ik):
   den(ie) = ((((-(spec(is)%z*spec(is)%dens))*beta)*spec(is)%stm)*e15(ie))
   for ie in 1..negrid:
      num = (((vnmult(1)*spec(is)%vnewk)*code_dt)*kperp2(ig,it,ik))     # real
      tmp = (num*aparnew(ig,it,ik)) / den(ie)                          # complex, component-wise
      for ixi in 1..nxi_lim: gle(ixi,ie) = gle(ixi,ie) + tmp*vpa_aj0_le(ixi,ie)
```
No `tunits` factor and no `skip` test here. (`-spec%z*spec%dens` parses as `-(z*dens)`; sign flips are exact.)

### 4.7 Lorentz solve `solfp_lorentz_le_layout` (C:3897-3980)

```
for each ile:
   if skip(is,ik): cycle ; if not active(it,ik): cycle
   nxi_scatt = JE - 1
   if is_wfb and special_wfb_lorentz: nxi_scatt = nxi_scatt - 1
   for ie in 1..negrid:
      if (JE == 2*cur_jend) and is_wfb and vpar_zero_mean:                     # trapped grid, wfb bounce point
          gle(cur_jend,ie) = (gle(cur_jend,ie) + gle(2*cur_jend,ie)) / 2.0
      gle(JE,ie) = 0.0                                                         # duplicate / spare point
      if is_wfb and special_wfb_lorentz:
          gwfb = gle(ng2+1,ie)
          gle(ng2+1:JE-2, ie) = gle(ng2+2:JE-1, ie)      # shift left (RHS evaluated first)
          gle(JE-1,ie) = 0.0
      delta(1) = gle(1,ie)
      for il in 1..nxi_scatt:
          delta(il+1) = gle(il+1,ie) - qle(il+1,ie)*delta(il)
      gle(nxi_scatt+1,ie) = delta(nxi_scatt+1)*betaale(nxi_scatt+1,ie)
      for il in nxi_scatt..1 step -1:
          gle(il,ie) = (delta(il) - c1le(il,ie)*gle(il+1,ie)) * betaale(il,ie)
      if is_wfb and special_wfb_lorentz:
          gle(ng2+2:JE-1, ie) = gle(ng2+1:JE-2, ie)      # shift right
          gle(ng2+1,ie) = gwfb
      if cur_jend /= 0: gle(2*cur_jend,ie) = gle(cur_jend,ie)          # bounce condition
```
Consequences: at every trapped bounce point the **isgn = 2 value** (ixi = cur_jend) is the one that is
scattered; the isgn = 1 input value is discarded (except at the wfb bounce point with
`vpar_zero_mean`, where the average is used), and after the solve both copies are equal.
On a trapped grid `nxi_scatt+1 = 2*cur_jend = te2+1`; `qle(te2+1)=betaale(te2+1)=0` so the last row
gives exactly 0 and the effective system size is `te2 = 2*cur_jend-1`.
**No-trapped grid**: `JE = 2*ng2+1` exceeds `nxi_lim = 2*ng2`, so `delta(2ng2+1)`, `qle(2ng2+1,·)`,
`betaale(2ng2+1,·)` are out of bounds (see §8.1). The intended semantics (and what the JAX port should
implement) is a `2*ng2` system with `gle(2ng2+1)=0`.

### 4.8 Lorentz conservation `conserve_lorentz_le_layout` (C:3256-3376)

All passes over `ile` with `active(it,ik)` (no `skip`):
```
if drag:                                   # pass 0
   v0y0 = 0
   for ie in 1..negrid: for ixi in 1..nxi_lim:
       v0y0 = v0y0 + ((vpa_aj0_le(ixi,ie)*gle(ixi,ie))*W(ie,is))*wxi(ixi,ig)
   gle = gle - z0le*v0y0
# pass 1 (conservative = .true.)
v1y1 = 0
for ie: nud = ((V(ie)*vnmult(1))*vnew_D(ik,ie,is))*W(ie,is)
        for ixi: v1y1 = v1y1 + ((((vpdiffle(ixi,ig)*nud)*aj0le(ixi,ie))*gle(ixi,ie))*wxi(ixi,ig))
#   (conservative = .false.: nud = (vnmult(1)*vnew_D(ik,ie,is))*W(ie,is);
#    v1y1 = v1y1 + (((nud*vpa_aj0_le)*gle)*wxi)  — C:3337-3344; GATED)
gle = gle - s0le*v1y1
# pass 2
v2y2 = 0
for ie: nud = (vnmult(1)*vnew_D(ik,ie,is))*W(ie,is)
        for ixi: v2y2 = v2y2 + (((nud*vperp_aj1le(ixi,ie))*gle(ixi,ie))*wxi(ixi,ig))
gle = gle - w0le*v2y2
```
Each pass uses the `gle` produced by the previous pass. Updates are over the full `(nxi_lim, negrid+1)`
block (padding entries of the vectors are 0).

### 4.9 Standard layout (`use_le_layout = .false.`) — differences only

`vspace_derivatives` LD:113-167 → `solfp1_standard_layout(g, g1=g_work, gc1, gc2, diagnostics)`
(C:2968-3032), all in g_lo with redistributes to e_lo / lz_lo:
* ediffuse (C:3986-4034): `gather(energy_map)` → per `ielo` (skip, active, `forbid(ig,il)` cycle) same
  sweep as §4.4 on `ged(1:negrid+1)` → `scatter`. Identical arithmetic per row.
* conserve_diffuse (C:3379-3538): moments via `integrate_moment_c34` (LG:1321-1422):
  `total(:,it,ik,is) += (W(ie,is)*wl(:,il))*(g(:,1,iglo)+g(:,2,iglo))` in **iglo order**, with
  integrands `gtmp = ((vnmult(2)*vnew_E)*aj0)*g`, `((vns*vpa)*aj0)*g1`, `((((vns*E)*al)*aj1)*g1)`
  (`vns = vnmult(2)*delvnew`), then `g1 = g - v0y0*bz0` etc. (only for active modes; others keep stale
  g_work but are never used).
* drag (C:3003-3027): `g += (((((vnmult(1)*vnewk)*code_dt)*vpa)*kperp2)*aparnew)*aj0 / den` (different
  association from §4.6).
* Lorentz (C:3631-3791): `gather(lambda_map)` (all points incl. forbidden) → identical per-`ilz` solve
  to §4.7 with `glz(:,ilz)` → `scatter`.
* conserve_lorentz (C:3093-3253): moments via `integrate_moment_c34` with integrands
  `(vpa*aj0)*g`, `(((vns*V)*vpdiff)*aj0)*g1` (conservative), `((((vns*E)*al)*aj1)*g1)`, `vns = vnmult(1)*vnew_D`.
Physics is identical; only summation order/association differ (verified 3.2e-16 after 40 steps on a
trapped grid, §8.1). **No-trapped grids: runtime out-of-bounds write `glz(2ng2+1,ilz)` aliases
`glz(1,ilz+1)` (zeroes the ξ≈1 point of the next row) — wrong physics (21 % error observed).**

### 4.10 `split_collisions = .true.`

`timeadv` skips collisions (DF:4083); after the second `timeadv`, `collisions_advance` (DF:4130-4167)
applies `vspace_derivatives` once, then recomputes fields from h: `g_adjust(gnew, phinew, bparnew, from)`
→ `calculate_potentials_from_nonadiabatic_dfn(gnew, phinew, aparnew, bparnew)` → `g_adjust(gnew, phinew,
bparnew, to)` (with the new fields). The response matrix then excludes collisions. Field-solve details
belong to the fields spec. `timesteps_between_collisions` is ignored (C:336-340).

### 4.11 Heating (`heating = .true.`) — evolution-neutral

Computes `c_rate` only (LD:72-110, C:3812-3894, C:3653-3705). Side effect: `g_adjust(gold=g, …)` is
applied and reversed, changing the *old* `g` by round-off; `g` is overwritten by `copy(gnew, g)` at the
start of the next step (fields_implicit.f90:337) and not otherwise read after collisions, so evolution
is unchanged. Do not port `c_rate` in phase 1.

---------------------------------------------------------------------------------------------------

## 5. Rebuild formulas for dt- / vnmult-dependent arrays

All rebuilds are deterministic functions of exported, dt-independent data plus `code_dt` and
`vnmult`. Recommended implementation: **NumPy on the host, float64, explicit loops along the
recurrence/summation axes, vectorised along independent axes**, then upload to device. Validate at the
initial dt against the exported arrays (expect bitwise equality if FMA contraction is disabled in the
reference build, else ≤ a few ulp; §8.5).

### 5.1 Lorentz matrices: `get_lorentz_matrix` (C:2392-2639) + elimination (C:2240-2285)

For each `(ig, it, ik, is)` and `ie in 1..negrid` (arrays `aa, bb, cc` of length `nxi_lim`, zeroed):
```
je = jend(ig)
if (je <= ng2) or (je == ng2+1 and special_wfb_lorentz):     # il_is_passing(je) .or. (il_is_wfb(je) .and. special)
     te = ng2;  te2 = 2*ng2;  teh = ng2                      # (je==0 ⇒ passing branch)
else:
     te = je;   te2 = 2*je-1; teh = je-1
if collision_model_switch == lorentz_test:
     vn = (vnmult(1)*abs(spec(is)%vnewk))*tunits(ik)
else:
     vhyp = vnewh(ig,it,ik,is) if hyper_colls else 0.0
     if hypermult: vn = (vnmult(1)*vnew(ik,ie,is))*(1 + vhyp)
     else:         vn = vnmult(1)*vnew(ik,ie,is) + vhyp
sq(x)  = sqrt(max(0.0, x))                                   # safe_sqrt C:569-573
s(il)  = sq(1.0 - bmag(ig)*al(il))
eeF(slb1) = ((((0.5*E(ie))*(1 + slb1**2))*kperp2(ig,it,ik))*cfac) / ((bmag(ig)*spec(is)%zstm)**2)
# --- lorentz_scheme = 'default' (C:2455-2554)
for il in 2..te-1:
   slb0 = s(il-1); slb1 = s(il); slb2 = s(il+1)
   slbl = (slb1 + slb0)*0.5 ; slbr = (slb1 + slb2)*0.5
   ee = eeF(slb1)
   cc(il) = (((2.0*vn)*code_dt)*(1 - slbr**2)) / (pw(il,ig)*(slb2 - slb1))
   aa(il) = (((2.0*vn)*code_dt)*(1 - slbl**2)) / (pw(il,ig)*(slb1 - slb0))
   bb(il) = (1.0 - (aa(il) + cc(il))) + ((ee*vn)*code_dt)
# xi = 1 boundary
slb1 = s(1); slb2 = s(2); slbr = (slb1 + slb2)*0.5 ; ee = eeF(slb1)
cc(1) = (((2.0*vn)*code_dt)*(1.0 - slbr**2)) / (pw(1,ig)*(slb2 - slb1))
aa(1) = 0.0
bb(1) = (1.0 - cc(1)) + ((ee*vn)*code_dt)
# xi = 0 boundary, il = te
slb0 = s(te-1)
if te == ng2: slb1 = s(te); slb2 = -slb1
else:         slb1 = 0.0;   slb2 = -slb0
slbl = (slb1 + slb0)*0.5 ; slbr = (slb1 + slb2)*0.5 ; ee = eeF(slb1)
deltaxi = pw(te,ig)
if te == je: deltaxi = 2*deltaxi                                   # vpar=0 point appears twice
if (not special_wfb_lorentz) and deltaxi == 0 and te == ng2+1:     # il_is_wfb(te)
     deltaxi = xi(ng2,ig) - xi(ng2+2,ig)                           # xi is indexed by ixi
cc(te) = (((2.0*vn)*code_dt)*(1 - slbr**2)) / (deltaxi*(slb2 - slb1))
aa(te) = (((2.0*vn)*code_dt)*(1 - slbl**2)) / (deltaxi*(slb1 - slb0))
bb(te) = (1.0 - (aa(te) + cc(te))) + ((ee*vn)*code_dt)
# symmetric fill (C:2630-2632), RHS read before write; indices teh..1 are untouched by the writes
aa(te+1:te2) = cc(teh:1:-1)
bb(te+1:te2) = bb(teh:1:-1)
cc(te+1:te2) = aa(teh:1:-1)
# hybrid electrons (C:2265-2266, 2366-2389), only if has_hybrid_electron_species(spec):
if spec(is)%type == hybrid and aky(ik) /= 0:
   L = ng2+2 ; U = 2*je-1-(ng2+1)
   aa(1:L) = 0; bb(1:L-1) = 1; cc(1:L-1) = 0; aa(U+1:) = 0; bb(U+1:) = 1; cc(U:) = 0
```
`x**2` is `x*x` (integer power). `1 - x` with integer 1 is `1.0 - x`. Unary minus placement never
changes the rounded result (sign-symmetric rounding).

Elimination (le layout, C:2268-2282):
```
c1le(:,ie) = cc(1:nxi_lim)
qle(1,ie) = 0.0
betaale(1,ie) = 1.0 / bb(1)
for ixi in 1..te2-1:
    qle(ixi+1,ie)     = aa(ixi+1) * betaale(ixi,ie)
    betaale(ixi+1,ie) = 1.0 / (bb(ixi+1) - qle(ixi+1,ie)*c1le(ixi,ie))
qle(te2+1:,ie) = 0.0 ; betaale(te2+1:,ie) = 0.0      # c1le beyond te2 is already 0 (cc zero-initialised)
```
(Standard layout C:2337-2351 is identical plus `c1(te2+1:)=0`.) `cc(te2) = aa(1) = 0` always.

`lorentz_scheme = 'old'` (C:2557-2626, GATED): interior `cc = (-vn*code_dt*(1-slbr**2))/((slbr-slbl)*(slb2-slb1))`,
`aa = (-vn*code_dt*(1-slbl**2))/((slbr-slbl)*(slb1-slb0))`; xi=1: `slb0=1`, `cc(1) = -vn*code_dt*(-1.0-slbr)/(slb2-slb1)`,
`bb(1) = 1.0-(aa(1)+cc(1)) + ee*vn*code_dt`; xi=0 as interior formula with the xi=0 `slb` choice; no `deltaxi`.

Heating-only coefficients `dd, hh` (C:2476-2481, 2498-2503, 2551-2554) need `vnc = vnmult(1)*vnew`,
`vnh = vhyp*vnc` (hypermult) or `vhyp`; not ported.

### 5.2 Energy-diffusion matrices: `get_ediffuse_matrix` (C:2041-2182) + elimination (C:1981-1999)

For each `(ig, it, ik, is)` and each `ixi in 1..nxi_lim` with `il = ixi_to_il(ixi,ig)` and
**not** `forbid(ig,il)` (forbidden rows stay 0; C:1988):
```
xe = speed_maxw ; WM = w_maxw
vn   = (vnmult(2)*spec(is)%vnewk)*tunits(ik)
slb1 = sq(1.0 - bmag(ig)*al(il))
eeE(ie) = ((((0.125*(1 - slb1**2))*vnew_s(ik,ie,is))*kperp2(ig,it,ik))*cfac) / ((bmag(ig)*spec(is)%zstm)**2)
# ediff_scheme = 'default'
for ie in 2..negrid-1:
    capgr = capg_half(ie) ; capgl = capg_half(ie-1)
    cc(ie) = -((((0.25*vn)*code_dt)*capgr) / (WM(ie)*(xe(ie+1) - xe(ie))))
    aa(ie) = -((((0.25*vn)*code_dt)*capgl) / (WM(ie)*(xe(ie) - xe(ie-1))))
    bb(ie) = (1.0 - (aa(ie) + cc(ie))) + (eeE(ie)*code_dt)
cc(1) = -((((0.25*vn)*code_dt)*capg_half(1)) / (WM(1)*(xe(2) - xe(1))))
aa(1) = 0.0
bb(1) = (1.0 - cc(1)) + (eeE(1)*code_dt)
cc(negrid) = 0.0
aa(negrid) = -((((0.25*vn)*code_dt)*capg_half(negrid-1)) / (WM(negrid)*(xe(negrid) - xe(negrid-1))))
bb(negrid) = (1.0 - aa(negrid)) + (eeE(negrid)*code_dt)
# elimination
ec1le(ixi,:) = cc
ebetaale(ixi,1) = 1.0 / bb(1)            # eqle(ixi,1) stays 0
for ie in 1..negrid-1:
    eqle(ixi,ie+1)     = aa(ie+1) * ebetaale(ixi,ie)
    ebetaale(ixi,ie+1) = 1.0 / (bb(ie+1) - eqle(ixi,ie+1)*ec1le(ixi,ie))
```
Notes: the classical (FLR) term `eeE*code_dt` is **not** scaled by `vnmult(2)` (unlike Lorentz,
which uses `ee*vn*code_dt`). `vnew_E` is **not** used by the matrix (only by conservation).
`ediff_scheme = 'old'` (C:2112-2161) uses `capg_old(xe1,xeB) = (0.5*exp(xe1**2-xeB**2)/xe1**2)*capg_kernel(xeB)/xeB`
and denominators `(xer-xel)*(…)`, `xer*(…)`, `(1.0-xel)*(…)` — GATED.
Standard layout (C:2016-2033): same per `ielo` (skips `forbid`).

### 5.3 Lorentz conservation vectors `init_lorentz_conserve` (C:691-1194)

Work in g-canonical real arrays (the Fortran uses complex arrays with zero imaginary part; real
arithmetic gives identical real parts). Define:
* `IM(f)` = `integrate_moment_r34(f)` with the exact iglo summation order of §5.5 → shape `(ig,it,ik,is)`.
* `SOLVE_L(v)`: `ctmp = gather_le(v)` (allowed points; others hold 0), apply the **runtime** Lorentz
  solve §4.7 (same skips, same wfb/duplicate handling), `v = scatter_le(ctmp)` (allowed points only).
  (Standard layout: `solfp_lorentz_standard_layout(v)`.)
* `ZH(v)`: `zero_out_passing_hybrid_electrons` (C:1684-1704): `v(:,:,iglo)=0` if
  `type==hybrid ∧ aky(ik)≠0 ∧ il ≤ ng2+1`. Identity if no hybrid species.
* `vns1(ik,ie,is) = vnmult(1)*vnew_D(ik,ie,is)`;
  `vns3(ik,ie,is) = ((vnmult(1)*spec(is)%vnewk)*tunits(ik))/e15(ie)` for electron-like species if `drag`, else 0.

Below `D(ig,it,ik,is)` denotes the moment just computed; broadcasting over `isgn, il, ie` is implicit;
conservative = .true. branches shown (non-conservative alternatives in brackets, GATED).
```
if drag:
   z0 = -(((((2.0*code_dt)*vns3)*vpdiff(ig,isgn,il))*V(ie))*aj0)        [ -((((2.0*code_dt)*vns3)*vpa)*aj0) ]
   z0 = ZH(z0); z0 = SOLVE_L(z0)
   D  = IM((vpa*aj0)*z0)                                    # v0z0
   z0 = z0 / (1.0 + D)
else:
   z0 = 0
du    = IM(((vns1*vpa)*vpdiff(ig,isgn,il))*V(ie))           [ IM(vns1*(vpa**2)) ]
duinv = where(abs(du) > epsilon(0.0), 1.0/du, du)           # C:856-859 (keeps du, i.e. 0, otherwise)
s0 = -(((((vns1*vpdiff)*V(ie))*aj0)*code_dt)*duinv)         [ -((((vns1*vpa)*aj0)*code_dt)*duinv) ]
s0 = ZH(s0); s0 = SOLVE_L(s0)
D  = IM((vpa*aj0)*s0)                                       # v0s0
s0 = s0 - D*z0
D  = IM(((((vns1*V(ie))*vpdiff)*aj0)*s0))                   [ IM((((vns1*vpa)*aj0)*s0)) ]   # v1s0
s0 = s0 / (1.0 + D)
w0 = -((((((((vns1*E(ie))*al(il))*aj1)*code_dt)*(smz**2))*kperp2(ig,it,ik))*duinv)/bmag(ig))
w0 = where(forbid(ig,il) and conserve_forbid_zero, 0.0, w0) # (aj1=0 there anyway ⇒ only sign of 0)
w0 = ZH(w0); w0 = SOLVE_L(w0)
D  = IM((vpa*aj0)*w0)                                       # v0w0
w0 = w0 - z0*D
D  = IM(((((vns1*V(ie))*vpdiff)*aj0)*w0))                   [ IM((((vns1*vpa)*aj0)*w0)) ]   # v1w1
w0 = w0 - s0*D
D  = IM(((((vns1*E(ie))*al(il))*aj1)*w0))                   # v2w2
w0 = w0 / (1.0 + D)
# le layout: s0le, w0le, z0le = gather_le(s0), gather_le(w0), gather_le(z0)  (forbidden/padding = 0)
```
`smz**2 = smz*smz`. In Fortran `duinv` and every `D` are full `(ig,it,ik,is)` arrays zero-initialised
(C:729-730); `kwork_filter` is false here. The sequence order of these steps is significant (each step
uses the latest vector).

### 5.4 Energy-diffusion conservation vectors `init_diffuse_conserve` (C:1198-1679)

`SOLVE_E(v)` = gather_le → runtime ediffuse solve §4.4 → scatter_le (standard: `solfp_ediffuse_standard_layout`).
`vns1(ik,ie,is) = vnmult(2)*delvnew(ik,ie,is)`.
```
du    = IM(((E(ie)*vnmult(2))*vnew_E(ik,ie,is)))           # integrand independent of ig,isgn,il except via wl
duinv = where(abs(du) > epsilon(0.0), 1/du, du)
bz0 = -((((code_dt*vnmult(2))*vnew_E)*aj0)*duinv) ; bz0 = 0 where forbid ∧ conserve_forbid_zero
bz0 = ZH(bz0); bz0 = SOLVE_E(bz0)
D   = IM(((vnmult(2)*vnew_E)*aj0)*bz0)                    # v0z0
bz0 = bz0 / (1.0 + D)
du    = IM(vns1*(vpa**2))                                  # second, different duinv
duinv = where(abs(du) > epsilon(0.0), 1/du, du)
bs0 = -((((vns1*vpa)*aj0)*code_dt)*duinv)
bs0 = ZH(bs0); bs0 = SOLVE_E(bs0)
D   = IM(((vnmult(2)*vnew_E)*aj0)*bs0)                    # v0s0
bs0 = bs0 - D*bz0
D   = IM(((vns1*vpa)*aj0)*bs0)                            # v1s0
bs0 = bs0 / (1.0 + D)
bw0 = -((((((((vns1*E(ie))*al(il))*aj1)*code_dt)*(smz**2))*kperp2)*duinv)/bmag) ; 0 where forbid ∧ cfz
bw0 = ZH(bw0); bw0 = SOLVE_E(bw0)
D   = IM(((vnmult(2)*vnew_E)*aj0)*bw0)                    # v0w0
bw0 = bw0 - bz0*D
D   = IM(((vns1*vpa)*aj0)*bw0)                            # v1w1
bw0 = bw0 - bs0*D
D   = IM((((vns1*E(ie))*al(il))*aj1)*bw0)                 # v2w2
bw0 = bw0 / (1.0 + D)
# le: bs0le, bw0le, bz0le = gather_le(...)
```
Integrand for the first `du` is `gtmp = energy(ie)*vnmult(2)*vnew_E(ik,ie,is)` broadcast over `ig, isgn`
(C:1251).

### 5.5 `integrate_moment_r34` summation order (LG:1427-1523), single process, 1 OpenMP thread

```
total(:,:,:,:) = 0
for iglo in g_lo%llim_proc..g_lo%ulim_proc (ascending):      # iglo order = layout string, first letter fastest
    (it,ik,il,ie,is) = idx(iglo)
    total(:,it,ik,is) = total(:,it,ik,is) + (W(ie,is)*wl(:,il))*(g(:,1,iglo) + g(:,2,iglo))
```
For a fixed `(it,ik,is)` the visiting order of `(il, ie)` is: **ie outer, il inner** for
`lxyes` (default), `lexys`, `yxles`, `lyxes`, `xyles`; **il outer, ie inner** for `yxels` (GL:917-930).
OpenMP `REDUCTION(+)` (LG:1464-1480) and MPI `sum_allreduce_sub` change the order — reference data must
be produced with 1 MPI rank and `OMP_NUM_THREADS=1`.

### 5.6 dt-independent precomputes (EXPORT; formulas for reference/verification)

* `init_le_bessel` (C:634-687): `aj0le(ixi,ie) = gather(aj0)`, `vperp_aj1le(ixi,ie) = (gather(aj1)*E(ie))*al(ixi_to_il(ixi,ig))` for `ie ≤ negrid`,
  `vpa_aj0_le = gather(vpa)*aj0le`. Forbidden/padding entries 0.
* `init_vpdiff` (C:4079-4155), conservative only:
```
for ig: je = jend(ig); te = ng2 if (je <= ng2 or (je==ng2+1 and special_wfb_lorentz)) else je
   for il in 2..te-1: slbl=(s(il)+s(il-1))*0.5; slbr=(s(il)+s(il+1))*0.5; vpdiff(ig,1,il) = (slbl**2 - slbr**2)/pw(il,ig)
   slbr = 0.5*(s(1)+s(2)); vpdiff(ig,1,1) = (1.0 - slbr**2)/pw(1,ig)
   slb0 = s(te-1); (slb1,slb2) = (s(te),-s(te)) if te==ng2 else (0.0,-slb0)
   slbl=(slb1+slb0)*0.5; slbr=(slb1+slb2)*0.5
   vpdiff(ig,1,te) = 0.0 if abs(pw(te,ig)) <= epsilon(0.0) else (slbl**2 - slbr**2)/pw(te,ig)
   vpdiff(ig,2,:) = -vpdiff(ig,1,:)                     # other il stay 0
vpdiffle(ixi,ig) = vpdiff(ig, ixi_to_isgn(ixi,ig), ixi_to_il(ixi,ig))
```
  (Note: no `deltaxi` doubling here, unlike the matrix.)
* `pitch_weights` (C:2199-2210): `setup_passing_lambda_grids(al, lw)` (LG:2160-2236; also recomputes
  `al(1:ng2)` identically) + `setup_trapped_lambda_grids_old_finite_difference(al, lw)` (LG:2617-2641),
  `pitch_weights = transpose(lw)`.
* `init_vnew` (C:1707-1876) — EXPORT (erf, exp, `**1.5`). `local_energy = 1 if const_v else E`;
  `hee = exp(-le)/sqrt(pi*le) + (1-0.5/le)*erf(sqrt(le))`; `hsg = hsg_func(sqrt(le))`
  (`0.5*erf(v)/v**2 - exp(-v**2)/(sqrt_pi*v)`, 0 if `|v|≤eps`, C:1879-1888);
  `vnew_s = vnewk/sqrt(le)*hsg*4*tunits`; electrons: `vnew = vnewk/le**1.5*(zeff+hee)*0.5*tunits`,
  `vnew_D = vnewk/le**1.5*hee*tunits`; ions: `vnew = vnewk/le**1.5*hee*0.5*tunits`, `vnew_D = 2*vnew`;
  `ei_coll_only`: `vnew_s=vnew_D=vnew_E=delvnew=0`, electrons `vnew = vnewk/le**1.5*zeff*0.5*tunits`, others 0;
  `conservative`: `vnew_E`, `delvnew` from flux differences at `v_{ie±1/2}` (C:1801-1836, uses `E` and
  `w_maxw`, **not** `local_energy`; **overrides the ei_coll_only zeros**); hyper:
  `vnewh(:,it,ik,is) = (nu_h*kperp2(:,it,ik)**2)/k4max`, `k4max = maxval(kperp2)**2`.

### 5.7 Dependency summary (what to rebuild when)

| trigger | rebuild |
|---|---|
| `code_dt` changes | §5.1 (c1le, betaale, qle), §5.2 (ec1le, ebetaale, eqle), §5.3 (z0le, s0le, w0le), §5.4 (bz0le, bs0le, bw0le); plus the field response matrix (fields spec), which contains the collision operator unless `split_collisions` |
| `vnmult(1)` changes | §5.1, §5.3 (runtime §4.6/§4.8 read vnmult directly) |
| `vnmult(2)` changes | §5.2, §5.4 (runtime §4.5 reads vnmult directly) |
| flow shear (`g_exb`) | nothing (kperp2 grid fixed) |
| restart | nothing beyond exporting the post-init `vnmult` |

---------------------------------------------------------------------------------------------------

## 6. Namelist / config options

`collisions_knobs` (C:144-375) unless stated.

| option | namelist | default | effect | decision |
|---|---|---|---|---|
| `collision_model` | collisions_knobs | `'default'` | `'default'`→full (Lorentz+ediffuse), `'lorentz'`, `'ediffuse'`, `'lorentz-test'`, `'none'`/`'collisionless'` (C:475-481, 541-550). Auto-disabled if all `vnew(:,1,:)≈0` (C:605-609) | SUPPORT all |
| `use_le_layout` | collisions_knobs | `.true.` | data path; physics identical | SUPPORT T; F → run the le algorithm (≤1e-15 deviation, verified 3e-16) except GATE on no-trapped grids |
| `conserve_moments` | collisions_knobs | `.true.` | restoring terms (§4.5, §4.8) | SUPPORT (GATE with no-trapped grid + Lorentz) |
| `conservative` | collisions_knobs | `.true.` | vpdiff-based conservative forms; F leaves vnew_E/delvnew endpoints uninitialised | SUPPORT T; GATE F |
| `resistivity` | collisions_knobs | `.true.` | enables `drag` (with beta>eps, nspec>1, fapar>0) | SUPPORT |
| `const_v` | collisions_knobs | `.false.` | energy-independent ν (init_vnew only) | SUPPORT (export) |
| `ei_coll_only` | collisions_knobs | `.false.` | only e-i collisions (init_vnew only; note conservative overrides vnew_E/delvnew) | SUPPORT (export) |
| `cfac` | collisions_knobs | `1.0` | FLR/classical diffusion factor in `ee` terms | SUPPORT |
| `hypermult` | collisions_knobs | `.false.` | `vn = vnmult*vnew*(1+vhyp)` vs `+vhyp` | SUPPORT |
| `nu_h` | species_parameters_N | `0.0` | hyper-collisions (`hyper_colls`), Lorentz only | SUPPORT |
| `vnewk` | species_parameters_N | `0.0` | collision frequency | SUPPORT |
| `type` (`'electron'`,`'hybrid'`) | species_parameters_N | `'ion'` | drag applies to electron-like; hybrid modifies matrices | SUPPORT electron/ion; GATE hybrid |
| `special_wfb_lorentz` | collisions_knobs | `.false.` | remove wfb from scattering at its bounce point | SUPPORT (spec'd §4.7, §5.1) |
| `vpar_zero_mean` | collisions_knobs | `.true.` | average the two vpar=0 wfb copies before the solve | SUPPORT |
| `conserve_forbid_zero` | collisions_knobs | `.true.` | zero w0/bz0/bw0 at forbidden points (numerically no effect: aj0=aj1=0 there) | SUPPORT |
| `lorentz_scheme` | collisions_knobs | `'default'` | `'old'` non-conservative matrix | SUPPORT default; GATE `'old'` |
| `ediff_scheme` | collisions_knobs | `'default'` | `'old'` non-conservative matrix | SUPPORT default; GATE `'old'` |
| `adjust` | collisions_knobs | `.true.` | apply g↔h transform around collisions | SUPPORT |
| `force_collisions` | collisions_knobs | `.false.` | disable the vnew≈0 skips/auto-disable | SUPPORT |
| `split_collisions` | collisions_knobs | `.false.` | collisions after 2nd timeadv + field re-solve; excluded from response | SUPPORT (needs fields spec) |
| `timesteps_between_collisions` | collisions_knobs | `1` | ignored by the code | SUPPORT (ignored) |
| `heating` | collisions_knobs | `.false.` | c_rate diagnostic only (diagnostics may turn it off) | SUPPORT (no-op; c_rate not ported) |
| `test` | collisions_knobs | `.false.` | redistribute self-checks | SUPPORT (ignored) |
| `vary_vnew` | collisions_knobs | `.false.` | adaptive vnmult via get_verr; stale response matrix | GATE |
| `ncheck` | collisions_knobs | `100` | vary_vnew cadence | SUPPORT (no effect unless vary_vnew) |
| `vnfac`, `vnslow` | collisions_knobs | `1.1`, `0.9` | vary_vnew factors | SUPPORT (no effect unless vary_vnew) |
| `etol`, `etola`, `ewindow`, `ewindowa` | collisions_knobs | `2e-2,2e-2,1e-2,1e-2` | vary_vnew thresholds | SUPPORT (no effect unless vary_vnew) |
| `beta` | knobs | `0.0` (run_parameters.f90:55) | drag denominator and enable (`beta > eps`) | SUPPORT |
| `fapar` | knobs | `0.0` (run_parameters.f90:80) | drag enable (`fapar > 0`) | SUPPORT |
| `zeff` | knobs | `1.0` | electron vnew | SUPPORT (export) |
| `delt` / `delt_option` | knobs | – | `code_dt`; `'check_restart'` restores dt and vnmult from restart | SUPPORT (rebuild §5 on dt change) |
| `eqzip_option` | knobs | `'none'` | `ieqzip` masks skip collisions for frozen modes | SUPPORT |
| `wstar_units` | knobs | `.false.` | `tunits = 2/aky`; note vnewh and drag lack tunits (literal port) | SUPPORT (literal) |
| `layout` | layouts_knobs | `'lxyes'` | iglo order ⇒ summation order of rebuild integrals §5.5 | SUPPORT (all 6 strings) |
| `trapped_particles`, `eps`(→ no trapped grid) | le_grids_knobs / theta_grid | `.true.` | no-trapped grid hits the OOB bug with Lorentz | SUPPORT trapped; GATE no-trapped ∧ Lorentz ∧ (conserve_moments ∨ ¬use_le_layout) |
| `radau_gauss_grid`, `new_trap_int`, `ngauss/npassing`, `negrid/nesub/nesuper`, `vcut`, `genquad`, `bouncefuzz` | le_grids_knobs | – | change grids/weights only | SUPPORT (all via export) |
| LOWFLOW build | compile flag | off | `get_and_set_gvpa` in le path (LD:85-87) | GATE |

---------------------------------------------------------------------------------------------------

## 7. GPU mapping notes

* **Independent (batch / vmap) axes**: `(is, ik, it, ig)` for every stage. Represent `gle` as an
  array `(nspec, naky, ntheta0, nθ, negrid+1, NX)` (or with ξ leading for scans).
* **Energy solve**: forward/backward recurrences over `ie` → two `lax.scan`s of length `negrid-1`,
  batched over `(is, ik, it, ig, ixi)`. Rows `ixi > max_xi(ig)` may be solved unmasked if the matrices are
  zero there (forbidden rows: `ebetaale = 0` ⇒ result 0 = input 0); the padding `ie = negrid+1` must be
  set to 0.
* **Lorentz solve**: recurrences over `ixi` → two `lax.scan`s of length `NX-1`, batched over
  `(is, ik, it, ig, ie)`. The variable system size `nxi_scatt(ig)+1` can be handled **without masks**
  provided (i) `qle, betaale, c1le` are zero beyond `te2` (true by construction), (ii) entries beyond
  the system are exactly 0 on entry (forbidden ξ and `gle(JE)` zeroed), and (iii) the
  `special_wfb_lorentz` shift is applied before/after. With these, the unmasked sweep reproduces the
  Fortran values for all valid entries (padded entries come out ±0). Keep explicit pre/post steps:
  wfb averaging (mask `is_wfb(ig)` on a trapped grid), `gle(JE)=0`, special shift, bounce copy
  `gle(2*cur_jend) = gle(cur_jend)` (gather/scatter with per-ig index).
* **Skip masks** (must be explicit `where(mask, new, old)` — matrices are not identity when skipped,
  e.g. hyper-collisions with vnewk = 0): `skip(is,ik)` for both solves; `active(it,ik)` (ieqzip) for all
  stages; electron mask for drag; `drag`, `has_lorentz`, `has_diffuse`, `conserve_moments` are static.
* **Conservation moments**: per `(is,ik,it,ig)` sums over `(ie, ixi)` — `jnp.sum`/`einsum` (tree
  reduction) changes round-off only at ~1e-16 relative to the term magnitudes; acceptable for the 1e-12
  target. For bitwise debugging provide a strict mode on CPU that accumulates sequentially
  (`ie` outer, `ixi` inner, starting from complex 0) — `lax.fori_loop` over 2·nlambda·negrid steps.
* **Gather/scatter g ↔ le**: static integer index maps (from `ixi_to_il/isgn`, `jend`, `forbid`) built
  once; scatter only allowed points; forbidden h passes through unchanged.
* **No linked chains / no communication**: collisions are local in θ and in (kx, ky). With flow shear
  nothing changes (index remap happens elsewhere).
* **Rebuild on dt change**: NumPy host code (rare event); reuse the device solve kernels for `SOLVE_L`,
  `SOLVE_E` inside §5.3/§5.4 only if they are bit-identical to the host version; simplest is a host
  implementation of the same sweeps.

---------------------------------------------------------------------------------------------------

## 8. Pitfalls

### 8.1 No-trapped-particle grids: out-of-bounds indexing (confirmed)

`nxi_lim = 2*ng2` but `solfp_lorentz_*` uses `JE = 2*ng2+1` (C:3730, 3917) on arrays sized `nxi_lim`:
* le path, runtime: `gle` has `nxi+1 = 2ng2+1` rows (OK), but `delta(nxi_lim)` (stack, C:3807),
  `qle/betaale(2ng2+1,ie,ile)` (heap reads aliasing `(1,ie+1,ile)`, past the end for the last `ile`)
  are out of bounds; result equals the intended 2ng2 system unless the stray write/reads hit something.
* le path, init (`init_lorentz_conserve`): `ctmp(nxi_lim,…)` (C:714) ⇒ `gle(2ng2+1,ie,ile)=0` zeroes
  `ctmp(1,ie+1,ile)` (the ξ≈1 input of the next energy) before it is solved ⇒ wrong z0/s0/w0, and
  observed heap corruption.
* standard path, runtime and init: `glz(nxi_lim,ilz)` (C:3648) ⇒ `glz(2ng2+1,ilz)=0` zeroes
  `glz(1,ilz+1)` every step (and writes one element past the array for the last `ilz`).
* conservation updates `gle(:,:,ile) - z0le(:,:,ile)*v` use non-conforming shapes `(2ng2+1)` vs `(2ng2)`
  in the no-trapped le path (harmless: only the spare row is affected).
Evidence (bin/gs2, s-alpha, ntheta=12, 1 ion species, vnewk=0.1, aky=0.4, delt=0.05, 40 steps,
`OMP_NUM_THREADS=1`):

| grid | use_le_layout | conserve_moments | result |
|---|---|---|---|
| eps=0.18 | T vs F | T | max rel diff φ = 3.2e-16 |
| eps=0.18 | T vs F | F | identical (0.0) |
| eps=0.0 | T vs F | F | max rel diff φ = 0.21, phi2 0.38 |
| eps=0.0 | T | T | crash: SIGTRAP in `libsystem_malloc` (heap corruption detected later in netCDF init) |

Port action: GATE as stated; recommend an upstream fix (size the ξ dimension of `gle/ctmp/glz/qle/
betaale/c1le/delta` as `max(nxi_lim, 2*ng2+1)` or use `JE = 2*ng2` when `cur_jend == 0`).

### 8.2 Bounce points and wfb
* isgn=1 at the bounce point is overwritten by the isgn=2 post-collision value (§4.7) on every trapped ig.
* `vpar_zero_mean` averaging only at θ where `jend = ng2+1` (wfb bounce, i.e. θ = ±π for standard geometry).
* `deltaxi` doubling at `te == je` in the matrix but **not** in `vpdiff`; wfb fix `deltaxi = xi(ng2)-xi(ng2+2)`
  triggers only when `pitch_weights(te,ig) == 0` exactly (e.g. `radau_gauss_grid = F`).
* `xi` is indexed by `ixi`, not `il` (`xi(ng2+2,ig)` is `-xi(ng2,ig)` at the wfb bounce point).
* No ttp-specific handling exists in collisions.
### 8.3 Forbidden region
le path never touches forbidden points (not gathered, not scattered). Standard path passes them through
the solves unchanged but **integrates** them (`wl = 0` there, but `0*NaN = NaN` if g holds NaN there).
`g_adjust` skips `il > jend(ig)` (uses `jend`, not `forbid`; identical except the exact `-bouncefuzz` tie).
### 8.4 Summation order / parallelism
Rebuild integrals follow g_lo iglo order (layout dependent, §5.5); OpenMP reductions and MPI
allreduce reorder them — generate references with 1 rank, 1 thread. Runtime le conservation sums are
per-`ile` sequential (thread-independent). `integrate_moment_lec` (heating) uses `SCHEDULE(guided)`
with a reduction — diagnostic only.
### 8.5 FMA contraction
The reference build is gfortran `-O3` on arm64 (Makefiles/Makefile.gnu-gfortran:24; `uname -m` = arm64),
where GCC contracts `a - b*c` into fused multiply-adds by default (`-ffp-contract=fast`). Expressions like
`delta(il+1) = gle(il+1) - qle*delta(il)` and `1/(bb - q*c)` are FMA candidates. XLA also fuses.
Bitwise equality therefore needs `-ffp-contract=off` in a validation build of GS2 and an XLA
configuration that does not contract (open question). 1e-12 agreement does not need it.
### 8.6 Mixed real/complex arithmetic
gfortran lowers `real*complex` to `(r*a, r*b)` and `complex/real` to `(a/r, b/r)`. JAX promotes the real
to complex: multiplication gives the same values (except signed zeros), but complex/complex division in
XLA may use a scaled algorithm — implement `complex/real` component-wise explicitly (drag `tmp`,
`z0/(1+D)` analogues in the rebuild are real-only).
### 8.7 vnmult
* Clamped to ≥ 1 at first allocation in `init_lorentz`/`init_ediffuse` (C:2225, 2298, 1971, 2006).
* `set_vnmult` is called after collision init in restart/override paths
  (gs2_init.fpp:322-329): runtime conservation uses the *new* vnmult while matrices were built with
  `vnm_init`. Consistent for dt resets (C:4169-4175 hack); possibly inconsistent for restarts from runs
  that used `vary_vnew`. Export `vnmult` after full init and check it equals the build value.
### 8.8 Units / literal quirks to replicate
* With `wstar_units`, `vnew*` include `tunits` but `vnewh` and the drag term do not.
* Ediffuse `ee*code_dt` is not multiplied by `vnmult(2)`; Lorentz `ee*vn*code_dt` is.
* `vns3` (z0) contains `tunits`, the runtime drag term does not.
* `ei_coll_only` + `conservative`: `vnew_E`/`delvnew` are recomputed non-zero for all species.
* `lorentz-test` uses `abs(vnewk)*tunits` and no hyper term, but the skip test still uses `vnew`.
### 8.9 Two applications per step
Collisions run in both `timeadv` calls (first with `phinew == phi`); the first result is discarded
(`gnew` is rebuilt from `g` in the second `invert_rhs`). Both must be applied to be consistent with the
response matrix (which includes collisions). The drag term reads the **module** `aparnew`.
### 8.10 `epsilon(0.0)` is double precision; `abs(du) > eps` test keeps `du` (= 0) otherwise.
### 8.11 No half steps: collisions always use the full `code_dt` (backward Euler), also with Strang-split nonlinearity.

---------------------------------------------------------------------------------------------------

## 9. Open questions

1. No-trapped grids with Lorentz collisions: should the port implement the *intended* semantics
   (2ng2 system, uncorrupted conservation vectors) — which cannot be validated against GS2 — or refuse
   (current recommendation for conserve_moments=T / use_le_layout=F)? Will upstream GS2 be patched
   first so that a reference exists?
2. The exporter needs access to PRIVATE module arrays (`c1le, betaale, qle, ec1le, …, s0le, …, vnew*,
   vnewh, pitch_weights, vpdiff(le), aj0le, vperp_aj1le, vpa_aj0_le`, flags `drag, has_lorentz,
   has_diffuse, conservative, …`). Is a small accessor (e.g. `public :: collisions_export` in
   collisions.fpp) acceptable in the export phase?
3. Can the validation build of GS2 use `-ffp-contract=off`, and is there an XLA setting to disable FMA
   contraction on the target GPU, if bitwise checks are wanted?
4. `capg_half` and `e15` must come from the same compiler/libm as GS2: is it acceptable to compute them
   in the exporter with a verbatim copy of `capg`, or should collisions.fpp expose them?
5. Should `use_le_layout = .false.` be implemented with its own exact summation order (integrate_moment
   based), or is the le algorithm (≤1e-15 deviation) acceptable?
6. `vary_vnew`: is adaptive collisionality needed? If so, the response matrix must also be rebuilt on
   every vnmult change (GS2 does not do this — replicate the inconsistency or fix it?).
7. Restarts from runs that used `vary_vnew` (vnmult ≠ 1 in the restart file): GS2 builds matrices with
   `vnm_init` and later calls `set_vnmult`; confirm which vnmult the matrices actually contain in that case.
8. Hybrid electrons (`type = 'hybrid'`): needed? The matrix modification `set_hzero_lorentz_collisions_matrix`
   has asymmetric index ranges (`bb(:L-1)` vs `aa(:L)`, `cc(U:)` vs `aa(U+1:)`); confirm intent if supported.
9. `pitch_weights` uses old finite-difference trapped weights even with `new_trap_int = T` while the
   conservation integrals use `wl`/`wxi` from `new_trap_int` — replicate literally (yes, via export), but
   confirm this is intended.
10. Edge case `1 - al*bmag == -bouncefuzz` exactly: `forbid` is false but `jend` does not count it — ignore?

---------------------------------------------------------------------------------------------------

## 10. Export checklist (for the exporter phase)

Static (always): `vnew, vnew_s, vnew_D, vnew_E, delvnew` `(naky,negrid,nspec)`; `vnewh` (if
`hyper_colls`); `pitch_weights`; `vpdiff`, `vpdiffle` (if conservative); `aj0le, vperp_aj1le, vpa_aj0_le`
(LE-canonical, negrid+1); helper `e15`, `capg_half`; all scalars/flags of §3.2; `vnmult` (post-init);
from other modules: `code_dt, tunits, kperp2, aky, ieqzip, beta, fapar, fphi, fbpar, zeff, bmag, al,
wl, wxi, xi, ixi_to_il, ixi_to_isgn, jend, forbid, w, energy_maxw, speed_maxw, w_maxw, ng2, nlambda,
negrid, nxi, aj0, aj1, vpa, vperp2, nonmaxw_corr, spec(:)%{vnewk,z,dens,stm,zstm,smz,zt,nu_h,type}, layout`.

dt-dependent (export at the initial dt for validation): `c1le, betaale, qle, ec1le, ebetaale, eqle`
(LE-canonical, negrid), `z0le, s0le, w0le, bz0le, bs0le, bw0le` (LE-canonical, negrid+1); for
`use_le_layout=F` additionally `c1/betaa/ql` (expand `ilz → (ig,it,ik,ie,is)`), `ec1/ebetaa/eql`
(expand `ielo → (ig,isgn,it,ik,il,is)`), `z0,s0,w0,bz0,bs0,bw0` (g-canonical).

Reference-run requirements: 1 MPI rank, `OMP_NUM_THREADS=1`, trapped grid (or a gated configuration
check), dump after `set_initial_field_and_dist_fn_values` so `vnmult` is final.
