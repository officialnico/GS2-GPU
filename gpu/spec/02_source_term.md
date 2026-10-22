# 02 — Source term (right-hand side) of the implicit GK equation

Status: contract for implementers of the JAX port. Read-only analysis of GS2 8.2.1
(`gpu-port` branch, commit `4d8c94bcf`). Every statement cites `file:line` of the
Fortran source in `/Users/nico/Workspace/gs2/src`.

Scope (this document):

* `dist_fn.fpp`: `get_source_term_on_local_domain`, `get_source_term`, `set_source`,
  `get_source_term_opt`, `set_source_opt`, the source part of `invert_rhs`
  (`sourcefac`, `g_work` handling), and how `invert_rhs_1` consumes the result.
* Coefficient initialisation used by the source: `init_vperp2`, `init_vpa_vpac`,
  `init_vpar`, `init_wdrift` (+ `wdrift_func`, `wcoriolis_func`), `init_wstar`,
  `init_bessel`, `init_source_term` (`ufac`, `source_coeffs_*`), the `a`/`b` part of
  `init_invert_rhs`, `hyper::D_res`.
* Explicit (Adams–Bashforth) history: `nonlinear_terms.fpp`
  `add_explicit_terms`, `add_explicit`, `nl_center`; `gs2_time.f90`
  `get_adams_bashforth_coefficients`, `time_history_available`, `update_time`.
  The bracket itself (`add_nl`) is out of scope (other spec).
* Flow shear as it touches the source (parallel-shear drive term; what the `exb_shear`
  remap does and does *not* touch).

Out of scope (other specs): the tridiagonal sweep / bounce / linked BCs in
`invert_rhs_1` + `apply_linked_boundary_conditions`, `r`/`ainv`/`g_h`, the field solve,
collisions, hyperviscosity `hyper_diff`, the nonlinear bracket `add_nl`, the
`exb_shear` remap algorithm in detail.

---------------------------------------------------------------------------------------

## 0. Conventions and build facts

### 0.1 Index conventions used in this document

* `ig ∈ [-ntgrid, ntgrid]` theta grid point. Python index `j = ig + ntgrid`,
  `T = 2*ntgrid+1` points. A *cell* `c` is the interval `[ig, ig+1]` for
  `ig ∈ [-ntgrid, ntgrid-1]`, Python index `c = ig + ntgrid ∈ [0, T-2]`, `N = T-1` cells.
* `isgn = 1` is `v_|| > 0` (rightward sweep), `isgn = 2` is `v_|| < 0`.
* Distribution-like arrays (`g`, `gexp_*`, `wdrift`, `vpar`, `vpac`, `vpa`, `a`, `b`,
  `r`, `ainv`, the source `g_work`), Fortran `(ig, isgn, iglo)`, are exported expanded as
  Fortran `(ig, isgn, it, ik, il, ie, is)` →
  **NumPy `(nspec, negrid, nlambda, naky, ntheta0, 2, T)`** (call this the G-shape,
  axes `S,E,L,Y,X,σ,T`).
* `(ig, iglo)` arrays (`aj0`, `aj1`, `vperp2`) → NumPy `(S, E, L, Y, X, T)`.
* Fields `phi, apar, bpar, phinew, aparnew, bparnew` Fortran `(ig, it, ik)` →
  NumPy `(naky, ntheta0, T)`.
* Velocity/theta arrays keep Fortran order; the NumPy view is the reversed index list
  (standard netCDF-from-Fortran behaviour). Given explicitly in §2.
* Cell arrays declared `(-ntgrid:ntgrid-1, …)` (`source_coeffs_*`) have length `N`
  on the theta axis.

### 0.2 Build facts that matter for round-off

* Double precision via `-fdefault-real-8 -fdefault-double-8`
  (`Makefiles/Makefile.gcc:18`, `Makefile:42` `DBLE ?= on`): every `real` *and every
  unsuffixed literal* (`0.5`, `2.0`, `3.`, `1.e8`, `epsilon(0.0)`) is binary64;
  `complex` is complex128. `epsilon(0.0) = 2.220446049250313e-16`.
* `-O3` (`Makefile:38` `OPT ?= on`, `Makefiles/Makefile.gcc:23`), gfortran on
  Apple arm64 (`bin/gs2` links `libgfortran.5`). GCC's default `-ffp-contract=fast`
  means **FMA contraction is likely**; see §6.1.
* The shipped binary is built **without OpenMP** (no `libgomp`; `otool -L bin/gs2`) and
  **without LOWFLOW** (`nm objs/dist_fn.o` has no `lowflow` symbols). Only the
  `#ifndef LOWFLOW` branches are specified below; LOWFLOW is GATED (§4).
* Complex arithmetic in gfortran follows `-fcx-fortran-rules` (default): complex
  multiplication is the naive `(ar*br - ai*bi, ar*bi + ai*br)` with no NaN recovery.
  GCC's complex-lowering pass tracks operands whose imaginary (or real) part is known
  to be zero (mixed real×complex, `zi*real`), and then emits only the needed real
  products. We mirror that with the helpers in §3.0.
* Fortran expression evaluation: binary operators of equal precedence associate
  left-to-right; unary minus has the precedence of binary `+/-` so `-a*b*c` is
  `-((a*b)*c)` and `-a*b - c*d` is `(-(a*b)) - (c*d)`. `x**2` with literal integer 2
  is `x*x`. `x**1.5` is a call to `pow`. gfortran does not reassociate without
  `-ffast-math`. All pseudo-code below is written fully parenthesised to fix the order.

---------------------------------------------------------------------------------------

## 1. Control flow / call graph

### 1.1 Initialisation (what computes the coefficients, and what is redone on dt change)

`gs2_init` levels (`gs2_init_level_list.inc`): `hyper`=19, `dist_fn_parameters`=21,
`nonlinear_terms`=24, `dist_fn_arrays`=26, `dist_fn_level_1`=27, `dist_fn_level_2`=28,
`override_timestep`=29, `collisions`=30, `dist_fn_level_3`=31, `fields_level_2`=32,
`set_initial_values`=34, `full`=35.

```
init_dist_fn_parameters                       dist_fn.fpp:849
  └ read_parameters                           dist_fn.fpp:1174   (knobs; fexp, bkdiff; tpdriftknob default;
                                                                  bkdiff forced equal; btor_slab→itor_over_B)
init_dist_fn_arrays                           dist_fn.fpp:885
  ├ init_gs2_time → adjust_time_norm          gs2_time.f90:45, 83   (wunits, tunits)
  ├ init_nonlinear_terms                      nonlinear_terms.fpp:287 (nonlin, split_nonlinear…)
  └ allocate_arrays                           dist_fn.fpp:3923   (g, gnew, g_work, g_h zeroed;
                                                                  gexp_1..3 allocated+zeroed iff
                                                                  nonlin .and. .not. split_nonlinear, :3970-3981;
                                                                  source_coeffs_* iff opt_source, :3951-3957;
                                                                  kx_shift/theta0_shift iff |g_exb*g_exbfac|>eps, :3983-3991)
init_dist_fn_level_1                          dist_fn.fpp:928
  ├ init_vperp2                               dist_fn.fpp:1560
  ├ init_vpa_vpac                             dist_fn.fpp:1592
  └ init_bc                                   dist_fn.fpp:2128   (links; not used by the source)
init_dist_fn_level_2                          dist_fn.fpp:939
  ├ init_bessel                               dist_fn.fpp:1756   (aj0, aj1; guarded by bessinit, :1768)
  └ init_fieldeq                              (field spec)
init_hyper                                    hyper.f90:359      (D_res; level 19, guarded by `initialized`)
init_dist_fn_level_3                          dist_fn.fpp:960
  ├ init_vpar                                 dist_fn.fpp:1681   ← dt
  ├ init_wdrift                               dist_fn.fpp:1405   ← dt
  ├ init_wstar                                dist_fn.fpp:1732   ← dt
  ├ init_collisions
  ├ init_invert_rhs                           dist_fn.fpp:1826   ← dt (a, b, r, ainv, then g_h)
  ├ init_hyper  (returns immediately: already initialised)
  └ init_source_term                          dist_fn.fpp:2159   (ufac once; source_coeffs_* ← dt, iff opt_source)
init_fields → fields_init_response (calls timeadv(...,0)) → nl_finish_init (initializing=.false.)
            → kwork_filter(1,1)=.true. if grid_option='box'                      fields.f90:354-377
```

**Time-step change** (`gs2_reinit::reset_time_step`, `gs2_reinit.f90:104-224`): goes down
to `override_timestep` (29) (`:166`), changes `code_dt` (`:174-198`), comes back to `full`
(`:216`). Therefore:

* Recomputed with the new `code_dt`: everything in `init_dist_fn_level_3`
  (`finish_dist_fn_level_3`, `dist_fn.fpp:1140-1154`, deallocates `vpar`, `wdrift`,
  `wdriftttp`, `wstar`; `a,b,r,ainv` are re-zeroed and recomputed in `init_invert_rhs`;
  `source_coeffs_*` recomputed in `init_source_term`), collisions, response matrix.
* **Not** recomputed: `ufac` (allocated once, `dist_fn.fpp:2180`), `aj0`, `aj1`, `vperp2`,
  `vpa`, `vpac` (dt independent anyway), **and `hyper::D_res`**, whose formula contains
  `code_dt` (`hyper.f90:403`) but whose module (level 19) is never re-initialised on a dt
  change (`init_hyper` returns at `hyper.f90:370`). `D_res` therefore stays at its value for
  the *initial* `code_dt` for the whole run. Replicate this (§6.7).
* Not touched: `gexp_1..3`, `code_dt_prev1/2`, `code_time`, `nonlinear_terms::istep_last`,
  `kx_shift`/`theta0_shift`.

### 1.2 One time step

```
gs2_main::evolve_equations, loop istep = istep_end+1 …        gs2_main.fpp:717
  advance(istep)                                              :730 → fields.f90:586
    advance_implicit(istep)                                   fields_implicit.f90:315
      [antenna]                                               :333
      if kx_shift/theta0_shift allocated: exb_shear(gnew, phinew, aparnew, bparnew, istep)   :335
      g ← gnew ; phi ← phinew ; apar ← aparnew ; bpar ← bparnew                              :337-338
      timeadv(phi, apar, bpar, phinew, aparnew, bparnew, istep)       # FIRST CALL  :341
           (here phinew == phi, aparnew == apar, bparnew == bpar bit-for-bit)
      if reset: return                                                              :343
      [apar_ext] ; getfield(...) ; phinew += phi ; …                                :345-356
      timeadv(phi, apar, bpar, phinew, aparnew, bparnew, istep, 1)    # SECOND CALL :365
      collisions_advance                                                            :369
  if reset (immediate CFL reset): reset_time_step; advance(istep) again              gs2_main.fpp:739-752
  update_time                                                 gs2_main.fpp:771 → gs2_time.f90:161
       code_dt_prev2 ← code_dt_prev1 ; code_dt_prev1 ← code_dt ; code_dt_old ← code_dt ;
       code_time ← code_time + code_dt ; user_time ← user_time + user_dt
  check_time_step → maybe reset_time_step (non-immediate path)  gs2_main.fpp:800-815
```

`fields_local::advance_local` (`fields_local.f90:4244-4342`) and
`fields_gf_local::advance_gf_local` follow the same pattern (exb_shear, copy, timeadv,
field solve, timeadv); in all of them `phinew==phi` on the first call.

`timeadv` (`dist_fn.fpp:4014-4126`):

```
if split_nonlinear .and. nonlin:   (GATED, other algorithm)                    :4046-4058
else: add_explicit_terms(gexp_1, gexp_2, gexp_3, phi, apar, bpar, istep, bkdiff(1))  :4060
if (reset .and. immediate_reset) return                                         :4066
invert_rhs(phi, apar, bpar, phinew, aparnew, bparnew, istep)                    :4068
hyper_diff, vspace_derivatives, enforce_parity, (split NL second half)          :4073-4111
istep_last = istep                                                              :4124
```

`invert_rhs` (`dist_fn.fpp:6152-6193`):

```
time = code_time                                                    :6166
sourcefac = (time > t0) ? source0*exp(-zi*omega0*time + gamma0*time)
                        : (0.5 - 0.5*cos(pi*time/t0))*exp(-zi*omega0*time + gamma0*time)   :6172-6176
g_work ← 0                                                          :6179
get_source_term_on_local_domain(g_work, phi, apar, bpar, phinew, aparnew, bparnew, istep, sourcefac)  :6181
invert_rhs_1(phinew, bparnew, source = g_work)                      :6184   (other spec)
if linked: apply_linked_boundary_conditions                         :6186-6191 (other spec)
```

`get_source_term_on_local_domain` (`dist_fn.fpp:5363-5420`):

```
include_explicit = (istep > 0) .and. nonlin .and. (.not. split_nonlinear)        :5382-5386
if include_explicit: ab = get_adams_bashforth_coefficients()                     :5388-5389
for iglo:                                                                        :5397
   if kwork_filter(it,ik): skip (g_work stays 0 for every ig, isgn)              :5401
   if opt_source: get_source_term_opt(..., source(:,:,iglo))                     :5407-5410  (GATED)
   else: for isgn=1,2: get_source_term(..., isgn, ..., source(:,isgn,iglo))      :5412-5416
```

Note: `ieqzip` is **not** applied here (only in `invert_rhs_1`, `dist_fn.fpp:5477`); the
source is computed for `ieqzip` modes and simply not used.

### 1.3 Other callers

* Response-matrix construction: `fields_implicit.f90:669` (and
  `fields_local.f90:2184`, `fields_gf_local.f90:1955`) call `timeadv(..., 0)` →
  `include_explicit = .false.` and `add_explicit_terms` does nothing for `istep=0`
  (`nonlinear_terms.fpp:424-426`). The source formulas are otherwise identical.
* Debug hook `dump_current_source_term` (`dist_fn.fpp:7850-7876`, **public**): zeroes
  `g_work`, calls `get_source_term_on_local_domain` with the given fields/istep/sourcefac
  and writes `g_work` with `gs2_save_for_restart(fileopt='.source')`. Ideal oracle for unit
  tests of this spec (§8).

### 1.4 Output contract consumed by `invert_rhs_1`

`g_work(ig, isgn, iglo)` complex, G-shape:

* Entries `ig ∈ [-ntgrid, ntgrid-1]` are set as specified in §3.
* `ig = ntgrid`: never written by `set_source` / matrix multiply / ttp code
  (comments at `dist_fn.fpp:4946-4955`), so it is the 0 from `zero_array(g_work)`
  (`:6179`) — except `il > lmax` or `homogeneous` where `source = 0.0` is assigned to the
  whole column anyway. **Always 0.**
* `kwork_filter(it,ik)` points: all zeros.
* Consumers (`invert_rhs_1`): `source(ig,2,iglo)` for all `ig ≤ ntgrid-1` (`:5596-5598`),
  `source(ig,1,iglo)` for `ig ≤ ntgrid-1` (`:5634-5645`), and `apply_parallel_boundary_conditions`
  reads `source(-ntgrid,1)` / `source(ntgrid-1,2)` (`:5753`, `:5784`, `:5810`, `:5826`). `invert_rhs_1`
  **overwrites** `source(ig,1,iglo) = gnew(ig+1,2,iglo)` at cells left of lower bounce
  points (`:5612-5621`); the JAX port should treat that as a local modification in the
  invert step, not as part of the source.

---------------------------------------------------------------------------------------

## 2. Data inventory

Legend: **dt** = depends on `code_dt`; **RT** = other runtime change. Decision:
EXPORT (dump from Fortran), RECOMPUTE (port formula), BOTH (export at initial dt for
verification *and* port formula so it can be rebuilt when dt changes).

Important practical note for the exporter: several of these are **private** module
variables of `dist_fn` (`wdrift`, `wdriftttp`, `vpar`, `a`, `b`, `r`, `ainv`, `ufac`,
`fexp`, `bkdiff`, `source_coeffs_*`, `omprimfac`, `mach`, `phi_ext`, `t0`, `omega0`,
`gamma0`, `source0`, `source_option_switch`, `opt_source`, `driftknob`, `tpdriftknob`)
and of `nonlinear_terms` (`zip`, `nl_forbid_force_zero`, `istep_last`). The public
`get_dist_fn_config()` returns **raw input values, not effective ones** (e.g.
`tpdriftknob=-9.9e9` sentinel, un-equalised `bakdif`, pre-override `g_exb`/`mach`/`opt_source`).
The export phase needs an accessor in `dist_fn` (and `nonlinear_terms`) that returns the
effective module values.

### 2.1 Runtime state (inputs to every call)

| Fortran name | module | decl / ranges | NumPy shape | notes | decision |
|---|---|---|---|---|---|
| `phi, apar, bpar` | `fields_arrays` | `complex(-ntgrid:ntgrid, ntheta0, naky)` | `(Y,X,T)` | fields at t^n (copied from `*new` at step start, `fields_implicit.f90:338`) | EXPORT initial state |
| `phinew, aparnew, bparnew` | `fields_arrays` | same | `(Y,X,T)` | 1st call: == `phi`; 2nd call: after field solve | EXPORT initial state |
| `g` | `dist_fn_arrays` | `complex(-ntgrid:ntgrid, 2, g_lo)` | G-shape | g^n (copy of `gnew`, `fields_implicit.f90:337`) | EXPORT initial state |
| `gexp_1, gexp_2, gexp_3` | `dist_fn_arrays` | `complex(-ntgrid:ntgrid, 2, g_lo)`; allocated iff `nonlin .and. .not. split_nonlinear` (`dist_fn.fpp:3970`) | G-shape | AB history: newest, previous, oldest NL source (already `nl_center`ed, i.e. ≈2× cell value) | EXPORT initial (zeros on fresh start); then evolved |
| `code_dt` | `gs2_time` | real | scalar | changes on reset | EXPORT; evolve |
| `code_dt_prev1, code_dt_prev2` | `gs2_time` | real, sentinel `dt_not_set=-1.0` (`gs2_time.f90:18-23`) | scalars | set by `update_time` / restart (`delt_option`) | EXPORT; evolve |
| `code_time` | `gs2_time` | real | scalar | `sourcefac` only | EXPORT; evolve (`code_time = code_time + code_dt`) |
| `code_dt_old` | `gs2_time` | real, init 0 | scalar | used by exb_shear (`gdt`) only | EXPORT; evolve |
| `istep` | — | integer | scalar | loop counter; explicit terms need `istep>0` | driver |
| `nonlinear_terms::istep_last` | private | integer, init 0 (`nonlinear_terms.fpp:32`) | scalar | history shuffle guard | RECOMPUTE (deterministic: = last istep advanced) |
| `kwork_filter` | `kt_grids` | `logical(ntheta0, naky)` | `(Y,X)` | `.false.` except `(1,1)=.true.` for box after `init_fields` (`fields.f90:376`) | EXPORT (value during time loop) |

### 2.2 Coefficients read by the source

| Fortran name | module | decl / ranges | NumPy shape | computed (file:line) | dt | RT | decision |
|---|---|---|---|---|---|---|---|
| `aj0` | `dist_fn_arrays` | `real(-ntgrid:ntgrid, g_lo)` | `(S,E,L,Y,X,T)` | `init_bessel` `dist_fn.fpp:1784-1793` (`spfunc::j0`, local spfunc on macOS) | no | no | EXPORT (Bessel impl. not portable bit-exactly) |
| `aj1` (= J1(x)/x) | `dist_fn_arrays` | same | same | same (`j1`) | no | no | EXPORT |
| `vperp2` | `dist_fn_arrays` | `real(-ntgrid:ntgrid, g_lo)` | `(S,E,L,Y,X,T)` (depends on is,ie,il,ig only) | `init_vperp2` `:1585` | no | no | EXPORT |
| `vpa` | `dist_fn_arrays` | `real(-ntgrid:ntgrid,2,g_lo)` | G-shape (depends on is,ie,il,ig,σ) | `init_vpa_vpac` `:1646-1652` | no | no | EXPORT (needed to recompute wdrift) |
| `vpac` | `dist_fn_arrays` | same | G-shape | `:1666-1675` | no | no | EXPORT |
| `vpar` | `dist_fn` (private) | `real(-ntgrid:ntgrid,2,g_lo)` | G-shape (no it dependence) | `init_vpar` `:1715-1725` | **yes** | dt reset | BOTH |
| `wdrift` | `dist_fn` (private) | `real(-ntgrid:ntgrid,2,g_lo)` (complex only with LOWFLOW, `:133-139`) | G-shape | `init_wdrift` `:1423-1491` + `wdrift_func` `:1500-1531` + `wcoriolis_func` `:1534-1557` | **yes** | dt reset; `mach` override | BOTH |
| `wdriftttp` | `dist_fn` (private) | same | G-shape (non-zero only where `is_ttp`) | `:1451-1480` | **yes** | dt reset | BOTH |
| `wstar` | `dist_fn_arrays` | `real(naky, negrid, nspec)` | Fortran `(Y,E,S)` → NumPy `(S,E,Y)` | `init_wstar` `:1746` | **yes** | dt reset | BOTH (trivial) |
| `a`, `b` | `dist_fn` (private) | `complex(-ntgrid:ntgrid,2,g_lo)`; entries `ig=ntgrid` are 0 | G-shape | `init_invert_rhs` `:1892-1897`, ttp override of `a` `:1919` | **yes** | dt reset | BOTH |
| `r`, `ainv` | `dist_fn` (private) | same | G-shape | `:1885-1891`, `:1899-1921` | **yes** | dt reset | BOTH (owned by invert spec; listed because they share the init loop) |
| `ufac` | `dist_fn` (private) | `real(negrid, nspec)` | Fortran `(E,S)` → NumPy `(S,E)` | `init_source_term` `:2180-2187` | no | never recomputed | EXPORT (contains `pow(x,1.5)`) |
| `D_res` | `hyper` | `real(ntheta0, naky)` | Fortran `(X,Y)` → NumPy `(Y,X)` | `init_hyper` `hyper.f90:400-410`; 0 otherwise (`hyper.f90:615`) | formula yes, **but frozen at initial dt** | — | EXPORT (do not rescale on dt change) |
| `source_coeffs_phim/phip/aparm/aparp` | `dist_fn` (private) | `complex(-ntgrid:ntgrid-1, 2, g_lo)`; `apar*` only if `has_apar` | `(S,E,L,Y,X,2,N)` | `init_source_term` `:2238-2259` | **yes** | dt reset | BOTH, only if `opt_source` (GATED) |
| `nonmaxw_corr` | `species` | `real(negrid, nspec)` | `(S,E)` | `species.f90:571-592` (1.0 for Maxwellian) | no | profile overrides | EXPORT |
| `dlnf0drho` | `species` | `real(negrid, nspec)` | `(S,E)` | `species.f90:587` (Maxwellian: `-fprim - (E-1.5)*tprim`) | no | profile overrides | EXPORT |
| `fexp` (= `fexpr`) | `dist_fn` (private) | `real(nspec)` | `(S,)` | `read_species_knobs` `:1398` | no | no | EXPORT |
| `bkdiff` (= `bakdif`) | `dist_fn` (private) | `real(nspec)`; forced `=bkdiff(1)` unless `mult_imp` (linear) (`:1330-1344`) | `(S,)` | `:1397` | no | no | EXPORT (effective) |
| `spec%zstm, tz, zt, stm, smz` | `species` | real per species | `(S,)` | `species.f90:476-480` | no | no | EXPORT |
| `wunits, tunits` | `gs2_time` | `real(naky)` | `(Y,)` | `adjust_time_norm` `gs2_time.f90:108-125` (`wunits=aky/2`, `tunits=1` unless `wstar_units`) | no | no | EXPORT |
| `aky, akx` | `kt_grids` | `real(naky)`, `real(ntheta0)` | `(Y,)`, `(X,)` | kt_grids | no | no | EXPORT |
| `theta0` | `kt_grids` | `real(ntheta0, naky)` | `(Y,X)` | kt_grids | no | no | EXPORT |
| `kperp2` | `kt_grids` | `real(-ntgrid:ntgrid, ntheta0, naky)` | `(Y,X,T)` | `kt_grids.f90:313-349` | no | no | EXPORT (only if aj0 recompute is ever attempted) |
| `bmag, gradpar, delthet, gbdrift, gbdrift0, cvdrift, cvdrift0, cdrift, cdrift0, itor_over_B` | `theta_grid` | `real(-ntgrid:ntgrid)`; `delthet(ntgrid)=0` (`theta_grid.f90:511`) | `(T,)` | geometry; drift knobs already applied (`theta_grid.f90:623-626`); `itor_over_B` overridden by `btor_slab` (`dist_fn.fpp:1369`) | no | no | EXPORT (effective values) |
| `shat` | `theta_grid` | real | scalar | | no | no | EXPORT |
| `energy` | `le_grids` | `real(negrid, nspec)` | `(S,E)` | | no | no | EXPORT |
| `al` | `le_grids` | `real(nlambda)` | `(L,)` | | no | no | EXPORT |
| `forbid` | `le_grids` | `logical(-ntgrid:ntgrid, nlambda)` | `(L,T)` | `le_grids.f90:2715-2738` | no | no | EXPORT |
| `is_ttp` | `le_grids` | `logical(-ntgrid:ntgrid, nlambda)` | `(L,T)` | `le_grids.f90:2907-2946` | no | no | EXPORT |
| `can_be_ttp` | `le_grids` | `logical(nlambda)` | `(L,)` | same | no | no | EXPORT |
| `is_bounce_point` | `le_grids` | `logical(-ntgrid:ntgrid, nlambda)` | `(L,T)` | `le_grids.f90:2746-2771` | no | no | EXPORT (needed for `ainv`; `is_lower_bounce_point` `:2775-2819` derived) |
| `jend` | `le_grids` | `integer(-ntgrid:ntgrid)` | `(T,)` | `le_grids.f90:2871-2903` | no | no | EXPORT (nl_center) |
| `ng2, nlambda, lmax` | `le_grids` | integers; `lmax = nlambda-1` if trapped grid else `nlambda` (`le_grids.f90:1897-1903`) | scalars | | no | no | EXPORT |
| `fphi, fapar, fbpar, has_phi, has_apar, has_bpar` | `run_parameters` | `fbpar` default `-1.0`; `beta=0` ⇒ `fapar=fbpar=0` (`run_parameters.f90:340-352`) | scalars | | no | no | EXPORT (effective) |
| knobs `driftknob, tpdriftknob(eff.), mach, g_exb, omprimfac, phi_ext, t0, omega0, gamma0, source0, source_option_switch, opt_source, nonlin, split_nonlinear, zip, nl_forbid_force_zero` | see §4 | | scalars | | no | `g_exb`,`mach` profile overrides; `opt_source` optimisation override (`dist_fn.fpp:1157-1172`) | EXPORT (effective) |
| AB coefficients | computed per call | `real(:)`, size 1..3 | `(1..3,)` | `gs2_time.f90:250-312` | **yes** | dt reset, step count | RECOMPUTE |
| `sourcefac` | local | complex | scalar | `dist_fn.fpp:6172-6176` | via `code_time` | each step | RECOMPUTE |

---------------------------------------------------------------------------------------

## 3. Exact computations

### 3.0 Arithmetic helpers (mirror GCC complex lowering)

All values are float64 / complex128 (`jax.config.update("jax_enable_x64", True)`).

```python
def rmul(z, r):        # Fortran  z*r  or  r*z  with r real, z complex  (componentwise)
    return lax.complex(z.real * r, z.imag * r)

def imul(w, z):        # Fortran  (zi*w)*z  with w real:  zi*w == (0, w), then (0,w)*z
    return lax.complex(-(w * z.imag), w * z.real)

def cmul(x, y):        # general complex*complex (both VARYING)
    return lax.complex(x.real*y.real - x.imag*y.imag, x.real*y.imag + x.imag*y.real)
```

Rationale: `zi` is the parameter `(0.0, 1.0)` (`externals/utils/constants.fpp:78`);
`zi*w` for real `w` has an exactly-zero real part and imaginary part `1.0*w = w`; a
complex with known-zero real part times a general complex costs two real products and
no addition. Using plain `1j*w*z` in JAX gives the same values except for signs of zero
(and NaN propagation) – acceptable, but the helpers make intent explicit and avoid
useless FMAs. Complex `+`/`-` are componentwise. A real added to a complex only touches
the real part.

### 3.1 Coefficient formulas (init)

Notation per element: `s=is, e=ie, l=il, y=ik, x=it, σ=isgn, ig`. All loops are over
independent elements; the order of the *outer* loops is irrelevant for values.

**vperp2** (`dist_fn.fpp:1584-1587`):
```
vperp2[s,e,l,y,x,ig] = forbid(ig,l) ? 0.0 : (bmag(ig)*al(l))*energy(e,s)
```

**vpa, vpac** (`:1646-1676`):
```
vpa[...,1,ig] = sqrt(energy(e,s)*max(0.0, 1.0 - al(l)*bmag(ig))) * vparknob
vpa[...,2,ig] = -vpa[...,1,ig]
if forbid(ig,l): vpa[...,1,ig] = vpa[...,2,ig] = 0.0
for ig in [-ntgrid, ntgrid-1]:
    if forbid(ig,l) or forbid(ig+1,l): vpac[...,σ,ig] = 0.0
    else: vpac[...,σ,ig] = 0.5*(vpa[...,σ,ig] + vpa[...,σ,ig+1])
vpac[...,σ,ntgrid] = 0.0
```
(Consequence: `vpac(:,2) == -vpac(:,1)` bit-for-bit.)

**vpar** (`:1715-1725`), for `ig ∈ [-ntgrid, ntgrid-1]`, both σ:
```
vpar[s,e,l,y,x,σ,ig] = ((((((zstm(s)*tunits(y))*code_dt)*0.5)/delthet(ig))
                         *(abs(gradpar(ig)) + abs(gradpar(ig+1))))*vpac[s,e,l,y,x,σ,ig])
vpar[..., ntgrid] = 0.0
```
(`vpar(:,2) == -vpar(:,1)` exactly.)

**wdrift_func** (`:1500-1531`) — at a grid point, real:
```
vperp2_l = vperp2[s,e,l,y,x,ig];  vpa2 = vpa[s,e,l,y,x,1,ig]*vpa[s,e,l,y,x,1,ig]
if aky(y) == 0.0:      # exact compare
    wd = ((((akx(x)/shat) * (cvdrift0(ig)*vpa2 + (gbdrift0(ig)*0.5)*vperp2_l)) * code_dt) / 2.0)
else:
    wd = ((((cvdrift(ig) + theta0(x,y)*cvdrift0(ig))*vpa2
            + ((gbdrift(ig) + theta0(x,y)*gbdrift0(ig))*0.5)*vperp2_l) * code_dt) * wunits(y))
```
Fortran text (zonal): `akx(it)/shat*(cvdrift0(ig)*vpa2_local + gbdrift0(ig)*0.5*vperp2_local)*code_dt/2.0`.

**wcoriolis_func** (`:1534-1557`):
```
if aky(y) == 0.0:
    wc = (((((mach*vpa[...,1,ig])*cdrift0(ig))*code_dt)*akx(x)) / ((2.0*shat)*stm(s)))
else:
    wc = (((((mach*vpa[...,1,ig])*(cdrift(ig) + theta0(x,y)*cdrift0(ig)))*code_dt)*wunits(y)) / stm(s))
```

**init_wdrift** (`:1405-1498`):
```
wdrift = 0 ; wdriftttp = 0
# pass 1 (grid values)
knob = tpdriftknob if (il_is_wfb(l) or il_is_trapped(l)) else driftknob   # wfb: l == ng2+1; trapped: l > ng2+1
for ig in [-ntgrid, ntgrid]:
    if forbid(ig,l): continue                     # stays 0
    wd = wdrift_func(ig) ; wc = wcoriolis_func(ig)
    wdrift[...,1,ig] = wd*knob + wc
    wdrift[...,2,ig] = wd*knob - wc
# pass 2 (ttp), only if can_be_ttp(l)
for ig in [-ntgrid, ntgrid]:
    if not is_ttp(ig,l): continue
    wd = wdrift_func(ig) ; wc = wcoriolis_func(ig)
    wdriftttp[...,1,ig] = wd*tpdriftknob + wc
    wdriftttp[...,2,ig] = wd*tpdriftknob - wc
# pass 3 (centre wdrift only; NOT wdriftttp)
for ig in [-ntgrid, ntgrid-1]:                    # ascending; reads ORIGINAL wdrift(ig+1)
    wdrift[...,σ,ig] = 0.5*(wdrift[...,σ,ig] + wdrift[...,σ,ig+1])
# wdrift[..., ntgrid] keeps its grid value
```
Vectorised: `wc_new[:, :-1] = 0.5*(w[:, :-1] + w[:, 1:])`, `wc_new[:, -1] = w[:, -1]`
(pure function of the pass-1 array; no recurrence). Forbidden grid points contribute 0 to
the average of a bounce cell.

**wstar** (`:1746`):
```
wstar[y,e,s] = -((code_dt*wunits(y))*dlnf0drho(e,s))
```

**ufac** (`:2183-2184`) — export; for reference:
```
ufac[e,s] = (2.0*uprim(s)) + (((uprim2(s)*(energy(e,s)**1.5))*sqrt(pi))/4.0)
```

**D_res** (`hyper.f90:385-410`) — export; for reference (only if
`hyper_option ∈ {visc_only,res_only,both}` and `D_hyperres > 0`, else 0):
```
akperp4_max = gridnorm ? (akx(ntheta0/2+1)**2 + aky(naky)**2)**nexp : 1.0
D_res[x,y] = ((D_hyperres*code_dt_INITIAL)*((aky(y)**2 + akx(x)**2)**nexp))/akperp4_max
```

**a, b** (`dist_fn.fpp:1866-1924`), for `ig ∈ [-ntgrid, ntgrid-1]`, both σ, with
`wd = wdrift[...,σ,ig]` (σ-dependent), `vp = vpar[...,1,ig]` (**always σ=1, positive**),
`bd = bkdiff(s)`, `fe = fexp(s)`, `tz = tz(s)`:
```
a = (1.0 + bd) + (fe*tz) * ( (-( (zi*wd)*(1.0+bd) )) - 2.0*vp )
b = (1.0 - bd) + (fe*tz) * ( (-( (zi*wd)*(1.0-bd) )) + 2.0*vp )
```
i.e. in components (exact up to signed zero):
```
a.re = (1.0+bd) + (fe*tz)*(-(2.0*vp))        a.im = (fe*tz)*(-(wd*(1.0+bd)))
b.re = (1.0-bd) + (fe*tz)*(2.0*vp)           b.im = (fe*tz)*(-(wd*(1.0-bd)))
```
Then, only `if grid_has_trapped_particles()` (`nlambda > ng2`) and `is_ttp(ig,l)`
(`:1917-1921`), with `wdttp = wdriftttp[...,σ,ig]`:
```
a = 1.0 - ((zi*fe)*tz)*wdttp        # a.re = 1.0, a.im = -((fe*tz)*wdttp)
```
`a` and `b` are **not** zeroed in forbidden cells (only `r`, `ainv` are, `:1901-1904`).
`a[..., ntgrid] = b[..., ntgrid] = 0` (never set). (For the invert spec, verbatim:
`ainv = 1/(1+bd + (1-fe)*tz*(zi*wd*(1+bd) + 2vp))`, `r = (1-bd + (1-fe)*tz*(zi*wd*(1-bd) - 2vp))*ainv`,
forbid → `r=ainv=0`; σ=1 and `is_lower_bounce_point(ig+1,l)` → `ainv = 1 + ainv`;
ttp → `ainv = 1/(1 + zi*(1-fe)*tz*wdttp)`, `r = 0`. Note the complex **division**.)

### 3.2 `sourcefac` (`dist_fn.fpp:6166-6176`)

```
time = code_time
arg  = (-((zi*omega0)*time)) + gamma0*time       # complex: (gamma0*time [+ -0], -(omega0*time))
if time > t0: sourcefac = source0 * exp(arg)
else:         sourcefac = (0.5 - 0.5*cos((pi*time)/t0)) * exp(arg)
```
Computed every call; only used if `source_option='phiext_full'` and `aky(ik) < epsilon(0.0)`.

### 3.3 Adams–Bashforth coefficients (`gs2_time.f90:234-312`)

```
order = 3 if code_dt_prev2 != -1.0 else (2 if code_dt_prev1 != -1.0 else 1)
dt = code_dt ; p1 = code_dt_prev1 ; p2 = code_dt_prev2
order 3:
  q  = dt + p1
  c1 = ( ((q*q)/p1) * (q/3. + p2/2.)  -  p1*(p1/3. + p2/2.) ) / (dt*(p1 + p2))
  c2 = -( (dt*(dt/3. + (p1+p2)/2.)) / (p1*p2) )
  c3 = (dt*(dt/3. + p1/2.)) / (p2*(p1 + p2))
order 2:
  c1 = (0.5*dt + p1)/p1
  c2 = -((0.5*dt)/p1)
order 1:
  c1 = 1.0
```
Fortran for c1 (verbatim, `:287-291`):
`(((this_dt + code_dt_prev1)**2/code_dt_prev1)*( (this_dt+code_dt_prev1)/3. + code_dt_prev2/2. ) - code_dt_prev1*(code_dt_prev1/3. + code_dt_prev2/2.) ) / (this_dt*(code_dt_prev1 + code_dt_prev2))`.
`(x)**2` → `x*x`; `/3.` is a true division. With equal steps this is *not* exactly
23/12, −4/3, 5/12 in floating point – use the formula.

### 3.4 `get_source_term` (default path, `opt_source=.false.`), per `(iglo, isgn)`

`dist_fn.fpp:4856-5130`. Inputs per call: `s,e,l,y,x,σ`. Local work vectors over the full
`ig ∈ [-ntgrid, ntgrid]`.

**Step A — time-decentred potentials** (`:4891-4896`), both over all `ig`:
```
fe = fexp(s) ; ofe = 1.0 - fexp(s)
phigavg[ig]  = ( ((fe*phi(ig,x,y) + ofe*phinew(ig,x,y)) * aj0[ig]) * fphi )
             + ( (((((fe*bpar(ig,x,y) + ofe*bparnew(ig,x,y)) * aj1[ig]) * fbpar) * 2.0) * vperp2[ig]) * tz(s) )
apargavg[ig] = ((fe*apar(ig,x,y) + ofe*aparnew(ig,x,y)) * aj0[ig]) * fapar
```
(each `real*complex` is `rmul`; `fe*phi + ofe*phinew` is componentwise.)
These are evaluated even when `fbpar`/`fapar` are 0 or `-1` (default `fbpar=-1.0` when
`beta≠0` without bpar; fields are then zero so the terms are ±0).

**Step B — select on `source_option_switch`** (`:4899-4927`):
```
full (1):        if l <= lmax: set_source()  else: source[all ig] = 0
phiext_full (5): if l <= lmax:
                     set_source()
                     if aky(y) < epsilon(0.0):
                         for ig in [-ntgrid, ntgrid-1]:
                             source[ig] = source[ig] - rmul( imul( ((wdrift[σ,ig]*nonmaxw(e,s))*2.0)*phi_ext, sourcefac ), aj0[ig] )
                 else: source[all ig] = 0
homogeneous (6): source[all ig] = 0
```
Fortran for the phiext term (`:4918-4920`):
`- zi*wdrift(:ntgrid-1,isgn,iglo)*nonmaxw_corr(ie,is)*2.0*phi_ext*sourcefac*aj0(:ntgrid-1,iglo)`
= `-(((((( zi*wdrift)*nonmaxw)*2.0)*phi_ext)*sourcefac)*aj0)`.

**set_source** (`:5017-5128`), for `ig ∈ [-ntgrid, ntgrid-1]`:
```
bd = bkdiff(s)
bdfac_p = 1.0 + bd*(3.0 - 2.0*real(σ))     # σ=1: 1+bd ; σ=2: 1-bd   (bd*(±1) exact)
bdfac_m = 1.0 - bd*(3.0 - 2.0*real(σ))     # σ=1: 1-bd ; σ=2: 1+bd

phi_p  = rmul(phigavg[ig+1], bdfac_p) + rmul(phigavg[ig], bdfac_m)
phi_m  = phigavg[ig+1] - phigavg[ig]
apar_p = rmul(apargavg[ig+1], bdfac_p) + rmul(apargavg[ig], bdfac_m)
apar_m = rmul(rmul(aparnew(ig+1)-apar(ig+1), bdfac_p), fapar)
       + rmul(rmul(aparnew(ig)  -apar(ig),   bdfac_m), fapar)
nm   = nonmaxw_corr(e,s)
vpr  = vpar[σ,ig] ; vpc = vpac[σ,ig] ; wdr = wdrift[σ,ig]

T1 = -( rmul(phi_m, (2.0*vpr)*nm) )
X  = rmul(apar_m, (aj0[ig+1] + aj0[ig])*0.5) + rmul(apar_p, D_res(x,y))
T2 = rmul(X, (zstm(s)*vpc)*nm)
T3 = imul(wdr*nm, phi_p)                     # ((zi*wdrift)*nm)*phi_p
A  = (T1 - T2) - T3

W  = ( wstar(y,e,s) + (((vpc*code_dt)*wunits(y))*ufac(e,s)) )
   - ( ((((((((2.0*omprimfac)*vpc)*nm)*code_dt)*wunits(y))*g_exb)*itor_over_B(ig)) / stm(s) )
P  = phi_p - rmul(rmul(apar_p, stm(s)), vpc)
source[ig] = A + imul(W, P)                  # zi*(W)*(P) = (zi*W)*P
```
Fortran (verbatim, `:5094-5102`):
```fortran
source(ig) = (-2.0*vpar(ig,isgn,iglo)*nonmaxw_corr(ie,is)*phi_m &
     -spec(is)%zstm*vpac(ig,isgn,iglo)*nonmaxw_corr(ie,is) &
     *((aj0(ig+1,iglo) + aj0(ig,iglo))*0.5*apar_m  &
     + D_res(it,ik)*apar_p) &
     -zi*wdrift(ig,isgn,iglo)*nonmaxw_corr(ie,is)*phi_p) &
     + zi*(wstar(ik,ie,is) &
     + vpac(ig,isgn,iglo)*code_dt*wunits(ik)*ufac(ie,is) &
     -2.0*omprimfac*vpac(ig,isgn,iglo)*nonmaxw_corr(ie,is)*code_dt*wunits(ik)*g_exb*itor_over_B(ig)/spec(is)%stm) &
     *(phi_p - apar_p*spec(is)%stm*vpac(ig,isgn,iglo))
```
Notes: `apar_m` is the plain (non-gyroaveraged) time difference, gyroaveraged by the
centred `aj0`; `D_res` multiplies the gyroaveraged `apar_p`. `g_exb` here is **not**
multiplied by `g_exbfac` and is **not** gated by `g_exb_start_time/timestep`.

Explicit terms inside `set_source` (`:5106-5127`), only if `include_explicit`:
```
k = size(ab)              # 1, 2 or 3
E = rmul(gexp_1[σ,ig], ab[0])
if k >= 2: E = E + rmul(gexp_2[σ,ig], ab[1])
if k == 3: E = E + rmul(gexp_3[σ,ig], ab[2])
source[ig] = source[ig] + rmul(E, 0.5*code_dt)
```
(`0.5*code_dt*(ab1*g1 + ab2*g2 + ab3*g3)` = `(0.5*code_dt)*((ab1*g1 + ab2*g2) + ab3*g3)`.
Padding missing coefficients with 0.0 and always summing 3 terms is exact up to signed
zeros: `x + 0.0*y == x` for finite `y`.)

**Step C — "matrix multiplication" with g^n** (`:4932-4944`), only if `l <= lmax`
(for all options including homogeneous), `ig ∈ [-ntgrid, ntgrid-1]`:
```
σ=1: source[ig] = (source[ig] + cmul(b[1,ig], g[1,ig])) + cmul(a[1,ig], g[1,ig+1])
σ=2: source[ig] = (source[ig] + cmul(a[2,ig], g[2,ig])) + cmul(b[2,ig], g[2,ig+1])
```
(Fortran writes `b*g` and `a*g`, i.e. coefficient first.)

**Step D — totally trapped particles** (`:4957-5012`). Condition:
`source_option ∈ {full, phiext_full}` and `nlambda > ng2` and `σ == 2` and `can_be_ttp(l)`.
For each `ig ∈ [-ntgrid, ntgrid-1]` with `is_ttp(ig,l)` (this **overwrites** whatever
Steps B–C produced at that point, including for `l > lmax`):
```
source[ig] = ( cmul(g[2,ig], a[2,ig]) - imul(wdriftttp[2,ig]*nm, phigavg[ig]) ) + imul(wstar(y,e,s), phigavg[ig])
if phiext_full and aky(y) < epsilon(0.0):
    source[ig] = source[ig] - rmul( imul( ((wdriftttp[2,ig]*nm)*2.0)*phi_ext, sourcefac ), aj0[ig] )
if include_explicit:
    source[ig] = source[ig] + rmul(E_ttp, 0.5*code_dt)    # same AB sum as above, σ=2
```
Fortran (`:4962-4971`): `source(ig) = g(ig,2,iglo)*a(ig,2,iglo) - zi*(wdriftttp(ig, 2, iglo))*nonmaxw_corr(ie,is)*phigavg(ig) + zi*wstar(ik,ie,is)*phigavg(ig)`
(note `g*a` order; grid-point `phigavg`, not centred; no bakdif; no ufac/g_exb term; no apar).
Explicit terms on a ttp point are therefore added exactly once (the `set_source` addition
is overwritten). ttp never occurs at `ig = ±ntgrid` (`le_grids.f90:2919`).

For σ=1 at ttp points the Steps B–C value is kept (it is not used by `invert_rhs_1`,
which sets `gnew(ig,1)=gnew(ig,2)` there, `dist_fn.fpp:5643`).

`il = nlambda` (> `lmax` on a trapped grid) therefore gets `source = 0` everywhere except
its ttp points, σ=2.

### 3.5 Per-element summary (vectorised reference)

```python
# canonical axes: S,E,L,Y,X,σ,T ; cells c = 0..T-2 ; slices lo=[..., :-1], hi=[..., 1:]
src = zeros(G, complex128)
pg  = phigavg            # (S,E,L,Y,X,T)
ag  = apargavg
base = set_source(...)   # (S,E,L,Y,X,2,N) per §3.4 (incl. explicit if include_explicit)
if option == PHIEXT: base = where(aky_zero[Y], base - phiext_term, base)
if option in (FULL, PHIEXT): body = where(l <= lmax, base, 0)
else:                        body = 0
mm = where(l <= lmax, matrix_terms(a, b, g), 0)          # §3.4 Step C (also for HOMOGENEOUS)
src[..., :-1] = body + mm        # NOTE: written as (body + term1) + term2 per element
ttp_mask = (option in (FULL,PHIEXT)) & trapped_grid & (σ==2) & can_be_ttp[L] & is_ttp[L, c]
src[..., :-1] = where(ttp_mask, ttp_value, src[..., :-1])
src = where(kwork_filter[Y,X], 0, src)
src[..., -1] = 0
```
Order inside one element must be `(((set_source) [ - phiext]) + b*g) + a*g`; do not
pre-sum `b*g + a*g`.

### 3.6 `get_source_term_opt` (`opt_source=.true.`, GATED; spec for later)

`dist_fn.fpp:5135-5356`, coefficients `init_source_term` `:2194-2264` (non-LOWFLOW branch `:2238-2259`). Differences vs §3.4
that change **results**, not just rounding:

* `apar_m = ((-apar(ig) - apar(ig+1)) + aparnew(ig)) + aparnew(ig+1)` and
  `apar_p = apargavg(ig) + apargavg(ig+1)`: **no bakdif weighting on apar** (`:5312-5315`),
  whereas the default path applies `bdfac` (GJW 2018 fix, `:5061-5071`).
* The `omprimfac·g_exb` term in `source_coeffs_phip/aparp` lacks `nonmaxw_corr`
  (`:2241-2245`, `:2252-2257`), while `set_source` includes it (`:5101`).
* phiext on ttp points lacks `nonmaxw_corr` (`:5258-5262`).
* `phigavg = (fac_now*phi + fac_next*phinew)*aj0` with `fac_now = fphi*fexp`,
  `fac_next = fphi*(1-fexp)`; bpar part `fac=(fbpar*2)*tz` …`*aj1*vperp2` (`:5162-5185`);
  only computed if `has_phi/has_bpar/has_apar` (different rounding).
* `source = ((phi_m*c_phim + phi_p*c_phip) + apar_m*c_aparm) + apar_p*c_aparp` with
  complex-typed coefficients (`:5309-5327`).

Coefficients (non-LOWFLOW):
```
c_phim  = -((2*vpar)*nm)                                                   (stored complex)
W_opt   = (wstar + ((vpac*code_dt)*wunits)*ufac) - ((((((2.0*omprimfac)*vpac)*code_dt)*wunits)*g_exb)*itor_over_B)/stm
c_phip  = (-((zi*wdrift)*nm)) + zi*W_opt
c_aparm = -(((((zstm*fapar)*vpac)*nm)*(aj0(ig+1)+aj0(ig)))*0.5)
c_aparp = (-((((fapar*D_res)*zstm)*vpac)*nm)) - ((((fapar*stm)*vpac)*zi)*W_opt)
```

### 3.7 Explicit (nonlinear) history bookkeeping

**add_explicit_terms** (`nonlinear_terms.fpp:407-440`), called once per `timeadv` call
(twice per step) with `bd = bkdiff(1)`:
```
if nonlin:
    if istep != 0: add_explicit(gexp_1, gexp_2, gexp_3, phi, apar, bpar, istep, bd, nl=.true.)
else:
    save_dt_cfl(1e8)          # linear: no history, gexp not even allocated
```

**add_explicit** (`nonlinear_terms.fpp:442-612`):
```
if initializing: save_dt_cfl(1e8); max_vel_components = 1e-8; return      # only before init_fields ends
if istep != istep_last:                                                     # once per step
    istep_last = istep
    gexp_3 ← gexp_2 ; gexp_2 ← gexp_1                                       # :505-506
    max_vel_components = epsilon
    gexp_1 ← add_nl(g, phi, apar, bpar)        # bracket at grid points, uses g^n and fields at t^n (= phinew)
    nl_center(gexp_1, bd)                      # → ≈2× upwinded cell value
    [CFL / order-based error → max_vel → maybe reset]                       # :529-583
    if reset: return                                                        # history already shuffled & recomputed
if zip: for ik != 1: g(:,:,iglo)=0 ; gexp_1(:,:,iglo)=0                     # GATED
istep_last = istep
```

**nl_center** (`nonlinear_terms.fpp:677-762`) — pure function of its input (the in-place
ascending loop only ever reads the not-yet-modified `ig+1`):
```
trapped = nlambda > ng2
gin = where(nl_forbid_force_zero & forbid[l, ig], 0, gin)      # whole column incl. ntgrid
for ig in [-ntgrid, ntgrid-1]:
    if is_ttp(ig,l):                      out[:,ig] = gin[:,ig]           # unchanged (no 2×)
    elif trapped and (l > jend(ig) or l > jend(ig+1)):  out[:,ig] = 0
    else:
        out[1,ig] = (1.+bd)*gin[1,ig+1] + (1.-bd)*gin[1,ig]   # rmul + rmul
        out[2,ig] = (1.-bd)*gin[2,ig+1] + (1.+bd)*gin[2,ig]
out[:,ntgrid] = gin[:,ntgrid]                                  # never used by the source
```
(The ttp grid value with the `0.5*code_dt` prefactor is consistent with the ttp equation
being written without the factor 2 that the cell equations carry; see §6.10.)

**State machine for the JAX driver** (carry: `gexp1,gexp2,gexp3, dt, p1, p2, time, dt_old, istep`):
```
start of step n (first timeadv):
    if nonlin and not split: (gexp3, gexp2, gexp1) = (gexp2, gexp1, nl_center(add_nl(g, phi)))
ab = AB(dt, p1, p2)          # identical for both half calls of the step
include_explicit = n > 0 and nonlin and not split
source(call 1, phinew=phi) ; field solve ; source(call 2, phinew) — same gexp, same ab
end of step: p2, p1, dt_old, time = p1, dt, dt, time + dt
```
* Fresh start: `gexp_* = 0`, `p1 = p2 = -1` → step 1 uses order 1 (`ab=[1.0]`), step 2
  order 2, step ≥3 order 3. With `delt_option='check_restart'` p1/p2 come from the restart
  file (`gs2_time.f90:56-62`) and the order can be 3 immediately; `gexp_*` are restored from
  the restart (`gs2_save.fpp:718-810`, `include_explicit_source_in_restart=.true.`).
* Immediate CFL reset (`immediate_reset=.true.`, default): the first `timeadv` returns
  after `add_explicit` (`dist_fn.fpp:4066`), `reset_time_step` rebuilds all dt-dependent
  coefficients (not `D_res`), the step is redone with the **same** `istep`: no second
  shuffle, `gexp_1` reused (NL term is dt independent), AB coefficients recomputed with
  the new `dt` and unchanged `p1,p2`.
* Non-immediate reset: the step completes with the old dt, `update_time` runs, then dt
  changes; the next step uses `this_dt = new dt`, `p1 = old dt`.
* Linear runs: none of this; `include_explicit=.false.`.

### 3.8 Flow shear touching the source

* Parallel-velocity-shear drive: the `-2·omprimfac·vpac·nm·code_dt·wunits·g_exb·itor_over_B/stm`
  term in `W` (§3.4). Active whenever `g_exb ≠ 0`, from step 1, independent of
  `g_exb_start_time/timestep` and of `g_exbfac`.
* Perpendicular shear is implemented by `exb_shear` (`dist_fn.fpp:4183-4600`), called on
  `gnew, phinew, aparnew, bparnew` at the start of `advance_implicit`
  (`fields_implicit.f90:335`). It updates `kx_shift`/`theta0_shift` with
  `gdt = 0.5*(code_dt + code_dt_old)` (`:4283`; `code_dt_old=0` on the first step, so the
  first increment is halved) and shifts `g`/fields by whole kx (or theta0) cells when
  `jump ≠ 0`. **It does not shift `gexp_1..3`**: after a jump the AB history stays
  attached to the old `it` indices. **It does not modify `kperp2`, `aj0`, `aj1`, `wdrift`**
  (no continuous kx in the coefficients). Detailed algorithm: flow-shear spec.

---------------------------------------------------------------------------------------

## 4. Options that change behaviour here

| name | namelist | default | effect on the source | decision |
|---|---|---|---|---|
| `source_option` | `source_knobs` | `'default'` (=`'full'`) | `full`: standard; `phiext_full`: extra `-i ω_d 2 phi_ext sourcefac J0` for `ky=0`; `homogeneous`: only `b*g + a*g` | SUPPORT (all three; trivial masks) |
| `phi_ext` | `source_knobs` | 0.0 | amplitude of phiext term | SUPPORT |
| `t0`, `omega0`, `gamma0`, `source0` | `source_knobs` | 100.0, 0.0, 0.0, 1.0 | `sourcefac` time dependence | SUPPORT |
| `fexpr` | `dist_fn_species_knobs` | 0.4 | time centring `fexp(is)` in `phigavg`, `apargavg`, `a`, `b` | SUPPORT (per species) |
| `bakdif` | `dist_fn_species_knobs` | 0.0 | spatial upwinding `bdfac_p/m`, `a`, `b`; `bkdiff(1)` in `nl_center` | SUPPORT |
| `mult_imp` | `dist_fn_knobs` | F | allows per-species `bakdif` (linear only; forced F if nonlinear, `dist_fn.fpp:1330`) | SUPPORT (export effective `bkdiff`) |
| `driftknob` | `dist_fn_knobs` | 1.0 | scales passing drifts | SUPPORT |
| `tpdriftknob` | `dist_fn_knobs` | -9.9e9 → `driftknob` (`:1286`) | scales trapped+wfb+ttp drifts | SUPPORT |
| `vparknob` | `dist_fn_knobs` | 1.0 | scales `vpa` (hence `vpac`, `vpar`, drifts) | SUPPORT (in exported `vpa`) |
| `mach` | `dist_fn_knobs` | 0.0 | Coriolis drift in `wdrift`, `wdriftttp` (profile-overridable) | SUPPORT |
| `g_exb` | `dist_fn_knobs` | 0.0 | parallel-shear drive in source; remap (other spec) | SUPPORT |
| `omprimfac` | `dist_fn_knobs` | 1.0 | multiplies parallel-shear drive | SUPPORT |
| `g_exbfac` | `dist_fn_knobs` | 1.0 | remap only; **not** in the source | SUPPORT |
| `g_exb_start_time`, `g_exb_start_timestep` | `dist_fn_knobs` | -1, -1 | gate remap only; **not** the source drive | SUPPORT |
| `btor_slab` | `dist_fn_knobs` | 0.0 | overrides `itor_over_B` (`:1369`) | SUPPORT (export effective `itor_over_B`) |
| `opt_source` | `dist_fn_knobs` (+ optimisation override) | F | alternative algorithm (§3.6), different results with apar+bakdif or non-Maxwellian+g_exb | GATE (v1): not used by standard studies; refuse if effective value is T |
| `lf_default`, `lf_decompose` | `dist_fn_knobs` | T, F | LOWFLOW only; no effect in this build | SUPPORT (no-op) |
| `gf_lo_integrate` | `dist_fn_knobs` | F | no effect on source | SUPPORT (no-op here) |
| `boundary_option`, `nonad_zero`, `zero_forbid`, `start_from_previous_solution`, `def_parity`, `boundary_off_grid`, `exponential_boundary` | `dist_fn_knobs` | — | no effect on the source (invert spec) | n/a here |
| `fphi`, `fapar`, `fbpar` | `knobs` | 1.0, 0.0, -1.0 | multiply potentials in `phigavg`, `apargavg`, `apar_m`; `beta=0` forces `fapar=fbpar=0` | SUPPORT (incl. non-unit values) |
| `beta` | `knobs` | 0.0 | through `fapar`/`fbpar` | SUPPORT |
| `delt` | `knobs` | 0.1 | `code_dt` | SUPPORT |
| `delt_option` | `knobs` | `'default'` | `'check_restart'` loads `code_dt`, `code_dt_prev1/2` (AB order) | SUPPORT (export prev dts) |
| `wstar_units` | `knobs` | F | `wunits=1`, `tunits=2/aky` (0 for ky=0: "garbage", `gs2_time.f90:108-120`) | SUPPORT via exported `wunits/tunits` (warn/gate if any `aky==0`) |
| `immediate_reset` | `knobs` | T | reset inside step vs at step end (§3.7) | SUPPORT |
| `eqzip_option` | `knobs` | `'none'` | not applied in source (invert only) | n/a here |
| `nonlinear_mode` | `nonlinear_terms_knobs` | `'default'` (=none) | enables explicit AB terms | SUPPORT |
| `split_nonlinear` | `nonlinear_terms_knobs` | F | NL advanced by separate RK scheme; `include_explicit=.false.` | GATE (different algorithm, other subsystem) |
| `zip` | `nonlinear_terms_knobs` | F | zeroes `g` and `gexp_1` for `ik≠1` inside `add_explicit` (`:599-609`) | GATE (expert/testing) |
| `nl_forbid_force_zero` | `nonlinear_terms_knobs` | T | zero NL input in forbidden region in `nl_center` | SUPPORT |
| `use_order_based_error`, `error_target`, `istep_error_start` | `nonlinear_terms_knobs` | F, 0.1, 30 | dt control using `gexp` history (`estimate_error`, `:620-672`) | GATE (v1; dt-control subsystem) |
| `cfl`, `use_cfl_limit`, `use_2d_cfl` | `nonlinear_terms_knobs` | 0.95, T, T | when resets happen | SUPPORT (dt-control spec) |
| `hyper_option`, `D_hyperres`, `nexp`, `gridnorm` | `hyper_knobs` | `'default'`, -10, 2, T | `D_res` in the apar source term | SUPPORT (export `D_res`, frozen) |
| `delt_adj`, `delt_cushion`, `dt0`, `delt_minimum` | `reinit_knobs` | 2.0, 1.5, delt, 1e-5 | dt values after reset ⇒ coefficient rebuild | SUPPORT |
| `f0type` | `species_parameters_*` | maxwellian | `nonmaxw_corr`, `dlnf0drho` | SUPPORT via export |
| `uprim`, `uprim2` | `species_parameters_*` | 0, 0 | `ufac` | SUPPORT via export |
| `bess_fac` | `species_parameters_*` | 1.0 | `aj0`, `aj1` | SUPPORT via export |
| `gbdriftknob`, `cvdriftknob` | theta_grid knobs | 1, 1 | baked into exported geometry | SUPPORT |
| `trapped_particles`, `bouncefuzz`, `ngauss/npassing` | `le_grids_knobs` | — | masks/`lmax` (exported) | SUPPORT |
| `grid_option` | `kt_grids_knobs` | — | `aky/akx/theta0`, `kwork_filter(1,1)` for box | SUPPORT (box/range/single/specified) |
| `field_option` | `fields_knobs` | `'default'`(implicit) | no effect on source formulas | SUPPORT |
| LOWFLOW compile flag | build | off | complex `wdrift`, `vparterm`, `wdfac`, `wstarfac`, explicit terms in linear runs | GATE (deprecated, not built) |
| `parameter_scan` | `parameter_scan_knobs` | off | runtime changes of `g_exb`, gradients, … | GATE |

---------------------------------------------------------------------------------------

## 5. GPU mapping notes

* **No recurrences** in the source. Every element `(s,e,l,y,x,σ,c)` depends only on
  grid values at `ig=c` and `ig=c+1` of the same `(s,e,l,y,x)` and of the fields at
  `(y,x)`. Implement as whole-array expressions (broadcast) or `vmap` over the leading
  axes; the theta stencil is `arr[..., :-1]` / `arr[..., 1:]`. Same for `nl_center`
  (after zeroing, a pure 2-point stencil) and the `init_wdrift` centring.
* Broadcasting: fields `(Y,X,T)` → `[None,None,None,:,:,(None,):]`; `fexp,bkdiff,zstm,tz,stm`
  `(S,)` → `[:,None,None,None,None,(None,)None]`; `wstar` NumPy `(S,E,Y)` →
  `[:, :, None, :, None, None, None]`; `nonmaxw_corr, ufac` `(S,E)` → `[:, :, None, None, None, None, None]`;
  `D_res` `(Y,X)` → `[None,None,None,:,:,None,None]`; `wunits` `(Y,)`; `itor_over_B[:-1]` `(N,)`;
  `forbid, is_ttp` `(L,T)` → `[None,None,:,None,None,None,:]`.
* Masks (all static per run unless dt changes nothing about them):
  `l <= lmax` (L axis), `is_ttp[l, c] & can_be_ttp[l] & (σ==2) & trapped_grid`,
  `kwork_filter[y,x]`, `aky[y] < eps` (phiext), `aky[y] == 0` (wdrift branch),
  `forbid` (nl_center, wdrift, vpa/vpac), `jend` cell test (nl_center),
  `il_is_wfb/il_is_trapped` (knob choice).
* Output padding: the source has `N = T-1` meaningful cells; store G-shape with the
  last theta entry 0, matching `g_work`.
* Two source evaluations per step (first with `phinew=phi`). Keep both; do not reuse the
  first source plus a linear correction (rounding would differ).
* Time history: keep `gexp1..3` in the `lax.scan` carry; the shuffle is a tuple
  rotation. The AB order (1,2,3) changes only during the first two steps after a fresh
  start → either use a fixed length-3 coefficient vector padded with 0.0 (exact, §3.4) or
  `lax.switch`. dt changes (resets) require rebuilding `vpar, wdrift, wdriftttp, wstar, a,
  b, r, ainv, g_h`, the response matrix and collisions: do them outside the jitted scan
  (host-driven chunks of constant dt).
* Memory: `wdrift`, `vpar`, `vpac`, `a`, `b` are G-shape. Several have redundant axes
  (`vpa`, `vpac`: no `y,x`; `vpar`: no `x`; `vperp2`: no `y,x`). Compressing them is
  exact (values identical by construction) but verify on import with `np.array_equal`
  against the broadcast.
* No irregular structures (linked chains, supercells) appear in the source; they only
  enter in `invert_rhs_1` / linked BCs.
* No MPI/OpenMP reductions in the source → result is independent of processor count and
  layout string (only the iglo ordering changes, which the canonical expansion removes).

---------------------------------------------------------------------------------------

## 6. Pitfalls

1. **FMA contraction.** gfortran `-O3` on aarch64 contracts `a*b+c` (default
   `-ffp-contract=fast`); XLA may or may not, and not at the same places. Expect last-bit
   differences in individual source elements (~1e-16 relative), not bitwise equality.
   The 1e-12 target is reachable; for bit-level debugging rebuild a reference GS2 with
   `F90OPTFLAGS += -ffp-contract=off` (and check what XLA does on the target device).
2. **Literal precision.** All literals are double (`-fdefault-real-8`); `epsilon(0.0)` is
   2.22e-16, not 1.19e-7.
3. **First call of each step has `phinew == phi` exactly.** Still evaluate
   `fexp*phi + (1-fexp)*phinew` (≠ `phi` in general by 1 ulp); `apar_m` is exactly 0.
4. **`vpar(ig,1)` (positive) is used for both σ in `a`, `b`**, while `wdrift(ig,σ)` is
   σ-dependent (Coriolis sign). In the source itself, `vpar(ig,σ)`, `vpac(ig,σ)` are
   signed.
5. **σ swaps:** σ=1 matrix term is `b*g(ig) + a*g(ig+1)`, σ=2 is `a*g(ig) + b*g(ig+1)`;
   `bdfac_p` is `1+bd` for σ=1 and `1-bd` for σ=2.
6. **`wdrift` is centred, `wdriftttp` is not**; `wdrift(ntgrid)` stays a grid value;
   forbidden grid values are 0 *before* centring. The `0.5*(w(ig)+w(ig+1))` reads the
   uncentred `w(ig+1)`.
7. **`D_res` is frozen at the initial dt** (hyper module is not re-initialised on dt
   change, §1.1). The response matrix rebuilt after a reset uses the same stale value, so
   the scheme is self-consistent; replicate, do not "fix".
8. **ttp overwrite.** At ttp points (σ=2) the full `set_source`+matrix value is replaced;
   explicit terms are then re-added once. For `il = nlambda > lmax`, everything except ttp
   points is 0. ttp source uses grid-point `phigavg` and `wstar` only.
9. **`kwork_filter` vs `ieqzip`.** The source honours `kwork_filter` (zeros) but not
   `ieqzip`. `kwork_filter(1,1)` is set true in box mode only after the response matrix
   is built (`fields.f90:376`); export the value seen in the time loop.
10. **Explicit-term normalisation.** `nl_center` returns ≈2× the centred NL value for normal
    cells but the unmodified grid value at ttp points; both are multiplied by
    `0.5*code_dt`. This matches the factor-2 convention of the cell equations
    (`a+b ≈ 2`) vs the ttp equation (`a ≈ 1`); replicate as is.
11. **History is not remapped by `exb_shear`** and is not touched by resets; the shuffle
    happens once per `istep` (guard `istep /= istep_last`), *before* a possible CFL
    reset. On a redo of the same `istep` the old `gexp_1` is reused with new AB
    coefficients.
12. **AB coefficients with equal steps are not the textbook rationals**; compute with the
    variable-step formula in the same order (`/3.`, `/2.` divisions, `q*q`).
13. **`g_exb` drive vs remap gating**: the source term uses raw `g_exb` at all times;
    `g_exbfac`, `g_exb_start_time`, `g_exb_start_timestep` only affect the remap.
14. **Effective vs input knobs**: `tpdriftknob` sentinel, `bkdiff` equalisation,
    `esv`, profile/optimisation overrides of `g_exb`, `mach`, `opt_source`; `itor_over_B`
    overridden by `btor_slab`; `fapar=fbpar=0` when `beta=0`; `fbpar` default is `-1.0`
    (has_bpar false) when `beta≠0` and `fbpar` not set. Export effective values.
15. **`aky == 0.0` exact compare** for the zonal drift branch; `aky < epsilon(0.0)` for
    phiext. With `shat == 0` exactly the zonal branch divides by zero (see §7).
16. **`source(ntgrid)` is 0**, and the source never wraps periodically in theta
    (periodicity / linking is imposed in the invert step).
17. **Signed zeros / NaN**: mixed real×complex and `zi*real` are lowered to real products
    by GCC; generic complex multiplication in JAX gives identical finite values except
    possibly the sign of zero. Forbidden-region values (`aj0=vpa=vperp2=0`) keep the source
    finite there, which matters because `ainv=0` multiplies it.
18. **Complex division** does not occur in the source, but does in `ainv` (and `beta1` in
    the invert step): gfortran uses Smith-style range reduction (`-fcx-fortran-rules`),
    XLA's algorithm may differ in the last bit. Prefer exporting `ainv`/`r` (invert spec).
19. **`wstar` uses `wunits(ik)`**, i.e. `aky/2` (normalisation factor 1/2,
    `gs2_time.f90:91-107`); `vpar` uses `tunits(ik)` (=1).
20. **ufac uses `energy**1.5`** (`pow`) and `sqrt(pi)`; export rather than recompute.
21. **The source routine is called with `istep=0` during response-matrix construction**;
    `include_explicit` must be false there even in nonlinear runs.

---------------------------------------------------------------------------------------

## 7. Open questions

1. FMA: which operations does the shipped gfortran build actually contract (inspect
   `objs/dist_fn.o` disassembly for `fmadd/fmsub` in `get_source_term`), and does XLA
   on the target GPU contract the same expressions? Decide whether the reference should
   be a `-ffp-contract=off` build.
2. `shat == 0` exactly with a `ky=0` mode: `wdrift_func`/`wcoriolis_func` divide by
   `shat` (`dist_fn.fpp:1519`, `:1550-1551`). Is `shat` ever exactly 0 in a supported run (zero-shear
   slab), and what does the exported `wdrift` contain (Inf/NaN)? Propose GATE on
   `shat == 0` with zonal modes until clarified.
3. `D_res` frozen at initial dt: intended GS2 behaviour or bug? (We replicate.) Same
   question for any other `hyper` quantity (hyper spec).
4. `opt_source` differences (§3.6) look like an unported 2018 bakdif fix and a missing
   `nonmaxw_corr`; confirm with GS2 developers before ever supporting it.
5. Whether to port the complex division used in `ainv` or have the exporter dump
   coefficient sets for every reachable dt (with `delt_adj = 2` the reachable steps are
   `dt0·2^-k`, capped at `code_dt_max`): decision for the invert/dt-control specs; it also
   determines whether `a`, `b`, `wdrift`, `vpar` must be ported (they are simple and
   recommended as RECOMPUTE regardless).
6. `nl_center` sits between the NL bracket and the source; confirm which spec owns it (it is
   fully specified here).
7. The `exb_shear` non-remapping of `gexp_*` (and of `kx`-continuous coefficients) is GS2
   behaviour; confirm it is to be replicated bit-for-bit (we assume yes).
8. `wstar_units = .true.` with `aky = 0` gives `tunits = 0` (GS2 prints a warning): gate or
   replicate?
9. `nonlinear_terms::istep_last`, `initializing` and the CFL state are private; the exporter
   cannot dump them directly. The JAX driver re-derives them (shuffle at every new step);
   confirm no exported run starts mid-step.

---------------------------------------------------------------------------------------

## 8. Verification plan and export checklist

Export (netCDF, double; complex as `_r/_i` pairs like `gs2_save.fpp:310-337`; logicals as
int8), after `init(full)` and before the first `advance`:

* Grids/geometry: `ntgrid, nspec, negrid, nlambda, naky, ntheta0, ng2, lmax`, `theta`,
  `bmag, gradpar, delthet, gbdrift, gbdrift0, cvdrift, cvdrift0, cdrift, cdrift0,
  itor_over_B, shat`, `aky, akx, theta0, kperp2`, `energy, al`, `forbid, is_ttp,
  can_be_ttp, is_bounce_point, jend`, `kwork_filter`.
* Species: `zstm, tz, zt, stm, smz`, `nonmaxw_corr, dlnf0drho`, `fexp, bkdiff`.
* Coefficients: `aj0, aj1, vperp2, vpa, vpac, vpar, wdrift, wdriftttp, wstar, ufac,
  D_res, a, b, r, ainv` (+ `source_coeffs_*` if ever supported).
* Scalars: all knobs of §4 (effective), `code_dt, code_dt_prev1, code_dt_prev2,
  code_time, code_dt_old, wunits, tunits`.
* State: `g, gnew, phi*, apar*, bpar*`, `gexp_1..3`.

Unit tests (oracle = `dump_current_source_term`, `dist_fn.fpp:7850`):

1. Linear ES, `phinew == phi` (first-call condition) and with a random `phinew`
   perturbation (second-call condition); compare `g_work` elementwise (rel. 1e-14).
2. Same with `bakdif ≠ 0`, `fexpr ≠ 0.5`, trapped particles (ttp points, `il = nlambda`).
3. EM (`beta>0`, `fapar=1`, `fbpar=1`) with `aparnew ≠ apar`; with `D_hyperres > 0`.
4. Nonlinear with non-zero `gexp_1..3` and `istep = 1, 2, ≥3` (AB orders 1–3), unequal
   `code_dt_prev1/2`.
5. `g_exb ≠ 0`, `mach ≠ 0`, `phiext_full` with `ky=0` and `code_time` on both sides of `t0`.
6. Recompute `vpar, wdrift, wdriftttp, wstar, a, b` in JAX from exported ingredients at
   the initial dt and at `dt/2` (Fortran reference: rerun the exporter with `delt/2`);
   require exact or ≤1-ulp agreement.
