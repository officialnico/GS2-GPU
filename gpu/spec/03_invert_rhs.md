# 03 — The implicit parallel solve (`invert_rhs`)

Status: contract for implementers of the JAX (float64/complex128) port. Read-only analysis of
GS2 8.2.1 (`gpu-port` branch, commit `4d8c94bcf`). Every statement cites `file:line`
(paths relative to `/Users/nico/Workspace/gs2`; `src/` unless stated).

Scope (this document):

* `dist_fn.fpp`: `invert_rhs` (6152-6194), `invert_rhs_1` (5440-5707),
  `apply_parallel_boundary_conditions` (5710-5849), `self_periodic` (5864-5906),
  `apply_linked_boundary_conditions` (5919-6149; this is what older comments call
  `invert_rhs_linked`, e.g. `dist_fn.fpp:2780`, `:5548`), `enforce_parity` (4796-4853) and
  `init_enforce_parity`/`get_parity_conn` (4603-4794), `init_invert_rhs` (1826-1929: `a, b, r,
  ainv`), `init_homogeneous_g` (1937-2125: `g_h`, `g_adj` h-part), `init_bc` (2128-2156),
  `init_connected_bc` (3513-3602) with `compute_jshift0` (2401-2436),
  `compute_itleft_and_itright` (2441-2513), `compute_connections` (2518-2565), `set_save_h`
  (2570-2604), `count_links*` (2608-2672), `setup_connected_bc_redistribute` (2759-3022),
  `setup_connected_bc_redistribute_mixed_wfb` (3028-3224),
  `setup_pass_incoming_boundary_to_connections` (3229-3508), `get_left/right_connection`
  (3837-3886), `get_leftmost_it`/`get_rightmost_it` (3889-3920), `find_left/rightmost_link`
  (7902-7964), `calculate_supercell_labels` (2355-2397).
* `dist_fn_arrays.fpp`: `get_adjust` (178-188), `set_h_zero` (246-294).
* `le_grids.f90`: `forbid`, `is_bounce_point`, `is_lower/upper_bounce_point`, `can_be_ttp`,
  `is_ttp`, `jend`, `il_is_wfb/passing/trapped`, `grid_has_trapped_particles`,
  `passing_wfb/trapped_wfb/mixed_wfb`, `is_passing_hybrid_electron`.
* `externals/utils/redistribute.f90`: semantics of `fill`/`scatter` used for the linked BCs.

Out of scope (owned by other specs, only interfaced here): the source `g_work`
(`02_source_term.md`, which also owns the formulas for `vpar`, `wdrift`, `wdriftttp`, `a`,
`b`), the field solve / response matrices, `hyper_diff`, collisions, `exb_shear` remap, the
time loop (`01_time_loop.md`).

----------------------------------------------------------------------------------------------

## 0. Conventions and build facts

### 0.1 Index conventions

* `ig ∈ [-ntgrid, ntgrid]`, `T = 2*ntgrid+1` grid points; Python index `j = ig + ntgrid`.
* `isgn`/σ: σ=1 is `v_|| > 0` (sweeps left→right, incoming boundary at `ig=-ntgrid`), σ=2 is
  `v_|| < 0` (sweeps right→left, incoming boundary at `ig=+ntgrid`).
* `iglo` is always expanded to `(it, ik, il, ie, is)` with `it_idx/ik_idx/il_idx/ie_idx/is_idx`.
  Distribution-like arrays (`gnew`, `g`, `g_work`(source), `g_h`, `a`, `b`, `r`, `ainv`,
  `wdrift`, `wdriftttp`, `vpar`, `vpac`), Fortran `(ig, isgn, iglo)`, are exported as Fortran
  `(ig, isgn, it, ik, il, ie, is)` → **NumPy `(nspec, negrid, nlambda, naky, ntheta0, 2, T)`**
  ("G-shape", axes `S,E,L,Y,X,σ,T`), identical to `02_source_term.md` §0.1.
* `(ig, iglo)` arrays (`vperp2`, `aj0`, `aj1`) → NumPy `(S,E,L,Y,X,T)`.
* Fields `phinew`, `bparnew` Fortran `(ig, it, ik)` → NumPy `(Y, X, T)`.
* `(it, ik)` arrays (`l_links`, `r_links`, `itleft`, `itright`, `kwork_filter`, `ieqzip`,
  `supercell_labels`, `n_cells`) Fortran `(ntheta0, naky)` → NumPy `(Y, X)`.
  `n_links(2, ntheta0, naky)` → NumPy `(Y, X, 2)`. `save_h(2, iglo)` → expanded Fortran
  `(2, it, ik, il, ie, is)` → NumPy `(S,E,L,Y,X,2)`.
* `(ig, il)` le_grids flags (`forbid`, `is_ttp`, `is_bounce_point`) Fortran
  `(-ntgrid:ntgrid, nlambda)` → NumPy `(L, T)`. `can_be_ttp(nlambda)` → `(L,)`.
* **Integer *values*** stored in `itleft/itright` are Fortran 1-based `it` indices, `-1` meaning
  "no connection"; keep them 1-based in the export and convert in Python (`it0 = v-1`).
* Pitch-angle classes (`le_grids.f90:443-491`): passing `il <= ng2`, wfb `il == ng2+1`,
  trapped `il > ng2+1` (wfb is *neither* passing nor trapped for these predicates). The wfb
  exists only if `nlambda > ng2` (`grid_has_trapped_particles()`).

### 0.2 Build facts that matter for round-off (verified on this machine)

* `bin/gs2 --build-config`: GCC 15.2.0, `-O3 -mcpu=apple-m1 -fdefault-real-8
  -fdefault-double-8`, MPI on, **OpenMP off**, **LOWFLOW off**, `SPFUNC=_SPF200X_`
  (`j0`/`j1` are `bessel_j0`/`bessel_j1` from the system libm, `externals/utils/spfunc.fpp:741-791`;
  `j1(x)` returns `J1(x)/x`, and `0.5` for `x == 0` exactly).
* All literals and `real` are binary64; `epsilon(0.0) = 2.220446049250313e-16`.
* **Complex division is Smith's algorithm (GCC "wide" method, inlined).** Verified: a test
  program built with the same flags plus `-ffp-contract=off` reproduces, bit for bit on 20 000
  random operands, the formula in §3.0 (`cdiv`); no `__divdc3` calls exist in `objs/dist_fn.o`.
* **FMA contraction is ON in the shipped build.** `otool -tv objs/dist_fn.o` shows 116
  `fmadd/fmsub/fnmsub` instructions inside `___dist_fn_MOD_invert_rhs.constprop.0` (which
  inlines `invert_rhs_1`), 90 in `apply_linked_boundary_conditions`, 4 in
  `init_homogeneous_g`, 106 in `___dist_fn_MOD_init_dist_fn_level_3` (which inlines
  `init_vpar`/`init_wdrift`/`init_invert_rhs`). With default flags gfortran's complex divide
  differs from plain Smith in ~19 % of random cases and complex multiply from the naive formula
  in ~29 %. See §6.1 for the strategy.
* Complex multiplication without contraction is the naive `(ar*br - ai*bi, ar*bi + ai*br)`
  (verified bitwise, 20 000 cases). Mixed real×complex and `zi*real` are lowered component-wise
  by GCC's complex-lowering pass (no products with the known-zero part); the component formulas
  of §3.3, §3.8 and §3.10 were checked bitwise against a `-ffp-contract=off` build (§A.2).
* Fortran evaluation order: equal-precedence binary operators associate left-to-right; unary
  minus binds like binary `+/-`, so `-a*b + c*d` is `(-(a*b)) + (c*d)`; gfortran does not
  reassociate without `-ffast-math`. All pseudo-code below is fully parenthesised. `-x + y` and
  `y - x` are bitwise identical in IEEE arithmetic (up to NaN payloads), so they are used
  interchangeably.
* MPI never enters the arithmetic of this subsystem: the only inter-process operations are
  `redistribute::fill`/`scatter`, which are pure copies (`externals/utils/redistribute.f90:7685-7739`,
  `6425-6485`, `6787-6845`). Results are independent of `nproc` and of the `layout` string.

----------------------------------------------------------------------------------------------

## 1. Control flow / call graph

### 1.1 Initialisation

`gs2_init` levels (`gs2_init_level_list.inc`): `le_grids`=18, `dist_fn_parameters`=21,
`dist_fn_arrays`=26, `dist_fn_level_1`=27, `dist_fn_level_2`=28, `override_timestep`=29,
`collisions`=30, `dist_fn_level_3`=31, `fields_level_2`=32.

```
init_le_grids                                   le_grids.f90:242
  read_parameters                               le_grids.f90:352-428  ng2 (:400-405), bouncefuzz (:370),
                                                                        passing_wfb/trapped_wfb/mixed_wfb (:413-427)
  set_grids (proc0) → lgridset                  le_grids.f90:1861, 2120
     calculate_forbidden_region                 le_grids.f90:2715   forbid
     calculate_bounce_points                    le_grids.f90:2746   is_bounce_point
     calculate_jend                             le_grids.f90:2871   jend
     calculate_ittp                             le_grids.f90:2907   can_be_ttp, is_ttp
  broadcast_results                             le_grids.f90:299-349
init_dist_fn_parameters → read_parameters       dist_fn.fpp:849, 1174
     boundary_option → boundary_option_switch   dist_fn.fpp:1290-1310 (default resolution, §4)
     tpdriftknob sentinel                       dist_fn.fpp:1286
     bkdiff forced equal unless mult_imp        dist_fn.fpp:1330-1344
init_dist_fn_arrays → allocate_arrays           dist_fn.fpp:885, 3923   g, gnew, g_work, g_h allocated+zeroed (3943-3949)
init_dist_fn_level_1                            dist_fn.fpp:928
  init_vperp2, init_vpa_vpac                    dist_fn.fpp:1560, 1592
  init_bc                                       dist_fn.fpp:2128
     l_links=r_links=n_links=0                  :2132-2137
     linked:  init_connected_bc                 :2141 → dist_fn.fpp:3513
                compute_jshift0                 :3532 → 2401
                compute_itleft_and_itright      :3535 → 2441
                compute_connections             :3538 → 2518
                set_save_h                      :3541 → 2570
                count_links                     :3543 → 2645  (l_links, r_links, n_links)
                n_links_max (+3 if mixed_wfb)   :3544-3547
                setup_fields_implicit_class_arrays :3551 → 2676 (M_class, N_class; fields spec)
                setup_connected_bc_redistribute :3555 → 2759  (links_p, links_h, no_connections)
                [return if no_connections]      :3568-3571
                setup_..._mixed_wfb (mixed_wfb) :3576 → 3028  (wfb_p, wfb_h)
                setup_pass_incoming_... (sfps)  :3588 → 3229  (incoming_links; only if start_from_previous_solution)
                allocate g_adj(n_links_max,2,g_lo), zero  :3599-3600
              init_enforce_parity (def_parity)  :2142 → 4603
     other:   i_class=1, M_class=naky*ntheta0, N_class=1   :2143-2150
     calculate_supercell_labels                 :2155 → 2355
init_dist_fn_level_2 (bessel, fieldeq)          dist_fn.fpp:939
init_dist_fn_level_3                            dist_fn.fpp:960        ← redone on every dt change
  init_vpar      (dt)                           :982  → 1681
  init_wdrift    (dt)                           :985  → 1405
  init_wstar     (dt)                           :988
  init_collisions                               :991
  init_invert_rhs (dt)                          :1003 → 1826   a, b, r, ainv
     init_homogeneous_g(g_h)                    :1928 → 1937   g_h; then fill(links_h/wfb_h, g_h, g_adj) :2120-2123
  init_hyper, init_source_term                  :1006, :1009
init_fields                                     fields.f90:353
  fields_init_response   (calls timeadv(...,istep=0) per unit field; see §1.3)
  if box: kwork_filter(1,1) = .true.            fields.f90:376
```

### 1.2 Per time step

```
gs2_main.fpp:730  advance(istep)                       fields.f90:586
  advance_implicit (fields_implicit.f90:315) | advance_local (fields_local.f90:4244) | advance_gf_local (fields_gf_local.f90:4239)
    [exb_shear remap of gnew/fields]                   (does not touch any invert_rhs coefficient)
    g = gnew ; phi = phinew ; ...                       fields_implicit.f90:337-338
    timeadv(phi,...,phinew,..., istep)                  fields_implicit.f90:341   ← call #1 (phinew == phi)
    getfield; phinew = phinew + phi ...                 fields_implicit.f90:352-356
    timeadv(phi,...,phinew,..., istep, diagnostics)     fields_implicit.f90:365   ← call #2 (corrected fields)
timeadv                                                dist_fn.fpp:4014
  explicit NL terms                                    :4046-4062
  if (reset .and. immediate_reset) return              :4066
  invert_rhs(phi, apar, bpar, phinew, aparnew, bparnew, istep)   :4068
      time/sourcefac                                   :6166-6176   (02_source_term §3.2)
      zero_array(g_work)                               :6179
      get_source_term_on_local_domain(g_work, ...)     :6181        (02_source_term)
      invert_rhs_1(phinew, bparnew, source=g_work)     :6184 → 5440
          per iglo: masks, set_h_zero, flags, init gnew, BCs, σ=2 sweep, bounce source,
          σ=1 sweep, self_periodic, trapped β, trapped-wfb β, zero_forbid      (§3.7)
          if def_parity: enforce_parity(parity_redist) :5706
      if linked: apply_linked_boundary_conditions(gnew, g_h, g_adj, phinew, bparnew)  :6188-6191 → 5919
          [start_from_previous_solution: fill(incoming_links) + difference]   :5951-5982
          fill(links_p, gnew, g_adj); if mixed_wfb fill(wfb_p, gnew, g_adj)   :5984-5985
          per iglo: add b0*g_h                                                 :5996-6146
  hyper_diff (if istep>0 or hyper_in_initialisation)   :4073-4080
  vspace_derivatives (if .not. split_collisions)       :4083-4087
  if def_parity: enforce_parity(parity_redist)         :4091
```

`invert_rhs` is therefore executed **twice per time step** with the same `istep`, the same
`g` and `phi/apar/bpar`, but different `phinew/bparnew` (the incoming BC values and the
self-periodic/linked `dadj` corrections use `phinew/bparnew`).

### 1.3 Response-matrix construction

`fields_implicit.f90:669` (`init_response_row`), `fields_local.f90:2184`,
`fields_gf_local.f90:1955` call `timeadv(..., 0)` (istep=0) with unit-impulse fields. The same
`invert_rhs` code runs; only the source differs (`include_explicit=.false.`). These matrices
depend on dt through `r/ainv/g_h` and must be rebuilt after a dt change (fields spec). Note
`kwork_filter` state during such rebuilds (§6.6).

### 1.4 dt change (reset)

`gs2_reinit.f90:104-225` (`reset_time_step`): `init(current_init, override_timestep)` (:166)
tears down levels ≥ 29, sets the new `code_dt`, then `init(current_init, full)` (:216) re-runs
`init_dist_fn_level_3` → `init_vpar`, `init_wdrift`, `init_invert_rhs` (+`init_homogeneous_g`
and the `links_h`/`wfb_h` fill into `g_adj`). Level-1 topology (`itleft`, `l_links`,
`connections`, `save_h`, redistributes, `g_adj` allocation) is **not** recomputed.
`finish_dist_fn_level_3` (:1140-1154) deallocates `vpar, wdrift, wdriftttp, wstar`; `a,b,r,ainv`
are re-zeroed (:1851-1852) and recomputed; `g_h` is re-zeroed inside `init_homogeneous_g`
(:1959).

----------------------------------------------------------------------------------------------

## 2. Data inventory

Decision legend: **EXPORT** (dump from Fortran), **RECOMPUTE** (port the formula, given in §3),
**BOTH** (export at the initial dt for validation *and* port the formula so it can be rebuilt on
dt change). "dt" = depends on `code_dt`.

### 2.1 Dynamic inputs/outputs of one `invert_rhs` call

| name | module | Fortran shape | NumPy | role | decision |
|---|---|---|---|---|---|
| `g_work` (= `source`) | `dist_fn_arrays` | `complex(-ntgrid:ntgrid,2,g_lo)` | G-shape | input from `02_source_term`; `ig=ntgrid` always 0; **modified in place** at `(lb-1, σ=1)` (`:5619`) | produced by source spec |
| `gnew` | `dist_fn_arrays` | same | G-shape | output; *pass-through* (unchanged) at `kwork_filter`/`ieqzip` points | state |
| `g` | `dist_fn_arrays` | same | G-shape | read only if `start_from_previous_solution` (`:5572`) | state |
| `phinew`, `bparnew` | `fields_arrays` | `complex(-ntgrid:ntgrid,ntheta0,naky)` | `(Y,X,T)` | BC values (`nonad_zero`), `dadj`, hybrid `set_h_zero` | state |
| `g_adj` | `dist_fn` (private) | `complex(n_links_max,2,g_lo)` | — | scratch: p-slots filled every call, h-slots at init | RECOMPUTE as gathers (§3.11); never export |

### 2.2 Coefficients (dt-dependent)

| name | module | Fortran shape | NumPy | computed in | dt | decision |
|---|---|---|---|---|---|---|
| `ainv` | `dist_fn` (private) | `complex(-ntgrid:ntgrid,2,g_lo)`; `ig=ntgrid` entries 0 | G-shape | `init_invert_rhs` `dist_fn.fpp:1885-1887`, zeroing `:1901-1904`, lower-bounce `+1` `:1910-1914`, ttp `:1917-1918` | yes | BOTH (§3.3) |
| `r` | `dist_fn` (private) | same | G-shape | `:1888-1891`, `:1901-1904`, `:1920` | yes | BOTH (§3.3) |
| `a`, `b` | `dist_fn` (private) | same | G-shape | `:1892-1897`, ttp `a` `:1919` (not used by invert_rhs_1; consumed by the source) | yes | BOTH (formula owned by `02_source_term.md` §3.1, restated §3.3) |
| `g_h` | `dist_fn` (private) | `complex(-ntgrid:ntgrid,2,g_lo)` | G-shape | `init_homogeneous_g` `:1937-2125` | yes (through `r`) | BOTH (§3.4) |
| `wdrift` | `dist_fn` (private) | `real(-ntgrid:ntgrid,2,g_lo)` | G-shape | `init_wdrift` `:1405-1498` | yes | BOTH (owned by 02 §3.1) |
| `wdriftttp` | `dist_fn` (private) | same | G-shape | `:1447-1480` (non-zero only at `is_ttp`) | yes | BOTH (owned by 02) |
| `vpar` | `dist_fn` (private) | same | G-shape (`vpar(:,2) == -vpar(:,1)`) | `init_vpar` `:1681-1729` | yes | BOTH (owned by 02). invert uses **only `vpar(ig,1,iglo)`** (`:1882`) |

### 2.3 Static (dt-independent) per-point data

| name | module | Fortran shape | NumPy | computed in | decision |
|---|---|---|---|---|---|
| `forbid` | `le_grids` | `logical(-ntgrid:ntgrid,nlambda)` | `(L,T)` | `le_grids.f90:2715-2738` | EXPORT |
| `is_bounce_point` | `le_grids` | same | `(L,T)` | `le_grids.f90:2746-2771` (only `il >= ng2+2`) | EXPORT |
| `is_lower_bounce_point(ig,il)` | `le_grids` (elemental fn) | — | `(L,T)` | `le_grids.f90:2775-2819` | RECOMPUTE from `is_bounce_point`,`forbid`,`trapped_wfb` (§3.1); EXPORT evaluated array for check |
| `is_upper_bounce_point(ig,il)` | `le_grids` (elemental fn) | — | `(L,T)` | `le_grids.f90:2823-2867` | RECOMPUTE (§3.1); EXPORT for check |
| `is_ttp` | `le_grids` | `logical(-ntgrid:ntgrid,nlambda)` | `(L,T)` | `le_grids.f90:2907-2946` | EXPORT |
| `can_be_ttp` | `le_grids` | `logical(nlambda)` | `(L,)` | same | EXPORT (only gates `wdriftttp`, upstream) |
| `ng2, nlambda, lmax, ntgrid, ntheta0, naky, negrid, nspec` | various | integers | scalars | | EXPORT |
| `passing_wfb, trapped_wfb, mixed_wfb` | `le_grids` | logical | scalars | `le_grids.f90:413-427` | EXPORT |
| `aky` | `kt_grids` | `real(naky)` | `(Y,)` | | EXPORT (exact `== 0.0` tests) |
| `kwork_filter` | `kt_grids` | `logical(ntheta0,naky)` | `(Y,X)` | `.false.`, then `(1,1)=.true.` in box (`fields.f90:376`) | EXPORT (value during the time loop) |
| `ieqzip` | `run_parameters` | `logical(ntheta0,naky)` | `(Y,X)` | `run_parameters.f90:259-272` | EXPORT (all-false unless `eqzip_option` set → GATE) |
| species type / `is_passing_hybrid_electron(is,ik,il)` | `le_grids.f90:4490-4499` | — | `(S,Y,L)` bool | `is_hybrid_electron_species(spec(is)) .and. aky(ik)/=0 .and. (il<=ng2+1)` | RECOMPUTE (export species type) |
| `vperp2` | `dist_fn_arrays` | `real(-ntgrid:ntgrid,g_lo)` | `(S,E,L,Y,X,T)` | `init_vperp2` `:1584-1587` | EXPORT (shared with 02) |
| `aj0`, `aj1` | `dist_fn_arrays` | same | same | `init_bessel` `:1784-1793` | EXPORT (libm Bessel) |
| `nonmaxw_corr` | `species` | `real(negrid,nspec)` | `(S,E)` | `species.f90:571-592` | EXPORT |
| `spec%zt, tz, zstm, smz, bess_fac, type` | `species` | per species | `(S,)` | | EXPORT |
| `fexp`, `bkdiff` | `dist_fn` (private) | `real(nspec)` | `(S,)` | `dist_fn.fpp:1397-1398`, forced equal `:1333-1344` | EXPORT (effective) |
| `fphi`, `fbpar` | `run_parameters` | real | scalar | `run_parameters.f90:314-352` (**`fbpar` default −1.0**, 0 if `beta==0`) | EXPORT (effective) |
| `delthet` | `theta_grid` | `real(-ntgrid:ntgrid)` | `(T,)` | | EXPORT (only for GATED `exponential_boundary`) |
| `wfb` (knob) | `dist_fn` | real | scalar | `dist_fn.fpp:1275`, default 1.0 (`:418`) | EXPORT |
| flags `boundary_option_switch, nonad_zero, zero_forbid, start_from_previous_solution, boundary_off_grid, exponential_boundary, exponential_boundary_factor, def_parity, even` | `dist_fn` | int/logical/real | scalars | `read_parameters` | EXPORT |

### 2.4 Linked-boundary topology (only `boundary_option_switch == boundary_option_linked`)

| name | module | Fortran shape | NumPy | computed in | decision |
|---|---|---|---|---|---|
| `itleft`, `itright` | `dist_fn` (public) | `integer(ntheta0,naky)`, values 1-based it or −1 | `(Y,X)` | `compute_itleft_and_itright` `:2441-2513` | EXPORT (formula §3.5 for reference) |
| `l_links`, `r_links` | `dist_fn` (public) | `integer(ntheta0,naky)` | `(Y,X)` | `count_links_one_side` `:2608-2640` | EXPORT + RECOMPUTE-check |
| `n_links` | `dist_fn` (private) | `integer(2,ntheta0,naky)` | `(Y,X,2)` | `count_links` `:2659-2671`: `2*links-1` or 0 | RECOMPUTE |
| `connections(iglo)%neighbor` | `dist_fn` (private) | `logical(g_lo)` | `(S,E,L,Y,X)`, depends on `(it,ik,il)` only | `compute_connections` `:2518-2565` | RECOMPUTE (§3.5); EXPORT for check |
| `save_h` | `dist_fn` (private) | `logical(2,g_lo)` | `(S,E,L,Y,X,2)`, depends on `(it,ik,il)` | `set_save_h` `:2570-2604` | RECOMPUTE (§3.5); EXPORT for check |
| `no_connections` | `dist_fn` (private) | logical (global) | scalar | `setup_connected_bc_redistribute` `:2870-2876` | EXPORT |
| `supercell_labels`, `n_cells` | `dist_fn` (public) | `integer(ntheta0,naky)` | `(Y,X)` | `calculate_supercell_labels` `:2355-2397` | EXPORT (convenience; chains can be rebuilt from `itright`) |
| mixed-wfb end factors `aj0l, aj1l, aj0r, aj1r, vperp2left, vperp2right` | locals of `apply_linked_boundary_conditions` | — | `(S,E,Y,X)` (il = ng2+1) | `:6019-6028` (**own argument association**, differs from `init_bessel`) | EXPORT (new export arrays, §3.11.4) |
| `links_p, links_h, wfb_p, wfb_h, incoming_links, parity_redist` | `dist_fn` (redist_type) | — | — | §3.5 | RECOMPUTE as static gather index tables |

**Accessibility warning.** `a, b, r, ainv, g_h, wdrift, wdriftttp, vpar, n_links, connections,
save_h, g_adj, no_connections, fexp, bkdiff` are **private** to module `dist_fn` (public list
`dist_fn.fpp:21-95`). The export program needs an accessor added to `dist_fn` (or must
re-evaluate the formulas — not recommended). Public helpers that can be reused:
`dump_homogeneous_solution` (`:7878-7892`, writes `g_h`-like data via `gs2_save_for_restart`),
`dump_current_source_term` (`:7850-7876`), `wdrift_func` (public function), `l_links`,
`r_links`, `itleft`, `itright`, `get_leftmost_it`, `supercell_labels`, `n_cells`.

----------------------------------------------------------------------------------------------

## 3. Exact computations

Notation per element: `s=is, e=ie, l=il, y=ik, x=it, σ=isgn`. All loops over `iglo` are over
independent elements (no cross-iglo data flow inside `invert_rhs_1`); only the loops over `ig`
are recurrences. Loop *order over iglo* is irrelevant for values.

### 3.0 Arithmetic helpers

Same names as `02_source_term.md` §3.0 (`rmul`, `imul`, `cmul`), plus `cdiv`. A "strict" variant
wraps each product in `lax.optimization_barrier` to forbid FMA contraction by XLA (§6.1); this
reproduces a `-ffp-contract=off` gfortran build bit for bit on XLA:CPU (verified, §A.2).

```python
from jax import lax
ob = lax.optimization_barrier            # strict mode; use identity for fast mode

def cmul(a, b):                          # Fortran complex*complex (both general)
    return lax.complex(ob(a.real*b.real) - ob(a.imag*b.imag),
                       ob(a.real*b.imag) + ob(a.imag*b.real))

def cdiv(a, b):                          # Fortran complex/complex: GCC wide (Smith) method
    ar, ai, br, bi = a.real, a.imag, b.real, b.imag
    t = jnp.abs(br) < jnp.abs(bi)        # strict '<': ties take the FALSE branch
    r1 = br/bi; d1 = ob(br*r1) + bi; x1 = ob(ar*r1) + ai; y1 = ob(ai*r1) - ar
    r2 = bi/br; d2 = ob(bi*r2) + br; x2 = ob(ai*r2) + ar; y2 = ai - ob(ar*r2)
    return lax.complex(jnp.where(t, x1/d1, x2/d2), jnp.where(t, y1/d1, y2/d2))

ONE = lax.complex(1.0, 0.0)              # Fortran 1.0 promoted to complex
```

* `1.0/z` (real numerator) is `cdiv(ONE, z)`: GCC keeps the general expansion with `ai = 0.0`
  (it may fold `1.0*ratio → ratio`, which is exact), so use `cdiv` with `ar=1, ai=0` literally;
  this also reproduces signs of zero.
* real ± complex touches only the real part; `1.0 - z` = `(1.0 - z.re, -z.im)`.
* real×complex and complex×real are component-wise (`rmul`); integer literals (`2*`, `3*`) are
  converted to double first.
* Unused `jnp.where` branches may produce inf/nan when `br` or `bi` is 0; they are discarded.

### 3.1 le_grids flags (definitions, for RECOMPUTE-checks)

```
forbid[l,ig]          = (1.0 - al(l)*bmag(ig)) < -bouncefuzz                    le_grids.f90:2727
is_bounce_point[l,ig] = (l >= ng2+2) and abs(1.0 - al(l)*bmag(ig)) <= bouncefuzz  :2760-2768 (stmt :2767)
is_lower_bounce_point(ig,l):                                                       :2775-2819
    if not is_bounce_point[l,ig]: False
    elif il_is_wfb(l) and not trapped_wfb: False      # never reached: wfb has no bounce points
    elif ig == -ntgrid: True
    elif ig ==  ntgrid: False
    else: forbid[l,ig-1]          # (forbid(ig-1) .or. (forbid(ig-1) .and. forbid(ig+1))) == forbid(ig-1)
is_upper_bounce_point(ig,l):                                                       :2823-2867
    same, with ig==-ntgrid → False, ig==ntgrid → True, else forbid[l,ig+1]
ittp(ig) = min{ l in [ng2+1, nlambda] : forbid[l,ig-1] and forbid[l,ig+1] and not forbid[l,ig] }
           for ig in [-ntgrid+1, ntgrid-1]; nlambda+1 otherwise / if none                 :2919-2930
is_ttp[l,ig] = (l >= ng2+1) and l >= ittp(ig) and not forbid[l,ig]                        :2940-2944
can_be_ttp[l] = (l >= ng2+1) and any_ig(l >= ittp(ig))                                    :2934-2936
grid_has_trapped_particles = nlambda > ng2                                                :443-446
```
`bouncefuzz` default `10*epsilon(0.0)` (`le_grids.f90:131`). All are `False` when
`nlambda == ng2`. A ttp point that is also flagged in `is_bounce_point` (the usual case: the
pitch angle reflecting exactly at that local `bmag` minimum) is simultaneously a lower and an
upper bounce point; a ttp point need not be a bounce point in general (the masks handle both).
The wfb is never forbidden, never a bounce point, never ttp. Export the arrays; the formulas are only
needed to understand the masks.

### 3.2 Upstream dt-dependent inputs (owned by `02_source_term.md` §3.1)

`vpar`, `wdrift` (centred), `wdriftttp` (grid value at ttp points, not centred). For dt
changes the JAX side rebuilds them with the formulas in 02 §3.1 (`code_dt` enters inside the
product chain, e.g. `(((zstm*tunits)*code_dt)*0.5)/delthet...` and `(Q*code_dt)*wunits`), so
export the dt-free inputs listed there (`vpac`, `vpa`, `vperp2`, geometry, `wunits`, `tunits`,
knobs). The invert step uses `wd = wdrift(ig,σ,iglo)`, `wdttp = wdriftttp(ig,σ,iglo)` and
`vp = vpar(ig,1,iglo)` (**σ=1 for both σ**, `:1882`) at `ig ∈ [-ntgrid, ntgrid-1]`.

### 3.3 `init_invert_rhs` — `a, b, r, ainv` (`dist_fn.fpp:1826-1929`)

Element-wise for every iglo, σ ∈ {1,2}, `ig ∈ [-ntgrid, ntgrid-1]` (loop nest iglo → isgn → ig,
`:1862-1925`; order irrelevant). Entries `ig = ntgrid` stay `0` (`zero_array`, `:1851-1852`).

Fortran (verbatim, `:1885-1897`):
```fortran
ainv = 1.0/(1.0 + bd + (1.0-fexp(is))*spec(is)%tz*(zi*wd*(1.0+bd) + 2.0*vp))
r    = (1.0 - bd + (1.0-fexp(is))*spec(is)%tz*(zi*wd*(1.0-bd) - 2.0*vp))*ainv
a    = 1.0 + bd + fexp(is)*spec(is)%tz*(-zi*wd*(1.0+bd) - 2.0*vp)
b    = 1.0 - bd + fexp(is)*spec(is)%tz*(-zi*wd*(1.0-bd) + 2.0*vp)
```
Component form (exact, incl. GCC's known-zero lowering; `zi*wd = (0, wd)`):
```
bd = bkdiff(s); fe = fexp(s); tz = tz(s)
opb = 1.0 + bd ; omb = 1.0 - bd
c1  = (1.0 - fe)*tz ; c2 = fe*tz
v2  = 2.0*vp
wp  = wd*opb ; wm = wd*omb
D    = complex(opb + c1*v2,  c1*wp)          # (1+bd) + c1*(2vp, wd(1+bd))
ainv = cdiv(ONE, D)
N    = complex(omb - c1*v2,  c1*wm)          # (1-bd) + c1*(-(2vp), wd(1-bd))  [x + (-y) == x - y]
r    = cmul(N, ainv)
a    = complex(opb - c2*v2, -(c2*wp))
b    = complex(omb + c2*v2, -(c2*wm))
```
Then, **only if `grid_has_trapped_particles()`** (`:1899`), in this order for the same element:
```
if forbid[l,ig] or forbid[l,ig+1]:          r = 0 ; ainv = 0           # :1901-1904 (a, b NOT zeroed)
if σ == 1 and is_lower_bounce_point(ig+1, l): ainv = complex(1.0 + ainv.re, ainv.im)   # :1910-1914
if is_ttp[l,ig]:                                                        # :1917-1921
    ainv = cdiv(ONE, complex(1.0, c1*wdttp))     # 1/(1 + zi*(1-fe)*tz*wdttp)
    a    = complex(1.0, -(c2*wdttp))             # 1 - zi*fe*tz*wdttp
    r    = 0
```
Facts used later: at an interior lower bounce point `lb`, `forbid(lb-1)` holds, so
`ainv(lb-1,1) = (1,0)` and `r(lb-1,1) = 0` exactly; at an interior upper bounce point `ub`
(`forbid(ub+1)`), `r(ub,σ) = ainv(ub,σ) = 0` unless `ub` is a ttp point (then the ttp values);
in forbidden cells `r = ainv = 0`.

### 3.4 `init_homogeneous_g` — `g_h` (`dist_fn.fpp:1937-2125`)

Called from `init_invert_rhs` with `force_passing` absent (`:1928`). Per iglo, independent:
```
g_h[...] = 0                                                            # :1959 (whole array)
use_pass_homog_boundary = boundary_option_switch in {self_periodic, linked}   # :1965-1967
for each iglo:
    use_pass_homog = use_pass_homog_boundary or aky(y) == 0.0           # :1986 (force_passing=False)
    ilmin = 1 if use_pass_homog else ng2+1                              # :1988-1992
    if l < ilmin: continue                       # passing with zero-BC & ky≠0 → g_h ≡ 0   :1999
    is_wfb = (l == ng2+1); is_passing = (l <= ng2); is_trapped = (l > ng2+1)
    g2[:, :] = 0  (complex, (T,2))                                       # :2027
    # --- incoming boundary values ---
    if use_pass_homog and is_passing: g_h[-nt,1] = 1.0 ; g_h[nt,2] = 1.0          # :2033-2036
    if   is_wfb and mixed_wfb:   g_h[-nt,1] = wfb ; g_h[nt,2] = wfb    (complex(wfb,0))  # :2038-2041
    elif is_wfb and trapped_wfb: g_h[nt,2] = 1.0                                     # :2042-2045
    elif is_wfb and passing_wfb: g_h[-nt,1] = 1.0 ; g_h[nt,2] = 1.0                  # :2046-2050
    if is_trapped:
        for ig in [-nt, nt]: if not is_ttp[l,ig] and is_upper_bounce_point(ig,l): g2[ig,2] = 1.0   # :2054-2060
    # --- σ=2 sweep (right → left) ---
    for ig = nt-1 down to -nt:                                                       # :2071-2073
        g_h[ig,2] = g2[ig,2] - cmul(g_h[ig+1,2], r[ig,2])      # Fortran: -g_h(ig+1,2)*r(ig,2) + g2(ig,2)
    # --- lower-bounce matching ---
    if is_trapped:
        for ig in [-nt, nt-1]:                                                       # :2080-2096
            if is_ttp[l,ig]: continue                    # NB tests ig, not ig+1
            if is_lower_bounce_point(ig+1, l): g2[ig+1,1] = g_h[ig+1,2]
    if is_wfb and trapped_wfb: g_h[-nt,1] = g_h[-nt,2]                               # :2102-2105
    # --- σ=1 sweep (left → right) ---
    for ig in [-nt, nt-1]:                                                           # :2108-2114
        if is_ttp[l,ig]: continue                        # g_h[ig+1,1] keeps its value (0)
        g_h[ig+1,1] = g2[ig+1,1] - cmul(g_h[ig,1], r[ig,1])
# Communication (linked only, after the loop)                                         # :2120-2123
if g_adj allocated and not no_connections:
    fill(links_h, g_h, g_adj) ; if mixed_wfb: fill(wfb_h, g_h, g_adj)
```
Notes:
* `g2` entries are complex; the σ=1 source `g2[lb,1] = g_h[lb,2]` is complex.
* For trapped `l` the σ=2 sweep yields `g_h[ub,2] = 1` (real part exactly 1) at each non-ttp
  upper bounce point (`r[ub,2] = 0`) and `0` in forbidden regions; the σ=1 sweep starts each well from
  `g_h[lb,1] = g_h[lb,2]` (because `r[lb-1,1] = 0`).
* The upper-bounce `g2` assignment is skipped at `ig = ntgrid` *in effect*: the σ=2 loop never
  uses `g2[nt,2]` (pitfall §6.9).
* At a ttp point `g_h = 0` for both σ (all contributions vanish exactly).
* `g_h` is what the JAX side stores; `g_adj`'s h-slots are just `g_h[nt,1]` / `g_h[-nt,2]` of
  other cells (§3.11) and need not be stored.

### 3.5 Linked-boundary topology (init, `boundary_option_linked` only)

All of this is static; export `itleft, itright, l_links, r_links, no_connections` and build the
index tables in Python at load time. Formulas for reference/verification:

**jshift0** (`:2401-2436`):
```
theta_extent = theta(ntgrid) - theta(-ntgrid)
if naky > 1 and ntheta0 > 1:                    spacing = theta0(2,2) - theta0(1,2)
elif naky == 1 and ntheta0 > 1 and aky(1) != 0: spacing = theta0(2,1) - theta0(1,1)
else:                                           spacing = theta_extent
jshift0 = nint(theta_extent/spacing)            # abort if |jshift0*spacing - extent| > 100*eps
```
In box mode `theta0(it,ik) = akx(it)/(aky(ik)*shat)` and `dkx = dky*2π(2nperiod-1)|shat|/jtwist`
(`kt_grids_box.f90:371-396`), so `jshift0 = jtwist*sign(shat)` — **negative for negative
shear**, which mirrors the link direction.

**itleft/itright** (`:2441-2513`), for every `(it, ik)` (1-based):
```
it0 = it - ntheta0 - 1 if it > (ntheta0+1)/2 else it - 1        # signed kx index
if ik == 1:
    if aky(1) != 0 and naky == 1: itl = it0 + jshift0 ; itr = it0 - jshift0
    else:                         itl = it0 ; itr = it0         # ky=0: self-connection
else:  itl = it0 + (ik-1)*jshift0 ; itr = it0 - (ik-1)*jshift0
remap(v) = v+1                 if 0 <= v < (ntheta0+1)/2
         = v+ntheta0+1         if (ntheta0+1)/2 < v+ntheta0+1 <= ntheta0
         = -1                  otherwise
itleft(it,ik) = remap(itl) ; itright(it,ik) = remap(itr)
```
So for `ik = 1` (ky = 0 in box) `itleft = itright = it` (self-periodic link).

**l_links / r_links** (`:2608-2640`): number of hops following `itleft` (resp. `itright`) until
the next index is `< 0` or equal to the current one (self-loop) — ky=0 cells have 0 links.
Abort if > 5000. `n_links(1) = 0 if l_links==0 else 2*l_links-1`; `n_links(2)` likewise with
`r_links` (`:2659-2671`). `n_links_max = maxval(n_links)`, `+3` if `mixed_wfb` and > 0
(`:3544-3547`) (= `2*Nmax` where `Nmax` is the longest chain).

**connections / neighbor** (`:2518-2565`), per `(it, ik, il)` (independent of ie, is):
```
excluded(l) = il_is_trapped(l) or (il_is_wfb(l) and trapped_wfb)
iglo_left  exists ⇔ not excluded(l) and itleft(it,ik)  >= 1
iglo_right exists ⇔ not excluded(l) and itright(it,ik) >= 1
neighbor = iglo_left exists or iglo_right exists
```
(For ky=0 cells `neighbor = .true.` because of the self-link.)

**save_h** (`:2570-2604`):
```
save_h[1] = (iglo_left exists)  and aky(y) != 0.0
save_h[2] = (iglo_right exists) and aky(y) != 0.0
if il_is_wfb(l) and neighbor and not passing_wfb and aky(y) != 0.0: save_h[1] = save_h[2] = True
```

**Chains (supercells).** For each `ik`, follow `itright` from every head `it` with
`l_links(it,ik) == 0`: `chain = [it, itright(it), itright(itright(it)), …]` until the next value
is `< 1` or equals the current one. Length `N = r_links(head)+1 = n_cells(it,ik)`; the cell at
0-based chain position `k` has `l_links = k`, `r_links = N-1-k`. Every `(it,ik)` belongs to
exactly one chain; ky=0 chains have `N=1`. Equivalent to `supercell_labels`
(`:2355-2397`, labels via `get_leftmost_it`). `get_leftmost_it(it,ik)` (`:3889-3903`) =
`chain[0]`, `get_rightmost_it` (`:3906-3920`) = `chain[N-1]`.

**Communication patterns** (the JAX side replaces them by gathers along the chain; all are
copies of single boundary values; `to`/`from` semantics: `fill(f, src, dst)` does
`dst(to_idx) = src(from_idx)`, `externals/utils/redistribute.f90:7685-7739`):

* `links_p` (`:2759-3022`): passing (and `passing_wfb`) only, `ncell > 1`. For a sender `S`
  and each `j = 1..r_links(S)`: `g_adj(j, 1, R_j) = gnew(ntgrid, 1, S)` where `R_j` is `j`
  hops to the right of `S`. For `j = 1..l_links(S)`: `g_adj(j, 2, L_j) = gnew(-ntgrid, 2, S)`,
  `L_j` = `j` hops to the left. So at receiver `R` (chain position `k`):
  `g_adj(j,1,R) = p1[k-j]` (j=1..k), `g_adj(j,2,R) = p2[k+j]` (j=1..N-1-k), with
  `p1[m] = gnew(nt,1,cell m)`, `p2[m] = gnew(-nt,2,cell m)`.
* `links_h` (same routine, "special restriction"): `g_adj(2*l_links(S)+j, 1, R_j) =
  g_h(ntgrid,1,S)` **only if `l_links(S) > 0`**; `g_adj(2*r_links(S)+j, 2, L_j) =
  g_h(-ntgrid,2,S)` only if `r_links(S) > 0`. At receiver position `k`:
  `g_adj(2k-j, 1, R) = h1[k-j]` for `j = 1..k-1`, `g_adj(2(N-1-k)-j, 2, R) = h2[k+j]` for
  `j = 1..N-2-k`, with `h1[m] = g_h(nt,1,cell m)`, `h2[m] = g_h(-nt,2,cell m)`. (The leftmost
  cell's `h1` / rightmost cell's `h2` are never needed.)
* `wfb_p`, `wfb_h` (`:3028-3224`, only `mixed_wfb`, `il == ng2+1`, `ncell > 1`): every cell
  `S` of the chain sends to **every** cell (itself included):
  `g_adj(r_links(S)+1, 1, ·) = gnew(nt,1,S)`, `g_adj(2N - r_links(S), 1, ·) = g_h(nt,1,S)`,
  `g_adj(l_links(S)+1, 2, ·) = gnew(-nt,2,S)`, `g_adj(2N - l_links(S), 2, ·) = g_h(-nt,2,S)`.
  In chain positions (S at position m): σ=1 p-slot `N-m` holds `p1[m]`, h-slot `N+1+m` holds
  `h1[m]`; σ=2 p-slot `m+1` holds `p2[m]`, h-slot `2N-m` holds `h2[m]`.
* `incoming_links` (`:3229-3508`, only if `start_from_previous_solution`): sends
  `gnew(-nt,1,S)` → `g_work(nt,1,left neighbour)` and `gnew(nt,2,S)` → `g_work(-nt,2,right
  neighbour)`; for `mixed_wfb` additionally wraps rightmost→leftmost and leftmost→rightmost.
* `no_connections` = no receiver anywhere (`max_allreduce(maxval(nn_to)) == 0`,
  `:2870-2876`); then `g_adj` is not allocated and `apply_linked_boundary_conditions` returns
  immediately (`:5948`).

### 3.6 `invert_rhs` (`dist_fn.fpp:6152-6194`)

```
sourcefac = ...                                    # :6166-6176 (02 §3.2)
g_work = 0                                         # :6179
g_work = source(phi, apar, bpar, phinew, aparnew, bparnew, istep, sourcefac)   # :6181 (02)
invert_rhs_1(phinew, bparnew, source=g_work)       # :6184  (§3.7)
if boundary_option_switch == linked:
    apply_linked_boundary_conditions(gnew, g_h, g_adj, phinew, bparnew)       # :6188-6191 (§3.11)
```

### 3.7 `invert_rhs_1` (`dist_fn.fpp:5440-5707`)

Per iglo (loop `:5466-5703`), completely independent across iglo. `nt = ntgrid`.

```
# ---- masks / early exits ----
if kwork_filter[y,x]: continue               # gnew unchanged (pass-through)        :5475
if ieqzip[y,x]:       continue               # gnew unchanged                       :5477
if is_passing_hybrid_electron(s,y,l):        # :5487-5490
    gnew[:, :] = set_h_zero(...)  (§3.9) ; continue

# ---- flags ----                                                                    :5503-5558
linked = (boundary_option_switch == linked)
use_pass_homog = (boundary_option_switch in {self_periodic, linked}) or aky(y) == 0.0
speriod_flag   = (aky(y) == 0.0)             # only defined/used when linked
is_wfb = (l == ng2+1) ; is_passing = (l <= ng2) ; is_trapped = (l > ng2+1)
is_passing_like = is_passing or (is_wfb and passing_wfb)
if linked:
    is_self_periodic = (speriod_flag and is_passing)
                    or (speriod_flag and is_wfb and not trapped_wfb)
                    or (is_wfb and not neighbor[x,y,l] and mixed_wfb)
else:
    is_self_periodic = (use_pass_homog and is_passing_like) or (is_wfb and mixed_wfb)

# ---- initial values ----                                                           :5563-5575
gnew[:, :] = g[:, :] if start_from_previous_solution else 0

# ---- incoming boundary values ----                                                 :5580-5583
if is_passing_like and not is_self_periodic:
    apply_parallel_boundary_conditions(..., is_left = (l_links[y,x] == 0),
                                            is_right = (r_links[y,x] == 0))      # §3.8
    # (non-linked: l_links = r_links = 0 everywhere → both ends)

# ---- σ=2 sweep, right → left ----                                                  :5596-5598
for ig = nt-1 down to -nt:
    gnew[ig,2] = cmul(ainv[ig,2], src[ig,2]) - cmul(gnew[ig+1,2], r[ig,2])
    # Fortran: -gnew(ig+1,2)*r(ig,2) + ainv(ig,2)*source(ig,2)

# ---- trapped: lower-bounce source (modifies the source in place) ----              :5604-5622
if is_trapped:
    for ig in [-nt, nt-1]:
        if is_ttp[l,ig]: continue
        if is_lower_bounce_point(ig+1, l): src[ig,1] = gnew[ig+1,2]

# ---- trapped wfb: bounce at the left end ----                                      :5628-5631
if is_wfb and trapped_wfb: gnew[-nt,1] = gnew[-nt,2]

# ---- σ=1 sweep, left → right ----                                                  :5634-5645
for ig in [-nt, nt-1]:
    if not is_ttp[l,ig]:
        t = cmul(ainv[ig,1], src[ig,1])                         # :5637
        gnew[ig+1,1] = t - cmul(gnew[ig,1], r[ig,1])           # :5638 (two statements, same value)
    else:
        if not forbid[l,ig]: gnew[ig,1] = gnew[ig,2]           # :5643 (always true for ttp)
        # NB gnew[ig+1,1] is NOT written: it keeps its initial value (0, or g if sfps)

# ---- periodicity ----                                                              :5651
if is_self_periodic: self_periodic(...)                                           # §3.10

# ---- trapped: add homogeneous solution well by well ----                           :5658-5671
if is_trapped:
    beta = 0
    for ig = nt-1 down to -nt:
        if forbid[l,ig] or is_ttp[l,ig]: beta = 0 ; continue
        elif forbid[l,ig+1]:
            beta = cdiv(gnew[ig,1] - gnew[ig,2], complex(1.0 - g_h[ig,1].re, -g_h[ig,1].im))
        gnew[ig,1] = gnew[ig,1] + cmul(beta, g_h[ig,1])
        gnew[ig,2] = gnew[ig,2] + cmul(beta, g_h[ig,2])
    # ig = nt is never corrected

# ---- trapped wfb: bounce at the right end ----                                     :5679-5688
if is_wfb and trapped_wfb:
    beta = cdiv(gnew[nt,1] - gnew[nt,2], 1.0 - g_h[nt,1])
    gnew[:, :] = gnew[:, :] + cmul(beta, g_h[:, :])            # all ig, both σ

# ---- optional zeroing ----                                                         :5696-5701
if zero_forbid: gnew[where forbid[l,:], :] = 0

# after the iglo loop:                                                               :5706
if def_parity: enforce_parity(parity_redist)                   # GATED (§3.12)
```

Behaviour summary by class (default `nonad_zero=.true.`, `wfbbc_option='mixed'`):

| class | linked, ky≠0, in a chain (N≥2) | linked, ky≠0, isolated (N=1) | linked, ky=0 | zero BC (ballooning), ky≠0 | zero BC, ky=0 | self-periodic |
|---|---|---|---|---|---|---|
| passing | BC only at chain ends (`l_links==0` left, `r_links==0` right); then linked fix (§3.11.2) | BC both ends | self_periodic | BC both ends | self_periodic | self_periodic |
| mixed wfb | **no BC**, cyclic linked fix over chain (§3.11.3) | self_periodic (with `dadj`) | self_periodic (with `dadj`) | self_periodic (with `dadj`) | self_periodic (with `dadj`) | self_periodic (with `dadj`) |
| passing wfb | as passing | BC both ends | self_periodic (`dadj=0`) | BC both ends | self_periodic | self_periodic |
| trapped wfb | trapped-wfb bounce at ±nt | same | same | same | same | same |
| trapped | bounce matching (β per well) | same | same | same | same | same |

### 3.8 `apply_parallel_boundary_conditions` (`dist_fn.fpp:5710-5849`, pure)

`adjust(ig) := get_adjust(ig, it, ik, ie, is, iglo, fphi, fbpar, phinew, bparnew)`
(`dist_fn_arrays.fpp:178-188`), verbatim
`2 * vperp2(ig,iglo) * aj1(ig,iglo) * bpar(ig,it,ik) * facbpar + spec(is)%zt * phi(ig,it,ik) * aj0(ig,iglo) * facphi * nonmaxw_corr(ie,is)`:
```
A = rmul(rmul(bparnew[y,x,ig], (2.0*vperp2[ig])*aj1[ig]), fbpar)          # (((2*vperp2)*aj1)*bpar)*fbpar
B = rmul(rmul(rmul(rmul(phinew[y,x,ig], zt(s)), aj0[ig]), fphi), nonmaxw_corr(e,s))   # (((zt*phi)*aj0)*fphi)*nonmaxw
adjust(ig) = A + B
```
(`rmul(z, c)` multiplies both components by `c`; the real prefactor `(2.0*vperp2)*aj1` is formed
first, then applied to `bpar`.)

```
if not (is_left or is_right): return                                       # :5722
if boundary_off_grid:            # GATED                                    :5723-5786
    if is_left:
        corr = (3*adjust(-nt) - adjust(-nt+1)) if nonad_zero else 0
        gnew[-nt,1] = cdiv(cmul(ainv[-nt,1], src[-nt,1]) - corr, complex(3.0 + r[-nt,1].re, r[-nt,1].im))
    if is_right:
        corr = (3*adjust(nt) - adjust(nt-1)) if nonad_zero else 0
        gnew[nt,2]  = cdiv(cmul(ainv[nt-1,2], src[nt-1,2]) - corr, complex(3.0 + r[nt-1,2].re, r[nt-1,2].im))
elif exponential_boundary:       # GATED                                    :5787-5832
    f = exponential_boundary_factor
    if is_left:
        corr = ((1 + f*delthet(-nt))*adjust(-nt) - adjust(-nt+1)) if nonad_zero else 0
        gnew[-nt,1] = cdiv(cmul(ainv[-nt,1], src[-nt,1]) - corr,
                           complex((r[-nt,1].re + 1) + f*delthet(-nt), r[-nt,1].im))
    if is_right:
        corr = ((1 + f*delthet(nt-1))*adjust(nt) - adjust(nt-1)) if nonad_zero else 0
        gnew[nt,2]  = cdiv(cmul(ainv[nt-1,2], src[nt-1,2]) - corr,
                           complex((r[nt-1,2].re + 1) + f*delthet(nt-1), r[nt-1,2].im))
elif nonad_zero:                 # DEFAULT: h = 0 for incoming particles      :5833-5844
    if is_left:  gnew[-nt,1] = -adjust(-nt)
    if is_right: gnew[nt,2]  = -adjust(nt)
else:                            # nonad_zero = F: leave gnew (0, or g if sfps) :5845-5847
    pass
```

### 3.9 `set_h_zero` (`dist_fn_arrays.fpp:246-294`) — passing hybrid electrons

```
facphi = -fphi ; facbpar = -fbpar
phi_factor  = (zt(s)*nonmaxw_corr(e,s))*facphi  if facphi  != 0 else 0.0
bpar_factor = 2.0*facbpar                        if facbpar != 0 else 0.0
for ig in [-nt, nt]:
    adj = rmul(bparnew[y,x,ig], (bpar_factor*vperp2[ig])*aj1[ig]) + rmul(rmul(phinew[y,x,ig], phi_factor), aj0[ig])
    gnew[ig,1] = adj ; gnew[ig,2] = adj
```
(Fortran: `bpar_factor*vperp2*aj1*bpar + phi_factor*phi*aj0`, left-to-right.)

### 3.10 `self_periodic` (`dist_fn.fpp:5864-5906`)

```
if is_wfb and nonad_zero and not passing_wfb:                                     # :5879
    adjl = rmul(rmul(bparnew[y,x,-nt], (2.0*vperp2[-nt])*aj1[-nt]), fbpar)
         + rmul(rmul(rmul(rmul(phinew[y,x,-nt], zt(s)), aj0[-nt]), nonmaxw_corr(e,s)), fphi)
    adjr = (same at +nt)
    dadj = adjl - adjr
    # NB order: ... *nonmaxw_corr*fphi  (get_adjust uses ... *facphi*nonmaxw_corr)
else:
    dadj = 0
if g_h[nt,1] != (1,0):            # exact complex /= 1. test                     :5893-5896
    beta = cdiv((gnew[nt,1] - gnew[-nt,1]) - dadj, 1.0 - g_h[nt,1])
    gnew[:,1] = gnew[:,1] + cmul(beta, g_h[:,1])
if g_h[-nt,2] != (1,0):                                                            :5898-5901
    beta = cdiv((gnew[-nt,2] - gnew[nt,2]) + dadj, 1.0 - g_h[-nt,2])
    gnew[:,2] = gnew[:,2] + cmul(beta, g_h[:,2])
```
In practice `dadj ≠ 0` only for the mixed wfb (`is_self_periodic` excludes trapped wfb; passing
wfb has `dadj = 0`).

### 3.11 `apply_linked_boundary_conditions` (`dist_fn.fpp:5919-6149`)

#### 3.11.1 Fortran control flow

```
if no_connections: return                                                         :5948
if start_from_previous_solution: (GATED, §3.12.1) p-values are differences        :5951-5982
else: fill(links_p, gnew, g_adj); if mixed_wfb: fill(wfb_p, gnew, g_adj)           :5984-5985
for each local iglo:                                                              :5996
    if not any(save_h[:,iglo]): continue
    if kwork_filter[y,x]: continue
    ncell = l_links + r_links + 1 ; if ncell == 1: continue
    if is_passing_hybrid_electron(s,y,l): continue
    if il_is_wfb(l) and mixed_wfb:  mixed-wfb cyclic algorithm (§3.11.3)
    else:                            passing open-chain algorithm (§3.11.2)
```
All `b0` values are computed from **pre-correction** `gnew` boundary values (the fills happen
before the loop); `ieqzip` is **not** checked here (§6.7).

#### 3.11.2 Passing (and passing wfb): open chain, zero incoming at the chain ends

For cell at chain position `k` (`l_links = k`), σ=1, executed iff `save_h[1] and k > 0`
(`:6121-6131`); Fortran:
```fortran
n = n_links(1,it,ik)                          ! = 2k-1
b0 = g_adj(1,1,iglo)
fac = 1.0
do i = 2, l_links(it,ik)
   fac = fac * g_adj(n+2-i,1,iglo)
   b0  = b0 + g_adj(i,1,iglo) * fac
end do
g_in(:,1,iglo) = g_in(:,1,iglo) + b0*g_h(:,1,iglo)
```
In chain terms (slot `i` ↔ `p1[k-i]`, slot `n+2-i = 2k+1-i` ↔ `h1[k-i+1]`):
```
b0  = p1[k-1]
fac = ONE
for i in 2..k:
    fac = cmul(fac, h1[k-i+1])
    b0  = b0 + cmul(p1[k-i], fac)
gnew[:,1] += cmul(b0, g_h[:,1])        # all ig
```
i.e. `b0 = p1[k-1] + p1[k-2]·h1[k-1] + p1[k-3]·(h1[k-1]h1[k-2]) + …`, accumulated **nearest
neighbour first, products built outward from the cell**. This is *not* the Horner/forward
recurrence `a_k = p_{k-1} + h_{k-1} a_{k-1}`; the two differ in rounding, so every cell's sum
must be formed independently (O(N²) per chain).

σ=2, iff `save_h[2] and R = N-1-k > 0` (`:6133-6143`), with `p2`, `h2`:
```
b0  = p2[k+1]
fac = ONE
for i in 2..R:
    fac = cmul(fac, h2[k+i-1])
    b0  = b0 + cmul(p2[k+i], fac)
gnew[:,2] += cmul(b0, g_h[:,2])
```
`fac = ONE` then `cmul(ONE, h)` is exactly `h` (also with FMA), but keep the multiply to stay
literal.

#### 3.11.3 Mixed wfb: cyclic chain (self-periodic over the supercell, on h)

For every mixed-wfb cell of a chain with `N ≥ 2` and `aky ≠ 0` (`save_h[1] = save_h[2] = True`
there), at position `k`. Fortran σ=1 (`:6049-6094`):
```fortran
facd = 1./(1. - product(g_adj(ncell+1:2*ncell,1,iglo)))
b0 = 0. ; fac = 1
do i = ncell-1, 1, -1
   fac = fac * g_adj(ncell+i+1,1,iglo)
   b0  = b0 + fac * g_adj(ncell+1-i,1,iglo)
end do
b0 = (b0 + g_adj(1,1,iglo)-dadj)*facd
do i = 1, l_links(it,ik)
   b0 = b0 * g_adj(ncell+i,1,iglo) + g_adj(ncell+1-i,1,iglo)
end do
g_in(:,1,iglo) = g_in(:,1,iglo) + b0*g_h(:,1,iglo)
```
In chain terms (σ=1 p-slot `N-m` ↔ `p1[m]`, h-slot `N+1+m` ↔ `h1[m]`):
```
P = ONE
for m in 0..N-1: P = cmul(P, h1[m])                   # PRODUCT intrinsic: left-to-right over slots N+1..2N
facd = cdiv(ONE, complex(1.0 - P.re, -P.im))
b0 = 0 ; fac = ONE
for i = N-1 down to 1:
    fac = cmul(fac, h1[i])
    b0  = b0 + cmul(fac, p1[i-1])
A0 = cmul((b0 + p1[N-1]) - dadj, facd)                # amplitude entering the leftmost cell
A  = A0
for m in 0..k-1: A = cmul(A, h1[m]) + p1[m]           # cascade to the right
gnew[:,1] += cmul(A, g_h[:,1])
```
σ=2 (`:6098-6114`) (p-slot `m+1` ↔ `p2[m]`, h-slot `2N-m` ↔ `h2[m]`):
```
P = ONE
for m = N-1 down to 0: P = cmul(P, h2[m])             # slots N+1..2N are cells N-1, …, 0
facd = cdiv(ONE, complex(1.0 - P.re, -P.im))
b0 = 0 ; fac = ONE
for i = N-1 down to 1:
    fac = cmul(fac, h2[N-1-i])
    b0  = b0 + cmul(fac, p2[N-i])
B0 = cmul((b0 + p2[0]) + dadj, facd)                  # amplitude entering the rightmost cell
B  = B0
for i in 1..(N-1-k): B = cmul(B, h2[N-i]) + p2[N-i]   # cascade to the left
gnew[:,2] += cmul(B, g_h[:,2])
```
Each cell of the chain recomputes `P`, `facd`, `b0`, `A0` from identical inputs with identical
operations, so the results are bitwise identical across the chain; the cascade *is* a forward
recurrence whose prefixes are shared → compute it once per chain and gather (`A` after `k`
steps for position `k`; `B` after `N-1-k` steps).

#### 3.11.4 `dadj` for the mixed wfb in a chain (`:6013-6047`)

Only `if nonad_zero`, else `dadj = 0`. With `itl = chain[0]`, `itr = chain[N-1]`, `l = ng2+1`:
```
vperp2left  = (bmag(-nt)*al(l))*energy(e,s)                  # :6019
vperp2right = (bmag(nt)*al(l))*energy(e,s)                   # :6020
argl = (bess_fac(s)*smz(s))*sqrt(((energy(e,s)*al(l))/bmag(-nt))*kperp2(-nt,itl,ik))   # :6022
argr = (bess_fac(s)*smz(s))*sqrt(((energy(e,s)*al(l))/bmag(nt))*kperp2(nt,itr,ik))     # :6023
aj0l = j0(argl) ; aj0r = j0(argr) ; aj1l = j1(argl) ; aj1r = j1(argr)                   # libm
adjl = rmul(rmul(bparnew[ik,itl,-nt], (2.0*vperp2left)*aj1l), fbpar)
     + rmul(rmul(rmul(rmul(phinew[ik,itl,-nt], zt(s)), aj0l), nonmaxw_corr(e,s)), fphi)
adjr = (same with right quantities at +nt, itr)
dadj = adjl - adjr
```
`argl` is associated differently from `init_bessel` (`arg = bess_fac*smz*sqrt(energy*al*kperp2/bmag)`,
`:1784-1786`), so `aj0l` need not equal `aj0(-nt, iglo_left)` bitwise. **EXPORT** new arrays
`wfb_aj0l, wfb_aj1l` (evaluated with `kperp2(-nt, it, ik)`) and `wfb_aj0r, wfb_aj1r` (with
`kperp2(nt, it, ik)`) of shape NumPy `(S,E,Y,X)` for every `it` (JAX gathers at `itl`/`itr`),
plus `vperp2left/right` `(S,E)` (or recompute: two exact products).

### 3.12 GATED paths (documented for completeness)

#### 3.12.1 `start_from_previous_solution` (`:5563-5572`, `:5951-5982`)
`gnew = g` initially (so incoming boundaries of internal cells and ttp-adjacent points start
non-zero). In the linked fix the p-values become `g_work(nt,1,S) = gnew(nt,1,S) -
gnew(-nt,1,right(S))` and `g_work(-nt,2,S) = gnew(-nt,2,S) - gnew(nt,2,left(S))` (mixed wfb
wraps rightmost↔leftmost), via `fill(incoming_links, gnew, g_work)` into a zeroed `g_work`
then the difference (`:5968-5971`), then `links_p`/`wfb_p` fills from `g_work`.

#### 3.12.2 `enforce_parity` (`:4796-4853`)
`mult = +1 if even else -1`. Ballooning (not linked): for `j = 1..nt`:
`gnew(-j,1) = mult*gnew(j,2)`, `gnew(j,1) = mult*gnew(-j,2)` (`ig = 0` untouched). Linked:
`scatter(parity_redist, gnew, gnew)` performs `gnew(-ig,1,iglo) = gnew(ig,2,conn(iglo))` for all
`ig` (including 0), where `conn` (`get_parity_conn`, `:4742-4794`) is the cell in the same chain
with `r_links(conn) = l_links(iglo)`; then `gnew(:,1,:) *= mult` if odd. Called at the end of
`invert_rhs_1` and again in `timeadv` (`:4091`).

#### 3.12.3 `boundary_off_grid`, `exponential_boundary` — formulas in §3.8.

----------------------------------------------------------------------------------------------

## 4. Configuration options

| option | namelist | default | effect here | decision |
|---|---|---|---|---|
| `boundary_option` | `dist_fn_knobs` | `'default'` → `linked` if `grid_option='box'`, else `zero`; `periodic` if `|shat| <= 1e-5` (`dist_fn.fpp:1295-1310`, `theta_grid.f90:38,674-681`) | selects zero / self-periodic / linked logic (§3.7, §3.4) | SUPPORT `default, zero, unconnected, self-periodic, periodic, kperiod=1, linked`; `alternate-zero` behaves exactly like `zero` (only `case default` branches) → SUPPORT as alias. `linked` with `grid_option /= 'box'` → GATE (GS2 itself warns "probably an error", `:592-599`; untested `ik=1, aky≠0` self-link path) |
| `nonad_zero` | `dist_fn_knobs` | `.true.` | incoming BC `h=0` (`gnew=-adjust`) and `dadj`; `.false.` leaves `gnew=0` and `dadj=0` | SUPPORT both |
| `wfbbc_option` | `le_grids_knobs` | `'default'` (= `mixed`) | `passing_wfb/trapped_wfb/mixed_wfb` (`le_grids.f90:413-427`); `dist_fn_knobs:wfbbc_option` must stay `'unset'` (abort, `:1281-1284`) | SUPPORT `default/mixed/passing/trapped` |
| `wfb` | `dist_fn_knobs` | `1.0` | `g_h` incoming value for mixed wfb (`:2040-2041`) | SUPPORT |
| `zero_forbid` | `dist_fn_knobs` | `.false.` | zero `gnew` in forbidden region (`:5696-5701`) | SUPPORT |
| `boundary_off_grid` | `dist_fn_knobs` | `.false.` | experimental half-cell BC (§3.8) | GATE (experimental; formula given) |
| `exponential_boundary` (+`exponential_boundary_factor`=1.0) | `dist_fn_knobs` | `.false.` | experimental decay BC (§3.8) | GATE (experimental; formula given) |
| `start_from_previous_solution` | `dist_fn_knobs` | `.false.` | `gnew=g` initial; needs `incoming_links` | GATE (documented experimental; extra communication pattern) |
| `def_parity`, `even` | `dist_fn_knobs` | `.false.`, `.true.` | parity projection (§3.12.2) | GATE (never tested nonlinearly per `check_dist_fn` `:684-690`; needs parity map) |
| `bakdif` | `dist_fn_species_knobs` | `0.0` | `bd` in §3.3 | SUPPORT |
| `fexpr` | `dist_fn_species_knobs` | `0.4` | `fe` in §3.3 | SUPPORT |
| `mult_imp` | `dist_fn_knobs` | `.false.` (forced `.false.` if nonlinear) | allows per-species `bkdiff` | SUPPORT (export effective `bkdiff`) |
| `driftknob`, `tpdriftknob`, `mach`, `vparknob` | `dist_fn_knobs` | 1, sentinel→`driftknob`, 0, 1 | enter `wdrift`/`wdriftttp`/`vpa` upstream (02) | SUPPORT (via 02 recipes) |
| `bouncefuzz`, `trapped_particles`, `ngauss`, `npassing`, `radau_gauss_grid` | `le_grids_knobs` | `10*eps`, `.true.`, 5, −1, `.true.` | determine `ng2`, `nlambda`, flags (exported) | SUPPORT |
| `fphi`, `fbpar`, `beta` | `knobs` / `parameters` | 1.0, **−1.0**, — | `adjust`, `dadj`, `set_h_zero` | SUPPORT (export effective values) |
| `eqzip_option` | `knobs` | `'none'` | `ieqzip` skips in `invert_rhs_1` but **not** in the linked fix | GATE (secondary/tertiary studies only; odd interaction §6.7) |
| species `type='hybrid_electron'` | `species_parameters_*` | — | `set_h_zero` for passing ky≠0 (§3.9); skipped in linked fix | SUPPORT |
| `grid_option`, `jtwist`, `nperiod`, `ntheta0`, `naky`, `shat` | `kt_grids_*`, `theta_grid_*` | | topology (exported `itleft/itright`), `kwork_filter(1,1)` in box | SUPPORT box/range/single |
| `delt` and dt changes (`reinit_knobs`) | `knobs` | | rebuild `vpar, wdrift, wdriftttp, a, b, r, ainv, g_h` (§3.3-3.4) | SUPPORT |
| `wstar_units` | `knobs` | `.false.` | enters `tunits/wunits` upstream (`tunits=0` for ky=0 → garbage, `gs2_time.f90:108-120`) | SUPPORT (exported units) |
| `g_exb` (flow shear) | `dist_fn_knobs` | 0 | no effect on invert coefficients (remap only) | SUPPORT |
| `field_option` | `fields_knobs` | `'default'` | invert_rhs identical for implicit/local/gf_local | SUPPORT |
| LOWFLOW (compile-time) | — | off | complex `wd` + `wstar_neo` terms (`:1873-1876`) | GATE |
| OpenMP (compile-time) | — | off | none (iglo loops independent) | SUPPORT (no effect) |

----------------------------------------------------------------------------------------------

## 5. GPU mapping notes

### 5.1 Layout

Keep the canonical batch `B = (S,E,L,Y,X)` flattened (C order) and put the theta axis **first**
for the sweeps: `ainv, r, src, gnew, g_h` as `(T, 2, B)` (or `(2, T, B)`), so `lax.scan` iterates
over the leading axis and every step is one vectorised complex op over `B` elements. The masks
depend on `(l, ig)` only (`forbid`, `is_ttp`, bounce flags) or `(y, x)` (`kwork_filter`,
`ieqzip`, `l_links==0`) and broadcast to `B`.

### 5.2 Independent vs sequential work

* `init_invert_rhs` (§3.3), BC values (§3.8), `set_h_zero`, `self_periodic`, trapped-wfb end
  correction, `zero_forbid`: fully element-wise → `vmap`/broadcast.
* σ=2 sweep: `lax.scan(reverse=True)` over `ig = nt-1 … -nt`, carry `c = gnew[ig+1,2]` (init
  = BC value or initial value at `nt`), output `gnew[ig,2] = cmul(ainv,src) - cmul(c,r)`.
  `2*ntgrid` sequential steps. Do **not** use an associative/parallel-prefix scan (different
  rounding).
* Lower-bounce source fix: element-wise `src1 = where(trapped & lbsrc[l,ig], gnew2[ig+1], src1)`
  with static `lbsrc[l,ig] = not is_ttp[l,ig] and is_lower_bounce_point(ig+1,l)`,
  `ig ∈ [-nt, nt-1]`.
* σ=1 sweep with the ttp branch: `lax.scan` over `ig = -nt … nt-1`, carry `c` = current stored
  value of `gnew[ig,1]`; per step:
  ```
  if not ttp[l,ig]: out[ig] = c ;                         c_next = cmul(ainv1[ig],src1[ig]) - cmul(c, r1[ig])
  else:             out[ig] = gnew2[ig] if not forbid[l,ig] else c ;  c_next = init1[ig+1]
  ```
  and `out[nt] = c_final`. `init1` is the initial `gnew[:,1]` (0, or `g` if sfps; the BC and the
  trapped-wfb copy only affect `ig=-nt`, which is the initial carry). Implement the branch with
  `jnp.where` on both outputs.
* Trapped β (§3.7, only `l > ng2+1`): replace the reverse loop by a **static gather**:
  precompute, by simulating the loop on the masks only,
  ```
  cur = -1
  for ig = nt-1 down to -nt:
      if forbid[l,ig] or is_ttp[l,ig]: cur = -1 ; upd[l,ig] = False
      else:
          if forbid[l,ig+1]: cur = ig
          top[l,ig] = cur ; upd[l,ig] = True
  upd[l,nt] = False
  ```
  Then `beta_at[ig] = cdiv(g1[ig]-g2[ig], 1-gh1[ig])` evaluated only where
  `ig` is a top (pre-update values — each point is updated exactly once, by its own well, and a
  top is updated only after its β is formed, so the parallel version is bitwise identical), and
  `gnew[ig,σ] += cmul(beta_at[top[ig]] or 0, g_h[ig,σ])` where `upd`. Apply the update with
  `beta = 0` where Fortran does (`upd` true but `top=-1`) to reproduce signed zeros, or skip it
  (values identical up to the sign of zero).
* `init_homogeneous_g`: two scans like the gnew sweeps (σ=2 reverse with `g2` source; σ=1 forward
  where a ttp step emits `out[ig] = c` and sets `c_next = 0`). Only on dt change.

### 5.3 Linked chains: padding and batching

Build at load time (NumPy, from exported `itright`, `l_links`, `r_links`, `aky`):
* `chains`: list of `(ik, [it_0, …, it_{N-1}])` with `N ≥ 2` and `aky(ik) ≠ 0`.
* Pad to `Nmax = max N` (≈ `2*ntheta0/|jtwist|`-ish at low ky, 1 at the highest ky): integer
  table `C[c, pos]` = flat `(ik, it)` index or a dummy, `valid[c,pos]`, `pos_of[(ik,it)]`,
  `chain_of[(ik,it)]`, `N_of[c]`. Optionally group chains by length class (`N_class`/`M_class`)
  to avoid padding waste; the arithmetic per valid cell must not depend on padding.
* Gather endpoint arrays for all `(c, pos, l, e, s)`:
  `p1 = gnew[nt,1]`, `p2 = gnew[-nt,2]`, `h1 = g_h[nt,1]`, `h2 = g_h[-nt,2]` (h-tables only
  change on dt change).
* Passing (§3.11.2): carry `(fac, b0)` per `(c, k, l, e, s)`; `b0 = p1[k-1]`, `fac = ONE`; then a
  `fori_loop` over `i = 2..Nmax-1` doing
  `fac' = cmul(fac, h1[k-i+1]); b0' = b0 + cmul(p1[k-i], fac')`, applied only where `i <= k`
  (`jnp.where`, never multiply-by-zero masking). Same for σ=2 with `i <= N-1-k`. Cost
  `O(Nmax² · B/ntheta0)`, negligible versus the sweeps.
* Mixed wfb (§3.11.3): per chain (only `l = ng2+1`): product scan for `P` (σ=1 over pos
  `0..N-1`; σ=2 over `N-1..0`; padded entries skipped via `where`), `facd`, the `b0` loop
  (`fori_loop` over `i = Nmax-1..1`, masked `i <= N-1`), `A0/B0`, then a forward `lax.scan` over
  positions for the cascade, gathered at `pos`.
* Scatter back: `gnew[(cells),:,σ] += cmul(b0_cell, g_h[(cells),:,σ])` for cells with
  `save_h[σ]` and the per-algorithm conditions; untouched cells unchanged.
* Multi-device: chains couple only `it` within the same `(ik, il, ie, is)`; shard over
  `(is, ie, il)` or `ik` to keep the linked fix local. Sharding over `it` would need an
  all-gather of the four endpoint arrays (that is what `links_p`/`wfb_p` do in Fortran).

### 5.4 Masks summary

`kwork_filter[y,x]`, `ieqzip[y,x]` (pass-through), hybrid `(s,y,l)` (replace by `set_h_zero`),
classes by `l` (passing / wfb / trapped), `forbid[l,ig]`, `is_ttp[l,ig]`, `lbsrc[l,ig]`,
`is_upper_bounce_point` (only in `g_h`), `top/upd[l,ig]`, `l_links==0`, `r_links==0`,
`neighbor[x,y,l]`, `save_h[σ,x,y,l]`, `aky==0`, `g_h[nt,1] != 1`, `g_h[-nt,2] != 1`.

----------------------------------------------------------------------------------------------

## 6. Pitfalls

6.1 **FMA contraction (verified).** The shipped `bin/gs2` contracts (`fmadd/fmsub/fnmsub` in
    `invert_rhs`, `apply_linked_boundary_conditions`, `init_homogeneous_g`; §0.2). XLA:CPU also
    contracts, but with a *different* pattern: on random inputs gfortran-FMA computes
    `re = fma(ar,br,-(ai*bi))`, XLA:CPU `re = fma(-ai,bi,ar*br)` (`im = fma(ar,bi,ai*br)` for
    both) and their complex products disagree in ~27 % of cases; their Smith divisions happened
    to agree on the test, but nothing guarantees this inside larger fused expressions. XLA:GPU
    is expected to contract as well (LLVM NVPTX / ptxas allow FMA fusion by default) —
    unverified, open question 2. Recommended validation setup:
    (a) build a reference GS2 (and the export program) with `F90FLAGS += -ffp-contract=off`;
    (b) run JAX with the strict helpers (§3.0), which match such a build bit for bit on XLA:CPU
    (§A.2); (c) production runs use fast helpers and are compared at ~1e-12 relative.
    Without (a)/(b), expect 1-ulp differences per operation that are then amplified by the
    divisions in 6.12.
6.2 **Smith division.** Complex `/` is Smith's algorithm with the strict `|br| < |bi|` test
    (§3.0). `numpy` complex division and naive `(a*conj(b))/|b|²` do **not** match. Occurs in
    `ainv` (twice), trapped β, trapped-wfb β, self-periodic β, mixed-wfb `facd`, and the GATED
    off-grid/exponential BCs.
6.3 **Operation orders differ between nearly identical expressions**: `get_adjust` ends
    `*facphi*nonmaxw_corr`, `self_periodic`/linked `dadj` end `*nonmaxw_corr*fphi`; `set_h_zero`
    uses pre-multiplied `phi_factor = (zt*nonmaxw)*(-fphi)` and `bpar_factor = 2*(-fbpar)`;
    the mixed-wfb `argl` associates as `((e*al)/bmag)*kperp2` while `init_bessel` uses
    `e*al*kperp2/bmag`. Port each literally; export the wfb end Bessel values.
6.4 **`vpar(ig,1)` for both σ** in `ainv/r/a/b` (`:1882`); `wdrift(ig,σ)` is σ-dependent
    (Coriolis sign). `wdriftttp` is the *uncentred* grid value, `wdrift` the centred one.
6.5 **Coefficient zeroing**: only `r`, `ainv` are zeroed where `forbid(ig) or forbid(ig+1)`;
    `a`, `b` are not. `ainv(lb-1,1) = 1 + 0`; ttp overrides `ainv`, `a`, `r` (not `b`). All four
    are 0 at `ig = ntgrid`. The zeroing block runs only if `nlambda > ng2`.
6.6 **`kwork_filter` pass-through.** In box mode `kwork_filter(1,1)=.true.` is set *after* the
    initial response matrix (`fields.f90:376`); during the time loop the `(kx=0,ky=0)` cell keeps
    its previous `gnew` (no zeroing, source is 0, linked fix skipped). On a dt-reset the flag stays
    `.true.` while the implicit response matrix is rebuilt (see `01_time_loop.md` §6.11).
6.7 **`ieqzip` is ignored by the linked fix**: a skipped cell keeps stale `gnew` but still
    receives `b0*g_h` and feeds its stale endpoints to its neighbours. Reason for GATE.
6.8 **Passing linked sum is not Horner** (§3.11.2) — each cell forms its own nearest-first sum;
    the mixed-wfb cascade *is* a forward recurrence (§3.11.3). Implementing either the other way
    breaks round-off agreement.
6.9 **Trapped β loop** starts at `ntgrid-1`: `ig = ntgrid` is never corrected, and an allowed
    region touching `ig = ntgrid` (bmax not at the domain ends) would get `β = 0`. Likewise
    `g2(ntgrid,2)` is never used. Replicate literally; do not "fix".
6.10 **ttp in the σ=1 sweep** writes `gnew(ig,1)`, not `gnew(ig+1,1)`; the next point keeps its
    initial value (0 / `g`). The source is modified in place at `lb-1` (`:5619`); treat the
    modified source as local to the invert step.
6.11 **Self-links for ky=0**: in linked mode `itleft = itright = it` for `ik=1`, giving
    `l_links = r_links = 0` but `neighbor = .true.` (so a ky=0 mixed wfb is self-periodic only via
    `speriod_flag`). `aky(ik) == 0.0` is an exact comparison everywhere.
6.12 **Ill-conditioned divisions**: `1 - g_h(ub,1)` (trapped β), `1 - g_h(nt,1)` (self-periodic,
    ≈ `1 - Π(-r)` → small for slow passing particles at ky=0), `1 - Π h` (mixed wfb `facd`).
    These amplify upstream last-bit differences; when comparing at 1e-12 compare the final
    `gnew`, and expect the worst agreement at small `|v_||`, deeply trapped near-bounce points.
6.13 **Signed zeros**: masked updates with `β = 0`, `-(x) + 0`, `1*h` etc. can flip the sign of
    zero relative to Fortran. Compare with `==`/tolerances, never bit patterns.
6.14 **`fbpar` default is −1.0** (`run_parameters.f90:82`) when `beta ≠ 0` and `fbpar` unset;
    then `has_bpar = .false.` and `bparnew ≡ 0`, but `fbpar = -1` still multiplies (producing
    `-0.0` terms). Export the effective value.
6.15 **Private module data**: see the accessibility warning in §2.4; the export program needs a
    hook in `dist_fn`.
6.16 **Two calls per step**: BC values and `dadj` use `phinew/bparnew`, which differ between the
    two `timeadv` calls of a step; never cache them across calls. `g_h`, `r`, `ainv` and the
    chain tables are reusable across calls and steps until dt changes.
6.17 **`reset .and. immediate_reset`** returns from `timeadv` *before* `invert_rhs`
    (`dist_fn.fpp:4066`), leaving `gnew` untouched.
6.18 **Periodic BCs are on g, except the mixed wfb**: `self_periodic` makes `g` periodic for
    passing/passing-wfb (`dadj = 0`) but `h` periodic for the mixed wfb when `nonad_zero`.
6.19 **Linked mixed-wfb end cells get no incoming BC at all** (not `-adjust`); the chain is closed
    cyclically. Passing end cells get the `-adjust` BC.
6.20 **No reductions** → no MPI/OpenMP ordering issues in this subsystem (§0.2). `PRODUCT` and
    the `b0` loops are sequential in the order written (no vectorised reduction without
    `-ffast-math`).

----------------------------------------------------------------------------------------------

## 7. Open questions

1. **Reference build**: may the validation reference (gs2 + export program) be built with
   `-ffp-contract=off`? Production numbers of the shipped binary will then differ from the
   reference at the ~1e-16 relative level per operation.
2. **XLA:GPU contraction**: does `lax.optimization_barrier` prevent FMA formation inside a fused
   GPU kernel as it does on XLA:CPU (§A.2), and what does strict mode cost? Needs a test on the
   target GPU.
3. **Export hook**: who adds the accessor for the private `dist_fn` arrays (a source change,
   outside this read-only phase), and should it also expose the wfb end Bessel factors
   (§3.11.4)?
4. **Upstream ownership**: confirm with `02_source_term.md` that the dt-free inputs needed to
   rebuild `vpar/wdrift/wdriftttp` bitwise are in the export list (prefix products before
   `code_dt` is applied).
5. **`kwork_filter` during dt-reset response rebuilds** (implicit fields): intended GS2
   behaviour or bug? The JAX port must replicate whichever the fields spec decides (§6.6).
6. **Geometries with bmax not at ±ntgrid** (trapped particles allowed at the domain end): any
   supported input where this happens? The Fortran algorithm silently skips the bounce matching
   there (§6.9).
7. **Linked BCs on non-box grids** (`ik = 1` with `aky ≠ 0` self-link): GATE proposed; confirm
   no target study needs it.
8. **Complex storage convention** in the netCDF export (two real variables vs trailing `ri`
   dimension) — must match the other specs.
9. Should `start_from_previous_solution`, `boundary_off_grid`, `exponential_boundary`,
   `def_parity` be supported later? All formulas are in §3.8/§3.12; the only extra machinery is
   the `incoming_links` gather and the parity map.

----------------------------------------------------------------------------------------------

## Appendix A — Verification performed for this spec

### A.1 Disassembly
`otool -tv objs/dist_fn.o`: `invert_rhs_1` and `init_invert_rhs` are inlined
(`___dist_fn_MOD_invert_rhs.constprop.0`, `___dist_fn_MOD_init_dist_fn_level_3`); FMA counts as
in §0.2; 16 `fabs` and 49 `fdiv` in `invert_rhs.constprop.0` (Smith's branch test), no
`__divdc3`/`__muldc3` library calls in the object.

### A.2 Arithmetic experiments (scratch programs, not part of the repo)
* gfortran 15.2 `-O3 -mcpu=apple-m1 -fdefault-real-8 -fdefault-double-8 -ffp-contract=off`:
  `a/b` and `1.0/b` equal the §3.0 Smith formula bitwise (20 000/20 000); `a*b` and `c - a*b`
  equal the naive formulas bitwise.
* Same flags without `-ffp-contract=off`: 3 783/20 000 `a/b` and 410/20 000 `1/b` differ from
  plain Smith; `a*b` matches `re=fma(ar,br,-(ai*bi))`, `im=fma(ar,bi,ai*br)`.
* JAX (x64, CPU): `x/y` equals gfortran-FMA Smith on this test; `x*y` uses
  `re=fma(-ai,bi,ar*br)`; `z - x*y` (real) is contracted; wrapping the product in
  `lax.optimization_barrier` removes the contraction; the strict `cmul`/`cdiv` of §3.0 match the
  `-ffp-contract=off` gfortran results bitwise (20 000/20 000).
* gfortran `-ffp-contract=off`, verbatim Fortran statements vs the component/helper forms of this
  spec (Python, no FMA), 20 000 random cases each, 0 mismatches: `ainv, r, a, b` and the ttp
  `ainv, a` (§3.3); `get_adjust` and the `dadj` ordering (§3.8, §3.10); the sweep step
  `-g*r + s*h` (§3.7); `beta = (g - s - adj)/(1.0 - h)` (§3.10).

### A.3 Proposed export variables (this subsystem)
G-shape complex: `ainv, r, a, b, g_h` (initial dt). G-shape real: `wdrift, wdriftttp, vpar`
(02). `(S,E,L,Y,X,T)`: `vperp2, aj0, aj1`. `(L,T)` bool: `forbid, is_bounce_point,
is_lower_bounce_point, is_upper_bounce_point, is_ttp`; `(L,)`: `can_be_ttp`. `(Y,X)` int:
`itleft, itright, l_links, r_links, supercell_labels, n_cells`; `(Y,X)` bool: `kwork_filter,
ieqzip`. `(S,E,L,Y,X,2)` bool: `save_h`; `(S,E,L,Y,X)` bool: `neighbor` (checks). `(S,E,Y,X)`:
`wfb_aj0l, wfb_aj1l, wfb_aj0r, wfb_aj1r`; `(S,E)`: `nonmaxw_corr, vperp2left, vperp2right`.
Scalars: `ntgrid, ntheta0, naky, nlambda, negrid, nspec, ng2, lmax, boundary_option_switch,
no_connections, passing_wfb, trapped_wfb, mixed_wfb, nonad_zero, zero_forbid,
start_from_previous_solution, boundary_off_grid, exponential_boundary,
exponential_boundary_factor, def_parity, even, wfb, fphi, fbpar`; per species `zt, tz, zstm,
smz, bess_fac, fexp, bkdiff, type`; `aky (Y,)`, `delthet (T,)`. For regression tests also dump
`gnew` after `invert_rhs_1` and after `apply_linked_boundary_conditions` for one step (both
calls), together with the `g_work` that was fed in.
