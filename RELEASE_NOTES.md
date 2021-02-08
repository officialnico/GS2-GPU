# GS2 v8.1

GS2 v8.1 is our first **feature release** using the new git + Bitbucket
workflow, and it's a big one! It's over **120 pull requests** - and over **1000
commits** - ahead of the current bugfix release 8.0.6.

v8.1 offers many new **features**, **enhancements** and **performance
improvements**, along with significant expansion of **testing and
documentation**, and lots of **bug fixes**!

This release enables **larger**, **faster**, more **robust simulations**,
without changing the familiar interface for users or developers. But since this
is a feature release, we have allowed some changes that **might require user
changes** to input or analysis scripts, or that **may alter results** slightly
- so do check out of the [breaking changes](#breaking-changes) section below.   

## Features and Enhancements

* **Smart defaults** - it’s now possible to change the default value of a
  variable, based on the input from another variable.  For example, we now set
  `write_fluxes` and `write_fluxes_by_mode` to true for nonlinear runs
  (`nonlin=true`) and to false for linear runs.
* New simpler option for **hyperviscosity**, `simple3D`, which has the form
  `D_hyper3D * (|kperp|/ max |kperp|)^P_hyper3D` that respects continuity along
  the field line, and does not require MPI communication.
* Users can now specify both **spatial grid resolution** `nx`, `ny` and
  **wavenumber resolution** `nakx`, `naky` (provided the anti-aliasing
  constraint is met), making it easier to select **sweetspots** (... but can
  still define only one of these and let GS2 choose the other).
* The **explicit source terms** are now stored in the output. This is the
  default, but can be turned off by setting
  `include_explicit_source_in_restart` in the `init_g_knobs` namelist to false.
* Users can now **specify a random seed** for initialization by setting `seed`
  in the `knobs` namelist. This makes it easier to exactly reproduce runs.
* Users can now control the **netcdf sync frequency** (the frequency at which
  updates are written to disk) by setting `nc_sync_freq` in
  `gs2_diagnostics_knobs`.
* GS2 now **saves the layout string** in output and checks for consistency in
  restarted runs.
* The new "run configuration" type allows **programmatic control of code
  inputs**, allowing better documentation and testing.

## Documentation

We've made significant improvements to our documentation:

* We've expanded our [**user
  manual**](https://gyrokinetics.gitlab.io/gs2/page/user_manual/index.html) and
  [**developer
  manual**](https://gyrokinetics.gitlab.io/gs2/page/developer_manual/index.html).
* We've expanded our **inline comments** and now have [**documentation of input
  variables**](https://gyrokinetics.gitlab.io/gs2/page/namelists/index.html). 
* Source code is parsed and put on the [**documentation
  website**](https://gyrokinetics.gitlab.io) that builds automatically for
  `master` and `next`.

## Testing and infrastructure

We've expanded our test suites and improved our building and continuous
integration frameworks.

* Makefile is now **tidier** and **more robust**. Git **hashes** for `utils`and
  `Makefiles` are now **stored** during the build process, making it easier to
  track down bugs.
* We've removed many instances of imported-but-unused variables from the source
  code, making the **build output quieter** and easier to read. 

## Performance improvements

8.1 is not really a performance release.
We have two big performance improvements coming up in future releases - an
improved field solve, and better data layout - but 8.1 still includes some
smaller improvements:  

* We've vastly **reduce creation and destruction of MPI communicators**,
  improving initialization performance, and enabling larger simulations which
  were previously failing due either to running out of communicators or to
  communicator fragmentation.
* We now **store response matrices for reuse** (for example, when the timestep
  halves and then doubles again). This can save orders-of-magnitude in
  reinitialization time when the timestep is changing frequently.
* We've improved **indexing** by removing "select" cases from tight loops.
* We've **removed** the **less performant** version of some routines (see
  [breaking changes](#breaking-changes)).
* We now **avoid unnecessary initialization** of unused data layouts.

## Bugfixes

Here are a selection of the more notable squashed bugs. We've:

* **tided** the source code, **improved code structure** and **inline
  comments**, and **removed unused code**
* fixed bugs in the `ceq` and `leq` **geometries**
* fixed the **initial theta grid** size for **eikcoefs**
* ensured **FFTs** are now **bitwise reproducible** by default
* added the missing **benchmark option** to the `build_gs2` script
* reduced the `bouncefuzz` parameter. `bouncefuzz` adds a small tolerance to
  the inequality determining whether a pitch angle is forbidden, to account for
  roundoff error in `lambda` and `B`. Reducing `bouncefuzz` makes the calculation
  more accurate, so any changes will actually be a bugfix.  This was introduced
  to prevent a numerical instability. 

## Breaking changes

Here are a list of changes that **might need users to alter their
input/analysis scripts**, or that **might change results**:
 
* `opt_redist_init` is no longer an option - we only support the optimal
  version of the redistribution initialization routine. **Users need** to
  remove `opt_redist_init` from `layout_knobs` namelist in input files.
* For consistency, the new diagnostics now uses lower case for `rplot`,
  `zplot`, `rprime` and `zprime`. This has the potential to **break analysis
  scripts**. 
* **Changes to defaults**: `write_fluxes` and `write_fluxes_by_mode` now
  default to `true` for nonlinear runs and `false` for nonlinear runs.
* Runs **now abort** if the user attempts to load a restart file with a
  different layout than the one set in the rest of the code / input file.
