# Breaking changes in v8.1

GS2 development policy is to minmize disruption to users by only introducing
breaking changes into major and minor releases (e.g. v8, v8.1, but not
v8.1.1). By "breaking change", we mean anything that will change the physics
results without the user changing anything in their input file. This typically
means changes to defaults, or improvements to physics/algorithms that aren't
fixing a definite bug.

The breaking changes for v8.1 are:

* [PR 353 Ensure maximum timestep stored in restart](https://bitbucket.org/gyrokinetics/gs2/pull-requests/353/). 
With this change, the maximum time step (input variable `dt0`) is stored in
the restart file and is restored when `delt_option = 'check_restart'`. The
value used can still be overridden by setting `dt0`.  This breaks backwards
compatibility with pre-8.1 restart files: now if `delt_max` is not found in
the restart file then we would fail to set the maximum time step or the time
step history from the restart file, and an error is thrown.
