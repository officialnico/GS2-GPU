Initialization and advance time benchmark
-----------------------------------------

This folder runs two **separate** benchmarks with the same executable.
The executable separately measures the time taken for

1. a complete gs2 initialization; and
2. `nstep` gs2 time steps (default: `nstep=100`).

The timings are stored in `timing` files whose roots are `time_initialization`
and `time_advance` respectively.

**Note.** To benchmark the initialization time only, set `nstep=0`.

These benchmarks are run using a single executable for efficiency. We cannot
measure the time for complete gs2 timesteps without first initializing the
code. This initialization can be costly for larger benchmarks, so it makes
sense to merge these two tests.
