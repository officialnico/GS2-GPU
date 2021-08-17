# GS2 Test Suites 

There are two test suites in GS2 - a legacy test suite and a pfunit test
suite. In addition there is a suite of performance benchmarks.

Please note that some of the tests unfortunately have hard-coded paths
to other files, for example input files or equilibrium files. This
means that test executables can generally only be run from the
location where they get built. If you use a machine that uses
cross-compilation, where executables built on the front-end can only
be run on the back-end (for example, MARCONI), the easiest thing to do
is to copy the whole GS2 directory to the work partition in your job
script, and run the tests as part of the job.

## Legacy test suite

The legacy test suite performs unit tests and integrated tests. These are found
in the folders `unit_tests`, `linear_tests`, `nonlinear_tests`, and subfolders.

To run these tests, do
```
$ make tests
```
or `make unit_tests`, `make linear_tests`, `make nonlinear_tests` to run a
subset of the tests.

To run the tests for a single subfolder only, set the `TESTS` variable, e.g.
```
$ make tests TESTS=fields_local
```
Multiple tests can be run using a space-separated string, e.g.
```
$ make tests TESTS="le_grids theta_grid"
```
Tests can also be _excluded_ by setting the `EXCLUDE_TESTS` variable.

### Performance benchmarks

The performance benchmarks are in the folder `benchmarks`.

To run these, do
```
$ make benchmarks
```
By default this runs the benchmarks with small resolutions, chosen so that the benchmarks
each take a few seconds on a small number of cores (<= 4). The input files for these
have the suffix `level_1`. To run larger benchmarks, do
```
$ make benchmarks BENCHMARK_LEVEL=2
```
These benchmarks are designed to take a few minutes each on 24 to 96 cores. These use the
input files that have the `level_2` suffix.

Benchmarks can be run for a range of processor counts by specifying a string in
`NTESTPROCS`, e.g.
```
$ make benchmarks NTESTPROCS='1 2 4'
```
Note: for this to work, the variable `TESTEXEC` must be set so that it can accept a
variable number of processors, e.g. the system makefile should define
```
TESTEXEC=mpirun -np $$(1)
```
or the command line should contain the definition
```
make benchmarks NTESTPROCS='1 2 4' TESTEXEC='mpirun -np $1'
```
Some benchmarks only make sense for a single processor and will ignore `NTESTPROCS`.

As with the test suite above, particular benchmarks can be selected or ignored by setting
the `TESTS` and `EXCLUDE_TESTS` variables.

## pFUnit tests

GS2 also has unit testing implemented in [pFUnit](https://github.com/Goddard-Fortran-Ecosystem/pFUnit).
These tests are in the folder `pfunit_tests` and subfolders. To run these
tests, set the variable `PFUNIT_DIR` to the location of your pFUnit install
(usually the value of the variable `PFUNIT`), and make the target `pfunit_tests`:
```
$ export PFUNIT_DIR=path/to/install/pfunit
$ make pfunit_tests
```
