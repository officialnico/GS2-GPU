# GS2 Test Suites 

There are two test suites in GS2 - a legacy test suite and a pfunit test
suite.

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

## pFUnit tests

GS2 also has unit testing implemented in [pFUnit](https://github.com/Goddard-Fortran-Ecosystem/pFUnit).
These tests are in the folder `pfunit_tests` and subfolders. To run these
tests, build GS2 with the variable `USE_PFUNIT=on` and make the target
`pfunit_tests`:
```
$ make USE_PFUNIT=on
$ make pfunit_tests USE_PFUNIT=on
```
