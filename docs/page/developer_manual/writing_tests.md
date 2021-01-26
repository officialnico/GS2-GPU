---
title: Writing Tests
---

[TOC]

There are two test suites in GS2 - a legacy test suite and a [pFUnit][pfunit] test
suite. The pFUnit test suite is a collection of _unit tests_. These try to test the
smallest units of code -- usually individual subroutines or functions. They also aim to be
very fast: "one-tenth of a second is _too slow_" is the usual motto. Unit tests should
also try to be independent of other tests: you don't want something going wrong in one
test to cause other tests to unnecessarily fail.

This document will take you through writing a new pFUnit test, detailing some of the
pitfalls and things to be aware of.

## Basic structure of a test

pFUnit can automatically detect and run new tests added to the test suite. In order to do
this, it uses a custom preprocessor to convert the `.pf` files to `.F90` which are then
compiled by the Fortran compiler as usual. We do also need to tell pFUnit about new files
however, and these are listed in `tests/pfunit_tests/testSuites.inc`.

A very simple pFUnit test for GS2 looks something like this:

`my_tests.pf`:

```
module test_basic_mod
  use PFUNIT_MOD
  implicit none

contains
  @test
  subroutine test_addition
    @assertEqual(4, 2 + 2)
  end subroutine test_addition
end module test_basic_mod
```

`testSuites.inc`:

```
ADD_TEST_SUITE(test_basic_mod_suite)
```

Note that the module name has an additional `_suite` appended in `testSuites.inc`.

Let's break down the test module line by line:

```
module test_basic_mod
```

The exact name isn't important, but in the GS2 tests, we prefer the naming convention
`test_<GS2 module name>_mod`; for example: `test_leq_mod`, or `test_gs2_time_mod`.

```
  use PFUNIT_MOD
```

After the `.pf` to `.F90` conversion step, the test module will use various functions and
types from the pFUnit module, so we `use` the entire module.

The all-caps is important here: we currently support using either pFUnit 3.x or 4.x, which
use different module names. We preprocess the resulting `.F90` to replace `PFUNIT_MOD`
with the correct module name.

```
  implicit none
  
contains
```

In this basic test module we just have the one test procedure, but we could declare
module-level parameters here. For example, you may wish to set a `parameter` defining a
tolerance for equalities.

```
  @test
  subroutine test_addition
```

The `@test` is a pFUnit directive declaring the immediately following subroutine as a test
to be run automatically. The subroutine name is not important, but we prefer the naming
convention of a `test_` prefix to distinguish the procedure from non-test procedures in
the same module.

```
    @assertEqual(4, 2 + 2)
```

This is the basic pFUnit testing directive that checks for equality between its two
arguments. If they don't compare equal, pFUnit will flag the test as having failed and
report it at the end of the run. Other asserts include `@assertTrue`, `@assertFalse`,
and `@assertSameShape`.

Importantly, pFUnit will stop running _this particular_ test if an `assertEqual` fails. If
the arguments are arrays, pFUnit will report only the first index where they differ.

The order of the arguments doesn't particularly matter, but when reporting errors pFUnit
assumes the **first** argument is the _expected_ value, and the **second** argument is the
_actual_ value.

Something to be aware of here is that the pFUnit preprocessor doesn't like this directive
split over multiple lines! Currently, we just have to live with long lines.

## Getting fancier

### Setup and Teardown

After writing a few tests, you may find that there's repeated code setting things up
first, and then cleaning up afterwards. pFUnit allows you to declare subroutines to be run
before and after each test in a module:

```
@before
subroutine setup
  ! Do something before every test
end subroutine setup

@after
subroutine teardown
  ! Do something after every test
end subroutine teardown
```

The subroutine names are conventional, but `setup` and `teardown` are preferred.

**Important!** If you mark one subroutine with either `@before` or `@after`, you **must**
also have the corresponding one, even if it's just an empty routine.

These `@before` and `@after` routines are run before/after each individual test, even if
that test fails. This can be very useful to, for example, call the `<module>_finish`
routine to cleanup a GS2 module, which helps make the tests more independent. Without the
`@after` routine, a failure in one test might result in a GS2 module being in a bad or
wrong state for the next test. For this reason, you may wish to call the `<module>_finish`
routine in both the setup and teardown routines.

### Tolerances

By default, pFUnit compares values exactly. This is usually not what you want when dealing
with floating point numbers. pFUnit can check equality within a tolerance, although
unfortunately, as of v4.1.9, only absolute tolerances (`abs(actual-expected) <=
tolerance`). Tolerances can be specified with the `tolerance` keyword to `@assertEqual`:

```
@assertEqual(3.14, pi, tolerance=1e-2)
```

### Message

When a test contains multiple asserts, it can be tricky to tell exactly which assert
caused a test to fail. Unfortunately, while pFUnit prints the values it's comparing, it
doesn't print the _expressions_. Instead, it has a mechanism to give a bit more
information with the `message` keyword:

```
@assertEqual(3.14, pi, message="comparing pi")
```

### MPI tests

pFUnit has the ability to _parameterise_ a test over the number of MPI ranks. That is, it
can rerun a test with different numbers of MPI processes. This can be useful to check that
a function works with different ranks.

Additionally, many routines often end up needing MPI to be set up because they contain
some code like `if (proc0) then ...`. For these routines to work correctly when called
from a pFUnit test, we essentially need to parameterise the test over one process.

The form of an MPI test is a little different:

```
@test(npes=[1, 2])
subroutine test_some_mpi(this)
  class(MPITestMethod), intent(inout) :: this
  ...
```

Here, we'll end up calling the test twice times, first with one rank, and second with
two. Next, our test subroutine now takes a single argument of `class(MPITestMethod)` --
note it is `intent(inout)`. Conventionally we name it `this`, but the exact name is not
important. This argument has some useful methods:

```
this%getMpiCommunicator()    ! Get the MPI_COMM for this test run
this%getNumProcesses()       ! Wrapper for MPI_Comm_size
this%getProcessRank()        ! Wrapper for MPI_Comm_rank
```

Parameterised tests, including MPI tests, **must** have setup/teardown routines. The setup
routine **must** call `mp::init_mp`, but the teardown **must not** call `mp::finish_mp`:

```
@before
subroutine setup(this)
  use mp, only : init_mp
  class (MpiTestMethod), intent(inout) :: this
  ! Make sure GS2 variables like proc0 are populated correctly:
  call init_mp(this%getMpiCommunicator())
end subroutine setup

@after
subroutine teardown(this)
  class (MpiTestMethod), intent(inout) :: this
  ! No call to finish_mp!
end subroutine teardown
```

Calling `finish_mp` will finalise MPI, which will cause problems for the next MPI test!

When using MPI tests for more than just getting `proc0` correct in GS2, you should use
`@mpiAssert*` in your tests.




[pfunit]: https://github.com/Goddard-Fortran-Ecosystem/pFUnit
