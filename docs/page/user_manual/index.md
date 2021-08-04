---
title: User manual
---

[TOC]

First steps
===========

So, you've never run GS2 before? This is the place to start. On this page we cover the
basics of how to get GS2 installed and running.

There are two different sorts of complexity in scientific computing; the complexity of the
physical problem in hand, and the additional complexity brought about by the realities of
software, including that of compiling, running and processing results. Most of the
complexity of GS2 derives from the former; hard work by many people means that the second
is much less of an issue. Now since the physical problem in hand, that of kinetic
turbulence in various magnetic geometries, is incredibly complicated, this still means
that to fully understand and use all the power of GS2 is a fiendishly tortuous
undertaking. But do not be discouraged. In simple physical situations, GS2 is simple to
use.

Installing GS2
==============

Downloading GS2
---------------

To download the source run:

``` shell
git clone --recurse-submodules https://bitbucket.org/gyrokinetics/gs2.git
```

or for git < v2.13:

```shell
git clone --recursive https://bitbucket.org/gyrokinetics/gs2.git
```

Alternatively, tarballs are available for the
[latest release](https://bitbucket.org/gyrokinetics/gs2/get/master.tar.gz)
and
[latest development branch](https://bitbucket.org/gyrokinetics/gs2/get/next.tar.gz)
along with other tags and branches on the
[Bitbucket downloads page](https://bitbucket.org/gyrokinetics/gs2/downloads/)

Please also see the note in the [README][readme] on sub-modules.

Building GS2
------------

### Dependencies

GS2 has the following dependencies, which should be installed before
building. On most HPC systems all of these will be installed and you can
skip this step.

- [FFTW3](http://www.fftw.org/) a popular open source package for
  performing fast Fourier transforms. Required for non-linear runs,
  possible to compile without these for linear-only simulations,
  though not recommended.
- [NetCDF](http://www.unidata.ucar.edu/software/netcdf/) a set of
  software libraries for data storage and retrieval.
- [HDF5](http://www.hdfgroup.org/HDF5/) - another set of libraries,
  needed if you want to install the parallel version of NetCDF.
- MPI is required to run on more than one processor. There are several
  MPI implementations, such as [Open MPI](http://www.open-mpi.org/).

### Compiling

There are two approaches to compiling GS2:

1. using the `build_gs2` script
2. using `make` directly

### Compilation using `build_gs2`

The level of difficulty of compiling depends on how much your system is supported.

To use the automated build system try:

``` shell
./build_gs2 -o ls
```

This will print a list of systems that are supported by the automated build system. This
configures your environment and then executes `make` automatically. If your system is
supported you can simply execute:

``` shell
./build_gs2 -s <system>
```

Where your `<system>` is replaced by your system name. If your system is not supported try

``` shell
./build_gs2 -o update
```

which will try and update the build system. Then try the previous instructions again.

If your system is still not supported, you can try the older process described below.
However, if your system is similar to one of the already supported systems, you could try
copying that system Makefile in the folder Makefiles and editing it to make it work.

Failing that, please [raise a support request][newissue].

### Compilation using `make`

`build_gs2` is essentially just a wrapper around the `Makefiles`, and they can still be
used manually.

To compile GS2, you must set two environment variables: `GK_SYSTEM` and `MAKEFLAGS`. Since
system-specific settings are written in `Makefiles/Makefile.$(GK_SYSTEM)`, you should
specify what system you are on. There are many supported systems -- look in `Makefiles` to
see if yours is already in there!

In order to load `Makefile.$(GK_SYSTEM)` from the `Makefiles/` directory, you must tell
`make` that it should search in that directory for files. This can be done with the
command line option

    make -I Makefiles

or, alternatively, use the `MAKEFLAGS` environment variable to give the `-I` option.

For example, to compile GS2 on the UK national supercomputer, ARCHER, you would run:

    make -I Makefiles GK_SYSTEM=archer

Instead of setting them on the command line, you can also set the `GK_SYSTEM` and
`MAKEFLAGS` environment variables in your `.cshrc` or `.bashrc`. Add the following:

(tcsh)

    setenv MAKEFLAGS -IMakefiles
    setenv GK_SYSTEM 'your system'

(bash)

    export MAKEFLAGS=-IMakefiles
    export GK_SYSTEM='your system'

After setting those environment variables, you just type `make` to get an executable. The
default target is automatically chosen depending on where you are making. Other common
targets are: `test_make`, `clean`, `distclean`, `tar`, and `depend`.

We have prepared some compile-time switches. To invoke them, add `SWITCH=value`
to your make command. For example, if you want to use the debug flags:

    make DEBUG=on

The following is the full list of switches the possible options in
`()` brackets and the default in `[]` square brackets:

 - `DEBUG` (defined [undefined]): turn on debug mode
 - `PROF` (`gprof`,`ipm` [undefined]): turns on profiling mode. It's possible to set
   `gprof` and `ipm` at the same time.  `ipm` works only on a few systems. Any other value
   is ignored.
 - `OPT` (defined,`aggressive` [defined]): turn on optimization
 - `STATIC` (defined, [undefined]): prevents linking with shared libraries
 - `DBLE` (defined [defined]): promotes precisions of `real` and `complex` to "double"
   (64-bit)
 - `USE_MPI` (defined [defined]): turns on distributed memory parallelization with MPI
 - `USE_SHMEM` (defined [undefined]): turns on SHMEM parallel communications on SGI (not working yet)
 - `USE_FFT` (`fftw`, `fftw3` [fftw]): which FFT library to use
 - `USE_NETCDF` (defined [defined]): uses NETCDF library or not
 - `USE_HDF5` (defined [undefined]): uses HDF5 library or not
 - `USE_MDSPLUS` (defined [undefind]): uses MDSplus for geo (gs2 only)
 - `USE_POSIX` (defined [undefined]): uses posix routines for argument processing
 - `USE_LOCAL_RAN` (`mt` [undefined]): uses local random number generator
 - `USE_LOCAL_SPFUNC` (defined [undefined]): use local special functions (this overrides
   compiler's intrinsic special functions and NAG library's)
 - `USE_NAGLIB` (spfunc [undefined]): uses NAG library. Currently, only used for special
   functions (NOTE: `NAG_PREC=dble` (or `sngl`) must be chosen appropriately for your NAG
   library)
 - `USE_GCOV` (defined [undefined]): enable coverage measurement
 - `NO_TEST_COLOURS` (defined [undefined]): do not use terminal colours when reporting
   test results
 - `KEEP_PROCESSED` (defined [undefined]): if defined then mark processed source files as
   `PRECIOUS` so they are not automatically deleted. This can be helpful for debugging the
   build process and can provide a small optimisation if generating dependency file and
   then recompiling.

Values that each switch takes are shown in brackets, where "defined" ("undefined") means
whether the switch is defined (or not), but not having that value. (Technically speaking,
the `ifdef` statement is used in the makefile.) Thus, note that `DEBUG=off` means `DEBUG`
is defined because this flag is just checked if it is defined or not. Be careful!

The default values in square brackets are those in the main makefile (`Makefile`). Some of
those default values are overwritten depending on hosts, compilers, and also environment
variables. This implies that, for some systems where `module` is provided, those default
values are set properly depending on what modules are loaded.

### More details in MACRO_DEFS

Here, the description of full `MACRO_DEFS` available are given. However,
`MACRO_DEFS` are basically set by above make switches and users need not
be aware of the exact `MACRO_DEFS`.

In the values column, the top line is the value chosen for the undefined
macro (`-D{variable}` is NOT invoked).  The second line corresponds to the
DEFAULT value when the variable is invoked but the value is not specified
(`-D{variable}` only).  The rest of choices will become available only when
the value is specified by `-D{variable}=_value_`.  Note that no space should
be placed within each definition (`-D{value} = _value_` does not work as one
might expect).

| variable        | values        | description                                      |
|-----------------|---------------|--------------------------------------------------|
| `FCOMPILER`     | `_XL_`        | Fortran compiler choices                         |
|                 | `_GFORTRAN_`  |                                                  |
|                 | `_G95_`       |                                                  |
|                 | `_NAG_`       |                                                  |
|                 | `_INTEL_`     |                                                  |
|                 | `_PGI_`       |                                                  |
|                 | `_PATHSCALE_` |                                                  |
|                 | `_LAHEY_`     |                                                  |
|                 |               |                                                  |
| `FFT`           | `_NONE_`      | We can choose to not compile fftw.               |
|                 | `_FFTW_`      | At the moment only fftw is the available         |
|                 | `_FFTW3_`     | fft option.                                      |
|                 |               |                                                  |
| `MPI`           | `_NONE_`      | Do not use MPI                                   |
|                 |               | Use of -DMPI makes mpi work.                     |
|                 |               |                                                  |
| `NETCDF`        | `_NONE_`      | Do not use NetCDF                                |
|                 |               | Use of -DNETCDF makes netcdf work.               |
|                 |               |                                                  |
| `HDF`           | `_NONE_`      | Do not use HDF5                                  |
|                 |               | Use of -DHDF makes hdf work.                     |
|                 |               |                                                  |
| `NO_SIZEOF`     | `_NONE_`      | This is for memory accounting.                   |
|                 |               | In most compilers, non-standard sizeof function  |
|                 |               | works fine as C function.                        |
|                 |               | For a compiler which does not support sizeof     |
|                 |               | function, one can use hard-wired fortran         |
|                 |               | implementation of size_of function instead       |
|                 |               | by specifying -DNO_SIZEOF.                       |
|                 |               |                                                  |
| `RANDOM`        | `_NONE_`      | Uses compiler's intrinsic random                 |
|                 | `_RANMT_`     | Uses Mersenne Twister 19937 random               |
|                 |               |                                                  |
| `SPFUNC`        | `_NONE_`      | Uses compiler's intrinsic special functions      |
|                 | `_LOCAL_`     | Uses local special functions                     |
|                 | `_NAG_`       | Uses NAG library's special functions             |
|                 |               | (NOTE: for compiles not having intrinsic         |
|                 |               | special functions like XL fortran, `_LOCAL_` or  |
|                 |               | `_NAG_` should be chosen. Current default for XL |
|                 |               | is `_LOCAL_`)                                    |
|                 |               |                                                  |
| `USE_L2E_MAP`   | `_NONE_`      | Affects collision.fpp.  When invoked, l2e_map    |
| (agk only)      |               | becomes available, which describes direct map    |
|                 |               | from lorentz to ediffuse layouts.                |
|                 |               | Requires additional memory, but might enhance    |
|                 |               | performance --- not much, though.                |

The standard `.F90` extension is not used since some file systems do not distinguish the
case difference of the file names, such as Windows or Mac OS. If you use Emacs, you may
add the following line to your `.emacs` file to make it recognize the `.fpp` extension as
F90 source code:

``` emacs-lisp
(setq auto-mode-alist (append '(("\\.fpp$" . f90-mode)) auto-mode-alist))
```

### Porting

Users of a new system need to prepare a new `Makefile.$(GK_SYSTEM)`. `GK_SYSTEM` can be
whatever you like, but we recommend to use the hostname of the machine. At its simplest,
you just need to set your compiler and to include your compiler's setting as written in
[Makefile.simple][simple_makefile]:

``` makefile
COMPILER = intel
include Makefile.$(COMPILER)
```

In order to build the parallelised version of GS2, you will need to override the `FC`
variable with the appropriate MPI wrapper for your system, for example `mpif90`. Do this
after the include of the compiler Makefile.

You may also need to set the location of the various libraries. Define `FFT_INC`,
`FFT_LIB`, `NETCDF_INC`, `NETCDF_LIB` in this file. You can also write other host-specific
configuration in this file.

### Local settings

If you wish to use an existing `GK_SYSTEM` Makefile but add or change options, you can use
a file named `Makefile.local` at the top-level. This file is always included if exists.

### Dependency

The dependency list is given in `Makefile.depend`, which is auto-generated by the depend
target using the perl script `fortdep`.

## Tests

Users should run the GS2 tests after installing, updating or developing the code. To run
the tests, run:

    make tests

Don't forget to include `-I Makefiles` and set `GK_SYSTEM` if you haven't set them as
environment variables!

If you have any problems running the tests, please report them by [creating an issue on
Bitbucket][newissue].

Tests are also run automatically on code that is pushed to the Bitbucket repository.

For more details on running the tests, see the [README][tests_readme] in the `tests` folder.

Running GS2
===========

Making an Input File
--------------------

A GS2 input file is a series of namelists in which various parameters
are specified:

```Fortran
&namelist
    parameter1 = value1 ! Fortran style comments are allowed
    parameter2 = value2
/
```

Documentation for the various parameters is [available here](../namelists). However, they
are so many and varied that making your own input file from scratch is not advised if
you're just starting out. It's much better to take an existing input file and tweak it to
suit your own needs. A list of example input files is here: [Example Input
Files](Example_Input_Files).

Take the input file, edit it as you need and save it as `run_name.in`, where `run_name` is
a name you have chosen for this simulation.

Running GS2
-----------

The command for running GS2 is simply:

    $ gs2 run_name.in

This will run GS2 on one processor, which, for all but the lowest resolutions, will take a
very long time. If you are on your own system with an MPI implementation installed, you'll
probably want to issue a command like

    $ mpirun -n <nprocs> gs2 run_name.in

where `nprocs` is the number of processors. If you are running on an HPC system, you will
need to consult their own manuals.

Restarting a simulation
-----------------------

By default, GS2 doesn't create the files needed to restart a simulation, as they can get
very large. In order to create them, you must first set the input variable
[gs2_diagnostics_knobs::save_for_restart][save_for_restart] to `T` (or `.true.`) in your
input file:

``` Fortran
&gs2_diagnostics_knobs
    save_for_restart = T
/
```

This will create one restart file per processor named like `<your input filename>.nc.<N>`
where `<N>` is the processor rank. There are two other inputs you can use to control where
these restart files get created and how they're named:

- you can use [init_g_knobs::restart_file][restart_file_option] to change the "base name"
  of the restart files. By default, this is the same as the input filename.
- you can change the directory where the restart files are created with
  [init_g_knobs::restart_dir][restart_dir_option]. `restart_dir` defaults to the directory
  where you run `gs2`. If you set this option, it's important to create the directory
  first, as GS2 won't do this for you. You can turn on
  [gs2_diagnostics_knobs::file_safety_check][file_safety_check_option] (default is on) to
  check during initialisation that `restart_dir` exists and that GS2 can write to the
  directory.

To then restart a simulation, set the following input variables:

- [knobs::delt_option][delt_option] to `check_restart`
- [init_g_knobs::ginit_option][ginit_option] to `restart`

``` Fortran
&knobs
    delt_option = "check_restart"
/
&init_g_knobs
    ginit_option = "restart"
/
```

Important: by default, GS2 will overwrite existing output files! You may wish to backup
these files before restarting the simulation. It is possible to get GS2 to append the
output to the existing files with [gs2_diagnostics_knobs::append_old][append_old]:

``` Fortran
&gs2_diagnostics_knobs
    append_old = T
/
```

Analysing the Results
---------------------

GS2 outputs a myriad collection of data and files. These are all dealt with in [Results
and Diagnostics][results_diagnostics].

A Warning
---------

We finish with a warning, which cannot be taken too seriously:

**Just because GS2 outputs a result, it does not follow that the result
is correct**

There are 1001 tweaks and adjustments that can be made to the resolution, the geometry,
the physical parameters and the diagnostics which may affect the final result. Like
Charles Babbage's adding machine, if you put wrong numbers in, if the resolution is
insufficient, you will not get right answers out. A rough guide to checking whether you
have sufficient resolution is provided in [Results and Diagnostics
page][results_diagnostics]. Please note that, unlike typical explicit algorithms, the
implicit algorithm can return results which look reasonable even when the time step is too
large. As such it is important to also check convergence in time step (usually only a
linear problem).

There is also a program, [[ingen]], which is built with GS2, which will check your input
file for inconsistencies. Although documentation for this is currently sparse, try running

    $ ingen run_name.in

If you're just starting out using GS2, try to find some well-known result that you can
reproduce, before going on to explore new territory.

Good luck!

[newissue]: https://bitbucket.org/gyrokinetics/gs2/issues?status=new&status=open
[readme]: https://bitbucket.org/gyrokinetics/gs2/src/6fbd2a6bc6ac5e5be0e0737f1860a0e0f78026ed/README.md?at=master&fileviewer=file-view-default
[simple_makefile]: https://bitbucket.org/gyrokinetics/makefiles/src/master/Makefile.simple
[tests_readme]: https://bitbucket.org/gyrokinetics/gs2/src/master/tests/README.md
[results_diagnostics]: |url|/page/user_manual/results_diagnostics.html
[append_old]: |url|/page/namelists/index.html#gs2_diagnostics_knobs-append_old
[save_for_restart]: |url|/page/namelists/index.html#gs2_diagnostics_knobs-save_for_restart
[delt_option]: |url|/page/namelists/index.html#knobs-delt_option
[ginit_option]: |url|/page/namelists/index.html#init_g_knobs-ginit_option
[restart_dir_option]: |url|/page/namelists/index.html#init_g_knobs-restart_dir
[restart_file_option]: |url|/page/namelists/index.html#init_g_knobs-restart_file
[file_safety_check_option]: |url|/page/namelists/index.html#gs2_diagnostics_knobs-file_safety_check
[use_old_diagnostics]: |url|/page/namelists/index.html#knobs-use_old_diagnostics

<!-- Local Variables: -->
<!-- mode: gfm -->
<!-- fill-column: 90 -->
<!-- End: -->
