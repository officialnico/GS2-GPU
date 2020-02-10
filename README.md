# GS2

[![gitlab pipeline status](https://gitlab.com/gyrokinetics/gs2/badges/master/pipeline.svg)](https://gitlab.com/gyrokinetics/gs2/commits/master)
[![Build Status](https://anvil.softeng-support.ac.uk/jenkins/buildStatus/icon?job=Gyrokinetics%2FMerge%2BBuild%2BTests-gfortran-7.2-coverage)](https://anvil.softeng-support.ac.uk/jenkins/job/Gyrokinetics/job/Merge+Build+Tests-gfortran-7.2-coverage/)
[![codecov](https://codecov.io/bb/gyrokinetics/gs2/branch/master/graph/badge.svg)](https://codecov.io/bb/gyrokinetics/gs2)

[TOC]

## Introduction

GS2 is a physics application, developed to study low-frequency
turbulence in magnetized plasma. It is typically used to assess the
microstability of plasmas produced in the laboratory and to calculate
key properties of the turbulence which results from instabilities. It
is also used to simulate turbulence in plasmas which occur in nature,
such as in astrophysical and magnetospheric systems.

To download the source run:

    git clone --recurse-submodules https://bitbucket.org/gyrokinetics/gs2.git

or for git < v2.13:

    git clone --recursive https://bitbucket.org/gyrokinetics/gs2.git

Alternatively, tarballs are available for the 
[latest release](https://bitbucket.org/gyrokinetics/gs2/get/master.tar.gz)
and 
[latest development branch](https://bitbucket.org/gyrokinetics/gs2/get/next.tar.gz)
along with other tags and branches on the 
[Bitbucket downloads page](https://bitbucket.org/gyrokinetics/gs2/downloads/)

## Dependencies

GS2 has no required external dependencies, but does have the following
optional ones:

- MPI
- NetCDF
- HDF5
- FFTW
- LAPACK
- NAG

## Developing GS2

Please see the [CONTRIBUTING.md](CONTRIBUTING.md) document for help
and guidelines on contributing to GS2.

## Help and Documentation

Auto-generated source code documentation is currently available at
[https://gyrokinetics.gitlab.io/gs2/](https://gyrokinetics.gitlab.io/gs2/).
Further details are available on the
[wiki](https://bitbucket.org/gyrokinetics/gs2/wiki/).

Please report any bugs by 
[raising an issue on Bitbucket](https://bitbucket.org/gyrokinetics/gs2/issues?status=open).

You can also get support by contacting users and developers on our 
[Slack channel](https://gyrokinetics.slack.com).

## Important note on sub-modules

We use git submodules to deal with the sub-projects `utils` and
`Makefiles` that are shared with other gyrokinetics
projects. Therefore, when cloning, you will find that the `utils` and
`Makefiles` directories are empty unless you clone using the
`--recursive` (git < v2.13) or `--recurse-submodules` (git >= v2.13)
options. alternatively, you can run the following after the initial
checkout:

    cd gs2
    make submodules

See this nice summary on stackoverflow for further details:
https://stackoverflow.com/questions/3796927/how-to-git-clone-including-submodules/4438292#4438292

## Compilation

There are two approaches to compiling GS2:

1. using the `build_gs2` script
2. using `make` directly

### Compilation using `build_gs2`

The level of difficulty of compiling depends on how much your system is
supported.

To use the automated build system try:

    ./build_gs2 -o ls

This will print a list of systems that are supported by the automated build
system. This configures your environment and then executes `make` automatically.
If your system is supported you can simply execute:

    ./build_gs2 -s <system>

Where your `<system>` is replaced by your system name. If your system is not
supported try

    ./build_gs2 -o update

which will try and update the build system. Then try the previous instructions
again.

If your system is still not supported, you can try the older process described
below and in more detail on the wiki:

http://gyrokinetics.sourceforge.net/wiki/index.php/GS2:_A_Guide_for_Beginners#Building_GS2

However, if your system is similar to one of the already supported systems, you could
try copying that system Makefile in the folder Makefiles and editing it to make it work.

Failing that, please raise a support request here:

https://bitbucket.org/gyrokinetics/gs2/issues

### Compilation using `make`

You must set two environment variables to compile agk/gs2: `GK_SYSTEM` and `MAKEFLAGS`.
Since system-specific settings are written in `Makefiles/Makefile.$(GK_SYSTEM)`,
you should specify what system you are on. There are many supported
systems -- look in `Makefiles` to see if yours is already in there!

In order to load `Makefile.$(GK_SYSTEM)` from the `Makefiles/`
directory, you must tell `make` that it should search in that
directory for files. This can be done with the command line option

    make -I Makefiles

or, alternatively, use the `MAKEFLAGS` environment variable to give
the `-I` option.

The `GK_SYSTEM` and `MAKEFLAGS` environment variables can be set in
your `.cshrc` or `.bashrc`. Add the following:

(tcsh)

    setenv MAKEFLAGS -IMakefiles
    setenv GK_SYSTEM 'your system'

(bash)

    export MAKEFLAGS=-IMakefiles
    export GK_SYSTEM='your system'

For utils and geo directories in gs2, do as follows,

    make -I.. -I../Makefiles

After setting those environment variables, you just type `make` to get
an executable. The default target is automatically chosen depending on
where you are making. Other common targets are: `test_make`, `clean`,
`distclean`, `tar`, and `depend`.

We have prepared some compile-time switches. To invoke them, add `SWITCH=value`
to your make command. For example, if you want to use the debug flags:

    make DEBUG=on

The following is the full list of switches the possible options in
`()` brackets and the default in `[]` square brackets:

 - `DEBUG` (defined [undefined]): turns on debug mode
 - `TEST` (defined [undefined]): turns on test mode (not working yet)
 - `PROF` (gprof,ipm [undefined]): turns on profiling mode
                   can set gprof and ipm at the same time
                   ipm is working on limited system. if not
                   supported this option is simply ignored
 - `OPT` (defined,aggressive [defined]): optimization
 - `STATIC` (defined, [undefined]): prevents linking with shared libraries
 - `DBLE` (defined [defined]): promotes precisions of real and complex
 - `USE_MPI` (defined [defined]): turns on distributed memory parallelization
                  with MPI
 - `USE_SHMEM` (defined [undefined]): turns on SHMEM parallel communications
                                  on SGI (not working yet)
 - `USE_FFT` (fftw [fftw]): FFT library to use
 - `USE_NETCDF` (defined [defined]): uses NETCDF library or no
 - `USE_HDF5` (defined [undefined]): uses HDF5 library or no
 - `USE_MDSPLUS` (defined [undefind]): uses MDSplus for geo (gs2 only)
 - `USE_C_INDEX` (defined [undefined]): uses function pointer
 - `USE_POSIX` (defined [undefined]): uses posix routines for argument processing
 - `USE_LOCAL_RAN` (mt [undefined]): uses local random number generator
 - `USE_LOCAL_SPFUNC` (defined [undefined]): use local special functions
     (this overrides compiler's intrinsic special functions and
     NAG library's)
 - `USE_NAGLIB` (spfunc [undefined]): uses NAG library
                          Currently, only used for special functions
                  (NOTE: NAG_PREC=dble (or sngl) must be
                  chosen appropriately for your NAG library)
 - `USE_GCOV` (defined [undefined]): generates gcov notes files
 - `NO_TEST_COLOURS` (defined [undefined]): do not use terminal colours when 
                  reporting test results
 - `KEEP_PROCESSED` (defined [undefined]): if defined then mark processed source files
                  as PRECIOUS so they are not automatically deleted. This can be helpful
                  for debugging the build process and can provide a small optimisation
                  if generating dependency file and then recompiling.

Values that each switch takes are shown in brackets, where "defined"
("undefined") means whether the switch is defined (or not), but not
having that value. (Technically speaking, `ifdef` statement is used in
the makefile.) Thus, note that `DEBUG=off` means `DEBUG` is defined
because this flag is just checked if it is defined or not. Be careful!

In the square bracket, default values in the main makefile (`Makefile`)
are given. Some of those default values are overwritten depending on
hosts, compilers, and also environment variables. This implies that,
for some systems where "module" is provided, those default values are
set properly depending on what modules are loaded.

### More details in CPPFLAGS

Here, the description of full `CPPFLAGS` available are given. However,
`CPPFLAGS` are basically set by above make switches and users need not
be aware of the exact `CPPFLAGS`.

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
| `USE_C_INDEX`   | `_NONE_`      | When invoked, function pointers in layouts_      |
|                 |               | indices.c is used.  A little bit of performance  |
|                 |               | upgrade.  For layouts_indices.c, there is        |
|                 |               | another CFLAG '-DNO_UNDERSCORE' for the          |
|                 |               | resolution of Fortran-C linking problem.         |
|                 |               | The CFLAG '-DUNDERSCORE=2' can also be used      |
|                 |               | when 2 underscores are needed.                   |
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

The standard .F90 extension is not used since some file systems do not
distinguish the case difference of the file names, eg. Windows or Mac.
In order to make your emacs recognize `.fpp` extension as a f90 source code,
you may add the following line to your .emacs file:

    (setq auto-mode-alist (append '(("\\.fpp$" . f90-mode)) auto-mode-alist))

### Porting

Users of a new system needs to prepare `Makefile.$(GK_SYSTEM)` in Makefiles/
directory. `GK_SYSTEM` can be whatever you like, but we recommend to use hostname.
In this file, most simply, you just needs to set your compiler
and to include your compiler's setting as written in Makefile.simple.
If you want to get parallelized version, override `FC` by `mpif90` etc after
loading `Makefile.$(COMPILER)`. You may also want to set libraries places.
Define `FFT_INC`, `FFT_LIB`, `NETCDF_INC`, `NETCDF_LIB` etc in this file.
Write other your host's specific setting if necessary.

Makefile.local is another recommended file to write your local settings.
This file is always included if exists.

### Dependency

The dependency list is given in Makefile.depend, which is auto-generated
by the depend target using the perl script "fortdep". The script is
distributed for developers convenience. Users do not need to reconstruct
it.

## Tests

Users should run the GS2 tests after installing, updating or developing the
code. Tests are also run automatically on code that is pushed to the
Bitbucket repository. Please report any problems or bugs by [creating an
issue on Bitbucket](https://bitbucket.org/gyrokinetics/gs2/issues?status=new&status=open).

For details on running the tests, see the README in the `tests` folder.

## Citing GS2

If you use GS2 in your work, please cite GS2 and acknowledge the grants that
support its development.

### Citation
GS2's DOI is [10.5281/zenodo.2551066](https://doi.org/10.5281/zenodo.2551066).
The DOI for the current stable version 8.0.3 is [10.5281/zenodo.3656334](https://doi.org/10.5281/zenodo.3656334).
The full citation for GS2 may be also found in CITATION.cff. You can convert the
CITATION.cff file into a Bibtex file as follows:

    pip3 install --user cffconvert
    cffconvert -if CITATION.cff -f bibtex -of CITATION.bib 

### Acknowledgement

Please acknowledge grants which support the development of GS2:

> "This work made use of computational support by CoSeC, the Computational
> Science Centre for Research Communities, through CCP Plasma ([EP/M022463/1](https://gow.epsrc.ukri.org/NGBOViewGrant.aspx?GrantRef=EP/M022463/1))
> and HEC Plasma ([EP/R029148/1](https://gow.epsrc.ukri.org/NGBOViewGrant.aspx?GrantRef=EP/R029148/1))."
