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
- python (for tests and documentation)

## Developing GS2

Please see the [CONTRIBUTING.md](CONTRIBUTING.md) document for help
and guidelines on contributing to GS2.

## Help and Documentation

Documentation, include auto-generated source code documentation is currently available at
[https://gyrokinetics.gitlab.io/gs2/](https://gyrokinetics.gitlab.io/gs2/).

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

## Building GS2

GS2 has been built on a variety of existing systems. To build on one of the
existing supported systems, set `GK_SYSTEM`, and provide the path to the
directory containing the makefiles. For example, to build GS2 on the UK national
supercomputer, ARCHER, run:

    make -I Makefiles GK_SYSTEM=archer

For further details on building GS2, including how to build it on a new system,
see the documentation in the [user manual][docs_building].

## Citing GS2

If you use GS2 in your work, please cite GS2 and acknowledge the grants that
support its development.

### Citation
GS2's DOI is [10.5281/zenodo.2551066](https://doi.org/10.5281/zenodo.2551066).
The DOI for the current stable version 8.0.5 is [10.5281/zenodo.4012004](https://doi.org/10.5281/zenodo.4012004).
The full citation for GS2 may be also found in CITATION.cff. You can convert the
CITATION.cff file into a Bibtex file as follows:

    pip3 install --user cffconvert
    cffconvert -if CITATION.cff -f bibtex -of CITATION.bib 

### Acknowledgement

Please acknowledge grants which support the development of GS2:

> "This work made use of computational support by CoSeC, the Computational
> Science Centre for Research Communities, through CCP Plasma ([EP/M022463/1](https://gow.epsrc.ukri.org/NGBOViewGrant.aspx?GrantRef=EP/M022463/1))
> and HEC Plasma ([EP/R029148/1](https://gow.epsrc.ukri.org/NGBOViewGrant.aspx?GrantRef=EP/R029148/1))."


[docs_building]: https://gyrokinetics.gitlab.io/gs2/page/user_manual/index.html
