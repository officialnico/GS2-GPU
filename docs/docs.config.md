---
project: GS2
project_bitbucket: https://bitbucket.org/gyrokinetics/gs2
project_download: https://bitbucket.org/gyrokinetics/gs2/get/master.tar.gz
summary:
    <p style="text-align: center">
        <img alt="GS2: A fast, flexible, physicist's toolkit for
            gyrokinetics" src="|media|/GS2_text.png" style="width: 35%;
            min-width: 44%">
        <img alt="GS2 logo"
            src="|media|/MAST_ETG_transparent_cropped.png"
            style="width: 55%; min-width: 55%">
    </p>
author: The GS2 team
author_description:
    Collection of developers from around the world.
    GS2 has been developed by dozens of developers over the years.
    For a complete list, please see <https://doi.org/10.5281/zenodo.2551066>
parallel: 16
src_dir: ../src
css: user.css
src_dir: ../externals/utils
output_dir: ./html
graph_dir: ./graphs
page_dir: ./page
exclude_dir: ../src/templates
media_dir: ./media
favicon: ./media/favicon.ico
md_extensions: markdown.extensions.toc
fpp_extensions: fpp
fpp_extensions: F90
predocmark: >
predocmark_alt: !
docmark: <
docmark_alt: #
display: public
display: protected
display: private
source: true
print_creation_date: true
search: false
graph: false
coloured_edges: true
extra_vartypes: MPI_Comm
                Vec
                Mat
                PetscErrorCode
                PetscInt
                PetscScalar
                PetscReal
                EPS
                EPSType
                EPSExtraction
                ST
---

GS2 is a physics application, developed to study low-frequency
turbulence in magnetized plasma. It is typically used to assess the
microstability of plasmas produced in the laboratory and to calculate
key properties of the turbulence which results from instabilities. It
is also used to simulate turbulence in plasmas which occur in nature,
such as in astrophysical and magnetospheric systems.

### GS2 is:

- open and accessible, including to new PhDs
- an excellent tool for doing physics and developing new improved
  physics calculations
- a great platform for trying out new algorithms
- highly optimised and designed to scale to use as many cores
  as possible

Getting started: [read the user manual](page/user_manual).

@note This documentation is currently a work in progress -- if you
find any errors please raise an
[issue][issues]

## Features

### Fast microstability assessment

Linear microstability growth rates are calculated on a wavenumber-by-wavenumber
basis with an implicit initial-value algorithm in the ballooning (or
"flux-tube") limit. Linear and quasilinear properties of the fastest growing (or
least damped) eigenmode at a given wavenumber may be calculated independently
(and therefore reasonably quickly).

### Fully gyrokinetic, nonlinear simulations

Nonlinear simulations of fully developed turbulence can be performed by users
with access to a parallel computer. All plasma species are treated on an equal,
gyrokinetic footing. Nonlinear simulations provide fluctuation spectra,
anomalous (turbulent) heating rates, and species-by-species anomalous transport
coefficients for particles, momentum and energy.

### Good performance

GS2 is a parallel code which scales well to large numbers of processors. There
are separate optimizations to work with for small (Beowulf-style) clusters and
large supercomputers.

### Flexible Simulation Geometry

Linear and nonlinear calculations may be carried out using a wide
range of assumptions, including:

-  Local slab
-  Local cylinder
-  Local torus
-  Vacuum magnetic dipole
-  High aspect ratio torus (analytic axisymmetric equilibrium)
-  Numerically generated, local equilibrium (Miller-style)
-  Axisymmetric, numerically generated equilibria with arbitrary poloidal shaping from
    - [EFIT](http://web.gat.com/efit)
    - [TRANSP](https://transp.pppl.gov/index.html)
    - [TOQ](http://web.gat.com/toq/)
    - [JSOLVER](http://w3.pppl.gov/rib/repositories/NTCC/catalog/Asset/jsolver.html)

### Efficient computational grid

Turbulent structures in gyrokinetics are highly elongated along the magnetic
field. GS2 uses field-line following (Clebsch) coordinates to resolve such
structures with maximal efficiency, in a flux tube of modest
extent. Pseudo-spectral algorithms are used in the spatial directions
perpendicular to the field, and for the two velocity space coordinate grids
(energy and pitch angle) for high accuracy on computable 5-D grids.


![GS2 scaling](|media|/Gs2_performance.png)


## Get the code

The project is hosted on [BitBucket][bitbucket]. To download the source run:

    git clone --recurse-submodules https://bitbucket.org/gyrokinetics/gs2.git

or download the latest stable version from [here][mastertarball].

## Dependencies

GS2 has no required external dependencies, but does have the following
optional ones:

- MPI
- NetCDF
- HDF5
- FFTW
- LAPACK
- NAG

## Building and running GS2

GS2 has been built on a variety of existing systems. To build on one of the
existing supported systems, set `GK_SYSTEM`, and provide the path to the
directory containing the makefiles. For example, to build GS2 on the UK national
supercomputer, ARCHER, run:

    make -I Makefiles GK_SYSTEM=archer

For further details on building GS2, including how to build it on a new system,
see the documentation in the [user manual](page/user_manual)

## Citing GS2

If you use GS2 please remember to cite the following publications and
the source code. For more details see CITATION.cff provided with GS2,
or see <https://doi.org/10.5281/zenodo.2551066>

## Documentation

This website has documentation built automatically using
[Ford][ford]. See below for details on individual source files,
modules, procedures, and data types:


[bitbucket]: https://bitbucket.org/gyrokinetics/gs2/
[ford]: https://github.com/Fortran-FOSS-Programmers/ford
[issues]: https://bitbucket.org/gyrokinetics/gs2/issues
[mastertarball]: https://bitbucket.org/gyrokinetics/gs2/get/master.tar.gz
[readme]: https://bitbucket.org/gyrokinetics/gs2/src/6fbd2a6bc6ac5e5be0e0737f1860a0e0f78026ed/README.md?at=master&fileviewer=file-view-default
