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
page_dir: ./pages
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

@note This documentation is currently a work in progress -- if you
find any errors please raise an
[issue][issues]

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

Please see the [README][readme] for details on building and running GS2.

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
