---
title: Running With OpenMP
---

# OpenMP

OpenMP is a convenient directive based approach to shared-memory work
sharing. In GS2 this is currently restricted to distributing loop
iterations over a team of threads. Only a small subset of loops are
parallelised with OpenMP, however these should mostly represent the
main areas of work in standard simulations.

## Compiling with OpenMP

To compile with OpenMP one can simply pass `USE_OPENMP=on` to make. We
also support the use of threaded FFTW routines (FFTW3 only) so if
building with FFTs then one may need to add an additional library to
`FFT_LIB`, `-lfftw3_omp`. Adding the following snippet to
Makefile.local (or the system specific makefile) should be sufficient:

```
ifdef USE_OPENMP
  FFT_LIB += -lfftw3_omp
endif
```

and one can then do

```
make USE_OPENMP=on gs2
```

to build GS2 with OpenMP support enabled.

## Running with OpenMP

The environment variable `OMP_NUM_THREADS` controls the maximum number
of threads to be used and this maximum is reported during GS2
initialisation. If unset this will default to the number of cores on
the system.

The optimal choice for `OMP_NUM_THREADS` depends both on the problem
size and characteristics of the machine. For example, on Archer2 there
are 16 cores per NUMA region (sharing main memory) but performance
often drops considerably when using larger than 4 threads as groups of
four cores share L3 cache.
