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
`FFT_LIB`, `-lfftw3_omp`. We attempt to handle this automatically in
the Makefile, but it may be necessary to modify this if using a
different library (e.g. mkl). One can then simply do:

```
make USE_OPENMP=on gs2
```

to build GS2 with OpenMP support enabled.

## Running with OpenMP

The environment variable `OMP_NUM_THREADS` controls the maximum number
of threads to be used and this maximum is reported during GS2
initialisation. If unset this will default to the number of cores on
the system. It is usually advisable to explicitly set
`OMP_NUM_THREADS` and whilst it's possible to over-subscribe a machine
(i.e. `nproc * omp_num_threads > ncpu`) this will generally
significantly harm performance.

The optimal choice for `OMP_NUM_THREADS` depends both on the problem
size and characteristics of the machine. For example, on Archer2 there
are 16 cores per NUMA region (sharing main memory) but performance
often drops considerably when using larger than 4 threads as groups of
four cores share L3 cache.

If is often recommended to run `ulimit -s unlimited` prior to
launching OpenMP enabled executables. In addition to this one may need
to set the
[OMP_STACKSIZE](https://www.openmp.org/spec-html/5.0/openmpse54.html)
environement variable to ensure that each thread has a sufficient
stack size. Failing to set this can lead to "stack smashing detected"
run time error messages. This is particularly important at higher
`OMP_NUM_THREADS` values when the local problem size is large.

## When should I use OpenMP?

There are generally two motivations for using OpenMP in GS2:

1. You need more memory per core than is available when splitting node
   memory amongst all cores. In this situation one can under-populate
   nodes and OpenMP may act to _mitigate_ the performance degradation
   that this would bring (assuming a fixed number of nodes). This is
   often more typical of low processor count runs, but can also hit
   high processor counts as not all memory consumption is distributed.
2. You've reached the scaling limit and MPI is now accounting for the
   dominant part of the run, particularly collectives. By fixing the
   number of MPI ranks but increasing the number of threads one can
   try to scale a little further by roughly fixing the MPI cost whilst
   continuing to distribute the work a little further. This generally
   won't give very efficient scaling due to Amdahl like behaviour
   coming from large fixed MPI costs. In practice MPI costs (in
   particular point-to-point calls in redistributes) can sometimes
   actually decrease with increases to the number of threads
   (presumably due to the potential to overlap communications) so it
   can be slightly more efficient than one might expect. It is often
   more effective to fix the total number of cores and reduce the
   number of MPI processors in this instance. It is important to note
   that recommended sweetspots refer to the number of MPI processes
   and not the total number of cores.
