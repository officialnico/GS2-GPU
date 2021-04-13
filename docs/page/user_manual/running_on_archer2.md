---
title: Running on ARCHER2
---

[ARCHER2][archer2] is the UK's top research supercomputer. There are a
few things to be aware of when compiling and running GS2 on ARCHER2
that we will list here.

## Compiling

There is a Makefile already written for ARCHER2, so you can set
`GK_SYSTEM=archer2` to use it. This Makefile uses the GNU compilers by
default, so you must load the GNU compiler suite in addition to the
other GS2 dependencies -- in contrast to previous systems, the
compiler suites use `module restore`:

```
module restore PrgEnv-gnu
module load cray-fftw
module load cray-netcdf
module load cray-hdf5   # The netCDF module needs HDF5
```

then we can compile GS2:

```
make -I Makefiles GK_SYSTEM=archer2 -j
```

## Running

ARCHER2 doesn't mount the home partition on the compute nodes, so
you will want to either make sure the compilation is done in the work
partition, or more simply, copy the `gs2` executable over.

In your job submission script, it is necessary to restore some
modules first to make sure the correct libraries are found. Note that
you only need to load `epcc-job-env` and restore `PrgEnv-gnu`, the
other libraries will be found automatically.

Here is an example job submission script for 1 minute on 1 node:

```shell
#!/bin/bash

# Job info
#SBATCH --job-name=gs2
#SBATCH --nodes=1
#SBATCH --tasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH --time=00:01:00

# Partition info
#SBATCH --account=<account>
#SBATCH --partition=standard

# Ensure simulation has access to libraries.
module load epcc-job-env
module restore /etc/cray-pe.d/PrgEnv-gnu

# Important: set OpenMP threads to 1, as we're not using OpenMP.
# This makes sure any libraries don't use threads
export OMP_NUM_THREADS=1

srun --distribution=block:block --hint=nomultithread gs2 <input file>
```

Note that the `--distribution=block:block --hint=nomultithread`
arguments to `srun` are _strongly_ [recommended by the ARCHER2
team][srun].

If you're developing GS2 and want to make sure your modified code
remains on the backed-up home partition, don't forget to copy the
`gs2` executable on to the work partition first. You can do this as
part of your job script.

### Short Experiments

Sometimes you'd like to run a short simulation and find you spend more
time waiting for it to get through the queue and start running than it
actually takes to run the simulation. In that case, you can try using
the [short queue][short] (which is called "Quality of Service" in
`SLURM` jargon). You can do that by adding the following to your
submission script:

```shell
# The following uses the short job partition, which can only be used for quick
jobs, but has a much sorter queue time
#SBATCH --reservation=shortqos
#SBATCH --qos=short
```

Note the queue limits in the documentation: 20 minutes on 8 nodes, and
only Monday to Friday.

[archer2]: https://www.archer2.ac.uk
[srun]: https://docs.archer2.ac.uk/user-guide/scheduler/#srun-launching-parallel-jobs
[short]: https://docs.archer2.ac.uk/user-guide/scheduler/#quality-of-service-qos

<!-- Local Variables: -->
<!-- mode: gfm -->
<!-- fill-column: 90 -->
<!-- End: -->
