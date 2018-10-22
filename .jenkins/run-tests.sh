#!/bin/bash

echo "Starting tests at $(date)"
export GK_SYSTEM=jenkins
export MAKEFLAGS=-IMakefiles DEBUG=on
export TESTEXEC='mpirun -np 6'

module load gcc/7.2.0 openmpi/1.10.0 fftw

# Manually installed netcdf/netcdff

export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:/home/jenkins/workspace/Gyrokinetics/test_job/software/gcc-7.2.0/netcdf/4.6.1/lib:/home/jenkins/workspace/Gyrokinetics/test_job/software/gcc-7.2.0/netcdf-fortran/4.4.4/lib

echo \$LD_LIBRARY_PATH

make tests
