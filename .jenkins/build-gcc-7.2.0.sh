#!/bin/bash

echo "Starting build at $(date)"
export GK_SYSTEM=jenkins
export MAKEFLAGS=-IMakefiles DEBUG=on

module load gcc/7.2.0 openmpi/1.10.0 
module load fftw

# Add path to manually installed netcdf/netcdff
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/jenkins/workspace/Gyrokinetics/test_job/software/gcc-7.2.0/netcdf/4.6.1/lib:/home/jenkins/workspace/Gyrokinetics/test_job/software/gcc-7.2.0/netcdf-fortran/4.4.4/lib:/opt/software/gcc-7.2.0/fftw/3.3.5/lib

cd Makefiles/
git checkout jenkins
cd ..
make depend
make test_make
make -j

