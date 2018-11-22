#!/bin/bash

export NETCDF_VERSION=4.6.1
export NETCDFF_VERSION=4.4.4
export MY_NETCDF_DIR=${WORKSPACE}/software/
export MY_NETCDF_DIR=software/${params.compiler}/netcdf/
export MY_NETCDF_DIR=\$NETCDF_VERSION
export MY_NETCDFF_DIR=${WORKSPACE}/software/${params.compiler}/netcdf-fortran/\$NETCDFF_VERSION
export MY_NETCDF_INC=\$MY_NETCDF_DIR/include
export MY_NETCDFF_INC=\$MY_NETCDFF_DIR/include
export MY_NETCDF_LIB=\$MY_NETCDF_DIR/lib
export MY_NETCDFF_LIB=\$MY_NETCDFF_DIR/lib
export GK_SYSTEM=${params.GK_SYSTEM}
export MAKEFLAGS=-IMakefiles
if [ params.opt == "debug" ] 
then
	echo "Compiling in DEBUG mode"
	export MAKEFLAGS=\$MAKEFLAGS DEBUG=on
fi				
export TESTEXEC='mpirun -np 6'
module load ${params.modules}
export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:\$MY_NETCDF_LIB:\$MY_NETCDFF_LIB:/opt/software/${params.compiler}/fftw/3.3.5/lib

make depend
make test_make
make -j
