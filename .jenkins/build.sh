#!/bin/bash

export NETCDF_VERSION=4.6.1
export NETCDFF_VERSION=4.4.4
export MY_NETCDF_DIR=${WORKSPACE}/software/$1/netcdf/\$NETCDF_VERSION
export MY_NETCDFF_DIR=${WORKSPACE}/software/$1/netcdf-fortran/\$NETCDFF_VERSION
export MY_NETCDF_INC=\$MY_NETCDF_DIR/include
export MY_NETCDFF_INC=\$MY_NETCDFF_DIR/include
export MY_NETCDF_LIB=\$MY_NETCDF_DIR/lib
export MY_NETCDFF_LIB=\$MY_NETCDFF_DIR/lib
export GK_SYSTEM=$2
export MAKEFLAGS=-IMakefiles
if [ params.opt == "debug" ] 
then
	echo "Compiling in DEBUG mode"
	export MAKEFLAGS=\$MAKEFLAGS DEBUG=on
fi				
export TESTEXEC='mpirun -np 6'
module load $3 $4 $5
export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:\$MY_NETCDF_LIB:\$MY_NETCDFF_LIB:/opt/software/$1/fftw/3.3.5/lib

make depend
make test_make
make -j
