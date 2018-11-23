#!/bin/bash
export MAKEFLAGS=-IMakefiles
if [ params.opt == "debug" ] 
then
	echo "Compiling in DEBUG mode"
	export MAKEFLAGS=\$MAKEFLAGS DEBUG=on
fi				
export TESTEXEC='mpirun -np 6'

make depend
make test_make
make -j
