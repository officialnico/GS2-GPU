#!/bin/bash

echo "Starting build and tests at $(date)"
git clone git@bitbucket.org:gyrokinetics/utils.git
git clone git@bitbucket.org:gyrokinetics/makefiles.git Makefiles
# Checkout the correct Makefiles branch
cd Makefiles
git checkout automated-testing
cd ..
# Checkout the utils branch which knows about git information
cd utils
git checkout git_info
cd ..
export GK_SYSTEM=bitbucket
export MAKEFLAGS=-IMakefiles DEBUG=on
export TESTEXEC='mpirun -np 8 --oversubscribe'
chmod +x tests/write_test_xml.sh
mpif90 --version
make depend
make -j
make tests
