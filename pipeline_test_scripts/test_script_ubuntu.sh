#!/bin/bash

echo "Starting build and tests at $(date)"
git submodule update --init --recursive
export GK_SYSTEM=jenkins
export MAKEFLAGS=-IMakefiles DEBUG=on
export TESTEXEC='mpirun -np 8'
chmod +x tests/write_test_xml.sh
mpif90 --version
make depend
make test_make
make -j
make tests
