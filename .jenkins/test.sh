#!/bin/bash

export MAKEFLAGS=-IMakefiles DEBUG=on
export TESTEXEC='mpirun -np 2'

make tests
cd ../
rm -rf gs2
