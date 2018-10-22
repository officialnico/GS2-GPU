#!/bin/bash

echo "Starting tests at $(date)"
export GK_SYSTEM=jenkins
export MAKEFLAGS=-IMakefiles DEBUG=on
export TESTEXEC='mpirun -np 8'
make tests
