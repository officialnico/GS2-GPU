#!/bin/bash

echo "Starting tests at $(date)"
export TESTEXEC='mpirun -np 8'
make tests
