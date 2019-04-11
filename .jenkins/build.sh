#!/bin/bash

export MAKEFLAGS=-IMakefiles DEBUG=on QUIET=

make depend
make test_make
make -j
