#!/bin/bash

export MAKEFLAGS=-IMakefiles DEBUG=on QUIET=
export EXCLUDE_TESTS="gs2_diagnostics_new cyclone_itg"

make depend
make test_make
make -j
