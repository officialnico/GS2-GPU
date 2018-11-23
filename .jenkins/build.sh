#!/bin/bash

export MAKEFLAGS=-IMakefiles DEBUG=on

make depend
make test_make
make -j
