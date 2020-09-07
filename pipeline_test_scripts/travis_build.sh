#!/bin/bash

set -e

build_gs2() {
    make test_make ${EXTRA_MAKE_FLAGS}
    make -j gs2 ${EXTRA_MAKE_FLAGS}
}

make_and_run_pfunit_tests() {
    make build_pfunit_library ${EXTRA_MAKE_FLAGS}
    make -j pfunit_tests ${EXTRA_MAKE_FLAGS}
}

make_and_run_integrated_tests() {
    make tests ${EXTRA_MAKE_FLAGS}
    pipeline_test_scripts/didTestsFail.sh ${PWD}/tests results_of_test.txt
}


if [[ "x${GK_TRAVIS_BUILD}" -ne "x" ]]
then
    time build_gs2
fi

if [[ "x${GK_TRAVIS_PFUNIT}" -ne "x" ]]
then
    time make_and_run_pfunit_tests
fi

if [[ "x${GK_TRAVIS_INTEGRATED}" -ne "x" ]]
then
    time make_and_run_integrated_tests
fi
