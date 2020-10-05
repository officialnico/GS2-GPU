#!/bin/bash

set -e

header() {
    echo ""
    echo "###############################################"
    echo " "${1:-"UNNAMED SECTION"}
    echo "###############################################"
    echo ""
}

build_gs2() {
    header "Build GS2"
    make test_make ${EXTRA_MAKE_FLAGS}
    make -j gs2 ${EXTRA_MAKE_FLAGS}
}

make_and_run_pfunit_tests() {
    header "Build and run pfunit tests"
    make build_pfunit_library ${EXTRA_MAKE_FLAGS}
    if [[ "x${GK_TRAVIS_PFUNIT_COVERAGE}" == "x" ]]
    then
        make -j pfunit_tests ${EXTRA_MAKE_FLAGS}
    else
        make -j pfunit_coverage ${EXTRA_MAKE_FLAGS}
    fi
}

make_and_run_integrated_tests() {
    header "Build and run integrated tests"
    make tests ${EXTRA_MAKE_FLAGS} ${EXTRA_TEST_FLAGS}
    pipeline_test_scripts/didTestsFail.sh ${PWD}/tests results_of_test.txt
}


if [[ "x${GK_TRAVIS_BUILD}" != "x" ]]
then
    time build_gs2
fi

if [[ "x${GK_TRAVIS_PFUNIT}" != "x" ]]
then
    time make_and_run_pfunit_tests
fi

if [[ "x${GK_TRAVIS_INTEGRATED}" != "x" ]]
then
    time make_and_run_integrated_tests
fi
