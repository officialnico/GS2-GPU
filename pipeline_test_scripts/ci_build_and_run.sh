#!/bin/bash

set -e
ulimit -s unlimited

header() {
    echo ""
    echo "###############################################"
    echo " "${1:-"UNNAMED SECTION"}
    echo "###############################################"
    echo ""
}

add_to_local_makefile() {
    header "Adding to Makefile.local"
    echo -e ${GK_MAKEFILE_LOCAL_EXTRA_FLAGS} >> Makefile.local
    cat Makefile.local
}

build_gs2() {
    header "Build ${GK_BUILD_TARGET}"
    make test_make ${EXTRA_MAKE_FLAGS}
    make -j ${GK_BUILD_TARGET} ${EXTRA_MAKE_FLAGS}
    ./bin/gs2 --version
    ./bin/gs2 --build-config
}

make_and_run_pfunit_tests() {
    header "Build and run pfunit tests"
    make build_pfunit_library ${EXTRA_MAKE_FLAGS} ${EXTRA_PFUNIT_FLAGS}
    if [[ "x${GK_CI_PFUNIT_COVERAGE}" == "x" ]]
    then
        make -j pfunit_tests ${EXTRA_MAKE_FLAGS} ${EXTRA_PFUNIT_FLAGS}
    else
        make -j pfunit_coverage ${EXTRA_MAKE_FLAGS} ${EXTRA_PFUNIT_FLAGS}
    fi
    pipeline_test_scripts/didPfunitTestsFail.sh ${PWD}/ pfunit.error
}

make_and_run_integrated_tests() {
    header "Build and run integrated tests"
    make tests ${EXTRA_MAKE_FLAGS} ${EXTRA_TEST_FLAGS}
    pipeline_test_scripts/didTestsFail.sh ${PWD}/tests results_of_test.txt
}


if [[ "x${GK_MAKEFILE_LOCAL_EXTRA_FLAGS}" != "x" ]]
then
    time add_to_local_makefile
fi

if [[ "x${GK_CI_BUILD}" != "x" ]]
then
    time build_gs2
fi

if [[ "x${GK_CI_PFUNIT}" != "x" ]]
then
    time make_and_run_pfunit_tests
fi

if [[ "x${GK_CI_INTEGRATED}" != "x" ]]
then
    time make_and_run_integrated_tests
fi
