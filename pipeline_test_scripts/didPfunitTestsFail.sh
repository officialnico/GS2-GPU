#!/bin/bash
DIR=${1:-"./"}
RESULTS_NAME=${2:-"pfunit.error"}

# First check that the error file exists
NTEST=$(find ${DIR} -name ${RESULTS_NAME} | wc -l)

if [ ${NTEST} == 0 ]
then
    echo "Found ${NTEST} tests when looking in ${DIR} for ${RESULTS_NAME}"
    echo "No test results found -- check compilation achieved."
    exit 1
fi

# Now know pfunit tests did run -- check if any failed
NFAIL=$(grep -E "\bEncountered 1 or more failures/errors during testing\b" $(find ${DIR} -name ${RESULTS_NAME}) | wc -l)
echo "Found ${NFAIL} failing tests"
exit ${NFAIL}
