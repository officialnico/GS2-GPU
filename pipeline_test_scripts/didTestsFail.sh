#!/bin/bash
DIR=${1:-"./tests"}
RESULTS_NAME=${2:-"results_of_test.txt"}

# First check that at least one test result available
NTEST=$(find ${DIR} -name ${RESULTS_NAME} | wc -l)

if [ ${NTEST} == 0 ]
then
   echo "No test results found -- check compilation achieved."
   exit -1
fi

# Now know tests did run -- check if any failed
NFAIL=$(grep -E "^FAILED:" $(find ${DIR} -name ${RESULTS_NAME}) | wc -l)
exit ${NFAIL}
