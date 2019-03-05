#!/bin/bash
DIR=${1:-"./tests"}
RESULTS_NAME=${2:-"results_of_test.txt"}

# First check that at least one test has passed/skipped
NPASS=$(grep -E "^(passed|skipped):" $(find ${DIR} -name ${RESULTS_NAME}) | wc -l)

if [ ${NPASS} == 0 ]
then
   echo "No passing/skipped tests found -- check compilation achieved."
   exit -1
fi

# Now know tests did run -- check if any failed
NFAIL=$(grep -E "^FAILED:" $(find ${DIR} -name ${RESULTS_NAME}) | wc -l)
exit ${NFAIL}
