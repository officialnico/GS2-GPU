#!/bin/bash
DIR=${1:-"./tests"}
RESULTS_NAME=${2:-"results_of_test.txt"}

NFAIL=$(grep -E "^FAILED:" $(find ${DIR} -name ${RESULTS_NAME}) | wc -l)
exit ${NFAIL}
