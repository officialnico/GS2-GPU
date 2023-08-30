#!/bin/bash

# A small script to automatically update the golden answers with whatever
# result is currently being produced.

EXE=${1}

# If the exectuable doesn't exist try running make
if [[ ! -x ${EXE} ]]
then
    make
fi

for INP in $(grep -l gamma_expected *.in)
do
    TTT=$(./${EXE} ${INP} | grep -A 1 'growth rate failed' | grep 'but got' | awk '{print $6}')
    if [[ -n ${TTT} ]]
    then
	sed --in-place 's/gamma_expected = .*/gamma_expected = '${TTT}'/g' ${INP}
    fi
done

for INP in $(grep -l gamma_expected *.in)
do
    TTT=$(./${EXE} ${INP} | grep -A 1 'frequency failed' | grep 'but got' | awk '{print $6}')
    if [[ -n ${TTT} ]]
    then
	sed --in-place 's/omega_expected = .*/omega_expected = '${TTT}'/g' ${INP}
    fi
done
