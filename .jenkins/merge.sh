#!/bin/bash

git fetch
git checkout ${sourceBranch}
git pull
git checkout origin/${targetBranch}
git checkout -b jobid-${BUILD_NUMBER}-${sourceBranch}-${targetBranch}
git merge ${sourceBranch}
rm -rf tests/unit_tests/collisions/

