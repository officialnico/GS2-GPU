#!/bin/bash

git fetch
git checkout ${sourceBranch}
git pull

# Until Jenkins branch merged with next, force the targetBranch to be
# jenkins, to ensure all necessary files exist
git checkout origin/${targetBranch}

git checkout -b jobid-${BUILD_NUMBER}-${sourceBranch}-${targetBranch}
git merge ${sourceBranch}
git submodule update
