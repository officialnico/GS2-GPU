#!/bin/bash
: '
 Script for merging the source branch into the target branch.
'

# Make sure branches are up to date by checking out and pulling
git fetch
git checkout ${sourceBranch}
git pull
git checkout ${targetBranch}
git pull

# Now on target branch. 
# Make test branch from this and attempt merge
git checkout -b jobid-${BUILD_NUMBER}-${sourceBranch}-${targetBranch}
git merge ${sourceBranch}

# Ensure submodules are the version that would be used after the PR is merged.
git submodule update --init --recursive
