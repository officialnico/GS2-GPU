#!/bin/bash

echo "Starting checkout at $(date)"
# Jenkins keeps old repos. Delete and start afresh
rm -rf gs2
git clone https://bitbucket.org/gyrokinetics/gs2.git"
git submodule update --init --recursive
git checkout jenkins
