#!/bin/bash

set -e

function short_help(){
    cat <<EOF
usage: ./$(basename $0) [options]

OPTIONS:
   -y          Proceed with setup
   -f          Force run even if gitlab CI not detected
EOF
}


header() {
    echo ""
    echo "###############################################"
    echo " "${1:-"UNNAMED SECTION"}
    echo "###############################################"
    echo ""
}

function add_intel_oneapi_repo() {
    header "Setup intel oneapi repository"
    dnf config-manager --add-repo https://yum.repos.intel.com/oneapi
    rpmkeys --import https://yum.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB
}

function install_intel_compiler() {
    header "Install intel compiler"
    dnf install -y intel-oneapi-compiler-fortran
}

# Handle options
SETUP_INTEL_ON_GITLAB_FORCE="OFF"
SETUP_INTEL_ON_GITLAB_YES="OFF"

#Parse arguments
while getopts "hfy" OPTIONS
do
  case ${OPTIONS} in
    h)
      short_help
      exit 0
      ;;
    f)
      SETUP_INTEL_ON_GITLAB_FORCE="ON"
      ;;
    y)
      SETUP_INTEL_ON_GITLAB_YES="ON"
      ;;
    *)
      short_help
      exit 1
      ;;
  esac
done

# Crude test to see if we're running on gitlab CI or not
if [ -z ${GITLAB_CI} ]
then
    # We haven't detected gitlab CI so exit unless we've forced execution
    if [ "${SETUP_INTEL_ON_GITLAB_FORCE}" != "ON" ]
    then
        echo "Not detected gitlab CI and not being forced so exiting"
        exit 0
    fi
fi

# If we've reached this point then we are on CI or being forced.
# Before setting things up check if we've asked to proceed. This
# is basically just a second level of "safety" to avoid users adding
# repos to their local machine by mistake etc.
if [ "${SETUP_INTEL_ON_GITLAB_YES}" != "ON" ]
then
    echo "Running on gitlab CI (or forced) but not asked to proceed. Pass -y if you are sure."
    exit 0
fi


time add_intel_oneapi_repo
time install_intel_compiler

# Write out a short file to source which should
# make the intel compiler available

cat << EOF > load_intel_compiler
source /opt/intel/oneapi/compiler/latest/env/vars.sh
. /opt/intel/oneapi/setvars.sh intel64
EOF
