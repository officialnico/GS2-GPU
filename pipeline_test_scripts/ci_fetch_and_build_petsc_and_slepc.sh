#!/bin/bash

set -e

function short_help(){
    cat <<EOF
usage: ./$(basename $0) [options]

OPTIONS:
   -h          Show this message and long help
   -d <VAL>    Specify the install destination (defaults ${DESTINATION})
   -p <VAL>    Specify the petsc verision (defaults ${PETSC_VERSION})
   -s <VAL>    Specify the slepc verision (defaults ${SLEPC_VERSION})
   -f          Forces download and build of packages even if already found
EOF
}

# Handle options
DESTINATION=${PWD}
PETSC_VERSION=3.15.2
SLEPC_VERSION=3.15.1
FORCE=""

#Parse arguments
while getopts "hfd:p:s:" OPTIONS
do
  case ${OPTIONS} in
    h)
      short_help
      exit 0
      ;;
    d)
      DESTINATION=${OPTARG}
      ;;
    p)
      PETSC_VERSION=${OPTARG}
      ;;
    s)
      SLEPC_VERSION=${OPTARG}
      ;;
    *)
      short_help
      exit 1
      ;;
  esac
done

PETSC_ORIGIN=https://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-${PETSC_VERSION}.tar.gz
SLEPC_ORIGIN=https://slepc.upv.es/download/distrib/slepc-${SLEPC_VERSION}.tar.gz

DESTINATION=$(readlink -f ${DESTINATION})

header() {
    echo ""
    echo "###############################################"
    echo " "${1:-"UNNAMED SECTION"}
    echo "###############################################"
    echo ""
}

setup_petsc() {
    header "Getting and building petsc version ${PETSC_VERSION} to ${PETSC_DESTINATION}"
    if [[ -e ${PETSC_DESTINATION} && ${FORCE} == "" ]]
    then
        header "${PETSC_DESTINATION} already exists -- assuming already built."
    else
        time download_and_unpack_petsc
        time configure_and_build_petsc
    fi
    export PETSC_DIR=${PETSC_DESTINATION}
    export PETSC_ARCH="" #Not strictly valid but needed for SLEPc
}

download_and_unpack_petsc() {
    download_petsc && unpack_petsc
}

download_petsc() {
    header "Downloading petsc"
    wget ${PETSC_ORIGIN}
}

unpack_petsc() {
    header "Unpacking petsc"
    tar -xzf petsc-${PETSC_VERSION}.tar.gz
}

configure_and_build_petsc() {
    configure_petsc && build_petsc
}

configure_petsc() {
    header "Configuring petsc"
    (cd petsc-${PETSC_VERSION} && ./configure --with-scalar-type=complex --prefix=${PETSC_DESTINATION})
}

build_petsc() {
    header "Building petsc"
    (cd petsc-${PETSC_VERSION} && (make && make install))
}

remove_petsc_share() {
    header "Removing petsc share"
    rm -rf ${PETSC_DESTINATION}/share
}

setup_slepc() {
    header "Getting and building slepc version ${SLEPC_VERSION} to ${SLEPC_DESTINATION}"
    if [[ -e ${SLEPC_DESTINATION} && ${FORCE} == "" ]]
    then
        header "${SLEPC_DESTINATION} already exists -- assuming already built."
    else
        time download_and_unpack_slepc
        time configure_and_build_slepc
    fi
    export SLEPC_DIR=${SLEPC_DESTINATION}
}

download_and_unpack_slepc() {
    download_slepc && unpack_slepc
}

download_slepc() {
    header "Downloading slepc"
    wget ${SLEPC_ORIGIN}
}

unpack_slepc() {
    header "Unpacking slepc"
    tar -xzf slepc-${SLEPC_VERSION}.tar.gz
}

configure_and_build_slepc() {
    configure_slepc && build_slepc
}

configure_slepc() {
    header "Configuring slepc"
    (cd slepc-${SLEPC_VERSION} && ./configure --prefix=${DESTINATION}/slepc)
}

build_slepc() {
    header "Building slepc"
    (cd slepc-${SLEPC_VERSION} && (make SLEPC_DIR=${PWD} && make SLEPC_DIR=${PWD} install))
}

remove_slepc_share() {
    header "Removing slepc share"
    rm -rf ${SLEPC_DESTINATION}/share
}


# Ensure the destination exists
mkdir -p ${DESTINATION}
PETSC_DESTINATION=${DESTINATION}/petsc
SLEPC_DESTINATION=${DESTINATION}/slepc
time (setup_petsc && setup_slepc)
time remove_petsc_share
time remove_slepc_share

# Write out a short file to source which should
# make petsc and slepc available for building
# and running
cat << EOF > load_petsc_slepc
export PETSC_DIR=${PETSC_DESTINATION}
export SLEPC_DIR=${SLEPC_DESTINATION}
export LD_LIBRARY_PATH=\${LD_LIBRARY_PATH}:\${PETSC_DIR}/lib:\${SLEPC_DIR}/lib
export CPATH=\${CPATH}:\${PETSC_DIR}/include:\${SLEPC_DIR}/include
EOF
