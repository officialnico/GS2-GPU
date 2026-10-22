#!/usr/bin/env bash
# Build the GS2 reference binaries bin/gs2 and bin/gs2_export.
#
# Contract: gpu/spec/00_decisions.md D2 and gpu/spec/09_export_infra.md §9.
#
# Usage (from anywhere):
#   gpu/scripts/build_reference.sh            # reference build: stock flags + -ffp-contract=off
#   gpu/scripts/build_reference.sh --stock    # stock build: exactly the flags of today's bin/gs2
#
# Options:
#   --stock        build with the original (stock) flags
#   --jobs N       parallel make jobs (default 8)
#   --clean        force a full rebuild even if the mode did not change
#   --gs2-only     build bin/gs2 only (no bin/gs2_export)
#   --dry-run      print the commands, do not run them
#
# Exact commands this script runs (repository root, single MPI rank builds, no OpenMP):
#
#   STOCK (the command that builds today's bin/gs2):
#     make -j8 -I Makefiles GK_SYSTEM=macosx-homebrew FFT_INC=-I/opt/homebrew/include \
#          FFT_LIB="-L/opt/homebrew/lib -lfftw3" NETCDF_INC="-I/opt/homebrew/include" \
#          NETCDF_LIB="-L/opt/homebrew/lib -lnetcdf -lnetcdff" depend
#     make -j8 -I Makefiles GK_SYSTEM=macosx-homebrew FFT_INC=-I/opt/homebrew/include \
#          FFT_LIB="-L/opt/homebrew/lib -lfftw3" NETCDF_INC="-I/opt/homebrew/include" \
#          NETCDF_LIB="-L/opt/homebrew/lib -lnetcdf -lnetcdff" gs2 gs2_export
#
#   REFERENCE (default): the same two commands with one extra command-line variable
#          F90OPTFLAGS="-O3 -ffp-contract=off"
#     F90OPTFLAGS is the Fortran optimisation-flag variable of the GS2 build: Makefile:149 defines
#     it empty, Makefiles/enable_features.make appends $(OPT_FLAG) (= -O3 for gnu-gfortran,
#     Makefiles/Makefile.gnu-gfortran:24) and Makefile:547 appends it to F90FLAGS (and hence to
#     LDFLAGS). Setting it on the command line therefore yields exactly the stock F90FLAGS with
#     "-O3" followed by "-ffp-contract=off"; every other flag (-fdefault-real-8, ..., and the
#     driver defaults such as -mcpu=apple-m1 reported by --build-config) is unchanged.
#
# Mode switches force a full rebuild: objects are shared in objs/ and include/, so the script
# records the mode in objs/.gs2_gpu_build_mode and runs `make clean` (objs/*.o, include/*.mod)
# whenever the requested mode differs from the recorded one (or nothing is recorded).
# Note: a plain `make` run by hand in between is not tracked; use --clean if in doubt.
#
# `make depend` is always run first: the GPU-port edits add the module gs2_export_probe to the
# use lists of physics modules and add the new exporter files (09 §9.2).
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

MODE=reference
JOBS=8
FORCE_CLEAN=0
TARGETS="gs2 gs2_export"
DRY=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --stock) MODE=stock ;;
    --reference) MODE=reference ;;
    --jobs) JOBS="$2"; shift ;;
    --jobs=*) JOBS="${1#--jobs=}" ;;
    --clean) FORCE_CLEAN=1 ;;
    --gs2-only) TARGETS="gs2" ;;
    --dry-run) DRY=1 ;;
    -h|--help) sed -n '2,40p' "$0"; exit 0 ;;
    *) echo "build_reference.sh: unknown option '$1'" >&2; exit 2 ;;
  esac
  shift
done

# The make command line of the stock build (do not edit without updating the header above).
MAKE_ARGS=(-I Makefiles GK_SYSTEM=macosx-homebrew
           FFT_INC=-I/opt/homebrew/include
           "FFT_LIB=-L/opt/homebrew/lib -lfftw3"
           "NETCDF_INC=-I/opt/homebrew/include"
           "NETCDF_LIB=-L/opt/homebrew/lib -lnetcdf -lnetcdff")
REF_FLAGS="-O3 -ffp-contract=off"
if [[ "${MODE}" == reference ]]; then
  MAKE_ARGS+=("F90OPTFLAGS=${REF_FLAGS}")
  MODE_STRING="reference F90OPTFLAGS=${REF_FLAGS}"
else
  MODE_STRING="stock"
fi

run() {
  printf '+'; printf ' %q' "$@"; printf '\n'
  if [[ "${DRY}" -eq 0 ]]; then "$@"; fi
}

cd "${ROOT}"
STAMP="${ROOT}/objs/.gs2_gpu_build_mode"
OLD_MODE="$(cat "${STAMP}" 2>/dev/null || echo unknown)"
echo "build_reference.sh: mode='${MODE_STRING}', recorded objects mode='${OLD_MODE}'"
if [[ "${FORCE_CLEAN}" -eq 1 || "${OLD_MODE}" != "${MODE_STRING}" ]]; then
  echo "build_reference.sh: full rebuild (mode change or --clean)"
  run make "${MAKE_ARGS[@]}" clean
  if [[ "${DRY}" -eq 0 ]]; then rm -f "${STAMP}" bin/gs2 bin/gs2_export; fi
fi

run make "${MAKE_ARGS[@]}" depend
run make -j"${JOBS}" "${MAKE_ARGS[@]}" ${TARGETS}

if [[ "${DRY}" -eq 1 ]]; then exit 0; fi
mkdir -p objs
echo "${MODE_STRING}" > "${STAMP}"

# Verify the flags that actually went into the binaries.
FLAGS_LINE="$(bin/gs2 --build-config | sed -n 's/^Build flags: //p')"
case "${FLAGS_LINE}" in
  *-ffp-contract=off*) BUILT=reference ;;
  *) BUILT=stock ;;
esac
if [[ "${BUILT}" != "${MODE}" ]]; then
  echo "build_reference.sh: ERROR: bin/gs2 reports a ${BUILT} build, expected ${MODE}" >&2
  exit 1
fi
if [[ " ${TARGETS} " == *" gs2_export "* ]]; then
  EFLAGS="$(bin/gs2_export --build-config | sed -n 's/^Build flags: //p')"
  if [[ "${EFLAGS}" != "${FLAGS_LINE}" ]]; then
    echo "build_reference.sh: ERROR: bin/gs2 and bin/gs2_export report different build flags" >&2
    exit 1
  fi
fi
echo "build_reference.sh: OK (${MODE}); build flags: ${FLAGS_LINE}"
