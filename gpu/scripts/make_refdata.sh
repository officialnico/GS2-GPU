#!/usr/bin/env bash
# Generate the reference data of the GS2 JAX port (gpu/spec/00_decisions.md D2, D9;
# gpu/spec/09_export_infra.md; gpu/spec/09a_exporter_notes.md).
#
# Usage:
#   gpu/scripts/make_refdata.sh [CASE ...]
#
#   CASE is a case name (gpu/cases/CASE.in) or a path to a GS2 input file. Without arguments every
#   gpu/cases/*.in present at run time is processed (cases added later are picked up on the next
#   run); if gpu/cases is missing or has no *.in, the GS2 Cyclone tests are used instead:
#   tests/linear_tests/cyclone_itg/cyclone_itg_low_res_base.in, tests/nonlinear_tests/cyclone_itg/*.in
#
# For every case, in gpu/refdata/<case>/ (old contents of the three run directories are replaced):
#   gs2/     full bin/gs2 run (all GS2 outputs, for end-to-end comparison; <case>.out.nc etc.)
#   export/  bin/gs2_export: <case>.gs2export.nc (static data, /state_t0, dt ladder /dt_ladder/k1..),
#            <case>.out.nc = GS2's own static part + istep-0 record (JAX appends from nout = 2), and
#            the GS2 restart file of the t = 0 state (*.gs2export_t0.0, canonical-map check 09 §12 T2)
#   trace/   bin/gs2_export traced for the first N steps: <case>.gs2trace.nc (all probes) and the
#            traced run's GS2 outputs; its own export file is compared with export/ and then removed
#   <case>.gs2export.nc, <case>.gs2trace.nc, <case>.out.nc   symlinks into the directories above
#   check_export.log   output of gpu/scripts/check_export.py;  STATUS   PASS or FAIL (one line)
#
# Environment:
#   GS2_REFDATA_NSTEP_TRACE   traced steps N (default 10)
#   GS2_REFDATA_DT_LADDER     dt-ladder rungs written to the export file (default 2)
#   GS2_REFDATA_COMPRESSION   netCDF deflate level of export/trace files (default 1; lossless)
#   GS2_REFDATA_DIR           output root (default gpu/refdata)
#   GS2_REFDATA_KEEP_TRACE_EXPORT=1   keep trace/<case>.gs2export.nc
#   GS2_REFDATA_BUILD=1       run gpu/scripts/build_reference.sh first
#   GS2_REFDATA_ALLOW_STOCK=1 accept binaries built with stock flags (default: refuse)
#   GS2_REFDATA_ONLY=export   regenerate only export/ (keep gs2/ and trace/), then re-run the checks
#   MPIRUN                    launcher (default "mpirun -np 1"); reference data are single-rank
#
# Prints a per-case PASS/FAIL summary; exit status 0 only if every case passes.
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
PY="${ROOT}/gpu/.venv/bin/python"
BIN="${ROOT}/bin"
NTRACE="${GS2_REFDATA_NSTEP_TRACE:-10}"
NLADDER="${GS2_REFDATA_DT_LADDER:-2}"
COMPRESSION="${GS2_REFDATA_COMPRESSION:-1}"
OUTROOT="${GS2_REFDATA_DIR:-${ROOT}/gpu/refdata}"
MPIRUN="${MPIRUN:-mpirun -np 1}"
ONLY="${GS2_REFDATA_ONLY:-}"
export OMP_NUM_THREADS=1
unset GK_FFTW3_WISDOM || true   # no shared FFTW wisdom (09 §11.14); every run uses a fresh directory

if [[ "${GS2_REFDATA_BUILD:-0}" == 1 ]]; then
  "${SCRIPT_DIR}/build_reference.sh" || { echo "make_refdata: build failed" >&2; exit 2; }
fi
for b in gs2 gs2_export; do
  [[ -x "${BIN}/${b}" ]] || { echo "make_refdata: ${BIN}/${b} missing; run gpu/scripts/build_reference.sh" >&2; exit 2; }
done
if ! "${BIN}/gs2" --build-config | grep -q -- "-ffp-contract=off"; then
  if [[ "${GS2_REFDATA_ALLOW_STOCK:-0}" != 1 ]]; then
    echo "make_refdata: bin/gs2 is not a reference build (-ffp-contract=off missing); run" \
         "gpu/scripts/build_reference.sh or set GS2_REFDATA_ALLOW_STOCK=1" >&2
    exit 2
  fi
  echo "make_refdata: WARNING: using stock-flag binaries"
fi

# ---------------------------------------------------------------- case list
INPUTS=()
if [[ $# -gt 0 ]]; then
  for a in "$@"; do
    if [[ -f "$a" ]]; then INPUTS+=("$(cd "$(dirname "$a")" && pwd)/$(basename "$a")")
    elif [[ -f "${ROOT}/gpu/cases/$a.in" ]]; then INPUTS+=("${ROOT}/gpu/cases/$a.in")
    else echo "make_refdata: no such case '$a'" >&2; exit 2; fi
  done
else
  shopt -s nullglob
  for f in "${ROOT}"/gpu/cases/*.in; do INPUTS+=("$f"); done
  if [[ ${#INPUTS[@]} -eq 0 ]]; then
    echo "make_refdata: gpu/cases has no *.in; using the GS2 Cyclone tests"
    INPUTS+=("${ROOT}/tests/linear_tests/cyclone_itg/cyclone_itg_low_res_base.in")
    for f in "${ROOT}"/tests/nonlinear_tests/cyclone_itg/*.in; do INPUTS+=("$f"); done
  fi
  shopt -u nullglob
fi
mkdir -p "${OUTROOT}"

# Copy an input file and the files it !includes into a run directory
stage_input() {  # $1 input file, $2 directory
  local inp="$1" dir="$2" f
  cp "${inp}" "${dir}/"
  while read -r f; do
    [[ -n "$f" && -f "$(dirname "${inp}")/$f" ]] && cp "$(dirname "${inp}")/$f" "${dir}/"
  done < <(sed -n 's/^[[:space:]]*!include[[:space:]]\{1,\}\([^[:space:]]*\).*/\1/p' "${inp}")
}

# Remove one run directory created by this script (refuse anything outside OUTROOT)
clean_dir() {
  local d="$1"
  case "$d" in
    "${OUTROOT}"/?*/?*) [[ -d "$d" ]] && rm -rf -- "$d" ;;
    *) echo "make_refdata: refusing to remove '$d'" >&2; exit 3 ;;
  esac
}

SUMMARY=()
NFAIL=0
T0=$(date +%s)
for inp in "${INPUTS[@]}"; do
  case="$(basename "${inp}" .in)"
  dir="${OUTROOT}/${case}"
  echo "================================================================ ${case}"
  mkdir -p "${dir}"
  if [[ "${ONLY}" == export ]]; then subs="export"; else subs="gs2 export trace"; fi
  for sub in ${subs}; do clean_dir "${dir}/${sub}"; mkdir -p "${dir}/${sub}"; stage_input "${inp}" "${dir}/${sub}"; done
  rm -f "${dir}/STATUS" "${dir}/check_export.log" "${dir}/${case}.gs2export.nc"
  [[ "${ONLY}" == export ]] || rm -f "${dir}/${case}.gs2trace.nc" "${dir}/${case}.out.nc"
  ts=$(date +%s)

  # 1. full reference run of bin/gs2
  if [[ "${ONLY}" == export ]]; then
    [[ -f "${dir}/gs2/${case}.out.nc" ]]; rc_gs2=$?
  else
    (cd "${dir}/gs2" && ${MPIRUN} "${BIN}/gs2" "${case}.in" > gs2.log 2>&1)
    rc_gs2=$?
  fi

  # 2. export file with the dt ladder (no stepping)
  cat > "${dir}/export/export_knobs.nml" <<EOF
&gs2_export_knobs
 dt_ladder_n = ${NLADDER}
 compression = ${COMPRESSION}
 write_restart_t0 = .true.
/
EOF
  (cd "${dir}/export" && ${MPIRUN} "${BIN}/gs2_export" "${case}.in" export_knobs.nml > gs2_export.log 2>&1)
  rc_exp=$?

  # 3. traced run: all probes for the first N steps
  [[ "${ONLY}" == export ]] || cat > "${dir}/trace/trace_knobs.nml" <<EOF
&gs2_export_knobs
 nstep_trace = ${NTRACE}
 trace_stages = 'all'
 compression = ${COMPRESSION}
/
EOF
  if [[ "${ONLY}" == export ]]; then
    [[ -f "${dir}/trace/${case}.gs2trace.nc" ]]; rc_tr=$?
  else
    (cd "${dir}/trace" && ${MPIRUN} "${BIN}/gs2_export" "${case}.in" trace_knobs.nml > gs2_export.log 2>&1)
    rc_tr=$?
  fi

  # 4. checks
  exp="${dir}/export/${case}.gs2export.nc"
  trc="${dir}/trace/${case}.gs2trace.nc"
  rc_chk=99
  if [[ ${rc_exp} -eq 0 && -f "${exp}" ]]; then
    args=("${exp}" --outnc "${dir}/export/${case}.out.nc")
    [[ ${rc_tr} -eq 0 && -f "${trc}" ]] && args+=(--trace "${trc}")
    [[ ${rc_tr} -eq 0 && -f "${dir}/trace/${case}.gs2export.nc" ]] && args+=(--compare-export "${dir}/trace/${case}.gs2export.nc")
    [[ ${rc_gs2} -eq 0 && -f "${dir}/gs2/${case}.out.nc" ]] && args+=(--ref-outnc "${dir}/gs2/${case}.out.nc")
    [[ ${rc_gs2} -eq 0 && ${rc_tr} -eq 0 ]] && args+=(--trace-outnc "${dir}/trace/${case}.out.nc")
    rst="$(ls "${dir}"/export/*.gs2export_t0.0 2>/dev/null | head -1)"
    [[ -n "${rst}" ]] && args+=(--restart "${rst}")
    "${PY}" "${SCRIPT_DIR}/check_export.py" "${args[@]}" > "${dir}/check_export.log" 2>&1
    rc_chk=$?
  fi
  [[ "${GS2_REFDATA_KEEP_TRACE_EXPORT:-0}" == 1 ]] || rm -f "${dir}/trace/${case}.gs2export.nc"
  [[ -f "${exp}" ]] && ln -sf "export/${case}.gs2export.nc" "${dir}/${case}.gs2export.nc"
  [[ -f "${trc}" ]] && ln -sf "trace/${case}.gs2trace.nc" "${dir}/${case}.gs2trace.nc"
  [[ -f "${dir}/gs2/${case}.out.nc" ]] && ln -sf "gs2/${case}.out.nc" "${dir}/${case}.out.nc"

  checks="$(grep -E '^check_export: [0-9]+ passed' "${dir}/check_export.log" 2>/dev/null | sed 's/^check_export: //')"
  gates="$(grep -ho 'gate_status = .*' "${dir}/export/gs2_export.log" 2>/dev/null | head -1 | sed 's/gate_status = //')"
  size="$(du -sh "${dir}" 2>/dev/null | cut -f1)"
  if [[ ${rc_gs2} -eq 0 && ${rc_exp} -eq 0 && ${rc_tr} -eq 0 && ${rc_chk} -eq 0 ]]; then
    status=PASS
  else
    status=FAIL
    NFAIL=$((NFAIL + 1))
    grep -E "FAIL" "${dir}/check_export.log" 2>/dev/null | head -10
    [[ ${rc_gs2} -ne 0 ]] && tail -5 "${dir}/gs2/gs2.log"
    [[ ${rc_exp} -ne 0 ]] && tail -5 "${dir}/export/gs2_export.log"
    [[ ${rc_tr} -ne 0 ]] && tail -5 "${dir}/trace/gs2_export.log"
  fi
  line="$(printf '%-34s %-4s gs2=%d export=%d trace=%d check=%d  [%s]  gates=%s  %s  %ds' "${case}" "${status}" \
          ${rc_gs2} ${rc_exp} ${rc_tr} ${rc_chk} "${checks:-no checks}" "${gates:-?}" "${size}" $(( $(date +%s) - ts )))"
  echo "${line}"
  echo "${line}" > "${dir}/STATUS"
  SUMMARY+=("${line}")
done

echo
echo "================================================================ make_refdata summary"
echo "binaries: $(${BIN}/gs2 --build-config | sed -n 's/^Build flags: .* -O3/-O3/p' | cut -c1-60)  N=${NTRACE} ladder=${NLADDER} compression=${COMPRESSION}"
for l in "${SUMMARY[@]}"; do echo "${l}"; done
echo "$(( ${#SUMMARY[@]} - NFAIL )) of ${#SUMMARY[@]} cases PASS; total $(( $(date +%s) - T0 ))s; $(du -sh "${OUTROOT}" | cut -f1) in ${OUTROOT}"
[[ ${NFAIL} -eq 0 ]]
