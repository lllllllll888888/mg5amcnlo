#!/usr/bin/env bash
#
# Regenerate parton-level events from the already-integrated NLO process -- i.e.
# use the "NLO gridpack" that ./generate.sh built. It reuses the stored MINT
# integration grids via
#       generate_events --only_generation -p
# so it SKIPS the (expensive) re-integration and only runs the event-generation
# step. Each call makes a fresh parton-level LHE run (run_02, run_03, ...) that
# the EW Sudakov reweight + shower pipeline can then consume.
#
# This is the NLO analogue of an LO gridpack's run.sh: aMC@NLO has no portable
# gridpack.tar.gz (gridpack=True is LO/MadEvent only), so the reusable object is
# the integrated process directory itself, driven through its own aMCatNLO.
#
# Prerequisite: ./generate.sh has been run once (process built AND integrated;
# the survey/integration grids must already exist under the process dir).
#
# Usage:  ./regenerate_events.sh [NEVENTS] [SEED] [RUNNAME]
#     NEVENTS  number of events to generate     (default: 10000)
#     SEED     random seed (0 = auto)           (default: 0)
#     RUNNAME  explicit run name                (default: MG auto-increments run_NN)
#
# Only NEVENTS/SEED are safe to vary here: the grids are fixed by the build, so
# changing cuts / accuracy / process settings would be inconsistent with them
# (you would need a fresh ./generate.sh for that).

set -eo pipefail

# ── env (mirror generate.sh) ──────────────────────────────────────────────────
MG5_DIR="$(cd "$(dirname "$0")" && pwd)"
HEPTOOLS="${MG5_DIR}/HEPTools"
export LD_LIBRARY_PATH="${HEPTOOLS}/lhapdf6_py3/lib:${HEPTOOLS}/pythia8/lib:${HEPTOOLS}/hepmc/lib:${HEPTOOLS}/zlib/lib:${HEPTOOLS}/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export PYTHONPATH="${HEPTOOLS}/lhapdf6_py3/lib64/python3.9/site-packages${PYTHONPATH:+:$PYTHONPATH}"
export LHAPDF_DATA_PATH="${HEPTOOLS}/lhapdf6_py3/share/LHAPDF${LHAPDF_DATA_PATH:+:$LHAPDF_DATA_PATH}"

PROC="PROCNLO_loop_qcd_qed_sm_Gmu_forSudakov_0"
PROC_DIR="${MG5_DIR}/${PROC}"

NEVENTS=${1:-10000}
SEED=${2:-0}
RUNNAME=${3:-}

if [[ ! -x "${PROC_DIR}/bin/aMCatNLO" ]]; then
  echo "ERROR: ${PROC_DIR}/bin/aMCatNLO not found." >&2
  echo "       Run ./generate.sh first to build + integrate the gridpack." >&2
  exit 1
fi

# Command file for the process-local aMCatNLO. No -f: the card editor reads the
# 'set' lines (nevents/seed) before regenerating from the stored grids.
CMD="$(mktemp "${TMPDIR:-/tmp}/regen_XXXXXX.cmd")"
trap 'rm -f "$CMD"' EXIT
{
  echo "generate_events --only_generation -p${RUNNAME:+ -n ${RUNNAME}}"
  echo "set nevents ${NEVENTS}"
  echo "set iseed ${SEED}"
} > "$CMD"

echo "=========================================="
echo "Regenerating from NLO gridpack (grid reuse)"
echo "  Process : ${PROC}"
echo "  Events  : ${NEVENTS}"
echo "  Seed    : ${SEED}"
echo "  Run name: ${RUNNAME:-<auto>}"
echo "=========================================="

( cd "$PROC_DIR" && ./bin/aMCatNLO < "$CMD" )

echo
echo "Done. New parton-level LHE under ${PROC_DIR}/Events/"
