#!/usr/bin/env bash

# Fixed-Order NLO QCD+EW Generation Script for Z+jets (single mixed run)
#
# Runs a single fixed-order [QCD QED] integration. The custom analysis
# (analysis_FO_QCDvsQCDEW.f) splits the per-PS-point amp_split vector
# into four correlated histogram slices, all filled at the same PS point:
#   |T@NLO     full NLO total       (no filter -- numerator of K_EW)
#   |T@LO      Born only            (ibody=3)
#   |T@NLOQCD  qed_squared=2 NLO    (ibody.ne.3 AND qed_squared==2)
#   |T@NLOEW   qed_squared>=4 NLO   (ibody.ne.3 AND qed_squared>=4)
# Slots are disjoint and exhaustive (i=2+i=3+i=4 == i=1). The K_EW
# ratio is post-processed as (numerator) / (i=2 + i=3), with full
# MC-stat correlation between numerator and denominator -- impossible
# from two independent [QCD] and [QED] runs. See the analysis file for
# the precise physical interpretation of each slot in the absence of
# a Born coupling restriction.
#
# This is the FO leg of the closure test:
#       (FxFx+EWSud)/FxFx  vs  (LO + delta_QCD + delta_EW)/(LO + delta_QCD)   [FO]
#                              = 1 + delta_EW / sigma_NLO_QCD
#
# NO Born restriction (no aS/aEW): MG5 generates all Born coupling
# combinations (alpha_s^N alpha^1, ..., alpha_s^0 alpha^(N+1) for
# Z+Njets) and all NLO QCD + NLO EW corrections on each. The four
# analysis slots therefore have a coarser interpretation than they
# would in the Born-restricted case -- see analysis_FO_QCDvsQCDEW.f
# for details. The slot algebra (i=2+i=3+i=4 = i=1) is preserved.
#
# WHY THE BARE loop_qcd_qed_sm_Gmu MODEL (NOT _forSudakov):
# The _forSudakov model variant carries the Denner-Pozzorini Sudakov
# machinery (used by the Sudakov+PS leg of the closure test) and does
# NOT correctly support the full NLO EW [QCD QED] expansion -- it is
# tuned for the Sudakov approximation only. The bare loop_qcd_qed_sm_Gmu
# is the right model for genuine FO NLO QCD+EW. Widths are zeroed via
# interactive `set WW/WZ/WT/WH 0.0` in the launch CMD below (the bare
# model is not assumed to ship with a restrict_no_widths.dat).
#
# Usage:
#   ./generate_FO.sh PROCESS
#
#   PROCESS : Zj   (p p > z j   [QCD QED])
#             Zjj  (p p > z j j [QCD QED])
#
# Examples:
#   ./generate_FO.sh Zj    # NLO mixed  for Z+1j
#   ./generate_FO.sh Zjj   # NLO mixed  for Z+2j
#
# Environment variables (override defaults):
#   SEED          random seed           (default: 0 = auto)
#   LHAID         LHAPDF ID            (default: 324900 = NNPDF31_nlo_as_0118_luxqed)
#   MODEL         UFO model name       (default: loop_qcd_qed_sm_Gmu)
#   REQ_ACC_FO    target FO accuracy   (default: 0.001)
#   FO_ANALYSIS   analysis .f file     (default: analysis_FO_QCDvsQCDEW.f)
#   PTJCUT        min jet pT in GeV    (default: 30)
#   JETRADIUS     jet radius           (default: 0.4)

# pipefail is needed: `MG5 ... 2>&1 | tee LOG` returns 0 when tee
# succeeds, which would mask MG5/calculate_xsect failures and let the
# script report success on a half-broken run.
set -eo pipefail

# The script uses relative paths (./bin/mg5_aMC, ./Template/...) and
# must be run from the MG5 root directory (mg5amcnlo_devewsudakov_on_valentin).
# Fail loud and early if we are not there.
if [ ! -x "./bin/mg5_aMC" ] || [ ! -d "./Template/NLO/FixedOrderAnalysis" ]; then
  echo "ERROR: this script must be run from the MG5 root directory." >&2
  echo "       expected ./bin/mg5_aMC and ./Template/NLO/ here, got" >&2
  echo "       \$PWD = $(pwd)" >&2
  exit 1
fi

# Ensure LHAPDF can find PDF sets in MG5's local data directory
MG5_LHAPDF_DATA="$(dirname "$0")/HEPTools/lhapdf6_py3/share/LHAPDF"
if [ -d "$MG5_LHAPDF_DATA" ]; then
  export LHAPDF_DATA_PATH="${MG5_LHAPDF_DATA}${LHAPDF_DATA_PATH:+:$LHAPDF_DATA_PATH}"
fi

#-----------------------------------------------------------------------------
# Parse arguments
#-----------------------------------------------------------------------------
PROCESS=${1:?"Usage: $0 PROCESS  (PROCESS=Zj|Zjj)"}

# No Born coupling restriction: MG5 generates the full mixed expansion
# at NLO including all Born types (QCD-Born, EW-Born, mixed Borns at
# Z+2j) and their NLO QCD + NLO EW corrections.
case "$PROCESS" in
  Zj)  GEN_LINE="generate p p > z j   [QCD QED]" ;;
  Zjj) GEN_LINE="generate p p > z j j [QCD QED]" ;;
  *)   echo "ERROR: PROCESS must be 'Zj' or 'Zjj', got '$PROCESS'"; exit 1 ;;
esac

#-----------------------------------------------------------------------------
# Defaults (overridable via environment)
#-----------------------------------------------------------------------------
SEED=${SEED:-0}
LHAID=${LHAID:-324900}
MODEL=${MODEL:-loop_qcd_qed_sm_Gmu}
REQ_ACC_FO=${REQ_ACC_FO:-0.001}
FO_ANALYSIS=${FO_ANALYSIS:-analysis_FO_QCDvsQCDEW.f}
PTJCUT=${PTJCUT:-30}
JETRADIUS=${JETRADIUS:-0.4}

MG_BIN="./bin/mg5_aMC"
PROC="PROC_FO_${PROCESS}_NLOmixed"
PROC_DIR="./${PROC}"
TAG="${PROCESS}_NLOmixed"

echo "=========================================="
echo "Fixed-Order NLO QCD+EW (mixed) Generation"
echo "=========================================="
echo "Process : $PROCESS  ($GEN_LINE)"
echo "Model   : $MODEL"
echo "Output  : $PROC"
echo "LHAPDF  : $LHAID"
echo "Seed    : $SEED"
echo "Accuracy: $REQ_ACC_FO"
echo "Analysis: $FO_ANALYSIS"
echo "Jet def : anti-kT, R=$JETRADIUS, pT>$PTJCUT GeV"
echo "=========================================="

#-----------------------------------------------------------------------------
# 1) Process generation
#    - define p = p a : include photon in proton (needed for NLO EW real
#      emission subtraction and photon-initiated contributions)
#    - define j = p   : photon included in jet definition
#-----------------------------------------------------------------------------
CMD1="run_fo_gen_${TAG}.cmd"
cat <<EOF > "$CMD1"
import model ${MODEL}
define p = p a
define j = p
${GEN_LINE}
output ${PROC}
EOF

echo ""
echo "Step 1: Creating process directory..."
"$MG_BIN" "$CMD1" 2>&1 | tee "fo_gen_${TAG}.log"

#-----------------------------------------------------------------------------
# 2) Configure the fixed-order analysis
#    Copy chosen analysis into the process FixedOrderAnalysis/ directory
#    and update the FO_analyse_card to point to it.
#-----------------------------------------------------------------------------
ANALYSIS_SRC="Template/NLO/FixedOrderAnalysis/${FO_ANALYSIS}"
if [ ! -f "$ANALYSIS_SRC" ]; then
  echo "WARNING: Analysis file '$ANALYSIS_SRC' not found."
  echo "         Falling back to analysis_HwU_general.f"
  FO_ANALYSIS="analysis_HwU_general.f"
  ANALYSIS_SRC="Template/NLO/FixedOrderAnalysis/${FO_ANALYSIS}"
fi

echo ""
echo "Step 2: Configuring analysis (${FO_ANALYSIS})..."
cp "$ANALYSIS_SRC" "${PROC_DIR}/FixedOrderAnalysis/${FO_ANALYSIS}"

# Update FO_analyse_card: replace the default analysis .o reference
ANALYSIS_OBJ="${FO_ANALYSIS%.f}.o"
sed -i "s|^FO_ANALYSE.*|FO_ANALYSE = ${ANALYSIS_OBJ}|" "${PROC_DIR}/Cards/FO_analyse_card.dat"

#-----------------------------------------------------------------------------
# 3) Resolve custom dummy_fct.f path (if present)
#-----------------------------------------------------------------------------
DUMMY_FCT_PATH=""
if [ -f dummy_fct.f ]; then
  DUMMY_FCT_PATH="$(readlink -f dummy_fct.f)"
  echo "         Will use custom_fcts: ${DUMMY_FCT_PATH}"
fi

#-----------------------------------------------------------------------------
# 4) Launch fixed-order run
#    Key differences from NLO+PS (generate.sh):
#    - --fixed_order instead of --parton
#    - No shower / merging settings (ickkw, parton_shower)
#    - Integration accuracy controlled by req_acc_FO
#    - Histograms filled on-the-fly during integration
#-----------------------------------------------------------------------------
CMD2="run_fo_launch_${TAG}.cmd"
cat <<EOF > "$CMD2"
launch ${PROC} --fixed_order

$([ -n "$DUMMY_FCT_PATH" ] && echo "set custom_fcts ${DUMMY_FCT_PATH}")

# Integration accuracy
set req_acc_FO ${REQ_ACC_FO}
set iseed ${SEED}

# Scale choice: HT/2 (same as FxFx runs)
set dynamical_scale_choice 3
set fixed_ren_scale False
set fixed_fac_scale False

# 5-flavour scheme: massless b-quark
set maxjetflavor 5
set mass 5 0.0

# Widths: zero W, Z, t, H to match the stable-boson Sudakov-limit
# prediction on the other leg. The bare loop_qcd_qed_sm_Gmu does not
# ship with a restrict_no_widths.dat so we have to do this interactively.
set WW 0.0
set WZ 0.0
set WT 0.0
set WH 0.0

# PDF: NNPDF31_nlo_as_0118_luxqed (includes photon PDF)
set pdlabel lhapdf
set lhaid ${LHAID}

# Jet definition: anti-kT algorithm
set jetalgo -1
set jetradius ${JETRADIUS}
set ptj ${PTJCUT}
set etaj -1

# Cluster photons into jets (as specified in plan.txt)
set gamma_is_j True

# Photon recombination with fermions (relevant for NLO EW real emission)
set lepphreco True
set quarkphreco True
EOF

echo ""
echo "Step 3: Running fixed-order integration..."
"$MG_BIN" "$CMD2" 2>&1 | tee "fo_launch_${TAG}.log"

#-----------------------------------------------------------------------------
# Check output
#-----------------------------------------------------------------------------
HWU_DIR="${PROC_DIR}/Events/run_01"
if ls "${HWU_DIR}"/*.HwU 1>/dev/null 2>&1; then
  echo ""
  echo "Fixed-order run complete."
  echo "  Histograms: ${HWU_DIR}/"
  ls "${HWU_DIR}"/*.HwU "${HWU_DIR}"/*.gnuplot 2>/dev/null | head -10
else
  echo ""
  echo "WARNING: No .HwU files found in ${HWU_DIR}/"
  echo "         Check fo_launch_${TAG}.log for errors."
fi