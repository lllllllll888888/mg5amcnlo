#!/usr/bin/env bash
# =============================================================================
# Parallel LHE Pipeline for Z+jets EW Sudakov Validation (v1)
# =============================================================================
# Processes LHE events through two parallel paths:
#   Path a: EW Sudakov reweight → shower (QCD+QED)  [NLO_QCD+PS+EWS]
#   Path b: shower only (QCD)                        [NLO_QCD+PS]
#
# No MadSpin — Z is kept stable throughout (production-level only).
# Designed for pp → Z + 0,1,2j @ NLO QCD with FxFx merging.
#
# Usage: ./parallel_pipeline_DYj.sh <input.lhe[.gz]> [options via env vars]
#
# Environment variables:
#   EVENTS_PER_CHUNK=50000    # Events per chunk (default: 50k)
#   MODE=both                 # both, path_a, path_b
#   MAX_CONCURRENT=24         # Max concurrent processes
#   DRY_RUN=false             # If true, show what would run
#   DEBUG_PATH_A_REWEIGHT=false  # If true, run a single input LHE through path_a reweight only
#   DEBUG_MAX_EVENTS=1000     # In debug path_a mode, build/use a reduced LHE with first N events (0 = all)
#   PROC_DIR=./pc16/madevent  # Path to the MG5 process directory
#   VERBOSE=false             # Enable debug logging
#   QCUT=30.0                 # FxFx merging scale (default: 30)
#   NJMAX=2                   # Max jet multiplicity for FxFx (default: 2)
#   CLEANUP_INTERMEDIATE=true # Clean temp files after each step
#   CLEANUP_PROCNLO=true      # Remove procnlo after worker completes
#   KEEP_LOGS=false           # Keep all logs during cleanup
#
# Examples:
#   ./parallel_pipeline_DYj.sh pc16/madevent/Events/run_pc16_seed1003/events.lhe.gz
#   QCUT=20 MODE=path_b ./parallel_pipeline_DYj.sh events.lhe.gz
#   EVENTS_PER_CHUNK=100000 QCUT=50 ./parallel_pipeline_DYj.sh events.lhe.gz
#   DEBUG_PATH_A_REWEIGHT=true ./parallel_pipeline_DYj.sh chunk_000.lhe
#   DEBUG_PATH_A_REWEIGHT=true DEBUG_MAX_EVENTS=1000 ./parallel_pipeline_DYj.sh events.lhe.gz
# =============================================================================

set -euo pipefail

# =============================================================================
# Configuration
# =============================================================================
INPUT_LHE="${1:?Usage: $0 <input.lhe[.gz]>}"
EVENTS_PER_CHUNK="${EVENTS_PER_CHUNK:-50000}"
MODE="${MODE:-both}"
NCPU=$(nproc 2>/dev/null || echo 32)
MAX_CONCURRENT="${MAX_CONCURRENT:-$NCPU}"
DRY_RUN="${DRY_RUN:-false}"
DEBUG_PATH_A_REWEIGHT="${DEBUG_PATH_A_REWEIGHT:-false}"
DEBUG_MAX_EVENTS="${DEBUG_MAX_EVENTS:-1000}"
VERBOSE="${VERBOSE:-false}"

# FxFx parameters
QCUT="${QCUT:-30.0}"
NJMAX="${NJMAX:-2}"

# EW Sudakov clustering threshold scan: space-separated list of ξ values for s_ij > ξ·M_W².
# Empty → legacy single ξ=1.5 (one weight pair 2001/2101 per event).
# Example: XI_SCAN="0.5 1.0 1.5 2.0 4.0" → five weight pairs (2001/2101 ... 2801/2901).
XI_SCAN="${XI_SCAN:-}"

# FxFx merging-scale variation list: space-separated qcut values (GeV) injected
# into Pythia as JetMatching:qCutList. Consumed by the locally patched
# Pythia83_hep driver to write per-event accept flags as named HepMC weights
# "FxFx_qCutAccept_<value>". Cross-section at qcut X is recovered as
#     sigma(X) = sum_evt  w_nominal[evt] * FxFx_qCutAccept_<X>[evt]
# Empty → variation pipeline disabled (single-qcut FxFx run, no extra weights).
QCUT_LIST="${QCUT_LIST:-20 30 40}"
# Build the dict entry. The key MUST be quoted because banner.py's dict parser
# does rsplit(':',1) -- without quotes it would split "JetMatching:qCutList" on
# its internal colon. Quotes get stripped at banner.py:1238-1239.
_qcut_list_entry=""
[[ -n "$QCUT_LIST" ]] && _qcut_list_entry=", \"JetMatching:qCutList\":${QCUT_LIST}"

# Disk optimization
CLEANUP_INTERMEDIATE="${CLEANUP_INTERMEDIATE:-true}"
CLEANUP_PROCNLO="${CLEANUP_PROCNLO:-true}"
KEEP_LOGS="${KEEP_LOGS:-false}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROC_DIR="${PROC_DIR:-$SCRIPT_DIR/PROCNLO_loop_qcd_qed_sm_Gmu_forSudakov_0}"

# ── HEPTools paths ──────────────────────────────────────────────────────────
HEPTOOLS="${SCRIPT_DIR}/HEPTools"
export LD_LIBRARY_PATH="${HEPTOOLS}/lhapdf6_py3/lib:${HEPTOOLS}/pythia8/lib:${HEPTOOLS}/hepmc/lib:${HEPTOOLS}/zlib/lib:${HEPTOOLS}/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export PYTHONPATH="${HEPTOOLS}/lhapdf6_py3/lib64/python3.9/site-packages${PYTHONPATH:+:$PYTHONPATH}"
export LHAPDF_DATA_PATH="${HEPTOOLS}/lhapdf6_py3/share/LHAPDF${LHAPDF_DATA_PATH:+:$LHAPDF_DATA_PATH}"

# Validate
[[ -f "$INPUT_LHE" ]] || { echo "ERROR: Input file not found: $INPUT_LHE"; exit 1; }
[[ -d "$PROC_DIR" ]] || { echo "ERROR: PROC_DIR not found: $PROC_DIR"; exit 1; }
[[ -f "$PROC_DIR/bin/aMCatNLO" ]] || { echo "ERROR: No bin/aMCatNLO in PROC_DIR: $PROC_DIR"; exit 1; }

# =============================================================================
# Embedded Card Templates
# =============================================================================

# Reweight card: enables EW Sudakov corrections
_xi_scan_line=""
[[ -n "$XI_SCAN" ]] && _xi_scan_line="change sudakov_xi $XI_SCAN"
REWEIGHT_CARD_TEMPLATE="#*************************************************************************
#                          Reweight Module                               *
#              EW Sudakov reweighting for Z+jets                         *
#*************************************************************************

change mode NLO

# Enable EW Sudakov reweighting
change include_sudakov True
${_xi_scan_line}

launch"

# Shower card generator — two variants: with and without QED shower
generate_shower_card() {
    local qed_shower="${1:-F}"  # T for path_a, F for path_b
    cat << SHOWERCARD_EOF
#***********************************************************************
#                        MadGraph5_aMC@NLO                             *
#                      shower_card.dat aMC@NLO                         *
#  Generated by parallel_pipeline_DYj.sh                               *
#***********************************************************************
nevents      = -1
nsplit_jobs  = 1
combine_td   = T
maxprint     = 2
maxerrs      = 0.1
rnd_seed     = 0
rnd_seed2    = 0
pdfcode      = 1
ue_enabled   = F
hadronize    = F
lambda_5     = -1.0
b_stable     = F
pi_stable    = T
wp_stable    = F
wm_stable    = F
z_stable     = T
h_stable     = F
tap_stable   = F
tam_stable   = F
mup_stable   = F
mum_stable   = F
b_mass       = -1.0
is_4lep      = F
is_bbar      = F
Qcut         = ${QCUT}
njmax        = ${NJMAX}
qed_shower   = ${qed_shower}
primordialkt = F
pythia8_options = {23:mayDecay = off${_qcut_list_entry}}
space_shower_me_corrections = F
time_shower_me_corrections  = T
time_shower_me_extended     = F
time_shower_me_after_first  = F
EXTRALIBS    = stdhep Fmcfio
EXTRAPATHS   = ../lib
INCLUDEPATHS =
ANALYSE      =
SHOWERCARD_EOF
}

# Write cards to a directory
# $1 = directory, $2 = "a" or "b" (determines QED shower)
write_cards_to_dir() {
    local cards_dir="$1"
    local path_id="${2:-a}"

    # Write reweight card (only meaningful for path_a, but harmless to have)
    echo "$REWEIGHT_CARD_TEMPLATE" > "$cards_dir/reweight_card.dat"

    # Write shower card
    if [[ "$path_id" == "a" ]]; then
        generate_shower_card "T" > "$cards_dir/shower_card.dat"
    else
        generate_shower_card "F" > "$cards_dir/shower_card.dat"
    fi

    log_debug "write_cards: Wrote cards to $cards_dir (path=$path_id, qed_shower=$([ "$path_id" = a ] && echo T || echo F))"
}

# =============================================================================
# Logging Functions
# =============================================================================
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }
log_error() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: $*" >&2; }
log_debug() { [[ "$VERBOSE" == "true" ]] && echo "[$(date '+%Y-%m-%d %H:%M:%S')] DEBUG: $*"; return 0; }
log_step() {
    local worker_id="$1" step="$2" status="$3"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [$worker_id] $step: $status"
}

get_elapsed() { local now; now=$(date +%s); echo $(($now - $1)); }

format_duration() {
    local s="$1"
    local h=$((s / 3600)) m=$(((s % 3600) / 60)) sec=$((s % 60))
    if [[ $h -gt 0 ]]; then printf "%dh %dm %ds" "$h" "$m" "$sec"
    elif [[ $m -gt 0 ]]; then printf "%dm %ds" "$m" "$sec"
    else printf "%ds" "$sec"; fi
}

verify_file() {
    local file="$1" context="$2"
    if [[ ! -f "$file" ]]; then log_error "$context: File not found: $file"; return 1; fi
    local size; size=$(stat -c%s "$file" 2>/dev/null || echo 0)
    if [[ "$size" -eq 0 ]]; then log_error "$context: File is empty: $file"; return 1; fi
    log_debug "$context: Verified $file (${size}B)"
    return 0
}

count_events() { grep -cE '<event[ >]' "$1" 2>/dev/null || echo 0; }

calculate_chunks() { echo $(( ($1 + EVENTS_PER_CHUNK - 1) / EVENTS_PER_CHUNK )); }

extract_first_events() {
    local input="$1" max_events="$2" output="$3"
    local outdir; outdir="$(dirname "$output")"
    mkdir -p "$outdir"

    if [[ "$max_events" -le 0 ]]; then
        cp "$input" "$output"
        return 0
    fi

    # Header = everything before the first event block.
    sed -n '1,/<event[ >]/{ /<event[ >]/!p }' "$input" > "$outdir/header.xml"

    awk -v max_events="$max_events" '
    BEGIN { in_event=0; event_count=0 }
    /<event[ >]/ {
        if (event_count >= max_events) exit
        in_event=1
        event_count++
    }
    in_event { print }
    /<\/event>/ { in_event=0 }
    ' "$input" > "$outdir/events.txt"

    cat "$outdir/header.xml" "$outdir/events.txt" > "$output"
    echo "</LesHouchesEvents>" >> "$output"
    rm -f "$outdir/header.xml" "$outdir/events.txt"
}

# =============================================================================
# Split LHE into Chunks
# =============================================================================
split_lhe() {
    local master="$1" nchunks="$2" outdir="$3"
    local start_time; start_time=$(date +%s)

    log "Splitting $master into $nchunks chunks..."
    mkdir -p "$outdir"

    # Extract header (everything before first <event>)
    sed -n '1,/<event[ >]/{ /<event[ >]/!p }' "$master" > "$outdir/header.xml"

    local nevents; nevents=$(count_events "$master")
    local epc=$(( (nevents + nchunks - 1) / nchunks ))
    log "  Total events: $nevents, ~$epc per chunk"

    awk -v outdir="$outdir" -v epc="$epc" -v nchunks="$nchunks" '
    BEGIN { chunk=0; count=0; in_event=0 }
    /<event[ >]/ {
        in_event=1
        if (count > 0 && count % epc == 0 && chunk < nchunks-1) chunk++
        count++
    }
    in_event { print >> (outdir "/events_" chunk ".txt") }
    /<\/event>/ { in_event=0 }
    ' "$master"

    for i in $(seq 0 $((nchunks-1))); do
        local events_file="$outdir/events_$i.txt"
        local chunk_file="$outdir/chunk_$(printf %03d $i).lhe"
        if [[ -f "$events_file" ]]; then
            cat "$outdir/header.xml" "$events_file" > "$chunk_file"
            echo "</LesHouchesEvents>" >> "$chunk_file"
            rm "$events_file"
            local nc; nc=$(count_events "$chunk_file")
            log "  Chunk $i: $nc events"
        fi
    done

    rm "$outdir/header.xml"
    log "  Split complete ($(format_duration "$(get_elapsed "$start_time")"))"
}

# =============================================================================
# Create Shadow PROCNLO
# =============================================================================
create_shadow_procnlo() {
    local shadow="$1" path_id="$2"
    local source="$PROC_DIR"

    log "Creating shadow PROCNLO for path_$path_id: $shadow ..."

    cp -r "$source" "$shadow"

    # Remove stale artifacts
    rm -rf "$shadow"/rw_me "$shadow"/rw_me_density_prod "$shadow"/rw_me_density_decay

    # Write cards for this path
    write_cards_to_dir "$shadow/Cards" "$path_id"
    log "  Wrote cards: shower (Qcut=$QCUT, njmax=$NJMAX, qed=$([ "$path_id" = a ] && echo T || echo F))"
    [[ "$path_id" == "a" ]] && log "  Wrote reweight card (include_sudakov=True)"

    # Fresh Events directory
    rm -rf "$shadow/Events"
    mkdir -p "$shadow/Events/run_01"

    log "  Shadow PROCNLO created"
}

# =============================================================================
# Create Per-Worker PROCNLO
# =============================================================================
create_worker_procnlo() {
    local worker_dir="$1" shadow="$2" worker_num="$3"
    local procnlo="$worker_dir/procnlo"

    mkdir -p "$worker_dir/output"
    cp -r "$shadow" "$procnlo"

    rm -rf "$procnlo/Events"
    mkdir -p "$procnlo/Events/run_01"
}

# =============================================================================
# Cleanup Functions
# =============================================================================
cleanup_step_intermediate() {
    local step="$1" worker_dir="$2" procnlo="$3"
    [[ "$CLEANUP_INTERMEDIATE" != "true" ]] && return 0

    case "$step" in
        reweight)
            rm -rf "$procnlo/SubProcesses/"*/GF*.0 2>/dev/null || true
            [[ "$KEEP_LOGS" != "true" ]] && rm -f "$procnlo/Events/run_01/reweight.log" 2>/dev/null || true
            ;;
        shower)
            rm -f "$procnlo/Events/run_01/events.lhe" 2>/dev/null || true
            rm -rf "$procnlo/MCatNLO/RUN_PYTHIA8_"*/Pythia8.hep 2>/dev/null || true
            ;;
    esac
    return 0
}

cleanup_worker_complete() {
    local worker_dir="$1"
    local procnlo="$worker_dir/procnlo"
    [[ "$CLEANUP_PROCNLO" != "true" ]] && return 0
    [[ ! -f "$worker_dir/final_output.txt" ]] && return 0

    local final_output; final_output=$(cat "$worker_dir/final_output.txt")
    [[ ! -f "$final_output" ]] && return 0

    rm -rf "$procnlo"
    if [[ "$KEEP_LOGS" != "true" ]]; then
        local output_file
        for output_file in "$worker_dir/output/"*.lhe; do
            [[ -e "$output_file" ]] || continue
            [[ "$output_file" == "$final_output" ]] && continue
            rm -f "$output_file"
        done
    fi
    return 0
}

# =============================================================================
# Processing Functions
# =============================================================================

# EW Sudakov reweight (path_a only)
process_reweight() {
    local input="$1" output="$2" procnlo="$3" logdir="$4"
    local start_time; start_time=$(date +%s)

    verify_file "$input" "reweight input" || return 1

    if [[ "$input" == *.gz ]]; then
        gunzip -c "$input" > "$procnlo/Events/run_01/events.lhe"
    else
        cp "$input" "$procnlo/Events/run_01/events.lhe"
    fi

    local mg_file; mg_file="$(realpath "$logdir")/reweight.mg"
    cat > "$mg_file" << EOF
reweight run_01
0
quit
EOF

    log_debug "reweight: Running aMCatNLO reweight"
    (cd "$procnlo" && ./bin/aMCatNLO < "$mg_file") &> "$logdir/reweight.log"
    local exit_code=$?

    if [[ $exit_code -ne 0 ]]; then
        log_error "reweight: aMCatNLO failed (exit $exit_code). See $logdir/reweight.log"
        return 1
    fi

    # MG5 can exit 0 even on failure — check log for actual errors
    if grep -q "interrupted with error" "$logdir/reweight.log"; then
        log_error "reweight: aMCatNLO reported error (exit 0 but failed). See $logdir/reweight.log"
        return 1
    fi

    if [[ ! -f "$procnlo/Events/run_01/events.lhe" ]]; then
        log_error "reweight: Output events.lhe not found after reweight"
        return 1
    fi

    cp "$procnlo/Events/run_01/events.lhe" "$output"
    verify_file "$output" "reweight output" || return 1

    log_debug "reweight: Complete ($(format_duration "$(get_elapsed "$start_time")"))"
    return 0
}

# Shower (both paths, but different shower_card)
process_shower() {
    local input="$1" output="$2" procnlo="$3" logdir="$4"
    local start_time; start_time=$(date +%s)

    verify_file "$input" "shower input" || return 1

    if [[ "$input" == *.gz ]]; then
        gunzip -c "$input" > "$procnlo/Events/run_01/events.lhe"
    else
        cp "$input" "$procnlo/Events/run_01/events.lhe"
    fi

    # Clean prior shower artifacts
    rm -rf "$procnlo/MCatNLO/RUN_PYTHIA8_"* 2>/dev/null || true
    rm -f "$procnlo/Events/run_01/"events_PYTHIA8_*.hepmc.gz 2>/dev/null || true
    # Clean leftover split files from reweight step
    rm -f "$procnlo/Events/run_01/"events.lhe_* 2>/dev/null || true

    local mg_file; mg_file="$(realpath "$logdir")/shower.mg"
    cat > "$mg_file" << EOF
shower run_01
0
quit
EOF

    log_debug "shower: Running aMCatNLO shower"
    (cd "$procnlo" && ./bin/aMCatNLO < "$mg_file") &> "$logdir/shower.log"
    local exit_code=$?

    if [[ $exit_code -ne 0 ]]; then
        log_error "shower: aMCatNLO failed (exit $exit_code). See $logdir/shower.log"
        return 1
    fi

    # Find output
    local hepmc=""
    hepmc=$(find "$procnlo/MCatNLO" -name "*.hepmc.gz" -type f 2>/dev/null | head -1)
    [[ -z "$hepmc" ]] && hepmc=$(find "$procnlo/Events/run_01" -name "events_PYTHIA8_*.hepmc.gz" -type f 2>/dev/null | head -1)

    if [[ -n "$hepmc" && -f "$hepmc" ]]; then
        cp "$hepmc" "$output"
        verify_file "$output" "shower output" || return 1
        log_debug "shower: Complete ($(format_duration "$(get_elapsed "$start_time")"))"
        return 0
    else
        log_error "shower: No HEPMC output found"
        tail -50 "$logdir/shower.log" 2>&1 | grep -iE "(error|fail|exception)" | head -10 | while read -r line; do log_error "  $line"; done
        return 1
    fi
}

# =============================================================================
# Worker Functions
# =============================================================================

run_step() {
    local step_name="$1" worker_dir="$2"
    shift 2
    echo "$step_name" > "$worker_dir/status"
    "$@"
}

# Path A: reweight (EW Sudakov) → shower (QCD+QED)
run_path_a_worker() {
    local id="$1" chunk="$2" worker_dir="$3"
    local procnlo="$worker_dir/procnlo"
    local logfile="$worker_dir/worker.log"
    local worker_start; worker_start=$(date +%s)

    {
        log_step "A-$id" "START" "reweight(EWS) → shower(QCD+QED)"
        log_step "A-$id" "INPUT" "$(basename "$chunk") ($(count_events "$chunk") events)"

        # Step 1: EW Sudakov reweight
        local step_start; step_start=$(date +%s)
        if ! run_step "reweight" "$worker_dir" \
            process_reweight "$chunk" "$worker_dir/output/reweighted.lhe" "$procnlo" "$worker_dir"; then
            log_step "A-$id" "FAILED" "reweight"
            echo "FAILED:reweight" > "$worker_dir/status"
            return 1
        fi
        log_step "A-$id" "DONE" "reweight ($(format_duration "$(get_elapsed "$step_start")"))"
        cleanup_step_intermediate "reweight" "$worker_dir" "$procnlo"

        # Step 2: Shower (QCD+QED)
        step_start=$(date +%s)
        if ! run_step "shower" "$worker_dir" \
            process_shower "$worker_dir/output/reweighted.lhe" "$worker_dir/output/reweighted.hepmc.gz" "$procnlo" "$worker_dir"; then
            log_step "A-$id" "FAILED" "shower"
            echo "FAILED:shower" > "$worker_dir/status"
            return 1
        fi
        log_step "A-$id" "DONE" "shower ($(format_duration "$(get_elapsed "$step_start")"))"
        cleanup_step_intermediate "shower" "$worker_dir" "$procnlo"

        # Success
        local total_elapsed; total_elapsed=$(get_elapsed "$worker_start")
        echo "$worker_dir/output/reweighted.hepmc.gz" > "$worker_dir/final_output.txt"
        echo "SUCCESS:$(format_duration "$total_elapsed")" > "$worker_dir/status"
        log_step "A-$id" "COMPLETE" "Total: $(format_duration "$total_elapsed")"
        cleanup_worker_complete "$worker_dir"
    } 2>&1 | tee "$logfile"
}

# Path B: shower only (QCD, no EW Sudakov)
run_path_b_worker() {
    local id="$1" chunk="$2" worker_dir="$3"
    local procnlo="$worker_dir/procnlo"
    local logfile="$worker_dir/worker.log"
    local worker_start; worker_start=$(date +%s)

    {
        log_step "B-$id" "START" "shower(QCD only)"
        log_step "B-$id" "INPUT" "$(basename "$chunk") ($(count_events "$chunk") events)"

        # Single step: Shower (QCD only, no reweight)
        local step_start; step_start=$(date +%s)
        if ! run_step "shower" "$worker_dir" \
            process_shower "$chunk" "$worker_dir/output/showered.hepmc.gz" "$procnlo" "$worker_dir"; then
            log_step "B-$id" "FAILED" "shower"
            echo "FAILED:shower" > "$worker_dir/status"
            return 1
        fi
        log_step "B-$id" "DONE" "shower ($(format_duration "$(get_elapsed "$step_start")"))"
        cleanup_step_intermediate "shower" "$worker_dir" "$procnlo"

        # Success
        local total_elapsed; total_elapsed=$(get_elapsed "$worker_start")
        echo "$worker_dir/output/showered.hepmc.gz" > "$worker_dir/final_output.txt"
        echo "SUCCESS:$(format_duration "$total_elapsed")" > "$worker_dir/status"
        log_step "B-$id" "COMPLETE" "Total: $(format_duration "$total_elapsed")"
        cleanup_worker_complete "$worker_dir"
    } 2>&1 | tee "$logfile"
}

# Path A debug: reweight only on a single input LHE
run_path_a_reweight_debug_worker() {
    local chunk="$1" worker_dir="$2"
    local procnlo="$worker_dir/procnlo"
    local logfile="$worker_dir/worker.log"
    local worker_start; worker_start=$(date +%s)

    {
        log_step "A-DBG" "START" "single-file reweight(EWS only)"
        log_step "A-DBG" "INPUT" "$(basename "$chunk") ($(count_events "$chunk") events)"

        local step_start; step_start=$(date +%s)
        if ! run_step "reweight" "$worker_dir" \
            process_reweight "$chunk" "$worker_dir/output/reweighted_debug.lhe" "$procnlo" "$worker_dir"; then
            log_step "A-DBG" "FAILED" "reweight"
            echo "FAILED:reweight" > "$worker_dir/status"
            return 1
        fi

        local total_elapsed; total_elapsed=$(get_elapsed "$worker_start")
        echo "$worker_dir/output/reweighted_debug.lhe" > "$worker_dir/final_output.txt"
        echo "SUCCESS:$(format_duration "$total_elapsed")" > "$worker_dir/status"
        log_step "A-DBG" "DONE" "reweight ($(format_duration "$(get_elapsed "$step_start")"))"
        log_step "A-DBG" "COMPLETE" "Total: $(format_duration "$total_elapsed")"
    } 2>&1 | tee "$logfile"
}

# Export for subshells
export -f log log_error log_debug log_step get_elapsed format_duration verify_file count_events run_step
export -f extract_first_events
export -f process_reweight process_shower
export -f run_path_a_worker run_path_b_worker run_path_a_reweight_debug_worker
export -f generate_shower_card write_cards_to_dir
export -f cleanup_step_intermediate cleanup_worker_complete
export SCRIPT_DIR PROC_DIR VERBOSE QCUT NJMAX QCUT_LIST _qcut_list_entry REWEIGHT_CARD_TEMPLATE
export CLEANUP_INTERMEDIATE CLEANUP_PROCNLO KEEP_LOGS DEBUG_MAX_EVENTS

# HEPTools paths already exported at script top

# =============================================================================
# Merge HEPMC Outputs
# =============================================================================
merge_hepmc() {
    local work_dir="$1" output="$2"
    local start_time; start_time=$(date +%s)

    log "Merging HEPMC files..."
    mkdir -p "$(dirname "$output")"

    local hepmc_files=()
    local failed_count=0

    for marker in "$work_dir"/workers/*/final_output.txt; do
        [[ -f "$marker" ]] || continue
        local hepmc; hepmc=$(cat "$marker")
        if [[ -f "$hepmc" ]]; then
            hepmc_files+=("$hepmc")
        else
            ((failed_count++))
            log_error "  Missing: $hepmc"
        fi
    done

    if [[ ${#hepmc_files[@]} -eq 0 ]]; then
        log_error "  No HEPMC files to merge"
        return 1
    fi

    local n_files=${#hepmc_files[@]}
    local merge_jobs=$MAX_CONCURRENT
    (( merge_jobs > n_files )) && merge_jobs=$n_files
    local per_group=$(( (n_files + merge_jobs - 1) / merge_jobs ))

    log "  Merging $n_files files ($merge_jobs shards)..."

    local shard_dir="$work_dir/merge_shards"
    mkdir -p "$shard_dir"

    local shard_pids=()
    local shard_files=()
    for (( g=0; g<merge_jobs; g++ )); do
        local start_idx=$((g * per_group))
        (( start_idx >= n_files )) && break

        local end_idx=$((start_idx + per_group))
        (( end_idx > n_files )) && end_idx=$n_files

        local group_files=("${hepmc_files[@]:$start_idx:$((end_idx - start_idx))}")
        local shard_file="$shard_dir/shard_$(printf '%03d' $g).hepmc.gz"
        shard_files+=("$shard_file")

        (
            local count=0 first_in_shard=true
            for hepmc in "${group_files[@]}"; do
                if [[ "$first_in_shard" == "true" && "$g" == "0" ]]; then
                    zcat "$hepmc" | sed '/^HepMC::IO_GenEvent-END_EVENT_LISTING$/d' | gzip >> "$shard_file"
                    first_in_shard=false
                else
                    zcat "$hepmc" | awk '/^E / { past=1 } past && !/^HepMC::IO_GenEvent-END/' | gzip >> "$shard_file"
                fi
                rm -f "$hepmc"
                count=$((count + 1))
            done
            log "    [shard $g] Merged $count files"
        ) &
        shard_pids+=($!)
    done

    local shard_failed=0
    for pid in "${shard_pids[@]}"; do
        wait "$pid" || shard_failed=$((shard_failed + 1))
    done

    if [[ $shard_failed -gt 0 ]]; then
        log_error "  $shard_failed shard(s) failed"
        return 1
    fi

    for shard in "${shard_files[@]}"; do
        cat "$shard" >> "$output"
        rm -f "$shard"
    done
    echo "HepMC::IO_GenEvent-END_EVENT_LISTING" | gzip >> "$output"
    rm -rf "$shard_dir"

    local size; size=$(stat -c%s "$output" 2>/dev/null || echo 0)
    log "  Merged $n_files files → $output ($(numfmt --to=iec-i --suffix=B "$size" 2>/dev/null || echo "${size}B")) in $(format_duration "$(get_elapsed "$start_time")")"
}

# =============================================================================
# Progress Monitoring
# =============================================================================
monitor_progress() {
    local work_dir="$1" total_workers="$2"
    local last_status=""

    while true; do
        local running=0 success=0 failed=0

        for status_file in "$work_dir"/*/workers/*/status; do
            [[ -f "$status_file" ]] || continue
            local status; status=$(cat "$status_file" 2>/dev/null)
            case "$status" in
                SUCCESS*) ((success++)) ;;
                FAILED*) ((failed++)) ;;
                *) ((running++)) ;;
            esac
        done

        local new_status="$success/$total_workers done, $running running, $failed failed"
        if [[ "$new_status" != "$last_status" ]]; then
            printf "\r[$(date '+%H:%M:%S')] Progress: %s   " "$new_status"
            last_status="$new_status"
        fi

        [[ $((success + failed)) -ge $total_workers ]] && break
        sleep 5
    done
    echo ""
}

# =============================================================================
# Summary Report
# =============================================================================
generate_summary() {
    local work_dir="$1"
    local report_file="$work_dir/summary_report.txt"

    {
        echo "==========================================================="
        echo " EW Sudakov Z+jets Pipeline — Summary Report"
        echo " Generated: $(date '+%Y-%m-%d %H:%M:%S')"
        echo "==========================================================="
        echo ""
        echo " Configuration:"
        echo "   Process:     Z+jets (stable Z, FxFx merged)"
        echo "   FxFx Qcut:     $QCUT"
        echo "   FxFx njmax:    $NJMAX"
        echo "   FxFx qCutList: ${QCUT_LIST:-<disabled>}"
        echo "   Input:       $INPUT_LHE"
        echo ""

        for path in path_a path_b; do
            [[ -d "$work_dir/$path/workers" ]] || continue

            if [[ "$path" == "path_a" ]]; then
                echo "--- Path A: EW Sudakov reweight → shower (QCD+QED) ---"
            else
                echo "--- Path B: shower only (QCD) ---"
            fi

            local sc=0 fc=0
            for sf in "$work_dir/$path/workers/"*/status; do
                [[ -f "$sf" ]] || continue
                local s; s=$(cat "$sf")
                if [[ "$s" == SUCCESS* ]]; then ((sc++)); echo "  $(basename "$(dirname "$sf")"): $s"
                elif [[ "$s" == FAILED* ]]; then ((fc++)); echo "  $(basename "$(dirname "$sf")"): $s"
                    local wl="$(dirname "$sf")/worker.log"
                    [[ -f "$wl" ]] && grep -i "error\|failed" "$wl" | tail -3 | sed 's/^/      /'
                fi
            done
            echo ""; echo "  Total: $sc success, $fc failed"; echo ""
        done

        echo "--- Output Files ---"
        for output in "$work_dir"/*/output/*.hepmc.gz; do
            [[ -f "$output" ]] || continue
            local size; size=$(stat -c%s "$output" 2>/dev/null || echo 0)
            echo "  $output ($(numfmt --to=iec-i --suffix=B "$size" 2>/dev/null || echo "${size}B"))"
        done
        echo ""
    } | tee "$report_file"

    log "Summary: $report_file"
}

# =============================================================================
# Main
# =============================================================================
cleanup() {
    log "Received interrupt, cleaning up..."
    [[ -n "${monitor_pid:-}" ]] && kill "$monitor_pid" 2>/dev/null
    jobs -p | xargs -r kill 2>/dev/null
    exit 130
}
trap cleanup SIGINT SIGTERM

run_debug_path_a_reweight() {
    local pipeline_start; pipeline_start=$(date +%s)
    local work_dir="./DYj_pipeline_debug_reweight_qcut${QCUT}_$(date +%Y%m%d_%H%M%S)"
    local full_input="$work_dir/path_a/input/full_input.lhe"
    local debug_input="$work_dir/path_a/input/debug_input.lhe"
    local worker_dir="$work_dir/path_a/workers/worker_debug"

    mkdir -p "$work_dir/path_a/input"

    log "Preparing single-file Path A reweight debug run..."
    if [[ "$INPUT_LHE" == *.gz ]]; then
        gunzip -c "$INPUT_LHE" > "$full_input"
    else
        cp "$INPUT_LHE" "$full_input"
    fi

    local total_events debug_events
    total_events=$(count_events "$full_input")

    if [[ "$DEBUG_MAX_EVENTS" -gt 0 && "$total_events" -gt "$DEBUG_MAX_EVENTS" ]]; then
        debug_input="$work_dir/path_a/input/debug_input_${DEBUG_MAX_EVENTS}ev.lhe"
        extract_first_events "$full_input" "$DEBUG_MAX_EVENTS" "$debug_input"
    else
        debug_input="$work_dir/path_a/input/debug_input.lhe"
        cp "$full_input" "$debug_input"
    fi

    debug_events=$(count_events "$debug_input")

    echo ""
    echo "==========================================================="
    echo " DYj Path A Reweight Debug"
    echo "==========================================================="
    echo " Input:            $INPUT_LHE"
    echo " Source events:    $total_events"
    echo " Debug max events: $DEBUG_MAX_EVENTS"
    echo " Debug input:      $debug_input"
    echo " Debug events:     $debug_events"
    echo " Mode:             path_a (reweight only)"
    echo " PROC_DIR:         $PROC_DIR"
    echo " FxFx Qcut:        $QCUT"
    echo " FxFx njmax:       $NJMAX"
    echo " FxFx qCutList:    ${QCUT_LIST:-<disabled>}"
    echo " Work directory:   $work_dir"
    echo "==========================================================="
    echo ""

    if [[ "$DRY_RUN" == "true" ]]; then
        echo "[DRY RUN] Would run one path_a reweight worker on $debug_input"
        return 0
    fi

    log "Creating shadow PROCNLO directory for path_a..."
    create_shadow_procnlo "$work_dir/path_a/shadow_procnlo" "a"

    log "Setting up debug worker directory..."
    create_worker_procnlo "$worker_dir" "$work_dir/path_a/shadow_procnlo" "0"

    if ! run_path_a_reweight_debug_worker "$debug_input" "$worker_dir"; then
        local total_failed; total_failed=$(get_elapsed "$pipeline_start")
        echo ""
        echo "==========================================================="
        echo " Reweight Debug Failed"
        echo "==========================================================="
        echo " Worker log: $worker_dir/worker.log"
        echo " Reweight log: $worker_dir/reweight.log"
        echo " Time: $(format_duration "$total_failed")"
        echo "==========================================================="
        return 1
    fi

    local total_elapsed; total_elapsed=$(get_elapsed "$pipeline_start")
    local output_file="$worker_dir/output/reweighted_debug.lhe"

    echo ""
    echo "==========================================================="
    echo " Reweight Debug Complete"
    echo "==========================================================="
    echo " Output LHE:   $output_file"
    echo " Worker log:   $worker_dir/worker.log"
    echo " Reweight log: $worker_dir/reweight.log"
    echo " PROCNLO:      $worker_dir/procnlo"
    echo " Total time:   $(format_duration "$total_elapsed")"
    echo "==========================================================="
}

main() {
    if [[ "$DEBUG_PATH_A_REWEIGHT" == "true" ]]; then
        MODE="path_a"
        run_debug_path_a_reweight
        return $?
    fi

    local pipeline_start; pipeline_start=$(date +%s)
    local work_dir="./DYj_pipeline_qcut${QCUT}_$(date +%Y%m%d_%H%M%S)"
    mkdir -p "$work_dir"/{master,path_a/chunks,path_b/chunks}

    # Decompress and analyze input
    log "Analyzing input..."
    if [[ "$INPUT_LHE" == *.gz ]]; then
        log "Decompressing (this may take a while for large files)..."
        gunzip -c "$INPUT_LHE" > "$work_dir/master/events.lhe"
    else
        cp "$INPUT_LHE" "$work_dir/master/events.lhe"
    fi

    local nevents nchunks
    nevents=$(count_events "$work_dir/master/events.lhe")
    nchunks=$(calculate_chunks "$nevents")

    echo ""
    echo "==========================================================="
    echo " EW Sudakov Z+jets Pipeline v1"
    echo "==========================================================="
    echo " Input:            $INPUT_LHE"
    echo " Total events:     $nevents"
    echo " Events per chunk: $EVENTS_PER_CHUNK"
    echo " Number of chunks: $nchunks"
    echo " Mode:             $MODE"
    echo " Max concurrent:   $MAX_CONCURRENT"
    echo " PROC_DIR:         $PROC_DIR"
    echo " FxFx Qcut:        $QCUT"
    echo " FxFx njmax:       $NJMAX"
    echo " FxFx qCutList:    ${QCUT_LIST:-<disabled>}"
    echo " Work directory:   $work_dir"
    echo "==========================================================="
    echo ""
    echo " Path A: reweight(EW Sudakov) → shower(QCD+QED)"
    echo " Path B: shower(QCD only)"
    echo ""

    if [[ "$DRY_RUN" == "true" ]]; then
        echo "[DRY RUN] Would process $nchunks chunks"
        [[ "$MODE" != "path_b" ]] && echo "  Path A: $nchunks workers"
        [[ "$MODE" != "path_a" ]] && echo "  Path B: $nchunks workers"
        return 0
    fi

    # Split
    log "Splitting into $nchunks chunks..."
    split_lhe "$work_dir/master/events.lhe" "$nchunks" "$work_dir/path_a/chunks"

    # Symlink chunks for path_b
    mkdir -p "$work_dir/path_b/chunks"
    for f in "$work_dir/path_a/chunks/"*.lhe; do
        ln -sf "$(realpath "$f")" "$work_dir/path_b/chunks/$(basename "$f")"
    done

    # Create shadow PROCNLOs (different shower cards for a vs b)
    log "Creating shadow PROCNLO directories..."
    [[ "$MODE" != "path_b" ]] && create_shadow_procnlo "$work_dir/path_a/shadow_procnlo" "a"
    [[ "$MODE" != "path_a" ]] && create_shadow_procnlo "$work_dir/path_b/shadow_procnlo" "b"

    # Create worker directories
    log "Setting up worker directories..."
    for i in $(seq 0 $((nchunks-1))); do
        local id; id=$(printf %03d "$i")
        [[ "$MODE" != "path_b" ]] && create_worker_procnlo "$work_dir/path_a/workers/worker_$id" "$work_dir/path_a/shadow_procnlo" "$i"
        [[ "$MODE" != "path_a" ]] && create_worker_procnlo "$work_dir/path_b/workers/worker_$id" "$work_dir/path_b/shadow_procnlo" "$i"
    done

    # Total workers
    local total_workers=0
    [[ "$MODE" == "both" ]] && total_workers=$((nchunks * 2))
    [[ "$MODE" == "path_a" ]] && total_workers=$nchunks
    [[ "$MODE" == "path_b" ]] && total_workers=$nchunks

    # Launch workers with concurrency limit
    log "Launching $total_workers workers (max $MAX_CONCURRENT concurrent)..."

    monitor_progress "$work_dir" "$total_workers" &
    monitor_pid=$!

    local active_pids=()
    for i in $(seq 0 $((nchunks-1))); do
        local id; id=$(printf %03d "$i")
        local chunk="$work_dir/path_a/chunks/chunk_$id.lhe"

        # Concurrency control
        while [[ ${#active_pids[@]} -ge $MAX_CONCURRENT ]]; do
            local new_pids=()
            for pid in "${active_pids[@]}"; do
                if kill -0 "$pid" 2>/dev/null; then
                    new_pids+=("$pid")
                else
                    wait "$pid" 2>/dev/null || true
                fi
            done
            active_pids=("${new_pids[@]}")
            sleep 1
        done

        if [[ "$MODE" != "path_b" ]]; then
            run_path_a_worker "$id" "$chunk" "$work_dir/path_a/workers/worker_$id" &
            active_pids+=($!)
        fi

        if [[ "$MODE" != "path_a" ]]; then
            run_path_b_worker "$id" "$chunk" "$work_dir/path_b/workers/worker_$id" &
            active_pids+=($!)
        fi
    done

    # Wait for all
    for pid in "${active_pids[@]}"; do
        wait "$pid" 2>/dev/null || true
    done

    kill "$monitor_pid" 2>/dev/null || true
    wait "$monitor_pid" 2>/dev/null || true

    # Count results
    local a_success=0 b_success=0 a_failed=0 b_failed=0
    [[ "$MODE" != "path_b" ]] && {
        a_success=$(grep -rl "SUCCESS" "$work_dir/path_a/workers/"*/status 2>/dev/null | wc -l || echo 0)
        a_failed=$(grep -rl "FAILED" "$work_dir/path_a/workers/"*/status 2>/dev/null | wc -l || echo 0)
    }
    [[ "$MODE" != "path_a" ]] && {
        b_success=$(grep -rl "SUCCESS" "$work_dir/path_b/workers/"*/status 2>/dev/null | wc -l || echo 0)
        b_failed=$(grep -rl "FAILED" "$work_dir/path_b/workers/"*/status 2>/dev/null | wc -l || echo 0)
    }

    echo ""
    log "Merging outputs..."
    mkdir -p "$work_dir/path_a/output" "$work_dir/path_b/output"
    [[ "$MODE" != "path_b" ]] && merge_hepmc "$work_dir/path_a" "$work_dir/path_a/output/Zj_FxFx_EWS_qcut${QCUT}.hepmc.gz"
    [[ "$MODE" != "path_a" ]] && merge_hepmc "$work_dir/path_b" "$work_dir/path_b/output/Zj_FxFx_noEWS_qcut${QCUT}.hepmc.gz"

    echo ""
    generate_summary "$work_dir"

    local total_elapsed; total_elapsed=$(get_elapsed "$pipeline_start")

    echo ""
    echo "==========================================================="
    echo " Pipeline Complete!"
    echo "==========================================================="
    [[ "$MODE" != "path_b" ]] && echo " Path A (EWS + QCD+QED shower): $a_success/$nchunks ok, $a_failed failed"
    [[ "$MODE" != "path_a" ]] && echo " Path B (QCD shower only):       $b_success/$nchunks ok, $b_failed failed"
    echo ""
    [[ "$MODE" != "path_b" ]] && echo " Output A: $work_dir/path_a/output/Zj_FxFx_EWS_qcut${QCUT}.hepmc.gz"
    [[ "$MODE" != "path_a" ]] && echo " Output B: $work_dir/path_b/output/Zj_FxFx_noEWS_qcut${QCUT}.hepmc.gz"
    echo " Report:   $work_dir/summary_report.txt"
    echo " Total time: $(format_duration "$total_elapsed")"
    echo "==========================================================="
}

main
