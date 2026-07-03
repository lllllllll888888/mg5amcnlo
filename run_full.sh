#!/usr/bin/env bash
# Convenience launcher for parallel_pipeline_DYj.sh on the standard input file.
# Usage:
#   ./run_full.sh                        # uses default input
#   ./run_full.sh /path/to/other.lhe.gz  # override input
# Writes a timestamped log next to this script and backgrounds the run.

set -euo pipefail
cd "$(dirname "$(readlink -f "$0")")"

INPUT="${1:-/data0/lenartj/MadGraphDevelopment/EWSudakov/pc12seed1001.lhe.gz}"
LOG="pipeline_full_$(date +%Y%m%d_%H%M%S).log"

if [[ ! -f "$INPUT" ]]; then
  echo "ERROR: input LHE not found: $INPUT" >&2
  exit 1
fi

nohup ./parallel_pipeline_DYj.sh "$INPUT" > "$LOG" 2>&1 &
PID=$!
echo "Launched."
echo "  PID:   $PID"
echo "  Input: $INPUT"
echo "  Log:   $(pwd)/$LOG"
echo
echo "To monitor:  tail -f $LOG"
echo "To check:    ps -fp $PID"
echo "To kill:     kill $PID  (or: kill -TERM $PID)"
