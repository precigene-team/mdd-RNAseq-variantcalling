#!/bin/bash
# [[_hpc/submit_notify.sh]]
# Goal: Run notification directly if target job is complete, otherwise submit notify.sbatch with dependency.

# Steps:
#   1. Parse and validate target job ID to monitor.
#   2. Query current execution status of target job.
#   3. If running or pending, submit notify.sbatch with --dependency=afterany constraint.
#   4. If already completed/failed, run notify.sbatch directly in the shell.

# Exit on error
set -e
set -o pipefail

# --- Project Root Discovery ---
export PROJECT_ROOT="${PROJECT_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"

# --- Input Parsing & Assertions ---
MONITOR_JOB="$1"
if [ -z "$MONITOR_JOB" ]; then
    echo "Assertion Error: No job ID provided to monitor." >&2
    echo "Usage: bash _hpc/submit_notify.sh <JOB_ID>" >&2
    exit 1
fi

# Assert that MONITOR_JOB is numeric
if [[ ! "$MONITOR_JOB" =~ ^[0-9]+$ ]]; then
    echo "Assertion Error: Monitored job ID must be numeric. Got: '$MONITOR_JOB'" >&2
    exit 1
fi

# Query the job state from sacct
exit_status=$(sacct -j "$MONITOR_JOB" --format=state -n | head -n 1 | xargs)

if [[ "$exit_status" == "RUNNING" || "$exit_status" == "PENDING" ]]; then
    echo "Job $MONITOR_JOB is still active ($exit_status). Submitting held notify.sbatch job..."
    sbatch --dependency=afterany:"$MONITOR_JOB" "$PROJECT_ROOT/_hpc/notify.sbatch" "$MONITOR_JOB"
else
    echo "Job $MONITOR_JOB has already completed/finished ($exit_status). Sending notification immediately..."
    bash "$PROJECT_ROOT/_hpc/notify.sbatch" "$MONITOR_JOB"
fi
