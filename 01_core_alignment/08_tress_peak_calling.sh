#!/bin/bash
# [[01_core_alignment/08_tress_peak_calling.sh]]
# Goal: Run TRESS peak calling for a single subject (supports Slurm Job Array and local execution).

# Steps:
#   1. Resolve subject name from Slurm array index (0-38) or first argument if local.
#   2. Check if output peaks.xls file already exists.
#   3. Check for local locks to prevent multiple processes overlapping on the same sample.
#   4. Execute Rscript 08_tress_peak_calling.R <subject_name>.

#SBATCH --job-name=tress_peak_calling
#SBATCH --output=_data/logs/tress_peak/%x_%A_%a.log
#SBATCH --error=_data/logs/tress_peak/%x_%A_%a.err
#SBATCH --array=0-38%2
#SBATCH --cpus-per-task=2
#SBATCH --mem=128G
#SBATCH --time=04:00:00

set -e
set -o pipefail

if [ -n "${SLURM_SUBMIT_DIR:-}" ]; then
    export PROJECT_ROOT="$SLURM_SUBMIT_DIR"
else
    export PROJECT_ROOT="${PROJECT_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
fi

# --- Safe Conda/Mamba Hook Initialization ---
eval "$(micromamba shell hook --shell bash 2>/dev/null || mamba shell hook --shell bash 2>/dev/null || conda shell.bash hook 2>/dev/null)"
micromamba activate mdd 2>/dev/null || mamba activate mdd 2>/dev/null || conda activate mdd

if [ -z "$CONDA_PREFIX" ]; then
    echo "Assertion Error: Mamba environment 'mdd' is not active." >&2
    exit 1
fi

CSV_FILE="$PROJECT_ROOT/_data/_metadata/sample_sheet.csv"
PEAKS_DIR="$PROJECT_ROOT/_data/peaks"
LOG_DIR="$PROJECT_ROOT/_data/logs/tress_peak"

mkdir -p "$PEAKS_DIR"
mkdir -p "$LOG_DIR"

# --- Resolve subject name ---
if [ -n "${SLURM_ARRAY_TASK_ID:-}" ]; then
    # Slurm Job Array Mode
    mapfile -t SUBJECTS < <(
      tail -n +2 "$CSV_FILE" | awk -F',' '{print $1}' | sed 's/_[12]$//' | \
      sed 's/_\(m6A\|input\)$//' | sort -u
    )
    SUBJECT_NAME="${SUBJECTS[$SLURM_ARRAY_TASK_ID]}"
    if [ -z "$SUBJECT_NAME" ]; then
        echo "Error: No subject resolved for SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID" >&2
        exit 1
    fi
else
    # Local Mode: argument must be passed
    SUBJECT_NAME="$1"
    if [ -z "$SUBJECT_NAME" ]; then
        echo "ERROR: SUBJECT_NAME argument is missing." >&2
        echo "Usage: bash 08_tress_peak_calling.sh <subject_name> [cpus] [memory]" >&2
        exit 1
    fi
fi

out_xls="$PEAKS_DIR/${SUBJECT_NAME}/${SUBJECT_NAME}_peaks.xls"
sample_log="$LOG_DIR/${SUBJECT_NAME}.peak.log"

# Register trap to clean up incomplete output directory in case of job cancellation (SIGTERM/SIGINT)
cleanup_leftovers() {
    echo ">>> [Signal Trap] Job interrupted. Cleaning up incomplete output directory for subject: $SUBJECT_NAME"
    rm -rf "$PEAKS_DIR/$SUBJECT_NAME"
    exit 1
}
trap cleanup_leftovers SIGTERM SIGINT

# Check if local workstation is processing this sample (lock file exists)
lock_file="$LOG_DIR/${SUBJECT_NAME}.local_lock"
if [ -f "$lock_file" ]; then
    echo "[Skip] $SUBJECT_NAME (Local workstation is processing this sample; locked)"
    exit 0
fi

if [ -s "$out_xls" ]; then
    echo "[Skip] $SUBJECT_NAME (TRESS peaks already exist)"
    # Clean up own log and error files in Slurm mode to prevent log clutter
    if [ -n "${SLURM_ARRAY_JOB_ID:-}" ] && [ -n "${SLURM_ARRAY_TASK_ID:-}" ]; then
        rm -f "$LOG_DIR/tress_peak_calling_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.log"
        rm -f "$LOG_DIR/tress_peak_calling_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err"
    fi
    exit 0
fi

echo ">>> Running TRESS peak calling for subject: $SUBJECT_NAME"

# Run R script
if [ -n "${SLURM_ARRAY_TASK_ID:-}" ]; then
    # In Slurm Array mode, print output to Slurm's log file directly
    Rscript "$PROJECT_ROOT/01_core_alignment/08_tress_peak_calling.R" "$SUBJECT_NAME"
else
    # In local mode, redirect output to the sample log file
    Rscript "$PROJECT_ROOT/01_core_alignment/08_tress_peak_calling.R" "$SUBJECT_NAME" > "$sample_log" 2>&1
fi
