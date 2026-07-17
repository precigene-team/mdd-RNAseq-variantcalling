#!/bin/bash
# [[01_core_alignment/10_tress_dmr_finder.sh]]
# Goal: Submit TRESS_DMRfinder jobs for male and female cohorts via Slurm,
#       or run locally by passing the sex argument directly.

# Steps:
#   1. Resolve project root (Slurm or local).
#   2. Activate mdd conda environment.
#   3. Determine sex from Slurm array index (0=male, 1=female) or argument.
#   4. Check if DMR output already exists (skip if so).
#   5. Run 10_tress_dmr_finder.R <sex>.

#SBATCH --job-name=tress_dmr_finder
#SBATCH --output=_data/logs/tress_dmr/%x_%A_%a.log
#SBATCH --error=_data/logs/tress_dmr/%x_%A_%a.err
#SBATCH --array=0-1
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --time=24:00:00

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

LOG_DIR="$PROJECT_ROOT/_data/logs/tress_dmr"
DMR_DIR="$PROJECT_ROOT/_data/dmr"
mkdir -p "$LOG_DIR"
mkdir -p "$DMR_DIR"

# --- Resolve sex ---
SEXES=("male" "female")

if [ -n "${SLURM_ARRAY_TASK_ID:-}" ]; then
    SEX="${SEXES[$SLURM_ARRAY_TASK_ID]}"
    if [ -z "$SEX" ]; then
        echo "Error: No sex resolved for SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID" >&2
        exit 1
    fi
else
    # Local mode
    SEX="${1:-}"
    if [ -z "$SEX" ]; then
        echo "ERROR: Sex argument required." >&2
        echo "Usage: bash 10_tress_dmr_finder.sh <male|female>" >&2
        exit 1
    fi
fi

echo ">>> Running TRESS_DMRfinder for sex: $SEX"

# --- Skip if DMR output already exists ---
OUT_DMR="$DMR_DIR/${SEX}/DMRs.RData"
if [ -s "$OUT_DMR" ]; then
    echo "[Skip] DMRs.RData already exists for $SEX — delete to rerun."
    exit 0
fi

# --- Run R script ---
if [ -n "${SLURM_ARRAY_TASK_ID:-}" ]; then
    Rscript "$PROJECT_ROOT/01_core_alignment/10_tress_dmr_finder.R" "$SEX"
else
    sample_log="$LOG_DIR/${SEX}_dmr_finder.log"
    Rscript "$PROJECT_ROOT/01_core_alignment/10_tress_dmr_finder.R" "$SEX" 2>&1 | tee "$sample_log"
fi

echo ">>> TRESS_DMRfinder complete for sex: $SEX"
