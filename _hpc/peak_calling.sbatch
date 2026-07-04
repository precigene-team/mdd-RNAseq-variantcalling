#!/bin/bash
# [[_hpc/peak_calling.sbatch]]
# Goal: Slurm submission script to execute MACS2 peak-calling on the HPC in parallel.
#SBATCH --job-name=MDD_Peak_Calling
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=_data/logs/peak_call/peak_call_log_%j.txt
#SBATCH --error=_data/logs/peak_call/peak_call_error_%j.txt

# Exit on error
set -e

# --- Project Root Discovery ---
find_project_root() {
    local curr="$(pwd)"
    while [ "$curr" != "/" ] && [ -n "$curr" ]; do
        if [ -f "$curr/run" ]; then
            echo "$curr"
            return 0
    fi
        curr="$(dirname "$curr")"
done
    echo "$(pwd)" # Fallback
}
export PROJECT_ROOT="${PROJECT_ROOT:-$(find_project_root)}"
# ------------------------------

CSV_FILE="$PROJECT_ROOT/_data/_metadata/sample_sheet.csv"
BAM_DIR="$PROJECT_ROOT/_data/bam"
PEAK_DIR="$PROJECT_ROOT/_data/peaks"
LOG_DIR="$PROJECT_ROOT/_data/logs/hpc_run"
PEAK_LOG_DIR="$PROJECT_ROOT/_data/logs/peak_call/samples"

mkdir -p "$PEAK_DIR"
mkdir -p "$LOG_DIR"
mkdir -p "$PEAK_LOG_DIR"

echo "========================================================="
echo "   HPC Peak Calling via MACS2"
echo "   Starting at: $(date)"
echo "========================================================="

# 1. Initialize Conda/Mamba
eval "$(micromamba shell hook --shell bash 2>/dev/null || mamba shell hook --shell bash 2>/dev/null || conda shell.bash hook 2>/dev/null)"
micromamba activate mdd 2>/dev/null || mamba activate mdd 2>/dev/null || conda activate mdd

# Set thread limits for underlying math libraries to align with 4-way parallelism (16 CPUs total)
export OMP_NUM_THREADS=4
export MKL_NUM_THREADS=4
export OPENBLAS_NUM_THREADS=4
export NUMEXPR_NUM_THREADS=4
export VECLIB_MAXIMUM_THREADS=4

# Check if macs2 is installed
if ! command -v macs2 &> /dev/null; then
    echo "ERROR: macs2 not found in active environment. Please install it with: mamba install -y macs2" >&2
    exit 1
fi

# Gated subjects exclusion list
is_gated() {
    local subj=$1
    if [[ "$subj" == "female_mdd_2" || "$subj" == "female_ctrl_6" || "$subj" == "male_ctrl_1" ]]; then
        return 0 # is gated
    fi
    return 1 # is not gated
}

# 2. Iterate through all unique subjects
echo "[INFO] Parsing sample sheet..."
subjects=$(awk -F',' 'NR>1 {print $1}' "$CSV_FILE" | sed -E 's/_(m6A|input)_1//' | sort -u)

for subj in $subjects; do
    if is_gated "$subj"; then
        echo "  [SKIP] Subject $subj is gated (low coverage)."
        continue
    fi

    # Define BAM and Peak output paths
    IP_BAM="$BAM_DIR/${subj}_m6A/${subj}_m6A_Aligned.sortedByCoord.out.bam"
    INPUT_BAM="$BAM_DIR/${subj}_input/${subj}_input_Aligned.sortedByCoord.out.bam"
    OUT_PEAKS="${PEAK_DIR}/${subj}/${subj}_peaks.narrowPeak"

    if [ -s "$OUT_PEAKS" ]; then
        echo "  [SKIP] Subject $subj (Peaks already exist)"
        continue
    fi

    if [ ! -f "$IP_BAM" ] || [ ! -f "$INPUT_BAM" ]; then
        echo "  [WARNING] BAM files for $subj not found. Skipping."
        echo "    IP: $IP_BAM"
        echo "    Input: $INPUT_BAM"
        continue
    fi

    echo ">>> Running MACS2 for subject: $subj"
    # Call peaks:
    # -t: IP BAM (treatment)
    # -c: Input BAM (control)
    # -f BAM: Treat alignments as single-end reads to bypass intron-splicing size calculation distortion
    # -g hs: Human genome size
    # -q 0.05: FDR threshold
    # --keep-dup all: Keep all duplicates to preserve quantitative epitranscriptomic signal on highly expressed genes
    # --nomodel --extsize 150: Force fixed 150bp extension to prevent intron-splicing distortion and speed up q-value lookup table computation
    macs2 callpeak \
        -t "$IP_BAM" \
        -c "$INPUT_BAM" \
        -f BAM \
        -g hs \
        -n "$subj" \
        --outdir "$PEAK_DIR/$subj" \
        --keep-dup all \
        --slocal 10000 \
        --llocal 50000 \
        -q 0.05 \
        --nomodel \
        --extsize 150 > "$PEAK_LOG_DIR/${subj}.macs2.log" 2>&1 &

    # Limit parallelism to 4 concurrent jobs
    # 4 jobs * 4 threads/job = 16 threads (matching --cpus-per-task=16)
    while [ $(jobs -rp | wc -l) -ge 4 ]; do
        sleep 5
    done
done

# Wait for all background peak-calling jobs to complete
wait

echo "========================================================="
echo "   MACS2 Peak Calling COMPLETED SUCCESSFULLY"
echo "   Ending at: $(date)"
echo "========================================================="
