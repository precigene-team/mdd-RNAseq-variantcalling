#!/bin/bash
# [[03_variant_calling/07_m6allele_asm.sh]]
# Goal: Run M6Allele AsmPeakDetection on preprocessed IP/Input BAM pairs.

# Steps:
#   1. Validate inputs (M6Allele.jar, split BAMs, BAI indexes, filtered VCFs, consensus peaks).
#   2. Resolve subject from Slurm array index (0-38) or run loop if run locally.
#   3. Execute java -jar M6Allele.jar AsmPeakDetection per subject.
#   4. Validate that the output asmPeak.txt results exist.

#SBATCH --job-name=m6allele_asm
#SBATCH --output=_data/logs/asm/%x_%A_%a.log
#SBATCH --error=_data/logs/asm/%x_%A_%a.err
#SBATCH --array=0-38
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00

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


# --- Paths ---
M6ALLELE_JAR="$PROJECT_ROOT/_data/tools/M6Allele/M6Allele.jar"
PEAKS_BED="$PROJECT_ROOT/_data/peaks/consensus_peaks.bed"
GTF_FILE="$PROJECT_ROOT/_data/genome/Homo_sapiens.GRCh38.115.gtf"
VCF_DIR="$PROJECT_ROOT/_data/variant/filtered_vcf"
SPLIT_DIR="$PROJECT_ROOT/_data/variant/split_n_cigar"
ASM_DIR="$PROJECT_ROOT/results/tables/12_asm"
LOG_DIR="$PROJECT_ROOT/_data/logs/asm"
SAMPLE_SHEET="$PROJECT_ROOT/_data/_metadata/sample_sheet_original.csv"

mkdir -p "$ASM_DIR" "$LOG_DIR"

# --- Validate M6Allele.jar exists ---
if [ ! -f "$M6ALLELE_JAR" ]; then
    echo "Assertion Error: M6Allele.jar was not found at: $M6ALLELE_JAR" >&2
    echo "Please clone the official repository to download it first:" >&2
    echo "  git clone https://github.com/RenLabBioinformatics/M6Allele.git $PROJECT_ROOT/_data/tools/M6Allele" >&2
    exit 1
fi

# --- Resolve subject(s) to process ---
mapfile -t SUBJECTS < <(
  tail -n +2 "$SAMPLE_SHEET" | awk -F',' '{print $1}' | sed 's/_[12]$//' | \
  sed 's/_\(m6A\|input\)$//' | sort -u
)

run_subject_asm() {
    local SUBJ="$1"
    echo "================================================================================"
    echo "   Processing ASM for: $SUBJ"
    echo "================================================================================"

    # Construct input paths from preprocessed split_n_cigar directory
    INPUT_BAM="$SPLIT_DIR/${SUBJ}_input.split.bam"
    IP_BAM="$SPLIT_DIR/${SUBJ}_m6A.split.bam"
    VCF_FILE="$VCF_DIR/${SUBJ}_input.filtered.vcf.gz"
    SUBJ_OUT_DIR="$ASM_DIR/${SUBJ}"

    # Resolve BAI files (can be named .split.bai or .split.bam.bai)
    INPUT_BAI="$SPLIT_DIR/${SUBJ}_input.split.bai"
    [ ! -f "$INPUT_BAI" ] && INPUT_BAI="${INPUT_BAM}.bai"
    
    IP_BAI="$SPLIT_DIR/${SUBJ}_m6A.split.bai"
    [ ! -f "$IP_BAI" ] && IP_BAI="${IP_BAM}.bai"

    # Input Assertions
    for f in "$INPUT_BAM" "$INPUT_BAI" "$IP_BAM" "$IP_BAI" "$VCF_FILE" "$PEAKS_BED" "$GTF_FILE"; do
        if [ ! -f "$f" ]; then
            echo "Assertion Error: Required file missing for $SUBJ: $f" >&2
            return 1
        fi
    done

    mkdir -p "$SUBJ_OUT_DIR"

    # Skip if outputs already exist and are non-empty
    if [ -s "$SUBJ_OUT_DIR/asmPeak.txt" ]; then
        echo "[SKIP] ASM analysis for $SUBJ already completed."
        return 0
    fi

    echo ">>> Running M6Allele.jar AsmPeakDetection..."
    java -jar "$M6ALLELE_JAR" AsmPeakDetection \
        -g "$GTF_FILE" \
        -bed "$PEAKS_BED" \
        -vcf "$VCF_FILE" \
        -inputBam "$INPUT_BAM" \
        -inputBai "$INPUT_BAI" \
        -ipBam "$IP_BAM" \
        -ipBai "$IP_BAI" \
        -o "$SUBJ_OUT_DIR"

    # Validate output
    if [ ! -s "$SUBJ_OUT_DIR/asmPeak.txt" ]; then
        echo "Assertion Error: M6Allele did not generate non-empty output for $SUBJ." >&2
        return 1
    fi
    echo ">>> Completed ASM for $SUBJ."
}

# --- Execution ---
if [ -n "${SLURM_ARRAY_TASK_ID:-}" ]; then
    # Array job mode
    SUBJECT="${SUBJECTS[$SLURM_ARRAY_TASK_ID]}"
    if [ -z "$SUBJECT" ]; then
        echo "Error: No subject resolved for SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID" >&2
        exit 1
    fi
    run_subject_asm "$SUBJECT"
else
    # Local loop mode
    echo "Running sequentially for all subjects locally..."
    for SUBJ in "${SUBJECTS[@]}"; do
        run_subject_asm "$SUBJ" || exit 1
    done
fi
