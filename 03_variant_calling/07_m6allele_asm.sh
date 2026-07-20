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
#SBATCH --mem=32G
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


# --- Paths ---
M6ALLELE_JAR="$PROJECT_ROOT/_data/tools/M6Allele/M6Allele.jar"
PEAKS_BED="$PROJECT_ROOT/_data/peaks/consensus_peaks.bed"
GTF_FILE="$PROJECT_ROOT/_data/genome/Homo_sapiens.GRCh38.115.gtf"
VCF_DIR="$PROJECT_ROOT/_data/variant/filtered_vcf"
SPLIT_DIR="$PROJECT_ROOT/_data/variant/split_n_cigar"
ASM_DIR="$PROJECT_ROOT/_data/asm"
LOG_DIR="$PROJECT_ROOT/_data/logs/asm"
SAMPLE_SHEET="$PROJECT_ROOT/_data/_metadata/sample_sheet.csv"

mkdir -p "$ASM_DIR" "$LOG_DIR"

# --- Compute Java heap (75% of Slurm allocation, default 8G locally) ---
if [ -n "${SLURM_MEM_PER_NODE:-}" ]; then
    JAVA_HEAP_MB=$(( SLURM_MEM_PER_NODE * 3 / 4 ))
else
    JAVA_HEAP_MB=8192
fi

# --- Auto-build M6Allele.jar if missing ---
if [ ! -f "$M6ALLELE_JAR" ]; then
    TOOLS_DIR="$(dirname "$M6ALLELE_JAR")"
    if [ ! -d "$TOOLS_DIR" ]; then
        echo ">>> Cloning M6Allele repository..."
        git clone https://github.com/RenLabBioinformatics/M6Allele.git "$TOOLS_DIR"
    fi
    echo ">>> Building M6Allele.jar with Maven..."
    cd "$TOOLS_DIR"
    mvn clean package -DskipTests -q
    BUILT_JAR=$(find target -name "M6Allele-*-jar-with-dependencies.jar" 2>/dev/null | head -1)
    if [ -z "$BUILT_JAR" ]; then
        echo "Error: Maven build completed but no fat JAR found in target/." >&2
        exit 1
    fi
    cp "$BUILT_JAR" "$M6ALLELE_JAR"
    echo ">>> M6Allele.jar built and placed at $M6ALLELE_JAR"
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
    ASM_TMP="$ASM_DIR/tmp"
    WORK_DIR="$ASM_TMP/${SUBJ}_work"

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

    mkdir -p "$WORK_DIR"

    # All outputs inside unique WORK_DIR to prevent cross-subject overwrites:
    # M6Allele writes aseRes.txt, peak_with_snp.txt, snp_location.txt to the
    # same directory as -o. If -o pointed to a shared dir, parallel array jobs
    # would corrupt each other.
    ASM_OUT="$WORK_DIR/${SUBJ}_asmPeak.txt"
    FIXED_VCF="$WORK_DIR/${SUBJ}_fixed.vcf"

    # Always clean up WORK_DIR on exit (skip, success, crash)
    trap 'rm -rf "$WORK_DIR"' EXIT

    # Skip if final output already exists
    FINAL_OUT="$ASM_DIR/${SUBJ}_asmPeak.txt"
    if [ -s "$FINAL_OUT" ]; then
        echo "[SKIP] ASM analysis for $SUBJ already completed."
        return 0
    fi

    # --- Filter VCF for M6Allele ---
    # M6Allele's FileInputStream does not handle gzipped files, and its
    # getMajorMinorHaplotype parser crashes (ArrayIndexOutOfBoundsException)
    # on multi-allelic SNPs. Filter to PASS, biallelic, heterozygous only.
    echo ">>> Filtering VCF to biallelic het SNPs for M6Allele..."
    bcftools view -f PASS -v snps -m2 -M2 -i 'GT="het"' -O v "$VCF_FILE" > "$FIXED_VCF"

    if [ ! -s "$FIXED_VCF" ]; then
        echo "[SKIP] No heterozygous biallelic PASS SNPs found for $SUBJ."
        touch "$FINAL_OUT"
        return 0
    fi

    echo ">>> Running M6Allele.jar AsmPeakDetection (workdir: $WORK_DIR)..."
    cd "$WORK_DIR"
    if java -Xmx${JAVA_HEAP_MB}M -jar "$M6ALLELE_JAR" AsmPeakDetection \
        -g "$GTF_FILE" \
        -bed "$PEAKS_BED" \
        -vcf "$FIXED_VCF" \
        -inputBam "$INPUT_BAM" \
        -inputBai "$INPUT_BAI" \
        -ipBam "$IP_BAM" \
        -ipBai "$IP_BAI" \
        -o "$ASM_OUT"; then
        if [ ! -s "$ASM_OUT" ]; then
            echo "Assertion Error: M6Allele did not generate non-empty output for $SUBJ." >&2
            return 1
        fi
        mv "$ASM_OUT" "$FINAL_OUT"
        echo ">>> Completed ASM for $SUBJ → $FINAL_OUT"
    else
        echo "[WARNING] M6Allele crashed for $SUBJ (internal Java error). Skipping." >&2
        touch "$FINAL_OUT"
    fi
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
