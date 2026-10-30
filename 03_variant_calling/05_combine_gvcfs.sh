#!/bin/bash
# [[03_variant_calling/05_combine_gvcfs.sh]]
# Goal: Combine per-sample GVCFs into a single cohort GVCF for joint genotyping.

#SBATCH --job-name=combine_gvcfs
#SBATCH --output=_data/logs/combine_gvcfs/%x_%j.log
#SBATCH --error=_data/logs/combine_gvcfs/%x_%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=48:00:00

set -e
set -o pipefail

if [ -n "${SLURM_SUBMIT_DIR:-}" ]; then
    export PROJECT_ROOT="$SLURM_SUBMIT_DIR"
else
    export PROJECT_ROOT="${PROJECT_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
fi

eval "$(micromamba shell hook --shell bash 2>/dev/null || mamba shell hook --shell bash 2>/dev/null || conda shell.bash hook 2>/dev/null)"
micromamba activate mdd 2>/dev/null || mamba activate mdd 2>/dev/null || conda activate mdd

RAW_VCF_DIR="$PROJECT_ROOT/_data/variant/raw_vcf"
GENOME_DIR="$PROJECT_ROOT/_data/genome"
FASTA_FILE="$GENOME_DIR/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
LOG_DIR="$PROJECT_ROOT/_data/logs/combine_gvcfs"
SAMPLE_SHEET="$PROJECT_ROOT/_data/_metadata/sample_sheet.csv"

mkdir -p "$LOG_DIR"

# Collect all input GVCFs (one per subject, from input fraction only)
mapfile -t SUBJECTS < <(
  tail -n +2 "$SAMPLE_SHEET" | awk -F',' '{print $1}' | sed 's/_[12]$//' | \
  sed 's/_\(m6A\|input\)$//' | sort -u
)

GVCFS=()
for subj in "${SUBJECTS[@]}"; do
    gvcf="$RAW_VCF_DIR/${subj}_input.g.vcf.gz"
    if [ ! -f "$gvcf" ]; then
        gvcf="$RAW_VCF_DIR/${subj}.g.vcf.gz"
    fi

    if [ -f "$gvcf" ]; then
        GVCFS+=("--variant $gvcf")
    else
        echo "[WARNING] Missing GVCF for $subj: checked ${subj}_input.g.vcf.gz and ${subj}.g.vcf.gz" >&2
    fi
done

if [ ${#GVCFS[@]} -eq 0 ]; then
    echo "Error: No GVCFs found in $RAW_VCF_DIR." >&2
    exit 1
fi

COMBINED_GVCF="$RAW_VCF_DIR/cohort_combined.g.vcf.gz"

if [ -f "$COMBINED_GVCF" ]; then
    echo "[SKIP] Combined GVCF already exists: $COMBINED_GVCF"
    exit 0
fi

echo ">>> Combining ${#GVCFS[@]} GVCFs..."
gatk CombineGVCFs \
    -R "$FASTA_FILE" \
    ${GVCFS[@]} \
    -O "$COMBINED_GVCF" \
    --tmp-dir "$LOG_DIR" 2>&1 | tee "$LOG_DIR/combine_gvcfs.log"

echo ">>> CombineGVCFs done → $COMBINED_GVCF"
