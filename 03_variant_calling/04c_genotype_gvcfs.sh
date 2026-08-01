#!/bin/bash
# [[03_variant_calling/04c_genotype_gvcfs.sh]]
# Goal: Joint-genotype the combined cohort GVCF to produce a multi-sample VCF.

#SBATCH --job-name=genotype_gvcfs
#SBATCH --output=_data/logs/genotype_gvcfs/%x_%j.log
#SBATCH --error=_data/logs/genotype_gvcfs/%x_%j.err
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
LOG_DIR="$PROJECT_ROOT/_data/logs/genotype_gvcfs"

mkdir -p "$LOG_DIR"

COMBINED_GVCF="$RAW_VCF_DIR/cohort_combined.g.vcf.gz"
JOINT_VCF="$RAW_VCF_DIR/cohort_joint.vcf.gz"

if [ ! -f "$COMBINED_GVCF" ]; then
    echo "Error: Combined GVCF not found at $COMBINED_GVCF. Run 04b_combine_gvcfs.sh first." >&2
    exit 1
fi

if [ -f "$JOINT_VCF" ]; then
    echo "[SKIP] Joint VCF already exists: $JOINT_VCF"
    exit 0
fi

echo ">>> Genotyping cohort GVCF..."
gatk GenotypeGVCFs \
    -R "$FASTA_FILE" \
    -V "$COMBINED_GVCF" \
    -O "$JOINT_VCF" \
    --tmp-dir "$LOG_DIR" 2>&1 | tee "$LOG_DIR/genotype_gvcfs.log"

echo ">>> GenotypeGVCFs done → $JOINT_VCF"
