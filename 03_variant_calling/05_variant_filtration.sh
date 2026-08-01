#!/bin/bash
# [[03_variant_calling/05_variant_filtration.sh]]
# Goal: Filter the joint-genotyped cohort VCF and split into per-sample VCFs.

#SBATCH --job-name=variant_filtration
#SBATCH --output=_data/logs/variant_filtration/%x_%j.log
#SBATCH --error=_data/logs/variant_filtration/%x_%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=12:00:00

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
OUTPUT_DIR="$PROJECT_ROOT/_data/variant/filtered_vcf"
RAW_SNPS_DIR="$PROJECT_ROOT/_data/variant/snps_raw"
LOG_DIR="$PROJECT_ROOT/_data/logs/variant_filtration"
GENOME_DIR="$PROJECT_ROOT/_data/genome"
FASTA_FILE="$GENOME_DIR/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
SAMPLE_SHEET="$PROJECT_ROOT/_data/_metadata/sample_sheet.csv"

mkdir -p "$OUTPUT_DIR" "$RAW_SNPS_DIR" "$LOG_DIR"

JOINT_VCF="$RAW_VCF_DIR/cohort_joint.vcf.gz"
COHORT_SNPS="$RAW_SNPS_DIR/cohort_snps.raw.vcf.gz"
COHORT_FILTERED="$RAW_VCF_DIR/cohort_filtered.vcf.gz"

if [ ! -f "$JOINT_VCF" ]; then
    echo "Error: Joint VCF not found at $JOINT_VCF. Run 04c_genotype_gvcfs.sh first." >&2
    exit 1
fi

# Skip if per-sample outputs already exist (check first subject as sentinel)
mapfile -t SUBJECTS < <(
  tail -n +2 "$SAMPLE_SHEET" | awk -F',' '{print $1}' | sed 's/_[12]$//' | \
  sed 's/_\(m6A\|input\)$//' | sort -u
)
SENTINEL="$OUTPUT_DIR/${SUBJECTS[0]}_input.filtered.vcf.gz"
if [ -f "$SENTINEL" ]; then
    echo "[SKIP] Per-sample filtered VCFs already exist."
    exit 0
fi

echo ">>> Selecting SNPs from joint VCF..."
gatk SelectVariants \
    -R "$FASTA_FILE" \
    -V "$JOINT_VCF" \
    -select-type SNP \
    --restrict-alleles-to BIALLELIC \
    -O "$COHORT_SNPS" 2>&1 | tee "$LOG_DIR/select_snps.log"

echo ">>> Hard-filtering SNPs..."
gatk VariantFiltration \
    -R "$FASTA_FILE" \
    -V "$COHORT_SNPS" \
    -O "$COHORT_FILTERED" \
    --window 35 \
    --cluster 3 \
    --filter-name "FS" \
    --filter-expression "FS > 30.0" \
    --filter-name "QD" \
    --filter-expression "QD < 2.0" 2>&1 | tee "$LOG_DIR/filter.log"

echo ">>> Splitting into per-sample VCFs..."
# Get list of samples in the joint VCF
bcftools query -l "$COHORT_FILTERED" > "$LOG_DIR/sample_list.txt"

while IFS= read -r sample; do
    out_vcf="$OUTPUT_DIR/${sample}.filtered.vcf.gz"
    if [ -f "$out_vcf" ]; then
        echo "  [SKIP] $sample already extracted"
        continue
    fi
    echo "  Extracting $sample (PASS, het, DP>=10, GQ>=20)..."
    bcftools view -s "$sample" -f PASS -g het -i 'FMT/DP >= 10 && FMT/GQ >= 20' -O z -o "$out_vcf" "$COHORT_FILTERED"
    bcftools index -t "$out_vcf"
done < "$LOG_DIR/sample_list.txt"

echo ">>> Variant filtration complete — per-sample VCFs in $OUTPUT_DIR"
