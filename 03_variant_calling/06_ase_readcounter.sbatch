#!/bin/bash
# [[03_variant_calling/06_ase_readcounter.sbatch]]
# Goal: Run GATK ASEReadCounter on IP BAMs at heterozygous SNP positions from Input VCFs.

# Steps:
#   1. Validate inputs (IP BAM, filtered VCF, reference FASTA).
#   2. Subset VCF to PASS-only, biallelic, heterozygous SNPs.
#   3. Run gatk ASEReadCounter on the IP BAM.
#   4. Assert output CSV is non-empty.

#SBATCH --job-name=ase_readcounter
#SBATCH --output=_data/logs/ase/%x_%A_%a.log
#SBATCH --error=_data/logs/ase/%x_%A_%a.err
#SBATCH --array=0-38
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
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

# --- Paths ---
REF_FA="$PROJECT_ROOT/_data/genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
VCF_DIR="$PROJECT_ROOT/_data/variant/filtered_vcf"
BAM_DIR="$PROJECT_ROOT/_data/bam"
ASE_DIR="$PROJECT_ROOT/_data/ase"
LOG_DIR="$PROJECT_ROOT/_data/logs/ase"
CSV_DIR="$PROJECT_ROOT/_data/metadata"
SAMPLE_SHEET="$PROJECT_ROOT/_data/_metadata/sample_sheet_original.csv"

mkdir -p "$ASE_DIR" "$LOG_DIR"

# --- Resolve subject from SLURM array index ---
# Build list of unique subjects (strip _m6A / _input suffix)
mapfile -t SUBJECTS < <(
  tail -n +2 "$SAMPLE_SHEET" | awk -F',' '{print $1}' | sed 's/_[12]$//' | \
  sed 's/_\(m6A\|input\)$//' | sort -u
)

SUBJECT="${SUBJECTS[$SLURM_ARRAY_TASK_ID]}"
if [ -z "$SUBJECT" ]; then
    echo "Error: No subject found for SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID" >&2
    exit 1
fi

echo ">>> Processing subject: $SUBJECT"

# --- Input paths ---
IP_BAM="$PROJECT_ROOT/_data/variant/split_n_cigar/${SUBJECT}_m6A.split.bam"
FILTERED_VCF="$VCF_DIR/${SUBJECT}_input.filtered.vcf.gz"
HET_VCF="$ASE_DIR/${SUBJECT}_het_snps.vcf.gz"
OUT_CSV="$ASE_DIR/${SUBJECT}_ip_ase_counts.csv"

# --- Input Assertions ---
if [ ! -f "$IP_BAM" ]; then
    echo "Assertion Error: IP BAM not found: $IP_BAM" >&2
    exit 1
fi
if [ ! -f "$FILTERED_VCF" ]; then
    echo "Assertion Error: Filtered VCF not found: $FILTERED_VCF" >&2
    exit 1
fi

# Skip if already done
if [ -s "$OUT_CSV" ]; then
    echo "[SKIP] $SUBJECT: ASE counts already exist."
    exit 0
fi

# --- Step 1: Subset VCF to PASS, biallelic, heterozygous SNPs ---
echo ">>> Step 1: Filtering VCF to heterozygous biallelic SNPs..."
gatk SelectVariants \
    -R "$REF_FA" \
    -V "$FILTERED_VCF" \
    --select-type-to-include SNP \
    --restrict-alleles-to BIALLELIC \
    --select "vc.getGenotype(0).isHet()" \
    --exclude-filtered \
    -O "$HET_VCF"

# Index if needed
if [ ! -f "${HET_VCF}.tbi" ]; then
    gatk IndexFeatureFile -I "$HET_VCF"
fi

HET_COUNT=$(gatk CountVariants -V "$HET_VCF" 2>/dev/null | tail -1)
echo "    Heterozygous biallelic SNPs: $HET_COUNT"

if [ "$HET_COUNT" -eq 0 ]; then
    echo "[WARNING] No heterozygous SNPs found for $SUBJECT. Skipping ASEReadCounter." >&2
    exit 0
fi

# --- Step 2: Run ASEReadCounter ---
echo ">>> Step 2: Running GATK ASEReadCounter on IP BAM..."
gatk ASEReadCounter \
    -R "$REF_FA" \
    -I "$IP_BAM" \
    -V "$HET_VCF" \
    -O "$OUT_CSV" \
    --min-mapping-quality 20 \
    --min-base-quality 20 \
    --disable-read-filter NotDuplicateReadFilter

# --- Output Assertion ---
if [ ! -s "$OUT_CSV" ]; then
    echo "Assertion Error: ASE output CSV for $SUBJECT is missing or empty." >&2
    exit 1
fi

SITE_COUNT=$(tail -n +2 "$OUT_CSV" | wc -l)
echo "    ASE sites counted: $SITE_COUNT"
echo ">>> Done: $SUBJECT → $OUT_CSV"
