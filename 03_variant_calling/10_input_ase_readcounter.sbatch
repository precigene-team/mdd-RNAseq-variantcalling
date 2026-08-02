#!/bin/bash
# [[03_variant_calling/10_input_ase_readcounter.sbatch]]
# Goal: Run GATK ASEReadCounter on INPUT BAMs at each subject's heterozygous SNP positions,
#       to obtain per-subject INPUT allelic counts for ASM group-comparison (Input-as-null).

# Why this exists:
#   The per-peak beta-binomial ASM model (docs/asm_analysis_methodology.md, §4) needs the
#   Input allelic fraction as the null/baseline per subject, not a fixed 0.5, to control
#   mapping/reference bias (Castel et al. 2015). 06_ase_readcounter.sbatch counted the IP
#   fraction only; this job counts the Input fraction at the SAME heterozygous SNP positions.

# How to run (on HPC, from the pipeline root):
#   sbatch 03_variant_calling/10_input_ase_readcounter.sbatch
#
# After completion, sync the new input counts to the local workstation:
#   rsync -avz --ignore-existing \
#       <hpc>:/mnt/22T/phongdinh/RNAseq_pipeline_depression/_data/ase/*_input_ase_counts.csv \
#       <local>:/home/mashxp/Documents/Bioinformatics/mdd_RNAseq_pipeline/_data/ase/

# Steps:
#   1. Validate inputs (Input BAM, existing het-SNP VCF, reference FASTA).
#   2. Run gatk ASEReadCounter on the Input BAM at the het-SNP positions.
#   3. Assert output CSV is non-empty.

#SBATCH --job-name=input_ase_readcounter
#SBATCH --output=_data/logs/ase/%x_%A_%a.log
#SBATCH --error=_data/logs/ase/%x_%A_%a.err
#SBATCH --array=0-38%5
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
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

# --- Paths ---
REF_FA="$PROJECT_ROOT/_data/genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
BAM_DIR="$PROJECT_ROOT/_data/variant/split_n_cigar"
ASE_DIR="$PROJECT_ROOT/_data/ase"
LOG_DIR="$PROJECT_ROOT/_data/logs/ase"
SAMPLE_SHEET="$PROJECT_ROOT/_data/_metadata/sample_sheet.csv"

mkdir -p "$ASE_DIR" "$LOG_DIR"

# --- Resolve subject from SLURM array index ---
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
INPUT_BAM="$BAM_DIR/${SUBJECT}_input.split.bam"
HET_VCF="$ASE_DIR/${SUBJECT}_het_snps.vcf.gz"
OUT_CSV="$ASE_DIR/${SUBJECT}_input_ase_counts.csv"

# --- Input Assertions ---
if [ ! -f "$INPUT_BAM" ]; then
    echo "Assertion Error: Input BAM not found: $INPUT_BAM" >&2
    exit 1
fi
if [ ! -f "$HET_VCF" ]; then
    echo "Assertion Error: Het-SNP VCF not found (run 06_ase_readcounter.sbatch first): $HET_VCF" >&2
    exit 1
fi

# Skip if already done
if [ -s "$OUT_CSV" ]; then
    echo "[SKIP] $SUBJECT: Input ASE counts already exist."
    exit 0
fi

# --- Step 1: Run ASEReadCounter on the INPUT BAM at the het-SNP positions ---
echo ">>> Running GATK ASEReadCounter on INPUT BAM..."
gatk ASEReadCounter \
    -R "$REF_FA" \
    -I "$INPUT_BAM" \
    -V "$HET_VCF" \
    -O "$OUT_CSV" \
    --min-mapping-quality 20 \
    --min-base-quality 20 \
    --disable-read-filter NotDuplicateReadFilter

# --- Output Assertion ---
if [ ! -s "$OUT_CSV" ]; then
    echo "Assertion Error: Input ASE output CSV for $SUBJECT is missing or empty." >&2
    exit 1
fi

SITE_COUNT=$(tail -n +2 "$OUT_CSV" | wc -l)
echo "    Input ASE sites counted: $SITE_COUNT"
echo ">>> Done: $SUBJECT → $OUT_CSV"
