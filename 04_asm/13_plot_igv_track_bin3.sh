#!/bin/bash
# [[04_asm/13_plot_igv_track_bin3.sh]]
# Goal: Render pyGenomeTracks and IGV visual proof of BIN3 m6A-associated SNP locus.

#SBATCH --job-name=igv_track_bin3
#SBATCH --output=_data/logs/igv_track_bin3/%x_%j.log
#SBATCH --error=_data/logs/igv_track_bin3/%x_%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00

# Steps:
#   1. Set up dynamic project paths, log directories, and active environment.
#   2. Extract bedGraph coverage for IP and Input BAMs if absent.
#   3. Render pyGenomeTracks multi-track figures (PNG and PDF).
#   4. Generate interactive igv-reports HTML alignment view centered on A>G SNP.

set -e
set -o pipefail

if [ -n "${SLURM_SUBMIT_DIR:-}" ]; then
    export PROJECT_ROOT="$SLURM_SUBMIT_DIR"
else
    export PROJECT_ROOT="${PROJECT_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
fi

cd "$PROJECT_ROOT"

eval "$(micromamba shell hook --shell bash 2>/dev/null || mamba shell hook --shell bash 2>/dev/null || conda shell.bash hook 2>/dev/null)"
micromamba activate mdd 2>/dev/null || mamba activate mdd 2>/dev/null || conda activate mdd 2>/dev/null

LOG_DIR="$PROJECT_ROOT/_data/logs/igv_track_bin3"
OUT_DIR="results/03_variant_calling/figures"
TRACK_DIR="results/03_variant_calling/tracks"
mkdir -p "$LOG_DIR" "$OUT_DIR" "$TRACK_DIR"

REGION="8:22619800-22621500"
SNP_LOCUS="8:22620952"
FASTA="_data/genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
GTF="_data/genome/Homo_sapiens.GRCh38.115.gtf"
IP_BAM="_data/bam/male_mdd_10_m6A/male_mdd_10_m6A_Aligned.sortedByCoord.out.bam"
INPUT_BAM="_data/bam/male_mdd_10_input/male_mdd_10_input_Aligned.sortedByCoord.out.bam"
IP_BG="$TRACK_DIR/male_mdd_10_m6A_bin3.bedgraph"
INPUT_BG="$TRACK_DIR/male_mdd_10_input_bin3.bedgraph"

# --- Step 1: bedGraph Extraction ---
extract_bg() {
  local bam="$1" out="$2"
  [[ -s "$out" ]] || samtools view -b "$bam" "$REGION" | bedtools genomecov -ibam stdin -bg > "$out"
}

extract_bg "$IP_BAM" "$IP_BG"
extract_bg "$INPUT_BAM" "$INPUT_BG"

# --- Step 2: pyGenomeTracks Rendering ---
echo ">>> Rendering pyGenomeTracks figure for BIN3 ($REGION)..."
INI="04_asm/13_bin3_tracks.ini"

pyGenomeTracks --tracks "$INI" --region "$REGION" --outFileName "$OUT_DIR/16_igv_track_bin3.png" --dpi 300
pyGenomeTracks --tracks "$INI" --region "$REGION" --outFileName "$OUT_DIR/16_igv_track_bin3.pdf"

# --- Step 3: Interactive IGV Read-Level HTML Report ---
SNP_BED="$TRACK_DIR/bin3_snp.bed"
echo -e "8\t22620951\t22620952\tBIN3_m6A_SNP_A_to_G" > "$SNP_BED"

if command -v create_report &>/dev/null; then
  echo ">>> Generating interactive IGV alignment report for BIN3 A>G SNP ($SNP_LOCUS)..."
  create_report "$SNP_BED" \
    --fasta "$FASTA" \
    --tracks "$IP_BAM" "$INPUT_BAM" "$GTF" \
    --flanking 200 \
    --output "$OUT_DIR/16_igv_report_bin3.html"
else
  echo ">>> Note: 'igv-reports' (create_report) not found in PATH yet."
fi

echo ">>> BIN3 track and alignment report generation complete."
