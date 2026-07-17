#!/bin/bash
# [[utils/ame_motif_enrichment.sh]]
# Goal: Extract TRESS consensus peak sequences and run AME to validate GGACH motif enrichment.

# Steps:
#   1. Validate environment and input files (consensus BED, reference FASTA).
#   2. Extract peak sequences using bedtools getfasta.
#   3. Build a 1st-order Markov background model from peak sequences.
#   4. Run AME with the GGACH motif matrix (--text mode; HTML template broken in v5.5.9).
#   5. Generate enrichment bar chart from the summary.
#   6. Clean up intermediate files.

set -e
set -o pipefail

export PROJECT_ROOT="${PROJECT_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"

# --- Environment Assertion ---
if [ -z "$CONDA_PREFIX" ]; then
    echo "Assertion Error: Mamba environment 'mdd' is not active." >&2
    exit 1
fi

# --- Paths ---
PEAKS_BED="$PROJECT_ROOT/_data/peaks/consensus_peaks.bed"
GENOME_FA="$PROJECT_ROOT/_data/genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
MOTIF_FILE="$PROJECT_ROOT/utils/ggach_motive.meme"
TMP_DIR="$PROJECT_ROOT/_data/tmp/ame_motif"
FIGURES_DIR="$PROJECT_ROOT/results/utils/figures"
TABLES_DIR="$PROJECT_ROOT/results/utils/tables"
PEAKS_FA="$TMP_DIR/consensus_peaks.fa"
BG_MODEL="$TMP_DIR/background.model"
LOG_DIR="$PROJECT_ROOT/_data/logs/ame_motif"

mkdir -p "$FIGURES_DIR" "$TABLES_DIR" "$LOG_DIR" "$TMP_DIR"

# --- Input Assertions ---
for f in "$PEAKS_BED" "$GENOME_FA" "$MOTIF_FILE"; do
    if [ ! -f "$f" ]; then
        echo "Assertion Error: Required input not found: $f" >&2
        exit 1
    fi
done

echo "================================================================================"
echo "   AME MOTIF ENRICHMENT: GGACH Validation on TRESS Consensus Peaks"
echo "================================================================================"
echo "Peaks BED:  $PEAKS_BED"
echo "Genome:     $GENOME_FA"
echo "Motif:      $MOTIF_FILE"
echo "Table:      $TABLES_DIR/ame_enrichment_summary.tsv"
echo "Figure:     $FIGURES_DIR/ame_enrichment.png"

# 1. Extract peak sequences
echo ">>> Step 1: Extracting peak sequences with bedtools getfasta..."
bedtools getfasta \
    -fi "$GENOME_FA" \
    -bed "$PEAKS_BED" \
    -fo "$PEAKS_FA" \
    -name 2> "$LOG_DIR/getfasta.log"

PEAK_COUNT=$(grep -c "^>" "$PEAKS_FA")
echo "    Extracted $PEAK_COUNT sequences."
if [ "$PEAK_COUNT" -eq 0 ]; then
    echo "Assertion Error: No sequences were extracted from peaks FASTA." >&2
    exit 1
fi

# 2. Build background Markov model
echo ">>> Step 2: Building 1st-order background Markov model..."
fasta-get-markov -m 1 "$PEAKS_FA" "$BG_MODEL" 2>> "$LOG_DIR/getfasta.log"

# 3. Run AME
echo ">>> Step 3: Running AME motif enrichment..."
ame \
    --text \
    --control --shuffle-- \
    --scoring avg \
    --method fisher \
    --bfile "$BG_MODEL" \
    "$PEAKS_FA" \
    "$MOTIF_FILE" \
    > "$TABLES_DIR/ame_enrichment_summary.tsv" \
    2>"$LOG_DIR/ame.log"

# --- Output Assertions ---
if [ ! -s "$TABLES_DIR/ame_enrichment_summary.tsv" ]; then
    echo "Assertion Error: AME enrichment summary was not generated or is empty." >&2
    exit 1
fi

# 4. Generate enrichment figure
echo ">>> Step 4: Generating enrichment figure..."
python3 "$PROJECT_ROOT/utils/ame_plot_enrichment.py" \
    "$TABLES_DIR/ame_enrichment_summary.tsv" \
    "$FIGURES_DIR/ame_enrichment.png" \
    2>> "$LOG_DIR/ame.log"

# --- Clean up intermediate files ---
echo ">>> Step 5: Cleaning up intermediate files..."
rm -rf "$TMP_DIR"

echo "================================================================================"
echo "   AME COMPLETE."
echo "   Summary:  $TABLES_DIR/ame_enrichment_summary.tsv"
echo "   Figure:   $FIGURES_DIR/ame_enrichment.png"
echo "================================================================================"
