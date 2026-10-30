#!/bin/bash
# [[01_core_alignment/10_multiqc.sh]]
# Goal: Run MultiQC to centrally aggregate all quality reports.

# Steps:
#   1. Validate inputs and environment (active mdd environment, config file existence).
#   2. Run MultiQC to aggregate all raw/trimmed QC, alignment, and quantification logs.
#   3. Assert that the multiqc HTML report is successfully generated and non-empty.

# Exit on error
set -e
set -o pipefail

# --- Project Root Discovery ---
export PROJECT_ROOT="${PROJECT_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
# ------------------------------

# --- Environment Assertion ---
if [ -z "$CONDA_PREFIX" ]; then
    echo "Assertion Error: Mamba environment 'mdd' is not active. Please activate it first." >&2
    exit 1
fi

QC_RAW_DIR="$PROJECT_ROOT/_data/qc_raw"
QC_TRIM_DIR="$PROJECT_ROOT/_data/qc_trim"
QUALIMAP_BAMQC_DIR="$PROJECT_ROOT/_data/qc_bam"
QUALIMAP_RNASEQ_DIR="$PROJECT_ROOT/_data/qc_rnaseq"
STAR_DIR="$PROJECT_ROOT/_data/bam"
FEATURECOUNTS_DIR="$PROJECT_ROOT/_data/featurecounts"

MULTIQC_OUT_DIR="$PROJECT_ROOT/_data/multiqc"
CONFIG_FILE="$PROJECT_ROOT/multiqc_config.yaml"

echo "================================================================================"
echo "   MULTIQC AGGREGATION: Running MultiQC"
echo "================================================================================"
echo "Config: $CONFIG_FILE"
echo "Output: $MULTIQC_OUT_DIR"

# --- Input Assertions ---
if [ ! -f "$CONFIG_FILE" ]; then
    echo "Assertion Error: MultiQC configuration not found at $CONFIG_FILE" >&2
    exit 1
fi

mkdir -p "$MULTIQC_OUT_DIR"

# Run MultiQC referencing the central config
multiqc \
    "$QC_RAW_DIR" \
    "$QC_TRIM_DIR" \
    "$QUALIMAP_BAMQC_DIR" \
    "$QUALIMAP_RNASEQ_DIR" \
    "$STAR_DIR" \
    "$FEATURECOUNTS_DIR" \
    -o "$MULTIQC_OUT_DIR" \
    -f \
    -n rna_seq_pipeline_summary \
    -c "$CONFIG_FILE"

# --- Output Assertions ---
if [ ! -s "$MULTIQC_OUT_DIR/rna_seq_pipeline_summary.html" ]; then
    echo "Assertion Error: MultiQC summary html was not generated or is empty." >&2
    exit 1
fi

echo "================================================================================"
echo "   MultiQC Complete. Summary html generated at:"
echo "   $MULTIQC_OUT_DIR/rna_seq_pipeline_summary.html"
echo "================================================================================"
