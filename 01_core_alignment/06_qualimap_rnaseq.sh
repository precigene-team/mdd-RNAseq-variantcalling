#!/bin/bash
# [[01_core_alignment/06_qualimap_rnaseq.sh]]
# Goal: Run Qualimap RNA-seq QC on aligned BAM files.

# Steps:
#   1. Validate inputs and environment (BAM folder and GTF file existence, active mdd environment).
#   2. Detect CPU resources and set parallelism.
#   3. Find all aligned BAM files to process.
#   4. Run Qualimap RNA-seq QC in parallel for each sample BAM file.
#   5. Assert that all expected Qualimap RNA-seq reports are successfully generated.

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

BAM_DIR="$PROJECT_ROOT/_data/bam"
QC_DIR="$PROJECT_ROOT/_data/qc_rnaseq"
GTF_FILE="$PROJECT_ROOT/_data/genome/Homo_sapiens.GRCh38.115.gtf"
LOG_DIR="$PROJECT_ROOT/_data/logs/qualimap_rnaseq"

mkdir -p "$QC_DIR"
mkdir -p "$LOG_DIR"

# --- Input Assertions ---
if [ ! -d "$BAM_DIR" ]; then
    echo "Assertion Error: BAM directory not found at $BAM_DIR" >&2
    exit 1
fi
if [ ! -f "$GTF_FILE" ]; then
    echo "Assertion Error: GTF annotation file not found at $GTF_FILE" >&2
    exit 1
fi

BAM_LIST=$(find "$BAM_DIR" -maxdepth 2 -name "*_Aligned.sortedByCoord.out.bam")
SAMPLE_COUNT=$(echo "$BAM_LIST" | wc -w)

if [ "$SAMPLE_COUNT" -eq 0 ]; then
    echo "Assertion Error: No sorted coordinate BAM files found in $BAM_DIR." >&2
    exit 1
fi

# --- Resource Auto-detection ---
TOTAL_CPUS=${SLURM_CPUS_PER_TASK:-32}
PARALLEL_SAMPLES=4

echo "================================================================================"
echo "   QUALIMAP RNA-SEQ QC: Starting Qualimap rnaseq"
echo "================================================================================"
echo "Resources: $PARALLEL_SAMPLES concurrent samples."

CURRENT_SAMPLE=0
for bam in $BAM_LIST; do
    [ -f "$bam" ] || continue
    sample_name=$(basename "$bam" _Aligned.sortedByCoord.out.bam)
    CURRENT_SAMPLE=$((CURRENT_SAMPLE + 1))
    
    # Check if Qualimap report already exists
    if [ -f "$QC_DIR/$sample_name/rnaseq_qc_results.txt" ]; then
        echo "[Skip] [$CURRENT_SAMPLE/$SAMPLE_COUNT] $sample_name (report already exists)"
        continue
    fi
    
    # Simple background job management
    while [ $(jobs -r | wc -l) -ge "$PARALLEL_SAMPLES" ]; do
        sleep 5
    done
    
    echo ">>> Submitting [$CURRENT_SAMPLE/$SAMPLE_COUNT]: $sample_name"
    SAMPLE_LOG="$LOG_DIR/${sample_name}.qualimap_rnaseq.log"
    
    # Run in background
    (
        mkdir -p "$QC_DIR/$sample_name"
        qualimap rnaseq \
            -bam "$bam" \
            -gtf "$GTF_FILE" \
            -outdir "$QC_DIR/$sample_name" \
            --java-mem-size=32G
    ) > "$SAMPLE_LOG" 2>&1 &
done

echo "Waiting for all Qualimap RNA-seq QC jobs to complete..."
wait

# --- Output Assertions ---
echo "Running post-execution assertions..."
MISSING_OUTPUTS=0
for bam in $BAM_LIST; do
    [ -f "$bam" ] || continue
    sample_name=$(basename "$bam" _Aligned.sortedByCoord.out.bam)
    
    HTML_OUT="$QC_DIR/$sample_name/qualimapReport.html"
    RESULTS_OUT="$QC_DIR/$sample_name/rnaseq_qc_results.txt"
    
    if [ ! -s "$HTML_OUT" ] || [ ! -s "$RESULTS_OUT" ]; then
        echo "Assertion Error: Qualimap RNA-seq QC output missing or empty for sample $sample_name" >&2
        MISSING_OUTPUTS=$((MISSING_OUTPUTS + 1))
    fi
done

if [ "$MISSING_OUTPUTS" -gt 0 ]; then
    echo "Assertion Error: $MISSING_OUTPUTS sample(s) failed Qualimap RNA-seq QC validation." >&2
    exit 1
fi

echo "All outputs successfully validated."
echo "=== Qualimap RNA-seq QC Complete ==="
