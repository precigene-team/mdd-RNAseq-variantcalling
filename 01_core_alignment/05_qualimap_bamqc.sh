#!/bin/bash
# [[01_core_alignment/05_qualimap_bamqc.sh]]
# Goal: Run Qualimap BAM QC on aligned BAM files.

# Steps:
#   1. Validate inputs and environment (BAM folder existence, active mdd environment).
#   2. Detect CPU resources and set parallelism.
#   3. Find all aligned BAM files to process.
#   4. Run Qualimap BAM QC in parallel for each sample BAM file.
#   5. Assert that all expected Qualimap reports are successfully generated.

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
QC_DIR="$PROJECT_ROOT/_data/qc_bam"
LOG_DIR="$PROJECT_ROOT/_data/logs/qualimap_bamqc"

mkdir -p "$QC_DIR"
mkdir -p "$LOG_DIR"

# --- Input Assertions ---
if [ ! -d "$BAM_DIR" ]; then
    echo "Assertion Error: BAM directory not found at $BAM_DIR" >&2
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
CPUS_PER_SAMPLE=$((TOTAL_CPUS / PARALLEL_SAMPLES))
[ "$CPUS_PER_SAMPLE" -lt 1 ] && CPUS_PER_SAMPLE=1

echo "================================================================================"
echo "   QUALIMAP BAM QC: Starting Qualimap bamqc"
echo "================================================================================"
echo "Resources: $PARALLEL_SAMPLES concurrent samples, $CPUS_PER_SAMPLE threads per sample."

CURRENT_SAMPLE=0
for bam in $BAM_LIST; do
    [ -f "$bam" ] || continue
    sample_name=$(basename "$bam" _Aligned.sortedByCoord.out.bam)
    CURRENT_SAMPLE=$((CURRENT_SAMPLE + 1))
    
    # Check if Qualimap report already exists
    if [ -f "$QC_DIR/$sample_name/genome_results.txt" ] && [ -f "$QC_DIR/$sample_name/qualimapReport.html" ]; then
        echo "[Skip] [$CURRENT_SAMPLE/$SAMPLE_COUNT] $sample_name (report already exists)"
        continue
    fi
    
    # Simple background job management
    while [ $(jobs -r | wc -l) -ge "$PARALLEL_SAMPLES" ]; do
        sleep 5
    done
    
    echo ">>> Submitting [$CURRENT_SAMPLE/$SAMPLE_COUNT]: $sample_name"
    SAMPLE_LOG="$LOG_DIR/${sample_name}.qualimap_bamqc.log"
    
    # Run in background
    (
        mkdir -p "$QC_DIR/$sample_name"
        qualimap bamqc \
            -bam "$bam" \
            -outdir "$QC_DIR/$sample_name" \
            -nt "$CPUS_PER_SAMPLE" \
            --java-mem-size=16G
    ) > "$SAMPLE_LOG" 2>&1 &
done

echo "Waiting for all Qualimap BAM QC jobs to complete..."
wait

# --- Output Assertions ---
echo "Running post-execution assertions..."
MISSING_OUTPUTS=0
for bam in $BAM_LIST; do
    [ -f "$bam" ] || continue
    sample_name=$(basename "$bam" _Aligned.sortedByCoord.out.bam)
    
    HTML_OUT="$QC_DIR/$sample_name/qualimapReport.html"
    RESULTS_OUT="$QC_DIR/$sample_name/genome_results.txt"
    
    if [ ! -s "$HTML_OUT" ] || [ ! -s "$RESULTS_OUT" ]; then
        echo "Assertion Error: Qualimap BAM QC output missing or empty for sample $sample_name" >&2
        MISSING_OUTPUTS=$((MISSING_OUTPUTS + 1))
    fi
done

if [ "$MISSING_OUTPUTS" -gt 0 ]; then
    echo "Assertion Error: $MISSING_OUTPUTS sample(s) failed Qualimap BAM QC validation." >&2
    exit 1
fi

echo "All outputs successfully validated."
echo "=== Qualimap BAM QC Complete ==="
