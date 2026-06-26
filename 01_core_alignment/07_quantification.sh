#!/bin/bash
# [[01_core_alignment/07_quantification.sh]]
# Goal: Run featureCounts on aligned BAM files.

# Steps:
#   1. Validate inputs and environment (BAM folder and GTF annotation existence, active mdd environment).
#   2. Detect CPU resources and set parallelism.
#   3. Find all aligned BAM files to process.
#   4. Run featureCounts in parallel for each BAM file to generate raw gene counts.
#   5. Assert that all expected count files and summary statistics are successfully generated.

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
COUNT_DIR="$PROJECT_ROOT/_data/featurecounts"
ANNOTATION_FILE="$PROJECT_ROOT/_data/genome/Homo_sapiens.GRCh38.115.gtf"
LOG_DIR="$PROJECT_ROOT/_data/logs/quantification"

mkdir -p "$COUNT_DIR"
mkdir -p "$LOG_DIR"

# --- Input Assertions ---
if [ ! -d "$BAM_DIR" ]; then
    echo "Assertion Error: BAM directory not found at $BAM_DIR" >&2
    exit 1
fi
if [ ! -f "$ANNOTATION_FILE" ]; then
    echo "Assertion Error: Annotation file not found at $ANNOTATION_FILE" >&2
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
PARALLEL_JOBS=4
CPUS_PER_JOB=$((TOTAL_CPUS / PARALLEL_JOBS))
[ "$CPUS_PER_JOB" -lt 1 ] && CPUS_PER_JOB=1

echo "================================================================================"
echo "   QUANTIFICATION: Running featureCounts"
echo "================================================================================"
echo "Annotation: $ANNOTATION_FILE"
echo "Output directory: $COUNT_DIR"
echo "Resources: $PARALLEL_JOBS concurrent samples, $CPUS_PER_JOB threads per sample."

CURRENT_SAMPLE=0
for bam in $BAM_LIST; do
    [ -f "$bam" ] || continue
    sample_name=$(basename "$bam" _Aligned.sortedByCoord.out.bam)
    CURRENT_SAMPLE=$((CURRENT_SAMPLE + 1))
    
    OUTPUT_FILE="$COUNT_DIR/${sample_name}.txt"
    
    # Check if counts file already exists and is non-empty
    if [ -s "$OUTPUT_FILE" ]; then
        echo "[Skip] [$CURRENT_SAMPLE/$SAMPLE_COUNT] $sample_name (counts already exist)"
        continue
    fi
    
    # Simple background job management
    while [ $(jobs -r | wc -l) -ge "$PARALLEL_JOBS" ]; do
        sleep 5
    done
    
    echo ">>> Submitting [$CURRENT_SAMPLE/$SAMPLE_COUNT]: $sample_name"
    SAMPLE_LOG="$LOG_DIR/${sample_name}.featurecounts.log"
    
    # Run in background
    featureCounts \
        -T "$CPUS_PER_JOB" \
        -p \
        -s 2 \
        -t exon \
        -g gene_id \
        -a "$ANNOTATION_FILE" \
        -o "$OUTPUT_FILE" \
        "$bam" > "$SAMPLE_LOG" 2>&1 &
        
done

echo "Waiting for all featureCounts jobs to complete..."
wait

# --- Output Assertions ---
echo "Running post-execution assertions..."
MISSING_OUTPUTS=0
for bam in $BAM_LIST; do
    [ -f "$bam" ] || continue
    sample_name=$(basename "$bam" _Aligned.sortedByCoord.out.bam)
    
    OUTPUT_FILE="$COUNT_DIR/${sample_name}.txt"
    SUMMARY_FILE="$COUNT_DIR/${sample_name}.txt.summary"
    
    if [ ! -s "$OUTPUT_FILE" ] || [ ! -s "$SUMMARY_FILE" ]; then
        echo "Assertion Error: featureCounts output missing or empty for sample $sample_name" >&2
        MISSING_OUTPUTS=$((MISSING_OUTPUTS + 1))
    fi
done

if [ "$MISSING_OUTPUTS" -gt 0 ]; then
    echo "Assertion Error: $MISSING_OUTPUTS sample(s) failed featureCounts validation." >&2
    exit 1
fi

echo "All outputs successfully validated."
echo "=== Quantification Complete ==="
