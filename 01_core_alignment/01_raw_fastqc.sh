#!/bin/bash
# [[01_core_alignment/01_raw_fastqc.sh]]
# Goal: Run FastQC on raw fastq files.

# Steps:
#   1. Validate inputs and environment (sample sheet existence, active mdd environment).
#   2. Detect CPU resources and set parallelism.
#   3. Parse sample names from sample sheet and check corresponding FASTQ files.
#   4. Run FastQC on each fastq.gz file in parallel.
#   5. Assert that all expected FastQC HTML and ZIP reports are successfully generated.

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

CSV_FILE="$PROJECT_ROOT/_data/_metadata/sample_sheet.csv"
FASTQ_DIR="$PROJECT_ROOT/_data/fastq"
# Fallback to sample_fastq if _data/fastq does not exist
[ ! -d "$FASTQ_DIR" ] && FASTQ_DIR="$PROJECT_ROOT/sample_fastq"

QC_RAW_DIR="$PROJECT_ROOT/_data/qc_raw"
LOG_DIR="$PROJECT_ROOT/_data/logs/raw_fastqc"

mkdir -p "$QC_RAW_DIR"
mkdir -p "$LOG_DIR"

# --- Input Assertions ---
if [ ! -f "$CSV_FILE" ]; then
    echo "Assertion Error: Sample sheet not found at $CSV_FILE" >&2
    exit 1
fi

# --- Resource Auto-detection ---
TOTAL_CPUS=${SLURM_CPUS_PER_TASK:-32}
PARALLEL_SAMPLES=8
[ "$PARALLEL_SAMPLES" -gt "$TOTAL_CPUS" ] && PARALLEL_SAMPLES="$TOTAL_CPUS"

echo "================================================================================"
echo "   RAW FASTQC: Running FastQC on raw FastQ files"
echo "================================================================================"
echo "Input directory: $FASTQ_DIR"
echo "Output directory: $QC_RAW_DIR"
echo "Resources: $PARALLEL_SAMPLES concurrent samples."

# Parse unique sample names from sample_sheet.csv (skip header, get File1 and File2 columns)
SAMPLES=$(awk -F',' 'NR>1 {print $1; print $2}' "$CSV_FILE" | grep -v '^$' | sort -u)
SAMPLE_COUNT=$(echo "$SAMPLES" | wc -w)

# Assert that we parsed at least one sample
if [ "$SAMPLE_COUNT" -eq 0 ]; then
    echo "Assertion Error: No sample names could be parsed from $CSV_FILE" >&2
    exit 1
fi

# Assert that at least one FASTQ file was found to process
FOUND_COUNT=0
for sample in $SAMPLES; do
    [ -n "$sample" ] || continue
    FASTQ_FILE="$FASTQ_DIR/${sample}.fastq.gz"
    [ -f "$FASTQ_FILE" ] && FOUND_COUNT=$((FOUND_COUNT + 1))
done

if [ "$FOUND_COUNT" -eq 0 ]; then
    echo "Assertion Error: No raw FASTQ files found in $FASTQ_DIR matching sample sheet entries." >&2
    exit 1
fi

CURRENT_SAMPLE=0
for sample in $SAMPLES; do
    [ -n "$sample" ] || continue
    CURRENT_SAMPLE=$((CURRENT_SAMPLE + 1))
    
    FASTQ_FILE="$FASTQ_DIR/${sample}.fastq.gz"
    
    # Check if fastq file exists
    if [ ! -f "$FASTQ_FILE" ]; then
        echo "[Warn] File not found: $FASTQ_FILE. Skipping."
        continue
    fi
    
    # Check if FastQC report already exists
    if [ -f "$QC_RAW_DIR/${sample}_fastqc.html" ] && [ -f "$QC_RAW_DIR/${sample}_fastqc.zip" ]; then
        echo "[Skip] [$CURRENT_SAMPLE/$SAMPLE_COUNT] $sample (report already exists)"
        continue
    fi
    
    # Simple background job management
    while [ $(jobs -r | wc -l) -ge "$PARALLEL_SAMPLES" ]; do
        sleep 2
    done
    
    echo ">>> Submitting [$CURRENT_SAMPLE/$SAMPLE_COUNT]: $sample"
    SAMPLE_LOG="$LOG_DIR/${sample}.fastqc.log"
    
    # Run in background
    fastqc -t 1 -o "$QC_RAW_DIR" "$FASTQ_FILE" > "$SAMPLE_LOG" 2>&1 &
done

echo "Waiting for all raw FastQC jobs to complete..."
wait

# --- Output Assertions ---
echo "Running post-execution assertions..."
MISSING_OUTPUTS=0
for sample in $SAMPLES; do
    [ -n "$sample" ] || continue
    FASTQ_FILE="$FASTQ_DIR/${sample}.fastq.gz"
    [ -f "$FASTQ_FILE" ] || continue # If input was missing, output won't exist
    
    HTML_OUT="$QC_RAW_DIR/${sample}_fastqc.html"
    ZIP_OUT="$QC_RAW_DIR/${sample}_fastqc.zip"
    
    if [ ! -s "$HTML_OUT" ] || [ ! -s "$ZIP_OUT" ]; then
        echo "Assertion Error: FastQC output missing or empty for sample $sample" >&2
        MISSING_OUTPUTS=$((MISSING_OUTPUTS + 1))
    fi
done

if [ "$MISSING_OUTPUTS" -gt 0 ]; then
    echo "Assertion Error: $MISSING_OUTPUTS sample(s) failed FastQC validation." >&2
    exit 1
fi

echo "All outputs successfully validated."
echo "=== Raw FastQC Complete ==="
