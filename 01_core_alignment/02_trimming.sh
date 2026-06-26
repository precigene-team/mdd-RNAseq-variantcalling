#!/bin/bash
# [[01_core_alignment/02_trimming.sh]]
# Goal: Run fastp adapter and quality trimming on paired-end FastQ files.

# Steps:
#   1. Validate inputs and environment (sample sheet existence, active mdd environment).
#   2. Detect CPU resources and set parallelism.
#   3. Parse sample metadata from sample sheet and match input FASTQ files.
#   4. Run fastp in parallel for paired-end quality and adapter trimming.
#   5. Assert that all expected trimmed paired FASTQ files are successfully generated.

# Exit on error
set -e
set -o pipefail

# --- Project Root Discovery ---
export PROJECT_ROOT="${PROJECT_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
# ------------------------------

CSV_FILE="$PROJECT_ROOT/_data/_metadata/sample_sheet.csv"
FASTQ_DIR="$PROJECT_ROOT/_data/fastq"
[ ! -d "$FASTQ_DIR" ] && FASTQ_DIR="$PROJECT_ROOT/sample_fastq"

OUTDIR="$PROJECT_ROOT/_data/trimmed_fastq"
QC_OUT_DIR="$PROJECT_ROOT/_data/qc_trim_fastp"
LOG_DIR="$PROJECT_ROOT/_data/logs/trimming"

mkdir -p "$OUTDIR"
mkdir -p "$QC_OUT_DIR"
mkdir -p "$LOG_DIR"

# --- Input Assertions ---
if [ ! -f "$CSV_FILE" ]; then
    echo "Assertion Error: Sample sheet not found at $CSV_FILE" >&2
    exit 1
fi

# --- Resource Auto-detection ---
TOTAL_CPUS=${SLURM_CPUS_PER_TASK:-32}
PARALLEL_SAMPLES=4
CPUS_PER_SAMPLE=$((TOTAL_CPUS / PARALLEL_SAMPLES))
[ "$CPUS_PER_SAMPLE" -lt 1 ] && CPUS_PER_SAMPLE=1

# --- fastp Check ---
if [ -z "$CONDA_PREFIX" ]; then
    echo "Error: Mamba environment 'mdd' is not active. Please activate it first." >&2
    exit 1
fi

if ! command -v fastp &> /dev/null; then
    echo "Error: fastp executable not found in $CONDA_PREFIX. Please install it first." >&2
    exit 1
fi

echo "================================================================================"
echo "   TRIMMING: Running fastp PE"
echo "================================================================================"
echo "Output directory: $OUTDIR"
echo "QC Reports directory: $QC_OUT_DIR"
echo "Resources: $PARALLEL_SAMPLES concurrent samples, $CPUS_PER_SAMPLE threads per sample."

# Parse sample sheet
# Accession is column 3, File1 is column 1, File2 is column 2
# Format: Accession File1 File2
SAMPLE_ROWS=$(awk -F',' 'NR>1 {print $3 "\t" $1 "\t" $2}' "$CSV_FILE" | grep -v '^$')
SAMPLE_COUNT=$(echo "$SAMPLE_ROWS" | wc -l)

# Assert that we parsed at least one sample row
if [ "$SAMPLE_COUNT" -eq 0 ]; then
    echo "Assertion Error: No sample rows parsed from $CSV_FILE" >&2
    exit 1
fi

# Assert that at least one paired fastq files exist
FOUND_INPUTS=0
while read -r row; do
    [ -n "$row" ] || continue
    file1=$(echo "$row" | cut -f2).fastq.gz
    file2=$(echo "$row" | cut -f3).fastq.gz
    if [ -f "$FASTQ_DIR/$file1" ] && [ -f "$FASTQ_DIR/$file2" ]; then
        FOUND_INPUTS=$((FOUND_INPUTS + 1))
    fi
done <<< "$SAMPLE_ROWS"

if [ "$FOUND_INPUTS" -eq 0 ]; then
    echo "Assertion Error: No matching paired FastQ files found in $FASTQ_DIR." >&2
    exit 1
fi

CURRENT_SAMPLE=0
while read -r row; do
    [ -n "$row" ] || continue
    CURRENT_SAMPLE=$((CURRENT_SAMPLE + 1))
    
    accession=$(echo "$row" | cut -f1)
    file1_base=$(echo "$row" | cut -f2)
    file2_base=$(echo "$row" | cut -f3)
    
    r1="$FASTQ_DIR/${file1_base}.fastq.gz"
    r2="$FASTQ_DIR/${file2_base}.fastq.gz"
    
    if [ ! -f "$r1" ] || [ ! -f "$r2" ]; then
        echo "[Warn] Input file missing for $accession: ${file1_base}.fastq.gz or ${file2_base}.fastq.gz. Skipping."
        continue
    fi
    
    # Name outputs with single extension to match 03_trimmed_fastqc.sh and 04_star_align.sh
    r1_paired="$OUTDIR/${file1_base}_paired.fastq.gz"
    r1_unpaired="$OUTDIR/${file1_base}_unpaired.fastq.gz"
    r2_paired="$OUTDIR/${file2_base}_paired.fastq.gz"
    r2_unpaired="$OUTDIR/${file2_base}_unpaired.fastq.gz"
    
    # fastp HTML and JSON reports
    report_html="$QC_OUT_DIR/${accession}_fastp.html"
    report_json="$QC_OUT_DIR/${accession}_fastp.json"
    
    # Check if trimmed outputs already exist
    if [ -f "$r1_paired" ] && [ -f "$r2_paired" ]; then
        echo "[Skip] [$CURRENT_SAMPLE/$SAMPLE_COUNT] $accession (trimmed files already exist)"
        continue
    fi
    
    # Simple background job management
    while [ $(jobs -r | wc -l) -ge "$PARALLEL_SAMPLES" ]; do
        sleep 5
    done
    
    echo ">>> Submitting [$CURRENT_SAMPLE/$SAMPLE_COUNT]: $accession"
    SAMPLE_LOG="$LOG_DIR/${accession}.fastp.log"
    
    # Run in background
    # Parameters matching Trimmomatic's sliding window behavior (4:20) and min length (36):
    # --detect_adapter_for_pe handles adapter auto-detection
    # --cut_front and --cut_tail act as sliding windows
    fastp \
        -w "$CPUS_PER_SAMPLE" \
        -i "$r1" -I "$r2" \
        -o "$r1_paired" -O "$r2_paired" \
        --unpaired1 "$r1_unpaired" --unpaired2 "$r2_unpaired" \
        -h "$report_html" -j "$report_json" \
        --detect_adapter_for_pe \
        --cut_front --cut_tail \
        --cut_window_size 4 --cut_mean_quality 20 \
        --length_required 36 > "$SAMPLE_LOG" 2>&1 &
        
done <<< "$SAMPLE_ROWS"

echo "Waiting for all fastp jobs to complete..."
wait

# --- Output Assertions ---
echo "Running post-execution assertions..."
MISSING_OUTPUTS=0
while read -r row; do
    [ -n "$row" ] || continue
    accession=$(echo "$row" | cut -f1)
    file1_base=$(echo "$row" | cut -f2)
    file2_base=$(echo "$row" | cut -f3)
    
    r1="$FASTQ_DIR/${file1_base}.fastq.gz"
    r2="$FASTQ_DIR/${file2_base}.fastq.gz"
    [ -f "$r1" ] && [ -f "$r2" ] || continue # If inputs were missing, outputs won't exist
    
    r1_paired="$OUTDIR/${file1_base}_paired.fastq.gz"
    r2_paired="$OUTDIR/${file2_base}_paired.fastq.gz"
    
    if [ ! -s "$r1_paired" ] || [ ! -s "$r2_paired" ]; then
        echo "Assertion Error: Trimmed output missing or empty for accession $accession (R1/R2 paired)" >&2
        MISSING_OUTPUTS=$((MISSING_OUTPUTS + 1))
    fi
done <<< "$SAMPLE_ROWS"

if [ "$MISSING_OUTPUTS" -gt 0 ]; then
    echo "Assertion Error: $MISSING_OUTPUTS accession(s) failed fastp output validation." >&2
    exit 1
fi

echo "All outputs successfully validated."
echo "=== Trimming Complete ==="
