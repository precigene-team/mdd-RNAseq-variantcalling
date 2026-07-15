#!/bin/bash
# [[01_core_alignment/04_trimmed_fastqc.sh]]
# Goal: Run FastQC on pre-alignment FastQ files — rrna_filtered for input samples,
#       trimmed_fastq for m6A samples — reflecting the actual reads fed to STAR.

# Steps:
#   1. Validate inputs and environment (sample sheet existence, active mdd environment).
#   2. Detect CPU resources and set parallelism.
#   3. Parse sample metadata (col 1, 2, 3, 9) from sample sheet.
#   4. For each sample, resolve the correct FASTQ source based on DataType:
#        - input  → _data/rrna_filtered/  (*_filtered.fastq.gz)
#        - m6A    → _data/trimmed_fastq/  (*_paired.fastq.gz)
#   5. Run FastQC on R1 and R2 in parallel.
#   6. Assert that all expected FastQC HTML and ZIP reports are generated.

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
TRIMMED_DIR="$PROJECT_ROOT/_data/trimmed_fastq"       # m6A samples
RRNA_FILTERED_DIR="$PROJECT_ROOT/_data/rrna_filtered" # input samples (rRNA-depleted)
QC_TRIM_DIR="$PROJECT_ROOT/_data/qc_trim"
LOG_DIR="$PROJECT_ROOT/_data/logs/trimmed_fastqc"

mkdir -p "$QC_TRIM_DIR"
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
echo "   TRIMMED FASTQC: Running FastQC on pre-alignment FastQ files"
echo "   (both fractions → rrna_filtered)"
echo "================================================================================"
echo "Resources: $PARALLEL_SAMPLES concurrent files."

# --- Parse sample sheet ---
# Col 1=file1_base, Col 2=file2_base, Col 3=accession, Col 9=DataType (m6A|input)
SAMPLE_ROWS=$(awk -F',' 'NR>1 {print $3 "\t" $1 "\t" $2 "\t" $9}' "$CSV_FILE" | grep -v '^$')
SAMPLE_COUNT=$(echo "$SAMPLE_ROWS" | wc -l)

if [ "$SAMPLE_COUNT" -eq 0 ]; then
    echo "Assertion Error: No sample rows parsed from $CSV_FILE" >&2
    exit 1
fi

CURRENT_SAMPLE=0
while read -r row; do
    [ -n "$row" ] || continue
    CURRENT_SAMPLE=$((CURRENT_SAMPLE + 1))

    accession=$(echo "$row" | cut -f1)
    file1_base=$(echo "$row" | cut -f2)
    file2_base=$(echo "$row" | cut -f3)
    data_type=$(echo "$row"  | cut -f4)

    # --- Resolve source directory and file suffix (both fractions are now rRNA-filtered) ---
    src_dir="$RRNA_FILTERED_DIR"
    suffix="_filtered.fastq.gz"

    r1="${src_dir}/${file1_base}${suffix}"
    r2="${src_dir}/${file2_base}${suffix}"

    if [ ! -f "$r1" ] && [ ! -f "$r2" ]; then
        echo "[Warn] [$CURRENT_SAMPLE/$SAMPLE_COUNT] $accession (type=$data_type): source files not found, skipping."
        continue
    fi

    # Simple background job management (FastQC is run per-file)
    for fastq_file in "$r1" "$r2"; do
        [ -f "$fastq_file" ] || continue

        base_name=$(basename "$fastq_file" .fastq.gz)
        html_out="$QC_TRIM_DIR/${base_name}_fastqc.html"
        zip_out="$QC_TRIM_DIR/${base_name}_fastqc.zip"

        if [ -f "$html_out" ] && [ -f "$zip_out" ]; then
            echo "[Skip] $base_name (report already exists)"
            continue
        fi

        while [ "$(jobs -r | wc -l)" -ge "$PARALLEL_SAMPLES" ]; do
            sleep 2
        done

        echo ">>> Submitting: $base_name (type=$data_type)"
        SAMPLE_LOG="$LOG_DIR/${base_name}.fastqc.log"
        fastqc -t 1 -o "$QC_TRIM_DIR" "$fastq_file" > "$SAMPLE_LOG" 2>&1 &
    done

done <<< "$SAMPLE_ROWS"

echo "Waiting for all FastQC jobs to complete..."
wait

# --- Output Assertions ---
echo "Running post-execution assertions..."
MISSING_OUTPUTS=0

while read -r row; do
    [ -n "$row" ] || continue
    accession=$(echo "$row" | cut -f1)
    file1_base=$(echo "$row" | cut -f2)
    file2_base=$(echo "$row" | cut -f3)
    data_type=$(echo "$row"  | cut -f4)

    src_dir="$RRNA_FILTERED_DIR"
    suffix="_filtered.fastq.gz"

    for base in "$file1_base" "$file2_base"; do
        src="${src_dir}/${base}${suffix}"
        [ -f "$src" ] || continue  # skip if source was missing

        html_out="$QC_TRIM_DIR/${base}${suffix%.fastq.gz}_fastqc.html"
        zip_out="$QC_TRIM_DIR/${base}${suffix%.fastq.gz}_fastqc.zip"

        if [ ! -s "$html_out" ] || [ ! -s "$zip_out" ]; then
            echo "Assertion Error: FastQC output missing or empty for $base ($accession)" >&2
            MISSING_OUTPUTS=$((MISSING_OUTPUTS + 1))
        fi
    done

done <<< "$SAMPLE_ROWS"

if [ "$MISSING_OUTPUTS" -gt 0 ]; then
    echo "Assertion Error: $MISSING_OUTPUTS file(s) failed FastQC validation." >&2
    exit 1
fi

echo "All outputs successfully validated."
echo "=== Trimmed FastQC Complete ==="
