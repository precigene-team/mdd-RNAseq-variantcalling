#!/bin/bash
# [[01_core_alignment/08_counts_qc.sh]]
# Goal: Compile sample gene count files and generate alignment/assignment QC summaries using bash/awk.

# Steps:
#   1. Validate inputs and environment (count directory existence, active mdd environment).
#   2. Find all individual count files and compile the raw count matrix (raw_counts.csv).
#   3. Find all summary statistic files and compile the Assign/QC summary table (counts_qc_summary.csv).
#   4. Assert that compiled files are valid, have correct dimensions, and are non-empty.
#   5. Run QC assertions on assignment rate (>= 50%) and read depth (>= 10M assigned reads).

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

COUNT_DIR="$PROJECT_ROOT/_data/featurecounts"
OUT_MATRIX="$COUNT_DIR/raw_counts.csv"
OUT_SUMMARY="$COUNT_DIR/counts_qc_summary.csv"

echo "================================================================================"
echo "   COUNTS AGGREGATION & QC: Running featureCounts summary aggregation           "
echo "================================================================================"

# --- Input Assertions ---
if [ ! -d "$COUNT_DIR" ]; then
    echo "Assertion Error: Count directory does not exist: $COUNT_DIR" >&2
    exit 1
fi

# Find all individual count files
FILES=( $(find "$COUNT_DIR" -maxdepth 1 -name "*.txt" | grep -v "raw_counts" | sort) )
if [ ${#FILES[@]} -eq 0 ]; then
    echo "Assertion Error: No featureCounts output .txt files found in $COUNT_DIR" >&2
    exit 1
fi

# Find all summary reports
SUM_FILES=( $(find "$COUNT_DIR" -maxdepth 1 -name "*.txt.summary" | sort) )
if [ ${#SUM_FILES[@]} -eq 0 ]; then
    echo "Assertion Error: No summary files found in $COUNT_DIR" >&2
    exit 1
fi

# Assert that counts files and summary files counts match
if [ ${#FILES[@]} -ne ${#SUM_FILES[@]} ]; then
    echo "Assertion Error: Mismatch between count files (${#FILES[@]}) and summary files (${#SUM_FILES[@]})." >&2
    exit 1
fi

echo "Found ${#FILES[@]} count files to compile."

# Create a workspace-local temporary directory for intermediate columns
TMP_DIR="$COUNT_DIR/tmp_merge"
mkdir -p "$TMP_DIR"

# Clean up temp files on exit
trap 'rm -rf "$TMP_DIR"' EXIT

# 1. Compile Count Matrix
# Extract the Geneid column from the first file (skipping metadata headers starting with #)
grep -v '^#' "${FILES[0]}" | cut -f1 > "$TMP_DIR/genes.txt"

# Extract the count column (column 7) from each file, replacing the path header with the sample name
for f in "${FILES[@]}"; do
    sample=$(basename "$f" .txt)
    # The first line of data (after skipping # header) is the featureCounts header line
    grep -v '^#' "$f" | awk -v s="$sample" 'NR==1 {print s; next} {print $7}' > "$TMP_DIR/${sample}.col"
done

# Paste all columns together to form the raw count matrix
paste "$TMP_DIR/genes.txt" $(for f in "${FILES[@]}"; do echo "$TMP_DIR/$(basename "$f" .txt).col"; done) | tr '\t' ',' > "$OUT_MATRIX"
echo "Saved aggregated count matrix to: $OUT_MATRIX"


# 2. Compile featureCounts summary reports
# Extract the status categories from the first summary file
cut -f1 "${SUM_FILES[0]}" > "$TMP_DIR/summary_categories.txt"

# Extract values from each summary file, replacing path headers with the sample name
for f in "${SUM_FILES[@]}"; do
    sample=$(basename "$f" .txt.summary)
    awk -v s="$sample" 'NR==1 {print s; next} {print $2}' "$f" > "$TMP_DIR/${sample}.sum_col"
done

# Paste categories and values together to form the QC summary table
paste "$TMP_DIR/summary_categories.txt" $(for f in "${SUM_FILES[@]}"; do echo "$TMP_DIR/$(basename "$f" .txt.summary).sum_col"; done) | tr '\t' ',' > "$OUT_SUMMARY"
echo "Saved featureCounts summary QC to: $OUT_SUMMARY"


# --- Matrix Dimension Assertions ---
if [ ! -s "$OUT_MATRIX" ]; then
    echo "Assertion Error: Count matrix file $OUT_MATRIX is missing or empty." >&2
    exit 1
fi
if [ ! -s "$OUT_SUMMARY" ]; then
    echo "Assertion Error: Summary QC file $OUT_SUMMARY is missing or empty." >&2
    exit 1
fi

EXPECTED_COLS=$(( ${#FILES[@]} + 1 ))
ACTUAL_MATRIX_COLS=$(head -n 1 "$OUT_MATRIX" | tr ',' '\n' | wc -l)
ACTUAL_SUMMARY_COLS=$(head -n 1 "$OUT_SUMMARY" | tr ',' '\n' | wc -l)

if [ "$ACTUAL_MATRIX_COLS" -ne "$EXPECTED_COLS" ]; then
    echo "Assertion Error: Count matrix has $ACTUAL_MATRIX_COLS columns, expected $EXPECTED_COLS." >&2
    exit 1
fi
if [ "$ACTUAL_SUMMARY_COLS" -ne "$EXPECTED_COLS" ]; then
    echo "Assertion Error: Summary table has $ACTUAL_SUMMARY_COLS columns, expected $EXPECTED_COLS." >&2
    exit 1
fi


# 3. Quality Checks using Awk
echo ""
echo "--- Count QC Quality Checks ---"
awk -F',' '
NR==1 {
    # Store sample names
    for (i=2; i<=NF; i++) {
        samples[i] = $i
        col_count++
    }
    next
}
{
    # Store values for each category
    category = $1
    for (i=2; i<=NF; i++) {
        data[category, i] = $i
    }
    # Track status categories
    cats[NR] = category
}
END {
    low_assign_count = 0
    low_depth_count = 0
    
    # Loop through each sample column
    for (i=2; i<=(col_count+1); i++) {
        sample = samples[i]
        assigned = data["Assigned", i]
        
        # Calculate total reads for the sample by summing all categories
        total = 0
        for (c in cats) {
            category = cats[c]
            total += data[category, i]
        }
        
        if (total > 0) {
            rate = (assigned / total) * 100
        } else {
            rate = 0
        }
        
        # Filter samples under standard quality thresholds
        if (rate < 50.0) {
            low_assign[low_assign_count++] = sprintf("  %s: %.2f%%", sample, rate)
        }
        if (assigned < 10000000) {
            low_depth[low_depth_count++] = sprintf("  %s: %'\''d reads", sample, assigned)
        }
    }
    
    # Print findings
    if (low_assign_count > 0) {
        print "WARNING: The following " low_assign_count " samples have assignment rate < 50%:"
        for (j=0; j<low_assign_count; j++) {
            print low_assign[j]
        }
    } else {
        print "SUCCESS: All samples have assignment rate >= 50%."
    }
    
    if (low_depth_count > 0) {
        print "WARNING: The following " low_depth_count " samples have < 10 million assigned reads:"
        for (j=0; j<low_depth_count; j++) {
            print low_depth[j]
        }
    } else {
        print "SUCCESS: All samples have >= 10 million assigned reads."
    }
}
' "$OUT_SUMMARY"

echo "================================================================================"
echo "   Counts Aggregation & QC Complete."
echo "================================================================================"
