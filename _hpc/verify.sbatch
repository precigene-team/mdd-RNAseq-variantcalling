#!/bin/bash
# [[_hpc/verify.sbatch]]
# Goal: Slurm submission script to verify MD5 checksums of downloaded fastq files.
#SBATCH --job-name=verify_fastq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=4:00:00
#SBATCH --output=_data/logs/verify/verify_%j.log
#SBATCH --error=_data/logs/verify/verify_%j.err

# Exit on error
set -e

# --- Project Root Discovery ---
find_project_root() {
    local curr="$(pwd)"
    while [ "$curr" != "/" ] && [ -n "$curr" ]; do
        if [ -f "$curr/run" ]; then
            echo "$curr"
            return 0
        fi
        curr="$(dirname "$curr")"
    done
    echo "$(pwd)" # Fallback
}
export PROJECT_ROOT="${PROJECT_ROOT:-$(find_project_root)}"

CSV_FILE="$PROJECT_ROOT/_data/_metadata/sample_sheet.csv"
OUTPUT_DIR="$PROJECT_ROOT/_data/fastq"
VERIFY_TXT="$PROJECT_ROOT/_data/_metadata/verify.txt"
mkdir -p "$(dirname "$VERIFY_TXT")"
touch "$VERIFY_TXT"

# Function to update verification status in verify.txt
update_verify_status() {
    local file=$1
    local status=$2
    local temp_file=$(mktemp)
    if [ -f "$VERIFY_TXT" ]; then
        grep -v "^${file}," "$VERIFY_TXT" > "$temp_file" || true
    fi
    echo "${file},${status}" >> "$temp_file"
    sort "$temp_file" > "$VERIFY_TXT"
    rm -f "$temp_file"
}

# Parse all unique accessions from column 3 of the sample sheet
ACCESSIONS=$(awk -F',' 'NR>1 {print $3}' "$CSV_FILE" | grep -v '^$' | sort -u)

echo "Starting cohort MD5 verification in pure Bash..."
echo "========================================================="

CORRUPTED_COUNT=0
MISSING_COUNT=0
OK_COUNT=0
CORRUPTED_LIST=""
MISSING_LIST=""

# Count accessions
ACC_ARRAY=($ACCESSIONS)
TOTAL_ACC=${#ACC_ARRAY[@]}
CURRENT_ACC=0

for acc in $ACCESSIONS; do
    CURRENT_ACC=$((CURRENT_ACC + 1))
    
    # Check if both files for this accession are already PASS in verify.txt
    file1="${acc}_1.fastq.gz"
    file2="${acc}_2.fastq.gz"
    
    status1=$(grep "^${file1}," "$VERIFY_TXT" | cut -d',' -f2 || true)
    status2=$(grep "^${file2}," "$VERIFY_TXT" | cut -d',' -f2 || true)
    
    if [ "$status1" = "PASS" ] && [ "$status2" = "PASS" ]; then
        echo "[$CURRENT_ACC/$TOTAL_ACC] Accession $acc: already verified (PASS)."
        OK_COUNT=$((OK_COUNT + 2))
        continue
    fi
    
    # If not both are PASS, fetch expected MD5 from ENA
    echo "[$CURRENT_ACC/$TOTAL_ACC] Querying ENA MD5s for $acc..."
    md5_line=$(curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${acc}&result=read_run&fields=fastq_md5&format=tsv" | tail -n 1 | awk '{print $2}')
    
    if [ -z "$md5_line" ] || [ "$md5_line" = "fastq_md5" ]; then
        echo "  Error: Could not fetch ENA MD5 for $acc." >&2
        continue
    fi
    
    for suffix in 1 2; do
        filename="${acc}_${suffix}.fastq.gz"
        filepath="$OUTPUT_DIR/$filename"
        
        # Check individual status
        status=$(grep "^${filename}," "$VERIFY_TXT" | cut -d',' -f2 || true)
        if [ "$status" = "PASS" ]; then
            echo "  $filename: PASS (skipped, already verified)"
            OK_COUNT=$((OK_COUNT + 1))
            continue
        fi
        
        if [ ! -f "$filepath" ]; then
            echo "  $filename: MISSING"
            MISSING_COUNT=$((MISSING_COUNT + 1))
            MISSING_LIST="$MISSING_LIST $filename"
            # Update status to FAILED/MISSING
            update_verify_status "$filename" "FAILED"
            continue
        fi
        
        if [ "$suffix" = 1 ]; then
            expected_md5=$(echo "$md5_line" | cut -d';' -f1)
        else
            expected_md5=$(echo "$md5_line" | cut -d';' -f2)
        fi
        expected_md5=$(echo "$expected_md5" | tr -d '\r\n[:space:]')
        
        # Calculate actual MD5
        actual_md5=$(md5sum "$filepath" | awk '{print $1}')
        
        if [ "$actual_md5" = "$expected_md5" ]; then
            echo "  $filename: PASS"
            OK_COUNT=$((OK_COUNT + 1))
            update_verify_status "$filename" "PASS"
        else
            echo "  $filename: CORRUPTED (Expected: $expected_md5, Got: $actual_md5)"
            CORRUPTED_COUNT=$((CORRUPTED_COUNT + 1))
            CORRUPTED_LIST="$CORRUPTED_LIST $filepath"
            update_verify_status "$filename" "FAILED"
        fi
    done
done

echo "========================================================="
echo "Verification Summary:"
echo "  - Valid files: $OK_COUNT/156"
echo "  - Missing files: $MISSING_COUNT/156"
echo "  - Corrupted files: $CORRUPTED_COUNT/156"

if [ $CORRUPTED_COUNT -gt 0 ]; then
    echo ""
    echo "Corrupted files list (copy-paste to delete):"
    for f in $CORRUPTED_LIST; do
        echo "rm $f"
    done
fi

if [ $MISSING_COUNT -gt 0 ]; then
    echo ""
    echo "Missing files list:"
    for f in $MISSING_LIST; do
        echo "$f"
    done
fi
