#!/bin/bash
# [[_hpc/generate_checksum_csv.sbatch]]
# Goal: Fetch MD5 checksums from ENA API for all accessions and produce
#       _data/_metadata/checksum.csv with columns: accession,filename,md5,url

# Steps:
#   1. Parse all unique accessions from sample_sheet.csv (col 3).
#   2. For each accession, query the ENA API for fastq_md5 (semicolon-delimited: R1;R2).
#   3. Extract the matching URL from all_ena_urls.txt.
#   4. Emit one row per file (R1 and R2) into checksum.csv.
#   5. Assert row count == 156 (78 accessions x 2 files each).

#SBATCH --job-name=MDD_Checksum_CSV
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=2:00:00
#SBATCH --output=_data/logs/checksum/checksum_csv_%j.log
#SBATCH --error=_data/logs/checksum/checksum_csv_%j.err

set -e
set -o pipefail

# --- Project Root Discovery ---
# Use SLURM_SUBMIT_DIR (set by SLURM to the dir where sbatch was called).
# Always submit this script from the project root: sbatch _hpc/generate_checksum_csv.sbatch
export PROJECT_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
# ------------------------------

CSV_FILE="$PROJECT_ROOT/_data/_metadata/sample_sheet.csv"
URL_FILE="$PROJECT_ROOT/_data/_metadata/all_ena_urls.txt"
OUT_FILE="$PROJECT_ROOT/_data/_metadata/checksum.csv"

# --- Input Assertions ---
if [ ! -f "$CSV_FILE" ]; then
    echo "Assertion Error: sample_sheet.csv not found at $CSV_FILE" >&2
    exit 1
fi
if [ ! -f "$URL_FILE" ]; then
    echo "Assertion Error: all_ena_urls.txt not found at $URL_FILE" >&2
    exit 1
fi

echo "========================================================="
echo "   Generating checksum.csv from ENA API"
echo "   Starting at: $(date)"
echo "========================================================="

# Write header
echo "accession,filename,md5,url" > "$OUT_FILE"

ACCESSIONS=$(awk -F',' 'NR>1 {print $3}' "$CSV_FILE" | grep -v '^$' | sort -u)
ACC_ARRAY=($ACCESSIONS)
TOTAL=${#ACC_ARRAY[@]}
CURRENT=0
ROW_COUNT=0
FAILED_LIST=""

for acc in $ACCESSIONS; do
    CURRENT=$((CURRENT + 1))
    echo "[$CURRENT/$TOTAL] Fetching MD5s for $acc..."

    # Query ENA API — same pattern as verify.sbatch
    md5_line=$(curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${acc}&result=read_run&fields=fastq_md5&format=tsv" \
                | tail -n 1 | awk '{print $2}')

    # --- Assertion: non-empty and not the header ---
    if [ -z "$md5_line" ] || [ "$md5_line" = "fastq_md5" ]; then
        echo "  [WARN] Could not fetch MD5 for $acc. Recording UNKNOWN." >&2
        FAILED_LIST="$FAILED_LIST $acc"
        md5_r1="UNKNOWN"
        md5_r2="UNKNOWN"
    else
        md5_r1=$(echo "$md5_line" | cut -d';' -f1 | tr -d '\r\n[:space:]')
        md5_r2=$(echo "$md5_line" | cut -d';' -f2 | tr -d '\r\n[:space:]')
    fi

    for suffix in 1 2; do
        filename="${acc}_${suffix}.fastq.gz"

        # Look up the URL from all_ena_urls.txt
        url=$(grep "/${filename}$" "$URL_FILE" || true)
        if [ -z "$url" ]; then
            echo "  [WARN] URL not found for $filename in all_ena_urls.txt" >&2
            url="UNKNOWN"
        fi

        if [ "$suffix" = "1" ]; then
            md5="$md5_r1"
        else
            md5="$md5_r2"
        fi

        echo "${acc},${filename},${md5},${url}" >> "$OUT_FILE"
        ROW_COUNT=$((ROW_COUNT + 1))
    done

    # Rate-limit: be polite to ENA API
    sleep 0.5
done

echo "========================================================="
echo "Output written to: $OUT_FILE"
echo "Total rows written: $ROW_COUNT (expected 156)"

# --- Output Assertion ---
if [ "$ROW_COUNT" -ne 156 ]; then
    echo "Assertion Error: Expected 156 rows, got $ROW_COUNT." >&2
    exit 1
fi

if [ -n "$FAILED_LIST" ]; then
    echo ""
    echo "[WARN] The following accessions returned UNKNOWN MD5s (check network or ENA availability):"
    for f in $FAILED_LIST; do
        echo "  $f"
    done
fi

echo "   Finished at: $(date)"
echo "========================================================="
