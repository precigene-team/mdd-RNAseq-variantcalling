#!/bin/bash
# [[_hpc/rename.sbatch]]
# Goal: Slurm submission script to rename verified raw fastq files to sample sheet names.
#SBATCH --job-name=rename_fastq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:30:00
#SBATCH --output=_data/logs/rename/rename_%j.log
#SBATCH --error=_data/logs/rename/rename_%j.err

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

# Parse all unique accessions from column 3 of the sample sheet
ACCESSIONS=$(awk -F',' 'NR>1 {print $3}' "$CSV_FILE" | grep -v '^$' | sort -u)

echo "Starting cohort renaming..."
echo "========================================================="

RENAME_COUNT=0
for acc in $ACCESSIONS; do
    for suffix in 1 2; do
        filename="${acc}_${suffix}.fastq.gz"
        filepath="$OUTPUT_DIR/$filename"
        
        # If the raw accession file exists, rename it
        if [ -f "$filepath" ]; then
            if [ "$suffix" = "1" ]; then
                target_name=$(awk -F',' -v acc="$acc" '$3 == acc {print $1}' "$CSV_FILE")
            else
                target_name=$(awk -F',' -v acc="$acc" '$3 == acc {print $2}' "$CSV_FILE")
            fi
            target_name=$(echo "$target_name" | tr -d '\r\n[:space:]')
            
            if [ -n "$target_name" ]; then
                echo "  Renaming $filename -> ${target_name}.fastq.gz"
                mv "$filepath" "$OUTPUT_DIR/${target_name}.fastq.gz"
                RENAME_COUNT=$((RENAME_COUNT + 1))
            else
                echo "  Warning: No matching sample sheet name found for $filename. Keeping raw name."
            fi
        fi
    done
done

echo "========================================================="
echo "Renaming complete: Renamed $RENAME_COUNT files."
