#!/bin/bash
# [[_hpc/backup.sbatch]]
# Goal: Backup _data/fastq to /mnt/22T/vuphong/mdd_backup using rsync.

# Steps:
#   1. Validate source directory exists.
#   2. Create destination directory if missing.
#   3. rsync source -> destination (incremental, checksum-verified, verbose).
#   4. Assert destination file count matches source.

#SBATCH --job-name=MDD_Backup
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --output=_data/logs/backup/backup_log_%j.txt
#SBATCH --error=_data/logs/backup/backup_error_%j.txt

# --- Project Root Discovery ---
# Use SLURM_SUBMIT_DIR (set by SLURM to the dir where sbatch was called).
# Always submit this script from the project root: sbatch _hpc/backup.sbatch
export PROJECT_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
# ------------------------------

SRC="$PROJECT_ROOT/_data/fastq"
DEST="/mnt/22T/vuphong/mdd_backup"
METADATA_DIR="$PROJECT_ROOT/_data/_metadata"
CHECKSUM_CSV="$METADATA_DIR/checksum.csv"
SAMPLE_SHEET_CSV="$METADATA_DIR/sample_sheet.csv"

echo "========================================================="
echo "   MDD Backup: _data/fastq + metadata CSVs -> $DEST"
echo "   Starting at: $(date)"
echo "========================================================="

# --- Step 1: Validate source ---
if [ ! -d "$SRC" ]; then
    echo "Assertion Error: Source directory not found: $SRC" >&2
    exit 1
fi

# --- Step 2: Create destination ---
mkdir -p "$DEST"
if [ ! -d "$DEST" ]; then
    echo "Assertion Error: Failed to create destination directory: $DEST" >&2
    exit 1
fi

echo "Source      : $SRC"
echo "Destination : $DEST"
echo "---------------------------------------------------------"

# --- Step 3: rsync (incremental, preserve metadata, verify checksums) ---
rsync \
    --archive \
    --checksum \
    --human-readable \
    --progress \
    --stats \
    --verbose \
    "$SRC/" "$DEST/"

RSYNC_EXIT=$?
if [ "$RSYNC_EXIT" -ne 0 ]; then
    echo "Assertion Error: rsync failed with exit code $RSYNC_EXIT" >&2
    exit 1
fi

# --- Step 4: Assert file counts match ---
SRC_COUNT=$(find "$SRC" -type f | wc -l)
DEST_COUNT=$(find "$DEST" -type f | wc -l)

echo "---------------------------------------------------------"
echo "Source file count      : $SRC_COUNT"
echo "Destination file count : $DEST_COUNT"

if [ "$SRC_COUNT" -ne "$DEST_COUNT" ]; then
    echo "Assertion Error: File count mismatch! Source=$SRC_COUNT, Dest=$DEST_COUNT" >&2
    exit 1
fi

# --- Step 5: Copy metadata CSVs ---
echo "---------------------------------------------------------"
echo "Copying metadata CSVs..."
for csv_file in "$CHECKSUM_CSV" "$SAMPLE_SHEET_CSV"; do
    if [ -f "$csv_file" ]; then
        cp "$csv_file" "$DEST/"
        echo "  Copied: $(basename "$csv_file") -> $DEST/"
    else
        echo "  [WARN] Not found, skipping: $csv_file" >&2
    fi
done

echo "========================================================="
echo "   Backup COMPLETED SUCCESSFULLY"
echo "   Finished at: $(date)"
echo "========================================================="
