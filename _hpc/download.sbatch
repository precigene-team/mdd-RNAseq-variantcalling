#!/bin/bash
# [[_hpc/download.sbatch]]
# Goal: Slurm submission script to download raw FASTQ files from ENA via Aspera (ascli).
#SBATCH --job-name=download_fastq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=48:00:00
#SBATCH --output=_data/logs/download_%j.log
#SBATCH --error=_data/logs/download_%j.err

# Steps:
#   1. Discover project root and set paths
#   2. Activate mamba env (mdd) and configure ascli / gem paths
#   3. Register SIGTERM trap so scancel cleanly kills ascp (which ignores SIGTERM)
#   4. Read ENA URL list; skip files already present in the output directory
#   5. Download remaining files via: ascli -Pera server download ... --to-folder=...

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

INPUT_FILE="$PROJECT_ROOT/_data/_metadata/all_ena_urls.txt"
OUTPUT_DIR="$PROJECT_ROOT/_data/fastq"

# 1. Initialize Mamba/Micromamba
eval "$(mamba shell hook --shell bash 2>/dev/null || micromamba shell hook --shell bash 2>/dev/null)"
mamba activate mdd 2>/dev/null || micromamba activate mdd 2>/dev/null
export GEM_PATH="$CONDA_PREFIX/share/rubygems"
export ASCLI_PROGRESS_BAR=yes

# Cleanup trap: ascp ignores SIGTERM by default, so we must explicitly kill
# all child processes when Slurm issues scancel. This ensures nothing lingers.
cleanup() {
    echo ">>> [SIGTERM] Job cancelled — terminating all transfer processes..."
    # Kill all children of this script (covers unbuffer, ascli, ascp)
    pkill -P $$ 2>/dev/null || true
    # Remove temp files if they exist
    [ -n "${TEMP_URLS:-}" ] && rm -f "$TEMP_URLS"
    exit 1
}
trap cleanup SIGTERM SIGINT

echo "================================================================================"
echo "   FASTQ DOWNLOADER (Slurm Job)"
echo "   Starting at: $(date)"
echo "================================================================================"
echo "Input file:  $INPUT_FILE"
echo "Output dir:  $OUTPUT_DIR"

if [ ! -f "$INPUT_FILE" ]; then
    echo "Error: Input URL file not found: $INPUT_FILE"
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

# --- Build list of URLs still needing download ---
TEMP_URLS=$(mktemp /tmp/ascli_urls_XXXXXX.txt)

while IFS= read -r url || [ -n "$url" ]; do
    [ -z "$url" ] && continue
    filename=$(basename "$url")

    # Skip if already fully downloaded (exact filename present)
    if [ -f "$OUTPUT_DIR/$filename" ]; then
        echo "[Skip] Fully downloaded but not yet renamed: $filename"
        continue
    fi

    echo "$url" >> "$TEMP_URLS"
done < "$INPUT_FILE"

URL_COUNT=$(wc -l < "$TEMP_URLS" | tr -d ' ')
echo "Found $URL_COUNT file(s) to download."
echo "-------------------------------------------------------------------------"

if [ "$URL_COUNT" -gt 0 ]; then

LOG_DIR="$PROJECT_ROOT/_data/logs/download"
mkdir -p "$LOG_DIR"

    # Check for download tool
    if ! command -v ascli &> /dev/null; then
        echo "Error: ascli (from the aspera-cli package) is not found in the PATH."
        echo "Please ensure the active mamba/conda environment has 'aspera-cli' installed."
        rm -f "$TEMP_URLS"
        exit 1
    fi

    echo ">>> Using ascli (Aspera) for download..."
    
    PARALLEL_DOWNLOADS=4
    CURRENT_SAMPLE=0

    while read -r url; do
        CURRENT_SAMPLE=$((CURRENT_SAMPLE + 1))
        filename=$(basename "$url")
        rel_path=$(echo "$url" | sed -E 's|^https?://ftp.sra.ebi.ac.uk/||')
        
        # Simple background job management
        while [ $(jobs -r | wc -l) -ge "$PARALLEL_DOWNLOADS" ]; do
            sleep 2
        done
        
        echo ">>> [$(date '+%Y-%m-%d %H:%M:%S')] Submitting download [$CURRENT_SAMPLE/$URL_COUNT]: $filename"
        SAMPLE_LOG="$LOG_DIR/${filename}.download.log"
        
        # Run download in background, capturing start/end time inside individual log file
        (
            echo ">>> [$(date '+%Y-%m-%d %H:%M:%S')] Starting download of $rel_path..."
            start_time=$(date +%s)
            
            ascli -Pera server download "$rel_path" --to-folder="$OUTPUT_DIR"
            
            end_time=$(date +%s)
            duration=$((end_time - start_time))
            echo ">>> [$(date '+%Y-%m-%d %H:%M:%S')] Finished download of $rel_path in ${duration}s"
        ) > "$SAMPLE_LOG" 2>&1 &
        
    done < "$TEMP_URLS"

    echo "Waiting for all background downloads to complete..."
    wait

fi

rm -f "$TEMP_URLS"
echo "-------------------------------------------------------------------------"
echo "=== FastQ Download Complete ==="
