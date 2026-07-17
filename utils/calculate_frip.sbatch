#!/bin/bash
# [[utils/calculate_frip.sbatch]]
# Goal: Calculate FRiP (Fraction of Reads in Peaks) scores for all subjects in the cohort using TRESS outputs.

# Steps:
#   1. Resolve project root and key paths.
#   2. Assert sample sheet exists.
#   3. Build BAM indexes in parallel if missing.
#   4. For each subject: count total mapped reads, count reads overlapping TRESS peaks, compute FRiP.
#   5. Write results TSV to results/utils/tables/.

#SBATCH --job-name=qc01_frip
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=_data/logs/frip/qc01_%j.txt
#SBATCH --error=_data/logs/frip/qc01_%j.err

set -e
set -o pipefail

# --- Project Root Discovery ---
if [ -z "$PROJECT_ROOT" ]; then
    curr_dir="$PWD"
    while [ "$curr_dir" != "/" ]; do
        if [ -f "$curr_dir/run" ] && [ -d "$curr_dir/01_core_alignment" ]; then
            export PROJECT_ROOT="$curr_dir"
            break
        fi
        curr_dir="$(dirname "$curr_dir")"
    done
    export PROJECT_ROOT="${PROJECT_ROOT:-$PWD}"
fi

# --- Safe Conda/Mamba Hook Initialization ---
eval "$(micromamba shell hook --shell bash 2>/dev/null || mamba shell hook --shell bash 2>/dev/null || conda shell.bash hook 2>/dev/null)"
micromamba activate mdd 2>/dev/null || mamba activate mdd 2>/dev/null || conda activate mdd

BAM_DIR="$PROJECT_ROOT/_data/bam"
PEAK_DIR="$PROJECT_ROOT/_data/peaks"
CSV_FILE="$PROJECT_ROOT/_data/_metadata/sample_sheet.csv"
OUTPUT_TABLE="$PROJECT_ROOT/results/utils/tables/frip_scores.tsv"
LOG_DIR="$PROJECT_ROOT/_data/logs/frip"

mkdir -p "$(dirname "$OUTPUT_TABLE")"
mkdir -p "$LOG_DIR"

if [ ! -f "$CSV_FILE" ]; then
    echo "[ERROR] Sample sheet not found at $CSV_FILE" >&2
    exit 1
fi

echo "=========================================================================="
echo "   FRiP Calculation started at: $(date)"
echo "   Sample sheet: $CSV_FILE"
echo "   Output table: $OUTPUT_TABLE"
echo "=========================================================================="

# Get unique subjects
subjects=$(awk -F',' 'NR>1 {print $1}' "$CSV_FILE" | sed -E 's/_(m6A|input)_1//' | sort -u)

# --- 1. Parallel Indexing Phase ---
echo ">>> Checking and building missing BAM indexes in parallel..."
max_jobs=4
job_count=0

for subj in $subjects; do
    # Skip gated samples
    if [[ "$subj" == "female_mdd_2" || "$subj" == "female_ctrl_6" || "$subj" == "male_ctrl_1" ]]; then
        continue
    fi

    IP_BAM="$BAM_DIR/${subj}_m6A/${subj}_m6A_Aligned.sortedByCoord.out.bam"

    if [ -f "$IP_BAM" ] && [ ! -s "${IP_BAM}.bai" ] && [ ! -s "${IP_BAM%.bam}.bai" ]; then
        echo "Indexing $IP_BAM in background..."
        samtools index -@ 1 "$IP_BAM" &
        job_count=$((job_count + 1))

        if [ "$job_count" -ge "$max_jobs" ]; then
            wait -n 2>/dev/null || wait
            job_count=$((job_count - 1))
        fi
    fi
done
wait
echo ">>> All BAM indexes are ready."

# Write header
echo -e "Subject\tTotal_Reads\tReads_in_Peaks\tFRiP_Score" | tee "$OUTPUT_TABLE"

for subj in $subjects; do
    # Skip gated samples (low coverage)
    if [[ "$subj" == "female_mdd_2" || "$subj" == "female_ctrl_6" || "$subj" == "male_ctrl_1" ]]; then
        continue
    fi

    IP_BAM="$BAM_DIR/${subj}_m6A/${subj}_m6A_Aligned.sortedByCoord.out.bam"
    PEAK_FILE="$PEAK_DIR/${subj}/${subj}_peaks.xls"

    if [ ! -f "$IP_BAM" ]; then
        echo -e "${subj}\t-\t-\t[NO_BAM]" | tee -a "$OUTPUT_TABLE"
        continue
    fi

    if [ ! -s "$PEAK_FILE" ]; then
        echo -e "${subj}\t-\t-\t[NO_PEAKS]" | tee -a "$OUTPUT_TABLE"
        continue
    fi

    total_reads=$(samtools view -@ "${SLURM_CPUS_PER_TASK:-4}" -c -F 260 "$IP_BAM")

    if [ "$total_reads" -eq 0 ]; then
        echo -e "${subj}\t0\t0\t0.0000" | tee -a "$OUTPUT_TABLE"
        continue
    fi

    # Check chr prefix convention
    if samtools view -H "$IP_BAM" | head -n 100 | grep -q "SN:chr"; then
        BAM_USES_CHR=true
    else
        BAM_USES_CHR=false
    fi

    # Convert TRESS XLS to clean BED
    TMP_BED=$(mktemp)
    awk -v use_chr="$BAM_USES_CHR" '
    NR>1 {
        gsub(/"/, "", $1)
        gsub(/"/, "", $2)
        gsub(/"/, "", $3)
        chrom = $1
        if (use_chr == "true" && chrom !~ /^chr/) chrom = "chr" chrom
        else if (use_chr == "false" && chrom ~ /^chr/) sub(/^chr/, "", chrom)
        if (chrom != "" && $2 != "" && $3 != "") print chrom "\t" $2 "\t" $3
    }' "$PEAK_FILE" > "$TMP_BED"

    if [ -s "$TMP_BED" ]; then
        reads_in_peaks=$(samtools view -@ "${SLURM_CPUS_PER_TASK:-4}" -c -F 260 -L "$TMP_BED" "$IP_BAM")
    else
        reads_in_peaks=0
    fi
    rm -f "$TMP_BED"

    frip=$(awk -v rip="$reads_in_peaks" -v tot="$total_reads" 'BEGIN { if (tot > 0) printf "%.4f", rip/tot; else print "0.0000" }')
    echo -e "${subj}\t${total_reads}\t${reads_in_peaks}\t${frip}" | tee -a "$OUTPUT_TABLE"
done

echo "=========================================================================="
echo "   FRiP Calculation completed at: $(date)"
echo "   Results: $OUTPUT_TABLE"
echo "=========================================================================="
