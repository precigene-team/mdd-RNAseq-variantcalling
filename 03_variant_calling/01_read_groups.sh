# [[03_variant_calling/01_read_groups.sh]]
# Goal: Add Read Group information to STAR-aligned coordinate-sorted BAM files using GATK.

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

BAM_DIR="$PROJECT_ROOT/_data/bam"
OUTPUT_DIR="$PROJECT_ROOT/_data/variant/readgroups"
LOG_DIR="$PROJECT_ROOT/_data/logs/readgroups"

mkdir -p "$OUTPUT_DIR"
mkdir -p "$LOG_DIR"

# Resource configuration (adjust for HPC node capacity)
PARALLEL_JOBS=4
CPUS_PER_JOB=2
MEM_PER_JOB="8G"

# Find all aligned coordinate-sorted BAM files
BAM_FILES=$(find "$BAM_DIR" -name "*_Aligned.sortedByCoord.out.bam" | sort)
BAM_COUNT=$(echo "$BAM_FILES" | grep -v '^$' | wc -l || echo 0)

if [ "$BAM_COUNT" -eq 0 ]; then
    echo "Warning: No coordinate-sorted BAM files found in $BAM_DIR."
    exit 0
fi

echo "[INFO] Found $BAM_COUNT BAM files to process."

CURRENT_SAMPLE=0
for bam in $BAM_FILES; do
    [ -f "$bam" ] || continue
    CURRENT_SAMPLE=$((CURRENT_SAMPLE + 1))
    
    sample_name=$(basename "$bam" _Aligned.sortedByCoord.out.bam)
    out_bam="$OUTPUT_DIR/${sample_name}.rg.bam"
    sample_log="$LOG_DIR/${sample_name}.rg.log"
    
    if [ -f "$out_bam" ]; then
        echo "[Skip] [$CURRENT_SAMPLE/$BAM_COUNT] $sample_name (Read group BAM already exists)"
        continue
    fi
    
    # Simple background job management
    while [ $(jobs -r | wc -l) -ge "$PARALLEL_JOBS" ]; do
        sleep 5
    done
    
    echo ">>> Running AddOrReplaceReadGroups [$CURRENT_SAMPLE/$BAM_COUNT]: $sample_name"
    
    (
        gatk --java-options "-Xmx${MEM_PER_JOB} -XX:ParallelGCThreads=${CPUS_PER_JOB}" AddOrReplaceReadGroups \
            -I "$bam" \
            -O "$out_bam" \
            -RGID "$sample_name" \
            -RGLB lib1 \
            -RGPL ILLUMINA \
            -RGPU unit1 \
            -RGSM "$sample_name" \
            --CREATE_INDEX true
    ) > "$sample_log" 2>&1 &
done

echo "Waiting for all AddOrReplaceReadGroups jobs to complete..."
wait
echo "[SUCCESS] AddOrReplaceReadGroups complete for all samples."
