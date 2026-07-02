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

# --- Resource Auto-detection ---
TOTAL_CPUS=${SLURM_CPUS_PER_TASK:-$(nproc 2>/dev/null || echo 16)}
if [ -n "$SLURM_MEM_PER_NODE" ]; then
    TOTAL_MEM_GB=$(($SLURM_MEM_PER_NODE / 1024))
else
    TOTAL_MEM_GB=$(free -g 2>/dev/null | awk '/^Mem:/{print $2}' || echo 128)
fi

# We target 4 CPUs per job for GATK tool execution
CPUS_PER_JOB=4
PARALLEL_JOBS=$((TOTAL_CPUS / CPUS_PER_JOB))
[ "$PARALLEL_JOBS" -lt 1 ] && PARALLEL_JOBS=1

# Ensure we have at least 8 GB RAM per job for JVM instances
MIN_MEM_PER_JOB=8
MEM_JOBS=$((TOTAL_MEM_GB / MIN_MEM_PER_JOB))
[ "$MEM_JOBS" -lt 1 ] && MEM_JOBS=1

if [ "$PARALLEL_JOBS" -gt "$MEM_JOBS" ]; then
    PARALLEL_JOBS=$MEM_JOBS
fi

# Recalculate requirements based on actual concurrent jobs
CPUS_PER_JOB=$((TOTAL_CPUS / PARALLEL_JOBS))
[ "$CPUS_PER_JOB" -lt 1 ] && CPUS_PER_JOB=1
MEM_PER_JOB_GB=$((TOTAL_MEM_GB / PARALLEL_JOBS))
[ "$MEM_PER_JOB_GB" -lt 4 ] && MEM_PER_JOB_GB=4
MEM_PER_JOB="${MEM_PER_JOB_GB}G"

echo "Resources allocated: $PARALLEL_JOBS concurrent jobs, each with $CPUS_PER_JOB threads and ${MEM_PER_JOB} memory."

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
