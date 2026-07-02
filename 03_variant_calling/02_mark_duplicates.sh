# [[03_variant_calling/02_mark_duplicates.sh]]
# Goal: Flag PCR duplication artifacts in BAM files using GATK MarkDuplicates.

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

INPUT_DIR="$PROJECT_ROOT/_data/variant/readgroups"
OUTPUT_DIR="$PROJECT_ROOT/_data/variant/duplicates"
LOG_DIR="$PROJECT_ROOT/_data/logs/duplicates"
TMP_DIR="$PROJECT_ROOT/_data/tmp"

mkdir -p "$OUTPUT_DIR"
mkdir -p "$LOG_DIR"
mkdir -p "$TMP_DIR"

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

# Find all read-grouped BAM files
BAM_FILES=$(find "$INPUT_DIR" -name "*.rg.bam" | sort)
BAM_COUNT=$(echo "$BAM_FILES" | grep -v '^$' | wc -l || echo 0)

if [ "$BAM_COUNT" -eq 0 ]; then
    echo "Warning: No read-grouped BAM files found in $INPUT_DIR."
    exit 0
fi

echo "[INFO] Found $BAM_COUNT BAM files to process."

CURRENT_SAMPLE=0
for bam in $BAM_FILES; do
    [ -f "$bam" ] || continue
    CURRENT_SAMPLE=$((CURRENT_SAMPLE + 1))
    
    sample_name=$(basename "$bam" .rg.bam)
    out_bam="$OUTPUT_DIR/${sample_name}.marked.bam"
    metrics_file="$OUTPUT_DIR/${sample_name}.metrics.txt"
    sample_log="$LOG_DIR/${sample_name}.dup.log"
    
    if [ -f "$out_bam" ]; then
        echo "[Skip] [$CURRENT_SAMPLE/$BAM_COUNT] $sample_name (Duplicate-marked BAM already exists)"
        continue
    fi
    
    # Simple background job management
    while [ $(jobs -r | wc -l) -ge "$PARALLEL_JOBS" ]; do
        sleep 5
    done
    
    echo ">>> Running MarkDuplicates [$CURRENT_SAMPLE/$BAM_COUNT]: $sample_name"
    
    (
        gatk --java-options "-Xmx${MEM_PER_JOB} -Djava.io.tmpdir=${TMP_DIR} -XX:ParallelGCThreads=${CPUS_PER_JOB}" MarkDuplicates \
            -I "$bam" \
            -O "$out_bam" \
            -M "$metrics_file" \
            --CREATE_INDEX true \
            --READ_NAME_REGEX null
    ) > "$sample_log" 2>&1 &
done

echo "Waiting for all MarkDuplicates jobs to complete..."
wait
echo "[SUCCESS] MarkDuplicates complete for all samples."
