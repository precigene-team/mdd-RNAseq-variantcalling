# [[03_variant_calling/03_split_n_cigar.sh]]
# Goal: Split splice-junction reads containing Ns in their CIGAR string using GATK SplitNCigarReads.

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

INPUT_DIR="$PROJECT_ROOT/_data/variant/duplicates"
OUTPUT_DIR="$PROJECT_ROOT/_data/variant/split_n_cigar"
LOG_DIR="$PROJECT_ROOT/_data/logs/split_n_cigar"
GENOME_DIR="$PROJECT_ROOT/_data/genome"
FASTA_FILE="$GENOME_DIR/Homo_sapiens.GRCh38.dna.primary_assembly.fa"

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

# --- Genome Reference Index Verification ---
if [ ! -f "$FASTA_FILE" ]; then
    echo "Assertion Error: Genome FASTA not found at $FASTA_FILE" >&2
    exit 1
fi

# Ensure samtools fai index exists
if [ ! -f "${FASTA_FILE}.fai" ]; then
    echo "[INFO] Genome FASTA index (.fai) not found. Building with samtools..."
    samtools faidx "$FASTA_FILE"
fi

# Ensure sequence dictionary exists
DICT_FILE="${FASTA_FILE%.fa}.dict"
if [ ! -f "$DICT_FILE" ]; then
    echo "[INFO] Sequence dictionary (.dict) not found. Building with GATK..."
    gatk CreateSequenceDictionary -R "$FASTA_FILE" -O "$DICT_FILE"
fi

# Find all duplicate-marked BAM files
BAM_FILES=$(find "$INPUT_DIR" -name "*.marked.bam" | sort)
BAM_COUNT=$(echo "$BAM_FILES" | grep -v '^$' | wc -l || echo 0)

if [ "$BAM_COUNT" -eq 0 ]; then
    echo "Warning: No duplicate-marked BAM files found in $INPUT_DIR."
    exit 0
fi

echo "[INFO] Found $BAM_COUNT BAM files to process."

CURRENT_SAMPLE=0
for bam in $BAM_FILES; do
    [ -f "$bam" ] || continue
    CURRENT_SAMPLE=$((CURRENT_SAMPLE + 1))
    
    sample_name=$(basename "$bam" .marked.bam)
    out_bam="$OUTPUT_DIR/${sample_name}.split.bam"
    sample_log="$LOG_DIR/${sample_name}.split.log"
    
    if [ -f "$out_bam" ]; then
        echo "[Skip] [$CURRENT_SAMPLE/$BAM_COUNT] $sample_name (Split CIGAR BAM already exists)"
        continue
    fi
    
    # Simple background job management
    while [ $(jobs -r | wc -l) -ge "$PARALLEL_JOBS" ]; do
        sleep 5
    done
    
    echo ">>> Running SplitNCigarReads [$CURRENT_SAMPLE/$BAM_COUNT]: $sample_name"
    
    (
        gatk --java-options "-Xmx${MEM_PER_JOB} -XX:ParallelGCThreads=${CPUS_PER_JOB}" SplitNCigarReads \
            -R "$FASTA_FILE" \
            -I "$bam" \
            -O "$out_bam" && \
        samtools index "$out_bam"
    ) > "$sample_log" 2>&1 &
done

echo "Waiting for all SplitNCigarReads jobs to complete..."
wait
echo "[SUCCESS] SplitNCigarReads complete for all samples."
