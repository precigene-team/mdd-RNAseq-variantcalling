# [[03_variant_calling/04_haplotype_caller.sh]]
# Goal: Run GATK HaplotypeCaller to call raw variants (SNPs and Indels) from processed BAM files.

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

INPUT_DIR="$PROJECT_ROOT/_data/variant/split_n_cigar"
OUTPUT_DIR="$PROJECT_ROOT/_data/variant/raw_vcf"
LOG_DIR="$PROJECT_ROOT/_data/logs/haplotype_caller"
GENOME_DIR="$PROJECT_ROOT/_data/genome"
FASTA_FILE="$GENOME_DIR/Homo_sapiens.GRCh38.dna.primary_assembly.fa"

mkdir -p "$OUTPUT_DIR"
mkdir -p "$LOG_DIR"

# Resource configuration (adjust for HPC node capacity)
PARALLEL_JOBS=4
CPUS_PER_JOB=2
MEM_PER_JOB="8G"

# --- Genome Reference Verification ---
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

# Find all split-CIGAR BAM files
BAM_FILES=$(find "$INPUT_DIR" -name "*.split.bam" | sort)
BAM_COUNT=$(echo "$BAM_FILES" | grep -v '^$' | wc -l || echo 0)

if [ "$BAM_COUNT" -eq 0 ]; then
    echo "Warning: No split-CIGAR BAM files found in $INPUT_DIR."
    exit 0
fi

echo "[INFO] Found $BAM_COUNT BAM files to process."

CURRENT_SAMPLE=0
for bam in $BAM_FILES; do
    [ -f "$bam" ] || continue
    CURRENT_SAMPLE=$((CURRENT_SAMPLE + 1))
    
    sample_name=$(basename "$bam" .split.bam)
    out_vcf="$OUTPUT_DIR/${sample_name}.raw.vcf.gz"
    sample_log="$LOG_DIR/${sample_name}.hc.log"
    
    if [ -f "$out_vcf" ]; then
        echo "[Skip] [$CURRENT_SAMPLE/$BAM_COUNT] $sample_name (Raw VCF already exists)"
        continue
    fi
    
    # Simple background job management
    while [ $(jobs -r | wc -l) -ge "$PARALLEL_JOBS" ]; do
        sleep 5
    done
    
    # Ensure input split BAM is indexed
    if [ ! -f "${bam}.bai" ]; then
        echo "[INFO] Index for $bam not found. Generating index..."
        samtools index "$bam"
    fi
    
    echo ">>> Running HaplotypeCaller [$CURRENT_SAMPLE/$BAM_COUNT]: $sample_name"
    
    (
        gatk --java-options "-Xmx${MEM_PER_JOB} -XX:ParallelGCThreads=${CPUS_PER_JOB}" HaplotypeCaller \
            -R "$FASTA_FILE" \
            -I "$bam" \
            -O "$out_vcf" \
            --dont-use-soft-clipped-bases \
            -stand-call-conf 20.0
    ) > "$sample_log" 2>&1 &
done

echo "Waiting for all HaplotypeCaller jobs to complete..."
wait
echo "[SUCCESS] HaplotypeCaller complete for all samples."
