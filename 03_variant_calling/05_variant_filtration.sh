# [[03_variant_calling/05_variant_filtration.sh]]
# Goal: Apply GATK standard hard-filtration filters for RNA-seq variants.

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

INPUT_DIR="$PROJECT_ROOT/_data/variant/raw_vcf"
RAW_SNPS_DIR="$PROJECT_ROOT/_data/variant/snps_raw"
OUTPUT_DIR="$PROJECT_ROOT/_data/variant/filtered_vcf"
LOG_DIR="$PROJECT_ROOT/_data/logs/variant_filtration"
GENOME_DIR="$PROJECT_ROOT/_data/genome"
FASTA_FILE="$GENOME_DIR/Homo_sapiens.GRCh38.dna.primary_assembly.fa"

mkdir -p "$RAW_SNPS_DIR"
mkdir -p "$OUTPUT_DIR"
mkdir -p "$LOG_DIR"

# Resource configuration (adjust for HPC node capacity)
PARALLEL_JOBS=4
CPUS_PER_JOB=2
MEM_PER_JOB="4G"

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

# Find all raw VCF files
VCF_FILES=$(find "$INPUT_DIR" -name "*.raw.vcf.gz" | sort)
VCF_COUNT=$(echo "$VCF_FILES" | grep -v '^$' | wc -l || echo 0)

if [ "$VCF_COUNT" -eq 0 ]; then
    echo "Warning: No raw VCF files found in $INPUT_DIR."
    exit 0
fi

echo "[INFO] Found $VCF_COUNT raw VCF files to filter."

CURRENT_SAMPLE=0
for vcf in $VCF_FILES; do
    [ -f "$vcf" ] || continue
    CURRENT_SAMPLE=$((CURRENT_SAMPLE + 1))
    
    sample_name=$(basename "$vcf" .raw.vcf.gz)
    out_vcf="$OUTPUT_DIR/${sample_name}.filtered.vcf.gz"
    sample_log="$LOG_DIR/${sample_name}.filter.log"
    
    if [ -f "$out_vcf" ]; then
        echo "[Skip] [$CURRENT_SAMPLE/$VCF_COUNT] $sample_name (Filtered VCF already exists)"
        continue
    fi
    
    # Simple background job management
    while [ $(jobs -r | wc -l) -ge "$PARALLEL_JOBS" ]; do
        sleep 5
    done
    
    raw_snps_vcf="$RAW_SNPS_DIR/${sample_name}.snps.raw.vcf.gz"
    
    echo ">>> Running SelectVariants & VariantFiltration [$CURRENT_SAMPLE/$VCF_COUNT]: $sample_name"
    
    (
        # 1. Select only SNPs from raw VCF
        gatk --java-options "-Xmx${MEM_PER_JOB} -XX:ParallelGCThreads=${CPUS_PER_JOB}" SelectVariants \
            -R "$FASTA_FILE" \
            -V "$vcf" \
            -select-type SNP \
            -O "$raw_snps_vcf" && \
            
        # 2. Hard-filter the SNPs
        gatk --java-options "-Xmx${MEM_PER_JOB} -XX:ParallelGCThreads=${CPUS_PER_JOB}" VariantFiltration \
            -R "$FASTA_FILE" \
            -V "$raw_snps_vcf" \
            -O "$out_vcf" \
            --window 35 \
            --cluster 3 \
            --filter-name "FS" \
            --filter-expression "FS > 30.0" \
            --filter-name "QD" \
            --filter-expression "QD < 2.0"
    ) > "$sample_log" 2>&1 &
done

echo "Waiting for all VariantFiltration jobs to complete..."
wait
echo "[SUCCESS] VariantFiltration complete for all samples."
