#!/bin/bash
# [[03_variant_calling/05_variant_filtration.sh]]
# Goal: Apply GATK standard hard-filtration filters for RNA-seq variants.

set -e
set -o pipefail

# Inputs
SAMPLE_NAME="$1"
CPUS="${2:-1}"
MEM="${3:-8G}"

if [ -z "$SAMPLE_NAME" ]; then
    echo "ERROR: SAMPLE_NAME argument is missing." >&2
    echo "Usage: bash 05_variant_filtration.sh <sample_name> [cpus] [memory]" >&2
    exit 1
fi

INPUT_DIR="$PROJECT_ROOT/_data/variant/raw_vcf"
RAW_SNPS_DIR="$PROJECT_ROOT/_data/variant/snps_raw"
OUTPUT_DIR="$PROJECT_ROOT/_data/variant/filtered_vcf"
LOG_DIR="$PROJECT_ROOT/_data/logs/variant_filtration"
GENOME_DIR="$PROJECT_ROOT/_data/genome"
FASTA_FILE="$GENOME_DIR/Homo_sapiens.GRCh38.dna.primary_assembly.fa"

input_vcf="$INPUT_DIR/${SAMPLE_NAME}.raw.vcf.gz"

if [ ! -f "$input_vcf" ]; then
    echo "ERROR: Raw VCF file for sample $SAMPLE_NAME not found at $input_vcf." >&2
    exit 1
fi

out_vcf="$OUTPUT_DIR/${SAMPLE_NAME}.filtered.vcf.gz"
raw_snps_vcf="$RAW_SNPS_DIR/${SAMPLE_NAME}.snps.raw.vcf.gz"
sample_log="$LOG_DIR/${SAMPLE_NAME}.filter.log"

# Check if local workstation is processing this sample (lock file exists)
lock_file="$LOG_DIR/${SAMPLE_NAME}.local_lock"
if [ -f "$lock_file" ]; then
    echo "[Skip] $SAMPLE_NAME (Local workstation is processing this sample; locked)"
    exit 0
fi

if [ -s "$out_vcf" ]; then
    echo "[Skip] $SAMPLE_NAME (Filtered VCF already exists)"
    exit 0
fi

echo ">>> Running SelectVariants & VariantFiltration: $SAMPLE_NAME"

(
    # 1. Select only SNPs from raw VCF
    gatk --java-options "-Xmx${MEM} -Dsamjdk.use_async_io_read_samtools=true -Dsamjdk.use_async_io_write_samtools=true -XX:ParallelGCThreads=${CPUS}" SelectVariants \
        -R "$FASTA_FILE" \
        -V "$input_vcf" \
        -select-type SNP \
        -O "$raw_snps_vcf" && \
        
    # 2. Hard-filter the SNPs
    gatk --java-options "-Xmx${MEM} -Dsamjdk.use_async_io_read_samtools=true -Dsamjdk.use_async_io_write_samtools=true -XX:ParallelGCThreads=${CPUS}" VariantFiltration \
        -R "$FASTA_FILE" \
        -V "$raw_snps_vcf" \
        -O "$out_vcf" \
        --window 35 \
        --cluster 3 \
        --filter-name "FS" \
        --filter-expression "FS > 30.0" \
        --filter-name "QD" \
        --filter-expression "QD < 2.0"
) > "$sample_log" 2>&1
