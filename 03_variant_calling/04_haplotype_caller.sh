#!/bin/bash
# [[03_variant_calling/04_haplotype_caller.sh]]
# Goal: Run GATK HaplotypeCaller to call raw variants (SNPs and Indels) from processed BAM files.

set -e
set -o pipefail

# Inputs
SAMPLE_NAME="$1"
CPUS="${2:-1}"
MEM="${3:-8G}"

if [ -z "$SAMPLE_NAME" ]; then
    echo "ERROR: SAMPLE_NAME argument is missing." >&2
    echo "Usage: bash 04_haplotype_caller.sh <sample_name> [cpus] [memory]" >&2
    exit 1
fi

INPUT_DIR="$PROJECT_ROOT/_data/variant/split_n_cigar"
OUTPUT_DIR="$PROJECT_ROOT/_data/variant/raw_vcf"
LOG_DIR="$PROJECT_ROOT/_data/logs/haplotype_caller"
GENOME_DIR="$PROJECT_ROOT/_data/genome"
FASTA_FILE="$GENOME_DIR/Homo_sapiens.GRCh38.dna.primary_assembly.fa"

input_bam="$INPUT_DIR/${SAMPLE_NAME}.split.bam"

if [ ! -f "$input_bam" ]; then
    echo "ERROR: Split-CIGAR BAM file for sample $SAMPLE_NAME not found at $input_bam." >&2
    exit 1
fi

out_vcf="$OUTPUT_DIR/${SAMPLE_NAME}.raw.vcf.gz"
sample_log="$LOG_DIR/${SAMPLE_NAME}.hc.log"

if [ -s "$out_vcf" ]; then
    echo "[Skip] $SAMPLE_NAME (Raw VCF already exists)"
    exit 0
fi

# Ensure input split BAM is indexed
if [ ! -f "${input_bam}.bai" ] && [ ! -f "${input_bam%.bam}.bai" ]; then
    echo "[INFO] Index for $input_bam not found. Generating index..." >> "$sample_log" 2>&1
    samtools index "$input_bam" >> "$sample_log" 2>&1
fi

echo ">>> Running HaplotypeCaller: $SAMPLE_NAME"
gatk --java-options "-Xmx${MEM} -Dsamjdk.use_async_io_read_samtools=true -Dsamjdk.use_async_io_write_samtools=true -XX:ParallelGCThreads=${CPUS}" HaplotypeCaller \
    -R "$FASTA_FILE" \
    -I "$input_bam" \
    -O "$out_vcf" \
    --dont-use-soft-clipped-bases \
    -stand-call-conf 20.0 > "$sample_log" 2>&1
