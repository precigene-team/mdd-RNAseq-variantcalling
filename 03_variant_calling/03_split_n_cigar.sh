#!/bin/bash
# [[03_variant_calling/03_split_n_cigar.sh]]
# Goal: Split splice-junction reads containing Ns in their CIGAR string using GATK SplitNCigarReads.

set -e
set -o pipefail

# Inputs
SAMPLE_NAME="$1"
CPUS="${2:-1}"
MEM="${3:-8G}"

if [ -z "$SAMPLE_NAME" ]; then
    echo "ERROR: SAMPLE_NAME argument is missing." >&2
    echo "Usage: bash 03_split_n_cigar.sh <sample_name> [cpus] [memory]" >&2
    exit 1
fi

INPUT_DIR="$PROJECT_ROOT/_data/variant/duplicates"
OUTPUT_DIR="$PROJECT_ROOT/_data/variant/split_n_cigar"
LOG_DIR="$PROJECT_ROOT/_data/logs/split_n_cigar"
GENOME_DIR="$PROJECT_ROOT/_data/genome"
FASTA_FILE="$GENOME_DIR/Homo_sapiens.GRCh38.dna.primary_assembly.fa"

input_bam="$INPUT_DIR/${SAMPLE_NAME}.marked.bam"

if [ ! -f "$input_bam" ]; then
    echo "ERROR: Duplicate-marked BAM file for sample $SAMPLE_NAME not found at $input_bam." >&2
    exit 1
fi

out_bam="$OUTPUT_DIR/${SAMPLE_NAME}.split.bam"
sample_log="$LOG_DIR/${SAMPLE_NAME}.split.log"

if [ -s "$out_bam" ] && ([ -s "${out_bam}.bai" ] || [ -s "${out_bam%.bam}.bai" ]); then
    echo "[Skip] $SAMPLE_NAME (Split CIGAR BAM already exists)"
    exit 0
fi

echo ">>> Running SplitNCigarReads: $SAMPLE_NAME"
gatk --java-options "-Xmx${MEM} -Dsamjdk.use_async_io_read_samtools=true -Dsamjdk.use_async_io_write_samtools=true -XX:+UseG1GC -XX:ParallelGCThreads=${CPUS}" SplitNCigarReads \
    -R "$FASTA_FILE" \
    -I "$input_bam" \
    -O "$out_bam" \
    --create-output-bam-index true > "$sample_log" 2>&1
