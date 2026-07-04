#!/bin/bash
# [[03_variant_calling/02_mark_duplicates.sh]]
# Goal: Flag PCR duplication artifacts in BAM files using GATK MarkDuplicates.

set -e
set -o pipefail

# Inputs
SAMPLE_NAME="$1"
CPUS="${2:-1}"
MEM="${3:-8G}"

if [ -z "$SAMPLE_NAME" ]; then
    echo "ERROR: SAMPLE_NAME argument is missing." >&2
    echo "Usage: bash 02_mark_duplicates.sh <sample_name> [cpus] [memory]" >&2
    exit 1
fi

INPUT_DIR="$PROJECT_ROOT/_data/variant/readgroups"
OUTPUT_DIR="$PROJECT_ROOT/_data/variant/duplicates"
LOG_DIR="$PROJECT_ROOT/_data/logs/duplicates"
TMP_DIR="${SLURM_TMPDIR:-$PROJECT_ROOT/_data/tmp}"

input_bam="$INPUT_DIR/${SAMPLE_NAME}.rg.bam"

if [ ! -f "$input_bam" ]; then
    echo "ERROR: Read-grouped BAM file for sample $SAMPLE_NAME not found at $input_bam." >&2
    exit 1
fi

out_bam="$OUTPUT_DIR/${SAMPLE_NAME}.marked.bam"
metrics_file="$OUTPUT_DIR/${SAMPLE_NAME}.metrics.txt"
sample_log="$LOG_DIR/${SAMPLE_NAME}.dup.log"

if [ -s "$out_bam" ] && ([ -s "${out_bam}.bai" ] || [ -s "${out_bam%.bam}.bai" ]); then
    echo "[Skip] $SAMPLE_NAME (Duplicate-marked BAM already exists)"
    exit 0
fi

echo ">>> Running MarkDuplicates: $SAMPLE_NAME"
gatk --java-options "-Xmx${MEM} -Djava.io.tmpdir=${TMP_DIR} -Dsamjdk.use_async_io_read_samtools=true -Dsamjdk.use_async_io_write_samtools=true -XX:ParallelGCThreads=${CPUS}" MarkDuplicates \
    -I "$input_bam" \
    -O "$out_bam" \
    -M "$metrics_file" \
    --CREATE_INDEX true \
    --READ_NAME_REGEX null > "$sample_log" 2>&1
