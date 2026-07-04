#!/bin/bash
# [[03_variant_calling/01_read_groups.sh]]
# Goal: Add Read Group information to STAR-aligned coordinate-sorted BAM files using GATK.

set -e
set -o pipefail

# Inputs
SAMPLE_NAME="$1"
CPUS="${2:-1}"
MEM="${3:-8G}"

if [ -z "$SAMPLE_NAME" ]; then
    echo "ERROR: SAMPLE_NAME argument is missing." >&2
    echo "Usage: bash 01_read_groups.sh <sample_name> [cpus] [memory]" >&2
    exit 1
fi

BAM_DIR="$PROJECT_ROOT/_data/bam"
OUTPUT_DIR="$PROJECT_ROOT/_data/variant/readgroups"
LOG_DIR="$PROJECT_ROOT/_data/logs/readgroups"

# Resolve the exact path of the input BAM file (it resides inside a sample-specific subdirectory)
input_bam=$(find "$BAM_DIR" -name "${SAMPLE_NAME}_Aligned.sortedByCoord.out.bam" | head -n 1)

if [ -z "$input_bam" ] || [ ! -f "$input_bam" ]; then
    echo "ERROR: Coordinate-sorted BAM file for sample $SAMPLE_NAME not found under $BAM_DIR." >&2
    exit 1
fi

out_bam="$OUTPUT_DIR/${SAMPLE_NAME}.rg.bam"
sample_log="$LOG_DIR/${SAMPLE_NAME}.rg.log"

if [ -s "$out_bam" ] && ([ -s "${out_bam}.bai" ] || [ -s "${out_bam%.bam}.bai" ]); then
    echo "[Skip] $SAMPLE_NAME (Read group BAM already exists)"
    exit 0
fi

echo ">>> Running AddOrReplaceReadGroups: $SAMPLE_NAME"
gatk --java-options "-Xmx${MEM} -Dsamjdk.use_async_io_read_samtools=true -Dsamjdk.use_async_io_write_samtools=true -XX:ParallelGCThreads=${CPUS}" AddOrReplaceReadGroups \
    -I "$input_bam" \
    -O "$out_bam" \
    -RGID "$SAMPLE_NAME" \
    -RGLB lib1 \
    -RGPL ILLUMINA \
    -RGPU unit1 \
    -RGSM "$SAMPLE_NAME" \
    --CREATE_INDEX true > "$sample_log" 2>&1
