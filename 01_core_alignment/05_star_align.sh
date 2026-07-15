#!/bin/bash
# [[01_core_alignment/05_star_align.sh]]
# Goal: Build STAR index (if missing) and map trimmed FastQ files.

# Steps:
#   1. Validate inputs and environment (sample sheet existence, genome FASTA and GTF existence, active mdd environment).
#   2. Detect CPU resources and set parallelism parameters.
#   3. Build STAR genome index if not already present.
#   4. Parse sample metadata from sample sheet and map input trimmed paired FASTQ files.
#   5. Run STAR alignment in parallel (2-pass mode) to generate sorted coordinate BAM files and gene count tables.
#   6. Assert that all expected BAM outputs are successfully generated and valid.

# Exit on error
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

CSV_FILE="$PROJECT_ROOT/_data/_metadata/sample_sheet.csv"
TRIMMED_DIR="$PROJECT_ROOT/_data/trimmed_fastq"        # m6A samples
RRNA_FILTERED_DIR="$PROJECT_ROOT/_data/rrna_filtered"  # input samples (rRNA-depleted)
OUTPUT_DIR="$PROJECT_ROOT/_data/bam"

GENOME_DIR="$PROJECT_ROOT/_data/genome"
STAR_INDEX="$GENOME_DIR/star_index"
FASTA_FILE="$GENOME_DIR/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
GTF_FILE="$GENOME_DIR/Homo_sapiens.GRCh38.115.gtf"

LOG_DIR="$PROJECT_ROOT/_data/logs/star"

mkdir -p "$OUTPUT_DIR"
mkdir -p "$LOG_DIR"

# --- Input Assertions ---
if [ ! -f "$CSV_FILE" ]; then
    echo "Assertion Error: Sample sheet not found at $CSV_FILE" >&2
    exit 1
fi
if [ ! -f "$FASTA_FILE" ]; then
    echo "Assertion Error: Genome FASTA not found at $FASTA_FILE" >&2
    exit 1
fi
if [ ! -f "$GTF_FILE" ]; then
    echo "Assertion Error: GTF annotation not found at $GTF_FILE" >&2
    exit 1
fi

# --- Resource Auto-detection ---
TOTAL_CPUS=${SLURM_CPUS_PER_TASK:-32}
PARALLEL_JOBS=2
CPUS_PER_JOB=$((TOTAL_CPUS / PARALLEL_JOBS))
[ "$CPUS_PER_JOB" -lt 1 ] && CPUS_PER_JOB=1

# RAM limit for STAR sorting
SORT_RAM=28000000000 # 28 GB per job

echo "================================================================================"
echo "   STAR ALIGNMENT: Mapping trimmed paired FastQ files"
echo "================================================================================"
echo "Genome index: $STAR_INDEX"
echo "BAM output: $OUTPUT_DIR"
echo "Resources: $PARALLEL_JOBS concurrent samples, $CPUS_PER_JOB threads per sample."

# 1. Genome Index generation if not present
if [ ! -d "$STAR_INDEX" ] || [ ! -f "$STAR_INDEX/Genome" ]; then
    echo "STAR index not found. Generating index..."
    mkdir -p "$STAR_INDEX"
    STAR \
        --runThreadN "$TOTAL_CPUS" \
        --runMode genomeGenerate \
        --genomeDir "$STAR_INDEX" \
        --genomeFastaFiles "$FASTA_FILE" \
        --sjdbGTFfile "$GTF_FILE" \
        --sjdbOverhang 149
    echo "STAR index generation complete."
fi

# Assert genome index exists
if [ ! -f "$STAR_INDEX/Genome" ] || [ ! -s "$STAR_INDEX/Genome" ]; then
    echo "Assertion Error: STAR index generation failed. Index file is missing or empty." >&2
    exit 1
fi

# 2. Alignment loop
# Col 1=file1_base, Col 2=file2_base, Col 3=accession, Col 9=DataType (m6A|input)
SAMPLE_ROWS=$(awk -F',' 'NR>1 {print $3 "\t" $1 "\t" $2 "\t" $9}' "$CSV_FILE" | grep -v '^$')
SAMPLE_COUNT=$(echo "$SAMPLE_ROWS" | wc -l)

# Assert that we parsed at least one sample row
if [ "$SAMPLE_COUNT" -eq 0 ]; then
    echo "Assertion Error: No sample rows parsed from $CSV_FILE" >&2
    exit 1
fi

# Assert that at least one paired fastq exists
FOUND_INPUTS=0
while read -r row; do
    [ -n "$row" ] || continue
    file1_base=$(echo "$row" | cut -f2)
    file2_base=$(echo "$row" | cut -f3)
    data_type=$(echo "$row" | cut -f4)

    r1_test="$RRNA_FILTERED_DIR/${file1_base}_filtered.fastq.gz"
    r2_test="$RRNA_FILTERED_DIR/${file2_base}_filtered.fastq.gz"

    if [ -f "$r1_test" ] && [ -f "$r2_test" ]; then
        FOUND_INPUTS=$((FOUND_INPUTS + 1))
    fi
done <<< "$SAMPLE_ROWS"

if [ "$FOUND_INPUTS" -eq 0 ]; then
    echo "Assertion Error: No matching pre-alignment FastQ files found in $TRIMMED_DIR or $RRNA_FILTERED_DIR." >&2
    exit 1
fi

CURRENT_SAMPLE=0
while read -r row; do
    [ -n "$row" ] || continue
    CURRENT_SAMPLE=$((CURRENT_SAMPLE + 1))
    
    accession=$(echo "$row" | cut -f1)
    file1_base=$(echo "$row" | cut -f2)
    file2_base=$(echo "$row" | cut -f3)
    data_type=$(echo "$row" | cut -f4)

    # --- Resolve source directory (both fractions are now rRNA-filtered) ---
    r1="$RRNA_FILTERED_DIR/${file1_base}_filtered.fastq.gz"
    r2="$RRNA_FILTERED_DIR/${file2_base}_filtered.fastq.gz"
    
    if [ ! -f "$r1" ] || [ ! -f "$r2" ]; then
        echo "[Warn] Missing input files for $accession: $r1 or $r2. Skipping."
        continue
    fi
    
    # Prefix naming logic: strips suffix '_1' if present
    sample_name=$(echo "$file1_base" | sed 's/_1$//')
    
    SAMPLE_OUT_DIR="$OUTPUT_DIR/${sample_name}"
    BAM_FILE="$SAMPLE_OUT_DIR/${sample_name}_Aligned.sortedByCoord.out.bam"
    
    # Check if BAM output already exists
    if [ -f "$BAM_FILE" ]; then
        echo "[Skip] [$CURRENT_SAMPLE/$SAMPLE_COUNT] $sample_name (BAM already exists)"
        continue
    fi
    
    # Simple background job management
    while [ $(jobs -r | wc -l) -ge "$PARALLEL_JOBS" ]; do
        sleep 5
    done
    
    echo ">>> Submitting [$CURRENT_SAMPLE/$SAMPLE_COUNT]: $sample_name"
    SAMPLE_LOG="$LOG_DIR/${sample_name}.star.log"
    
    # Run STAR in background
    (
        mkdir -p "$SAMPLE_OUT_DIR"
        STAR \
            --runThreadN "$CPUS_PER_JOB" \
            --genomeDir "$STAR_INDEX" \
            --readFilesIn "$r1" "$r2" \
            --readFilesCommand gunzip -c \
            --outFileNamePrefix "$SAMPLE_OUT_DIR/${sample_name}_" \
            --outSAMtype BAM SortedByCoordinate \
            --limitBAMsortRAM "$SORT_RAM" \
            --outFilterType BySJout \
            --outFilterMultimapNmax 20 \
            --alignSJoverhangMin 8 \
            --alignSJDBoverhangMin 1 \
            --quantMode GeneCounts
    ) > "$SAMPLE_LOG" 2>&1 &
        
done <<< "$SAMPLE_ROWS"

echo "Waiting for all STAR alignment jobs to complete..."
wait

# --- Output Assertions ---
echo "Running post-execution assertions..."
MISSING_OUTPUTS=0
while read -r row; do
    [ -n "$row" ] || continue
    file1_base=$(echo "$row" | cut -f2)
    file2_base=$(echo "$row" | cut -f3)
    data_type=$(echo "$row" | cut -f4)

    r1="$RRNA_FILTERED_DIR/${file1_base}_filtered.fastq.gz"
    r2="$RRNA_FILTERED_DIR/${file2_base}_filtered.fastq.gz"
    [ -f "$r1" ] && [ -f "$r2" ] || continue # If inputs were missing, outputs won't exist
    
    sample_name=$(echo "$file1_base" | sed 's/_1$//')
    SAMPLE_OUT_DIR="$OUTPUT_DIR/${sample_name}"
    BAM_FILE="$SAMPLE_OUT_DIR/${sample_name}_Aligned.sortedByCoord.out.bam"
    
    if [ ! -s "$BAM_FILE" ]; then
        echo "Assertion Error: STAR alignment output BAM file missing or empty for sample $sample_name" >&2
        MISSING_OUTPUTS=$((MISSING_OUTPUTS + 1))
    fi
done <<< "$SAMPLE_ROWS"

if [ "$MISSING_OUTPUTS" -gt 0 ]; then
    echo "Assertion Error: $MISSING_OUTPUTS sample(s) failed STAR alignment validation." >&2
    exit 1
fi

echo "All outputs successfully validated."
echo "=== STAR Alignment Complete ==="
