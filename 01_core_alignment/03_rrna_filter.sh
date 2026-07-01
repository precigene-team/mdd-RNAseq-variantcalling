#!/bin/bash
# [[01_core_alignment/03_rrna_filter.sh]]
# Goal: Remove rRNA reads from trimmed paired-end FastQ files using STAR pre-mapping.

# Steps:
#   1. Validate inputs and environment (sample sheet, active mdd env, STAR binary).
#   2. Check or build a small STAR index for the rRNA databases.
#   3. Parse sample metadata from sample sheet.
#   4. For each sample, run STAR mapping to rRNA index, saving only UNMAPPED reads.
#   5. Compress and rename unmapped reads to _data/rrna_filtered/ directory.
#   6. Clean up temporary mapping logs and folders.

# Exit on error
set -e
set -o pipefail

# --- Project Root Discovery ---
export PROJECT_ROOT="${PROJECT_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
# ------------------------------

CSV_FILE="$PROJECT_ROOT/_data/_metadata/sample_sheet.csv"
TRIMMED_DIR="$PROJECT_ROOT/_data/trimmed_fastq"
OUTDIR="$PROJECT_ROOT/_data/rrna_filtered"
LOG_DIR="$PROJECT_ROOT/_data/logs/rrna_filter"

GENOME_DIR="$PROJECT_ROOT/_data/genome"
RRNA_INDEX="$GENOME_DIR/rrna_star_index"
DB_FILE="$PROJECT_ROOT/_data/rrna_databases/smr_v4.3_default_db.fasta"

mkdir -p "$OUTDIR"
mkdir -p "$LOG_DIR"
mkdir -p "$RRNA_INDEX"

# --- Input Assertions ---
if [ ! -f "$CSV_FILE" ]; then
    echo "Assertion Error: Sample sheet not found at $CSV_FILE" >&2
    exit 1
fi

# --- Environment Assertions ---
if [ -z "$CONDA_PREFIX" ]; then
    echo "Error: Mamba environment 'mdd' is not active. Please activate it first." >&2
    exit 1
fi

if ! command -v STAR &> /dev/null; then
    echo "Error: STAR aligner not found in path." >&2
    exit 1
fi

# --- Paths to Reference Genome ---
FASTA_FILE="$GENOME_DIR/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
GTF_FILE="$GENOME_DIR/Homo_sapiens.GRCh38.115.gtf"
DB_FILE="$PROJECT_ROOT/_data/genome/human_rrna.fasta"

# --- Download and Parse Ensembl + SortMeRNA human rRNAs ---
if [ ! -f "$DB_FILE" ] || [ ! -s "$DB_FILE" ]; then
    echo "Custom human rRNA reference FASTA not found at $DB_FILE. Building database..."
    
    # 1. Extract human sequences from SortMeRNA database
    SMR_DB="$PROJECT_ROOT/_data/rrna_databases/smr_v4.3_default_db.fasta"
    if [ -f "$SMR_DB" ]; then
        echo "Extracting human sequences from SortMeRNA database..."
        awk '
            BEGIN { print_seq = 0 }
            /^>/ {
                if ($0 ~ /Homo sapiens/) {
                    print_seq = 1
                    print $0
                } else {
                    print_seq = 0
                }
                next
            }
            {
                if (print_seq == 1) {
                    print $0
                }
            }
        ' "$SMR_DB" > "$DB_FILE"
    fi

    # 2. Extract MT-rRNA, 5S, 5.8S from Ensembl
    NCRNA_FILE="$GENOME_DIR/Homo_sapiens.GRCh38.ncrna.fa.gz"
    if [ ! -f "$NCRNA_FILE" ] || [ ! -s "$NCRNA_FILE" ]; then
        echo "Downloading Ensembl 115 ncRNA FASTA..."
        curl -L -o "$NCRNA_FILE" "https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz"
    fi
    
    echo "Appending Ensembl Mt-rRNA and 5S/5.8S sequences..."
    zcat "$NCRNA_FILE" | awk '
        BEGIN { print_seq = 0 }
        /^>/ {
            if ($0 ~ /gene_biotype:rRNA/ || $0 ~ /gene_biotype:Mt_rRNA/ || $0 ~ /gene_biotype:rRNA_pseudogene/) {
                print_seq = 1
                print $0
            } else {
                print_seq = 0
            }
            next
        }
        {
            if (print_seq == 1) {
                print $0
            }
        }
    ' >> "$DB_FILE"
fi

# --- Check or Build STAR rRNA Index ---
if [ ! -f "$RRNA_INDEX/Genome" ] || [ ! -s "$RRNA_INDEX/Genome" ]; then
    echo "rRNA STAR index not found. Building index..."
    # Build with scaled parameters for a small genome (rRNA database is ~15-20MB)
    STAR \
        --runThreadN "${SLURM_CPUS_PER_TASK:-32}" \
        --runMode genomeGenerate \
        --genomeDir "$RRNA_INDEX" \
        --genomeFastaFiles "$DB_FILE" \
        --genomeSAindexNbases 7 \
        --limitGenomeGenerateRAM 60000000000 \
        --genomeChrBinNbits 11
    echo "rRNA STAR index generation complete."
fi

# Assert index exists
if [ ! -f "$RRNA_INDEX/Genome" ]; then
    echo "Assertion Error: rRNA STAR index building failed." >&2
    exit 1
fi

echo "================================================================================"
echo "   rRNA FILTER: STAR pre-mapping depletion"
echo "================================================================================"
echo "Input trimmed FASTQ  : $TRIMMED_DIR"
echo "Output filtered FASTQ : $OUTDIR"
echo "rRNA Genome index     : $RRNA_INDEX"

# --- Resource Auto-detection ---
TOTAL_CPUS=${SLURM_CPUS_PER_TASK:-32}
# STAR is highly multithreaded; run 2 samples in parallel
PARALLEL_SAMPLES=2
CPUS_PER_SAMPLE=$((TOTAL_CPUS / PARALLEL_SAMPLES))
[ "$CPUS_PER_SAMPLE" -lt 1 ] && CPUS_PER_SAMPLE=1

# RAM limit for STAR sorting
SORT_RAM=16000000000 # 16 GB per job

echo "Resources: $PARALLEL_SAMPLES concurrent samples, $CPUS_PER_SAMPLE threads each."

# --- Parse sample sheet ---
# Col 1=file1_base, Col 2=file2_base, Col 3=accession, Col 9=DataType (m6A|RNA-seq)
SAMPLE_ROWS=$(awk -F',' 'NR>1 {print $3 "\t" $1 "\t" $2 "\t" $9}' "$CSV_FILE" | grep -v '^$')
SAMPLE_COUNT=$(echo "$SAMPLE_ROWS" | wc -l)

if [ "$SAMPLE_COUNT" -eq 0 ]; then
    echo "Assertion Error: No sample rows parsed from $CSV_FILE" >&2
    exit 1
fi

# Track background job PIDs for robust parallelization control
pids=()
CURRENT_SAMPLE=0

while read -r row; do
    [ -n "$row" ] || continue
    CURRENT_SAMPLE=$((CURRENT_SAMPLE + 1))

    accession=$(echo "$row" | cut -f1)
    file1_base=$(echo "$row" | cut -f2)
    file2_base=$(echo "$row" | cut -f3)
    data_type=$(echo "$row" | cut -f4)

    # --- Only process RNA-seq (input) samples ---
    if [ "$data_type" != "RNA-seq" ]; then
        echo "[Skip] [$CURRENT_SAMPLE/$SAMPLE_COUNT] $accession (type=$data_type, rRNA filter not needed)"
        continue
    fi

    r1_in="$TRIMMED_DIR/${file1_base}_paired.fastq.gz"
    r2_in="$TRIMMED_DIR/${file2_base}_paired.fastq.gz"

    r1_out="$OUTDIR/${file1_base}_filtered.fastq.gz"
    r2_out="$OUTDIR/${file2_base}_filtered.fastq.gz"

    # Skip if both outputs already exist and are non-empty
    if [ -s "$r1_out" ] && [ -s "$r2_out" ]; then
        echo "[Skip] [$CURRENT_SAMPLE/$SAMPLE_COUNT] $accession (filtered files already exist)"
        continue
    fi

    # Skip if trimmed inputs are missing
    if [ ! -f "$r1_in" ] || [ ! -f "$r2_in" ]; then
        echo "[Warn] Trimmed input missing for $accession — skipping." >&2
        continue
    fi

    # Manage active background jobs using explicit PIDs
    while true; do
        active_pids=()
        for pid in "${pids[@]}"; do
            if kill -0 "$pid" 2>/dev/null; then
                active_pids+=("$pid")
            fi
        done
        pids=("${active_pids[@]}")
        
        if [ "${#pids[@]}" -lt "$PARALLEL_SAMPLES" ]; then
            break
        fi
        sleep 5
    done

    echo ">>> Submitting [$CURRENT_SAMPLE/$SAMPLE_COUNT]: $accession"
    SAMPLE_LOG="$LOG_DIR/${accession}.rrna_star.log"
    SAMPLE_WORK_DIR="$LOG_DIR/work_${accession}"

    (
        mkdir -p "$SAMPLE_WORK_DIR"
        
        # Align to rRNA, output only unmapped reads in Fastx format
        STAR \
            --runThreadN "$CPUS_PER_SAMPLE" \
            --genomeDir "$RRNA_INDEX" \
            --readFilesIn "$r1_in" "$r2_in" \
            --readFilesCommand gunzip -c \
            --outFileNamePrefix "$SAMPLE_WORK_DIR/${accession}_" \
            --outReadsUnmapped Fastx \
            --outSAMtype None \
            --outFilterMismatchNmax 5 \
            --outFilterMismatchNoverLmax 0.05 \
            --outFilterMultimapNmax 10 \
            --scoreDelOpen -10000 --scoreInsOpen -10000 \
            --alignIntronMax 1 \
            --seedPerWindowNmax 10 \
            --alignMatesGapMax 1000 \
            --limitBAMsortRAM "$SORT_RAM"
            
        # STAR outputs: Unmapped.out.mate1 and Unmapped.out.mate2
        mate1="$SAMPLE_WORK_DIR/${accession}_Unmapped.out.mate1"
        mate2="$SAMPLE_WORK_DIR/${accession}_Unmapped.out.mate2"
        
        if [ -s "$mate1" ] && [ -s "$mate2" ]; then
            # Compress and move unmapped reads
            gzip -c "$mate1" > "$r1_out"
            gzip -c "$mate2" > "$r2_out"
            
            # Clean up temp files
            rm -rf "$SAMPLE_WORK_DIR"
            echo "    Done: $accession (rRNA filtered successfully)"
        else
            echo "Assertion Error: STAR mapping failed or output unmapped files missing for $accession" >&2
            exit 1
        fi
    ) > "$SAMPLE_LOG" 2>&1 &
    pids+=("$!")

done <<< "$SAMPLE_ROWS"

echo "Waiting for all STAR rRNA filter jobs to complete..."
wait

# --- Output Assertions ---
echo "Running post-execution assertions..."
MISSING_OUTPUTS=0

while read -r row; do
    [ -n "$row" ] || continue
    accession=$(echo "$row" | cut -f1)
    file1_base=$(echo "$row" | cut -f2)
    file2_base=$(echo "$row" | cut -f3)
    data_type=$(echo "$row" | cut -f4)

    [ "$data_type" = "RNA-seq" ] || continue

    r1_in="$TRIMMED_DIR/${file1_base}_paired.fastq.gz"
    r2_in="$TRIMMED_DIR/${file2_base}_paired.fastq.gz"
    [ -f "$r1_in" ] && [ -f "$r2_in" ] || continue

    r1_out="$OUTDIR/${file1_base}_filtered.fastq.gz"
    r2_out="$OUTDIR/${file2_base}_filtered.fastq.gz"

    if [ ! -s "$r1_out" ] || [ ! -s "$r2_out" ]; then
        echo "Assertion Error: Filtered output missing or empty for $accession" >&2
        MISSING_OUTPUTS=$((MISSING_OUTPUTS + 1))
    fi
done <<< "$SAMPLE_ROWS"

if [ "$MISSING_OUTPUTS" -gt 0 ]; then
    echo "Assertion Error: $MISSING_OUTPUTS accession(s) failed rRNA filter validation." >&2
    exit 1
fi

echo "All outputs successfully validated."
echo "=== rRNA Filter Complete ==="
