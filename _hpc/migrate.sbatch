#!/bin/bash
#SBATCH --job-name=migrate_rna
#SBATCH --output=/mnt/22T/phongdinh/RNAseq_pipeline_depression/migrate_%j.out
#SBATCH --error=/mnt/22T/phongdinh/RNAseq_pipeline_depression/migrate_%j.err
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=12:00:00

set -e

DEST="/mnt/22T/phongdinh/RNAseq_pipeline_depression"
SRC="/mnt/22T/kienvu"

echo "=== Starting Migration to $DEST ==="

echo "1. Creating directory structure..."
mkdir -p "$DEST/_data/bam"
mkdir -p "$DEST/_data/genome"
mkdir -p "$DEST/_data/featurecounts"
mkdir -p "$DEST/_data/qc_raw"
mkdir -p "$DEST/_data/qc_trim"
mkdir -p "$DEST/_data/_metadata"
mkdir -p "$DEST/_data/RData"
mkdir -p "$DEST/01_core_alignment"
mkdir -p "$DEST/02_dge_analysis"
mkdir -p "$DEST/03_variant_calling"
mkdir -p "$DEST/utils"

echo "2. Copying metadata..."
cp "$SRC/sample_sheet.csv" "$DEST/_data/_metadata/"
cp "$SRC"/*.txt "$DEST/_data/_metadata/" 2>/dev/null || true

echo "3. Copying pipeline scripts..."
# Alignment stage
cp -r "$SRC/src/fastqc" "$SRC/src/trim" "$SRC/src/star" "$SRC/src/qualimap" "$SRC/src/count" "$DEST/01_core_alignment/" 2>/dev/null || true
# DGE analysis stage
cp "$SRC/main/"*.R "$DEST/02_dge_analysis/" 2>/dev/null || true
# Variant calling stage (03)
cp -r "$SRC/src/gatk"/* "$DEST/03_variant_calling/" 2>/dev/null || true
# Utilities
cp "$SRC/helper.py" "$SRC/read_distribution.py" "$DEST/utils/" 2>/dev/null || true
cp -r "$SRC/src/utils"/* "$DEST/utils/" 2>/dev/null || true
# Config
cp "$SRC/main/config.yaml" "$DEST/02_dge_analysis/config.yaml" 2>/dev/null || true

echo "4. Copying reference genome and annotations..."
cp "$SRC/star/ref/Homo_sapiens.GRCh38.115.gtf" "$DEST/_data/genome/"
cp "$SRC/star/ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa" "$DEST/_data/genome/"
cp "$SRC/star/ref/Homo_sapiens.GRCh38.115.bed" "$DEST/_data/genome/"
cp "$SRC/star/ref/genes.genePred" "$DEST/_data/genome/"
cp "$SRC/star/ref"/smr_*.fasta "$DEST/_data/genome/" 2>/dev/null || true

echo "5. Copying pre-computed feature counts..."
cp -r "$SRC/counts"/* "$DEST/_data/featurecounts/" 2>/dev/null || true

echo "6. Copying FastQC HTML reports..."
cp -r "$SRC/QC_1"/* "$DEST/_data/qc_raw/" 2>/dev/null || true
cp -r "$SRC/QC_2"/* "$DEST/_data/qc_trim/" 2>/dev/null || true

echo "7. Copying aligned BAM files (~582 GB)..."
cp -r "$SRC/star/aligned_bam"/*.bam "$DEST/_data/bam/"

echo "=== Migration Completed Successfully ==="
