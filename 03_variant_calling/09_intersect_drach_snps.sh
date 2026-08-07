#!/bin/bash
# [[03_variant_calling/09_intersect_drach_snps.sh]]
# Goal: Intersect GATK cohort SNPs with DRACH motif coordinates in TRESS peaks
#       to identify candidate motif-disrupting variants.

# Steps:
#   1. Validate inputs (DRACH BED, merged cohort VCF).
#   2. Convert cohort VCF to SNP-only BED for intersection.
#   3. Intersect SNP BED with DRACH BED using bedtools.
#   4. For each overlapping SNP, determine which DRACH position it hits
#      and annotate whether the ALT allele breaks the DRACH pattern.
#   5. Save annotated disrupting variants table and count summary.

set -e
set -o pipefail

export PROJECT_ROOT="${PROJECT_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"

if [ -z "$CONDA_PREFIX" ]; then
    echo "Assertion Error: Mamba environment 'mdd' is not active." >&2
    exit 1
fi

# --- Paths ---
DRACH_BED="$PROJECT_ROOT/_data/peaks/drach_sites.bed"
VCF_DIR="$PROJECT_ROOT/_data/variant/filtered_vcf"
OUT_TABLES="$PROJECT_ROOT/results/03_variant_calling/tables"
OUT_FIGURES="$PROJECT_ROOT/results/03_variant_calling/figures"
LOG_DIR="$PROJECT_ROOT/_data/logs/drach"
TMP_DIR="$PROJECT_ROOT/_data/logs/drach/tmp"
DRACH_ANNO="$OUT_TABLES/09_drach_disrupting_snps.csv"

mkdir -p "$LOG_DIR" "$TMP_DIR" "$OUT_TABLES" "$OUT_FIGURES"

# --- Input Assertions ---
if [ ! -f "$DRACH_BED" ]; then
    echo "Assertion Error: DRACH sites BED not found: $DRACH_BED" >&2
    echo "Run 08_find_drach_sites.py first." >&2
    exit 1
fi

DRACH_COUNT=$(grep -v "^#" "$DRACH_BED" | wc -l)
echo ">>> DRACH sites loaded: $DRACH_COUNT"

# --- Step 1: Build per-subject SNP BEDs and merge ---
echo ">>> Step 1: Building per-subject SNP BED files from filtered VCFs..."
MERGED_SNP_BED="$TMP_DIR/all_snps.bed"
> "$MERGED_SNP_BED"

for vcf in "$VCF_DIR"/*.filtered.vcf.gz; do
    [ -f "$vcf" ] || continue
    subj=$(basename "$vcf" .filtered.vcf.gz)
    # Extract PASS-only, biallelic, heterozygous (0/1) SNPs → BED (chr, pos-1, pos, subj, REF, ALT)
    bcftools view -f PASS,. -v snps -m2 -M2 -g het "$vcf" | \
    bcftools query -f '%CHROM\t%POS0\t%POS\t'"$subj"'\t%REF\t%ALT\n' >> "$MERGED_SNP_BED"
done

SNP_COUNT=$(wc -l < "$MERGED_SNP_BED")
echo "    Total per-subject SNP sites: $SNP_COUNT"

if [ "$SNP_COUNT" -eq 0 ]; then
    echo "Assertion Error: No SNPs found in filtered VCFs. Run 05_variant_filtration.sh first." >&2
    exit 1
fi

# --- Step 2: Intersect SNPs with DRACH sites ---
echo ">>> Step 2: Intersecting SNPs with DRACH sites..."
INTERSECT_BED="$TMP_DIR/snp_drach_intersect.bed"

bedtools intersect \
    -a "$MERGED_SNP_BED" \
    -b <(grep -v "^#" "$DRACH_BED") \
    -wa -wb \
    > "$INTERSECT_BED"

OVERLAP_COUNT=$(wc -l < "$INTERSECT_BED")
echo "    SNP × DRACH overlaps: $OVERLAP_COUNT"

if [ "$OVERLAP_COUNT" -eq 0 ]; then
    echo "[WARNING] No SNPs overlap any DRACH site. Output will be empty."
    echo "chrom,snp_pos,subject,ref,alt,drach_start,drach_end,drach_motif,strand,drach_pos_hit,disrupts_drach" \
        > "$DRACH_ANNO"
    exit 0
fi

# --- Step 3: Annotate disruption type ---
echo ">>> Step 3: Annotating motif disruption type (delegating to 10_annotate_drach_disruption.py)..."
python3 "$PROJECT_ROOT/03_variant_calling/09_annotate_drach_disruption.py"


# --- Output Assertion ---
if [ ! -s "$DRACH_ANNO" ]; then
    echo "Assertion Error: Final annotated DRACH disruption table is empty or missing." >&2
    exit 1
fi

TOTAL=$(tail -n +2 "$DRACH_ANNO" | wc -l)
DISRUPTING=$(tail -n +2 "$DRACH_ANNO" | awk -F',' '$11=="DISRUPTS_DRACH"' | wc -l)
echo "================================================================================"
echo "   DRACH SNP Intersection Complete"
echo "   Total overlapping SNPs: $TOTAL"
echo "   Motif-disrupting variants: $DISRUPTING"
echo "   Output: $DRACH_ANNO"
echo "================================================================================"
