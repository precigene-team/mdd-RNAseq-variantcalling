# [[04_asm/03_intersect_drach_genotypes.sh]]
# Goal: Intersect all PASS biallelic genotypes with DRACH motif sites and annotate strand-aware disruption.
#
# Steps:
#   1. Validate inputs (DRACH BED from 01, 39 cohort VCFs).
#   2. Emit one BED row per donor x variant for all PASS biallelic SNPs.
#   3. bedtools intersect against DRACH sites (-wa -wb).
#   4. Gate on 02, then delegate strand-aware annotation to 04.
#
# Notes: no -g het. Under the gene-dosage model both 0/1 and 1/1 are valid m6A-SNPs.

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
OUT_TABLES="$PROJECT_ROOT/results/04_asm/tables"
OUT_FIGURES="$PROJECT_ROOT/results/04_asm/figures"
LOG_DIR="$PROJECT_ROOT/_data/logs/drach"
TMP_DIR="$PROJECT_ROOT/_data/logs/drach/tmp"
DRACH_ANNO="$OUT_TABLES/04_drach_variant_overlaps.csv"

mkdir -p "$LOG_DIR" "$TMP_DIR" "$OUT_TABLES" "$OUT_FIGURES"

# --- Input Assertions ---
if [ ! -f "$DRACH_BED" ]; then
    echo "Assertion Error: DRACH sites BED not found: $DRACH_BED" >&2
    echo "Run 01_find_drach_sites.py first." >&2
    exit 1
fi

DRACH_COUNT=$(grep -vc "^#" "$DRACH_BED")
echo ">>> DRACH sites loaded: $DRACH_COUNT"

VCF_COUNT=$(find "$VCF_DIR" -name "*_input.filtered.vcf.gz" | wc -l)
if [ "$VCF_COUNT" -ne 39 ]; then
    echo "Assertion Error: expected 39 Input VCFs, found $VCF_COUNT" >&2
    exit 1
fi
echo ">>> Cohort VCFs: $VCF_COUNT donors"

# --- Step 1: Build donor x variant BED ---
echo ">>> Step 1: Building donor x variant BED from all PASS biallelic SNPs..."
MERGED_SNP_BED="$TMP_DIR/all_snps.bed"
> "$MERGED_SNP_BED"

for vcf in "$VCF_DIR"/*_input.filtered.vcf.gz; do
    [ -f "$vcf" ] || continue
    donor=$(basename "$vcf" _input.filtered.vcf.gz)   # donor-level key
    # PASS + biallelic SNPs, ALL genotypes (no -g het)
    bcftools view -f PASS -v snps -m2 -M2 "$vcf" | \
    bcftools query -f '%CHROM\t%POS0\t%POS\t'"$donor"'\t%REF\t%ALT\n' >> "$MERGED_SNP_BED"
done

SNP_COUNT=$(wc -l < "$MERGED_SNP_BED")
echo "    Total donor x variant rows: $SNP_COUNT"

if [ "$SNP_COUNT" -eq 0 ]; then
    echo "Assertion Error: No SNPs found in filtered VCFs. Run 05_variant_filtration.sh first." >&2
    exit 1
fi

# --- Step 2: Intersect with DRACH sites ---
echo ">>> Step 2: Intersecting variants with DRACH sites..."
INTERSECT_BED="$TMP_DIR/snp_drach_intersect.bed"

bedtools intersect \
    -a "$MERGED_SNP_BED" \
    -b <(grep -v "^#" "$DRACH_BED") \
    -wa -wb \
    > "$INTERSECT_BED"

OVERLAP_COUNT=$(wc -l < "$INTERSECT_BED")
echo "    donor x variant x motif overlaps: $OVERLAP_COUNT"

if [ "$OVERLAP_COUNT" -eq 0 ]; then
    echo "Assertion Error: no variant overlaps any DRACH motif." >&2
    exit 1
fi

# --- Step 3: Strand-aware disruption annotation ---
echo ">>> Step 3: Annotating strand-aware motif disruption..."
python3 "$PROJECT_ROOT/04_asm/02_verify_drach_reference.py"
Rscript "$PROJECT_ROOT/04_asm/04_annotate_strand_aware_disruption.R"

# --- Output Assertion ---
if [ ! -s "$DRACH_ANNO" ]; then
    echo "Assertion Error: Final annotated DRACH overlap table is empty or missing." >&2
    exit 1
fi

echo "================================================================================"
echo "   DRACH genotype intersection complete"
echo "   Output: $DRACH_ANNO"
echo "================================================================================"