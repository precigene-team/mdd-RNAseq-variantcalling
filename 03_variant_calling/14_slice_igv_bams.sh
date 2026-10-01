#!/bin/bash
# [[03_variant_calling/14_slice_igv_bams.sh]]
# Goal: Extract +/-100bp BAM slices around ZNF419 and CCDC30 candidate SNPs for IGV visualization.

#SBATCH --job-name=slice_igv_bams
#SBATCH --output=_data/logs/igv_slices/%x_%j.log
#SBATCH --error=_data/logs/igv_slices/%x_%j.err
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00

# Steps:
#   1. Set up dynamic project paths and environment.
#   2. Define target loci coordinates with +/-100bp forward and backward margins.
#   3. Resolve input BAM files for 2 shared MDD carrier subjects (male_mdd_10, female_mdd_10).
#   4. Extract coordinate-sorted, indexed BAM slices via samtools.
#   5. Generate complementary BED annotation files for IGV locus markers.

set -euo pipefail

# --- Step 1: Dynamic Pathing & Setup ---
if [ -n "${SLURM_SUBMIT_DIR:-}" ]; then
    export PROJECT_ROOT="$SLURM_SUBMIT_DIR"
else
    export PROJECT_ROOT="${PROJECT_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
fi

cd "$PROJECT_ROOT"

LOG_DIR="_data/logs/igv_slices"
OUT_DIR="results/03_variant_calling/tracks/igv_slices"
mkdir -p "$LOG_DIR" "$OUT_DIR"

# 2 Common MDD carrier subjects sharing both ZNF419 and CCDC30 SNPs
MDD_SUBJECTS=("male_mdd_10_input" "female_mdd_10_input")

# --- Step 2: Define Genomic Loci (+/-100 bp window) ---
# ZNF419: chr19:57492978 (G>C, rs2074076) -> window: [57492878, 57493078]
# CCDC30: chr1:42497129 (A>T, rs3748841)  -> window: [42497029, 42497229]
declare -A LOCI_CHR
declare -A LOCI_START
declare -A LOCI_END
declare -A LOCI_POS
declare -A LOCI_REF_ALT

LOCI_CHR["ZNF419"]="19"
LOCI_POS["ZNF419"]=57492978
LOCI_START["ZNF419"]=$((57492978 - 100))
LOCI_END["ZNF419"]=$((57492978 + 100))
LOCI_REF_ALT["ZNF419"]="G_to_C_rs2074076"

LOCI_CHR["CCDC30"]="1"
LOCI_POS["CCDC30"]=42497129
LOCI_START["CCDC30"]=$((42497129 - 100))
LOCI_END["CCDC30"]=$((42497129 + 100))
LOCI_REF_ALT["CCDC30"]="A_to_T_rs3748841"

# Create locus marker BED for IGV
BED_FILE="$OUT_DIR/candidate_snps_pm100bp.bed"
> "$BED_FILE"

echo "========================================================================"
echo "    EXTRACTING IGV BAM SLICES FOR ZNF419 & CCDC30 (+/-100 bp MARGIN)    "
echo "========================================================================"

# --- Step 3: BAM Slicing ---
find_bam() {
    local sample="$1"
    local candidates=(
        "_data/variant/split_n_cigar/${sample}.split.bam"
        "_data/variant/duplicates/${sample}.dedup.bam"
        "_data/variant/readgroups/${sample}.rg.bam"
        "_data/bam/${sample}/${sample}_Aligned.sortedByCoord.out.bam"
    )
    for b in "${candidates[@]}"; do
        if [ -s "$b" ]; then
            echo "$b"
            return 0
        fi
    done
    return 1
}

for gene in "ZNF419" "CCDC30"; do
    c_raw="${LOCI_CHR[$gene]}"
    p="${LOCI_POS[$gene]}"
    s="${LOCI_START[$gene]}"
    e="${LOCI_END[$gene]}"
    ann="${gene}_${LOCI_REF_ALT[$gene]}"

    echo -e "${c_raw}\t$((p - 1))\t${p}\t${ann}" >> "$BED_FILE"

    echo ""
    echo ">>> Locus: ${gene} (Chr ${c_raw}:${p}, Window: ${s}-${e})"

    for sample in "${MDD_SUBJECTS[@]}"; do
        in_bam=$(find_bam "$sample" || true)

        if [ -z "$in_bam" ]; then
            echo "  [WARNING] No BAM found for sample ${sample} under _data/. Skipping."
            continue
        fi

        # Detect chromosome naming convention in input BAM header
        if samtools idxstats "$in_bam" 2>/dev/null | grep -q "^chr${c_raw}\b"; then
            region="chr${c_raw}:${s}-${e}"
        else
            region="${c_raw}:${s}-${e}"
        fi

        out_slice="$OUT_DIR/${sample}_${gene}_slice_pm100bp.bam"

        echo "  [Processing] ${sample} -> ${region}"
        samtools view -b -h "$in_bam" "$region" > "$out_slice"
        samtools index "$out_slice"
        
        n_reads=$(samtools view -c "$out_slice" || echo 0)
        echo "  [SUCCESS] Extracted ${n_reads} reads -> ${out_slice}"
    done
done

echo ""
echo "========================================================================"
echo "[COMPLETE] Slices and BED generated in: ${OUT_DIR}"
echo "To view in IGV Desktop or igv-reports:"
echo "  Genome: GRCh38"
echo "  Loci:   19:57492878-57493078 (ZNF419) and 1:42497029-42497229 (CCDC30)"
echo "  BED:    ${BED_FILE}"
echo "========================================================================"
