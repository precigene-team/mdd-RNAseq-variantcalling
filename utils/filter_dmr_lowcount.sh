#!/bin/bash
# Filter DMR CSVs to approximate paper's lowcount=15 threshold.
#
# TRESS lowcount filters on TOTAL count across all samples.
# baseMean is the per-sample average, so:
#   threshold = lowcount / N_samples
# Male: 20 subjects × 2 fractions (IP+Input) = 40 samples → 15/40 = 0.375
# Female: 19 subjects × 2 fractions = 38 samples → 15/38 ≈ 0.395

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BASE="$SCRIPT_DIR/results/01_core_alignment/tables"
OUTDIR="$SCRIPT_DIR/results/utils/tables"

mkdir -p "$OUTDIR"

for sex in male female; do
    infile="$BASE/10_dmr_${sex}_full.csv"
    outfile="$OUTDIR/10_dmr_${sex}_filtered.csv"

    [ "$sex" = "male" ] && n=40 || n=38
    thresh=$(awk "BEGIN {printf \"%.6f\", 15/$n}")

    total=$(tail -n +2 "$infile" | wc -l)

    # Copy header + filter
    awk -F',' -v t="$thresh" \
        'NR==1{print; next} $5+0>=t && $9+0!="NA" && $9+0<0.05 && sqrt($6*$6)>=0.263' \
        "$infile" > "$outfile"

    sig=$(awk -F',' 'NR>1 && $9+0!="NA" && $9+0<0.05 && sqrt($6*$6)>=0.263' "$infile" | wc -l)
    lc15=$(tail -n +2 "$outfile" | wc -l)

    printf "%-7s N=%-2s lowcount=15 (baseMean≥%.3f): total %5d → sig %5d → after filter: %5d\n" \
        "$sex" "$n" "$thresh" "$total" "$sig" "$lc15"
done
