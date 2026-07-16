#!/bin/bash
# [[utils/run_qc03_audit.sh]]
# Goal: Wrapper to run the full cohort-wide QC-03 audit (gene-level vs peak-level)
#       across all contrasts that have significant results.

# Steps:
#   1. Resolve project root and key paths.
#   2. Assert required inputs exist (GTF, peaks dir, Python script).
#   3. Run cohort audit for each contrast that produced significant DMGs/DEGs.
#   4. Print a consolidated summary of all audit outputs.

set -e
set -o pipefail

export PROJECT_ROOT="${PROJECT_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"

# --- Paths ---
SCRIPT="$PROJECT_ROOT/utils/qc03_audit_gene_vs_peak.py"
GTF="$PROJECT_ROOT/_data/genome/Homo_sapiens.GRCh38.115.gtf"
PEAKS_DIR="$PROJECT_ROOT/_data/peaks"
TABLES="$PROJECT_ROOT/results/02_dge_analysis/tables"
OUT_DIR="$PROJECT_ROOT/results/utils/tables"
mkdir -p "$OUT_DIR"

PYTHON="${CONDA_PREFIX}/bin/python3"
if [ ! -x "$PYTHON" ]; then
    PYTHON=$(which python3)
fi

# --- Assertions ---
if [ ! -f "$SCRIPT" ]; then
    echo "ERROR: audit script not found at $SCRIPT" >&2; exit 1
fi
if [ ! -f "$GTF" ]; then
    echo "ERROR: GTF not found at $GTF" >&2; exit 1
fi
if [ ! -d "$PEAKS_DIR" ]; then
    echo "ERROR: Peaks directory not found at $PEAKS_DIR" >&2; exit 1
fi

echo "================================================================================"
echo "   QC-03 Cohort Audit — all contrasts"
echo "================================================================================"
echo "  Python  : $PYTHON"
echo "  Script  : $SCRIPT"
echo "  GTF     : $GTF"
echo "  Peaks   : $PEAKS_DIR"
echo "  Output  : $OUT_DIR"
echo "================================================================================"

# Helper: run one audit and report whether the DMG file is empty
run_audit() {
    local label="$1"
    local dmg_male="$2"
    local dmg_female="$3"
    local out_base="$4"

    echo ""
    echo ">>> Contrast: $label"

    # Count DMGs to give context
    local n_male=$(( $(wc -l < "$dmg_male") - 1 ))
    local n_female=$(( $(wc -l < "$dmg_female") - 1 ))
    echo "    Male DMGs   : $n_male  ($dmg_male)"
    echo "    Female DMGs : $n_female  ($dmg_female)"

    if [ "$n_male" -eq 0 ] && [ "$n_female" -eq 0 ]; then
        echo "    [SKIP] Both DMG files are empty — no sensitivity check possible."
        echo "           Exonic mapping rates will still be calculated."
    fi

    "$PYTHON" "$SCRIPT" \
        --peaks-dir "$PEAKS_DIR" \
        --dmg-male  "$dmg_male" \
        --dmg-female "$dmg_female" \
        --gtf "$GTF" \
        --out "${out_base}_summary.txt"

    echo "    Saved: ${out_base}_summary.txt"
    echo "    Saved: ${out_base}_summary.csv"
}

# --- Run all contrasts ---

# 1. Methylation interaction (IP/Input relative methylation — primary TRESS target)
run_audit \
    "Methylation Interaction (Female vs Male)" \
    "$TABLES/04_meth_male_interaction_sig.csv" \
    "$TABLES/04_meth_female_interaction_sig.csv" \
    "$OUT_DIR/qc03_meth_interaction"

# 2. IP fraction baseline DGE (m6A-enriched transcript abundance)
run_audit \
    "IP Baseline DGE" \
    "$TABLES/04_deg_male_ip_sig.csv" \
    "$TABLES/04_deg_female_ip_sig.csv" \
    "$OUT_DIR/qc03_ip_baseline"

# 3. Input fraction baseline DGE (steady-state transcription)
run_audit \
    "Input Baseline DGE" \
    "$TABLES/04_deg_male_input_sig.csv" \
    "$TABLES/04_deg_female_input_sig.csv" \
    "$OUT_DIR/qc03_input_baseline"

# 4. Sex interaction on Input (opposing DGE across sexes)
run_audit \
    "Sex Interaction (Input)" \
    "$TABLES/04_deg_sex_interaction_input_sig.csv" \
    "$TABLES/04_deg_sex_interaction_input_sig.csv" \
    "$OUT_DIR/qc03_sex_interaction_input"

# 5. Methylation sex interaction (opposing methylation shifts across sexes)
run_audit \
    "Methylation Sex Interaction" \
    "$TABLES/04_meth_sex_interaction_sig.csv" \
    "$TABLES/04_meth_sex_interaction_sig.csv" \
    "$OUT_DIR/qc03_meth_sex_interaction"

echo ""
echo "================================================================================"
echo "   All audits complete."
echo "   Summary files written to: $OUT_DIR/qc03_*.txt"
echo "   Detailed CSVs written to:  $OUT_DIR/qc03_*.csv"
echo "================================================================================"
