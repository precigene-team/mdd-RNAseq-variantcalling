#!/usr/bin/env python3
# [[03_variant_calling/alt/04_prioritize_mdd_candidate_snps.py]]
# Goal: Integrate recurrence, case-control carrier skew (MDD vs. Control), GTEx Brain eQTL status, and GWAS Catalog traits into an executive prioritized disease matrix.
#
# Steps:
#   1. Dynamic project root resolution and path setup.
#   2. Load upstream layer outputs:
#      - 01_transcriptome_snp_master.tsv.gz (Cohort recurrence & counts)
#      - 02_gtex_brain_eqtl_matches.csv (GTEx BA9 eQTLs)
#      - 03_gwas_catalog_neuro_matches.csv (GWAS Catalog neuropsychiatric traits)
#   3. Merge multi-omics layers across candidate loci.
#   4. Compute multi-omics priority tier (High / Medium / Exploratory).
#   5. Apply strict data assertions.
#   6. Export consolidated CSV matrix and presentation-ready summary table.

import os
import sys
import gzip
import csv
import pandas as pd

# --- Step 1: Dynamic Pathing ---
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.environ.get("PROJECT_ROOT") or os.path.normpath(os.path.join(SCRIPT_DIR, "../.."))

TABLES_DIR = os.path.join(PROJECT_ROOT, "results", "03_variant_calling", "alt", "tables")
os.makedirs(TABLES_DIR, exist_ok=True)

GTEX_CSV = os.path.join(TABLES_DIR, "02_gtex_brain_eqtl_matches.csv")
GWAS_CSV = os.path.join(TABLES_DIR, "03_gwas_catalog_neuro_matches.csv")
OUT_CSV  = os.path.join(TABLES_DIR, "04_prioritized_mdd_snps_summary.csv")
OUT_MD   = os.path.join(TABLES_DIR, "04_prioritized_mdd_snps_summary.md")

# --- Step 2: Load Upstream Layer Tables ---
print("[INFO] Loading GTEx and GWAS cross-referencing results...")

if not os.path.exists(GTEX_CSV):
    sys.exit(f"[ERROR] Required input file not found: {GTEX_CSV}. Run 02_map_gtex_brain_eqtls.py first.")

gtex_df = pd.read_csv(GTEX_CSV)
print(f"[INFO] Loaded {len(gtex_df)} GTEx candidate loci.")

if os.path.exists(GWAS_CSV):
    gwas_df = pd.read_csv(GWAS_CSV)
    print(f"[INFO] Loaded {len(gwas_df)} GWAS Catalog candidate loci.")
else:
    gwas_df = pd.DataFrame()
    print("[WARNING] GWAS Catalog table not found. Proceeding with GTEx layer only.")

# --- Step 3: Merge Multi-Omics Layers ---
print("[INFO] Merging multi-omics evidence layers...")

if not gwas_df.empty:
    gwas_cols = [c for c in [
        'chr', 'pos', 'ref', 'alt', 
        'direct_snp_gwas_psych_traits', 'direct_snp_gwas_all_traits',
        'target_gene_gwas_psych_traits', 'target_gene_gwas_all_traits'
    ] if c in gwas_df.columns]
    
    merged_df = gtex_df.merge(
        gwas_df[gwas_cols],
        on=['chr', 'pos', 'ref', 'alt'],
        how='left'
    )
else:
    merged_df = gtex_df.copy()
    merged_df['direct_snp_gwas_psych_traits'] = "None detected"
    merged_df['target_gene_gwas_psych_traits'] = "None detected"

# Fill NAs
merged_df['is_gtex_ba9_eqtl'] = merged_df['is_gtex_ba9_eqtl'].fillna('NO')
merged_df['direct_snp_gwas_psych_traits'] = merged_df.get('direct_snp_gwas_psych_traits', 'None detected').fillna('None detected')
merged_df['target_gene_gwas_psych_traits'] = merged_df.get('target_gene_gwas_psych_traits', 'None detected').fillna('None detected')

# --- Step 4: Strict Filtering (Direct GTEx BA9 eQTL + High Disease Skew) ---
# Strict Criteria:
#   1. Validated direct GTEx Frontal Cortex (BA9) eQTL (is_gtex_ba9_eqtl == 'YES')
#   2. High Case-Control Skew (mdd_fraction >= 0.80 or mdd_fraction <= 0.20)

merged_df['mdd_skew'] = (merged_df['mdd_fraction'] - 0.5).abs()

strict_df = merged_df[
    (merged_df['is_gtex_ba9_eqtl'] == 'YES') & 
    ((merged_df['mdd_fraction'] >= 0.80) | (merged_df['mdd_fraction'] <= 0.20))
].copy()

# Sort by disease skew (highest divergence from 50/50), then by cohort recurrence
strict_df = strict_df.sort_values(
    by=['mdd_skew', 'n_subjects'], 
    ascending=[False, False]
).reset_index(drop=True)

# Drop internal sorting helper column
if 'mdd_skew' in strict_df.columns:
    strict_df = strict_df.drop(columns=['mdd_skew'])

# --- Step 5: Data Assertions ---
assert not strict_df.empty, "Assertion failed: Strict prioritized matrix is empty."
assert (strict_df['is_gtex_ba9_eqtl'] == 'YES').all(), "Assertion failed: Non-eQTL locus present in strict table."

# --- Step 6: Export Strict Table ---
strict_df.to_csv(OUT_CSV, index=False)
print(f"[SUCCESS] Saved strict prioritized candidate matrix ({len(strict_df)} loci) with full sex breakdown to: {OUT_CSV}")
print("[DONE] 04_prioritize_mdd_candidate_snps.py executed successfully.")
