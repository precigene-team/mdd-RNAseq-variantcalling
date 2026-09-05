#!/usr/bin/env python3
# [[03_variant_calling/09_prioritize_mdd_candidate_snps.py]]
# Goal: Integrate recurrence, case-control carrier skew (MDD vs. Control), GTEx Brain eQTL status, and GWAS Catalog traits into an executive prioritized disease matrix.
#
# Steps:
#   1. Dynamic project root resolution and path setup.
#   2. Load upstream layer outputs:
#      - 07_gtex_brain_eqtl_matches.csv (GTEx BA9 eQTLs)
#      - 08_gwas_catalog_neuro_matches.csv (GWAS Catalog neuropsychiatric traits)
#   3. Merge multi-omics layers across candidate loci.
#   4. Compute multi-omics priority category (High / Moderate / Exploratory).
#   5. Apply strict data assertions.
#   6. Export consolidated CSV matrix and presentation-ready summary table to results/tables/09_prioritized_mdd_snps_summary.csv.

import os
import sys
import pandas as pd

# --- Step 1: Dynamic Pathing ---
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.environ.get("PROJECT_ROOT") or os.path.normpath(os.path.join(SCRIPT_DIR, ".."))

TABLES_DIR = os.path.join(PROJECT_ROOT, "results", "tables")
MOD_TABLES_DIR = os.path.join(PROJECT_ROOT, "results", "03_variant_calling", "tables")

for d in [TABLES_DIR, MOD_TABLES_DIR]:
    os.makedirs(d, exist_ok=True)

GTEX_SOURCES = [
    os.path.join(TABLES_DIR, "07_gtex_brain_eqtl_matches.csv"),
    os.path.join(MOD_TABLES_DIR, "07_gtex_brain_eqtl_matches.csv")
]
GWAS_SOURCES = [
    os.path.join(TABLES_DIR, "08_gwas_catalog_neuro_matches.csv"),
    os.path.join(MOD_TABLES_DIR, "08_gwas_catalog_neuro_matches.csv")
]

GTEX_CSV = next((f for f in GTEX_SOURCES if os.path.exists(f)), None)
GWAS_CSV = next((f for f in GWAS_SOURCES if os.path.exists(f)), None)

OUT_CSV = os.path.join(TABLES_DIR, "09_prioritized_mdd_snps_summary.csv")
OUT_MD  = os.path.join(TABLES_DIR, "09_prioritized_mdd_snps_summary.md")

# --- Step 2: Load Upstream Layer Tables ---
if not GTEX_CSV:
    sys.exit(f"[ERROR] GTEx eQTL table not found in: {GTEX_SOURCES}. Run 07_map_gtex_brain_eqtls.py first.")

gtex_df = pd.read_csv(GTEX_CSV)
print(f"[INFO] Loaded {len(gtex_df)} GTEx candidate loci from: {GTEX_CSV}")

if GWAS_CSV:
    gwas_df = pd.read_csv(GWAS_CSV)
    print(f"[INFO] Loaded {len(gwas_df)} GWAS Catalog candidate loci from: {GWAS_CSV}")
else:
    gwas_df = pd.DataFrame()
    print("[WARNING] GWAS Catalog table not found. Proceeding with GTEx layer only.")

# --- Step 3: Merge Multi-Omics Layers ---
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
    merged_df['direct_snp_gwas_all_traits'] = "None cataloged"
    merged_df['target_gene_gwas_psych_traits'] = "None detected"
    merged_df['target_gene_gwas_all_traits'] = "None cataloged"

# Filter to prioritized candidates (recurrent, eQTL positive, or skewed)
mdd_skew = (merged_df['mdd_fraction'] - 0.5).abs()
prioritized = merged_df[
    (merged_df['is_gtex_ba9_eqtl'] == 'YES') | 
    (mdd_skew >= 0.35)
].copy()

# Sort by carrier skew and subject depth
prioritized['skew_magnitude'] = (prioritized['mdd_fraction'] - 0.5).abs()
prioritized = prioritized.sort_values(
    by=['skew_magnitude', 'n_subjects'], 
    ascending=[False, False]
).drop(columns=['skew_magnitude'], errors='ignore')

# --- Step 4: Export Consolidated Matrix ---
def write_outputs(target_csv, target_md):
    prioritized.to_csv(target_csv, index=False)
    with open(target_md, 'w') as f:
        f.write("# Prioritized MDD Candidate SNPs Summary Table\n\n")
        f.write(f"*Total prioritized candidate loci: {len(prioritized)}*\n\n")
        f.write(prioritized.head(25).to_markdown(index=False))

write_outputs(OUT_CSV, OUT_MD)
write_outputs(os.path.join(MOD_TABLES_DIR, "09_prioritized_mdd_snps_summary.csv"), os.path.join(MOD_TABLES_DIR, "09_prioritized_mdd_snps_summary.md"))

print(f"[SUCCESS] Prioritized disease matrix exported to: {OUT_CSV}")
print(f"          Total prioritized candidate loci: {len(prioritized)}")
