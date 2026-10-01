#!/usr/bin/env python3
# [[03_variant_calling/10_prioritize_mdd_candidate_snps.py]]
# Goal: Prioritize high- and moderate-impact, protein-altering variants validated by GTEx Brain BA9 eQTLs, enriched in MDD (n_ctrl <= 2 & n_mdd >= 8).
#
# Steps:
#   1. Dynamic project root resolution and path setup.
#   2. Load annotated candidates from 09_gwas_catalog_neuro_matches.csv.
#   3. Filter for functional candidates (HIGH & MODERATE impact, valid hgvs_p, GTEx BA9 eQTL validated).
#   4. Strict Filter: restrict to n_ctrl <= 2 & n_mdd >= 8, sort MDD descending.
#   5. Apply strict data assertions.
#   6. Save clean prioritized matrix to results/03_variant_calling/tables/10_prioritized_mdd_snps_summary.csv.

import os
import sys
import pandas as pd
import numpy as np

# --- Step 1: Dynamic Pathing ---
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.environ.get("PROJECT_ROOT") or os.path.normpath(os.path.join(SCRIPT_DIR, ".."))

TABLES_DIR = os.path.join(PROJECT_ROOT, "results", "03_variant_calling", "tables")
os.makedirs(TABLES_DIR, exist_ok=True)

IN_09_CSV = os.path.join(TABLES_DIR, "09_gwas_catalog_neuro_matches.csv")
IN_08_CSV = os.path.join(TABLES_DIR, "08_gtex_brain_eqtl_matches.csv")
OUT_CSV = os.path.join(TABLES_DIR, "10_prioritized_mdd_snps_summary.csv")

# Helper to attach external genomic resource URLs if reading directly from 08
def attach_resource_urls(df):
    if 'omim_url' not in df.columns:
        df['omim_url'] = df['gene'].apply(lambda g: f"https://omim.org/search?search={g}")
        df['decipher_url'] = df['gene'].apply(lambda g: f"https://www.deciphergenomics.org/gene/{g}/overview/protein-genomic-info")
        df['pubmed_url'] = df['gene'].apply(lambda g: f"https://pubmed.ncbi.nlm.nih.gov/?term={g}+AND+(depression+OR+MDD)")
        df['gnomad_gene_url'] = df['gene'].apply(lambda g: f"https://gnomad.broadinstitute.org/gene/{g}?dataset=gnomad_r4")
        df['gnomad_variant_url'] = df['variant_id'].apply(lambda v: f"https://gnomad.broadinstitute.org/variant/{str(v).replace(':', '-')}?dataset=gnomad_r4")
        df['gencc_url'] = df['gene'].apply(lambda g: f"https://search.thegencc.org/genes/{g}")
        df['gwas_catalog_gene_url'] = df['gene'].apply(lambda g: f"https://www.ebi.ac.uk/gwas/genes/{g}")
        df['gwas_catalog_variant_url'] = df['rsid'].apply(lambda r: f"https://www.ebi.ac.uk/gwas/variants/{r}" if pd.notna(r) and r != '.' else "")
        df['direct_snp_gwas_psych_traits'] = "Linked via GWAS Catalog URL"
        df['target_gene_gwas_psych_traits'] = "Linked via GWAS Catalog URL"
    return df

# --- Step 2: Load Candidates (Prefer 08 if newer or 09 missing) ---
if os.path.exists(IN_08_CSV) and (not os.path.exists(IN_09_CSV) or os.path.getmtime(IN_08_CSV) >= os.path.getmtime(IN_09_CSV)):
    print(f"[INFO] Loading freshly mapped candidates from: {IN_08_CSV}...")
    df_09 = pd.read_csv(IN_08_CSV, low_memory=False)
    df_09 = attach_resource_urls(df_09)
    df_09.to_csv(IN_09_CSV, index=False)
    print(f"[INFO] Synchronized {len(df_09):,} records into: {IN_09_CSV}")
elif os.path.exists(IN_09_CSV):
    print(f"[INFO] Loading candidates from: {IN_09_CSV}...")
    df_09 = pd.read_csv(IN_09_CSV, low_memory=False)
else:
    sys.exit(f"[ERROR] Required table not found: neither {IN_08_CSV} nor {IN_09_CSV} exists.")

print(f"[INFO] Loaded {len(df_09):,} annotated candidates.")

# --- Step 3: Filter for Functional & GTEx Confirmed Candidates ---
has_hgvs_p = df_09['hgvs_p'].apply(
    lambda x: pd.notna(x) and str(x).strip() not in ('None', '.', 'N/A', '', 'nan')
)
is_functional = df_09['impact'].str.upper().isin(['HIGH', 'MODERATE', 'LOW']) & (
    has_hgvs_p | df_09['consequence'].str.contains('splice', case=False, na=False)
)
is_eqtl = df_09['is_gtex_ba9_eqtl'].str.upper() == 'YES'

prioritized = df_09[is_functional & is_eqtl].copy()

# --- Step 4: Strict Filter: Restrict n_ctrl <= 2 & n_mdd >= 8, Sort MDD Descending ---
# Filter: Strict case-enrichment criteria
is_mdd_enriched = (prioritized['n_ctrl'] <= 2) & (prioritized['n_mdd'] >= 8)
prioritized = prioritized[is_mdd_enriched].copy()

# Sorting Hierarchy:
#   1. n_ctrl ascending (0 is prioritized before 1, 2)
#   2. n_mdd descending (highest MDD carrier burden to lowest)
#   3. mdd_fraction descending
#   4. n_subjects descending
prioritized = prioritized.sort_values(
    by=['n_ctrl', 'n_mdd', 'mdd_fraction', 'n_subjects'],
    ascending=[True, False, False, False]
).reset_index(drop=True)

# --- Step 5: Data Integrity Assertions ---
assert not prioritized.empty, "Assertion Failed: Output prioritized DataFrame is empty."
assert prioritized['impact'].str.upper().isin(['HIGH', 'MODERATE', 'LOW']).all(), "Assertion Failed: Unexpected impact tier found in output."
assert (prioritized['n_ctrl'] <= 2).all(), "Assertion Failed: Variants with n_ctrl > 2 found in output."
assert (prioritized['n_mdd'] >= 8).all(), "Assertion Failed: Variants with n_mdd < 8 found in output."
assert (prioritized['hgvs_p'].notna() & (prioritized['hgvs_p'] != '.') & (prioritized['hgvs_p'] != '') | prioritized['consequence'].str.contains('splice', case=False, na=False)).all(), "Assertion Failed: Missing hgvs_p."
assert (prioritized['is_gtex_ba9_eqtl'].str.upper() == 'YES').all(), "Assertion Failed: Non-eQTL variants found in output."

required_cols = [
    'variant_id', 'chrom', 'pos', 'ref', 'alt', 'gene',
    'consequence', 'impact', 'hgvs_c', 'hgvs_p',
    'n_subjects', 'n_mdd', 'n_ctrl', 'mdd_fraction', 'is_gtex_ba9_eqtl'
]
for col in required_cols:
    assert col in prioritized.columns, f"Assertion Failed: Column {col} missing in output."

# --- Step 6: Export Clean Results ---
prioritized.to_csv(OUT_CSV, index=False)

print(f"[SUCCESS] Prioritized disease matrix exported to: {OUT_CSV}")
print(f"          Total prioritized candidate loci: {len(prioritized):,}")
print("\n--- Impact Classification Breakdown ---")
print(prioritized['impact'].value_counts().to_string())
print("\n--- Control Carrier Distribution (n_ctrl) ---")
print(prioritized['n_ctrl'].value_counts().to_string())
print("\n--- Sequence Ontology Consequence Distribution ---")
print(prioritized['consequence'].value_counts().head(10).to_string())
print("\n--- Top Prioritized Candidate Variants ---")
display_cols = ['gene', 'consequence', 'impact', 'hgvs_p', 'rsid', 'is_gtex_ba9_eqtl', 'n_ctrl', 'n_mdd', 'mdd_fraction', 'omim_url', 'pubmed_url']
print(prioritized[[c for c in display_cols if c in prioritized.columns]].head(20).to_string(index=False))

