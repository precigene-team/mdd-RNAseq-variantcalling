#!/usr/bin/env python3
# [[03_variant_calling/10_prioritize_mdd_candidate_snps.py]]
# Goal: Take the annotated SNPs from 09, keep the protein-altering candidates that are
#       MDD-enriched (n_ctrl <= 2 & n_mdd >= 8), and attach the external evidence URLs
#       for manual review. This is a terminal table: no script consumes it.
#
# Steps:
#   1. Dynamic project root resolution and path setup.
#   2. Load annotated candidates from 09_functional_annotation.csv.
#   3. Filter for functional candidates (HIGH & MODERATE & LOW impact, valid hgvs_p, coding/splice).
#   4. Strict filter: restrict to n_ctrl <= 2 & n_mdd >= 8, sort MDD descending.
#   5. Build the external resource URLs by string formatting (OMIM, DECIPHER, PubMed,
#      gnomAD gene + variant, GENCC, GWAS Catalog gene). No network call is made.
#   6. Apply strict data assertions.
#   7. Save the candidate table to results/03_variant_calling/tables/10_mdd_candidate_snps_urls.csv.
#
# No rsID, and why: 09 does not carry rsid and no script in this repo produces one, so the
# GWAS Catalog link here is gene-level. The variant-level Catalog lookup requires rsID, and
# the Catalog's position-based REST endpoint (singleNucleotidePolymorphisms/chromosome/
# position/) returns 410 Gone. gnomAD's variant URL is built from chrom:pos:ref:alt and
# needs no rsID, so per-variant depth is reachable through gnomAD.

import os
import sys
import pandas as pd
import numpy as np

# --- Step 1: Dynamic Pathing ---
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.environ.get("PROJECT_ROOT") or os.path.normpath(os.path.join(SCRIPT_DIR, ".."))

TABLES_DIR = os.path.join(PROJECT_ROOT, "results", "03_variant_calling", "tables")
os.makedirs(TABLES_DIR, exist_ok=True)

IN_ANNOT_CSV = os.path.join(TABLES_DIR, "09_functional_annotation.csv")
OUT_CSV = os.path.join(TABLES_DIR, "10_mdd_candidate_snps_urls.csv")

# Builds the external resource URLs by string formatting from columns 09 already writes:
# gene symbol for the gene-level links, chrom:pos:ref:alt for gnomAD. No network call is
# made, so 10 has no API dependency and no rate limits.
def attach_resource_urls(df):
    df['omim_url'] = df['gene'].apply(lambda g: f"https://omim.org/search?search={g}")
    df['decipher_url'] = df['gene'].apply(lambda g: f"https://www.deciphergenomics.org/gene/{g}/overview/protein-genomic-info")
    df['pubmed_url'] = df['gene'].apply(lambda g: f"https://pubmed.ncbi.nlm.nih.gov/?term={g}+AND+(depression+OR+MDD)")
    df['gnomad_gene_url'] = df['gene'].apply(lambda g: f"https://gnomad.broadinstitute.org/gene/{g}?dataset=gnomad_r4")
    df['gnomad_variant_url'] = df['variant_id'].apply(lambda v: f"https://gnomad.broadinstitute.org/variant/{str(v).replace(':', '-')}?dataset=gnomad_r4")
    df['gencc_url'] = df['gene'].apply(lambda g: f"https://search.thegencc.org/genes/{g}")
    df['gwas_catalog_gene_url'] = df['gene'].apply(lambda g: f"https://www.ebi.ac.uk/gwas/genes/{g}")
    return df

# --- Step 2: Load Candidates ---
if not os.path.exists(IN_ANNOT_CSV):
    sys.exit(f"[ERROR] Required table not found: {IN_ANNOT_CSV}. Run 09_functional_annotation.py first.")
print(f"[INFO] Loading candidates from: {IN_ANNOT_CSV}...")
df_cand = pd.read_csv(IN_ANNOT_CSV, low_memory=False)

print(f"[INFO] Loaded {len(df_cand):,} annotated candidates.")

# --- Step 3: Filter for Functional Candidates ---
has_hgvs_p = df_cand['hgvs_p'].apply(
    lambda x: pd.notna(x) and str(x).strip() not in ('None', '.', 'N/A', '', 'nan')
)
is_functional = df_cand['impact'].str.upper().isin(['HIGH', 'MODERATE', 'LOW']) & (
    has_hgvs_p | df_cand['consequence'].str.contains('splice', case=False, na=False)
)

prioritized = df_cand[is_functional].copy()

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

# --- Step 5: Attach External Evidence URLs ---
prioritized = attach_resource_urls(prioritized)

# --- Step 6: Data Integrity Assertions ---
assert not prioritized.empty, "Assertion Failed: Output prioritized DataFrame is empty."
assert prioritized['impact'].str.upper().isin(['HIGH', 'MODERATE', 'LOW']).all(), "Assertion Failed: Unexpected impact tier found in output."
assert (prioritized['n_ctrl'] <= 2).all(), "Assertion Failed: Variants with n_ctrl > 2 found in output."
assert (prioritized['n_mdd'] >= 8).all(), "Assertion Failed: Variants with n_mdd < 8 found in output."
assert (prioritized['hgvs_p'].notna() & (prioritized['hgvs_p'] != '.') & (prioritized['hgvs_p'] != '') | prioritized['consequence'].str.contains('splice', case=False, na=False)).all(), "Assertion Failed: Missing hgvs_p."

required_cols = [
    'variant_id', 'chrom', 'pos', 'ref', 'alt', 'gene',
    'consequence', 'impact', 'hgvs_c', 'hgvs_p',
    'n_subjects', 'n_mdd', 'n_ctrl', 'mdd_fraction',
    'omim_url', 'decipher_url', 'pubmed_url', 'gnomad_gene_url',
    'gnomad_variant_url', 'gencc_url', 'gwas_catalog_gene_url'
]
for col in required_cols:
    assert col in prioritized.columns, f"Assertion Failed: Column {col} missing in output."

# --- Step 7: Export Clean Results ---
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
display_cols = ['gene', 'consequence', 'impact', 'hgvs_p', 'n_ctrl', 'n_mdd', 'mdd_fraction', 'omim_url', 'gwas_catalog_gene_url', 'gnomad_variant_url']
print(prioritized[[c for c in display_cols if c in prioritized.columns]].head(20).to_string(index=False))

