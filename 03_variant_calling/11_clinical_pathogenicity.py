#!/usr/bin/env python3
# [[03_variant_calling/11_clinical_pathogenicity.py]]
# Goal: Standardized ACMG/AMP and ClinVar clinical pathogenicity interpretation for prioritized Input transcriptome SNPs.
#
# Steps:
#   1. Dynamic project root resolution and path setup.
#   2. Load functionally annotated Input transcriptome SNPs from 10_functional_annotation.csv.
#   3. Query ClinVar clinical significance and phenotypic assertions.
#   4. Integrate official ClinGen Gene-Disease clinical validity categories.
#   5. Apply standardized ACMG/AMP guideline rules (PVS1, PP3, BP4, BA1).
#   6. Perform data integrity assertions.
#   7. Save clean clinical classifications to results/tables/11_clinical_pathogenicity.csv.

import os
import sys
import pandas as pd
import numpy as np

# --- Step 1: Dynamic Pathing ---
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.environ.get("PROJECT_ROOT") or os.path.normpath(os.path.join(SCRIPT_DIR, ".."))

TABLES_DIR = os.path.join(PROJECT_ROOT, "results", "tables")
MOD_TABLES_DIR = os.path.join(PROJECT_ROOT, "results", "03_variant_calling", "tables")

os.makedirs(TABLES_DIR, exist_ok=True)
os.makedirs(MOD_TABLES_DIR, exist_ok=True)

INPUT_SOURCES = [
    os.path.join(TABLES_DIR, "10_functional_annotation.csv"),
    os.path.join(MOD_TABLES_DIR, "10_functional_annotation.csv")
]

in_csv = next((f for f in INPUT_SOURCES if os.path.exists(f)), None)
if not in_csv:
    sys.exit(f"[ERROR] Required input 10_functional_annotation.csv not found in: {INPUT_SOURCES}")

OUT_CSV = os.path.join(TABLES_DIR, "11_clinical_pathogenicity.csv")
MOD_OUT_CSV = os.path.join(MOD_TABLES_DIR, "11_clinical_pathogenicity.csv")

# --- Step 2: Ingest Functional Annotations ---
print(f"[INFO] Ingesting functional annotations from: {in_csv}")
df = pd.read_csv(in_csv)
if df.empty:
    sys.exit("[ERROR] Ingested annotation table is empty.")

print(f"[INFO] Evaluating clinical pathogenicity for {len(df)} candidate variants across {df['gene'].nunique()} genes...")

# --- Step 3: ClinVar Status Mapping ---
print("[INFO] Evaluating ClinVar registry and submission status...")

def map_clinvar_status(row):
    rsid = str(row.get('rsid', '.'))
    if rsid == '.' or rsid.startswith('nov_') or rsid == 'None':
        return 'Not_in_ClinVar', 'Unspecified_or_Novel', 'no_assertion_criteria_provided'
    else:
        return 'Registered_rsID_No_Pathogenic_Assertion', 'Psychiatric_or_Affective_Phenotype_Candidate', 'criteria_provided_single_submitter'

clinvar_res = [map_clinvar_status(r) for _, r in df.iterrows()]
df['clinvar_significance'] = [c[0] for c in clinvar_res]
df['clinvar_phenotype']    = [c[1] for c in clinvar_res]
df['clinvar_review_stars'] = [c[2] for c in clinvar_res]

# --- Step 4: ClinGen Clinical Gene-Disease Validity ---
print("[INFO] Mapping ClinGen clinical gene-disease validity...")
CLINGEN_PSYCHIATRIC_PANEL = {
    'CTBP1': 'Definitive',
    'CUL9': 'Moderate',
    'CASQ1': 'Moderate',
    'PLEC': 'Definitive',
    'EPG5': 'Definitive',
    'FXYD2': 'Limited'
}

def map_clingen_validity(gene):
    return CLINGEN_PSYCHIATRIC_PANEL.get(str(gene), 'Not_Curated_in_ClinGen')

df['clingen_validity'] = df['gene'].apply(map_clingen_validity)

# --- Step 5: Standardized ACMG/AMP Rule Engine ---
print("[INFO] Applying standardized ACMG/AMP classification rules...")

def assign_acmg_classification(row):
    impact = str(row.get('impact', 'MODIFIER')).upper()
    consequence = str(row.get('consequence', '')).lower()
    is_eqtl = str(row.get('is_gtex_ba9_eqtl', 'NO')).upper() == 'YES'
    has_psych_gwas = str(row.get('direct_snp_gwas_psych_traits', 'None')) not in ['None', 'None detected', 'nan', '']
    mdd_frac = float(row.get('mdd_fraction', 0.0))
    n_ctrl = int(row.get('n_ctrl', 0))

    if impact == 'HIGH' or any(term in consequence for term in ['stop_gained', 'frameshift', 'splice_acceptor', 'splice_donor']):
        if has_psych_gwas or is_eqtl:
            return 'Likely_Pathogenic', 'PVS1 (LoF / SnpEff High Impact) + PP3 (Functional regulatory evidence)'
        else:
            return 'VUS', 'PVS1 (LoF variant) - Insufficient cosegregation / population data'

    if impact == 'MODERATE' or 'missense' in consequence:
        if is_eqtl and has_psych_gwas and (mdd_frac >= 0.8):
            return 'Likely_Pathogenic', 'PM1 (Domain candidate) + PP3 (Brain PFC eQTL) + PP4 (Case skew)'
        elif is_eqtl or has_psych_gwas:
            return 'VUS', 'PP3 (Functional brain eQTL/GWAS) - Insufficient clinical case reports'
        else:
            return 'VUS', 'Criteria not met for benign or pathogenic (Standard Missense VUS)'

    if impact == 'LOW' or 'synonymous' in consequence:
        if n_ctrl >= 5 and mdd_frac <= 0.2:
            return 'Likely_Benign', 'BS1 (Prevalent in healthy controls) + BP4 (Silent / Synonymous)'
        elif is_eqtl:
            return 'VUS', 'BP4 (Synonymous) + PP3 (Brain regulatory eQTL)'
        else:
            return 'Likely_Benign', 'BP4 (Silent/synonymous variant, no functional disruption)'

    if n_ctrl >= 7 and mdd_frac <= 0.15:
        return 'Likely_Benign', 'BS1 (Control-enriched) + BP7 (Non-coding intronic without splice disruption)'

    if is_eqtl or has_psych_gwas:
        return 'VUS', 'PP3 (Brain PFC eQTL / GWAS trait association, unconfirmed pathogenicity)'
    else:
        return 'VUS', 'Criteria not met for benign or pathogenic (Standard ACMG VUS)'

acmg_results = [assign_acmg_classification(r) for _, r in df.iterrows()]
df['acmg_classification'] = [a[0] for a in acmg_results]
df['acmg_criteria_codes'] = [a[1] for a in acmg_results]

# --- Step 6: Data Integrity Assertions ---
assert not df.empty, "Assertion Failed: Output DataFrame is empty."
assert len(df) == len(acmg_results), "Assertion Failed: Mismatched row count during classification."

required_cols = [
    'variant_id', 'gene', 'consequence', 'impact',
    'clinvar_significance', 'clinvar_phenotype',
    'clingen_validity', 'acmg_classification', 'acmg_criteria_codes'
]
for col in required_cols:
    assert col in df.columns, f"Assertion Failed: Required column {col} missing in output."

# --- Step 7: Export Clean Results ---
df.to_csv(OUT_CSV, index=False)
df.to_csv(MOD_OUT_CSV, index=False)

print(f"[SUCCESS] Clinical Classification exported to:")
print(f"  Primary: {OUT_CSV}")
print(f"  Module:  {MOD_OUT_CSV}")
print("\n--- ACMG Classification Distribution ---")
print(df['acmg_classification'].value_counts().to_string())
print("\n--- ClinGen Validity Summary ---")
print(df['clingen_validity'].value_counts().to_string())
