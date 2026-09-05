#!/usr/bin/env python3
# [[03_variant_calling/12_pharmacogenomics_actionability.py]]
# Goal: Official CPIC catalog-wide pharmacogenomics actionability mapping for prioritized Input transcriptome SNPs.
#
# Steps:
#   1. Dynamic project root resolution and path setup.
#   2. Load clinically classified Input transcriptome SNPs from 11_clinical_pathogenicity.csv.
#   3. Ingest official CPIC reference database (_data/references/cpic_pairs.json).
#   4. Perform catalog-wide unbiased join against official CPIC gene-drug pairs.
#   5. Assign ClinPGx Clinical Actionability Categories (Actionable, Informative, Unclassified).
#   6. Perform data integrity assertions.
#   7. Save clean pharmacogenomic results to results/tables/12_pharmacogenomics_actionability.csv.

import os
import sys
import json
import pandas as pd
import numpy as np

# --- Step 1: Dynamic Pathing ---
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.environ.get("PROJECT_ROOT") or os.path.normpath(os.path.join(SCRIPT_DIR, ".."))

TABLES_DIR = os.path.join(PROJECT_ROOT, "results", "tables")
MOD_TABLES_DIR = os.path.join(PROJECT_ROOT, "results", "03_variant_calling", "tables")
REFS_DIR = os.path.join(PROJECT_ROOT, "_data", "references")

os.makedirs(TABLES_DIR, exist_ok=True)
os.makedirs(MOD_TABLES_DIR, exist_ok=True)

INPUT_SOURCES = [
    os.path.join(TABLES_DIR, "11_clinical_pathogenicity.csv"),
    os.path.join(MOD_TABLES_DIR, "11_clinical_pathogenicity.csv")
]

in_csv = next((f for f in INPUT_SOURCES if os.path.exists(f)), None)
if not in_csv:
    sys.exit(f"[ERROR] Required input 11_clinical_pathogenicity.csv not found in: {INPUT_SOURCES}")

CPIC_JSON = os.path.join(REFS_DIR, "cpic_pairs.json")
if not os.path.exists(CPIC_JSON):
    sys.exit(f"[ERROR] Official CPIC reference database not found: {CPIC_JSON}")

OUT_CSV = os.path.join(TABLES_DIR, "12_pharmacogenomics_actionability.csv")
MOD_OUT_CSV = os.path.join(MOD_TABLES_DIR, "12_pharmacogenomics_actionability.csv")

# --- Step 2: Ingest Classified Variants ---
print(f"[INFO] Ingesting clinical classifications from: {in_csv}")
df = pd.read_csv(in_csv)
if df.empty:
    sys.exit("[ERROR] Ingested clinical classification table is empty.")

print(f"[INFO] Evaluating pharmacogenomics for {len(df)} candidate variants across {df['gene'].nunique()} genes...")

# --- Step 3: Load Official CPIC Database ---
print(f"[INFO] Ingesting official CPIC database from: {CPIC_JSON}")
with open(CPIC_JSON, "r") as f:
    cpic_raw = json.load(f)

cpic_map = {}
for item in cpic_raw:
    sym = item.get("genesymbol")
    if not sym:
        continue
    if sym not in cpic_map:
        cpic_map[sym] = {
            'drugs': set(),
            'cpic_levels': set(),
            'clinpgx_levels': set(),
            'guidelines': set(),
            'pgx_testing': set()
        }
    if item.get("drugname"):
        cpic_map[sym]['drugs'].add(item["drugname"])
    if item.get("cpiclevel"):
        cpic_map[sym]['cpic_levels'].add(item["cpiclevel"])
    if item.get("clinpgxlevel"):
        cpic_map[sym]['clinpgx_levels'].add(item["clinpgxlevel"])
    if item.get("guidelinename"):
        cpic_map[sym]['guidelines'].add(item["guidelinename"])
    if item.get("pgxtesting"):
        cpic_map[sym]['pgx_testing'].add(item["pgxtesting"])

print(f"[INFO] Successfully indexed {len(cpic_map)} authoritative CPIC pharmacogenes.")

# --- Step 4: Catalog-Wide Unbiased Join ---
print("[INFO] Matching cohort candidate genes against official CPIC catalog...")

pgx_records = []
for _, r in df.iterrows():
    gene = str(r.get('gene', 'Unknown'))
    cpic_info = cpic_map.get(gene)
    
    if cpic_info:
        drugs_str = "; ".join(sorted(cpic_info['drugs']))
        cpic_lvl_str = "; ".join(sorted(cpic_info['cpic_levels']))
        cpgx_lvl_str = "; ".join(sorted(cpic_info['clinpgx_levels']))
        guide_str = "; ".join(sorted(cpic_info['guidelines']))
        testing_str = "; ".join(sorted(cpic_info['pgx_testing']))
        
        if any(lvl in ['A', 'B'] for lvl in cpic_info['cpic_levels']):
            action_cat = "Level_A_Actionable"
        elif any(lvl in ['C', 'D'] for lvl in cpic_info['cpic_levels']):
            action_cat = "Level_B_Informative"
        else:
            action_cat = "Level_B_Informative"
            
        pgx_records.append({
            'cpic_level': cpic_lvl_str,
            'clinpgx_level': cpgx_lvl_str,
            'drug_indications': drugs_str,
            'pgx_testing_status': testing_str,
            'guideline_name': guide_str,
            'actionability_category': action_cat
        })
    else:
        pgx_records.append({
            'cpic_level': 'None',
            'clinpgx_level': 'None',
            'drug_indications': 'None',
            'pgx_testing_status': 'Unclassified',
            'guideline_name': 'No_CPIC_Guideline',
            'actionability_category': 'Unclassified'
        })

df_pgx = pd.DataFrame(pgx_records)
for col in df_pgx.columns:
    df[col] = df_pgx[col]

# --- Step 5: Data Integrity Assertions ---
assert not df.empty, "Assertion Failed: Output DataFrame is empty."
assert len(df) == len(df_pgx), "Assertion Failed: Mismatched row count after CPIC merge."

required_cols = [
    'variant_id', 'gene', 'consequence', 'impact', 'acmg_classification',
    'cpic_level', 'clinpgx_level', 'drug_indications', 'actionability_category'
]
for col in required_cols:
    assert col in df.columns, f"Assertion Failed: Column {col} missing in output."

# --- Step 6: Export Clean Results ---
df.to_csv(OUT_CSV, index=False)
df.to_csv(MOD_OUT_CSV, index=False)

print(f"[SUCCESS] Pharmacogenomics Actionability exported to:")
print(f"  Primary: {OUT_CSV}")
print(f"  Module:  {MOD_OUT_CSV}")
print("\n--- CPIC Actionability Distribution in Cohort ---")
print(df['actionability_category'].value_counts().to_string())
