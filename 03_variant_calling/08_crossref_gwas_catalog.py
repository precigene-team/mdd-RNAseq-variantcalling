#!/usr/bin/env python3
# [[03_variant_calling/08_crossref_gwas_catalog.py]]
# Goal: Sequential dual-tier GWAS cross-referencing: (1) Direct SNP rsID associations and (2) Target gene-level disease traits.
#
# Steps:
#   1. Dynamic project root resolution and path setup.
#   2. Load prioritized candidates from 07_gtex_brain_eqtl_matches.csv.
#   3. Query direct variant rsIDs sequentially from GWAS Catalog.
#   4. Query unique target genes sequentially from GWAS Catalog with psychiatric keyword filtering.
#   5. Map results to all candidate loci and export results/tables/08_gwas_catalog_neuro_matches.csv.

import os
import sys
import json
import time
import urllib.request
import urllib.parse
import pandas as pd

# --- Step 1: Dynamic Pathing ---
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.environ.get("PROJECT_ROOT") or os.path.normpath(os.path.join(SCRIPT_DIR, ".."))

TABLES_DIR = os.path.join(PROJECT_ROOT, "results", "tables")
MOD_TABLES_DIR = os.path.join(PROJECT_ROOT, "results", "03_variant_calling", "tables")

for d in [TABLES_DIR, MOD_TABLES_DIR]:
    os.makedirs(d, exist_ok=True)

INPUT_SOURCES = [
    os.path.join(TABLES_DIR, "07_gtex_brain_eqtl_matches.csv"),
    os.path.join(MOD_TABLES_DIR, "07_gtex_brain_eqtl_matches.csv")
]

IN_GTEX_CSV = next((f for f in INPUT_SOURCES if os.path.exists(f)), None)
OUT_CSV = os.path.join(TABLES_DIR, "08_gwas_catalog_neuro_matches.csv")

if not IN_GTEX_CSV or not os.path.exists(IN_GTEX_CSV):
    sys.exit(f"[ERROR] Required input file not found in: {INPUT_SOURCES}. Run 07_map_gtex_brain_eqtls.py first.")

var_df = pd.read_csv(IN_GTEX_CSV)
print(f"[INFO] Loaded {len(var_df)} candidate loci from: {IN_GTEX_CSV}")

PSYCH_KEYWORDS = [
    'depress', 'mood', 'suicid', 'psych', 
    'cognit', 'neuro', 'brain', 'trauma',
    'anxiety', 'intellect', 'behavior', 'major depressive', 'stress'
]

# --- Step 2: Direct Variant-Level GWAS Query (by rsID) ---
if 'rsid' in var_df.columns:
    unique_rsids = [r for r in var_df['rsid'].dropna().unique() if str(r).startswith('rs')]
else:
    unique_rsids = []

print(f"[INFO] Querying {len(unique_rsids)} unique dbSNP rsIDs sequentially...")

def fetch_single_rsid(rsid):
    url = f"https://www.ebi.ac.uk/gwas/rest/api/singleNucleotidePolymorphisms/{rsid}"
    req = urllib.request.Request(url, headers={'Accept': 'application/json', 'User-Agent': 'MDD-RNAseq-Pipeline/1.0'})
    try:
        with urllib.request.urlopen(req, timeout=10) as response:
            if response.status == 200:
                data = json.loads(response.read().decode('utf-8'))
                traits = []
                for study in data.get('_embedded', {}).get('studies', []):
                    trait = study.get('diseaseTrait', {}).get('trait')
                    if trait and trait not in traits:
                        traits.append(trait)
                return traits
    except Exception:
        pass
    return []

rsid_trait_map = {}
for i, rs in enumerate(unique_rsids, 1):
    traits = fetch_single_rsid(rs)
    rsid_trait_map[rs] = traits
    time.sleep(0.1)
    if i % 10 == 0 or i == len(unique_rsids):
        print(f"  [GWAS rsID Query] {i}/{len(unique_rsids)} processed.")

# --- Step 3: Target Gene-Level GWAS Query ---
if 'target_gene_symbol' in var_df.columns:
    unique_genes = [g for g in var_df['target_gene_symbol'].dropna().unique() if g not in ('None', 'Unknown', '')]
else:
    unique_genes = []

print(f"[INFO] Querying {len(unique_genes)} unique target genes sequentially...")

def fetch_single_gene(gene):
    url = f"https://www.ebi.ac.uk/gwas/rest/api/genes/{gene}"
    req = urllib.request.Request(url, headers={'Accept': 'application/json', 'User-Agent': 'MDD-RNAseq-Pipeline/1.0'})
    try:
        with urllib.request.urlopen(req, timeout=10) as response:
            if response.status == 200:
                data = json.loads(response.read().decode('utf-8'))
                studies = data.get('_embedded', {}).get('studies', [])
                traits = []
                for s in studies:
                    t = s.get('diseaseTrait', {}).get('trait')
                    if t and t not in traits:
                        traits.append(t)
                return traits
    except Exception:
        pass
    return []

gene_trait_map = {}
for i, g in enumerate(unique_genes, 1):
    traits = fetch_single_gene(g)
    gene_trait_map[g] = traits
    time.sleep(0.1)
    if i % 10 == 0 or i == len(unique_genes):
        print(f"  [GWAS Gene Query] {i}/{len(unique_genes)} processed.")

# --- Step 4: Cross-Reference Mapping ---
print("[INFO] Mapping traits back to candidate variants...")

def filter_psych(traits):
    if not traits:
        return []
    return [t for t in traits if any(k in t.lower() for k in PSYCH_KEYWORDS)]

direct_all, direct_psych = [], []
gene_all, gene_psych = [], []

for _, row in var_df.iterrows():
    rs = str(row.get('rsid', ''))
    g = str(row.get('target_gene_symbol', ''))
    
    r_traits = rsid_trait_map.get(rs, [])
    r_psych = filter_psych(r_traits)
    direct_all.append("; ".join(r_traits[:5]) if r_traits else "None cataloged")
    direct_psych.append("; ".join(r_psych) if r_psych else "None detected")
    
    g_traits = gene_trait_map.get(g, [])
    g_psych = filter_psych(g_traits)
    gene_all.append("; ".join(g_traits[:5]) if g_traits else "None cataloged")
    gene_psych.append("; ".join(g_psych) if g_psych else "None detected")

var_df['direct_snp_gwas_psych_traits'] = direct_psych
var_df['direct_snp_gwas_all_traits'] = direct_all
var_df['target_gene_gwas_psych_traits'] = gene_psych
var_df['target_gene_gwas_all_traits'] = gene_all

# --- Step 5: Export Annotated Table ---
def write_df(target_path):
    var_df.to_csv(target_path, index=False)

write_df(OUT_CSV)
write_df(os.path.join(MOD_TABLES_DIR, "08_gwas_catalog_neuro_matches.csv"))

print(f"[SUCCESS] Exported GWAS Catalog results to: {OUT_CSV}")
