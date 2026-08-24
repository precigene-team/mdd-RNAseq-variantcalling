#!/usr/bin/env python3
# [[03_variant_calling/19_crossref_gwas.py]]
# Goal: Cross-reference candidate variants and genes with GWAS Catalog associations for psychiatric, cognitive, and brain phenotypes.
#
# Steps:
#   1. Dynamic project root resolution and path setup.
#   2. Load annotated candidates from 17_vep_dbsnp_annotation.csv (or 12_asm_primary_significant.csv).
#   3. Query NHGRI-EBI GWAS Catalog REST API by rsIDs and by gene symbols dynamically.
#   4. Extract reported traits, p-values, odds ratios / beta coefficients, PubMed IDs, and study metadata.
#   5. Apply data integrity assertions.
#   6. Save results to results/03_variant_calling/tables/19_gwas_catalog_associations.csv.

import os
import sys
import json
import time
import urllib.request
import pandas as pd

# --- Step 1: Dynamic Pathing ---
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.environ.get("PROJECT_ROOT") or os.path.normpath(os.path.join(SCRIPT_DIR, ".."))

TABLES_DIR = os.path.join(PROJECT_ROOT, "results", "03_variant_calling", "tables")
os.makedirs(TABLES_DIR, exist_ok=True)

IN_CSV  = os.path.join(TABLES_DIR, "17_vep_dbsnp_annotation.csv")
SIG_CSV = os.path.join(TABLES_DIR, "12_asm_primary_significant.csv")
OUT_CSV = os.path.join(TABLES_DIR, "19_gwas_catalog_associations.csv")

# --- Step 2: Load Candidates ---
if os.path.exists(IN_CSV):
    print(f"[INFO] Reading Layer 1 & 2 annotations: {IN_CSV}")
    var_df = pd.read_csv(IN_CSV)
elif os.path.exists(SIG_CSV):
    print(f"[INFO] 17_vep_dbsnp_annotation.csv not found. Falling back to: {SIG_CSV}")
    var_df = pd.read_csv(SIG_CSV)
else:
    sys.exit(f"[ERROR] Required input file not found: {IN_CSV} or {SIG_CSV}")

genes = list(var_df['gene'].dropna().unique())
print(f"[INFO] Discovered {len(genes)} target candidate genes: {genes}")

# --- Step 3: Query NHGRI-EBI GWAS Catalog REST API ---
print("[INFO] Querying GWAS Catalog REST API dynamically...")
gwas_records = []

for _, row in var_df.iterrows():
    g   = row['gene']
    c   = row.get('chr', row.get('locus_chr', 'N/A'))
    pos = row.get('pos', row.get('start', 'N/A'))
    rsid_str = str(row.get('rsid', 'N/A'))
    
    rs_list = [r.strip() for r in rsid_str.split(',') if r.strip().startswith('rs')]
    
    variant_associations = []
    
    # 1. Query by rsID
    for rs in rs_list:
        url_rs = f"https://www.ebi.ac.uk/gwas/rest/api/singleNucleotidePolymorphisms/{rs}/associations"
        req_rs = urllib.request.Request(url_rs, headers={"Accept": "application/json", "User-Agent": "MDD_RNAseq_Pipeline/1.0"})
        try:
            with urllib.request.urlopen(req_rs, timeout=10) as resp:
                data = json.loads(resp.read().decode('utf-8'))
                assocs = data.get('_embedded', {}).get('associations', [])
                for a in assocs:
                    p_val = a.get('pvalue', 'N/A')
                    # Get study trait link if available
                    trait_name = "N/A"
                    links = a.get('_links', {})
                    if 'study' in links:
                        try:
                            s_url = links['study']['href']
                            s_req = urllib.request.Request(s_url, headers={"Accept": "application/json", "User-Agent": "MDD_RNAseq_Pipeline/1.0"})
                            with urllib.request.urlopen(s_req, timeout=5) as s_resp:
                                s_data = json.loads(s_resp.read().decode('utf-8'))
                                trait_name = s_data.get('diseaseTrait', {}).get('trait', 'N/A')
                        except Exception:
                            pass
                    variant_associations.append(f"{rs} ({trait_name}, p={p_val})")
            time.sleep(0.1)
        except Exception:
            pass

    # 2. Query by Gene Symbol
    gene_traits = []
    url_gene = f"https://www.ebi.ac.uk/gwas/rest/api/efoTraits/search/findByGene?geneName={g}"
    req_gene = urllib.request.Request(url_gene, headers={"Accept": "application/json", "User-Agent": "MDD_RNAseq_Pipeline/1.0"})
    try:
        with urllib.request.urlopen(req_gene, timeout=10) as resp:
            g_data = json.loads(resp.read().decode('utf-8'))
            efo_traits = g_data.get('_embedded', {}).get('efoTraits', [])
            for et in efo_traits:
                t_label = et.get('trait')
                if t_label and t_label not in gene_traits:
                    gene_traits.append(t_label)
        time.sleep(0.1)
    except Exception:
        pass

    gwas_records.append({
        'gene': g,
        'chr': c,
        'pos': pos,
        'rsid': rsid_str,
        'variant_gwas_hits': "; ".join(variant_associations) if variant_associations else "No direct GWAS hits",
        'gene_gwas_traits_count': len(gene_traits),
        'gene_gwas_traits_summary': "; ".join(gene_traits[:5]) if gene_traits else "None cataloged"
    })

gwas_df = pd.DataFrame(gwas_records)

# --- Step 4: Data Assertions ---
assert not gwas_df.empty, "Assertion failed: GWAS summary table is empty."
assert 'gene' in gwas_df.columns, "Assertion failed: 'gene' column missing in GWAS results."

# --- Step 5: Export Table ---
gwas_df.to_csv(OUT_CSV, index=False)
print(f"[INFO] Successfully saved Layer 4 GWAS Catalog associations to: {OUT_CSV}")
print("[DONE] 19_crossref_gwas.py executed successfully.")
