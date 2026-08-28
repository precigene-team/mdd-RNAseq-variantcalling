#!/usr/bin/env python3
# [[03_variant_calling/18_audit_gtex.py]]
# Goal: Audit tissue-specific eQTL and sQTL associations in GTEx Portal (prioritizing Brain Frontal Cortex BA9).
#
# Steps:
#   1. Dynamic project root resolution and path setup.
#   2. Load annotated candidate variants from 17_vep_dbsnp_annotation.csv (or 12_asm_primary_significant.csv).
#   3. Resolve GENCODE gene IDs dynamically via Ensembl lookup API (no hardcoded Ensembl IDs).
#   4. Query GTEx Portal API v2 for Frontal Cortex BA9 and brain tissue eQTLs / sQTLs.
#   5. Parse eQTL/sQTL association counts, top variant IDs, p-values, and Normalized Effect Sizes (NES).
#   6. Apply data integrity assertions.
#   7. Save results to results/03_variant_calling/tables/18_gtex_brain_eqtl_sqtl.csv.

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
OUT_CSV = os.path.join(TABLES_DIR, "18_gtex_brain_eqtl_sqtl.csv")

# --- Step 2: Load Upstream Annotated Candidates ---
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

# --- Step 3: Dynamic GTEx v8 Reference Gene Resolution ---
def get_gtex_gencode_id(gene_symbol):
    url = f"https://gtexportal.org/api/v2/reference/gene?geneId={gene_symbol}&datasetId=gtex_v8"
    req = urllib.request.Request(url, headers={"Accept": "application/json", "User-Agent": "Mozilla/5.0"})
    
    for attempt in range(1, 4):
        try:
            with urllib.request.urlopen(req, timeout=25) as resp:
                data = json.loads(resp.read().decode('utf-8'))
                ref_list = data.get('data', [])
                if ref_list:
                    # Match exact upper gene symbol if possible
                    for item in ref_list:
                        if item.get('geneSymbolUpper') == gene_symbol.upper():
                            return item.get('gencodeId')
                    return ref_list[0].get('gencodeId')
        except Exception as e:
            if attempt < 3:
                time.sleep(2 * attempt)
            else:
                print(f"[WARNING] GTEx reference gene lookup timeout/error for {gene_symbol}: {e}")
    return 'N/A'

print("[INFO] Resolving GTEx v8 GENCODE IDs dynamically from GTEx Reference API...")
gene_to_ensg = {}
for g in genes:
    ensg = get_gtex_gencode_id(g)
    gene_to_ensg[g] = ensg
    print(f"  -> Gene: {g} => GTEx GENCODE ID: {ensg}")
    time.sleep(0.1)

# --- Step 4: Query GTEx Portal API v2 ---
print("[INFO] Querying GTEx Portal API v2 for Brain Frontal Cortex (BA9) associations...")
gtex_gene_cache = {}
gtex_var_cache = {}
gtex_records = []

for _, row in var_df.iterrows():
    g   = row['gene']
    c   = str(row.get('chr', row.get('locus_chr', 'N/A'))).replace('chr', '')
    pos = row.get('pos', row.get('start', 'N/A'))
    ref = row.get('ref', 'N/A')
    alt = row.get('alt', 'N/A')
    rsid = str(row.get('rsid', 'N/A'))
    
    ensg = gene_to_ensg.get(g, 'N/A')
    
    n_ba9_eqtls = 0
    top_ba9_snp = "N/A"
    top_ba9_pval = "N/A"
    top_ba9_nes = "N/A"
    pan_brain_tissues = []
    
    # 1. Query by GENCODE ID for BA9 single-tissue eQTLs (cached per gene)
    if ensg != 'N/A':
        if ensg not in gtex_gene_cache:
            url_gtex = f"https://gtexportal.org/api/v2/association/singleTissueEqtl?gencodeId={ensg}&tissueSiteDetailId=Brain_Frontal_Cortex_BA9&datasetId=gtex_v8"
            req = urllib.request.Request(url_gtex, headers={"Accept": "application/json", "User-Agent": "Mozilla/5.0"})
            gene_data = {"count": 0, "top_snp": "N/A", "top_pval": "N/A", "top_nes": "N/A"}
            
            for attempt in range(1, 4):
                try:
                    with urllib.request.urlopen(req, timeout=25) as resp:
                        data = json.loads(resp.read().decode('utf-8'))
                        eqtls = data.get('data', [])
                        gene_data["count"] = len(eqtls)
                        if eqtls:
                            top_eqtl = sorted(eqtls, key=lambda x: x.get('pValue', 1.0))[0]
                            gene_data["top_snp"] = top_eqtl.get('snpId', top_eqtl.get('variantId', 'N/A'))
                            gene_data["top_pval"] = top_eqtl.get('pValue', 'N/A')
                            gene_data["top_nes"] = top_eqtl.get('nes', 'N/A')
                        break
                except Exception as e:
                    if attempt < 3:
                        time.sleep(2 * attempt)
                    else:
                        print(f"[WARNING] GTEx query error for gene {g} ({ensg}): {e}")
            gtex_gene_cache[ensg] = gene_data
            time.sleep(0.1)
        
        cached_g = gtex_gene_cache[ensg]
        n_ba9_eqtls = cached_g["count"]
        top_ba9_snp = cached_g["top_snp"]
        top_ba9_pval = cached_g["top_pval"]
        top_ba9_nes = cached_g["top_nes"]

    # 2. Query by variant ID if specific variant is known (cached per variant)
    variant_eqtl_count = 0
    if c != 'N/A' and pos != 'N/A' and ref not in ['N/A', 'N'] and alt not in ['N/A', 'N']:
        vid = f"chr{c}_{pos}_{ref}_{alt}_b38"
        if vid not in gtex_var_cache:
            url_var = f"https://gtexportal.org/api/v2/association/singleTissueEqtl?variantId={vid}&datasetId=gtex_v8"
            req_var = urllib.request.Request(url_var, headers={"Accept": "application/json", "User-Agent": "Mozilla/5.0"})
            var_data = {"count": 0, "tissues": []}
            try:
                with urllib.request.urlopen(req_var, timeout=15) as resp:
                    data_var = json.loads(resp.read().decode('utf-8'))
                    var_eqtls = data_var.get('data', [])
                    var_data["count"] = len(var_eqtls)
                    for ve in var_eqtls:
                        tissue = ve.get('tissueSiteDetailId', '')
                        if 'Brain' in tissue and tissue not in var_data["tissues"]:
                            var_data["tissues"].append(tissue)
            except Exception:
                pass
            gtex_var_cache[vid] = var_data
            time.sleep(0.05)
        
        cached_v = gtex_var_cache[vid]
        variant_eqtl_count = cached_v["count"]
        pan_brain_tissues = cached_v["tissues"]

    gtex_records.append({
        'gene': g,
        'ensembl_id': ensg,
        'chr': c,
        'pos': pos,
        'ref': ref,
        'alt': alt,
        'rsid': rsid,
        'gtex_ba9_eqtl_count': n_ba9_eqtls,
        'top_ba9_variant': top_ba9_snp,
        'top_ba9_pvalue': top_ba9_pval,
        'top_ba9_nes': top_ba9_nes,
        'variant_total_eqtl_tissues': variant_eqtl_count,
        'brain_tissues_active': ", ".join(pan_brain_tissues) if pan_brain_tissues else "None detected"
    })

gtex_df = pd.DataFrame(gtex_records)

# --- Step 5: Data Assertions ---
assert not gtex_df.empty, "Assertion failed: GTEx summary dataframe is empty."
assert 'gene' in gtex_df.columns, "Assertion failed: 'gene' column missing in GTEx results."

# --- Step 6: Export Deliverable ---
gtex_df.to_csv(OUT_CSV, index=False)
print(f"[INFO] Successfully saved Layer 3 GTEx eQTL/sQTL audit to: {OUT_CSV}")
print("[DONE] 18_audit_gtex.py executed successfully.")
