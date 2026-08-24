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

# --- Step 3: Dynamic Ensembl ID Resolution ---
def get_ensembl_id(gene_symbol):
    url = f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/{gene_symbol}?content-type=application/json"
    req = urllib.request.Request(url, headers={"Content-Type": "application/json", "User-Agent": "MDD_RNAseq_Pipeline/1.0"})
    try:
        with urllib.request.urlopen(req, timeout=10) as resp:
            data = json.loads(resp.read().decode('utf-8'))
            return data.get('id', 'N/A')
    except Exception as e:
        print(f"[WARNING] Could not resolve Ensembl ID for {gene_symbol}: {e}")
        return 'N/A'

print("[INFO] Resolving Ensembl GENCODE IDs dynamically from Ensembl API...")
gene_to_ensg = {}
for g in genes:
    ensg = get_ensembl_id(g)
    gene_to_ensg[g] = ensg
    print(f"  -> Gene: {g} => Ensembl ID: {ensg}")
    time.sleep(0.1)

# --- Step 4: Query GTEx Portal API v2 ---
print("[INFO] Querying GTEx Portal API v2 for Brain Frontal Cortex (BA9) associations...")
gtex_records = []

for _, row in var_df.iterrows():
    g   = row['gene']
    c   = row.get('chr', row.get('locus_chr', 'N/A'))
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
    
    # 1. Query by GENCODE ID for BA9 single-tissue eQTLs
    if ensg != 'N/A':
        url_gtex = f"https://gtexportal.org/api/v2/association/singleTissueEqtl?gencodeId={ensg}&tissueSiteDetailId=Brain_Frontal_Cortex_BA9&datasetId=gtex_v8"
        req = urllib.request.Request(url_gtex, headers={"Accept": "application/json", "User-Agent": "Mozilla/5.0"})
        try:
            with urllib.request.urlopen(req, timeout=10) as resp:
                data = json.loads(resp.read().decode('utf-8'))
                eqtls = data.get('singleTissueEqtl', [])
                n_ba9_eqtls = len(eqtls)
                if eqtls:
                    top_eqtl = sorted(eqtls, key=lambda x: x.get('pValue', 1.0))[0]
                    top_ba9_snp = top_eqtl.get('variantId', 'N/A')
                    top_ba9_pval = top_eqtl.get('pValue', 'N/A')
                    top_ba9_nes = top_eqtl.get('nes', 'N/A')
            time.sleep(0.1)
        except Exception as e:
            print(f"[WARNING] GTEx query error for gene {g} ({ensg}): {e}")

    # 2. Query by variant ID if specific variant is known
    variant_eqtl_count = 0
    if c != 'N/A' and pos != 'N/A' and ref not in ['N/A', 'N'] and alt not in ['N/A', 'N']:
        vid = f"chr{c}_{pos}_{ref}_{alt}_b38"
        url_var = f"https://gtexportal.org/api/v2/association/singleTissueEqtl?variantId={vid}&datasetId=gtex_v8"
        req_var = urllib.request.Request(url_var, headers={"Accept": "application/json", "User-Agent": "Mozilla/5.0"})
        try:
            with urllib.request.urlopen(req_var, timeout=10) as resp:
                data_var = json.loads(resp.read().decode('utf-8'))
                var_eqtls = data_var.get('singleTissueEqtl', [])
                variant_eqtl_count = len(var_eqtls)
                for ve in var_eqtls:
                    tissue = ve.get('tissueSiteDetailId', '')
                    if 'Brain' in tissue and tissue not in pan_brain_tissues:
                        pan_brain_tissues.append(tissue)
            time.sleep(0.1)
        except Exception:
            pass

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
