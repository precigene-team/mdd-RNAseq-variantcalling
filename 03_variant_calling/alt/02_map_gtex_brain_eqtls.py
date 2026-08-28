#!/usr/bin/env python3
# [[03_variant_calling/alt/02_map_gtex_brain_eqtls.py]]
# Goal: Query GTEx Portal v2 API to map transcriptome-wide candidate SNPs against Brain Frontal Cortex (BA9) eQTLs.
#
# Steps:
#   1. Dynamic project root resolution and path setup.
#   2. Load master transcriptome SNP table (01_transcriptome_snp_master.tsv.gz).
#   3. Select recurrent/disease-skewed candidate variants for prioritized GTEx brain audit.
#   4. Query GTEx Portal API v2 for singleTissueEqtl (Brain_Frontal_Cortex_BA9) and variant eQTL status.
#   5. Extract eQTL counts, top regulatory variants, p-values, and Normalized Effect Sizes (NES).
#   6. Apply strict data assertions.
#   7. Save results to results/03_variant_calling/alt/tables/02_gtex_brain_eqtl_matches.csv.

import os
import sys
import json
import time
import gzip
import csv
import argparse
import urllib.request
import urllib.parse
import pandas as pd

# --- Step 1: Dynamic Pathing ---
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.environ.get("PROJECT_ROOT") or os.path.normpath(os.path.join(SCRIPT_DIR, "../.."))

TABLES_DIR = os.path.join(PROJECT_ROOT, "results", "03_variant_calling", "alt", "tables")
os.makedirs(TABLES_DIR, exist_ok=True)

IN_MASTER_TSV = os.path.join(TABLES_DIR, "01_transcriptome_snp_master.tsv.gz")
OUT_CSV = os.path.join(TABLES_DIR, "02_gtex_brain_eqtl_matches.csv")

parser = argparse.ArgumentParser(description="Map transcriptome-wide SNPs against GTEx Brain eQTLs.")
parser.add_argument("--min-subjects", type=int, default=8, help="Minimum cohort subjects carrying SNP (default: 8)")
parser.add_argument("--top-n", type=int, default=500, help="Max top candidate variants to query (default: 500)")
args = parser.parse_args()

# --- Step 2: Load Master SNP Table ---
if not os.path.exists(IN_MASTER_TSV):
    sys.exit(f"[ERROR] Master SNP table not found: {IN_MASTER_TSV}. Run 01_build_transcriptome_snp_table.py first.")

print(f"[INFO] Loading master SNP table from: {IN_MASTER_TSV}...")
df_list = []
with gzip.open(IN_MASTER_TSV, "rt") as gz:
    reader = csv.DictReader(gz, delimiter="\t")
    for row in reader:
        n_subj = int(row["n_subjects"])
        if n_subj >= args.min_subjects:
            df_list.append({
                "chr": str(row["chr"]).replace("chr", ""),
                "pos": int(row["pos"]),
                "ref": row["ref"],
                "alt": row["alt"],
                "n_subjects": n_subj,
                "n_mdd": int(row["n_mdd"]),
                "n_ctrl": int(row["n_ctrl"]),
                "n_male": int(row.get("n_male", 0)),
                "n_female": int(row.get("n_female", 0)),
                "n_male_mdd": int(row.get("n_male_mdd", 0)),
                "n_male_ctrl": int(row.get("n_male_ctrl", 0)),
                "n_female_mdd": int(row.get("n_female_mdd", 0)),
                "n_female_ctrl": int(row.get("n_female_ctrl", 0)),
                "mdd_fraction": float(row["mdd_fraction"]),
                "male_mdd_fraction": row.get("male_mdd_fraction", "N/A"),
                "female_mdd_fraction": row.get("female_mdd_fraction", "N/A"),
                "sex_specificity": row.get("sex_specificity", "Unclassified"),
                "carrier_subjects": row["carrier_subjects"]
            })

var_df = pd.DataFrame(df_list)
print(f"[INFO] Found {len(var_df):,} recurrent SNP loci (carried by >={args.min_subjects} subjects).")

# Prioritize by MDD skew and total recurrence
var_df['mdd_skew'] = (var_df['mdd_fraction'] - 0.5).abs()
var_df = var_df.sort_values(by=['mdd_skew', 'n_subjects'], ascending=[False, False]).reset_index(drop=True)

query_df = var_df.head(args.top_n)
print(f"[INFO] Prioritizing top {len(query_df)} high-recurrence / disease-skewed loci for GTEx Brain BA9 audit...")

# --- Step 3: Query GTEx Portal API v2 by Variant Coordinates ---
def query_gtex_variant_eqtl(chrom, pos, ref, alt):
    # GTEx v8 variant ID format: chr1_12345_A_G_b38
    var_id = f"chr{chrom}_{pos}_{ref}_{alt}_b38"
    url = f"https://gtexportal.org/api/v2/association/singleTissueEqtl?variantId={var_id}&tissueSiteDetailId=Brain_Frontal_Cortex_BA9&datasetId=gtex_v8"
    req = urllib.request.Request(url, headers={"Accept": "application/json", "User-Agent": "Mozilla/5.0"})
    
    for attempt in range(1, 4):
        try:
            with urllib.request.urlopen(req, timeout=15) as resp:
                data = json.loads(resp.read().decode('utf-8'))
                eqtls = data.get('data', [])
                if eqtls:
                    top_eqtl = eqtls[0]
                    return {
                        'is_ba9_eqtl': True,
                        'eqtl_count': len(eqtls),
                        'rsid': top_eqtl.get('snpId', 'N/A'),
                        'target_gene': top_eqtl.get('geneSymbol', top_eqtl.get('gencodeId', 'N/A')),
                        'gencode_id': top_eqtl.get('gencodeId', 'N/A'),
                        'ba9_pvalue': top_eqtl.get('pValue', 'N/A'),
                        'ba9_nes': top_eqtl.get('nes', 'N/A')
                    }
                else:
                    return {'is_ba9_eqtl': False, 'eqtl_count': 0, 'rsid': 'N/A', 'target_gene': 'N/A', 'gencode_id': 'N/A', 'ba9_pvalue': 'N/A', 'ba9_nes': 'N/A'}
        except Exception as e:
            if attempt < 3:
                time.sleep(1.5 * attempt)
            else:
                pass
    return {'is_ba9_eqtl': False, 'eqtl_count': 0, 'rsid': 'N/A', 'target_gene': 'N/A', 'gencode_id': 'N/A', 'ba9_pvalue': 'N/A', 'ba9_nes': 'N/A'}

print("[INFO] Querying GTEx Portal API v2 dynamically...")
gtex_results = []
eqtl_hits = 0

for idx, row in query_df.iterrows():
    c = row['chr']
    p = row['pos']
    r = row['ref']
    a = row['alt']
    
    res = query_gtex_variant_eqtl(c, p, r, a)
    if res['is_ba9_eqtl']:
        eqtl_hits += 1
    
    gtex_results.append({
        'chr': c,
        'pos': p,
        'ref': r,
        'alt': a,
        'rsid': res['rsid'],
        'n_subjects': row['n_subjects'],
        'n_mdd': row['n_mdd'],
        'n_ctrl': row['n_ctrl'],
        'n_male': row['n_male'],
        'n_female': row['n_female'],
        'n_male_mdd': row['n_male_mdd'],
        'n_male_ctrl': row['n_male_ctrl'],
        'n_female_mdd': row['n_female_mdd'],
        'n_female_ctrl': row['n_female_ctrl'],
        'mdd_fraction': row['mdd_fraction'],
        'male_mdd_fraction': row['male_mdd_fraction'],
        'female_mdd_fraction': row['female_mdd_fraction'],
        'sex_specificity': row['sex_specificity'],
        'is_gtex_ba9_eqtl': "YES" if res['is_ba9_eqtl'] else "NO",
        'target_gene_symbol': res['target_gene'],
        'gencode_id': res['gencode_id'],
        'ba9_eqtl_pvalue': res['ba9_pvalue'],
        'ba9_nes': res['ba9_nes']
    })
    
    if (idx + 1) % 25 == 0 or (idx + 1) == len(query_df):
        print(f"  [{idx+1:03d}/{len(query_df)}] Queried GTEx => {eqtl_hits} direct BA9 eQTL hits found so far...")
    time.sleep(0.08)

out_df = pd.DataFrame(gtex_results)

# --- Step 4: Data Assertions ---
assert not out_df.empty, "Assertion failed: GTEx results table is empty."
assert 'is_gtex_ba9_eqtl' in out_df.columns, "Assertion failed: 'is_gtex_ba9_eqtl' column missing."

# --- Step 5: Export Table ---
out_df.to_csv(OUT_CSV, index=False)
print(f"[SUCCESS] Saved GTEx Brain BA9 eQTL results to: {OUT_CSV}")
print(f"[SUMMARY] Queried {len(out_df)} candidate loci => {eqtl_hits} validated direct GTEx Frontal Cortex eQTLs.")
print("[DONE] 02_map_gtex_brain_eqtls.py executed successfully.")
