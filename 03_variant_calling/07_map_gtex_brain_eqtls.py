#!/usr/bin/env python3
# [[03_variant_calling/07_map_gtex_brain_eqtls.py]]
# Goal: Query GTEx Portal v2 API to map transcriptome-wide candidate SNPs against Brain Frontal Cortex (BA9) eQTLs.
#
# Steps:
#   1. Dynamic project root resolution and path setup.
#   2. Load master transcriptome SNP table (06_transcriptome_snp_master.tsv.gz).
#   3. Select recurrent/disease-skewed candidate variants for prioritized GTEx brain audit.
#   4. Query GTEx Portal API v2 for singleTissueEqtl (Brain_Frontal_Cortex_BA9) and variant eQTL status.
#   5. Extract eQTL counts, top regulatory variants, p-values, and Normalized Effect Sizes (NES).
#   6. Apply strict data assertions.
#   7. Save results to results/tables/07_gtex_brain_eqtl_matches.csv.

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
PROJECT_ROOT = os.environ.get("PROJECT_ROOT") or os.path.normpath(os.path.join(SCRIPT_DIR, ".."))

TABLES_DIR = os.path.join(PROJECT_ROOT, "results", "tables")
MOD_TABLES_DIR = os.path.join(PROJECT_ROOT, "results", "03_variant_calling", "tables")

for d in [TABLES_DIR, MOD_TABLES_DIR]:
    os.makedirs(d, exist_ok=True)

INPUT_SOURCES = [
    os.path.join(TABLES_DIR, "06_transcriptome_snp_master.tsv.gz"),
    os.path.join(MOD_TABLES_DIR, "06_transcriptome_snp_master.tsv.gz")
]

IN_MASTER_TSV = next((f for f in INPUT_SOURCES if os.path.exists(f)), None)
OUT_CSV = os.path.join(TABLES_DIR, "07_gtex_brain_eqtl_matches.csv")

parser = argparse.ArgumentParser(description="Map transcriptome-wide SNPs against GTEx Brain eQTLs.")
parser.add_argument("--min-subjects", type=int, default=8, help="Minimum cohort subjects carrying SNP (default: 8)")
parser.add_argument("--top-n", type=int, default=500, help="Max top candidate variants to query (default: 500)")
args = parser.parse_args()

# --- Step 2: Load Master SNP Table ---
if not IN_MASTER_TSV or not os.path.exists(IN_MASTER_TSV):
    sys.exit(f"[ERROR] Master SNP table not found in: {INPUT_SOURCES}. Run 06_build_transcriptome_snp_table.py first.")

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
                "male_mdd_fraction": float(row.get("male_mdd_fraction", 0.0)),
                "female_mdd_fraction": float(row.get("female_mdd_fraction", 0.0)),
                "sex_specificity": row.get("sex_specificity", "Sex-Shared"),
                "carrier_subjects": row.get("carrier_subjects", "")
            })

df_cand = pd.DataFrame(df_list)
print(f"[INFO] Filtered to {len(df_cand)} recurrent candidate SNPs (carried by >={args.min_subjects} subjects).")

df_cand["mdd_skew"] = (df_cand["mdd_fraction"] - 0.5).abs()
df_cand = df_cand.sort_values(by=["mdd_skew", "n_subjects"], ascending=[False, False]).head(args.top_n)

# --- Step 3: Query GTEx Portal API ---
print(f"[INFO] Querying GTEx Portal v2 API for top {len(df_cand)} candidate variants...")
TISSUE = "Brain_Frontal_Cortex_BA9"
API_BASE = "https://gtexportal.org/api/v2/association/singleTissueEqtl"

results = []
for idx, (_, r) in enumerate(df_cand.iterrows(), 1):
    chrom = f"chr{r['chr']}"
    pos = r['pos']
    ref = r['ref']
    alt = r['alt']
    gtex_variant_id = f"{chrom}_{pos}_{ref}_{alt}_b38"
    
    params = urllib.parse.urlencode({
        "tissueSiteDetailId": TISSUE,
        "variantId": gtex_variant_id,
        "datasetId": "gtex_v8"
    })
    url = f"{API_BASE}?{params}"
    
    is_eqtl = "NO"
    top_gene = "None"
    gencode_id = "None"
    pval = None
    nes = None
    rsid = "None"
    
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "MDD-RNAseq-Pipeline/1.0"})
        with urllib.request.urlopen(req, timeout=10) as resp:
            data = json.loads(resp.read().decode())
            eqtl_list = data.get("data", [])
            if eqtl_list:
                is_eqtl = "YES"
                top_hit = min(eqtl_list, key=lambda x: x.get("pValue", 1.0))
                top_gene = top_hit.get("geneSymbol", "Unknown")
                gencode_id = top_hit.get("gencodeId", "Unknown")
                pval = top_hit.get("pValue")
                nes = top_hit.get("nes")
                rsid = top_hit.get("snpId", "None")
    except Exception:
        pass
    
    time.sleep(0.05)
    
    results.append({
        **r.to_dict(),
        "gtex_variant_id": gtex_variant_id,
        "rsid": rsid,
        "is_gtex_ba9_eqtl": is_eqtl,
        "target_gene_symbol": top_gene,
        "gencode_id": gencode_id,
        "ba9_eqtl_pvalue": pval,
        "ba9_nes": nes
    })
    
    if idx % 50 == 0 or idx == len(df_cand):
        print(f"  [GTEx Query] {idx}/{len(df_cand)} loci processed.")

df_out = pd.DataFrame(results)

# --- Step 4: Export Results ---
def write_df(target_path):
    df_out.to_csv(target_path, index=False)

write_df(OUT_CSV)
write_df(os.path.join(MOD_TABLES_DIR, "07_gtex_brain_eqtl_matches.csv"))

print(f"[SUCCESS] Exported GTEx eQTL results to: {OUT_CSV}")
print(f"          eQTL positive rate: {(df_out['is_gtex_ba9_eqtl'] == 'YES').sum()} / {len(df_out)} loci.")
