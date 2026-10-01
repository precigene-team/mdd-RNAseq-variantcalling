#!/usr/bin/env python3
# [[03_variant_calling/08_map_gtex_brain_eqtls.py]]
# Goal: Map functionally annotated candidate transcriptome SNPs against GTEx Brain Frontal Cortex (BA9) eQTLs.
#
# Steps:
#   1. Dynamic project root resolution and path setup.
#   2. Load functional candidates from 07_functional_annotation.csv.
#   3. Prioritize functional candidates (HIGH/MODERATE, valid hgvs_p, n_ctrl in {0, 1}).
#   4. Query GTEx Portal v2 API for top candidates.
#   5. Apply strict data assertions.
#   6. Save annotated results to results/03_variant_calling/tables/08_gtex_brain_eqtl_matches.csv.

import os
import sys
import argparse
import urllib.request
import urllib.parse
import json
import time
import pandas as pd
import numpy as np

# --- Step 1: Dynamic Pathing ---
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.environ.get("PROJECT_ROOT") or os.path.normpath(os.path.join(SCRIPT_DIR, ".."))

TABLES_DIR = os.path.join(PROJECT_ROOT, "results", "03_variant_calling", "tables")
os.makedirs(TABLES_DIR, exist_ok=True)

IN_FUNC_CSV = os.path.join(TABLES_DIR, "07_functional_annotation.csv")
OUT_CSV = os.path.join(TABLES_DIR, "08_gtex_brain_eqtl_matches.csv")

parser = argparse.ArgumentParser(description="Map functional SNPs against GTEx Brain eQTLs.")
parser.add_argument("--top-n", type=int, default=500, help="Max candidate variants to query (default: 500)")
args = parser.parse_args()

# --- Step 2: Load Functional Annotation Table ---
if not os.path.exists(IN_FUNC_CSV):
    sys.exit(f"[ERROR] Functional annotation table not found: {IN_FUNC_CSV}. Run 07_functional_annotation.py first.")

print(f"[INFO] Loading annotated functional variants from: {IN_FUNC_CSV}...")
df_func = pd.read_csv(IN_FUNC_CSV, low_memory=False)
print(f"[INFO] Loaded {len(df_func):,} total annotated recurrent SNPs.")

# --- Step 3: Prioritize Candidates for GTEx Screening ---
# Pre-filter: HIGH, MODERATE & LOW impact (protein-altering and synonymous/splice), and strict: n_ctrl <= 2 & n_mdd >= 8
df_func['has_hgvs_p'] = df_func['hgvs_p'].apply(
    lambda x: 1 if pd.notna(x) and str(x).strip() not in ('None', '.', 'N/A', '', 'nan') else 0
)
is_functional = df_func['impact'].str.upper().isin(['HIGH', 'MODERATE', 'LOW']) & (
    (df_func['has_hgvs_p'] == 1) | df_func['consequence'].str.contains('splice', case=False, na=False)
)
is_mdd_enriched = (df_func['n_ctrl'] <= 2) & (df_func['n_mdd'] >= 8)

df_filtered = df_func[is_functional & is_mdd_enriched].copy()

# Sort: n_ctrl ascending (0 prioritized before 1, 2), n_mdd descending, mdd_fraction descending
df_cand = df_filtered.sort_values(
    by=['n_ctrl', 'n_mdd', 'mdd_fraction', 'n_subjects'],
    ascending=[True, False, False, False]
).head(args.top_n).copy()

print(f"[INFO] Selected {len(df_cand):,} candidates (HIGH/MODERATE/LOW, n_ctrl <= 2 & n_mdd >= 8) for GTEx eQTL mapping.")
print(f"       Control carrier breakdown: {df_cand['n_ctrl'].value_counts().to_dict()}")
print(f"       Impact breakdown: {df_cand['impact'].value_counts().to_dict()}")
print(f"       Annotated coding/splice loci: {len(df_cand)}")
print(f"       Top MDD carrier counts: {df_cand['n_mdd'].head(5).tolist()}")

# --- Step 4: Query GTEx Portal API ---
print(f"[INFO] Querying GTEx Portal v2 API for candidate variants...")
TISSUE = "Brain_Frontal_Cortex_BA9"
API_BASE = "https://gtexportal.org/api/v2/association/singleTissueEqtl"

results = []
for idx, (_, r) in enumerate(df_cand.iterrows(), 1):
    chrom = f"chr{r['chrom']}"
    pos = int(r['pos'])
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
    top_gene = r.get('gene', 'None')
    gencode_id = r.get('gene_id', 'None')
    pval = np.nan
    nes = np.nan
    rsid = "."
    
    req = urllib.request.Request(url, headers={"Accept": "application/json", "User-Agent": "Mozilla/5.0"})
    try:
        with urllib.request.urlopen(req, timeout=5) as response:
            if response.status == 200:
                data = json.loads(response.read().decode('utf-8'))
                eqtls = data.get("data", [])
                if eqtls:
                    is_eqtl = "YES"
                    best_eqtl = min(eqtls, key=lambda x: x.get("pValue", 1.0))
                    pval = best_eqtl.get("pValue")
                    nes = best_eqtl.get("nes")
                    top_gene = best_eqtl.get("geneSymbol", top_gene)
                    gencode_id = best_eqtl.get("gencodeId", gencode_id)
                    rsid = best_eqtl.get("snpId", ".")
                    print(f"  [{idx:03d}/{len(df_cand)}] eQTL HIT: {top_gene:12s} ({rsid:11s}) p={pval:.2e}, NES={nes:+.2f}", flush=True)
    except Exception:
        pass
    
    rec = dict(r)
    rec.update({
        "gtex_variant_id": gtex_variant_id,
        "rsid": rsid if rsid != "." else r.get("rsid", "."),
        "is_gtex_ba9_eqtl": is_eqtl,
        "target_gene_symbol": top_gene,
        "gencode_id": gencode_id,
        "ba9_eqtl_pvalue": pval,
        "ba9_nes": nes
    })
    results.append(rec)
    
    if idx % 25 == 0 or idx == len(df_cand):
        n_hits = sum(1 for x in results if x['is_gtex_ba9_eqtl'] == 'YES')
        print(f"  [GTEx Query] {idx:03d}/{len(df_cand)} loci processed ({n_hits} eQTL hits found so far)...", flush=True)
    time.sleep(0.01)

df_out = pd.DataFrame(results)
df_out = df_out.drop(columns=['impact_rank', 'has_hgvs_p', 'mdd_skew'], errors='ignore')

# --- Step 5: Data Integrity Assertions ---
assert not df_out.empty, "Assertion Failed: Output table is empty."
assert len(df_out) == len(df_cand), f"Assertion Failed: Row count mismatch (input {len(df_cand)} vs output {len(df_out)})"

# --- Step 6: Export Results ---
df_out.to_csv(OUT_CSV, index=False)

print(f"[SUCCESS] Exported GTEx eQTL results to: {OUT_CSV}")
print(f"          eQTL positive rate: {(df_out['is_gtex_ba9_eqtl'] == 'YES').sum()} / {len(df_out)} loci.")
