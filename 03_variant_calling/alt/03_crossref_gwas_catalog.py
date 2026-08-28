#!/usr/bin/env python3
# [[03_variant_calling/alt/03_crossref_gwas_catalog.py]]
# Goal: Sequential dual-tier GWAS cross-referencing: (1) Direct SNP rsID associations and (2) Target gene-level disease traits.
#
# Steps:
#   1. Dynamic project root resolution and path setup.
#   2. Load prioritized candidates from 02_gtex_brain_eqtl_matches.csv.
#   3. Query direct variant rsIDs sequentially from GWAS Catalog.
#   4. Query unique target genes sequentially from GWAS Catalog with psychiatric keyword filtering.
#   5. Map results to all candidate loci and export results/03_variant_calling/alt/tables/03_gwas_catalog_neuro_matches.csv.

import os
import sys
import json
import time
import urllib.request
import urllib.parse
import pandas as pd

# --- Step 1: Dynamic Pathing ---
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.environ.get("PROJECT_ROOT") or os.path.normpath(os.path.join(SCRIPT_DIR, "../.."))

TABLES_DIR = os.path.join(PROJECT_ROOT, "results", "03_variant_calling", "alt", "tables")
os.makedirs(TABLES_DIR, exist_ok=True)

IN_GTEX_CSV = os.path.join(TABLES_DIR, "02_gtex_brain_eqtl_matches.csv")
OUT_CSV = os.path.join(TABLES_DIR, "03_gwas_catalog_neuro_matches.csv")

if not os.path.exists(IN_GTEX_CSV):
    sys.exit(f"[ERROR] Required input file not found: {IN_GTEX_CSV}. Run 02_map_gtex_brain_eqtls.py first.")

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
    url = f"https://www.ebi.ac.uk/gwas/rest/api/singleNucleotidePolymorphisms/{rsid}/studies"
    req = urllib.request.Request(url, headers={"Accept": "application/json", "User-Agent": "Mozilla/5.0"})
    all_traits = []
    psych_traits = []
    try:
        with urllib.request.urlopen(req, timeout=5) as resp:
            data = json.loads(resp.read().decode('utf-8'))
            studies = data.get('_embedded', {}).get('studies', [])
            for st in studies:
                tr = st.get('diseaseTrait', {}).get('trait')
                if tr and tr not in all_traits:
                    all_traits.append(tr)
                    if any(k in tr.lower() for k in PSYCH_KEYWORDS):
                        psych_traits.append(tr)
    except Exception:
        pass
    
    direct_all = "; ".join(all_traits[:3]) if all_traits else "None cataloged"
    direct_psych = "; ".join(psych_traits[:4]) if psych_traits else "None detected"
    return direct_all, direct_psych

snp_cache = {}
for idx, rs in enumerate(unique_rsids, 1):
    d_all, d_psych = fetch_single_rsid(rs)
    snp_cache[rs] = (d_all, d_psych)
    if d_psych != "None detected":
        print(f"  [{idx:02d}/{len(unique_rsids)}] Direct SNP {rs:12s} => {d_psych}")
    time.sleep(0.04)

# --- Step 3: Gene-Level GWAS Query ---
unique_genes = [g for g in var_df['target_gene_symbol'].dropna().unique() if g not in ['N/A', '', 'None']]
print(f"[INFO] Querying {len(unique_genes)} unique candidate genes sequentially...")

def fetch_single_gene(gene_symbol):
    url = f"https://www.ebi.ac.uk/gwas/rest/api/singleNucleotidePolymorphisms/search/findByGene?geneName={urllib.parse.quote(gene_symbol)}&size=50"
    req = urllib.request.Request(url, headers={"Accept": "application/json", "User-Agent": "Mozilla/5.0"})
    all_traits = []
    psych_traits = []
    try:
        with urllib.request.urlopen(req, timeout=6) as resp:
            data = json.loads(resp.read().decode('utf-8'))
            snps = data.get('_embedded', {}).get('singleNucleotidePolymorphisms', [])
            for s in snps[:20]:
                st_href = s.get('_links', {}).get('studies', {}).get('href')
                if st_href:
                    try:
                        st_req = urllib.request.Request(st_href, headers={"Accept": "application/json", "User-Agent": "Mozilla/5.0"})
                        with urllib.request.urlopen(st_req, timeout=4) as st_resp:
                            st_data = json.loads(st_resp.read().decode('utf-8'))
                            for st in st_data.get('_embedded', {}).get('studies', []):
                                tr = st.get('diseaseTrait', {}).get('trait')
                                if tr and tr not in all_traits:
                                    all_traits.append(tr)
                                    if any(k in tr.lower() for k in PSYCH_KEYWORDS):
                                        psych_traits.append(tr)
                    except Exception:
                        pass
    except Exception:
        pass
    
    # Sort psych_traits so core depression/mood/suicidality terms appear first
    def trait_priority(t):
        t_low = t.lower()
        if any(k in t_low for k in ['depress', 'mood', 'suicid', 'affective']):
            return 0
        return 1
    
    sorted_psych = sorted(psych_traits, key=trait_priority)
    psych_summary = "; ".join(sorted_psych[:4]) if sorted_psych else "None detected"
    general_summary = "; ".join(all_traits[:3]) if all_traits else "None cataloged"
    return len(all_traits), len(sorted_psych), psych_summary, general_summary

gene_cache = {}
for idx, g in enumerate(unique_genes, 1):
    tot_cnt, psych_cnt, psych_sum, gen_sum = fetch_single_gene(g)
    gene_cache[g] = {
        'total_traits': tot_cnt,
        'psych_traits_count': psych_cnt,
        'psych_summary': psych_sum,
        'general_summary': gen_sum
    }
    if psych_cnt > 0:
        print(f"  [{idx:02d}/{len(unique_genes)}] Candidate Gene {g:12s} => {psych_sum}")
    time.sleep(0.04)

# --- Step 4: Map Results to Candidate Table ---
print("[INFO] Consolidating sequential annotations to master candidate table...")
results = []
for _, row in var_df.iterrows():
    rsid = row.get('rsid', 'N/A')
    d_all, d_psych = snp_cache.get(rsid, ("None cataloged", "None detected"))
    
    gene = row.get('target_gene_symbol', 'N/A')
    g_cached = gene_cache.get(gene, {
        'total_traits': 0,
        'psych_traits_count': 0,
        'psych_summary': 'None detected',
        'general_summary': 'None cataloged'
    })
    
    rec = dict(row)
    rec['direct_snp_gwas_psych_traits'] = d_psych
    rec['direct_snp_gwas_all_traits'] = d_all
    rec['target_gene_gwas_psych_traits'] = g_cached['psych_summary']
    rec['target_gene_gwas_all_traits'] = g_cached['general_summary']
    results.append(rec)

out_df = pd.DataFrame(results)

# --- Step 5: Export Table ---
out_df.to_csv(OUT_CSV, index=False)
print(f"[SUCCESS] Saved dual-layer GWAS results ({len(out_df)} loci) sequentially to: {OUT_CSV}")
print("[DONE] 03_crossref_gwas_catalog.py executed successfully.")
