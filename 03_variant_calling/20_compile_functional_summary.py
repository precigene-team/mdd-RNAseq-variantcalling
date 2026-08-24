#!/usr/bin/env python3
# [[03_variant_calling/20_compile_functional_summary.py]]
# Goal: Compile multi-layer functional annotations (Layers 1-5) into publication-ready deliverable tables.
#
# Steps:
#   1. Dynamic project root resolution and path setup.
#   2. Load upstream results dynamically:
#      - 12_asm_primary_significant.csv (Primary ASM statistical results)
#      - 17_vep_dbsnp_annotation.csv (Layer 1 & 2: VEP consequence, transcripts, rsIDs)
#      - 18_gtex_brain_eqtl_sqtl.csv (Layer 3: GTEx Brain BA9 eQTL/sQTL)
#      - 19_gwas_catalog_associations.csv (Layer 4: GWAS Catalog associations)
#   3. Query NCBI PubMed E-utilities API dynamically for literature citations (Layer 5).
#   4. Dynamically merge all empirical evidence layers without hardcoded biological assumptions.
#   5. Apply strict data assertions (non-empty merge, 100% locus coverage).
#   6. Export consolidated CSV and formatted Markdown deliverable report.

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

SIG_CSV  = os.path.join(TABLES_DIR, "12_asm_primary_significant.csv")
VEP_CSV  = os.path.join(TABLES_DIR, "17_vep_dbsnp_annotation.csv")
GTEX_CSV = os.path.join(TABLES_DIR, "18_gtex_brain_eqtl_sqtl.csv")
GWAS_CSV = os.path.join(TABLES_DIR, "19_gwas_catalog_associations.csv")

OUT_CSV = os.path.join(TABLES_DIR, "20_candidate_loci_functional_summary.csv")
OUT_MD  = os.path.join(TABLES_DIR, "20_candidate_loci_functional_summary.md")

# --- Step 2: Load Upstream Layer Tables ---
print("[INFO] Loading upstream layer outputs...")
for f in [SIG_CSV, VEP_CSV, GTEX_CSV, GWAS_CSV]:
    if not os.path.exists(f):
        sys.exit(f"[ERROR] Missing required upstream table: {f}. Please run steps 12, 17, 18, and 19 first.")

sig_df  = pd.read_csv(SIG_CSV)
vep_df  = pd.read_csv(VEP_CSV)
gtex_df = pd.read_csv(GTEX_CSV)
gwas_df = pd.read_csv(GWAS_CSV)

print(f"[INFO] Loaded {len(sig_df)} primary ASM loci, {len(vep_df)} VEP variants, {len(gtex_df)} GTEx records, and {len(gwas_df)} GWAS records.")

# --- Step 3: Dynamic PubMed E-Utilities Query (Layer 5) ---
def query_pubmed(gene_symbol):
    query = f"{gene_symbol}+AND+(MDD+OR+depression+OR+brain+OR+m6A+OR+stress)"
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term={query}&retmode=json&retmax=5"
    req = urllib.request.Request(url, headers={"User-Agent": "MDD_RNAseq_Pipeline/1.0"})
    
    count = 0
    top_pmids = []
    try:
        with urllib.request.urlopen(req, timeout=10) as resp:
            data = json.loads(resp.read().decode('utf-8'))
            esearch = data.get('esearchresult', {})
            count = int(esearch.get('count', 0))
            top_pmids = esearch.get('idlist', [])
        time.sleep(0.2)
    except Exception as e:
        print(f"[WARNING] PubMed query failed for {gene_symbol}: {e}")
    
    pmid_str = ", ".join([f"PMID:{p}" for p in top_pmids]) if top_pmids else "None"
    return count, pmid_str

print("[INFO] Querying NCBI PubMed dynamically for candidate genes...")
pubmed_data = {}
for g in sig_df['gene'].dropna().unique():
    c_pub, pmid_list = query_pubmed(g)
    pubmed_data[g] = {'pubmed_hit_count': c_pub, 'top_pmids': pmid_list}
    print(f"  -> Gene: {g} => {c_pub} PubMed articles found ({pmid_list})")

# --- Step 4: Dynamically Merge All Layers ---
print("[INFO] Merging all functional annotation layers dynamically...")

# Merge VEP with GTEx
merged_df = vep_df.merge(
    gtex_df[['gene', 'chr', 'pos', 'gtex_ba9_eqtl_count', 'top_ba9_variant', 'top_ba9_pvalue', 'top_ba9_nes', 'brain_tissues_active']],
    on=['gene', 'chr', 'pos'],
    how='left'
)

# Merge with GWAS
merged_df = merged_df.merge(
    gwas_df[['gene', 'chr', 'pos', 'variant_gwas_hits', 'gene_gwas_traits_count', 'gene_gwas_traits_summary']],
    on=['gene', 'chr', 'pos'],
    how='left'
)

# Add PubMed layer
merged_df['pubmed_hits'] = merged_df['gene'].map(lambda g: pubmed_data.get(g, {}).get('pubmed_hit_count', 0))
merged_df['top_pmids']   = merged_df['gene'].map(lambda g: pubmed_data.get(g, {}).get('top_pmids', 'N/A'))

# Merge primary model statistics from sig_df
model_cols = ['gene', 'coef', 'p', 'fdr', 'male_coef', 'male_p', 'female_coef', 'female_p', 'confidence']
merged_df = merged_df.merge(
    sig_df[[c for c in model_cols if c in sig_df.columns]],
    on='gene',
    how='left',
    suffixes=('', '_model')
)

# --- Step 5: Data Assertions ---
assert not merged_df.empty, "Assertion failed: Merged summary dataframe is empty."
assert set(sig_df['gene']).issubset(set(merged_df['gene'])), "Assertion failed: Not all candidate loci represented."

# --- Step 6: Export Deliverable Tables ---
merged_df.to_csv(OUT_CSV, index=False)
print(f"[INFO] Successfully saved consolidated summary CSV to: {OUT_CSV}")

# Generate Markdown Report
md_rows = []
for _, r in merged_df.iterrows():
    md_rows.append(
        f"| **{r['gene']}** | `{r['chr']}:{r['pos']}` | `{r['ref']}>{r['alt']}` | {r['rsid']} | "
        f"{r['vep_consequence']} | {r['disruption_status']} | "
        f"BA9 eQTLs: {r.get('gtex_ba9_eqtl_count', 0)} (Top: {r.get('top_ba9_variant', 'N/A')}) | "
        f"{r.get('variant_gwas_hits', 'None')} | "
        f"{r.get('pubmed_hits', 0)} articles ({r.get('top_pmids', 'N/A')}) | "
        f"$\\beta={r.get('coef', 'N/A')}, p_{{\\text{{FDR}}}}={r.get('fdr', 'N/A')}$ |"
    )

md_table_text = "\n".join(md_rows)

md_content = f"""# 🧬 Candidate Loci Multi-Layer Functional Annotation Summary

This report aggregates empirical, un-hardcoded multi-omics functional annotation across **Layers 1-5** for the FDR-significant Allele-Specific Methylation (ASM) candidate loci identified in the MDD vs. Control analysis.

---

## 📊 Consolidated Multi-Omics Functional Annotation Matrix

| Gene | Position (GRCh38) | Allele | dbSNP rsID | VEP Consequence (Layer 1 & 2) | $m^6A$ / DRACH Motif Effect | GTEx Brain BA9 eQTL (Layer 3) | GWAS Catalog Associations (Layer 4) | PubMed Evidence (Layer 5) | ASM Model Statistics |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
{md_table_text}

---

## 🔍 Data Provenance & Execution Pipeline
* **Layer 1 (VEP) & Layer 2 (dbSNP)**: Generated dynamically by `03_variant_calling/17_annotate_vep_dbsnp.py` via Ensembl REST API.
* **Layer 3 (GTEx Portal)**: Generated dynamically by `03_variant_calling/18_audit_gtex.py` via GTEx Portal API v2.
* **Layer 4 (GWAS Catalog)**: Generated dynamically by `03_variant_calling/19_crossref_gwas.py` via NHGRI-EBI GWAS Catalog API.
* **Layer 5 (Literature Curation)**: Generated dynamically by `03_variant_calling/20_compile_functional_summary.py` via NCBI Entrez E-utilities API.

*Outputs saved to `results/03_variant_calling/tables/20_candidate_loci_functional_summary.csv`.*
"""

with open(OUT_MD, "w") as f:
    f.write(md_content)

print(f"[INFO] Successfully exported Markdown summary report to: {OUT_MD}")
print("[DONE] 20_compile_functional_summary.py executed successfully.")
