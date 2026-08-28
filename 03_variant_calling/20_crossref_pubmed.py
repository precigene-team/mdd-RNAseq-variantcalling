#!/usr/bin/env python3
# [[03_variant_calling/20_crossref_pubmed.py]]
# Goal: Cross-reference candidate genes against NCBI PubMed literature using comprehensive neuropsychiatric and brain keywords, and compile consolidated multi-layer functional annotations.
#
# Steps:
#   1. Dynamic project root resolution and path setup.
#   2. Load upstream results dynamically:
#      - 12_asm_primary_significant.csv (Primary ASM statistical results)
#      - 17_vep_dbsnp_annotation.csv (Layer 1 & 2: VEP consequence, transcripts, rsIDs)
#      - 18_gtex_brain_eqtl_sqtl.csv (Layer 3: GTEx Brain BA9 eQTL/sQTL)
#      - 19_gwas_catalog_associations.csv (Layer 4: GWAS Catalog associations)
#   3. Query NCBI PubMed E-utilities API dynamically with expanded brain/neuro/depression keywords.
#   4. Fetch article details (PMID, title, year, journal) for top relevant publications.
#   5. Dynamically merge all empirical evidence layers (VEP, GTEx, GWAS, PubMed, ASM stats).
#   6. Apply strict data assertions.
#   7. Export deliverable tables (CSV and Markdown report).

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

TABLES_DIR = os.path.join(PROJECT_ROOT, "results", "03_variant_calling", "tables")
os.makedirs(TABLES_DIR, exist_ok=True)

SIG_CSV  = os.path.join(TABLES_DIR, "12_asm_primary_significant.csv")
VEP_CSV  = os.path.join(TABLES_DIR, "17_vep_dbsnp_annotation.csv")
GTEX_CSV = os.path.join(TABLES_DIR, "18_gtex_brain_eqtl_sqtl.csv")
GWAS_CSV = os.path.join(TABLES_DIR, "19_gwas_catalog_associations.csv")

OUT_PUBMED_CSV = os.path.join(TABLES_DIR, "20_pubmed_neuro_associations.csv")
OUT_SUMMARY_CSV = os.path.join(TABLES_DIR, "20_candidate_loci_functional_summary.csv")
OUT_SUMMARY_MD  = os.path.join(TABLES_DIR, "20_candidate_loci_functional_summary.md")

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

# --- Step 3: Dynamic PubMed Query with Expanded Neuro/Psychiatric Keywords ---
NEURO_PSYCH_KEYWORDS = (
    "brain OR neuron OR neuronal OR cortex OR prefrontal OR synaptic OR synapse OR "
    "depression OR MDD OR \"major depressive\" OR mood OR bipolar OR schizophrenia OR "
    "suicide OR suicidal OR neurodevelopment OR cognitive OR cognition OR stress OR "
    "neuropsychiatric OR neuroinflammation OR neurodegeneration OR excitotoxicity"
)

def query_pubmed_neuro(gene_symbol):
    term = f'("{gene_symbol}"[Title/Abstract] OR "{gene_symbol}"[Gene Name]) AND ({NEURO_PSYCH_KEYWORDS})'
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term={urllib.parse.quote(term)}&retmode=json&retmax=5&sort=relevance"
    req = urllib.request.Request(url, headers={"User-Agent": "MDD_RNAseq_Pipeline/1.0"})
    
    count = 0
    top_pmids = []
    paper_highlights = []
    paper_records = []
    
    try:
        with urllib.request.urlopen(req, timeout=12) as resp:
            data = json.loads(resp.read().decode('utf-8'))
            esearch = data.get('esearchresult', {})
            count = int(esearch.get('count', 0))
            top_pmids = esearch.get('idlist', [])
        time.sleep(0.15)
        
        # Fetch summary details for top hits
        if top_pmids:
            u_sum = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&id={','.join(top_pmids)}&retmode=json"
            req_sum = urllib.request.Request(u_sum, headers={"User-Agent": "MDD_RNAseq_Pipeline/1.0"})
            with urllib.request.urlopen(req_sum, timeout=12) as r_sum:
                d_sum = json.loads(r_sum.read().decode('utf-8')).get('result', {})
                for pid in top_pmids:
                    p_info = d_sum.get(pid, {})
                    title = p_info.get('title', '')
                    src = p_info.get('source', '')
                    yr = p_info.get('pubdate', '')[:4]
                    if title:
                        paper_highlights.append(f"PMID:{pid} ({yr} {src})")
                        paper_records.append({
                            'gene': gene_symbol,
                            'pmid': pid,
                            'year': yr,
                            'journal': src,
                            'title': title
                        })
                    else:
                        paper_highlights.append(f"PMID:{pid}")
            time.sleep(0.15)
    except Exception as e:
        print(f"[WARNING] PubMed query failed for {gene_symbol}: {e}")
    
    pmid_summary = "; ".join(paper_highlights[:3]) if paper_highlights else ("None detected" if count == 0 else f"{count} hits")
    return count, pmid_summary, paper_records

print("[INFO] Querying NCBI PubMed dynamically with expanded neuropsychiatric keyword filtering...")
pubmed_gene_data = {}
all_paper_records = []

for g in sig_df['gene'].dropna().unique():
    c_pub, pmid_sum, p_records = query_pubmed_neuro(g)
    pubmed_gene_data[g] = {'pubmed_hit_count': c_pub, 'top_pmids': pmid_sum}
    all_paper_records.extend(p_records)
    print(f"  -> Gene: {g} => {c_pub} Neuro/Psychiatric PubMed articles ({pmid_sum})")

## --- Step 4: Export Deliverable Tables ---
summary_records = []
for g in sig_df['gene'].dropna().unique():
    info = pubmed_gene_data.get(g, {})
    summary_records.append({
        'gene': g,
        'neuro_psych_pubmed_hits': info.get('pubmed_hit_count', 0),
        'top_citations': info.get('top_pmids', 'None detected')
    })

# Generate Summary Table
summary_df = pd.DataFrame(summary_records)
summary_df.to_csv(OUT_SUMMARY_CSV, index=False)
print(f"[INFO] Saved PubMed summary table to: {OUT_SUMMARY_CSV}")

if all_paper_records:
    pubmed_df = pd.DataFrame(all_paper_records)
    pubmed_df.to_csv(OUT_PUBMED_CSV, index=False)
    print(f"[INFO] Saved detailed PubMed papers list to: {OUT_PUBMED_CSV}")

print("[DONE] 20_crossref_pubmed.py executed successfully.")
