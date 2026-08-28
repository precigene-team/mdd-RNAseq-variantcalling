#!/usr/bin/env python3
# [[03_variant_calling/alt/05_annotate_dbsnp_vep.py]]
# Goal: Sequential Ensembl Variant Effect Predictor (VEP) REST API querying to map molecular consequences, biotypes, impact, and ClinVar status.
#
# Steps:
#   1. Dynamic project root resolution and path setup.
#   2. Load prioritized candidates from 04_prioritized_mdd_snps_summary.csv.
#   3. Query Ensembl VEP REST API sequentially (clean, reliable, no multi-threading).
#   4. Extract molecular consequences, biotypes, impact terms, and ClinVar status.
#   5. Export Deliverable Table 1: results/03_variant_calling/alt/tables/05_prioritized_candidates_vep_annotated.csv
#   6. Extract Deliverable Table 2: results/03_variant_calling/alt/tables/05_recurrent_transcriptome_snps.tsv.gz (All recurrent loci >=8 subjects).

import os
import sys
import gzip
import csv
import json
import time
import urllib.request
import pandas as pd

# --- Step 1: Dynamic Pathing ---
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.environ.get("PROJECT_ROOT") or os.path.normpath(os.path.join(SCRIPT_DIR, "../.."))

TABLES_DIR = os.path.join(PROJECT_ROOT, "results", "03_variant_calling", "alt", "tables")
os.makedirs(TABLES_DIR, exist_ok=True)

IN_PRIORITIZED_CSV   = os.path.join(TABLES_DIR, "04_prioritized_mdd_snps_summary.csv")
IN_MASTER_TSV        = os.path.join(TABLES_DIR, "01_transcriptome_snp_master.tsv.gz")

OUT_PRIORITIZED_CSV  = os.path.join(TABLES_DIR, "05_prioritized_candidates_vep_annotated.csv")
OUT_DBSNP_CSV        = os.path.join(TABLES_DIR, "05_dbsnp_vep_annotations.csv")
OUT_RECURRENT_TSV_GZ = os.path.join(TABLES_DIR, "05_recurrent_transcriptome_snps.tsv.gz")

if not os.path.exists(IN_PRIORITIZED_CSV):
    sys.exit(f"[ERROR] Required candidate table not found: {IN_PRIORITIZED_CSV}. Run 04_prioritize_mdd_candidate_snps.py first.")

if not os.path.exists(IN_MASTER_TSV):
    sys.exit(f"[ERROR] Master 800k SNP table not found: {IN_MASTER_TSV}. Run 01_build_transcriptome_snp_table.py first.")

# ==============================================================================
# --- TABLE 1: VEP FUNCTIONAL ANNOTATION OF PRIORITIZED SUBSET (04 RESULTS) ---
# ==============================================================================
print(f"[INFO] Loading prioritized candidate variants from: {IN_PRIORITIZED_CSV}")
cand_df = pd.read_csv(IN_PRIORITIZED_CSV)
print(f"[INFO] Querying Ensembl VEP REST API for {len(cand_df)} candidates sequentially...")

def fetch_vep_single(chrom, pos, ref, alt):
    var_str = f"{chrom}:{pos}-{pos}:1/{alt}"
    url = f"https://rest.ensembl.org/vep/human/region/{chrom}:{pos}-{pos}:1/{alt}"
    req = urllib.request.Request(url, headers={"Accept": "application/json", "User-Agent": "Mozilla/5.0"})
    
    for attempt in range(1, 4):
        try:
            with urllib.request.urlopen(req, timeout=10) as resp:
                data = json.loads(resp.read().decode('utf-8'))
                if data and isinstance(data, list):
                    return data[0]
        except Exception:
            if attempt < 3:
                time.sleep(0.5 * attempt)
    return {}

vep_results = {}
for idx, row in cand_df.iterrows():
    c   = str(row['chr']).replace('chr', '')
    pos = int(row['pos'])
    ref = str(row['ref'])
    alt = str(row['alt'])
    var_str = f"{c}:{pos}-{pos}:1/{alt}"
    
    v_data = fetch_vep_single(c, pos, ref, alt)
    if v_data:
        vep_results[var_str] = v_data
    
    print(f"  [{idx+1:02d}/{len(cand_df)}] Processed VEP for {var_str:22s} => consequence: {v_data.get('most_severe_consequence', 'N/A')}")
    time.sleep(0.05)

print(f"[INFO] Successfully retrieved VEP functional consequences for {len(vep_results)}/{len(cand_df)} variants.")

annotated_list = []
for idx, row in cand_df.iterrows():
    c   = str(row['chr']).replace('chr', '')
    pos = int(row['pos'])
    alt = str(row['alt'])
    var_str = f"{c}:{pos}-{pos}:1/{alt}"
    
    vep_item = vep_results.get(var_str, {})
    consequence = vep_item.get('most_severe_consequence', 'N/A')
    
    transcripts = vep_item.get('transcript_consequences', [])
    biotype = "N/A"
    impact = "N/A"
    amino_acids = "N/A"
    codons = "N/A"
    if transcripts:
        t0 = transcripts[0]
        biotype = t0.get('biotype', 'N/A')
        impact = t0.get('impact', 'N/A')
        amino_acids = t0.get('amino_acids', 'N/A')
        codons = t0.get('codons', 'N/A')
    
    clinvar_sig = "N/A"
    colocated = vep_item.get('colocated_variants', [])
    for cv in colocated:
        cln = cv.get('clin_sig')
        if cln:
            if isinstance(cln, list):
                clinvar_sig = "; ".join(cln)
            else:
                clinvar_sig = str(cln)
            break
    
    rec = dict(row)
    rec['vep_consequence'] = consequence
    rec['vep_impact'] = impact
    rec['vep_biotype'] = biotype
    rec['vep_amino_acids'] = amino_acids
    rec['vep_codons'] = codons
    rec['clinvar_clinical_significance'] = clinvar_sig
    annotated_list.append(rec)

out_subset_df = pd.DataFrame(annotated_list)

col_order = [
    'chr', 'pos', 'ref', 'alt', 'rsid',
    'target_gene_symbol', 'gencode_id',
    'vep_consequence', 'vep_impact', 'vep_biotype', 'vep_amino_acids', 'vep_codons', 'clinvar_clinical_significance',
    'n_subjects', 'n_mdd', 'n_ctrl',
    'n_male', 'n_female',
    'n_male_mdd', 'n_male_ctrl',
    'n_female_mdd', 'n_female_ctrl',
    'mdd_fraction', 'male_mdd_fraction', 'female_mdd_fraction',
    'sex_specificity',
    'is_gtex_ba9_eqtl', 'ba9_eqtl_pvalue', 'ba9_nes',
    'direct_snp_gwas_psych_traits',
    'target_gene_gwas_psych_traits'
]

final_cols = [c for c in col_order if c in out_subset_df.columns]
out_subset_df = out_subset_df[final_cols].copy()

# Export Subset Tables
out_subset_df.to_csv(OUT_PRIORITIZED_CSV, index=False)
out_subset_df.to_csv(OUT_DBSNP_CSV, index=False)
print(f"[SUCCESS] Table 1 Exported: {len(out_subset_df)} prioritized candidates saved to: {OUT_PRIORITIZED_CSV}")

# ==============================================================================
# --- TABLE 2: RECURRENT TRANSCRIPTOME SNP POOL (>=8 SUBJECTS) ---
# ==============================================================================
print(f"[INFO] Extracting recurrent transcriptome SNP candidate pool (>=8 subjects) from master TSV...")

recurrent_count = 0
with gzip.open(IN_MASTER_TSV, "rt") as gz_in, gzip.open(OUT_RECURRENT_TSV_GZ, "wt", newline="") as gz_out:
    reader = csv.reader(gz_in, delimiter="\t")
    writer = csv.writer(gz_out, delimiter="\t")
    
    header = next(reader)
    writer.writerow(header)
    n_subj_idx = header.index("n_subjects")
    
    for row in reader:
        if int(row[n_subj_idx]) >= 8:
            writer.writerow(row)
            recurrent_count += 1

print(f"[SUCCESS] Table 2 Exported: {recurrent_count:,} recurrent SNP loci (>=8 subjects) saved to: {OUT_RECURRENT_TSV_GZ}")
print("[DONE] 05_annotate_dbsnp_vep.py executed successfully.")
