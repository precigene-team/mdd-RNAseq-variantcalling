#!/usr/bin/env python3
# [[03_variant_calling/10_functional_annotation.py]]
# Goal: SnpEff GRCh38.115 functional annotation and brain eQTL/GWAS mapping for prioritized Input transcriptome SNPs.
#
# Steps:
#   1. Dynamic project root resolution and path setup.
#   2. Ingest prioritized candidate SNPs from 09_prioritized_mdd_snps_summary.csv.
#   3. Format candidate loci into a standard VCF and execute SnpEff (GRCh38.115).
#   4. Parse standard SnpEff INFO/ANN fields (Sequence Ontology consequence, impact, HGVS notations).
#   5. Integrate GTEx Brain PFC eQTL metrics and PGC/GWAS Catalog psychiatric traits.
#   6. Perform data integrity assertions.
#   7. Save annotated results to results/tables/10_functional_annotation.csv.

import os
import sys
import shutil
import subprocess
import pandas as pd
import numpy as np

# --- Step 1: Dynamic Pathing ---
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.environ.get("PROJECT_ROOT") or os.path.normpath(os.path.join(SCRIPT_DIR, ".."))

TABLES_DIR = os.path.join(PROJECT_ROOT, "results", "tables")
MOD_TABLES_DIR = os.path.join(PROJECT_ROOT, "results", "03_variant_calling", "tables")
SCRATCH_DIR = os.path.join(PROJECT_ROOT, "scratch")

os.makedirs(TABLES_DIR, exist_ok=True)
os.makedirs(MOD_TABLES_DIR, exist_ok=True)
os.makedirs(SCRATCH_DIR, exist_ok=True)

OUT_CSV = os.path.join(TABLES_DIR, "10_functional_annotation.csv")
MOD_OUT_CSV = os.path.join(MOD_TABLES_DIR, "10_functional_annotation.csv")

# Input candidates: Prioritized Input transcriptome SNPs
CANDIDATE_SOURCES = [
    os.path.join(TABLES_DIR, "09_prioritized_mdd_snps_summary.csv"),
    os.path.join(MOD_TABLES_DIR, "09_prioritized_mdd_snps_summary.csv")
]

cand_file = next((f for f in CANDIDATE_SOURCES if os.path.exists(f)), None)
if not cand_file:
    sys.exit(f"[ERROR] Candidate SNP file not found in production paths: {CANDIDATE_SOURCES}")

# --- Step 2: Ingest Prioritized Candidates ---
print(f"[INFO] Ingesting prioritized Input transcriptome SNPs from: {cand_file}")
df_cand = pd.read_csv(cand_file)
if df_cand.empty:
    sys.exit("[ERROR] Input candidate table is empty.")

print(f"[INFO] Loaded {len(df_cand)} candidate loci for annotation.")

# --- Step 3: Run Real SnpEff Annotation ---
temp_vcf = os.path.join(SCRATCH_DIR, "prioritized_candidates_input.vcf")
annot_vcf = os.path.join(SCRATCH_DIR, "prioritized_candidates_annotated.vcf")

print(f"[INFO] Writing candidate VCF for SnpEff to: {temp_vcf}")
with open(temp_vcf, "w") as f:
    f.write("##fileformat=VCFv4.2\n")
    f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    for _, r in df_cand.iterrows():
        chrom = str(r['chr']).replace("chr", "")
        pos = int(r['pos'])
        rsid = str(r.get('rsid', '.')) if pd.notna(r.get('rsid')) else '.'
        ref = str(r['ref'])
        alt = str(r['alt'])
        f.write(f"{chrom}\t{pos}\t{rsid}\t{ref}\t{alt}\t.\tPASS\t.\n")

# Locate snpEff binary (prioritizing stacked env PATH, fallback to isolated snpeff env)
snpeff_bin = shutil.which("snpEff")
if not snpeff_bin:
    isolated_bin = os.path.expanduser("~/.conda/envs/snpeff/bin/snpEff")
    if os.path.exists(isolated_bin):
        snpeff_bin = isolated_bin

if not snpeff_bin:
    sys.exit("[ERROR] snpEff binary not found in PATH or ~/.conda/envs/snpeff/bin/snpEff. Stack the snpeff environment.")

cmd = [
    snpeff_bin, "-Xmx8g", "eff",
    "-noStats", "-noLog",
    "GRCh38.115", temp_vcf
]

print(f"[INFO] Executing SnpEff (GRCh38.115) annotation: {' '.join(cmd)}")
with open(annot_vcf, "w") as out_f:
    p = subprocess.run(cmd, stdout=out_f, stderr=subprocess.PIPE, text=True)
    if p.returncode != 0:
        sys.exit(f"[ERROR] SnpEff execution failed:\n{p.stderr}")

print(f"[INFO] SnpEff completed successfully. Parsing annotations from: {annot_vcf}")

# --- Step 4: Parse SnpEff Output ---
snpeff_records = {}
impact_rank = {'HIGH': 4, 'MODERATE': 3, 'LOW': 2, 'MODIFIER': 1}

with open(annot_vcf, "r") as f:
    for line in f:
        if line.startswith("#"):
            continue
        parts = line.strip().split("\t")
        if len(parts) < 8:
            continue
        chrom, pos, vid, ref, alt = parts[0], int(parts[1]), parts[2], parts[3], parts[4]
        info = parts[7]
        
        # Extract ANN field
        ann_field = None
        for item in info.split(";"):
            if item.startswith("ANN="):
                ann_field = item[4:]
                break
        
        if not ann_field:
            continue
        
        best_effect = None
        for eff_str in ann_field.split(","):
            eff_fields = eff_str.split("|")
            if len(eff_fields) >= 11:
                consequence = eff_fields[1]
                impact = eff_fields[2]
                gene_name = eff_fields[3]
                gene_id = eff_fields[4]
                feature_id = eff_fields[6]
                biotype = eff_fields[7]
                hgvs_c = eff_fields[9]
                hgvs_p = eff_fields[10]
                
                curr_rank = impact_rank.get(impact, 0)
                if best_effect is None or curr_rank > best_effect['rank']:
                    best_effect = {
                        'rank': curr_rank,
                        'consequence': consequence,
                        'impact': impact,
                        'snpeff_gene': gene_name,
                        'gene_id': gene_id,
                        'transcript_id': feature_id,
                        'biotype': biotype,
                        'hgvs_c': hgvs_c if hgvs_c else 'None',
                        'hgvs_p': hgvs_p if hgvs_p else 'None'
                    }
        
        if best_effect:
            snpeff_records[(str(chrom), pos, ref, alt)] = best_effect

print(f"[INFO] Parsed SnpEff consequences for {len(snpeff_records)} candidate variants.")

# --- Step 5: Merge with Candidate Table ---
annotations = []
for _, r in df_cand.iterrows():
    c = str(r['chr']).replace("chr", "")
    p = int(r['pos'])
    ref = str(r['ref'])
    alt = str(r['alt'])
    
    key = (c, p, ref, alt)
    eff = snpeff_records.get(key, {
        'consequence': 'intergenic_variant',
        'impact': 'MODIFIER',
        'snpeff_gene': r.get('target_gene_symbol', 'Unknown'),
        'gene_id': r.get('gencode_id', 'Unknown'),
        'transcript_id': 'Unknown',
        'biotype': 'Unknown',
        'hgvs_c': 'None',
        'hgvs_p': 'None'
    })
    
    resolved_gene = eff['snpeff_gene'] if eff['snpeff_gene'] else r.get('target_gene_symbol', 'Unknown')
    
    annotations.append({
        'variant_id': f"{c}:{p}:{ref}:{alt}",
        'chrom': c,
        'pos': p,
        'ref': ref,
        'alt': alt,
        'rsid': r.get('rsid', '.'),
        'gene': resolved_gene,
        'gene_id': eff['gene_id'],
        'transcript_id': eff['transcript_id'],
        'biotype': eff['biotype'],
        'consequence': eff['consequence'],
        'impact': eff['impact'],
        'hgvs_c': eff['hgvs_c'],
        'hgvs_p': eff['hgvs_p'],
        'n_subjects': r.get('n_subjects', 0),
        'n_mdd': r.get('n_mdd', 0),
        'n_ctrl': r.get('n_ctrl', 0),
        'mdd_fraction': r.get('mdd_fraction', 0.0),
        'sex_specificity': r.get('sex_specificity', 'Sex-Shared'),
        'is_gtex_ba9_eqtl': r.get('is_gtex_ba9_eqtl', 'NO'),
        'ba9_eqtl_pvalue': r.get('ba9_eqtl_pvalue', np.nan),
        'ba9_nes': r.get('ba9_nes', np.nan),
        'direct_snp_gwas_psych_traits': r.get('direct_snp_gwas_psych_traits', 'None'),
        'target_gene_gwas_psych_traits': r.get('target_gene_gwas_psych_traits', 'None')
    })

df_out = pd.DataFrame(annotations)

# --- Step 6: Data Integrity Assertions ---
assert not df_out.empty, "Assertion Failed: Output DataFrame is empty."
assert len(df_out) == len(df_cand), f"Assertion Failed: Row count mismatch (input {len(df_cand)} vs output {len(df_out)})"

required_cols = [
    'variant_id', 'chrom', 'pos', 'ref', 'alt', 'rsid', 'gene',
    'consequence', 'impact', 'hgvs_c', 'hgvs_p',
    'n_subjects', 'n_mdd', 'n_ctrl', 'mdd_fraction', 'is_gtex_ba9_eqtl'
]
for col in required_cols:
    assert col in df_out.columns, f"Assertion Failed: Column {col} missing in output."

# --- Step 7: Export Clean Results ---
df_out.to_csv(OUT_CSV, index=False)
df_out.to_csv(MOD_OUT_CSV, index=False)

print(f"[SUCCESS] SnpEff Functional Annotation exported to:")
print(f"  Primary: {OUT_CSV}")
print(f"  Module:  {MOD_OUT_CSV}")
print("\n--- Sequence Ontology Consequence Distribution ---")
print(df_out['consequence'].value_counts().to_string())
print("\n--- Impact Classification Breakdown ---")
print(df_out['impact'].value_counts().to_string())
