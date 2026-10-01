#!/usr/bin/env python3
# [[03_variant_calling/07_functional_annotation.py]]
# Goal: SnpEff GRCh38.115 functional annotation across recurrent transcriptome-wide SNPs (>=8 subjects).
#
# Steps:
#   1. Dynamic project root resolution and path setup.
#   2. Ingest recurrent SNPs (>=8 subjects) from 06_transcriptome_snp_master.tsv.gz.
#   3. Format candidate loci into a standard VCF and execute SnpEff (GRCh38.115).
#   4. Parse standard SnpEff INFO/ANN fields (Sequence Ontology consequence, impact, HGVS notations).
#   5. Prioritize transcript consequences with explicit protein annotations (hgvs_p).
#   6. Perform data integrity assertions.
#   7. Save annotated results to results/03_variant_calling/tables/07_functional_annotation.csv.

import os
import sys
import gzip
import csv
import shutil
import subprocess
import pandas as pd
import numpy as np

# --- Step 1: Dynamic Pathing ---
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.environ.get("PROJECT_ROOT") or os.path.normpath(os.path.join(SCRIPT_DIR, ".."))

TABLES_DIR = os.path.join(PROJECT_ROOT, "results", "03_variant_calling", "tables")
SCRATCH_DIR = os.path.join(PROJECT_ROOT, "scratch")

os.makedirs(TABLES_DIR, exist_ok=True)
os.makedirs(SCRATCH_DIR, exist_ok=True)

IN_MASTER_TSV = os.path.join(TABLES_DIR, "06_transcriptome_snp_master.tsv.gz")
OUT_CSV = os.path.join(TABLES_DIR, "07_functional_annotation.csv")

if not os.path.exists(IN_MASTER_TSV):
    sys.exit(f"[ERROR] Master SNP table not found: {IN_MASTER_TSV}. Run 06_build_transcriptome_snp_table.py first.")

# --- Step 2: Ingest Recurrent SNPs ---
print(f"[INFO] Ingesting recurrent SNPs (>=8 subjects) from: {IN_MASTER_TSV}")
recurrent_snps = []
with gzip.open(IN_MASTER_TSV, "rt") as gz:
    reader = csv.DictReader(gz, delimiter="\t")
    for r in reader:
        n_subj = int(r["n_subjects"])
        if n_subj >= 8:
            recurrent_snps.append({
                "chrom": str(r["chr"]).replace("chr", ""),
                "pos": int(r["pos"]),
                "ref": str(r["ref"]),
                "alt": str(r["alt"]),
                "n_subjects": n_subj,
                "n_het": int(r.get("n_het", n_subj)),
                "n_hom_alt": int(r.get("n_hom_alt", 0)),
                "n_mdd": int(r["n_mdd"]),
                "n_ctrl": int(r["n_ctrl"]),
                "n_male": int(r.get("n_male", 0)),
                "n_female": int(r.get("n_female", 0)),
                "n_male_mdd": int(r.get("n_male_mdd", 0)),
                "n_male_ctrl": int(r.get("n_male_ctrl", 0)),
                "n_female_mdd": int(r.get("n_female_mdd", 0)),
                "n_female_ctrl": int(r.get("n_female_ctrl", 0)),
                "mdd_fraction": float(r["mdd_fraction"]),
                "male_mdd_fraction": float(r.get("male_mdd_fraction", 0.0)),
                "female_mdd_fraction": float(r.get("female_mdd_fraction", 0.0)),
                "sex_specificity": r.get("sex_specificity", "Sex-Shared"),
                "carrier_subjects": r.get("carrier_subjects", "")
            })

df_recurrent = pd.DataFrame(recurrent_snps)
if df_recurrent.empty:
    sys.exit("[ERROR] No recurrent SNPs found meeting threshold (>=8 subjects).")

print(f"[INFO] Loaded {len(df_recurrent):,} recurrent candidate loci for SnpEff annotation.")

# --- Step 3: Run Real SnpEff Annotation ---
temp_vcf = os.path.join(SCRATCH_DIR, "recurrent_candidates_input.vcf")
annot_vcf = os.path.join(SCRATCH_DIR, "recurrent_candidates_annotated.vcf")

print(f"[INFO] Writing candidate VCF for SnpEff to: {temp_vcf}")
with open(temp_vcf, "w") as f:
    f.write("##fileformat=VCFv4.2\n")
    f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    for _, r in df_recurrent.iterrows():
        f.write(f"{r['chrom']}\t{r['pos']}\t.\t{r['ref']}\t{r['alt']}\t.\tPASS\t.\n")

# Locate snpEff binary with Java 21 environment binding
isolated_dir = os.path.expanduser("~/.conda/envs/snpeff/bin")
isolated_bin = os.path.join(isolated_dir, "snpEff")

run_env = os.environ.copy()
if os.path.exists(isolated_bin):
    snpeff_bin = isolated_bin
    run_env["PATH"] = f"{isolated_dir}:{run_env.get('PATH', '')}"
    run_env["JAVA_HOME"] = os.path.expanduser("~/.conda/envs/snpeff")
else:
    snpeff_bin = shutil.which("snpEff")

if not snpeff_bin:
    sys.exit("[ERROR] snpEff binary not found in PATH or ~/.conda/envs/snpeff/bin/snpEff. Stack the snpeff environment.")

cmd = [
    snpeff_bin, "-Xmx8g", "eff",
    "-noStats", "-noLog",
    "GRCh38.115", temp_vcf
]

print(f"[INFO] Executing SnpEff (GRCh38.115) annotation: {' '.join(cmd)}")
with open(annot_vcf, "w") as out_f:
    p = subprocess.run(cmd, stdout=out_f, stderr=subprocess.PIPE, text=True, env=run_env)
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
        chrom, pos, ref, alt = parts[0], int(parts[1]), parts[3], parts[4]
        info = parts[7]
        
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
                has_p = bool(hgvs_p and hgvs_p not in ('.', 'None', '', 'N/A'))
                if (best_effect is None or 
                    curr_rank > best_effect['rank'] or 
                    (curr_rank == best_effect['rank'] and has_p and not best_effect['has_p'])):
                    best_effect = {
                        'rank': curr_rank,
                        'has_p': has_p,
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

print(f"[INFO] Parsed SnpEff consequences for {len(snpeff_records):,} candidate variants.")

# --- Step 5: Merge Annotations with Recurrent Loci ---
annotations = []
for _, r in df_recurrent.iterrows():
    c = str(r['chrom']).replace("chr", "")
    p = int(r['pos'])
    ref = str(r['ref'])
    alt = str(r['alt'])
    
    key = (c, p, ref, alt)
    eff = snpeff_records.get(key, {
        'consequence': 'intergenic_variant',
        'impact': 'MODIFIER',
        'snpeff_gene': 'Unknown',
        'gene_id': 'Unknown',
        'transcript_id': 'Unknown',
        'biotype': 'Unknown',
        'hgvs_c': 'None',
        'hgvs_p': 'None'
    })
    
    annotations.append({
        'variant_id': f"{c}:{p}:{ref}:{alt}",
        'chrom': c,
        'pos': p,
        'ref': ref,
        'alt': alt,
        'gene': eff['snpeff_gene'] if eff['snpeff_gene'] else 'Unknown',
        'gene_id': eff['gene_id'],
        'transcript_id': eff['transcript_id'],
        'biotype': eff['biotype'],
        'consequence': eff['consequence'],
        'impact': eff['impact'],
        'hgvs_c': eff['hgvs_c'],
        'hgvs_p': eff['hgvs_p'],
        'n_subjects': r['n_subjects'],
        'n_het': r['n_het'],
        'n_hom_alt': r['n_hom_alt'],
        'n_mdd': r['n_mdd'],
        'n_ctrl': r['n_ctrl'],
        'n_male': r['n_male'],
        'n_female': r['n_female'],
        'n_male_mdd': r['n_male_mdd'],
        'n_male_ctrl': r['n_male_ctrl'],
        'n_female_mdd': r['n_female_mdd'],
        'n_female_ctrl': r['n_female_ctrl'],
        'mdd_fraction': r['mdd_fraction'],
        'male_mdd_fraction': r['male_mdd_fraction'],
        'female_mdd_fraction': r['female_mdd_fraction'],
        'sex_specificity': r['sex_specificity'],
        'carrier_subjects': r['carrier_subjects']
    })

df_out = pd.DataFrame(annotations)

# --- Step 6: Multi-Tier Priority Sorting ---
impact_order = {'HIGH': 4, 'MODERATE': 3, 'LOW': 2, 'MODIFIER': 1}
df_out['impact_rank'] = df_out['impact'].map(lambda x: impact_order.get(str(x).upper(), 0))
df_out['has_hgvs_p'] = df_out['hgvs_p'].apply(
    lambda x: 1 if pd.notna(x) and str(x) not in ('None', '.', 'N/A', '') else 0
)
df_out['mdd_skew'] = (df_out['mdd_fraction'] - 0.5).abs()

df_out = df_out.sort_values(
    by=['impact_rank', 'has_hgvs_p', 'mdd_skew', 'n_subjects'],
    ascending=[False, False, False, False]
).reset_index(drop=True)

df_out = df_out.drop(columns=['impact_rank', 'has_hgvs_p', 'mdd_skew'])

# --- Step 7: Data Integrity Assertions ---
assert not df_out.empty, "Assertion Failed: Output DataFrame is empty."
assert len(df_out) == len(df_recurrent), f"Assertion Failed: Row count mismatch (input {len(df_recurrent)} vs output {len(df_out)})"

required_cols = [
    'variant_id', 'chrom', 'pos', 'ref', 'alt', 'gene',
    'consequence', 'impact', 'hgvs_c', 'hgvs_p',
    'n_subjects', 'n_mdd', 'n_ctrl', 'mdd_fraction'
]
for col in required_cols:
    assert col in df_out.columns, f"Assertion Failed: Column {col} missing in output."

# --- Step 8: Export Clean Results ---
df_out.to_csv(OUT_CSV, index=False)

print(f"[SUCCESS] Prioritized SnpEff Functional Annotation exported to: {OUT_CSV}")
print("\n--- Sequence Ontology Consequence Distribution ---")
print(df_out['consequence'].value_counts().head(15).to_string())
print("\n--- Impact Classification Breakdown ---")
print(df_out['impact'].value_counts().to_string())
print(f"\n--- Top 10 Prioritized Functional Variants ---")
print(df_out[['gene', 'consequence', 'impact', 'hgvs_p', 'n_subjects', 'mdd_fraction']].head(10).to_string(index=False))
