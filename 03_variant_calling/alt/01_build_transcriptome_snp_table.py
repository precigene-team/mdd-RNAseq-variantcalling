#!/usr/bin/env python3
# [[03_variant_calling/alt/01_build_transcriptome_snp_table.py]]
# Goal: Build the unified transcriptome-wide heterozygous SNP cohort master table across 39 subjects without m6A peak filtering.
#
# Steps:
#   1. Dynamic project root resolution and path setup.
#   2. Locate all 39 per-sample filtered VCF files (_data/ase/*_het_snps.vcf.gz).
#   3. Parse high-confidence biallelic heterozygous SNPs (PASS, 0/1, DP>=10, GQ>=20).
#   4. Aggregate cohort-wide carrier metrics (total subjects, MDD count, Control count, sex breakdown).
#   5. Apply strict data assertions (non-empty, valid chromosomes, coordinate bounds).
#   6. Export compressed master TSV and cohort summary statistics.

import os
import sys
import glob
import gzip
import csv
from collections import defaultdict

# --- Step 1: Dynamic Pathing ---
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.environ.get("PROJECT_ROOT") or os.path.normpath(os.path.join(SCRIPT_DIR, "../.."))

ASE_DIR = os.path.join(PROJECT_ROOT, "_data", "ase")
OUT_DIR = os.path.join(PROJECT_ROOT, "results", "03_variant_calling", "alt", "tables")
os.makedirs(OUT_DIR, exist_ok=True)

OUT_MASTER_TSV = os.path.join(OUT_DIR, "01_transcriptome_snp_master.tsv.gz")
OUT_SUMMARY_CSV = os.path.join(OUT_DIR, "01_transcriptome_snp_summary.csv")

# --- Step 2: Discover VCF Files ---
vcf_files = sorted(glob.glob(os.path.join(ASE_DIR, "*_het_snps.vcf.gz")))
if len(vcf_files) != 39:
    sys.exit(f"[ERROR] Expected 39 per-sample VCF files in {ASE_DIR}, found {len(vcf_files)}.")

print(f"[INFO] Discovered {len(vcf_files)} per-sample VCF files across the cohort.")

# --- Step 3: Parse and Aggregate SNPs ---
print("[INFO] Parsing and aggregating transcriptome-wide heterozygous SNPs...")

# locus -> dict of metrics
# key: (chr, pos, ref, alt)
locus_data = defaultdict(lambda: {
    'male_mdd': 0, 'male_ctrl': 0,
    'female_mdd': 0, 'female_ctrl': 0,
    'subjects': []
})

total_observations = 0

for idx, vcf_path in enumerate(vcf_files, 1):
    fname = os.path.basename(vcf_path)
    subject = fname.replace("_het_snps.vcf.gz", "")
    
    # Parse metadata from subject ID
    condition = "mdd" if "_mdd_" in subject else "ctrl"
    sex = "male" if subject.startswith("male_") else "female"
    
    sample_calls = 0
    with gzip.open(vcf_path, "rt") as gz:
        for line in gz:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 10:
                continue
            
            chrom, pos, var_id, ref, alt, qual, flt, info, fmt, sample_val = parts[:10]
            
            # Require standard canonical chromosomes or filter non-standard
            if flt not in ["PASS", "."]:
                continue
            
            fmt_keys = fmt.split(":")
            fmt_vals = sample_val.split(":")
            fmt_dict = dict(zip(fmt_keys, fmt_vals))
            
            gt = fmt_dict.get("GT", "")
            dp = int(fmt_dict.get("DP", 0)) if fmt_dict.get("DP", "").isdigit() else 0
            gq = int(fmt_dict.get("GQ", 0)) if fmt_dict.get("GQ", "").isdigit() else 0
            
            # Biallelic het SNP QC threshold
            if gt in ["0/1", "0|1", "1|0"] and dp >= 10 and gq >= 20 and len(ref) == 1 and len(alt) == 1:
                key = (chrom, int(pos), ref, alt)
                rec = locus_data[key]
                if sex == "male":
                    if condition == "mdd":
                        rec['male_mdd'] += 1
                    else:
                        rec['male_ctrl'] += 1
                else:
                    if condition == "mdd":
                        rec['female_mdd'] += 1
                    else:
                        rec['female_ctrl'] += 1
                
                rec['subjects'].append(subject)
                sample_calls += 1
                total_observations += 1
    
    print(f"  [{idx:02d}/39] Processed {subject:20s}: {sample_calls:,} pass het-SNPs")

n_unique_loci = len(locus_data)
print(f"[INFO] Cohort aggregation complete: {total_observations:,} observations across {n_unique_loci:,} unique SNP loci.")

# --- Step 4: Apply Strict Data Assertions ---
assert n_unique_loci > 500000, f"Assertion failed: Unique loci count ({n_unique_loci}) is below expected minimum (500k)."
assert total_observations > 1000000, f"Assertion failed: Total observations ({total_observations}) is below expected minimum (1M)."

# --- Step 5: Export Master Table ---
print(f"[INFO] Writing compressed master TSV to: {OUT_MASTER_TSV}...")

# Sort loci by chromosome and position
def chrom_sort_key(k):
    c = str(k[0]).replace("chr", "")
    if c == "X": val = 23
    elif c == "Y": val = 24
    elif c == "M" or c == "MT": val = 25
    elif c.isdigit(): val = int(c)
    else: val = 99
    return (val, k[1], k[2], k[3])

sorted_keys = sorted(locus_data.keys(), key=chrom_sort_key)

with gzip.open(OUT_MASTER_TSV, "wt", newline="") as gz_out:
    writer = csv.writer(gz_out, delimiter="\t")
    writer.writerow([
        "chr", "pos", "ref", "alt", 
        "n_subjects", "n_mdd", "n_ctrl", 
        "n_male", "n_female",
        "n_male_mdd", "n_male_ctrl",
        "n_female_mdd", "n_female_ctrl",
        "mdd_fraction", "male_mdd_fraction", "female_mdd_fraction",
        "sex_specificity", "carrier_subjects"
    ])
    
    for k in sorted_keys:
        d = locus_data[k]
        n_mdd = d['male_mdd'] + d['female_mdd']
        n_ctrl = d['male_ctrl'] + d['female_ctrl']
        n_subj = n_mdd + n_ctrl
        n_male = d['male_mdd'] + d['male_ctrl']
        n_female = d['female_mdd'] + d['female_ctrl']
        
        mdd_frac = round(n_mdd / n_subj, 4) if n_subj > 0 else 0.0
        male_mdd_frac = round(d['male_mdd'] / n_male, 4) if n_male > 0 else "N/A"
        female_mdd_frac = round(d['female_mdd'] / n_female, 4) if n_female > 0 else "N/A"
        
        if n_male == 0 and n_female > 0:
            sex_spec = "Female-Specific (100% Female)"
        elif n_female == 0 and n_male > 0:
            sex_spec = "Male-Specific (100% Male)"
        elif n_male > 0 and n_female > 0:
            if male_mdd_frac != "N/A" and female_mdd_frac != "N/A":
                if male_mdd_frac >= 0.80 and female_mdd_frac <= 0.30:
                    sex_spec = "Male-Driven MDD"
                elif female_mdd_frac >= 0.80 and male_mdd_frac <= 0.30:
                    sex_spec = "Female-Driven MDD"
                else:
                    sex_spec = "Sex-Shared"
            else:
                sex_spec = "Sex-Shared"
        else:
            sex_spec = "Unclassified"
        
        writer.writerow([
            k[0], k[1], k[2], k[3],
            n_subj, n_mdd, n_ctrl,
            n_male, n_female,
            d['male_mdd'], d['male_ctrl'],
            d['female_mdd'], d['female_ctrl'],
            mdd_frac, male_mdd_frac, female_mdd_frac,
            sex_spec, ",".join(d['subjects'])
        ])

print(f"[SUCCESS] Exported {len(sorted_keys):,} SNP loci to {OUT_MASTER_TSV}.")

# --- Step 6: Export Cohort Summary Metrics ---
recurrence_bins = defaultdict(int)
chr_counts = defaultdict(int)

for k, d in locus_data.items():
    n_subj = d['male_mdd'] + d['male_ctrl'] + d['female_mdd'] + d['female_ctrl']
    chr_counts[k[0]] += 1
    if n_subj == 1:
        recurrence_bins["Private (1 Subject)"] += 1
    elif 2 <= n_subj <= 4:
        recurrence_bins["Low Recurrence (2-4 Subjects)"] += 1
    elif 5 <= n_subj <= 7:
        recurrence_bins["Moderate Recurrence (5-7 Subjects)"] += 1
    elif n_subj >= 8:
        recurrence_bins["High Recurrence (>=8 Subjects)"] += 1

summary_rows = [
    {"metric": "total_subjects", "value": 39, "note": "20 MDD; 19 Control (20 Male; 19 Female)"},
    {"metric": "total_het_observations", "value": total_observations, "note": "Sum of per-subject pass het-SNPs"},
    {"metric": "unique_snp_loci", "value": n_unique_loci, "note": "Distinct (chr; pos; ref; alt) loci"},
    {"metric": "mean_het_snps_per_subject", "value": round(total_observations / 39, 1), "note": "Average het calls per individual"},
    {"metric": "loci_private_1_subject", "value": recurrence_bins["Private (1 Subject)"], "note": f"{recurrence_bins['Private (1 Subject)']/n_unique_loci*100:.2f}% of all loci"},
    {"metric": "loci_low_recurrence_2_4", "value": recurrence_bins["Low Recurrence (2-4 Subjects)"], "note": f"{recurrence_bins['Low Recurrence (2-4 Subjects)']/n_unique_loci*100:.2f}% of all loci"},
    {"metric": "loci_moderate_recurrence_5_7", "value": recurrence_bins["Moderate Recurrence (5-7 Subjects)"], "note": f"{recurrence_bins['Moderate Recurrence (5-7 Subjects)']/n_unique_loci*100:.2f}% of all loci"},
    {"metric": "loci_high_recurrence_gte_8", "value": recurrence_bins["High Recurrence (>=8 Subjects)"], "note": f"{recurrence_bins['High Recurrence (>=8 Subjects)']/n_unique_loci*100:.2f}% of all loci"}
]

with open(OUT_SUMMARY_CSV, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=["metric", "value", "note"])
    writer.writeheader()
    writer.writerows(summary_rows)

print(f"[SUCCESS] Exported summary statistics to: {OUT_SUMMARY_CSV}")
print("[DONE] 01_build_transcriptome_snp_table.py executed successfully.")
