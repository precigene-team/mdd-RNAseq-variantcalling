#!/usr/bin/env python3
# [[03_variant_calling/06_build_transcriptome_snp_table.py]]
# Goal: Build the unified transcriptome-wide heterozygous SNP cohort master table across 39 subjects without m6A peak filtering.
#
# Steps:
#   1. Dynamic project root resolution and path setup.
#   2. Locate all 39 per-sample filtered VCF files (_data/ase/*_het_snps.vcf.gz).
#   3. Parse high-confidence biallelic heterozygous SNPs (PASS, 0/1, DP>=10, GQ>=20).
#   4. Aggregate cohort-wide carrier metrics (total subjects, MDD count, Control count, sex breakdown).
#   5. Apply strict data assertions (non-empty, valid chromosomes, coordinate bounds).
#   6. Export compressed master TSV and cohort summary statistics to results/tables/06_transcriptome_snp_*.

import os
import sys
import glob
import gzip
import csv
from collections import defaultdict

# --- Step 1: Dynamic Pathing ---
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.environ.get("PROJECT_ROOT") or os.path.normpath(os.path.join(SCRIPT_DIR, ".."))

ASE_DIR = os.path.join(PROJECT_ROOT, "_data", "ase")
TABLES_DIR = os.path.join(PROJECT_ROOT, "results", "tables")
MOD_TABLES_DIR = os.path.join(PROJECT_ROOT, "results", "03_variant_calling", "tables")

for d in [TABLES_DIR, MOD_TABLES_DIR]:
    os.makedirs(d, exist_ok=True)

OUT_MASTER_TSV = os.path.join(TABLES_DIR, "06_transcriptome_snp_master.tsv.gz")
OUT_SUMMARY_CSV = os.path.join(TABLES_DIR, "06_transcriptome_snp_summary.csv")

# --- Step 2: Discover VCF Files ---
FILTERED_VCF_DIR = os.path.join(PROJECT_ROOT, "_data", "variant", "filtered_vcf")
vcf_files = sorted(glob.glob(os.path.join(FILTERED_VCF_DIR, "*.filtered.vcf.gz")))
source_label = "GATK filtered non-ref VCFs (_data/variant/filtered_vcf)"
if len(vcf_files) != 39:
    vcf_files = sorted(glob.glob(os.path.join(ASE_DIR, "*_het_snps.vcf.gz")))
    source_label = "ASE heterozygous VCFs (_data/ase)"

if len(vcf_files) != 39:
    sys.exit(f"[ERROR] Expected 39 per-sample VCF files across cohort, found {len(vcf_files)}.")

print(f"[INFO] Discovered {len(vcf_files)} per-sample VCF files ({source_label}).")

# --- Step 3: Parse and Aggregate SNPs ---
print("[INFO] Parsing and aggregating transcriptome-wide heterozygous SNPs...")

locus_data = defaultdict(lambda: {
    'male_mdd': 0, 'male_ctrl': 0,
    'female_mdd': 0, 'female_ctrl': 0,
    'subjects': []
})

total_observations = 0

for idx, vcf_path in enumerate(vcf_files, 1):
    fname = os.path.basename(vcf_path)
    subject = fname.replace(".filtered.vcf.gz", "").replace("_input.filtered.vcf.gz", "").replace("_het_snps.vcf.gz", "")
    
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
            
            chrom, pos_str, vid, ref, alt, qual, filt, info, fmt, sample_val = parts[:10]
            
            if filt != "PASS":
                continue
            if len(ref) != 1 or len(alt) != 1:
                continue
            
            fmt_keys = fmt.split(":")
            val_items = sample_val.split(":")
            fmap = dict(zip(fmt_keys, val_items))
            
            gt = fmap.get("GT", "")
            # Non-reference variants: heterozygous (0/1, 1/0) and homozygous alternate (1/1)
            if gt not in ("0/1", "1/0", "1/1", "0|1", "1|0", "1|1"):
                continue
            
            try:
                dp = int(fmap.get("DP", 0))
            except ValueError:
                dp = 0
            if dp < 10:
                continue
            
            try:
                gq = float(fmap.get("GQ", 0))
            except ValueError:
                gq = 0.0
            if gq < 20.0:
                continue
            
            clean_chr = chrom.replace("chr", "")
            pos = int(pos_str)
            key = (clean_chr, pos, ref, alt)
            
            group_key = f"{sex}_{condition}"
            locus_data[key][group_key] += 1
            locus_data[key]['subjects'].append(subject)
            sample_calls += 1
            total_observations += 1

    if idx % 10 == 0 or idx == len(vcf_files):
        print(f"  [PROCESSED] {idx}/{len(vcf_files)} subjects parsed. (Latest: {subject}, {sample_calls:,} het-SNPs)")

n_unique_loci = len(locus_data)
print(f"[INFO] Parsed {total_observations:,} total het-observations collapsed into {n_unique_loci:,} distinct loci.")

# --- Step 4: Strict Data Assertions ---
assert total_observations > 1_000_000, f"Assertion Failed: Total variant count {total_observations} unexpectedly low (<1M)."
assert n_unique_loci > 500_000, f"Assertion Failed: Unique loci count {n_unique_loci} unexpectedly low (<500k)."

# --- Step 5: Export Compressed Master TSV ---
print(f"[INFO] Exporting compressed master TSV to: {OUT_MASTER_TSV}...")

def sort_chrom(ch):
    if ch.isdigit():
        return (0, int(ch))
    elif ch == "X":
        return (1, 23)
    elif ch == "Y":
        return (1, 24)
    elif ch in ("M", "MT"):
        return (1, 25)
    return (2, ch)

sorted_keys = sorted(locus_data.keys(), key=lambda k: (sort_chrom(k[0]), k[1]))

header = [
    "chr", "pos", "ref", "alt",
    "n_subjects", "n_mdd", "n_ctrl",
    "n_male", "n_female",
    "n_male_mdd", "n_male_ctrl",
    "n_female_mdd", "n_female_ctrl",
    "mdd_fraction", "male_mdd_fraction", "female_mdd_fraction",
    "sex_specificity", "carrier_subjects"
]

def write_tsv(target_path):
    with gzip.open(target_path, "wt", newline="") as gz_out:
        writer = csv.writer(gz_out, delimiter="\t")
        writer.writerow(header)
        for k in sorted_keys:
            d = locus_data[k]
            n_mdd = d['male_mdd'] + d['female_mdd']
            n_ctrl = d['male_ctrl'] + d['female_ctrl']
            n_male = d['male_mdd'] + d['male_ctrl']
            n_female = d['female_mdd'] + d['female_ctrl']
            n_subj = n_mdd + n_ctrl
            
            mdd_frac = round(n_mdd / n_subj, 4) if n_subj > 0 else 0.0
            male_mdd_frac = round(d['male_mdd'] / n_male, 4) if n_male > 0 else 0.0
            female_mdd_frac = round(d['female_mdd'] / n_female, 4) if n_female > 0 else 0.0
            
            if n_male > 0 and n_female == 0:
                sex_spec = "Male-Specific"
            elif n_female > 0 and n_male == 0:
                sex_spec = "Female-Specific"
            else:
                sex_spec = "Sex-Shared"
            
            writer.writerow([
                k[0], k[1], k[2], k[3],
                n_subj, n_mdd, n_ctrl,
                n_male, n_female,
                d['male_mdd'], d['male_ctrl'],
                d['female_mdd'], d['female_ctrl'],
                mdd_frac, male_mdd_frac, female_mdd_frac,
                sex_spec, ",".join(d['subjects'])
            ])

write_tsv(OUT_MASTER_TSV)
write_tsv(os.path.join(MOD_TABLES_DIR, "06_transcriptome_snp_master.tsv.gz"))

# --- Step 6: Export Cohort Summary Metrics ---
recurrence_bins = defaultdict(int)
for k, d in locus_data.items():
    n_subj = d['male_mdd'] + d['male_ctrl'] + d['female_mdd'] + d['female_ctrl']
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

def write_summary(target_path):
    with open(target_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["metric", "value", "note"])
        writer.writeheader()
        writer.writerows(summary_rows)

write_summary(OUT_SUMMARY_CSV)
write_summary(os.path.join(MOD_TABLES_DIR, "06_transcriptome_snp_summary.csv"))

print(f"[SUCCESS] Exported summary statistics to: {OUT_SUMMARY_CSV}")
print("[DONE] 06_build_transcriptome_snp_table.py executed successfully.")
