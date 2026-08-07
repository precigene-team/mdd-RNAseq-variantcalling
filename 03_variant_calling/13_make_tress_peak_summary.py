#!/usr/bin/env python3
# [[03_variant_calling/13_make_tress_peak_summary.py]]
# Goal: Generate TRESS peak-count summary tables from the ASM input matrix.

# Input : _data/RData/03_variant_calling/asm_matrix.tsv
#         (one row per subject x consensus peak; n_peaks == 1 per row)
# Output: results/03_variant_calling/tables/tress_peak_summary.csv
#         results/03_variant_calling/tables/tress_peak_per_subject.csv

import csv
import os
import statistics
import sys
from collections import defaultdict

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.environ.get("PROJECT_ROOT") or os.path.normpath(
    os.path.join(SCRIPT_DIR, "..")
)
MATRIX = os.path.join(PROJECT_ROOT, "_data", "RData", "03_variant_calling", "asm_matrix.tsv")
OUT_TABLES = os.path.join(PROJECT_ROOT, "results", "03_variant_calling", "tables")
OUT_SUMMARY = os.path.join(OUT_TABLES, "tress_peak_summary.csv")
OUT_PER_SUBJECT = os.path.join(OUT_TABLES, "tress_peak_per_subject.csv")

if not os.path.isfile(MATRIX):
    print(f"Assertion Error: ASM matrix not found: {MATRIX}", file=sys.stderr)
    sys.exit(1)


def cohort_of(condition, sex):
    return f"{sex}_{condition}"


def main():
    os.makedirs(OUT_TABLES, exist_ok=True)

    with open(MATRIX, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")

        subject_peaks = defaultdict(int)
        subject_meta = {}
        subject_genes = defaultdict(set)
        cohort_subjects = defaultdict(set)
        cohort_loci = defaultdict(set)
        cohort_genes = defaultdict(set)
        all_loci = set()
        all_genes = set()

        for row in reader:
            subject = row["subject"]
            condition = row["condition"]
            sex = row["sex"]
            cohort = cohort_of(condition, sex)

            if subject not in subject_meta:
                subject_meta[subject] = (cohort, sex, condition)
            subject_peaks[subject] += 1
            subject_genes[subject].add(row["gene"])
            cohort_subjects[cohort].add(subject)
            cohort_loci[cohort].add(row["locus_id"])
            cohort_genes[cohort].add(row["gene"])
            all_loci.add(row["locus_id"])
            all_genes.add(row["gene"])

    # --- Per-subject table ---
    with open(OUT_PER_SUBJECT, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["subject", "cohort", "sex", "group", "n_peaks", "n_loci", "n_genes"])
        for subject in sorted(subject_peaks):
            cohort, sex, group = subject_meta[subject]
            writer.writerow([
                subject, cohort, sex, group,
                subject_peaks[subject], subject_peaks[subject], len(subject_genes[subject]),
            ])

    # --- Cohort summary table ---
    cohort_order = ["female_ctrl", "female_mdd", "male_ctrl", "male_mdd"]
    with open(OUT_SUMMARY, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow([
            "cohort", "n_subjects", "mean_peaks_subject", "sd_peaks",
            "median_peaks", "min_peaks", "max_peaks", "distinct_loci", "distinct_genes",
        ])
        for cohort in cohort_order:
            members = sorted(cohort_subjects[cohort])
            vals = [subject_peaks[s] for s in members]
            sd = statistics.stdev(vals) if len(vals) > 1 else 0.0
            writer.writerow([
                cohort,
                len(members),
                round(statistics.mean(vals), 1),
                sd,
                statistics.median(vals),
                min(vals),
                max(vals),
                len(cohort_loci[cohort]),
                len(cohort_genes[cohort]),
            ])
        writer.writerow(["All", len(subject_peaks), "", "", "", "", "", len(all_loci), len(all_genes)])

    print(f"  Per-subject table -> {OUT_PER_SUBJECT}")
    print(f"  Cohort summary   -> {OUT_SUMMARY}")


if __name__ == "__main__":
    main()