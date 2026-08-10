#!/usr/bin/env python3
# [[03_variant_calling/09_annotate_drach_disruption.py]]
# Goal: Annotate which overlapping SNPs disrupt the DRACH motif pattern.

# Steps:
#   1. Load the bedtools SNP × DRACH intersection BED.
#   2. Compute 0-based motif offset for each overlapping SNP.
#   3. Check whether the ALT allele is allowed at that IUPAC position.
#   4. Print summary statistics.
#   5. Write the fully annotated table to results/03_variant_calling/tables/09_drach_disrupting_snps.csv.

import csv
import os
import sys
import pandas as pd

# --- Dynamic Pathing ---
script_dir   = os.path.dirname(os.path.abspath(__file__))
project_root = os.environ.get("PROJECT_ROOT") or os.path.normpath(
    os.path.join(script_dir, "..")
)

intersect = os.path.join(project_root, "_data", "logs", "drach", "tmp", "snp_drach_intersect.bed")
out_dir   = os.path.join(project_root, "results", "03_variant_calling", "tables")
os.makedirs(out_dir, exist_ok=True)

out_table     = os.path.join(out_dir, "09_drach_disrupting_snps.csv")
summary_table = os.path.join(out_dir, "09_drach_summary.csv")

# --- Input Assertion ---
if not os.path.isfile(intersect):
    print(f"Assertion Error: Intersection BED not found: {intersect}", file=sys.stderr)
    print("Run 09_intersect_drach_snps.sh first (Steps 1-2).", file=sys.stderr)
    sys.exit(1)

# Column layout from bedtools intersect -wa -wb:
# 0: snp_chr, 1: snp_start(0-based), 2: snp_pos(1-based), 3: subject, 4: ref, 5: alt
# 6: drach_chr, 7: drach_start, 8: drach_end, 9: drach_name, 10: score, 11: strand, 12: methyl_pos
cols = [
    "snp_chr", "snp_start", "snp_pos", "subject", "ref", "alt",
    "drach_chr", "drach_start", "drach_end", "drach_name", "score", "strand", "methyl_pos"
]

df = pd.read_csv(intersect, sep="\t", header=None, names=cols)
df["snp_pos"]     = df["snp_pos"].astype(int)
df["drach_start"] = df["drach_start"].astype(int)
df["drach_end"]   = df["drach_end"].astype(int)

# Assertion: data loaded
if df.empty:
    print("Assertion Error: Intersection dataframe is empty.", file=sys.stderr)
    sys.exit(1)

# --- IUPAC Definition ---
IUPAC_ALLOWED = {
    0: set("AGT"), # D = A/G/T (non-C)
    1: set("AG"),  # R = A/G (Purine)
    2: set("A"),   # A = A (central m6A site)
    3: set("C"),   # C = C (Cytosine)
    4: set("ACT"), # H = A/C/T (non-G)
}

def is_allowed(offset: int, allele: str) -> bool:
    allowed = IUPAC_ALLOWED.get(offset)
    if allowed is None:
        return False
    return allele.upper() in allowed

# Calculate motif offset:
df["motif_offset"] = df.apply(
    lambda r: (r["snp_pos"] - 1) - r["drach_start"]
    if r["strand"] == "+"
    else r["drach_end"] - r["snp_pos"],
    axis=1
)

df["ref_allowed"] = df.apply(lambda r: is_allowed(r["motif_offset"], r["ref"]), axis=1)
df["alt_allowed"] = df.apply(lambda r: is_allowed(r["motif_offset"], r["alt"]), axis=1)

df["disruption_status"] = df["alt_allowed"].map(
    {True: "PRESERVES_DRACH", False: "DISRUPTS_DRACH"}
)

df["ov_key"] = (
    df["snp_chr"].astype(str) + ":" +
    df["snp_pos"].astype(str) + ":" +
    df["ref"] + ":" + df["alt"] + "__" +
    df["drach_name"]
)

# Output fully annotated table
df.to_csv(out_table, index=False)

subjects        = df["subject"].nunique()
rows            = len(df)
distinct_snps   = df.groupby(["snp_chr", "snp_pos", "ref", "alt"]).ngroups
distinct_motifs = df["drach_name"].nunique()
distinct_overlap= df["ov_key"].nunique()

disrupting  = df.loc[df["disruption_status"] == "DISRUPTS_DRACH", "ov_key"].nunique()
preserving  = df.loc[df["disruption_status"] == "PRESERVES_DRACH", "ov_key"].nunique()
both        = len(df.groupby("ov_key")["disruption_status"].nunique().loc[lambda s: s == 2])

summary_rows = [
    ("cohort_subjects", subjects, "number of unique subjects carrying >=1 overlapping SNP"),
    ("subject_x_snp_x_motif_het_rows", rows, "total per-subject het SNP x DRACH motif rows"),
    ("distinct_snp_sites_het", distinct_snps, "distinct (chr:pos:ref:alt) heterozygous sites"),
    ("distinct_drach_motif_instances", distinct_motifs, "distinct DRACH motif instances (bed_names)"),
    ("distinct_snp_x_motif_overlaps", distinct_overlap, "deduplicated (chr:pos:ref:alt + drach_name)"),
    ("disrupts_drach_overlaps", disrupting, "deduplicated overlaps where SNP DISRUPTS the DRACH motif"),
    ("preserves_drach_overlaps", preserving, "deduplicated overlaps where SNP PRESERVES the motif"),
    ("overlaps_both", both, "overlaps where a SNP disrupts one motif and preserves another"),
]

with open(summary_table, "w", newline="") as fh:
    writer = csv.writer(fh)
    writer.writerow(["metric", "value", "note"])
    writer.writerows(summary_rows)

print(f"  Annotated table -> {out_table}")
print(f"  Cohort summary table -> {summary_table}")

if not os.path.isfile(out_table) or os.path.getsize(out_table) == 0:
    print(f"Assertion Error: Output table empty or not written: {out_table}", file=sys.stderr)
    sys.exit(1)
