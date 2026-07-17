#!/usr/bin/env python3
# [[03_variant_calling/10_annotate_drach_disruption.py]]
# Goal: Annotate which overlapping SNPs disrupt the DRACH motif pattern.

# Steps:
#   1. Load the bedtools SNP × DRACH intersection BED.
#   2. Compute 0-based motif offset for each overlapping SNP.
#   3. Check whether the ALT allele is allowed at that IUPAC position.
#   4. Print summary statistics.
#   5. Write the fully annotated table to results/tables/13_drach_disrupting_snps.csv.

import os
import sys
import pandas as pd

# --- Dynamic Pathing ---
script_dir   = os.path.dirname(os.path.abspath(__file__))
project_root = os.environ.get("PROJECT_ROOT") or os.path.normpath(
    os.path.join(script_dir, "..")
)

intersect = os.path.join(project_root, "_data", "logs", "drach", "tmp", "snp_drach_intersect.bed")
out_table = os.path.join(project_root, "results", "tables", "13_drach_disrupting_snps.csv")

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
    print("Assertion Error: Intersection BED is empty. No SNPs overlap DRACH sites.")
    sys.exit(1)

# --- Compute motif offset (0-indexed position within the 5-mer) ---
# snp_pos is 1-based; drach_start is 0-based → offset = (snp_pos - 1) - drach_start
df["motif_offset"] = (df["snp_pos"] - 1) - df["drach_start"]

# DRACH IUPAC: D=[AGT], R=[AG], A=[A], C=[C], H=[ACT]
DRACH_IUPAC = [
    frozenset("AGT"),   # position 0: D
    frozenset("AG"),    # position 1: R
    frozenset("A"),     # position 2: A ← m6A methylation site
    frozenset("C"),     # position 3: C
    frozenset("ACT"),   # position 4: H
]

def check_disruption(row):
    alt = str(row["alt"]).upper()
    if len(alt) != 1:
        return "INDEL_NOT_SCORED"
    offset = int(row["motif_offset"])
    if offset < 0 or offset >= 5:
        return "OUT_OF_RANGE"
    return "DISRUPTS_DRACH" if alt not in DRACH_IUPAC[offset] else "PRESERVES_DRACH"

df["disruption_status"] = df.apply(check_disruption, axis=1)

# --- Summary ---
total       = len(df)
disrupting  = (df["disruption_status"] == "DISRUPTS_DRACH").sum()
preserving  = (df["disruption_status"] == "PRESERVES_DRACH").sum()
methyl_hits = (df["motif_offset"] == 2).sum()
methyl_disrupt = ((df["motif_offset"] == 2) & (df["disruption_status"] == "DISRUPTS_DRACH")).sum()

print(f"  Total SNP-DRACH overlaps:              {total}")
print(f"  Motif-disrupting variants:             {disrupting}")
print(f"  Motif-preserving variants:             {preserving}")
print(f"  SNPs at methylation A (offset=2):      {methyl_hits}")
print(f"  Disrupting at methylation A:           {methyl_disrupt}")

# --- Write output ---
os.makedirs(os.path.dirname(out_table), exist_ok=True)
out_cols = [
    "snp_chr", "snp_pos", "subject", "ref", "alt",
    "drach_start", "drach_end", "drach_name", "strand",
    "motif_offset", "disruption_status", "methyl_pos"
]
df[out_cols].to_csv(out_table, index=False)

# --- Output Assertion ---
if not os.path.isfile(out_table) or os.path.getsize(out_table) == 0:
    print(f"Assertion Error: Output table empty or not written: {out_table}", file=sys.stderr)
    sys.exit(1)

print(f"  Annotated table → {out_table}")
