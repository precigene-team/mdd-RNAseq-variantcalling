#!/usr/bin/env python3
# [[03_variant_calling/11_build_burden_matrix.py]]
# Goal: Build the donor x variant carrier matrix consumed by the burden steps.
#       DATA MUNGING ONLY -- no statistics, no network access.
#
# Inputs:
#   _data/variant/filtered_vcf/*_input.filtered.vcf.gz -> per-sample hard GT
#   results/03_variant_calling/tables/09_functional_annotation.csv -> SnpEff impact
#   _data/_metadata/sample_sheet.csv -> donor, condition, sex
#
# Filters applied here (single derivation point, so every later figure and test
# sees identical carriage):
#   1. impact in {HIGH, MODERATE, LOW}. LOW/synonymous retained: synonymous
#      variants are not inert -- they alter mRNA secondary structure, codon
#      optimality, translation kinetics and RNA-binding/m6A motifs. The HAX1
#      core-ablation locus recovered by 04_asm/07 is itself synonymous/LOW, so
#      dropping LOW would discard it.
#   2. Exclude unplaced contigs (GL*/KI*/JH*/chrUn*).
#   3. No carrier-count filter is applied; see the note below IMPACT_KEEP.
#
# Steps:
#   1. Load donor metadata from the sample sheet.
#   2. Filter the annotation table to the qualifying panel.
#   3. Parse GT per sample; carrier = >=1 alt allele; missing calls are dropped,
#      never scored as reference.
#   4. Emit tidy long form plus the per-donor totals used as the GLM offset.
#
# Outputs:
#   11_burden_matrix.tsv.gz  donor x variant carriers
#   11_donor_totals.csv      per-donor totals; total_variants_detected is the
#                            offset consumed by 17
#   11_build_summary.csv     row counts and exclusions
#
# Usage:
#   python3 03_variant_calling/11_build_burden_matrix.py

from __future__ import annotations

import csv
import gzip
import os
import re
import sys
from collections import defaultdict

import pandas as pd

# --- Dynamic pathing (mirrors the convention used across the pipeline) ---
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.environ.get("PROJECT_ROOT") or os.path.normpath(os.path.join(SCRIPT_DIR, ".."))

VCF_DIR = os.path.join(PROJECT_ROOT, "_data", "variant", "filtered_vcf")
META_FILE = os.path.join(PROJECT_ROOT, "_data", "_metadata", "sample_sheet.csv")
ANNOT_FILE = os.path.join(PROJECT_ROOT, "results", "03_variant_calling", "tables",
                          "09_functional_annotation.csv")
TABLE_DIR = os.path.join(PROJECT_ROOT, "results", "03_variant_calling", "tables")

OUT_MATRIX = os.path.join(TABLE_DIR, "11_burden_matrix.tsv.gz")
OUT_TOTALS = os.path.join(TABLE_DIR, "11_donor_totals.csv")
OUT_SUMMARY = os.path.join(TABLE_DIR, "11_build_summary.csv")

IMPACT_KEEP = ("HIGH", "MODERATE", "LOW")

# No carrier-count filter. An earlier version excluded rows at
# (n_subjects=8, n_mdd=4, n_ctrl=4) as placeholder padding; that was wrong --
# 8 is the minimum of a smooth distribution, not a cap, and (4,4) is the modal
# group split. Rationale and evidence: docs/polygenic_burden_methodology.md s4.
#
# Real constraint, upstream: this table holds only the "High Recurrence (>=8
# Subjects)" tier, so private/low-recurrence variants are absent by construction.

UNPLACED_PREFIXES = ("GL", "KI", "JH", "chrUn")


def log(msg: str) -> None:
    print(f"[11] {msg}", flush=True)


def load_metadata() -> dict[str, dict]:
    """donor -> {condition, sex}. Sample names carry an _input/_m6A suffix."""
    donors: dict[str, dict] = {}
    with open(META_FILE) as fh:
        for row in csv.DictReader(fh):
            sample = re.sub(r"_[12]$", "", row["File1"])
            donor = re.sub(r"_(input|m6A)$", "", sample)
            donors[donor] = {
                "condition": row["SampleType"].strip(),
                "sex": row["Gender"].strip(),
            }
    return donors


def parse_gt(sample_field: str) -> tuple[int, int] | None:
    """Return (ref_alleles, alt_alleles) for a hard genotype, else None.

    Takes the WHOLE sample field (e.g. '0/1:13,7:20:99:208,0,458') and splits
    off the GT subfield first. Handles 0/1, 1/1, 0|1, 1|1, mixed phasing and
    mixed depth. Missing calls ('.', './.') return None so they are never
    scored as reference -- an uncalled site is not evidence of homozygous ref.
    """
    if not sample_field:
        return None
    gt = sample_field.split(":")[0]
    if not gt:
        return None
    alleles = gt.replace("|", "/").split("/")
    if not alleles or any(a in (".", "") for a in alleles):
        return None
    try:
        idx = [int(a) for a in alleles]
    except ValueError:
        return None
    n_alt = sum(1 for a in idx if a > 0)
    return len(idx) - n_alt, n_alt


def load_annotation() -> dict[str, dict]:
    """variant_id -> annotation, keeping only qualifying impact classes."""
    csv.field_size_limit(10**7)
    keep: dict[str, dict] = {}
    stats = {"total": 0, "impact": 0, "unplaced": 0, "no_alt": 0}

    with open(ANNOT_FILE) as fh:
        for row in csv.DictReader(fh):
            stats["total"] += 1
            vid = row["variant_id"]
            chrom = (row["chrom"] or "").strip()

            ns = int(row["n_subjects"] or 0)

            if chrom.startswith(UNPLACED_PREFIXES):
                stats["unplaced"] += 1
                continue
            if row["impact"] not in IMPACT_KEEP:
                continue
            if ns > 0 and (int(row["n_het"] or 0) + int(row["n_hom_alt"] or 0)) == 0:
                stats["no_alt"] += 1
                continue

            stats["impact"] += 1
            keep[vid] = {
                "chrom": chrom,
                "pos": int(row["pos"]),
                "gene": row["gene"],
                "impact": row["impact"],
                "consequence": row["consequence"],
                "n_subjects": ns,
            }
    return keep, stats


def parse_vcf(path: str, annot: dict[str, dict]) -> tuple[list[tuple], int, int, int]:
    """Return (carrier rows, panel sites seen, off-panel records, total records).

    `total records` is the donor's FULL detected variant load and is the value
    used as the GLM offset in 17. It must NOT be the panel count: the panel
    count is the outcome being modelled, so using it as its own offset would
    be circular. Total records scale with the donor's detection depth across
    the whole transcriptome and are independent of any specific gene set.
    """
    rows: list[tuple] = []
    seen = skipped = total = 0
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            total += 1
            f = line.rstrip("\n").split("\t")
            vid = f"{f[0]}:{f[1]}:{f[3]}:{f[4]}"
            a = annot.get(vid)
            if a is None:
                skipped += 1
                continue
            seen += 1
            parsed = parse_gt(f[9])
            if parsed is None:
                continue
            _ref_n, alt_n = parsed
            if alt_n == 0:
                continue  # reference at this site contributes no burden
            gene = a["gene"] or ""
            if not gene:
                continue
            rows.append((gene, a["impact"], vid, a["chrom"], a["pos"], alt_n))
    return rows, seen, skipped, total


def main() -> int:
    log("=" * 78)
    log("BUILD BURDEN MATRIX (donor x variant carriers)")
    log("=" * 78)

    if not os.path.isdir(VCF_DIR):
        log(f"FATAL: VCF directory not found: {VCF_DIR}")
        return 1
    if not os.path.isfile(ANNOT_FILE):
        log(f"FATAL: annotation table not found: {ANNOT_FILE}")
        return 1
    os.makedirs(TABLE_DIR, exist_ok=True)

    donors = load_metadata()
    log(f"donors from sample sheet: {len(donors)}")
    n_mdd = sum(1 for d in donors.values() if d["condition"] == "mdd")
    log(f"  MDD {n_mdd} / Ctrl {len(donors) - n_mdd}")

    annot, astats = load_annotation()
    log(f"annotation rows: {astats['total']:,}")
    log(f"  qualifying impact {IMPACT_KEEP}: {astats['impact']:,}")
    log(f"  excluded unplaced contigs: {astats['unplaced']:,}")
    log(f"  excluded no-alt rows: {astats['no_alt']:,}")

    vcf_files = sorted(f for f in os.listdir(VCF_DIR) if f.endswith("_input.filtered.vcf.gz"))
    log(f"per-sample VCFs: {len(vcf_files)}")

    all_rows: list[tuple] = []
    donor_seen: dict[str, int] = {}
    donor_panel: dict[str, int] = {}
    donor_total: dict[str, int] = {}

    for fn in vcf_files:
        donor = fn.replace("_input.filtered.vcf.gz", "")
        if donor not in donors:
            log(f"  WARNING: {donor} absent from sample sheet; skipped")
            continue
        path = os.path.join(VCF_DIR, fn)
        rows, seen, _skipped, total = parse_vcf(path, annot)
        all_rows.extend((donor,) + r for r in rows)
        donor_seen[donor] = seen
        donor_panel[donor] = len(rows)
        donor_total[donor] = total
        log(f"  {donor:<16} total={total:>6,}  panel_carriers={len(rows):>6,}  "
            f"panel_share={100 * len(rows) / total:4.1f}%")

    if not all_rows:
        log("FATAL: no carrier rows parsed -- check VCF and annotation join")
        return 1

    df = pd.DataFrame(all_rows, columns=["donor", "gene", "impact", "variant_id",
                                         "chrom", "pos", "alt_alleles"])
    df["condition"] = df["donor"].map(lambda d: donors[d]["condition"])
    df["sex"] = df["donor"].map(lambda d: donors[d]["sex"])

    df.to_csv(OUT_MATRIX, sep="\t", index=False, compression="gzip")

    # total_variants_detected is the donor's FULL detected load (all VCF
    # records, on and off panel). It is the GLM offset consumed by 17 and must
    # be independent of the panel outcome -- see parse_vcf.
    totals = pd.DataFrame([
        {
            "donor": d,
            "condition": meta["condition"],
            "sex": meta["sex"],
            "total_variants_detected": donor_total.get(d, 0),
            "panel_variants_carried": donor_panel.get(d, 0),
        }
        for d, meta in sorted(donors.items())
    ])
    totals.to_csv(OUT_TOTALS, index=False)

    summary = pd.DataFrame([
        ("donors", len(donors)),
        ("mdd", n_mdd),
        ("ctrl", len(donors) - n_mdd),
        ("annotation_rows_total", astats["total"]),
        ("annotation_qualifying_impact", astats["impact"]),
        ("annotation_excluded_unplaced", astats["unplaced"]),
        ("annotation_excluded_no_alt", astats["no_alt"]),
        ("vcf_files", len(donor_seen)),
        ("variants_in_panel_union", len(df["variant_id"].unique())),
        ("genes_in_panel_union", len(df["gene"].unique())),
        ("carrier_rows", len(df)),
        ("burden_total_mdd", int((df["condition"] == "mdd").sum())),
        ("burden_total_ctrl", int((df["condition"] == "ctrl").sum())),
        ("mean_total_variants_mdd", round(float(totals.loc[totals.condition == "mdd", "total_variants_detected"].mean()), 1)),
        ("mean_total_variants_ctrl", round(float(totals.loc[totals.condition == "ctrl", "total_variants_detected"].mean()), 1)),
        ("mean_panel_burden_mdd", round(float(totals.loc[totals.condition == "mdd", "panel_variants_carried"].mean()), 1)),
        ("mean_panel_burden_ctrl", round(float(totals.loc[totals.condition == "ctrl", "panel_variants_carried"].mean()), 1)),
    ], columns=["metric", "value"])
    summary.to_csv(OUT_SUMMARY, index=False)

    log("-" * 78)
    for _, r in summary.iterrows():
        v = float(r["value"])
        log(f"  {r['metric']:<32} {int(v):>12,}" if float(v).is_integer() else f"  {r['metric']:<32} {v:>12,.1f}")
    log("-" * 78)
    log(f"matrix -> {OUT_MATRIX}")
    log(f"totals -> {OUT_TOTALS}")
    log("next: 12_fetch_opentargets_panel.py, then 14_gene_burden_test.py")
    return 0


if __name__ == "__main__":
    sys.exit(main())