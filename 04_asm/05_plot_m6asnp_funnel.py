#!/usr/bin/env python3
# [[04_asm/05_plot_m6asnp_funnel.py]]
# Goal: Render the m6A-SNP discovery funnel and emit the core-m6A site table consumed by steps 06 and 07.
#
# Steps:
#   1. Derive all six funnel levels from source data (master locus table, consensus peak BED, DRACH overlap table).
#   2. Assert each level is a strict subset of the one above it.
#   3. Write the level table, the core-m6A site table, and the recurrence distribution of the pool.
#   4. Render the pyramid.
#
# Notes:
#   - This is the single derivation point for the funnel. Steps 06 and 07 read
#     05_core_m6a_sites.csv from here rather than re-deriving it, so no script can
#     disagree about which sites are in the core set.
#   - LEVEL 6 IS NOT "m6A SITES". Motif presence does not establish methylation:
#     direct-RNA work finds only ~5% of available DRACH motifs are methylated in
#     vivo. These are adenosines in a DRACH context inside an IP-enriched region.
#     See docs/04_asm_literature.md section 2.1.
#   - No recurrence threshold is applied. Filtering on donor count before the
#     mechanistic steps would condition the m6A set on an arbitrary cutoff, so
#     recurrence is reported as a distribution and left for a later decision.
#   - Peak containment is half-open: the BED is 0-based [start,end) and a 1-based
#     SNP position p sits at 0-based p-1, so the test is start < p <= end.
#     Treating the BED end as inclusive admits variants one base past every peak.

import os
import bisect

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import pandas as pd

N_MDD, N_CTRL, N_DONORS = 20, 19, 39

# --- 1. Dynamic Pathing ---
def get_project_root():
    curr = os.getcwd()
    while curr != "/" and not os.path.exists(os.path.join(curr, "run")):
        curr = os.path.dirname(curr)
    return curr

project_root = get_project_root()
TABLE_DIR = os.path.join(project_root, "results", "04_asm", "tables")
FIG_DIR = os.path.join(project_root, "results", "04_asm", "figures")
os.makedirs(TABLE_DIR, exist_ok=True)
os.makedirs(FIG_DIR, exist_ok=True)

MASTER = os.path.join(project_root, "results", "03_variant_calling", "tables",
                      "06_transcriptome_snp_master.tsv.gz")
PEAKS = os.path.join(project_root, "_data", "peaks", "consensus_peaks.bed")
OVERLAPS = os.path.join(TABLE_DIR, "04_drach_variant_overlaps.csv")

out_png = os.path.join(FIG_DIR, "05_m6asnp_funnel.png")
out_pdf = os.path.join(FIG_DIR, "05_m6asnp_funnel.pdf")

for p in (MASTER, PEAKS, OVERLAPS):
    if not os.path.exists(p):
        raise SystemExit(f"Assertion Error: missing input: {p}\n"
                         "Run steps 01-04 first.")

# --- 2. Expressed loci ---
# carrier_subjects holds a comma-joined donor list on all 855k rows and dominates
# read time, so it is not read here. Read only the columns the funnel needs.
NEEDED = ["chr", "pos", "ref", "alt", "n_subjects", "n_het", "n_hom_alt",
          "n_mdd", "n_ctrl", "n_male", "n_female", "n_male_mdd", "n_female_mdd",
          "n_male_ctrl", "n_female_ctrl"]

print("[STEP 1] Reading expressed-locus master table ...")
master = pd.read_csv(MASTER, sep="\t", usecols=NEEDED, low_memory=False)
master["locus"] = (master["chr"].astype(str) + ":" + master["pos"].astype(str)
                   + ":" + master["ref"].astype(str) + ":" + master["alt"].astype(str))
master = master.drop_duplicates("locus").reset_index(drop=True)
n_obs = int(master["n_subjects"].sum())
n_distinct = len(master)
print(f"      donor x variant observations : {n_obs:,}")
print(f"      distinct expressed loci       : {n_distinct:,}")

# --- 3. Within consensus m6A peaks ---
peaks = pd.read_csv(PEAKS, sep="\t", header=None, usecols=[0, 1, 2],
                    names=["chrom", "start", "end"])
peaks["chrom"] = peaks["chrom"].astype(str)
index = {}
for chrom, grp in peaks.groupby("chrom", sort=False):
    spans = sorted(zip(grp["start"].astype(int), grp["end"].astype(int)))
    index[chrom] = (np.array([s for s, _ in spans]),
                    np.array([e for _, e in spans]))
print(f"      consensus peaks              : {len(peaks):,}")

chrom = master["chr"].astype(str).values
pos = master["pos"].astype(int).values
in_peak = np.zeros(len(master), dtype=bool)
for c, idx in master.groupby(master["chr"].astype(str)).indices.items():
    entry = index.get(c)
    if entry is None:
        continue
    starts, ends = entry
    i = np.searchsorted(starts, pos[idx] - 1, side="right") - 1
    ok = i >= 0
    ii = np.where(ok, i, 0)
    ok &= starts[ii] < pos[idx]
    ok &= pos[idx] <= ends[ii]
    in_peak[idx] = ok

master["in_peak"] = in_peak
in_peak_loci = set(master.loc[in_peak, "locus"])
print(f"      expressed loci in peaks      : {len(in_peak_loci):,}")

# --- 4. DRACH overlap, disruption, core ablation ---
print("[STEP 2] Annotating DRACH overlap and disruption ...")
ov = pd.read_csv(OVERLAPS, low_memory=False)
# The overlap table is per donor, so a variant carried by several donors
# contributes several rows. Collapse to one row per (variant, motif instance).
per_overlap = ov.drop_duplicates("overlap_key")
motif_sites = set(per_overlap["snp_id"])
disrupt_sites = set(per_overlap.loc[
    per_overlap["disruption_status"] == "DISRUPTS_DRACH", "snp_id"])
core_motifs = per_overlap.loc[
    (per_overlap["disruption_status"] == "DISRUPTS_DRACH")
    & (per_overlap["is_core_m6a"].astype(str).str.upper() == "TRUE")]
core_loci = set(core_motifs["snp_id"])

set_drach = in_peak_loci & motif_sites
set_disrupt = in_peak_loci & disrupt_sites
set_core = in_peak_loci & core_loci
print(f"      ... overlapping DRACH         : {len(set_drach):,}")
print(f"      ... disrupting DRACH          : {len(set_disrupt):,}")
print(f"      ... DRACH-core adenine        : {len(set_core):,}")

# --- 5. Assemble levels, assert strict nesting ---
levels = [
    (1, "Donor x Variant Observations", n_obs,
     "PASS biallelic non-reference calls, 39 donors"),
    (2, "Distinct Expressed Loci", n_distinct,
     "deduplicated chr:pos:ref:alt"),
    (3, "Within Consensus m6A Peaks", len(in_peak_loci),
     "overlaps a TRESS IP-vs-Input region"),
    (4, "Overlapping a DRACH Motif", len(set_drach),
     "writer recognition sequence"),
    (5, "Disrupts DRACH", len(set_disrupt),
     "abolishes writer competence at >=1 instance"),
    (6, "DRACH-Core Adenine (offset 2)", len(set_core),
     "catalytic A removed, within an IP peak"),
]

for (_, _, prev, _), (_, _, cur, _) in zip(levels, levels[1:]):
    if cur > prev:
        raise SystemExit(
            f"Assertion Error: level {cur} exceeds the level above it ({prev}); "
            "the funnel must be strictly nested.")

lvl_df = pd.DataFrame(levels, columns=["level", "label", "value", "filter"])
lvl_df["dropped"] = lvl_df["value"].shift(1) - lvl_df["value"]
lvl_df["pct_of_prev"] = 100 * lvl_df["value"] / lvl_df["value"].shift(1)
lvl_df["pct_of_full"] = 100 * lvl_df["value"] / lvl_df["value"].iloc[0]
lvl_df.to_csv(os.path.join(TABLE_DIR, "05_funnel_levels.csv"), index=False)

# --- 6. Core site table (handed to steps 06 and 07) ---
core = master[master["locus"].isin(set_core)].copy()
if core["locus"].nunique() != len(set_core):
    raise SystemExit("Assertion Error: core site count disagrees with the master table")

mech = (core_motifs.sort_values(["snp_id", "motif_start0"])
        .drop_duplicates("snp_id")
        [["snp_id", "strand", "motif_offset", "motif_seq", "motif_start0",
          "methyl_pos_1based", "tx_ref", "tx_alt"]]
        .rename(columns={"snp_id": "locus", "strand": "motif_strand",
                         "motif_offset": "core_motif_offset",
                         "motif_seq": "core_motif_seq",
                         "motif_start0": "core_motif_start0",
                         "methyl_pos_1based": "core_methyl_pos",
                         "tx_ref": "core_ref", "tx_alt": "core_alt"}))
core = core.merge(mech, on="locus", how="left")
core["mdd_fraction"] = core["n_mdd"] / (core["n_mdd"] + core["n_ctrl"])
core["case_excess"] = core["mdd_fraction"] - (N_MDD / (N_MDD + N_CTRL))

core_cols = ["locus", "n_mdd", "n_ctrl", "n_subjects", "n_het", "n_hom_alt",
             "mdd_fraction", "case_excess", "n_male", "n_female", "n_male_mdd",
             "n_female_mdd", "n_male_ctrl", "n_female_ctrl", "motif_strand",
             "core_motif_offset", "core_motif_seq", "core_methyl_pos",
             "core_ref", "core_alt", "core_motif_start0"]
core = core.sort_values(["n_mdd", "n_ctrl"], ascending=[False, True])
core[[c for c in core_cols if c in core.columns]].to_csv(
    os.path.join(TABLE_DIR, "05_core_m6a_sites.csv"), index=False)

# --- 7. Deferred: recurrence distribution of the core pool ---
rec_core = core["n_subjects"]
rec_tab = (pd.cut(rec_core, bins=[0, 1, 2, 3, 5, 8, 12, 20, 39],
                  labels=["1", "2", "3", "4-5", "6-8", "9-12", "13-20", "21-39"])
           .value_counts().sort_index().rename_axis("n_carrier_donors")
           .reset_index(name="n_sites"))
rec_tab["pct"] = 100 * rec_tab["n_sites"] / rec_tab["n_sites"].sum()
rec_tab.to_csv(os.path.join(TABLE_DIR, "05_core_pool_recurrence.csv"), index=False)

# --- 8. Symmetric pyramid ---
plt.rcParams['font.sans-serif'] = 'DejaVu Sans'
n_lv = len(levels)
top_y, bot_y = 88.0, 3.0
slot = (top_y - bot_y) / n_lv
tier_h = slot * 0.64
gap = slot - tier_h

BLUE, AMBER, TEAL = "#ebf5fb", "#fef9e7", "#e8f8f5"
EC_BLUE, EC_AMBER, EC_TEAL = "#1f4e79", "#f39c12", "#17a589"
TC_BLUE, TC_AMBER, TC_TEAL = "#1b4f72", "#d35400", "#117864"

# widths taper linearly, independent of magnitude, so the pyramid reads as a
# funnel while the numbers stay honest in the labels
w_top, w_bot = 94.0, 46.0
widths = [w_top + (w_bot - w_top) * (i / (n_lv - 1)) for i in range(n_lv)]

fig, ax = plt.subplots(figsize=(9.6, 9.2), dpi=300, facecolor="#ffffff")
ax.set_facecolor('#ffffff')
ax.set_xlim(0, 100)
ax.set_ylim(0, 100)
ax.axis('off')

for i, (num, label, value, _filt) in enumerate(levels):
    y = top_y - slot * i - gap
    w = widths[i]
    x = (100.0 - w) / 2.0

    if i == n_lv - 1:
        box = patches.FancyBboxPatch(
            (x, y), w, tier_h, boxstyle="round,pad=0.2",
            facecolor="#fdecec", edgecolor="#b91c1c", lw=2.8, zorder=2)
        ax.add_patch(box)
        badge_h = tier_h * 0.34
        badge = patches.FancyBboxPatch(
            (x + 2.0, y + tier_h - badge_h - 1.0), w - 4.0, badge_h,
            boxstyle="round,pad=0.05", facecolor="#b91c1c", edgecolor="none", zorder=3)
        ax.add_patch(badge)
        ax.text(50.0, y + tier_h - badge_h / 2.0 - 1.0, label, fontsize=10.0,
                fontweight="bold", color="#ffffff", ha="center", va="center", zorder=4)
        ax.text(50.0, y + tier_h * 0.26, f"{value:,} Sites", fontsize=13.0,
                fontweight="bold", color="#991b1b", ha="center", va="center", zorder=4)
    else:
        if i < 2:
            fc, ec, tc = BLUE, EC_BLUE, TC_BLUE
        elif i < 4:
            fc, ec, tc = AMBER, EC_AMBER, TC_AMBER
        else:
            fc, ec, tc = TEAL, EC_TEAL, TC_TEAL
        box = patches.FancyBboxPatch(
            (x, y), w, tier_h, boxstyle="round,pad=0.2",
            facecolor=fc, edgecolor=ec, lw=1.5, zorder=2)
        ax.add_patch(box)
        ax.text(50.0, y + tier_h * 0.72, label, fontsize=10.0, fontweight="bold",
                color=tc, ha="center", va="center", zorder=3)
        ax.text(50.0, y + tier_h * 0.28, f"{value:,}", fontsize=12.5,
                fontweight="bold", color="#2c3e50", ha="center", va="center", zorder=3)

    # Arrow spanning level i -> i+1, labelled by the criterion that PRODUCES
    # level i+1 (i.e. the next level's own definition) plus the drop count.
    if i < n_lv - 1:
        nxt_y = top_y - slot * (i + 1) - gap + tier_h
        y_mid = (y + nxt_y) / 2.0
        dropped = value - levels[i + 1][2]
        ax.annotate("", xy=(50.0, nxt_y + 0.3), xytext=(50.0, y - 0.3),
                    arrowprops=dict(arrowstyle="->", color="#34495e", lw=2.0), zorder=4)
        ax.text(48.5, y_mid, levels[i + 1][3], fontsize=8.0, fontweight="bold",
                color="#2c3e50", ha="right", va="center", zorder=5)
        ax.text(51.5, y_mid, f"- {dropped:,}", fontsize=8.0, fontweight="bold",
                color="#c0392b", ha="left", va="center", zorder=5)

ax.text(50.0, 99.5, "m6A-SNP Discovery Funnel", fontsize=15, fontweight="bold",
        color="#111827", ha="center", va="top")
ax.text(50.0, 96.2,
        "39 postmortem vmPFC donors (20 MDD / 19 Ctrl)   |   every level derived from source data",
        fontsize=8.5, color="#4b5563", ha="center", va="top")
ax.text(50.0, 93.4,
        "Level 6 = adenosines in a DRACH context within IP-enriched regions; motif presence does not establish methylation",
        fontsize=7.4, color="#9ca3af", ha="center", va="top", style="italic")

plt.savefig(out_png, dpi=300, bbox_inches='tight', facecolor='#ffffff')
plt.savefig(out_pdf, bbox_inches='tight', facecolor='#ffffff')
plt.close()

print("=" * 78)
print("   FUNNEL")
for num, label, value, _filt in levels:
    print(f"   {num}. {label:<34} {value:>12,}")
print("=" * 78)
print("   [DEFERRED] Core-pool recurrence distribution (not a funnel level)")
for _, r in rec_tab.iterrows():
    print(f"      {str(r['n_carrier_donors']):>6} donors : {int(r['n_sites']):>6,}"
          f"  ({r['pct']:.1f}%)")
print("=" * 78)
print(f"   Levels     -> {os.path.join(TABLE_DIR, '05_funnel_levels.csv')}")
print(f"   Core sites -> {os.path.join(TABLE_DIR, '05_core_m6a_sites.csv')}")
print(f"   Recurrence -> {os.path.join(TABLE_DIR, '05_core_pool_recurrence.csv')}")
print(f"   Figure     -> {out_png}")
print("=" * 78)
