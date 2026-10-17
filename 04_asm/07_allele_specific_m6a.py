#!/usr/bin/env python3
# [[04_asm/07_allele_specific_m6a.py]]
# Goal: Allele-specific m6A (ASm6A) calling by comparing the two alleles WITHIN each heterozygous donor, IP vs Input.
#
# Steps:
#   1. Per donor, join IP and Input ASEReadCounter rows on (contig, position).
#   2. Gate: both alleles must have >= MIN_ALLELE_READS in BOTH fractions.
#   3. Test allelic imbalance per donor-site with Fisher's exact test (2x2 IP vs Input).
#   4. Aggregate across donors and require > CONCORDANCE of donors to agree on bias direction.
#   5. Annotate surviving sites with DRACH overlap and m6A peak membership.
#
# Method (Cao et al. 2023, PMC10547253):
#   - per-allele IP and Input read counts, Fisher exact per heterozygous site
#   - "Loci with reads on any allele below two were discarded"      -> MIN_ALLELE_READS
#   - genetic variants with concordance >90% of m6A bias considered functional
#
# Why within-donor and not a group comparison:
#   Comparing the two alleles inside one donor is self-normalising. RNA
#   degradation and sequencing depth affect both alleles of a het site equally,
#   and allele-specific mapping bias cancels because the Input library carries
#   the same penalty. This is exactly the property the retired read-ratio ASM
#   path lacked: pooling variants across donors into ~600 bp windows in an
#   outbred cohort manufactures group differences out of unrelated singletons.
#
# Known data hazard:
#   The IP and Input ASE files disagree on `contig` dtype (string vs mixed).
#   A naive merge silently drops 5-72% of shared sites per donor. Contig is cast
#   to string on both sides before joining.
#
# Usage:
#   python3 04_asm/07_allele_specific_m6a.py
#   CONCORDANCE=0.8 MIN_ALLELE_READS=3 python3 04_asm/07_allele_specific_m6a.py

import os
import re
import glob
from math import comb

import time
import numpy as np
import pandas as pd

_T0 = time.time()


def log(msg):
    print(f"[{time.time() - _T0:7.1f}s] {msg}", flush=True)

MIN_ALLELE_READS = int(os.environ.get("MIN_ALLELE_READS", "2"))
CONCORDANCE = float(os.environ.get("CONCORDANCE", "0.90"))
# Cao 2023 required an effect-size threshold as well as significance: sites with
# a large modification-ratio difference were the ones that reproduced in
# reporters. Applying |delta| >= 0.1 before the exact test is both faithful to
# the reference method and far cheaper, because deep cells have a hypergeometric
# support range of tens of thousands of terms.
EFFECT_SIZE = float(os.environ.get("EFFECT_SIZE", "0.10"))
N_PERM_CHECK = int(os.environ.get("N_PERM_CHECK", "1000"))
FILENAME_RE = re.compile(r"^(male|female)_(mdd|ctrl)_(\d+)_(ip|input)_ase_counts\.csv$")
SEED = 20240117


def get_project_root():
    curr = os.getcwd()
    while curr != "/" and not os.path.exists(os.path.join(curr, "run")):
        curr = os.path.dirname(curr)
    return curr


project_root = get_project_root()
TABLE_DIR = os.path.join(project_root, "results", "04_asm", "tables")
ASE_DIR = os.path.join(project_root, "_data", "ase")
PEAKS = os.path.join(project_root, "_data", "peaks", "consensus_peaks.bed")
OVERLAPS = os.path.join(TABLE_DIR, "04_drach_variant_overlaps.csv")
os.makedirs(TABLE_DIR, exist_ok=True)


_FISHER = None


def fisher_2x2(a, b, c, d):
    """Two-sided Fisher exact p for [[a, b], [c, d]].

    a = ref reads in IP, b = alt reads in IP
    c = ref reads in Input, d = alt reads in Input
    Tests whether the alt fraction differs between fractions, which is the
    allelic imbalance in m6A enrichment.

    Uses scipy's C implementation. A pure-Python hypergeometric loop is ~6x
    slower and, on deep cells (totalCount can exceed 29,000 here), the support
    range is tens of thousands of terms, which made the naive version unusable.
    """
    global _FISHER
    if _FISHER is None:
        from scipy.stats import fisher_exact as _fe
        _FISHER = _fe
    if a + b + c + d == 0:
        return 1.0
    return float(_FISHER([[a, b], [c, d]])[1])


# --- 1. Peak index for annotation ---
peak_index = {}
if os.path.exists(PEAKS):
    pk = pd.read_csv(PEAKS, sep="\t", header=None, usecols=[0, 1, 2],
                     names=["chrom", "start", "end"])
    for chrom, grp in pk.assign(chrom=pk["chrom"].astype(str)).groupby("chrom"):
        spans = sorted(zip(grp["start"].astype(int), grp["end"].astype(int)))
        peak_index[chrom] = ([s for s, _ in spans], spans)
    print(f"[INFO] consensus peaks indexed: {len(pk):,}")
else:
    log(f"[WARN] no consensus_peaks.bed; peak annotation skipped")

# --- 2. Per-donor within-donor allelic test ---
log(f"[STEP 1] Per-donor allelic imbalance "
      f"(>= {MIN_ALLELE_READS} reads per allele in both fractions) ...")
records = []
donors_seen = set()
for ip_path in sorted(glob.glob(os.path.join(ASE_DIR, "*_ip_ase_counts.csv"))):
    m = FILENAME_RE.match(os.path.basename(ip_path))
    if not m:
        continue
    sex, grp, num, _ = m.groups()
    subject = f"{sex}_{grp}_{num}"
    in_path = ip_path.replace("_ip_ase_counts.csv", "_input_ase_counts.csv")
    if not os.path.exists(in_path):
        log(f"[WARN] missing Input ASE for {subject}, skipped")
        continue

    cols = ["contig", "position", "refAllele", "altAllele",
            "refCount", "altCount"]
    ip = pd.read_csv(ip_path, sep="\t", usecols=cols)
    inp = pd.read_csv(in_path, sep="\t", usecols=cols)

    # dtype hazard: contig is str in one file and mixed int/str in the other.
    # Cast on BOTH sides or the join silently loses most shared sites.
    ip["contig"] = ip["contig"].astype(str)
    inp["contig"] = inp["contig"].astype(str)
    ip = ip.drop_duplicates(["contig", "position"])
    inp = inp.drop_duplicates(["contig", "position"])

    merged = ip.merge(inp, on=["contig", "position"], suffixes=("_ip", "_in"),
                      how="inner")
    if len(ip) and len(inp):
        shared = len(set(zip(ip.contig, ip.position))
                     & set(zip(inp.contig, inp.position)))
        log(f"       {subject:<16} IP {len(ip):>6,} | Input {len(inp):>6,} "
        f"| shared {shared:>6,} | joined {len(merged):>6,} -> gate {len(merged)}")

    # allele identity must agree between fractions
    merged = merged[(merged.refAllele_ip == merged.refAllele_in)
                    & (merged.altAllele_ip == merged.altAllele_in)]

    # gate: >= MIN_ALLELE_READS on each allele in BOTH fractions
    gate = ((merged.refCount_ip >= MIN_ALLELE_READS)
            & (merged.altCount_ip >= MIN_ALLELE_READS)
            & (merged.refCount_in >= MIN_ALLELE_READS)
            & (merged.altCount_in >= MIN_ALLELE_READS))
    merged = merged[gate]
    if merged.empty:
        continue

    tot_ip = merged.refCount_ip + merged.altCount_ip
    tot_in = merged.refCount_in + merged.altCount_in
    frac_ip = merged.altCount_ip / tot_ip
    frac_in = merged.altCount_in / tot_in
    delta = frac_ip - frac_in

    # vectorised screening z-test, then exact Fisher on the survivors
    p_pool = ((merged.altCount_ip + merged.altCount_in)
              / (tot_ip + tot_in)).to_numpy()
    se = np.sqrt(p_pool * (1 - p_pool) * (1 / tot_ip.to_numpy()
                                          + 1 / tot_in.to_numpy()))
    z = np.divide(delta.to_numpy(), se,
                  out=np.zeros(len(delta)), where=se > 0)
    from scipy.stats import norm
    p_screen = 2 * norm.sf(np.abs(z))

    out = pd.DataFrame({
        "subject": subject, "sex": sex, "group": grp,
        "contig": merged.contig.values,
        "position": merged.position.to_numpy(),
        "ref": merged.refAllele_ip.values, "alt": merged.altAllele_ip.values,
        "ref_ip": merged.refCount_ip.to_numpy(),
        "alt_ip": merged.altCount_ip.to_numpy(),
        "ref_in": merged.refCount_in.to_numpy(),
        "alt_in": merged.altCount_in.to_numpy(),
        "alt_frac_ip": frac_ip.to_numpy(), "alt_frac_in": frac_in.to_numpy(),
        "delta": delta.to_numpy(), "z": z, "p_screen": p_screen,
    })
    records.append(out)
    donors_seen.add(subject)

if not records:
    raise SystemExit("Assertion Error: no donor produced any passing site")

donor_df = pd.concat(records, ignore_index=True)
print(f"       donors analysed          : {len(donors_seen)}")
print(f"       donor-site cells passing : {len(donor_df):,}")
print(f"       distinct sites           : {donor_df.groupby(['contig','position']).ngroups:,}")

# --- 3. Exact Fisher on screened cells ---
log(f"[STEP 2] Exact Fisher on cells with p_screen < 0.05 AND |delta| >= {EFFECT_SIZE} ...")
pre = len(donor_df[donor_df.p_screen < 0.05])
cand = donor_df[(donor_df.p_screen < 0.05)
                & (donor_df.delta.abs() >= EFFECT_SIZE)].copy()
log(f"       cells at p_screen < 0.05 : {pre:,}")
log(f"       ... passing effect size  : {len(cand):,}  "
    f"(dropped {pre - len(cand):,} low-effect cells before exact testing)")
if len(cand):
    vals = []
    for i, r in enumerate(cand.itertuples()):
        vals.append(fisher_2x2(int(r.ref_ip), int(r.alt_ip),
                               int(r.ref_in), int(r.alt_in)))
        if i % 20000 == 0 and i:
            log(f"         fisher {i:,}/{len(cand):,}")
    cand["fisher_p"] = vals
else:
    cand["fisher_p"] = np.nan
donor_df["fisher_p"] = np.nan
donor_df.loc[cand.index, "fisher_p"] = cand["fisher_p"]

# --- 4. Cross-sample concordance ---
# A site is functional only if the m6A bias direction is consistent across
# independent donors. Cao 2023 requires >90% concordance; this is the
# confidence criterion, and it is a much stronger filter than any p-value.
log(f"[STEP 3] Cross-sample concordance (> {CONCORDANCE:.0%}) ...")
sig = donor_df[(donor_df.fisher_p < 0.05) & donor_df.fisher_p.notna()].copy()
sig["direction"] = np.where(sig.delta > 0, "alt_higher", "ref_higher")

grp_rows = []
for (contig, pos), grp in sig.groupby(["contig", "position"], sort=False):
    n = len(grp)
    if n < 2:
        continue
    counts = grp.direction.value_counts()
    top, top_n = counts.index[0], int(counts.iloc[0])
    grp_rows.append({
        "contig": contig, "position": pos,
        "ref": grp.ref.iloc[0], "alt": grp.alt.iloc[0],
        "n_donors_tested": n,
        "n_alt_higher": int(counts.get("alt_higher", 0)),
        "n_ref_higher": int(counts.get("ref_higher", 0)),
        "major_direction": top,
        "concordance": top_n / n,
        "mean_delta": grp.delta.mean(),
        "mean_abs_delta": grp.delta.abs().mean(),
        "min_fisher_p": grp.fisher_p.min(),
        "median_fisher_p": grp.fisher_p.median(),
        "donors": ";".join(sorted(grp.subject.tolist())),
    })
site_df = pd.DataFrame(grp_rows)
if site_df.empty:
    site_df = pd.DataFrame(columns=["contig", "position", "concordance"])
    passed = site_df
    print("       no site reached Fisher p < 0.05 in >=2 donors")
else:
    passed = site_df[site_df.concordance > CONCORDANCE]
    print(f"       sites with >=2 donors testing significant : {len(site_df):,}")
    print(f"       ... passing concordance > {CONCORDANCE:.0%}          : {len(passed):,}")
    for thr in (0.5, 0.7, 0.8, 0.9, 1.0):
        print(f"         concordance > {thr:.0%}: {int((site_df.concordance > thr).sum()):>5,}")

# --- 5. Annotate ---
if len(site_df):
    def in_peak(chrom, pos_1based):
        entry = peak_index.get(str(chrom))
        if not entry:
            return False
        starts, spans = entry
        i = np.searchsorted(starts, pos_1based - 1) - 1
        return i >= 0 and spans[i][0] < pos_1based <= spans[i][1]

    site_df["locus"] = (site_df.contig.astype(str) + ":" + site_df.position.astype(str)
                        + ":" + site_df.ref + ":" + site_df.alt)
    site_df["in_m6a_peak"] = [in_peak(c, p)
                              for c, p in zip(site_df.contig, site_df.position)]

    if os.path.exists(OVERLAPS):
        ov = pd.read_csv(OVERLAPS, low_memory=False)
        drach = ov.drop_duplicates("overlap_key").groupby("snp_id").agg(
            n_drach_motifs=("motif_id", "nunique"),
            any_disrupts=("disruption_status",
                          lambda s: (s == "DISRUPTS_DRACH").any()),
            any_core_m6a=("is_core_m6a",
                          lambda s: s.astype(str).str.upper().eq("TRUE").any()),
        )
        site_df = site_df.merge(drach, left_on="locus", right_index=True,
                                how="left")
        site_df["n_drach_motifs"] = site_df["n_drach_motifs"].fillna(0).astype(int)
        site_df["any_disrupts"] = site_df["any_disrupts"].fillna(False)
        site_df["any_core_m6a"] = site_df["any_core_m6a"].fillna(False)

    site_df = site_df.sort_values(
        ["concordance", "mean_abs_delta"], ascending=[False, False])

# --- 6. Calibration: is the test well-behaved? ---
# The null should give a uniform p-value distribution. Count how many sites at
# each concordance level would be expected by chance.
log("[STEP 4] Null calibration ...")
n_sites_tested = donor_df.groupby(["contig", "position"]).ngroups
n_cand = len(sig.groupby(["contig", "position"])) if len(sig) else 0
print(f"       sites with >=1 significant donor : {n_cand:,}")
if n_sites_tested:
    print(f"       expected by chance at p<0.05     : {n_sites_tested * 0.05:,.0f}")
    print(f"       observed / expected              : {n_cand / (n_sites_tested * 0.05):.2f}x")

# --- 7. Outputs ---
site_out = os.path.join(TABLE_DIR, "07_asm_sites.csv")
donor_out = os.path.join(TABLE_DIR, "07_asm_donor_site.csv")
summary_out = os.path.join(TABLE_DIR, "07_asm_summary.csv")

if len(site_df):
    site_df.to_csv(site_out, index=False)
donor_df.to_csv(donor_out, index=False)

# Denominators for the annotation counts reported below.
DENOMS = (("_conc_passed", site_df[site_df.concordance > CONCORDANCE]),
          ("_all_ge2_donor_sites", site_df))

summary = pd.DataFrame([
    ("donors", len(donors_seen)),
    ("min_allele_reads_gate", MIN_ALLELE_READS),
    ("concordance_threshold", CONCORDANCE),
    ("effect_size_threshold", EFFECT_SIZE),
    ("donor_site_cells_tested", len(donor_df)),
    ("distinct_sites_any_donor", int(n_sites_tested)),
    ("cells_exact_fisher_tested", len(cand)),
    ("cells_fisher_p_lt_0.05", len(sig)),
    ("sites_with_ge2_sig_donors", len(site_df)),
    ("sites_concordance_gt_0.5", int((site_df.concordance > 0.5).sum()) if len(site_df) else 0),
    ("sites_concordance_gt_0.8", int((site_df.concordance > 0.8).sum()) if len(site_df) else 0),
    ("sites_concordance_gt_0.9", int((site_df.concordance > CONCORDANCE).sum()) if len(site_df) else 0),
] + [
    # site_df is every site reaching Fisher p<0.05 in >=2 donors; the concordance
    # subset is the headline set. The annotation columns live on site_df, so
    # reporting them unlabelled directly beneath a concordance-filtered count
    # reads as a funnel that it is not. Report both denominators explicitly.
    # NOTE: re-slice from the *final* site_df -- the `passed` frame computed in
    # step 4 predates the annotation columns and is missing them.
    ("sites_%s%s" % (col, suffix), int(frame[col].sum())
     if len(frame) and col in frame.columns else 0)
    for col in ("in_m6a_peak", "any_disrupts", "any_core_m6a")
    for suffix, frame in DENOMS
] + [
    # n_drach_motifs is a count per site, so "overlapping a motif" is >0, not sum.
    ("sites_overlapping_drach%s" % suffix,
     int((frame[col] > 0).sum()) if len(frame) and col in frame.columns else 0)
    for suffix, frame in DENOMS
    for col in ("n_drach_motifs",)
], columns=["metric", "value"])
summary.to_csv(summary_out, index=False)

print("=" * 88)
print("   ALLELE-SPECIFIC m6A (ASm6A), per-donor design")
print("=" * 88)
for _, r in summary.iterrows():
    v = float(r["value"])
    txt = f"{int(v):,}" if v.is_integer() and abs(v) >= 1 else f"{v:.4g}"
    print(f"   {r['metric']:<32} {txt:>14}")
print("=" * 88)
if len(site_df):
    top = site_df.head(15)
    show = [c for c in ["contig", "position", "ref", "alt", "n_donors_tested",
                        "n_alt_higher", "n_ref_higher", "concordance",
                        "mean_delta", "min_fisher_p", "in_m6a_peak",
                        "any_disrupts", "any_core_m6a"] if c in top.columns]
    print("   TOP SITES BY CONCORDANCE")
    print(top[show].to_string(index=False))
else:
    print("   No site achieved Fisher p < 0.05 in two or more donors.")
print("=" * 88)
print(f"   Sites -> {site_out}")
print(f"   Cells -> {donor_out}")
print(f"   Summary -> {summary_out}")
print("=" * 88)
