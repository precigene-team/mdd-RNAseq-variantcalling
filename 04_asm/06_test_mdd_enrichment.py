#!/usr/bin/env python3
# [[04_asm/06_test_mdd_enrichment.py]]
# Goal: Test whether the core-m6A site set is enriched for MDD carriers, calibrated by permuting donor labels.
#
# Steps:
#   1. Reconstruct the donor x site carrier matrix for the core-m6A set from the master table.
#   2. Diagnose the coverage confound: compare per-donor variant call rates between groups.
#   3. Test the pooled MDD carrier excess against a donor-label permutation null.
#   4. Contrast the core set against the non-peak expressed-locus background, which is the interpretable comparison.
#   5. Stratify by sex and report.
#
# Rationale (docs/04_asm_literature.md):
#   - Per-site selection on n_mdd is circular: n_mdd IS the case-control skew.
#     A pooled test asks the question without selecting on the outcome.
#   - m6A-QTL work (Zhang 2020 Nat Genet) calibrates empirically rather than
#     trusting a nominal p-value across many related tests.
#   - Donor-label permutation is the correct null here. It preserves each donor's
#     total call count, so a donor with deeper coverage stays a high-count donor.
#     That is essential here because group call rates are NOT equal (see step 2).
#
# Usage:
#   python3 04_asm/06_test_mdd_enrichment.py
#   N_PERM=50000 python3 04_asm/06_test_mdd_enrichment.py

import os
import re
import collections

import numpy as np
import pandas as pd

N_PERM = int(os.environ.get("N_PERM", "20000"))
SEED = 20240117          # fixed so the reported p is reproducible
DONOR_RE = re.compile(r"^(male|female)_(mdd|ctrl)_(\d+)_input$")


def get_project_root():
    curr = os.getcwd()
    while curr != "/" and not os.path.exists(os.path.join(curr, "run")):
        curr = os.path.dirname(curr)
    return curr


project_root = get_project_root()
TABLE_DIR = os.path.join(project_root, "results", "04_asm", "tables")
CORE = os.path.join(TABLE_DIR, "05_core_m6a_sites.csv")
MASTER = os.path.join(project_root, "results", "03_variant_calling", "tables",
                      "06_transcriptome_snp_master.tsv.gz")

if not os.path.exists(CORE):
    raise SystemExit(f"Assertion Error: missing input: {CORE}\n"
                     "Run 05_plot_m6asnp_funnel.py first.")


# --- 1. Donor x site carrier matrix ---
print("[STEP 1] Reconstructing the donor x site carrier matrix ...")
core = pd.read_csv(CORE, low_memory=False)
core_loci = set(core["locus"])

master = pd.read_csv(MASTER, sep="\t",
                     usecols=["chr", "pos", "ref", "alt", "n_subjects", "n_mdd",
                              "n_ctrl", "carrier_subjects"],
                     low_memory=False)
master["locus"] = (master["chr"].astype(str) + ":" + master["pos"].astype(str)
                   + ":" + master["ref"].astype(str) + ":" + master["alt"].astype(str))
master = master.drop_duplicates("locus")

donors, sex_of, grp_of = [], {}, {}
carrier_sets = {}
donor_calls = collections.Counter()

for locus, subj, n_subj, n_mdd, n_ctrl in zip(
        master["locus"], master["carrier_subjects"].fillna(""),
        master["n_subjects"], master["n_mdd"], master["n_ctrl"]):
    who = []
    for d in subj.split(","):
        mo = DONOR_RE.match(d.strip())
        if not mo:
            continue
        key = f"{mo.group(1)}_{mo.group(2)}_{mo.group(3)}"
        if key not in grp_of:
            donors.append(key)
            sex_of[key] = mo.group(1)
            grp_of[key] = mo.group(2)
        who.append(key)
        donor_calls[key] += 1
    # consistency gate: the listed donors must match the stored counts
    if len(who) != n_subj or sum(grp_of[k] == "mdd" for k in who) != n_mdd \
            or sum(grp_of[k] == "ctrl" for k in who) != n_ctrl:
        raise SystemExit(f"Assertion Error: carrier list disagrees with counts at {locus}")
    if locus in core_loci:
        carrier_sets[locus] = who

donors = sorted(donors)
is_mdd = np.array([grp_of[d] == "mdd" for d in donors])
is_male = np.array([sex_of[d] == "male" for d in donors])
N_DONORS = len(donors)
N_MDD = int(is_mdd.sum())
N_CTRL = int((~is_mdd).sum())
print(f"      donors: {N_DONORS} ({N_MDD} MDD / {N_CTRL} Ctrl)")
print(f"      core sites with a carrier matrix: {len(carrier_sets):,}")

sites = sorted(carrier_sets)
site_of_donor = collections.defaultdict(set)
for i, loc in enumerate(sites):
    for d in carrier_sets[loc]:
        site_of_donor[d].add(i)
# X[d, i] = 1 if donor d carries core site i
X = np.zeros((N_DONORS, len(sites)), dtype=np.float64)
for di, d in enumerate(donors):
    for si in site_of_donor.get(d, ()):
        X[di, si] = 1.0
print(f"      matrix X: {X.shape[0]} donors x {X.shape[1]} sites, "
      f"{int(X.sum()):,} carrier observations")


# --- 2. Coverage confound diagnostic ---
print("[STEP 2] Diagnosing the per-donor call-rate confound ...")
calls = np.array([donor_calls[d] for d in donors], dtype=float)
mdd_calls, ctrl_calls = calls[is_mdd], calls[~is_mdd]


def welch(a, b):
    na, nb = len(a), len(b)
    va, vb = a.var(ddof=1), b.var(ddof=1)
    t = (a.mean() - b.mean()) / np.sqrt(va / na + vb / nb)
    df = (va / na + vb / nb) ** 2 / ((va / na) ** 2 / (na - 1) + (vb / nb) ** 2 / (nb - 1))
    # two-sided p from the t distribution via a normal approximation is adequate
    # here; the ratio, not the p, is the quantity of interest.
    from math import erfc, sqrt
    p = erfc(abs(t) / sqrt(2))
    return t, p


t_all, p_all = welch(mdd_calls, ctrl_calls)
ratio_all = mdd_calls.mean() / ctrl_calls.mean()
print(f"      per-donor calls  MDD {mdd_calls.mean():,.0f}  Ctrl {ctrl_calls.mean():,.0f}"
      f"  ratio {ratio_all:.3f}")
print(f"      Welch t={t_all:+.2f} (normal approx p={p_all:.3g})")
sex_rows = []
for sx, nm in ((is_male, "male"), (~is_male, "female")):
    a, b = calls[is_mdd & sx], calls[~is_mdd & sx]
    tt, pp = welch(a, b)
    sex_rows.append({"stratum": nm, "mdd_donors": len(a), "ctrl_donors": len(b),
                     "mdd_mean_calls": a.mean(), "ctrl_mean_calls": b.mean(),
                     "call_rate_ratio": a.mean() / b.mean(), "welch_t": tt,
                     "approx_p": pp})
    print(f"      {nm:<6} MDD {a.mean():>8,.0f}  Ctrl {b.mean():>8,.0f}"
          f"  ratio {a.mean()/b.mean():.3f}")
call_confound = pd.DataFrame(sex_rows)

# Expected pooled MDD fraction under a call-rate-only null: weight each donor by
# its total call count rather than treating donors as exchangeable.
expected_mdd_frac = calls[is_mdd].sum() / calls.sum()
observed_mdd_frac = X[is_mdd].sum() / X.sum()
print(f"      MDD share of donors            : {N_MDD/N_DONORS:.4f}")
print(f"      MDD share of variant CALLS     : {calls[is_mdd].sum()/calls.sum():.4f}"
      "   <- the confound")
print(f"      core pooled MDD fraction       : {observed_mdd_frac:.4f}")


# --- 3. Donor-label permutation ---
# Null: donor labels are exchangeable within sex. Each donor keeps its own total
# call count, so the coverage imbalance is preserved under the null and cannot
# masquerade as signal.
print(f"[STEP 3] Permuting donor labels within sex ({N_PERM:,} draws) ...")
rng = np.random.default_rng(SEED)
base_frac_early = N_MDD / N_DONORS
# Statistic: excess MDD carriers relative to the donor base rate, signed so that
# positive means MDD-enriched.
obs_stat = (X[is_mdd].sum() * (1 - base_frac_early)
            - X[~is_mdd].sum() * base_frac_early)

perm = np.empty(N_PERM, dtype=np.float64)
base_frac = base_frac_early
for b in range(N_PERM):
    perm_mask = np.zeros(N_DONORS, dtype=bool)
    perm_mask[is_male] = is_mdd[is_male]
    perm_mask[is_male] = rng.permutation(is_mdd[is_male])
    perm_mask[~is_male] = is_mdd[~is_male]
    perm_mask[~is_male] = rng.permutation(is_mdd[~is_male])
    a, c = X[perm_mask].sum(), X[~perm_mask].sum()
    perm[b] = a * (1 - base_frac) - c * base_frac

p_perm = float((np.sum(np.abs(perm) >= abs(obs_stat)) + 1) / (N_PERM + 1))
null_sd = perm.std(ddof=1)
z_perm = obs_stat / null_sd
print(f"      observed statistic  : {obs_stat:+.1f}")
print(f"      permutation null sd : {null_sd:.1f}")
print(f"      z (permutation sd)  : {z_perm:+.2f}")
print(f"      permutation p       : {p_perm:.4f}")


# --- 4. Core vs background contrast ---
print("[STEP 4] Contrasting core against the non-peak expressed-locus background ...")
bg_loci = set(master["locus"]) - core_loci
bg_carrier_sets, bg_donor_calls = {}, collections.Counter()
for locus, subj in zip(master["locus"], master["carrier_subjects"].fillna("")):
    if locus not in bg_loci:
        continue
    who = []
    for d in subj.split(","):
        mo = DONOR_RE.match(d.strip())
        if not mo:
            continue
        key = f"{mo.group(1)}_{mo.group(2)}_{mo.group(3)}"
        bg_donor_calls[key] += 1
        who.append(key)
    bg_carrier_sets[locus] = who

bg_sites = list(bg_carrier_sets)
Xb = np.zeros((N_DONORS, len(bg_sites)), dtype=np.float64)
bg_site_of_donor = collections.defaultdict(set)
for i, loc in enumerate(bg_sites):
    for d in bg_carrier_sets[loc]:
        bg_site_of_donor[d].add(i)
for di, d in enumerate(donors):
    for si in bg_site_of_donor.get(d, ()):
        Xb[di, si] = 1.0
bg_frac = Xb[is_mdd].sum() / Xb.sum()
print(f"      background sites {len(bg_sites):,} | pooled MDD fraction {bg_frac:.4f}")

# The interpretable comparison: core versus background, both on the same scale.
n_c, N_c = int(X.sum()), int(X.size)
n_b, N_b = int(Xb.sum()), int(Xb.size)
se_diff = np.sqrt((n_c / N_c) * (1 - n_c / N_c) / n_c
                  + (n_b / N_b) * (1 - n_b / N_b) / n_b)
z_diff = (observed_mdd_frac - bg_frac) / se_diff
from math import erfc, sqrt
p_diff = erfc(abs(z_diff) / sqrt(2))
print(f"      core   {observed_mdd_frac:.4f}  vs  background {bg_frac:.4f}")
print(f"      difference z = {z_diff:+.2f}  (p={p_diff:.3g})")


# --- 5. Sex strata on the core set ---
strata = []
for nm, mask in (("male", is_male), ("female", ~is_male)):
    a, b = X[is_mdd & mask].sum(), X[~is_mdd & mask].sum()
    nm_m, nm_c = int((is_mdd & mask).sum()), int((~is_mdd & mask).sum())
    strata.append({"stratum": nm, "carriers": int(a + b), "mdd_carriers": int(a),
                   "ctrl_carriers": int(b), "mdd_fraction": a / (a + b),
                   "donor_base_rate": nm_m / (nm_m + nm_c),
                   "call_rate_base_rate": (calls[is_mdd & mask].sum()
                                           / (calls[is_mdd & mask].sum()
                                              + calls[~is_mdd & mask].sum()))})
strata.append({"stratum": "both",
               "carriers": int(X.sum()), "mdd_carriers": int(X[is_mdd].sum()),
               "ctrl_carriers": int(X[~is_mdd].sum()),
               "mdd_fraction": observed_mdd_frac,
               "donor_base_rate": N_MDD / N_DONORS,
               "call_rate_base_rate": expected_mdd_frac})
strata_df = pd.DataFrame(strata)

summary = pd.DataFrame([
    ("core_sites", len(sites)),
    ("donors", N_DONORS),
    ("core_carriers", int(X.sum())),
    ("core_mdd_carriers", int(X[is_mdd].sum())),
    ("core_mdd_fraction", round(observed_mdd_frac, 4)),
    ("donor_base_rate", round(N_MDD / N_DONORS, 4)),
    ("call_rate_base_rate", round(expected_mdd_frac, 4)),
    ("mdd_over_ctrl_call_rate_ratio", round(ratio_all, 3)),
    ("permutation_statistic", round(obs_stat, 1)),
    ("permutation_z", round(z_perm, 2)),
    ("permutation_p_two_sided", round(p_perm, 4)),
    ("n_permutations", N_PERM),
    ("background_sites", len(bg_sites)),
    ("background_mdd_fraction", round(bg_frac, 4)),
    ("core_minus_background_z", round(z_diff, 2)),
], columns=["metric", "value"])

strata_df.to_csv(os.path.join(TABLE_DIR, "06_mdd_enrichment_strata.csv"), index=False)
call_confound.to_csv(os.path.join(TABLE_DIR, "06_call_rate_confound.csv"), index=False)
summary.to_csv(os.path.join(TABLE_DIR, "06_mdd_enrichment_summary.csv"), index=False)

# --- 6. Report ---
print("=" * 92)
print("   MDD ENRICHMENT IN THE CORE-m6A SET (pooled, permutation-calibrated)")
print("=" * 92)
for _, r in summary.iterrows():
    v = r["value"]
    txt = f"{int(v):,}" if float(v).is_integer() and abs(v) >= 1 else f"{float(v):.4g}"
    print(f"   {r['metric']:<32} {txt:>14}")
print("-" * 92)
print("   SEX STRATA (compare mdd_fraction against BOTH nulls)")
print(strata_df.to_string(index=False))
print("=" * 92)
print(f"   Strata   -> {os.path.join(TABLE_DIR, '06_mdd_enrichment_strata.csv')}")
print(f"   Confound -> {os.path.join(TABLE_DIR, '06_call_rate_confound.csv')}")
print(f"   Summary  -> {os.path.join(TABLE_DIR, '06_mdd_enrichment_summary.csv')}")
print("=" * 92)
