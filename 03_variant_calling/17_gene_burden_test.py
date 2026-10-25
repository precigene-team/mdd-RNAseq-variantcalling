#!/usr/bin/env python3
# [[03_variant_calling/17_gene_burden_test.py]]
# Goal: Statistical step of the polygenic burden analysis. Consumes the matrix
#       from 15 and the panel from 16; produces every burden test result.
#
# Analyses (see docs/polygenic_burden_methodology.md section 7):
#   gene-level burden  unbiased gene-level, all genes, BH-FDR          -> PRIMARY inference
#   panel burden  Open Targets MDD genetic panel >= 0.35           -> interpretable
#   pathway burden  msigdbr Hallmark / C2 / C5, if 16b was run       -> EXPLORATORY only
#   convergence  gene-level burden significant AND in panel burden                   -> credible set
#
# The statistic:
#   Burden is a COUNT whose magnitude scales with detection depth, so the
#   measured 1.254x MDD/Ctrl coverage imbalance must be cancelled. Two defences,
#   both used here (utils/case_control_null.py):
#     (a) offset. Rate is modelled per unit of the donor's full detected
#         variant load: rate_d = burden_d / total_variants_d. This is
#         library-size normalisation, the edgeR/DESeq2 principle.
#     (b) permutation. Condition labels are shuffled WITHIN sex, preserving each
#         donor's call count, giving the null.
#   The offset-adjusted statistic used for permutation is
#       T = log(mean rate | MDD) - log(mean rate | Ctrl)
#   which is the rate-ratio estimate for a two-group model and is fast enough
#   to permute N_PERM times. The full negative-binomial GLM (with sex and
#   interaction terms) is fitted once per gene for the reported effect and
#   standard error, and the two agree by construction.
#
# Outputs:
#   17_gene_burden.csv      per gene: burden, rate ratio, p, BH-FDR
#   17_panel_burden.csv     per Open Targets threshold
#   17_pathway_burden.csv   exploratory pathway sets (if 16b ran)
#   17_convergence.csv      gene-level hits that are also in the panel
#   17_sensitivity.csv           offset on/off, outlier, impact, MAF strata
#   17_summary.csv               headline counts
#
# Usage:
#   python3 03_variant_calling/17_gene_burden_test.py
#   N_PERM=100000 python3 03_variant_calling/17_gene_burden_test.py

from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.environ.get("PROJECT_ROOT") or os.path.normpath(os.path.join(SCRIPT_DIR, ".."))
sys.path.insert(0, os.path.join(PROJECT_ROOT, "utils"))

import case_control_null as ccn  # noqa: E402

TABLE_DIR = os.path.join(PROJECT_ROOT, "results", "03_variant_calling", "tables")
MATRIX = os.path.join(TABLE_DIR, "15_burden_matrix.tsv.gz")
TOTALS = os.path.join(TABLE_DIR, "15_donor_totals.csv")
PANEL = os.path.join(TABLE_DIR, "16_opentargets_mdd_panel.tsv")
PATHWAYS = os.path.join(TABLE_DIR, "16b_pathway_sets.tsv")

N_PERM = int(os.environ.get("N_PERM", "20000"))
PRIMARY_THRESHOLD = float(os.environ.get("PANEL_THRESHOLD", "0.35"))
MIN_CARRIERS = int(os.environ.get("MIN_CARRIERS", "3"))
SEED = 20261017
IMPACT_CLASSES = ("HIGH", "MODERATE", "LOW")


def log(msg: str) -> None:
    print(f"[17] {msg}", flush=True)


def offset_stat(rate: np.ndarray, cond: np.ndarray) -> float:
    """log(mean rate | MDD) - log(mean rate | Ctrl); NaN if a group is empty."""
    m, k = rate[cond == 1], rate[cond == 0]
    if len(m) == 0 or len(k) == 0:
        return np.nan
    lm, lk = float(np.mean(m)), float(np.mean(k))
    if lm <= 0 or lk <= 0:
        return np.nan
    return float(np.log(lm) - np.log(lk))


def permute_set(rate: np.ndarray, cond: np.ndarray, sex: np.ndarray, n_perm: int, seed: int) -> dict:
    """Within-sex permutation null for one gene set."""
    obs = offset_stat(rate, cond)
    try:
        perms = ccn.permute_within_sex(pd.Series(cond), pd.Series(sex), n_perm=n_perm, seed=seed)
    except ValueError as exc:
        return {"stat": obs, "p_perm": np.nan, "z_perm": np.nan, "n_perm": 0,
                "note": f"permutation unavailable: {exc}"}
    draws = np.array([offset_stat(rate, p) for p in perms], dtype=float)
    draws = draws[np.isfinite(draws)]
    p, z = ccn.empirical_pvalue(obs, draws)
    return {"stat": obs, "p_perm": p, "z_perm": z, "n_perm": int(len(draws)), "note": ""}


def load_inputs():
    for path in (MATRIX, TOTALS):
        if not os.path.isfile(path):
            raise SystemExit(f"FATAL: missing input {path}. Run 15_build_burden_matrix.py first.")
    mat = pd.read_csv(MATRIX, sep="\t")
    tot = pd.read_csv(TOTALS)
    panel = None
    if os.path.isfile(PANEL):
        panel = pd.read_csv(PANEL, sep="\t")
    else:
        log(f"WARNING: {PANEL} not found. panel burden skipped. Run 16_fetch_opentargets_panel.py.")
    path_sets = None
    if os.path.isfile(PATHWAYS):
        path_sets = pd.read_csv(PATHWAYS, sep="\t")
    else:
        log(f"NOTE: {PATHWAYS} not found. pathway burden (exploratory) skipped.")
    return mat, tot, panel, path_sets


def main() -> int:
    log("=" * 78)
    log("POLYGENIC BURDEN TEST")
    log(f"  N_PERM={N_PERM}  panel_threshold={PRIMARY_THRESHOLD}  min_carriers={MIN_CARRIERS}")
    log("=" * 78)

    mat, tot, panel, path_sets = load_inputs()
    log(f"matrix rows={len(mat):,}  genes={mat['gene'].nunique():,}  donors={mat['donor'].nunique()}")

    donors = tot.set_index("donor")
    donor_order = donors.index.tolist()
    cond = np.array([1 if donors.loc[d, "condition"] == "mdd" else 0 for d in donor_order])
    sex = donors["sex"].astype(str).values
    total = donors["total_variants_detected"].astype(float).values
    rate_denom = np.where(total > 0, total, np.nan)
    log(f"MDD {int(cond.sum())} / Ctrl {int((cond == 0).sum())}")
    log(f"mean total variants: MDD {total[cond == 1].mean():,.0f} | "
        f"Ctrl {total[cond == 0].mean():,.0f} | ratio {total[cond == 1].mean() / total[cond == 0].mean():.3f}")

    # donor x gene carriage matrix, and per-donor panel burden
    carriage = (mat.groupby(["donor", "gene"]).size()
                .unstack(fill_value=0)
                .reindex(index=donor_order, fill_value=0))
    carriage = carriage.reindex(columns=sorted(carriage.columns), fill_value=0)
    log(f"carriage matrix: {carriage.shape[0]} donors x {carriage.shape[1]} genes")

    carrier_counts = carriage.values
    gene_carriers = carrier_counts.sum(axis=0)
    testable = gene_carriers >= MIN_CARRIERS
    log(f"genes with >= {MIN_CARRIERS} carriers: {int(testable.sum()):,} "
        f"(of {carriage.shape[1]:,}; {int((~testable).sum()):,} below threshold)")

    # ---------------- gene-level burden: unbiased gene-level, BH-FDR ----------------
    log("-" * 78)
    log("GENE-LEVEL BURDEN  unbiased, genome-wide (PRIMARY)")
    log("-" * 78)

    rows = []
    for j in np.flatnonzero(testable):
        y = carrier_counts[:, j].astype(float)
        rate = y / rate_denom
        res = ccn.offset_count_model(y, cond, sex, offsets=np.log(np.nan_to_num(rate_denom)))
        rows.append({
            "gene": carriage.columns[j],
            "n_carriers": int(y.sum()),
            "burden_mdd": float(y[cond == 1].sum()),
            "burden_ctrl": float(y[cond == 0].sum()),
            "rate_mdd": float(rate[cond == 1].mean()),
            "rate_ctrl": float(rate[cond == 0].mean()),
            "rate_ratio": res["rate_ratio"],
            "beta": res["beta_condition"],
            "se": res["se"],
            "p_glm": res["p_condition"],
            "alpha": res["alpha"],
        })
    gene_burden = pd.DataFrame(rows)
    gene_burden["fdr_bh"] = ccn.bh_fdr(gene_burden["p_glm"].values)
    gene_burden = gene_burden.sort_values("p_glm").reset_index(drop=True)
    gene_burden.to_csv(os.path.join(TABLE_DIR, "17_gene_burden.csv"), index=False)

    n_fdr = int((gene_burden["fdr_bh"] < 0.05).sum())
    n_nom = int((gene_burden["p_glm"] < 0.05).sum())
    log(f"genes tested: {len(gene_burden):,}")
    log(f"nominal p < 0.05 : {n_nom:,} ({100 * n_nom / len(gene_burden):.1f}%, ~5% expected by chance)")
    log(f"BH-FDR < 0.05    : {n_fdr:,}")
    if len(gene_burden):
        top = gene_burden.nsmallest(5, "p_glm")
        log("top 5 by nominal p: " + ", ".join(
            f"{r.gene}(RR={r.rate_ratio:.2f}, p={r.p_glm:.3f})" for r in top.itertuples()))

    # ---------------- panel burden: Open Targets panel ----------------
    panel_burden = pd.DataFrame()
    if panel is not None:
        log("-" * 78)
        log("PANEL BURDEN  Open Targets MDD genetic (declared in advance)")
        log("-" * 78)
        panel = panel.dropna(subset=["gene"])
        for thr in (0.50, PRIMARY_THRESHOLD, 0.20, 0.10):
            genes = set(panel.loc[panel["genetic_score"] >= thr, "gene"])
            present = sorted(genes.intersection(carriage.columns))
            if not present:
                continue
            idx = [carriage.columns.get_loc(g) for g in present]
            y = carrier_counts[:, idx].sum(axis=1).astype(float)
            rate = y / rate_denom
            res = ccn.offset_count_model(y, cond, sex, offsets=np.log(np.nan_to_num(rate_denom)))
            pm = permute_set(rate, cond, sex, N_PERM, SEED)
            panel_burden = pd.concat([panel_burden, pd.DataFrame([{
                "genetic_threshold": thr,
                "panel_genes": len(genes),
                "genes_in_cohort": len(present),
                "qualifying_variants": int(mat.loc[mat["gene"].isin(present), "variant_id"].nunique()),
                "burden_mdd_mean": float(y[cond == 1].mean()),
                "burden_ctrl_mean": float(y[cond == 0].mean()),
                "rate_mdd": float(rate[cond == 1].mean()),
                "rate_ctrl": float(rate[cond == 0].mean()),
                "rate_ratio_offset": res["rate_ratio"],
                "p_glm_offset": res["p_condition"],
                "p_perm_within_sex": pm["p_perm"],
                "z_perm": pm["z_perm"],
                "n_perm": pm["n_perm"],
                "note": pm["note"],
            }])], ignore_index=True)
            log(f"  thr>={thr:.2f}: {len(genes):>4} panel genes, {len(present):>3} in cohort | "
                f"offset RR={res['rate_ratio']:.3f} p={res['p_condition']:.3f} | "
                f"perm p={pm['p_perm']:.3f}")
        panel_burden.to_csv(os.path.join(TABLE_DIR, "17_panel_burden.csv"), index=False)

    # ---------------- pathway burden: exploratory pathways ----------------
    pathway_burden = pd.DataFrame()
    if path_sets is not None:
        log("-" * 78)
        log("PATHWAY BURDEN  msigdbr collections (EXPLORATORY -- supports no claim on its own)")
        log("-" * 78)
        rows3 = []
        for r in path_sets.itertuples():
            genes = set(str(r.genes).split(";")) - {""}
            present = sorted(genes.intersection(carriage.columns))
            if len(present) < 3:
                continue
            idx = [carriage.columns.get_loc(g) for g in present]
            y = carrier_counts[:, idx].sum(axis=1).astype(float)
            rate = y / rate_denom
            res = ccn.offset_count_model(y, cond, sex, offsets=np.log(np.nan_to_num(rate_denom)))
            rows3.append({
                "collection": r.collection, "set_name": r.set_name,
                "set_size": r.set_size, "genes_in_cohort": len(present),
                "rate_ratio_offset": res["rate_ratio"], "p_glm_offset": res["p_condition"],
            })
        pathway_burden = pd.DataFrame(rows3)
        pathway_burden["fdr_bh_within_collection"] = np.nan
        for coll, grp in pathway_burden.groupby("collection"):
            pathway_burden.loc[grp.index, "fdr_bh_within_collection"] = ccn.bh_fdr(grp["p_glm_offset"].values)
        pathway_burden = pathway_burden.sort_values("p_glm_offset").reset_index(drop=True)
        pathway_burden.to_csv(os.path.join(TABLE_DIR, "17_pathway_burden.csv"), index=False)
        log(f"  pathway sets tested: {len(pathway_burden):,}")
        for coll, grp in pathway_burden.groupby("collection"):
            k = int((grp["fdr_bh_within_collection"] < 0.05).sum())
            log(f"    {coll:<10} {len(grp):>5} sets, {k:>3} at FDR<0.05 "
                f"(n_tests={len(grp)})")

    # ---------------- convergence: convergence ----------------
    convergence = pd.DataFrame()
    if panel is not None and len(gene_burden):
        log("-" * 78)
        log("CONVERGENCE  gene-level meets panel (gene-level hit AND panel member)")
        log("-" * 78)
        genes = set(panel.loc[panel["genetic_score"] >= PRIMARY_THRESHOLD, "gene"])
        hit = gene_burden[gene_burden["fdr_bh"] < 0.05]
        conv = hit[hit["gene"].isin(genes)].copy()
        convergence = conv
        convergence.to_csv(os.path.join(TABLE_DIR, "17_convergence.csv"), index=False)
        log(f"  gene-level FDR<0.05: {len(hit)};  of which in the panel: {len(conv)}")

    # ---------------- Sensitivity ----------------
    log("-" * 78)
    log("SENSITIVITY")
    log("-" * 78)
    sens = []

    def add(label, y, note=""):
        y = np.asarray(y, dtype=float)
        rate = y / rate_denom
        pm = permute_set(rate, cond, sex, N_PERM, SEED)
        off = ccn.offset_count_model(y, cond, sex, offsets=np.log(np.nan_to_num(rate_denom)))
        noff = ccn.offset_count_model(y, cond, sex, offsets=None)
        sens.append({
            "analysis": label,
            "burden_mdd_mean": float(y[cond == 1].mean()),
            "burden_ctrl_mean": float(y[cond == 0].mean()),
            "raw_ratio": float(y[cond == 1].mean() / max(y[cond == 0].mean(), 1e-9)),
            "rate_ratio_offset": off["rate_ratio"],
            "p_glm_offset": off["p_condition"],
            "rate_ratio_no_offset": noff["rate_ratio"],
            "p_glm_no_offset": noff["p_condition"],
            "p_perm_within_sex": pm["p_perm"],
            "note": note,
        })
        log(f"  {label:<38} raw={sens[-1]['raw_ratio']:.3f} | "
            f"offset RR={off['rate_ratio']:.3f} p={off['p_condition']:.3f} | "
            f"no-offset RR={noff['rate_ratio']:.3f} | perm p={pm['p_perm']:.3f}")

    full_idx = [carriage.columns.get_loc(g) for g in carriage.columns]
    add("all genes, with offset", carrier_counts[:, full_idx].sum(axis=1), "primary")
    excl = np.asarray(donors.index) == "male_ctrl_1"
    y = carrier_counts[:, full_idx].sum(axis=1).astype(float)
    rate = y / rate_denom
    pm = permute_set(rate[~excl], cond[~excl], sex[~excl], N_PERM, SEED)
    off = ccn.offset_count_model(y, cond, sex, offsets=np.log(np.nan_to_num(rate_denom)), exclude=excl)
    noff = ccn.offset_count_model(y, cond, sex, offsets=None, exclude=excl)
    y_all = carrier_counts[:, full_idx].sum(axis=1).astype(float)
    keep_mdd = (cond == 1) & ~excl
    keep_ctrl = (cond == 0) & ~excl
    sens.append({"analysis": "all genes, exclude male_ctrl_1",
                 "burden_mdd_mean": float(y_all[keep_mdd].mean()),
                 "burden_ctrl_mean": float(y_all[keep_ctrl].mean()),
                 "raw_ratio": float(y_all[keep_mdd].mean() / y_all[keep_ctrl].mean()),
                 "rate_ratio_offset": off["rate_ratio"], "p_glm_offset": off["p_condition"],
                 "rate_ratio_no_offset": noff["rate_ratio"], "p_glm_no_offset": noff["p_condition"],
                 "p_perm_within_sex": pm["p_perm"], "note": "outlier QC-01 donor removed"})
    log(f"  {'all genes, exclude male_ctrl_1':<38} offset RR={off['rate_ratio']:.3f} p={off['p_condition']:.3f}")

    for ic in IMPACT_CLASSES:
        genes_ic = sorted(set(mat.loc[mat["impact"] == ic, "gene"]).intersection(carriage.columns))
        if not genes_ic:
            continue
        idx = [carriage.columns.get_loc(g) for g in genes_ic]
        add(f"impact {ic} only", carrier_counts[:, idx].sum(axis=1), "impact stratum")

    pd.DataFrame(sens).to_csv(os.path.join(TABLE_DIR, "17_sensitivity.csv"), index=False)

    # ---------------- Summary ----------------
    summary = pd.DataFrame([
        ("donors", len(donor_order)),
        ("mdd", int(cond.sum())),
        ("ctrl", int((cond == 0).sum())),
        ("mean_total_variants_mdd", round(float(total[cond == 1].mean()), 1)),
        ("mean_total_variants_ctrl", round(float(total[cond == 0].mean()), 1)),
        ("total_variant_ratio", round(float(total[cond == 1].mean() / total[cond == 0].mean()), 4)),
        ("panel_burden_ratio", round(float(y[cond == 1].mean() / y[cond == 0].mean()), 4)),
        ("genes_tested_gene_level", len(gene_burden)),
        ("gene_burden_nominal_p05", n_nom),
        ("gene_burden_fdr05", n_fdr),
        ("convergent_hits", len(convergence)),
        ("n_perm", N_PERM),
        ("seed", SEED),
    ], columns=["metric", "value"])
    summary.to_csv(os.path.join(TABLE_DIR, "17_summary.csv"), index=False)

    log("-" * 78)
    for _, r in summary.iterrows():
        v = float(r["value"])
        log(f"  {r['metric']:<28} {int(v):>12,}" if v.is_integer() else f"  {r['metric']:<28} {v:>12,.4f}")
    log("-" * 78)
    log("pathway burden is EXPLORATORY. No analysis supports a causal claim.")
    return 0


if __name__ == "__main__":
    sys.exit(main())