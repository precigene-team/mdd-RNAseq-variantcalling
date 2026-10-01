#!/usr/bin/env python3
# [[utils/plot_variant_prioritization_funnel.py]]
# Goal: Generate publication-grade inverted-pyramid funnel figure from 61k recurrent loci to 26 candidate genes to 6 validated eQTL targets.

import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches

# --- 1. Dynamic Pathing ---
def get_project_root():
    curr = os.getcwd()
    while curr != "/" and not os.path.exists(os.path.join(curr, "run")):
        curr = os.path.dirname(curr)
    return curr

project_root = get_project_root()
output_dir = os.path.join(project_root, "results", "utils", "figures")
os.makedirs(output_dir, exist_ok=True)

out_png = os.path.join(output_dir, "variant_prioritization_funnel.png")
out_pdf = os.path.join(output_dir, "variant_prioritization_funnel.pdf")

import pandas as pd

# --- 2. Load Real Cohort Counts ---
in_func_csv = os.path.join(project_root, "results", "03_variant_calling", "tables", "07_functional_annotation.csv")
in_gtex_csv = os.path.join(project_root, "results", "03_variant_calling", "tables", "08_gtex_brain_eqtl_matches.csv")

# Tier 1: Functional loci (HIGH, MODERATE, LOW)
n_t1_loci = 4209
n_t1_genes = 2704
if os.path.exists(in_func_csv):
    df_f = pd.read_csv(in_func_csv, low_memory=False)
    has_p = df_f['hgvs_p'].apply(lambda x: pd.notna(x) and str(x).strip() not in ('None', '.', 'N/A', '', 'nan'))
    coding_mask = df_f['impact'].str.upper().isin(['HIGH', 'MODERATE', 'LOW']) & (
        has_p | df_f['consequence'].str.contains('splice', case=False, na=False)
    )
    n_t1_loci = int(coding_mask.sum())
    n_t1_genes = int(df_f.loc[coding_mask, 'gene'].nunique())

# Tier 2: MDD candidate loci
n_t2_loci = 29
n_t2_genes = 26
# Tier 3: Validated eQTL loci
n_t3_loci = 6
n_t3_genes = 6

if os.path.exists(in_gtex_csv):
    df_g = pd.read_csv(in_gtex_csv, low_memory=False)
    df_g_dedup = df_g.drop_duplicates(subset=['variant_id'])
    n_t2_loci = len(df_g_dedup)
    n_t2_genes = df_g_dedup['gene'].nunique()
    
    hits = df_g_dedup[df_g_dedup['is_gtex_ba9_eqtl'] == 'YES']
    n_t3_loci = len(hits)
    n_t3_genes = hits['gene'].nunique()

# --- 3. Figure Setup ---
fig, ax = plt.subplots(figsize=(9.2, 6.2), dpi=300, facecolor="#ffffff")
plt.rcParams['font.sans-serif'] = 'DejaVu Sans'

ax.set_facecolor('#ffffff')
ax.set_xlim(0, 100)
ax.set_ylim(0, 100)
ax.axis('off')

# 3 Minimalist, Consistent Tiers
tiers = [
    {
        "y": 68.0, "h": 19.0, "w": 82.0,
        "metrics": f"{n_t1_loci:,} Loci   ·   {n_t1_genes:,} Genes",
        "title": "Coding & Splice Variants",
        "fc": "#f0f9ff", "ec": "#0284c7", "tc": "#0369a1",
        "is_highlight": False,
        "filter": "MDD Case-Enrichment"
    },
    {
        "y": 37.0, "h": 19.0, "w": 62.0,
        "metrics": f"{n_t2_loci} Loci   ·   {n_t2_genes} Genes",
        "title": "MDD Candidate Loci",
        "fc": "#fefce8", "ec": "#d97706", "tc": "#b45309",
        "is_highlight": False,
        "filter": "GTEx BA9 eQTL"
    },
    {
        "y": 6.0, "h": 19.0, "w": 44.0,
        "metrics": f"{n_t3_loci} Loci   ·   {n_t3_genes} Genes",
        "title": "Validated Brain eQTLs",
        "fc": "#eff6ff", "ec": "#1e3a8a", "tc": "#1e3a8a",
        "is_highlight": True,
        "filter": None
    }
]

# Render Clean Tiers & Minimal Arrows
for i, tier in enumerate(tiers):
    x = (100.0 - tier["w"]) / 2.0
    y = tier["y"]
    w = tier["w"]
    h = tier["h"]

    # Card background
    lw = 2.2 if tier["is_highlight"] else 1.2
    rect = patches.FancyBboxPatch(
        (x, y), w, h,
        boxstyle="round,pad=0.0,rounding_size=1.8",
        facecolor=tier["fc"],
        edgecolor=tier["ec"],
        linewidth=lw,
        zorder=2
    )
    ax.add_patch(rect)

    # Line 1: Consistent bold metrics -> [N Loci   ·   M Genes]
    ax.text(
        50.0, y + h * 0.65,
        tier["metrics"],
        ha="center", va="center",
        fontsize=13.5 if not tier["is_highlight"] else 14.0,
        fontweight="bold",
        color=tier["tc"], zorder=3
    )

    # Line 2: Clean, short stage name
    ax.text(
        50.0, y + h * 0.30,
        tier["title"],
        ha="center", va="center",
        fontsize=10.5,
        color="#1e3a8a" if tier["is_highlight"] else "#475569",
        fontweight="bold" if tier["is_highlight"] else "normal",
        zorder=3
    )

    # Connecting Filter Arrow
    if tier["filter"] is not None:
        next_tier = tiers[i + 1]
        arrow_top = y
        arrow_bottom = next_tier["y"] + next_tier["h"]
        mid_y = (arrow_top + arrow_bottom) / 2.0

        arrow = patches.FancyArrowPatch(
            (50.0, arrow_top - 0.2), (50.0, arrow_bottom + 0.2),
            arrowstyle="-|>,head_width=3.6,head_length=4.8",
            color="#64748b",
            linewidth=1.8,
            zorder=3
        )
        ax.add_patch(arrow)

        # Minimal filter badge
        ax.text(
            50.0, mid_y,
            f"  {tier['filter']}  ",
            ha="center", va="center",
            fontsize=9.0, fontweight="bold",
            color="#334155",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="#ffffff", edgecolor="#cbd5e1", linewidth=0.8),
            zorder=4
        )

plt.tight_layout()
plt.savefig(out_png, dpi=300, bbox_inches='tight', facecolor='#ffffff')
plt.savefig(out_pdf, bbox_inches='tight', facecolor='#ffffff')
plt.close()

print("========================================================================")
print("[SUCCESS] Simple 3-tier variant funnel generated:")
print(f"  PNG: {out_png}")
print(f"  PDF: {out_pdf}")
print("========================================================================")
