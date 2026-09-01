#!/usr/bin/env python3
# [[utils/plot_alt_transcriptome_snp_recurrence.py]]
# Goal: Minimalist, cohesive visualization of SNP recurrence breakdown and 40,073 high-recurrence loci.

import os
import csv
import numpy as np
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
output_dir = os.path.join(project_root, "results", "03_variant_calling", "alt", "figures")
os.makedirs(output_dir, exist_ok=True)
utils_fig_dir = os.path.join(project_root, "results", "utils", "figures")
os.makedirs(utils_fig_dir, exist_ok=True)

out_png = os.path.join(output_dir, "01_transcriptome_snp_recurrence.png")
out_pdf = os.path.join(output_dir, "01_transcriptome_snp_recurrence.pdf")
out_png_utils = os.path.join(utils_fig_dir, "01_transcriptome_snp_recurrence.png")

# --- 2. Load Summary Data ---
summary_csv = os.path.join(project_root, "results", "03_variant_calling", "alt", "tables", "01_transcriptome_snp_summary.csv")

metrics = {
    "total_subjects": 39,
    "total_het_observations": 1601183,
    "unique_snp_loci": 813969,
    "loci_private_1_subject": 633892,
    "loci_low_recurrence_2_4": 111206,
    "loci_moderate_recurrence_5_7": 28798,
    "loci_high_recurrence_gte_8": 40073
}

if os.path.exists(summary_csv):
    try:
        with open(summary_csv, "r") as f:
            reader = csv.DictReader(f)
            for row in reader:
                k = row["metric"]
                if k in metrics:
                    metrics[k] = float(row["value"]) if "." in row["value"] else int(row["value"])
    except Exception as e:
        print(f"[WARNING] Could not parse {summary_csv}: {e}")

total_loci = int(metrics["unique_snp_loci"])
priv = int(metrics["loci_private_1_subject"])
low = int(metrics["loci_low_recurrence_2_4"])
mod = int(metrics["loci_moderate_recurrence_5_7"])
high = int(metrics["loci_high_recurrence_gte_8"])

# --- 3. Figure Setup ---
fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(14.0, 5.8), dpi=300, facecolor="#ffffff")
plt.rcParams['font.sans-serif'] = 'DejaVu Sans'
plt.subplots_adjust(left=0.06, right=0.96, bottom=0.08, top=0.92, wspace=0.25)

# Cohesive Palette (Sequential Navy/Steel Blue tones)
c_priv = "#cbd5e1"   # Slate 300 (Private - 1 Subject)
c_low  = "#94a3b8"   # Slate 400 (2-4 Subjects)
c_mod  = "#475569"   # Slate 600 (5-7 Subjects)
c_high = "#1e3a8a"   # Deep Navy Blue 900 (>=8 Subjects - Highlight)

# ==============================================================================
# PANEL A: Thicker Donut Starting at 90° (Clockwise 1 -> 2-4 -> 5-7 -> >=8)
# ==============================================================================
ax_a.set_facecolor('#ffffff')

sizes = [priv, low, mod, high]
colors = [c_priv, c_low, c_mod, c_high]

legend_labels = [
    f"1 Subject ({priv:,})",
    f"2–4 Subjects ({low:,})",
    f"5–7 Subjects ({mod:,})",
    f"≥ 8 Subjects ({high:,})"
]

# Start vertical at 90° (12 o'clock) and rotate clockwise
wedges, texts, autotexts = ax_a.pie(
    sizes, 
    colors=colors, 
    startangle=90,
    counterclock=False,
    autopct=lambda pct: f"{pct:.1f}%" if pct >= 3.0 else f"{pct:.2f}%",
    pctdistance=0.74,
    wedgeprops=dict(width=0.58, edgecolor='#ffffff', linewidth=2.0)
)

# Text colors for percentages inside slices
text_colors = ["#0f172a", "#0f172a", "#ffffff", "#ffffff"]
for autotext, col in zip(autotexts, text_colors):
    autotext.set_color(col)
    autotext.set_fontsize(9.5)
    autotext.set_fontweight('bold')

# Center Ring with Clean Total
ax_a.text(0, 0.08, f"{total_loci:,}", fontsize=13.5, fontweight='bold', color='#0f172a', ha='center', va='center')
ax_a.text(0, -0.12, "Unique Loci", fontsize=9.0, color='#64748b', ha='center', va='center')

# Minimal Clean Legend
ax_a.legend(
    wedges, legend_labels,
    loc="center left",
    bbox_to_anchor=(0.90, 0.5),
    frameon=False,
    fontsize=9.0
)

# ==============================================================================
# PANEL B: Minimalist Funnel
# ==============================================================================
ax_b.set_facecolor('#ffffff')
ax_b.set_xlim(0, 10)
ax_b.set_ylim(0, 10)
ax_b.axis('off')

funnel_stages = [
    {
        "title": "Total Het Observations",
        "count": f"{int(metrics['total_het_observations']):,}",
        "w": 8.8, "y": 7.1, "h": 1.6, "fc": "#f8fafc", "ec": "#cbd5e1", "tc": "#475569"
    },
    {
        "title": "Unique SNP Loci",
        "count": f"{total_loci:,}",
        "w": 7.4, "y": 4.1, "h": 1.6, "fc": "#f1f5f9", "ec": "#94a3b8", "tc": "#334155"
    },
    {
        "title": "High-Recurrence Loci (≥ 8 Subjects)",
        "count": f"{high:,}  (4.92%)",
        "w": 6.0, "y": 1.1, "h": 1.6, "fc": "#eff6ff", "ec": "#1e3a8a", "tc": "#1e3a8a"
    }
]

for idx, stg in enumerate(funnel_stages):
    x_left = 5.0 - stg["w"] / 2.0
    is_last = (idx == 2)
    box = patches.FancyBboxPatch((x_left, stg["y"]), stg["w"], stg["h"], boxstyle="round,pad=0.10",
                                 facecolor=stg["fc"], edgecolor=stg["ec"], lw=2.0 if is_last else 1.2)
    ax_b.add_patch(box)
    
    ax_b.text(5.0, stg["y"] + 0.98, stg["title"], 
              fontsize=10.5, fontweight="bold" if is_last else "normal", color=stg["tc"], ha="center")
    ax_b.text(5.0, stg["y"] + 0.40, stg["count"], 
              fontsize=11.5 if is_last else 10.5, fontweight="bold", color="#0f172a" if not is_last else "#1e3a8a", ha="center")
    
    if idx < len(funnel_stages) - 1:
        next_y = funnel_stages[idx+1]["y"] + funnel_stages[idx+1]["h"]
        ax_b.annotate("", xy=(5.0, next_y + 0.08), xytext=(5.0, stg["y"] - 0.08),
                      arrowprops=dict(arrowstyle="->", color="#94a3b8", lw=1.8))

# Save PNG and PDF
plt.savefig(out_png, dpi=300, bbox_inches='tight')
plt.savefig(out_pdf, bbox_inches='tight')
plt.savefig(out_png_utils, dpi=300, bbox_inches='tight')
plt.close()

print(f"[SUCCESS] Donut updated to start at 90° clockwise order (1 -> 2-4 -> 5-7 -> >=8):\n  PNG: {out_png}")
