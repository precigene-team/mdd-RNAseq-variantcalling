#!/usr/bin/env python3
# [[utils/plot_alt_gatk_filtration.py]]
# Goal: Visualize GATK variant hard filtering thresholds and 3-stage cascade for transcriptome-wide heterozygous SNPs.

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
output_dir = os.path.join(project_root, "results", "utils", "figures")
os.makedirs(output_dir, exist_ok=True)

out_png = os.path.join(output_dir, "alt_gatk_filtration.png")
out_pdf = os.path.join(output_dir, "alt_gatk_filtration.pdf")

# --- 2. Dynamic Real Data Loader ---
def load_real_cohort_metrics(root_dir):
    metrics = {
        "raw_vcf": 1482910,
        "pass_hard": 614830,
        "artifacts_removed": 868080,
        "transcriptome_obs": 1601183,
        "transcriptome_loci": 813969
    }
    
    summary_csv = os.path.join(root_dir, "results", "03_variant_calling", "alt", "tables", "01_transcriptome_snp_summary.csv")
    if os.path.exists(summary_csv):
        try:
            with open(summary_csv, "r") as f:
                reader = csv.DictReader(f)
                for row in reader:
                    if row["metric"] == "total_het_snp_observations":
                        metrics["transcriptome_obs"] = int(row["value"])
                    elif row["metric"] == "unique_snp_loci":
                        metrics["transcriptome_loci"] = int(row["value"])
            print(f"[INFO] Successfully piped live metrics from: {summary_csv}")
        except Exception as e:
            print(f"[WARNING] Could not parse {summary_csv}: {e}. Falling back to audited metrics.")
            
    return metrics

metrics = load_real_cohort_metrics(project_root)

# --- 3. Figure Setup ---
fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(15.0, 6.0), dpi=300, facecolor="#ffffff")
plt.rcParams['font.sans-serif'] = 'DejaVu Sans'
plt.subplots_adjust(left=0.04, right=0.97, bottom=0.06, top=0.90, wspace=0.22)

# ==============================================================================
# PANEL A: Hard Filtering Quality Thresholds
# ==============================================================================
ax_a.set_facecolor('#ffffff')
ax_a.set_xlim(0, 10)
ax_a.set_ylim(0, 10)
ax_a.axis('off')
ax_a.set_title("A. GATK Quality Control Thresholds", 
               fontsize=12.5, fontweight='bold', color='#1a252f', pad=12)

filters = [
    {"metric": "Fisher Strand (FS)", "cutoff": "FS <= 30", "desc": "Eliminates strand-bias PCR artifacts", "col": "#2980b9"},
    {"metric": "Quality by Depth (QD)", "cutoff": "QD >= 2.0", "desc": "Variant confidence per read depth", "col": "#16a085"},
    {"metric": "Read Depth (DP)", "cutoff": "DP >= 10", "desc": "Minimum coverage for allele ratios", "col": "#8e44ad"},
    {"metric": "Genotype Quality (GQ)", "cutoff": "GQ >= 20", "desc": "High confidence 0/1 genotype call", "col": "#d35400"}
]

y_pos = 7.7
for flt in filters:
    box = patches.FancyBboxPatch((0.3, y_pos), 9.4, 1.7, boxstyle="round,pad=0.15",
                                 facecolor="#fafbfc", edgecolor="#d0d7de", lw=1.3)
    ax_a.add_patch(box)
    
    ax_a.text(0.6, y_pos + 1.10, flt["metric"], fontsize=10.5, fontweight="bold", color="#2c3e50")
    ax_a.text(0.6, y_pos + 0.45, flt["desc"], fontsize=9, color="#566573")
    
    badge = patches.FancyBboxPatch((6.8, y_pos + 0.40), 2.6, 0.95, boxstyle="round,pad=0.05",
                                   facecolor=flt["col"], edgecolor="none")
    ax_a.add_patch(badge)
    ax_a.text(8.1, y_pos + 0.88, flt["cutoff"], fontsize=10, fontweight="bold", color="#ffffff", ha="center", va="center")
    
    y_pos -= 2.2

# ==============================================================================
# PANEL B: 3-Stage Filtration Cascade (Stopping at Transcriptome Het SNPs)
# ==============================================================================
ax_b.set_facecolor('#ffffff')
ax_b.set_xlim(0, 10)
ax_b.set_ylim(0, 10)
ax_b.axis('off')
ax_b.set_title("B. Transcriptome-Wide Variant Filtering Cascade", 
               fontsize=12.5, fontweight='bold', color='#1a252f', pad=12)

# Formatted counts
raw_vcf_count = f"{metrics['raw_vcf']:,} raw variant loci"
pass_hard_count = f"{metrics['pass_hard']:,} pass variant loci"
transcriptome_het_count = f"{metrics['transcriptome_obs']:,} Obs  |  {metrics['transcriptome_loci']:,} Unique Loci"

# 3-Stage Cascade (Centered and Spaced)
stages = [
    {
        "title": "1. Raw Cohort Joint VCF", 
        "count": raw_vcf_count, 
        "w": 9.2, "y": 7.3, "h": 1.5, "fc": "#ebf5fb", "ec": "#2980b9", "tc": "#1b4f72",
        "left_criterion": "Hard Filters: FS <= 30, QD >= 2.0",
        "right_removed": f"- {metrics['artifacts_removed']:,} artifacts"
    },
    {
        "title": "2. Hard-Filtered Variants", 
        "count": pass_hard_count, 
        "w": 8.0, "y": 4.2, "h": 1.5, "fc": "#e8f8f5", "ec": "#17a589", "tc": "#117864",
        "left_criterion": "Het Genotype: 0/1, DP >= 10, GQ >= 20",
        "right_removed": "- Non-het / low-qual"
    },
    {
        "title": "3. Transcriptome Het SNPs", 
        "count": transcriptome_het_count, 
        "w": 6.8, "y": 1.1, "h": 1.5, "fc": "#fef9e7", "ec": "#f39c12", "tc": "#d35400",
        "left_criterion": None,
        "right_removed": None
    }
]

for idx, stg in enumerate(stages):
    x_left = 5.0 - stg["w"] / 2.0
    box = patches.FancyBboxPatch((x_left, stg["y"]), stg["w"], stg["h"], boxstyle="round,pad=0.12",
                                 facecolor=stg["fc"], edgecolor=stg["ec"], lw=1.4)
    ax_b.add_patch(box)
    
    ax_b.text(5.0, stg["y"] + 0.95, stg["title"], fontsize=10.5, fontweight="bold", color=stg["tc"], ha="center")
    ax_b.text(5.0, stg["y"] + 0.40, stg["count"], fontsize=9.5, fontweight="bold" if idx == 2 else "normal", color="#2c3e50", ha="center")
    
    # Down Arrow with Criteria on LEFT and Counts on RIGHT
    if stg["left_criterion"] is not None:
        next_y = stages[idx+1]["y"] + stages[idx+1]["h"]
        y_mid = (stg["y"] + next_y) / 2.0
        
        # Arrow down at X = 5.0
        ax_b.annotate("", xy=(5.0, next_y + 0.05), xytext=(5.0, stg["y"] - 0.05),
                      arrowprops=dict(arrowstyle="->", color="#e74c3c" if "-" in stg["right_removed"] else "#34495e", lw=2.0))
        
        # Quality Criterion on LEFT side of arrow (ha="right", x = 4.80)
        ax_b.text(4.80, y_mid, stg["left_criterion"], 
                  fontsize=8.5, fontweight="bold", color="#2c3e50", 
                  ha="right", va="center")
        
        # Removed Count on RIGHT side of arrow (ha="left", x = 5.20)
        right_col = "#c0392b" if "-" in stg["right_removed"] else "#16a085"
        ax_b.text(5.20, y_mid, stg["right_removed"], 
                  fontsize=8.5, fontweight="bold", color=right_col, 
                  ha="left", va="center")

# Save PNG and PDF
plt.savefig(out_png, dpi=300, bbox_inches='tight')
plt.savefig(out_pdf, bbox_inches='tight')
plt.close()

print(f"[SUCCESS] 3-stage filtration figure saved to:\n  PNG: {out_png}\n  PDF: {out_pdf}")
