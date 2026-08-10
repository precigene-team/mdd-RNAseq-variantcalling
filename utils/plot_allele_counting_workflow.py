#!/usr/bin/env python3
# [[utils/plot_allele_counting_workflow.py]]
# Goal: Visualize allele counting workflow with LaTeX formula rendering for Slide 9.

import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches

# Enable LaTeX mathtext formatting in matplotlib
plt.rcParams['mathtext.fontset'] = 'dejavusans'

# --- 1. Dynamic Pathing ---
def get_project_root():
    curr = os.getcwd()
    while curr != "/" and not os.path.exists(os.path.join(curr, "run")):
        curr = os.path.dirname(curr)
    return curr

project_root = get_project_root()
output_dir = os.path.join(project_root, "results", "utils", "figures")
os.makedirs(output_dir, exist_ok=True)

out_png = os.path.join(output_dir, "allele_counting_workflow.png")
out_pdf = os.path.join(output_dir, "allele_counting_workflow.pdf")

# --- 2. Figure Setup ---
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15.5, 5.6), dpi=300, facecolor="#ffffff")
plt.rcParams['font.sans-serif'] = 'DejaVu Sans'
plt.subplots_adjust(left=0.04, right=0.96, bottom=0.06, top=0.90, wspace=0.22)

# ==============================================================================
# PANEL 1: Spatial Intersection (Enhanced with Peak Curve & Clean Genomic Axis)
# ==============================================================================
ax1.set_facecolor('#ffffff')
ax1.set_xlim(0, 10)
ax1.set_ylim(0, 10)
ax1.axis('off')
ax1.set_title("1. Spatial Intersection", fontsize=12.5, fontweight='bold', color='#1a252f', pad=12)

# Outer box
box1 = patches.FancyBboxPatch((0.3, 0.4), 9.4, 8.8, boxstyle="round,pad=0.15",
                              facecolor="#fafbfc", edgecolor="#d0d7de", lw=1.2)
ax1.add_patch(box1)

ax1.text(5.0, 8.5, "TRESS Peak Calling & Intersect", fontsize=10, fontweight="bold", color="#d35400", ha="center")

# Draw Peak Curve
x_peak = np.linspace(1.0, 9.0, 300)
y_curve = 5.8 + 2.2 * np.exp(-((x_peak - 5.0) / 1.4)**2)

ax1.plot(x_peak, y_curve, color="#d35400", lw=2.2, zorder=3)
ax1.fill_between(x_peak, 5.8, y_curve, color="#f39c12", alpha=0.25, zorder=2)

# Highlighted TRESS Peak Interval Shading
ax1.axvspan(3.2, 6.8, ymin=0.42, ymax=0.82, color="#fef9e7", alpha=0.5, zorder=1)
ax1.plot([3.2, 6.8], [5.8, 5.8], color="#d35400", lw=2, zorder=3)
ax1.text(5.0, 6.1, "m6A Peak Region (350 bp)", fontsize=8.5, fontweight="bold", color="#d35400", ha="center")

# Clean Vertical Projection lines dropping from peak to genomic axis
ax1.plot([3.2, 3.2], [4.2, 5.8], color="#f39c12", linestyle="--", lw=1.2, zorder=2)
ax1.plot([6.8, 6.8], [4.2, 5.8], color="#f39c12", linestyle="--", lw=1.2, zorder=2)

# Genomic Axis Track (Solid line at y = 4.2)
ax1.plot([0.8, 9.2], [4.2, 4.2], color="#34495e", lw=2.2, zorder=3)

for tick_x in [1.5, 3.2, 5.0, 6.8, 8.5]:
    ax1.plot([tick_x, tick_x], [4.1, 4.3], color="#34495e", lw=1.5, zorder=3)

ax1.text(5.0, 3.5, "Genomic Axis (chr12: 52,140,000 - 52,141,000)", 
         fontsize=8.5, fontweight="bold", color="#2c3e50", ha="center")

# Het SNPs on Genomic Axis
ax1.plot(4.0, 4.2, marker="o", markersize=9, color="#27ae60", zorder=5)
ax1.plot(6.0, 4.2, marker="o", markersize=9, color="#27ae60", zorder=5)
ax1.text(5.0, 4.5, "Het SNPs (Inside Peak)", fontsize=8, fontweight="bold", color="#27ae60", ha="center")

ax1.plot(1.5, 4.2, marker="o", markersize=8, color="#95a5a6", zorder=5)
ax1.plot(8.5, 4.2, marker="o", markersize=8, color="#95a5a6", zorder=5)

# Bottom Callout Box
res_badge1 = patches.FancyBboxPatch((0.8, 0.7), 8.4, 2.1, boxstyle="round,pad=0.1",
                                    facecolor="#e8f8f5", edgecolor="#17a589", lw=1.4)
ax1.add_patch(res_badge1)
ax1.text(5.0, 2.2, "bedtools intersect", fontsize=9.5, fontweight="bold", color="#117864", ha="center")
ax1.text(5.0, 1.3, "Isolates Candidate ASm6A Sites\nOnly retains Het SNPs sitting inside TRESS m6A peaks", 
         fontsize=8.5, color="#148f77", ha="center", va="center")

# ==============================================================================
# PANEL 2: Pileup Allele Counting (ASEReadCounter)
# ==============================================================================
ax2.set_facecolor('#ffffff')
ax2.set_xlim(0, 10)
ax2.set_ylim(0, 10)
ax2.axis('off')
ax2.set_title("2. Pileup Allele Counting", fontsize=12.5, fontweight='bold', color='#1a252f', pad=12)

box2 = patches.FancyBboxPatch((0.3, 0.4), 9.4, 8.8, boxstyle="round,pad=0.15",
                              facecolor="#fafbfc", edgecolor="#d0d7de", lw=1.2)
ax2.add_patch(box2)

ax2.text(5.0, 8.5, "GATK ASEReadCounter", fontsize=10, fontweight="bold", color="#1b4f72", ha="center")

qc_badge = patches.FancyBboxPatch((1.5, 7.3), 7.0, 0.85, boxstyle="round,pad=0.05",
                                  facecolor="#ebf5fb", edgecolor="#2980b9", lw=1.2)
ax2.add_patch(qc_badge)
ax2.text(5.0, 7.7, "Quality Filters: MapQ >= 20, BaseQ >= 20", fontsize=8.5, fontweight="bold", color="#1b4f72", ha="center", va="center")

# Input Fraction Pileup Box
b_in = patches.Rectangle((0.8, 3.4), 4.0, 3.4, facecolor="#ffffff", edgecolor="#16a085", lw=1.2)
ax2.add_patch(b_in)
ax2.text(2.8, 6.4, "Input BAM (Total)", fontsize=9, fontweight="bold", color="#16a085", ha="center")
ax2.text(2.8, 5.5, "REF (A): 15 reads", fontsize=8.5, color="#27ae60", ha="center")
ax2.text(2.8, 4.7, "ALT (G): 14 reads", fontsize=8.5, color="#d35400", ha="center")
ax2.text(2.8, 3.8, "Ratio ~ 50:50", fontsize=8.5, fontweight="bold", color="#566573", ha="center")

# IP Fraction Pileup Box
b_ip = patches.Rectangle((5.2, 3.4), 4.0, 3.4, facecolor="#ffffff", edgecolor="#1f4e79", lw=1.2)
ax2.add_patch(b_ip)
ax2.text(7.2, 6.4, "IP BAM (m6A Pool)", fontsize=9, fontweight="bold", color="#1f4e79", ha="center")
ax2.text(7.2, 5.5, "REF (A):  4 reads", fontsize=8.5, color="#27ae60", ha="center")
ax2.text(7.2, 4.7, "ALT (G): 28 reads", fontsize=8.5, color="#d35400", ha="center")
ax2.text(7.2, 3.8, "Ratio = 87.5% ALT", fontsize=8.5, fontweight="bold", color="#d35400", ha="center")

# Result Box
res_badge2 = patches.FancyBboxPatch((0.8, 0.7), 8.4, 2.1, boxstyle="round,pad=0.1",
                                    facecolor="#fef9e7", edgecolor="#f39c12", lw=1.4)
ax2.add_patch(res_badge2)
ax2.text(5.0, 2.2, "Allele Count Extraction", fontsize=9.5, fontweight="bold", color="#d35400", ha="center")
ax2.text(5.0, 1.3, "Extracts REF and ALT counts per subject\nPaired Input baseline & IP methylated reads", 
         fontsize=8.5, color="#7e5109", ha="center", va="center")

# ==============================================================================
# PANEL 3: Allelic Ratio Testing (Rendered with LaTeX Math Formula)
# ==============================================================================
ax3.set_facecolor('#ffffff')
ax3.set_xlim(0, 10)
ax3.set_ylim(0, 10)
ax3.axis('off')
ax3.set_title("3. Allelic Ratio Testing", fontsize=12.5, fontweight='bold', color='#1a252f', pad=12)

box3 = patches.FancyBboxPatch((0.3, 0.4), 9.4, 8.8, boxstyle="round,pad=0.15",
                              facecolor="#fafbfc", edgecolor="#d0d7de", lw=1.2)
ax3.add_patch(box3)

# Formula Box with rendered LaTeX math
ax3.text(5.0, 8.5, "ALT Allele Fraction Formula", fontsize=10, fontweight="bold", color="#6c3483", ha="center")
f_box = patches.FancyBboxPatch((1.2, 6.7), 7.6, 1.4, boxstyle="round,pad=0.1",
                               facecolor="#f5eef8", edgecolor="#8e44ad", lw=1.4)
ax3.add_patch(f_box)

# Render crisp LaTeX formula
latex_formula = r"$\mathrm{ALT}_{\mathrm{af}} = \frac{\mathrm{ALT}}{\mathrm{REF} + \mathrm{ALT}}$"
ax3.text(5.0, 7.4, latex_formula, fontsize=14, color="#6c3483", ha="center", va="center")

# Statistical Test Box
ax3.text(5.0, 6.1, "IP vs. Input Binomial Test", fontsize=9.5, fontweight="bold", color="#2c3e50", ha="center")

test_box = patches.Rectangle((0.8, 3.4), 8.4, 2.3, facecolor="#ffffff", edgecolor="#bdc3c7", lw=1.0)
ax3.add_patch(test_box)
ax3.text(5.0, 5.1, "Binomial Test on IP vs. Input Ratios", fontsize=8.5, fontweight="bold", color="#2c3e50", ha="center")

latex_q = r"$\text{Significance: } q < 0.05$"
ax3.text(5.0, 4.3, latex_q, fontsize=9, color="#27ae60", ha="center")

latex_eff = r"$\text{Effect Size: } |\mathrm{IP}_{\mathrm{af}} - 0.5| \geq 0.2$"
ax3.text(5.0, 3.7, latex_eff, fontsize=9, color="#d35400", ha="center")

# Prioritization Box
res_badge3 = patches.FancyBboxPatch((0.8, 0.7), 8.4, 2.1, boxstyle="round,pad=0.1",
                                    facecolor="#f5eef8", edgecolor="#8e44ad", lw=1.4)
ax3.add_patch(res_badge3)
ax3.text(5.0, 2.2, "Candidate ASm6A Prioritization", fontsize=9.5, fontweight="bold", color="#6c3483", ha="center")
ax3.text(5.0, 1.3, "Filters loci with significant allelic imbalance\nPrioritizes high effect size m6A methylation sites", 
         fontsize=8.5, color="#512e5f", ha="center", va="center")

# Save PNG and PDF
plt.savefig(out_png, dpi=300, bbox_inches='tight')
plt.savefig(out_pdf, bbox_inches='tight')
plt.close()

print(f"[SUCCESS] Allele counting workflow figure with LaTeX math saved to:\n  PNG: {out_png}\n  PDF: {out_pdf}")
