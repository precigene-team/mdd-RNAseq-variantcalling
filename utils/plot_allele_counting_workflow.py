#!/usr/bin/env python3
# [[utils/plot_allele_counting_workflow.py]]
# Goal: Minimalist, clean 3-step ASm6A workflow diagram with worked numerical example in Card 3 based on Card 2 counts for Slide 9.

import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches

# Enable Matplotlib text formatting
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
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15.5, 5.2), dpi=300, facecolor="#ffffff")
plt.rcParams['font.sans-serif'] = 'DejaVu Sans'
plt.subplots_adjust(left=0.03, right=0.97, bottom=0.06, top=0.88, wspace=0.22)

# ==============================================================================
# CARD 1: 1. Spatial Intersection
# ==============================================================================
ax1.set_facecolor('#ffffff')
ax1.set_xlim(0, 10)
ax1.set_ylim(0, 10)
ax1.axis('off')
ax1.set_title("1. Spatial Intersection", fontsize=13, fontweight='bold', color='#1a252f', pad=12)

# Outer Card
card1 = patches.FancyBboxPatch((0.4, 0.4), 9.2, 9.0, boxstyle="round,pad=0.15",
                              facecolor="#ffffff", edgecolor="#d0d7de", lw=1.4)
ax1.add_patch(card1)

# INPUT Pill
in_pill1 = patches.FancyBboxPatch((1.0, 8.0), 8.0, 0.85, boxstyle="round,pad=0.05",
                                  facecolor="#eef2f7", edgecolor="#5c6b73", lw=1.2)
ax1.add_patch(in_pill1)
ax1.text(5.0, 8.42, "INPUT: Het SNPs (VCF) + m6A Peaks (BED)", fontsize=8.5, fontweight="bold", color="#2c3e50", ha="center", va="center")

# Peak Curve Visual
x_peak = np.linspace(1.5, 8.5, 200)
y_curve = 5.0 + 2.0 * np.exp(-((x_peak - 5.0) / 1.3)**2)

ax1.plot(x_peak, y_curve, color="#d35400", lw=2.2, zorder=3)
ax1.fill_between(x_peak, 5.0, y_curve, color="#f39c12", alpha=0.25, zorder=2)

# Genomic Axis Line & SNP
ax1.plot([1.0, 9.0], [4.0, 4.0], color="#34495e", lw=2.2, zorder=3)
ax1.plot(5.0, 4.0, marker="o", markersize=10, color="#27ae60", zorder=5)
ax1.text(5.0, 4.4, "Het SNP inside Peak", fontsize=8.5, fontweight="bold", color="#27ae60", ha="center")

# PROCESS Text
ax1.text(5.0, 2.6, "bedtools intersect", fontsize=9, fontweight="bold", color="#117864", ha="center",
         bbox=dict(boxstyle="round,pad=0.3", fc="#e8f8f5", ec="#17a589", lw=1.0))

# OUTPUT Pill
out_pill1 = patches.FancyBboxPatch((1.0, 0.8), 8.0, 0.9, boxstyle="round,pad=0.05",
                                   facecolor="#d4efdf", edgecolor="#27ae60", lw=1.4)
ax1.add_patch(out_pill1)
ax1.text(5.0, 1.25, "OUTPUT: Candidate m6A-SNPs", fontsize=9, fontweight="bold", color="#1e8449", ha="center", va="center")

# ==============================================================================
# CARD 2: 2. Pileup Allele Counting
# ==============================================================================
ax2.set_facecolor('#ffffff')
ax2.set_xlim(0, 10)
ax2.set_ylim(0, 10)
ax2.axis('off')
ax2.set_title("2. Pileup Allele Counting", fontsize=13, fontweight='bold', color='#1a252f', pad=12)

card2 = patches.FancyBboxPatch((0.4, 0.4), 9.2, 9.0, boxstyle="round,pad=0.15",
                              facecolor="#ffffff", edgecolor="#d0d7de", lw=1.4)
ax2.add_patch(card2)

# INPUT Pill
in_pill2 = patches.FancyBboxPatch((1.0, 8.0), 8.0, 0.85, boxstyle="round,pad=0.05",
                                  facecolor="#eef2f7", edgecolor="#5c6b73", lw=1.2)
ax2.add_patch(in_pill2)
ax2.text(5.0, 8.42, "INPUT: Paired Input & IP BAMs", fontsize=8.5, fontweight="bold", color="#2c3e50", ha="center", va="center")

# Input Box
b_in = patches.FancyBboxPatch((0.9, 4.0), 3.9, 3.4, boxstyle="round,pad=0.05", facecolor="#f4f6f7", edgecolor="#16a085", lw=1.2)
ax2.add_patch(b_in)
ax2.text(2.85, 6.7, "Input (Total)", fontsize=9, fontweight="bold", color="#16a085", ha="center")
ax2.text(2.85, 5.6, "REF: 15 reads\nALT: 14 reads", fontsize=8, color="#2c3e50", ha="center")
ax2.text(2.85, 4.5, "48.3% ALT", fontsize=8.5, fontweight="bold", color="#566573", ha="center")

# IP Box
b_ip = patches.FancyBboxPatch((5.2, 4.0), 3.9, 3.4, boxstyle="round,pad=0.05", facecolor="#fef9e7", edgecolor="#d35400", lw=1.2)
ax2.add_patch(b_ip)
ax2.text(7.15, 6.7, "IP (m6A Pool)", fontsize=9, fontweight="bold", color="#d35400", ha="center")
ax2.text(7.15, 5.6, "REF:  4 reads\nALT: 28 reads", fontsize=8, color="#2c3e50", ha="center")
ax2.text(7.15, 4.5, "87.5% ALT", fontsize=8.5, fontweight="bold", color="#d35400", ha="center")

# PROCESS Text
ax2.text(5.0, 2.6, "GATK ASEReadCounter", fontsize=9, fontweight="bold", color="#1b4f72", ha="center",
         bbox=dict(boxstyle="round,pad=0.3", fc="#ebf5fb", ec="#2980b9", lw=1.0))

# OUTPUT Pill
out_pill2 = patches.FancyBboxPatch((1.0, 0.8), 8.0, 0.9, boxstyle="round,pad=0.05",
                                   facecolor="#fdebd0", edgecolor="#d35400", lw=1.4)
ax2.add_patch(out_pill2)
ax2.text(5.0, 1.25, "OUTPUT: Paired Count Matrix", fontsize=9, fontweight="bold", color="#a04000", ha="center", va="center")

# ==============================================================================
# CARD 3: 3. Allelic Ratio Testing (Worked Example from Card 2)
# ==============================================================================
ax3.set_facecolor('#ffffff')
ax3.set_xlim(0, 10)
ax3.set_ylim(0, 10)
ax3.axis('off')
ax3.set_title("3. Allelic Ratio Testing", fontsize=13, fontweight='bold', color='#1a252f', pad=12)

card3 = patches.FancyBboxPatch((0.4, 0.4), 9.2, 9.0, boxstyle="round,pad=0.15",
                              facecolor="#ffffff", edgecolor="#d0d7de", lw=1.4)
ax3.add_patch(card3)

# INPUT Pill
in_pill3 = patches.FancyBboxPatch((1.0, 8.0), 8.0, 0.85, boxstyle="round,pad=0.05",
                                  facecolor="#eef2f7", edgecolor="#5c6b73", lw=1.2)
ax3.add_patch(in_pill3)
ax3.text(5.0, 8.42, "INPUT: Allele Count Matrix", fontsize=8.5, fontweight="bold", color="#2c3e50", ha="center", va="center")

# LaTeX Worked Example Box
f_box = patches.FancyBboxPatch((0.9, 4.0), 8.2, 3.5, boxstyle="round,pad=0.08",
                               facecolor="#f5eef8", edgecolor="#8e44ad", lw=1.2)
ax3.add_patch(f_box)

ax3.text(5.0, 7.0, "Example Calculation:", fontsize=9, fontweight="bold", color="#6c3483", ha="center")

# Formula & Calculation
latex_calc = r"$\mathrm{IP}_{\mathrm{af}} = \frac{28}{4 + 28} = \mathbf{0.875\ (87.5\%)}$"
ax3.text(5.0, 6.0, latex_calc, fontsize=11, color="#6c3483", ha="center", va="center")

latex_skew = r"$\Delta\mathrm{Skew} = |0.875 - 0.483| = \mathbf{39.2\% > 20\%}$"
ax3.text(5.0, 4.8, latex_skew, fontsize=10, color="#d35400", ha="center", va="center")

# PROCESS Text
ax3.text(5.0, 2.6, "Binomial Test (q < 0.05)", fontsize=9, fontweight="bold", color="#27ae60", ha="center",
         bbox=dict(boxstyle="round,pad=0.3", fc="#e8f8f5", ec="#27ae60", lw=1.0))

# OUTPUT Pill
out_pill3 = patches.FancyBboxPatch((1.0, 0.8), 8.0, 0.9, boxstyle="round,pad=0.05",
                                   facecolor="#e8daef", edgecolor="#7d3c98", lw=1.4)
ax3.add_patch(out_pill3)
ax3.text(5.0, 1.25, "OUTPUT: Candidate ASm6A Loci", fontsize=9, fontweight="bold", color="#512e5f", ha="center", va="center")

# ==============================================================================
# CONNECTING ARROWS BETWEEN CARDS
# ==============================================================================
fig.text(0.352, 0.50, "➔", fontsize=22, fontweight="bold", color="#34495e", ha="center", va="center")
fig.text(0.665, 0.50, "➔", fontsize=22, fontweight="bold", color="#34495e", ha="center", va="center")

# Save PNG and PDF
plt.savefig(out_png, dpi=300, bbox_inches='tight')
plt.savefig(out_pdf, bbox_inches='tight')
plt.close()

print(f"[SUCCESS] Worked calculation example figure saved to:\n  PNG: {out_png}\n  PDF: {out_pdf}")
