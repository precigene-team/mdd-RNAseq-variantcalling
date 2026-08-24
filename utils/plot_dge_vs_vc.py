#!/usr/bin/env python3
# [[utils/plot_dge_vs_vc.py]]
# Goal: Visualize conceptual comparison between Differential Gene Expression (DGE) vs Variant-Calling ASm6A Analysis.

import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches

# Enable Matplotlib LaTeX formatting
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

out_png = os.path.join(output_dir, "dge_vs_variant_calling.png")
out_pdf = os.path.join(output_dir, "dge_vs_variant_calling.pdf")

# --- 2. Figure Setup ---
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14.0, 6.2), dpi=300, facecolor="#ffffff")
plt.rcParams['font.sans-serif'] = 'DejaVu Sans'
plt.subplots_adjust(left=0.06, right=0.94, bottom=0.12, top=0.86, wspace=0.28)

# Figure Title
fig.suptitle("Comparison: Differential Gene Expression (DGE) vs. Variant-Calling (ASm6A)", 
             fontsize=14, fontweight="bold", color="#1a252f", y=0.96)

# ==============================================================================
# PANEL 1: Differential Gene Expression (DGE Focus)
# ==============================================================================
ax1.set_facecolor('#ffffff')

categories = ['Control', 'MDD']
total_counts = [100, 100]

bars1 = ax1.bar(categories, total_counts, color='#95a5a6', edgecolor='#34495e', width=0.42, linewidth=1.2)

ax1.set_ylim(0, 140)
ax1.set_ylabel("Gene-Level Read Depth", fontsize=11, fontweight='bold', color='#2c3e50')
ax1.set_title("1. Differential Gene Expression (DGE / DESeq2)", fontsize=12, fontweight='bold', color='#1a252f', pad=12)

# Value Labels on Bars
for bar in bars1:
    yval = bar.get_height()
    ax1.text(bar.get_x() + bar.get_width()/2.0, yval + 4, f"100 reads\n(Total)", ha='center', va='bottom', fontsize=9.5, fontweight='bold', color='#2c3e50')

# Spine styling
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.spines['left'].set_color('#bdc3c7')
ax1.spines['bottom'].set_color('#bdc3c7')
ax1.tick_params(axis='x', labelsize=11)

# Method Details Box for Panel 1
ax1.text(0.5, -0.22, "• Focus: Overall transcript abundance (log2FC)\n• Tool: DESeq2 (VST normalization)\n• Limitation: Aggregates maternal & paternal reads;\n  blind to opposing allelic shifts (log2FC = 0.0, p = 1.0)",
         fontsize=9, color='#2c3e50', ha='center', va='top', transform=ax1.transAxes,
         bbox=dict(boxstyle="round,pad=0.5", fc="#f4f6f7", ec="#bdc3c7", lw=1.3))

# ==============================================================================
# PANEL 2: Variant-Calling & Allele-Specific Methylation (ASm6A Focus)
# ==============================================================================
ax2.set_facecolor('#ffffff')

ref_counts = [50, 10]
alt_counts = [50, 90]

b_alt = ax2.bar(categories, alt_counts, color='#e74c3c', edgecolor='#c0392b', width=0.42, linewidth=1.2, label='ALT Allele (G)')
b_ref = ax2.bar(categories, ref_counts, bottom=alt_counts, color='#2980b9', edgecolor='#1b4f72', width=0.42, linewidth=1.2, label='REF Allele (A)')

ax2.set_ylim(0, 140)
ax2.set_ylabel("Allele-Specific Read Depth", fontsize=11, fontweight='bold', color='#2c3e50')
ax2.set_title("2. Variant-Calling & ASm6A (GATK + glmmTMB)", fontsize=12, fontweight='bold', color='#1a252f', pad=12)

# Value Labels inside Stacked Bars
ax2.text(0, 25, "ALT: 50\n(50%)", ha='center', va='center', fontsize=9, fontweight='bold', color='#ffffff')
ax2.text(0, 75, "REF: 50\n(50%)", ha='center', va='center', fontsize=9, fontweight='bold', color='#ffffff')

ax2.text(1, 45, "ALT: 90\n(90% Imbalance!)", ha='center', va='center', fontsize=9, fontweight='bold', color='#ffffff')
ax2.text(1, 95, "REF: 10", ha='center', va='center', fontsize=8, fontweight='bold', color='#ffffff')

ax2.legend(loc='upper right', frameon=False, fontsize=9.5)

# Spine styling
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax2.spines['left'].set_color('#bdc3c7')
ax2.spines['bottom'].set_color('#bdc3c7')
ax2.tick_params(axis='x', labelsize=11)

# Method Details Box for Panel 2
ax2.text(0.5, -0.22, "• Focus: Single-allele resolution at Het SNPs\n• Tool: GATK ASEReadCounter + glmmTMB model\n• Advantage: Uncovers hidden cis-regulatory methylation shifts\n  (FDR < 0.05, Delta Skew = 40% in BIN3, MGAT5, IER2)",
         fontsize=9, color='#145a32', ha='center', va='top', transform=ax2.transAxes,
         bbox=dict(boxstyle="round,pad=0.5", fc="#d4efdf", ec="#27ae60", lw=1.3))

# Save PNG and PDF
plt.savefig(out_png, dpi=300, bbox_inches='tight')
plt.savefig(out_pdf, bbox_inches='tight')
plt.close()

print(f"[SUCCESS] DGE vs Variant Calling figure saved to:\n  PNG: {out_png}\n  PDF: {out_pdf}")
