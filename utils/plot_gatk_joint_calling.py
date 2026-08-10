#!/usr/bin/env python3
# [[utils/plot_gatk_joint_calling.py]]
# Goal: Visualize GATK HaplotypeCaller GVCF mode and Cohort Joint Genotyping for Slide 5.

import os
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

out_png = os.path.join(output_dir, "gatk_joint_calling.png")
out_pdf = os.path.join(output_dir, "gatk_joint_calling.pdf")

# --- 2. Figure Setup ---
fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(14, 5.6), dpi=300, facecolor="#ffffff")
plt.rcParams['font.sans-serif'] = 'DejaVu Sans'
plt.subplots_adjust(left=0.04, right=0.97, bottom=0.08, top=0.88, wspace=0.18)

# ==============================================================================
# PANEL A: GVCF vs Regular Mode
# ==============================================================================
ax_a.set_facecolor('#ffffff')
ax_a.set_xlim(0, 10)
ax_a.set_ylim(0, 10)
ax_a.axis('off')
ax_a.set_title("A. HaplotypeCaller: Regular vs. GVCF Mode (-ERC GVCF)", 
               fontsize=12, fontweight='bold', color='#1a252f', pad=12)

# Box for Regular Mode
b_reg = patches.FancyBboxPatch((0.4, 5.2), 9.2, 4.0, boxstyle="round,pad=0.2",
                               facecolor="#fafbfc", edgecolor="#bdc3c7", lw=1.2)
ax_a.add_patch(b_reg)
ax_a.text(0.8, 8.6, "Regular VCF Mode", fontsize=10.5, fontweight="bold", color="#7f8c8d")
ax_a.text(0.8, 8.1, "Emits ONLY positions with variants (drops reference sites)", fontsize=8.5, color="#566573")

# Track for Regular Mode
ax_a.plot([1.0, 9.0], [6.5, 6.5], color="#bdc3c7", lw=2, linestyle=":")
ax_a.add_patch(patches.Rectangle((4.5, 6.1), 1.0, 0.8, facecolor="#e74c3c"))
ax_a.text(5.0, 6.5, "SNP", fontsize=8.5, color="#ffffff", ha="center", va="center", fontweight="bold")
ax_a.text(2.5, 5.7, "Missing Data / No Record", fontsize=8, color="#95a5a6", fontstyle="italic")
ax_a.text(7.5, 5.7, "Missing Data / No Record", fontsize=8, color="#95a5a6", fontstyle="italic")

# Box for GVCF Mode
b_gvcf = patches.FancyBboxPatch((0.4, 0.4), 9.2, 4.2, boxstyle="round,pad=0.2",
                                facecolor="#fef9e7", edgecolor="#f39c12", lw=1.5)
ax_a.add_patch(b_gvcf)
ax_a.text(0.8, 4.0, "GVCF Mode (-ERC GVCF)  [OUR CHOICE]", fontsize=10.5, fontweight="bold", color="#d35400")
ax_a.text(0.8, 3.5, "Emits variants AND confident reference blocks (<NON_REF>)", fontsize=8.5, color="#7e5109")

# Track for GVCF Mode
ax_a.add_patch(patches.Rectangle((1.0, 1.8), 3.2, 0.8, facecolor="#27ae60", alpha=0.85))
ax_a.text(2.6, 2.2, "Ref Block (0/0)", fontsize=8.5, color="#ffffff", ha="center", va="center", fontweight="bold")

ax_a.add_patch(patches.Rectangle((4.5, 1.8), 1.0, 0.8, facecolor="#e74c3c"))
ax_a.text(5.0, 2.2, "SNP", fontsize=8.5, color="#ffffff", ha="center", va="center", fontweight="bold")

ax_a.add_patch(patches.Rectangle((5.8, 1.8), 3.2, 0.8, facecolor="#27ae60", alpha=0.85))
ax_a.text(7.4, 2.2, "Ref Block (0/0)", fontsize=8.5, color="#ffffff", ha="center", va="center", fontweight="bold")

ax_a.text(5.0, 0.9, "Guarantees complete coverage info across all 39 subjects", 
          fontsize=9, fontweight="bold", color="#d35400", ha="center")

# ==============================================================================
# PANEL B: 39 GVCFs -> 1 Cohort VCF Merge (Clipping Fixed)
# ==============================================================================
ax_b.set_facecolor('#ffffff')
ax_b.set_xlim(0, 10.2)
ax_b.set_ylim(0, 10)
ax_b.axis('off')
ax_b.set_title("B. Joint Genotyping: 39 GVCFs -> 1 Cohort VCF Matrix", 
               fontsize=12, fontweight='bold', color='#1a252f', pad=12)

# Individual GVCF stacks
ax_b.text(1.5, 8.8, "39 Individual GVCFs", fontsize=10, fontweight="bold", color="#2c3e50", ha="center")

gvcf_colors = ["#2980b9", "#16a085", "#8e44ad", "#d35400"]
for i in range(4):
    y_pos = 7.5 - i * 0.7
    lbl = f"Sample_{i+1:02d}.g.vcf" if i < 3 else "... Sample_39.g.vcf"
    box = patches.FancyBboxPatch((0.1, y_pos), 2.8, 0.55, boxstyle="round,pad=0.05",
                                 facecolor=gvcf_colors[i], edgecolor="none", alpha=0.85)
    ax_b.add_patch(box)
    ax_b.text(1.5, y_pos + 0.275, lbl, fontsize=8.5, color="#ffffff", ha="center", va="center", fontweight="bold")

# Merge Arrow
ax_b.annotate("", xy=(4.4, 6.2), xytext=(3.1, 6.2),
            arrowprops=dict(arrowstyle="->", color="#34495e", lw=2.5))
ax_b.text(3.75, 6.6, "CombineGVCFs\nGenotypeGVCFs", fontsize=8, fontweight="bold", color="#2c3e50", ha="center")

# Final Cohort VCF Matrix Header
ax_b.text(7.25, 8.8, "Cohort VCF Matrix", fontsize=10, fontweight="bold", color="#117864", ha="center")

# Table Container Box (Well within 10.2 right boundary)
matrix_box = patches.FancyBboxPatch((4.5, 0.5), 5.3, 7.8, boxstyle="round,pad=0.15",
                                   facecolor="#e8f8f5", edgecolor="#17a589", lw=1.5)
ax_b.add_patch(matrix_box)

# Column Coordinates
x_pos = 4.8
x_s01 = 6.4
x_s02 = 7.3
x_s03 = 8.2
x_dot = 9.1

# Table Column Headers
ax_b.text(x_pos, 7.7, "POS", fontsize=8.5, fontweight="bold", color="#117864")
ax_b.text(x_s01, 7.7, "S01", fontsize=8.5, fontweight="bold", color="#117864", ha="center")
ax_b.text(x_s02, 7.7, "S02", fontsize=8.5, fontweight="bold", color="#117864", ha="center")
ax_b.text(x_s03, 7.7, "S03", fontsize=8.5, fontweight="bold", color="#117864", ha="center")
ax_b.text(x_dot, 7.7, "...", fontsize=8.5, fontweight="bold", color="#117864", ha="center")

ax_b.plot([4.7, 9.6], [7.45, 7.45], color="#17a589", lw=1.2)

# Matrix Rows
rows_data = [
    ("chr1:1042", "0/0", "0/1", "0/0"),
    ("chr1:1058", "0/1", "0/1", "1/1"),
    ("chr2:4012", "0/0", "0/0", "0/1"),
    ("chr3:8901", "1/1", "0/1", "0/0"),
    ("...", "...", "...", "...")
]

y_row = 6.7
for pos, s1, s2, s3 in rows_data:
    ax_b.text(x_pos, y_row, pos, fontsize=8, color="#2c3e50")
    
    c1 = "#d35400" if s1 == "0/1" else "#2c3e50"
    c2 = "#d35400" if s2 == "0/1" else "#2c3e50"
    c3 = "#d35400" if s3 == "0/1" else "#2c3e50"
    
    ax_b.text(x_s01, y_row, s1, fontsize=8, fontweight="bold" if s1 == "0/1" else "normal", color=c1, ha="center")
    ax_b.text(x_s02, y_row, s2, fontsize=8, fontweight="bold" if s2 == "0/1" else "normal", color=c2, ha="center")
    ax_b.text(x_s03, y_row, s3, fontsize=8, fontweight="bold" if s3 == "0/1" else "normal", color=c3, ha="center")
    ax_b.text(x_dot, y_row, "...", fontsize=8, color="#7f8c8d", ha="center")
    y_row -= 0.95

# Bottom Legend inside table container
ax_b.text(7.15, 0.9, "Genotypes: 0/0 (Ref), 0/1 (Het), 1/1 (Alt)", 
          fontsize=8.5, fontweight="bold", color="#117864", ha="center")

# Save PNG and PDF with tight bounding box
plt.savefig(out_png, dpi=300, bbox_inches='tight')
plt.savefig(out_pdf, bbox_inches='tight')
plt.close()

print(f"[SUCCESS] Fixed clipping figure saved to:\n  PNG: {out_png}\n  PDF: {out_pdf}")
