#!/usr/bin/env python3
# [[utils/plot_gatk_preprocessing.py]]
# Goal: Visualize GATK RNA-seq preprocessing steps (Read Grouping, MarkDuplicates, SplitNCigarReads)
# + HaplotypeCaller variant calling for Slides 4-5.

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

out_png = os.path.join(output_dir, "gatk_preprocessing.png")
out_pdf = os.path.join(output_dir, "gatk_preprocessing.pdf")

# --- 2. Figure Setup (2x2 Grid) ---
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(14, 9.2), dpi=300, facecolor="#ffffff")
plt.rcParams['font.sans-serif'] = 'DejaVu Sans'
plt.subplots_adjust(left=0.04, right=0.96, bottom=0.06, top=0.90, wspace=0.22, hspace=0.28)

def setup_axes(ax):
    ax.set_facecolor('#fafbfc')
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.axis('off')
    return patches.FancyBboxPatch((0.4, 0.4), 9.2, 8.8, boxstyle="round,pad=0.2",
                                  facecolor="#ffffff", edgecolor="#d0d7de", lw=1.2)

# ==============================================================================
# STEP 1: Read Grouping (AddOrReplaceReadGroups)
# ==============================================================================
box1 = setup_axes(ax1)
ax1.add_patch(box1)
ax1.set_title("1. Read Grouping\n(AddOrReplaceReadGroups)", fontsize=11.5, fontweight='bold', color='#1a252f')

ax1.text(5.0, 8.2, "Assign Metadata Tags (RG)", fontsize=10, fontweight='bold', color='#2c3e50', ha='center')

# Read 1 (Input)
r1 = patches.Rectangle((1.2, 5.8), 5.5, 0.8, facecolor="#16a085", alpha=0.85)
ax1.add_patch(r1)
ax1.text(3.95, 6.2, "Read #1042 (RNA-seq)", fontsize=9, color="#ffffff", ha="center", va="center", fontweight="bold")
tag1 = patches.FancyBboxPatch((7.0, 5.8), 2.2, 0.8, boxstyle="round,pad=0.05", facecolor="#e8f8f5", edgecolor="#16a085", lw=1.2)
ax1.add_patch(tag1)
ax1.text(8.1, 6.2, "RG:Input_01", fontsize=8.5, fontweight="bold", color="#117864", ha="center", va="center")

# Read 2 (IP)
r2 = patches.Rectangle((1.2, 4.2), 5.5, 0.8, facecolor="#1f4e79", alpha=0.85)
ax1.add_patch(r2)
ax1.text(3.95, 4.6, "Read #1043 (MeRIP)", fontsize=9, color="#ffffff", ha="center", va="center", fontweight="bold")
tag2 = patches.FancyBboxPatch((7.0, 4.2), 2.2, 0.8, boxstyle="round,pad=0.05", facecolor="#ebf5fb", edgecolor="#1f4e79", lw=1.2)
ax1.add_patch(tag2)
ax1.text(8.1, 4.6, "RG:IP_01", fontsize=8.5, fontweight="bold", color="#1b4f72", ha="center", va="center")

ax1.text(5.0, 1.8, "Tracks sample origin & keeps\nInput and IP paired without mixing", 
         fontsize=9, color="#566573", ha="center", va="center", bbox=dict(boxstyle="round,pad=0.4", fc="#f2f4f4", ec="none"))

# ==============================================================================
# STEP 2: MarkDuplicates
# ==============================================================================
box2 = setup_axes(ax2)
ax2.add_patch(box2)
ax2.set_title("2. Mark Duplicates\n(MarkDuplicates)", fontsize=11.5, fontweight='bold', color='#1a252f')

ax2.text(5.0, 8.2, "Flag PCR Copies (Do NOT Delete)", fontsize=10, fontweight='bold', color='#2c3e50', ha='center')

r_orig = patches.Rectangle((1.2, 6.2), 5.5, 0.7, facecolor="#2980b9", alpha=0.9)
ax2.add_patch(r_orig)
ax2.text(3.95, 6.55, "Unique Read (Original)", fontsize=8.5, color="#ffffff", ha="center", va="center", fontweight="bold")

r_dup1 = patches.Rectangle((1.2, 5.0), 5.5, 0.7, facecolor="#95a5a6", alpha=0.7)
ax2.add_patch(r_dup1)
ax2.text(3.95, 5.35, "PCR Duplicate #1", fontsize=8.5, color="#ffffff", ha="center", va="center")
dup_tag1 = patches.Rectangle((7.0, 5.0), 2.0, 0.7, facecolor="#f39c12")
ax2.add_patch(dup_tag1)
ax2.text(8.0, 5.35, "FLAG:1024", fontsize=8, color="#ffffff", ha="center", va="center", fontweight="bold")

r_dup2 = patches.Rectangle((1.2, 3.8), 5.5, 0.7, facecolor="#95a5a6", alpha=0.7)
ax2.add_patch(r_dup2)
ax2.text(3.95, 4.15, "PCR Duplicate #2", fontsize=8.5, color="#ffffff", ha="center", va="center")
dup_tag2 = patches.Rectangle((7.0, 3.8), 2.0, 0.7, facecolor="#f39c12")
ax2.add_patch(dup_tag2)
ax2.text(8.0, 4.15, "FLAG:1024", fontsize=8, color="#ffffff", ha="center", va="center", fontweight="bold")

ax2.text(5.0, 1.8, "Marked, not deleted!\nDeletions would distort IP:Input ratios", 
         fontsize=9, fontweight="bold", color="#d35400", ha="center", va="center", 
         bbox=dict(boxstyle="round,pad=0.4", fc="#fef9e7", ec="#f39c12", lw=1.0))

# ==============================================================================
# STEP 3: SplitNCigarReads
# ==============================================================================
box3 = setup_axes(ax3)
ax3.add_patch(box3)
ax3.set_title("3. Split Intron Gaps\n(SplitNCigarReads)", fontsize=11.5, fontweight='bold', color='#1a252f')

ax3.text(5.0, 8.2, "Split Reads at Intron 'N' CIGAR Gaps", fontsize=10, fontweight='bold', color='#2c3e50', ha='center')

ax3.text(1.0, 6.7, "Before:", fontsize=9, fontweight="bold", color="#7f8c8d")
r_ex1 = patches.Rectangle((2.2, 6.4), 2.2, 0.7, facecolor="#8e44ad", alpha=0.85)
ax3.add_patch(r_ex1)
ax3.text(3.3, 6.75, "Exon 1", fontsize=8, color="#ffffff", ha="center", va="center", fontweight="bold")

ax3.plot([4.4, 6.6], [6.75, 6.75], color="#e74c3c", lw=2, linestyle=":")
ax3.text(5.5, 7.1, "Intron Gap (N)", fontsize=8, color="#c0392b", ha="center", fontweight="bold")

r_ex2 = patches.Rectangle((6.6, 6.4), 2.2, 0.7, facecolor="#8e44ad", alpha=0.85)
ax3.add_patch(r_ex2)
ax3.text(7.7, 6.75, "Exon 2", fontsize=8, color="#ffffff", ha="center", va="center", fontweight="bold")

ax3.annotate("", xy=(5.5, 5.0), xytext=(5.5, 5.9), arrowprops=dict(arrowstyle="->", color="#34495e", lw=1.8))

ax3.text(1.0, 4.1, "After:", fontsize=9, fontweight="bold", color="#27ae60")
r_split1 = patches.Rectangle((2.2, 3.8), 2.2, 0.7, facecolor="#27ae60", alpha=0.9)
ax3.add_patch(r_split1)
ax3.text(3.3, 4.15, "Piece 1", fontsize=8, color="#ffffff", ha="center", va="center", fontweight="bold")

r_split2 = patches.Rectangle((6.6, 3.8), 2.2, 0.7, facecolor="#27ae60", alpha=0.9)
ax3.add_patch(r_split2)
ax3.text(7.7, 4.15, "Piece 2", fontsize=8, color="#ffffff", ha="center", va="center", fontweight="bold")

ax3.text(5.0, 1.8, "Prevents mistaking intron gaps\nfor false genomic deletions", 
         fontsize=9, color="#566573", ha="center", va="center", bbox=dict(boxstyle="round,pad=0.4", fc="#f2f4f4", ec="none"))

# ==============================================================================
# STEP 4: HaplotypeCaller (GVCF mode)
# ==============================================================================
box4 = setup_axes(ax4)
ax4.add_patch(box4)
ax4.set_title("4. Haplotype Caller\n(GVCF mode)", fontsize=11.5, fontweight='bold', color='#1a252f')

ax4.text(5.0, 8.2, "Active Region: Reassemble Reads & Score Genotypes", fontsize=10, fontweight='bold', color='#2c3e50', ha='center')

# Reference sequence track with SNP column
ref_bases = ["C", "A", "G", "A", "C", "U"]
bx, bw, bgap, bh = 1.6, 0.9, 0.15, 0.7
snp_idx = 3
ref_y = 6.6
for i, letter in enumerate(ref_bases):
    x = bx + i * (bw + bgap)
    if i == snp_idx:
        face, edge, tc = "#fdecea", "#e74c3c", "#c0392b"
    else:
        face, edge, tc = "#e8f8f5", "#bdc3c7", "#1e8449"
    ax4.add_patch(patches.FancyBboxPatch((x, ref_y), bw, bh, boxstyle="round,pad=0.04",
                 facecolor=face, edgecolor=edge, lw=1.4))
    ax4.text(x + bw/2, ref_y + bh/2, letter, fontsize=8.5, fontweight="bold", color=tc, ha="center", va="center")
snp_x = bx + snp_idx * (bw + bgap) + bw/2
ax4.text(snp_x, ref_y + bh + 0.45, "SNP A\u2192G", fontsize=7.5, fontweight="bold", color="#c0392b", ha="center")

# Aligned reads (REF = blue, ALT = red at the SNP column)
reads = [
    ("G", 2.2, 6.0, "#fdecea", "#c0392b"),  # ALT
    ("A", 1.8, 5.5, "#ebf5fb", "#1b4f72"),  # REF
    ("G", 2.4, 5.0, "#fdecea", "#c0392b"),  # ALT
    ("G", 1.7, 4.5, "#fdecea", "#c0392b"),  # ALT
    ("A", 2.0, 4.0, "#ebf5fb", "#1b4f72"),  # REF
]
for letter, xw, y, face, tc in reads:
    ax4.add_patch(patches.Rectangle((xw, y), 3.2, 0.42, facecolor="#ecf0f1", edgecolor="#bdc3c7", lw=0.8))
    ax4.add_patch(patches.Rectangle((snp_x - 0.22, y + 0.01), 0.44, 0.40, facecolor=face, edgecolor=tc, lw=1.0))
    ax4.text(snp_x, y + 0.21, letter, fontsize=7.5, fontweight="bold", color=tc, ha="center", va="center")

# Arrow to genotype call
ax4.annotate("", xy=(5.0, 2.9), xytext=(5.0, 3.9), arrowprops=dict(arrowstyle="->", color="#34495e", lw=1.8))

# Genotype call box
gen = patches.FancyBboxPatch((2.2, 2.0), 5.6, 0.9, boxstyle="round,pad=0.08",
                             facecolor="#f5eef8", edgecolor="#8e44ad", lw=1.5)
ax4.add_patch(gen)
ax4.text(5.0, 2.45, "Genotype: 0/1 (Heterozygous)", fontsize=9, fontweight="bold", color="#6c3483", ha="center", va="center")

ax4.text(5.0, 1.1, "3\u00d7 G (ALT) vs 2\u00d7 A (REF)\nGVCF emits every position (variants + reference blocks)\nfor joint calling across the cohort", 
         fontsize=8, color="#566573", ha="center", va="center")

# Save PNG and PDF
plt.savefig(out_png, dpi=300, bbox_inches='tight')
plt.savefig(out_pdf, bbox_inches='tight')
plt.close()

print(f"[SUCCESS] Preprocessing + HaplotypeCaller figure saved to:\n  PNG: {out_png}\n  PDF: {out_pdf}")
