#!/usr/bin/env python3
# [[utils/plot_input_only_rationale.py]]
# Goal: Visual diagram explaining Input-only variant calling with bolded main workflow steps.

import os
import sys
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

out_png = os.path.join(output_dir, "input_only_rationale.png")
out_pdf = os.path.join(output_dir, "input_only_rationale.pdf")

# --- 2. Create Presentation Visual Diagram ---
fig, axes = plt.subplots(2, 1, figsize=(13.5, 9.8), dpi=300, facecolor="#ffffff")
plt.rcParams['font.sans-serif'] = 'DejaVu Sans'
plt.subplots_adjust(left=0.06, right=0.96, bottom=0.06, top=0.92, hspace=0.25)

fig.suptitle("Why Call Variants Exclusively from Input RNA-seq?", 
             fontsize=15.0, fontweight="bold", color="#1a252f", y=0.97)

# -----------------------------------------------------------------------------
# 1. LOCALIZED COVERAGE BIAS (Genomic Track Cartoon with Generous Margins)
# -----------------------------------------------------------------------------
ax1 = axes[0]
ax1.set_facecolor("#ffffff")
ax1.axis("off")
ax1.set_xlim(0, 100)
ax1.set_ylim(0, 11)

ax1.text(0, 10.2, "1. Localized Coverage Bias (IP vs. Input Read Distribution)", 
         fontsize=12.0, fontweight="bold", color="#1a252f")

# Exon & Intron Transcript Model with generous margins
ax1.plot([10, 90], [5.0, 5.0], color="#7f8c8d", lw=2.2, zorder=1) # Intron line
rect_exon1 = patches.Rectangle((14, 4.2), 20, 1.6, facecolor="#34495e", edgecolor="#2c3e50", zorder=2)
rect_exon2 = patches.Rectangle((45, 4.2), 20, 1.6, facecolor="#34495e", edgecolor="#2c3e50", zorder=2)
rect_exon3 = patches.Rectangle((76, 4.2), 14, 1.6, facecolor="#34495e", edgecolor="#2c3e50", zorder=2)
ax1.add_patch(rect_exon1)
ax1.add_patch(rect_exon2)
ax1.add_patch(rect_exon3)

ax1.text(24, 3.2, "Exon 1 (Unmethylated)", fontsize=9.0, fontweight="bold", color="#2c3e50", ha="center")
ax1.text(55, 3.2, "Exon 2 (m6A Methylated Peak)", fontsize=9.0, fontweight="bold", color="#c0392b", ha="center")
ax1.text(83, 3.2, "Exon 3 (Unmethylated)", fontsize=9.0, fontweight="bold", color="#2c3e50", ha="center")

# m6A Peak star on Exon 2
ax1.plot(55, 6.2, marker="*", markersize=15, color="#e74c3c", zorder=5)
ax1.text(55, 6.8, "m6A Peak Site", fontsize=9.0, fontweight="bold", color="#e74c3c", ha="center")

# Input Reads (Uniform Green Coverage across all exons)
x_in = np.array([14, 25, 35, 45, 55, 65, 76, 85, 90])
y_in = np.array([7.8, 8.1, 7.9, 8.2, 8.4, 8.2, 8.0, 7.8, 7.6])
ax1.fill_between(x_in, 7.6, y_in, color="#27ae60", alpha=0.35, zorder=3)
ax1.plot(x_in, y_in, color="#27ae60", lw=2.2, zorder=4)
ax1.text(12, 8.2, "Input Reads:", fontsize=9.5, fontweight="bold", color="#27ae60", ha="right")
ax1.text(35, 8.9, "[UNBIASED] Uniform coverage across all exons -> Detects variants everywhere", 
         fontsize=8.5, fontweight="bold", color="#1e8449")

# IP Reads (Localized Red Peak around Exon 2 only)
x_ip = np.array([14, 35, 42, 50, 55, 60, 68, 76, 90])
y_ip = np.array([1.5, 1.5, 1.6, 2.7, 3.1, 2.6, 1.6, 1.5, 1.5])
ax1.fill_between(x_ip, 1.5, y_ip, color="#e74c3c", alpha=0.35, zorder=3)
ax1.plot(x_ip, y_ip, color="#e74c3c", lw=2.2, zorder=4)
ax1.text(12, 2.0, "IP Reads:", fontsize=9.5, fontweight="bold", color="#e74c3c", ha="right")
ax1.text(55, 0.4, "[BIASED COVERAGE] Zero reads on Exon 1 & Exon 3 -> Misses all unmethylated variants!", 
         fontsize=9.0, fontweight="bold", color="#c0392b", ha="center")

# -----------------------------------------------------------------------------
# 2. CIRCULAR SELECTION BIAS (Bolded Main Workflow Steps)
# -----------------------------------------------------------------------------
ax2 = axes[1]
ax2.set_facecolor("#ffffff")
ax2.axis("off")
ax2.set_xlim(0, 100)
ax2.set_ylim(0, 11)

ax2.text(0, 10.2, "2. Circular Selection Bias vs. Independent Input Guardrail", 
         fontsize=12.0, fontweight="bold", color="#1a252f")

# =============================================================================
# LEFT SECTION: FLAWED (Calling SNPs from IP Reads)
# =============================================================================
box_left = patches.FancyBboxPatch((1, 0.5), 47, 9.2, boxstyle="round,pad=0.2,rounding_size=0.3",
                                  linewidth=1.8, edgecolor="#e74c3c", facecolor="#fdf2e9")
ax2.add_patch(box_left)
ax2.text(24.5, 9.1, "FLAWED: Variant Calling from IP Reads", fontsize=10.5, fontweight="bold", color="#c0392b", ha="center")

# Sub-panel 1A: Workflow Steps (Bolded Main Steps)
sub_1a = patches.FancyBboxPatch((2.5, 3.2), 21, 5.2, boxstyle="round,pad=0.1,rounding_size=0.2",
                                 linewidth=1.0, edgecolor="#d98880", facecolor="#ffffff")
ax2.add_patch(sub_1a)
ax2.text(13, 7.9, "Workflow Steps", fontsize=9.0, fontweight="bold", color="#78281f", ha="center")

ax2.text(3.5, 7.1, "Step 1. Read Alignment", fontsize=8.5, fontweight="bold", color="#1a252f")
ax2.text(3.5, 6.5, "   IP Antibody Capture", fontsize=7.8, color="#566573")

ax2.text(3.5, 5.6, "Step 2. Variant Calling", fontsize=8.5, fontweight="bold", color="#1a252f")
ax2.text(3.5, 5.0, "   Call from IP BAM", fontsize=7.8, color="#566573")

ax2.text(3.5, 4.1, "Step 3. ASM Testing", fontsize=8.5, fontweight="bold", color="#1a252f")
ax2.text(3.5, 3.5, "   Evaluate IP Imbalance", fontsize=7.8, color="#566573")

# Sub-panel 1B: Concrete Read Example
sub_1b = patches.FancyBboxPatch((25, 3.2), 21.5, 5.2, boxstyle="round,pad=0.1,rounding_size=0.2",
                                 linewidth=1.0, edgecolor="#d98880", facecolor="#fadbd8")
ax2.add_patch(sub_1b)
ax2.text(35.75, 7.9, "Read Example", fontsize=9.0, fontweight="bold", color="#78281f", ha="center")

ax2.text(26, 6.8, "IP Counts:\n90 ALT (G) | 0 REF (A)", fontsize=8.0, fontweight="bold", color="#922b21")
ax2.text(26, 5.3, "Genotype Result:\nMiscalled 1/1 (G/G)", fontsize=8.0, fontweight="bold", color="#c0392b")
ax2.text(26, 3.8, "Test Result:\nCircular Bias / Fail", fontsize=8.0, fontweight="bold", color="#78281f")

# Left Bottom Summary Badge
rect_left_res = patches.FancyBboxPatch((2.5, 0.9), 44, 1.8, boxstyle="round,pad=0.1,rounding_size=0.15",
                                        linewidth=1.0, edgecolor="#c0392b", facecolor="#f5b7b1")
ax2.add_patch(rect_left_res)
ax2.text(24.5, 1.8, "CIRCULAR SELECTION BIAS", fontsize=9.0, fontweight="bold", color="#78281f", ha="center")
ax2.text(24.5, 1.2, "Site discovered BECAUSE it was methylated (False Positives)", fontsize=7.8, color="#78281f", ha="center")


# =============================================================================
# RIGHT SECTION: CORRECT (Input-First Guardrail Workflow)
# =============================================================================
box_right = patches.FancyBboxPatch((51, 0.5), 47, 9.2, boxstyle="round,pad=0.2,rounding_size=0.3",
                                   linewidth=1.8, edgecolor="#27ae60", facecolor="#eafaf1")
ax2.add_patch(box_right)
ax2.text(74.5, 9.1, "CORRECT: Input-First Guardrail Workflow", fontsize=10.5, fontweight="bold", color="#1e8449", ha="center")

# Sub-panel 2A: Workflow Steps (Bolded Main Steps)
sub_2a = patches.FancyBboxPatch((52.5, 3.2), 21, 5.2, boxstyle="round,pad=0.1,rounding_size=0.2",
                                 linewidth=1.0, edgecolor="#82e0aa", facecolor="#ffffff")
ax2.add_patch(sub_2a)
ax2.text(63, 7.9, "Workflow Steps", fontsize=9.0, fontweight="bold", color="#145a32", ha="center")

ax2.text(53.5, 7.1, "Step 1. Read Alignment", fontsize=8.5, fontweight="bold", color="#1a252f")
ax2.text(53.5, 6.5, "   Input Baseline", fontsize=7.8, color="#566573")

ax2.text(53.5, 5.6, "Step 2. Variant Calling", fontsize=8.5, fontweight="bold", color="#1a252f")
ax2.text(53.5, 5.0, "   Call from Input BAM", fontsize=7.8, color="#566573")

ax2.text(53.5, 4.1, "Step 3. ASM Testing", fontsize=8.5, fontweight="bold", color="#1a252f")
ax2.text(53.5, 3.5, "   Quantify IP at 0/1", fontsize=7.8, color="#566573")

# Sub-panel 2B: Concrete Read Example
sub_2b = patches.FancyBboxPatch((75, 3.2), 21.5, 5.2, boxstyle="round,pad=0.1,rounding_size=0.2",
                                 linewidth=1.0, edgecolor="#82e0aa", facecolor="#d4efdf")
ax2.add_patch(sub_2b)
ax2.text(85.75, 7.9, "Read Example", fontsize=9.0, fontweight="bold", color="#145a32", ha="center")

ax2.text(76, 6.8, "Input Baseline:\n50 REF (A) | 50 ALT (G)", fontsize=8.0, fontweight="bold", color="#145a32")
ax2.text(76, 5.3, "Genotype Result:\nTrue Heterozygous 0/1", fontsize=8.0, fontweight="bold", color="#1e8449")
ax2.text(76, 3.8, "IP Read Counts:\n10 REF (A) | 90 ALT (G)", fontsize=8.0, fontweight="bold", color="#145a32")

# Right Bottom Summary Badge
rect_right_res = patches.FancyBboxPatch((52.5, 0.9), 44, 1.8, boxstyle="round,pad=0.1,rounding_size=0.15",
                                         linewidth=1.0, edgecolor="#27ae60", facecolor="#abebc6")
ax2.add_patch(rect_right_res)
ax2.text(74.5, 1.8, "UNBIASED ALLELE-SPECIFIC TEST", fontsize=9.0, fontweight="bold", color="#145a32", ha="center")
ax2.text(74.5, 1.2, "Variant discovery is 100% independent of m6A methylation", fontsize=7.8, color="#145a32", ha="center")

# Save PNG and PDF
plt.savefig(out_png, dpi=300, bbox_inches='tight')
plt.savefig(out_pdf, bbox_inches='tight')
plt.close()

print(f"[SUCCESS] Updated Input-only rationale diagram with bold main steps saved to:\n  PNG: {out_png}\n  PDF: {out_pdf}")
