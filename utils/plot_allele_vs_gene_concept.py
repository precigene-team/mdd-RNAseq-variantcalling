#!/usr/bin/env python3
# [[utils/plot_allele_vs_gene_concept.py]]
# Goal: Visualize Gene vs Allele concept with dual arrows pointing to BOTH Maternal and Paternal chromosomes.

import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches

# Enable LaTeX mathtext
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

out_png = os.path.join(output_dir, "allele_vs_gene_concept.png")
out_pdf = os.path.join(output_dir, "allele_vs_gene_concept.pdf")

# Color Palette:
MATERNAL_COLOR = "#e91e63" # Pink
PATERNAL_COLOR = "#2980b9" # Blue
GENE_COLOR     = "#d35400" # Orange
VARIANT_COLOR  = "#8e44ad" # Purple

# --- 2. Figure Setup (2x2 Grid) ---
fig = plt.figure(figsize=(15.5, 8.2), dpi=300, facecolor="#ffffff")
plt.rcParams['font.sans-serif'] = 'DejaVu Sans'

gs = fig.add_gridspec(2, 2, width_ratios=[1.0, 1.3], height_ratios=[1.0, 1.0], 
                       wspace=0.18, hspace=0.25,
                       left=0.04, right=0.96, bottom=0.05, top=0.93)

ax_tl = fig.add_subplot(gs[0, 0]) # Top Left: Case 1 Chromosomes
ax_bl = fig.add_subplot(gs[1, 0]) # Bottom Left: Case 2 Chromosomes
ax_tr = fig.add_subplot(gs[0, 1]) # Top Right: Case 1 Magnified Track
ax_br = fig.add_subplot(gs[1, 1]) # Bottom Right: Case 2 Magnified Track

# Figure Main Title
fig.suptitle("Gene & Allele Concept: Homozygous vs. Heterozygous ASm6A Loci", 
             fontsize=14, fontweight="bold", color="#1a252f", y=0.97)

# Helper function to draw Upright X-shaped Chromosome pair
def draw_upright_chromosomes(ax, center_x, center_y, title_lbl, scale=0.85, show_variant_site=True):
    arm_len = 20.0 * scale
    w = 3.8 * scale
    
    # Left Homolog (Maternal - Pink)
    cx1 = center_x - 14.0 * scale
    cy1 = center_y
    # Sister Chromatid 1 (Left Arm Pair - Faded Alpha = 0.35)
    ax.plot([cx1, cx1 - arm_len*0.45], [cy1, cy1 + arm_len], color=MATERNAL_COLOR, lw=w*2, solid_capstyle='round', alpha=0.35, zorder=2)
    ax.plot([cx1, cx1 - arm_len*0.45], [cy1, cy1 - arm_len], color=MATERNAL_COLOR, lw=w*2, solid_capstyle='round', alpha=0.35, zorder=2)
    
    # Sister Chromatid 2 (Right Arm Pair - Primary Solid Alpha = 1.0)
    ax.plot([cx1, cx1 + arm_len*0.45], [cy1, cy1 + arm_len], color=MATERNAL_COLOR, lw=w*2, solid_capstyle='round', zorder=3)
    ax.plot([cx1, cx1 + arm_len*0.45], [cy1, cy1 - arm_len], color=MATERNAL_COLOR, lw=w*2, solid_capstyle='round', zorder=3)
    
    # Paint Gene Region RED on BOTH Chromatid Arms (Left Faded, Right Solid)
    g1_tl_s = (cx1 - arm_len*0.18, cy1 + arm_len*0.40)
    g1_tl_e = (cx1 - arm_len*0.38, cy1 + arm_len*0.84)
    ax.plot([g1_tl_s[0], g1_tl_e[0]], [g1_tl_s[1], g1_tl_e[1]], color="#e74c3c", lw=w*2.5, solid_capstyle='round', alpha=0.35, zorder=4)

    g1_tr_s = (cx1 + arm_len*0.18, cy1 + arm_len*0.40)
    g1_tr_e = (cx1 + arm_len*0.38, cy1 + arm_len*0.84)
    ax.plot([g1_tr_s[0], g1_tr_e[0]], [g1_tr_s[1], g1_tr_e[1]], color="#e74c3c", lw=w*2.5, solid_capstyle='round', zorder=4)
    
    snp1 = None
    if show_variant_site:
        snp1 = ((g1_tr_s[0] + g1_tr_e[0])/2.0, (g1_tr_s[1] + g1_tr_e[1])/2.0)
        ax.plot(snp1[0], snp1[1], marker="o", markersize=9*scale, color=VARIANT_COLOR, markeredgecolor="#ffffff", markeredgewidth=1.2, zorder=6)
    
    ax.plot(cx1, cy1, marker="o", markersize=9*scale, color="#c2185b", zorder=5)

    # Right Homolog (Paternal - Blue)
    cx2 = center_x + 14.0 * scale
    cy2 = center_y
    # Sister Chromatid 1 (Left Arm Pair - Faded Alpha = 0.35)
    ax.plot([cx2, cx2 - arm_len*0.45], [cy2, cy2 + arm_len], color=PATERNAL_COLOR, lw=w*2, solid_capstyle='round', alpha=0.35, zorder=2)
    ax.plot([cx2, cx2 - arm_len*0.45], [cy2, cy2 - arm_len], color=PATERNAL_COLOR, lw=w*2, solid_capstyle='round', alpha=0.35, zorder=2)

    # Sister Chromatid 2 (Right Arm Pair - Primary Solid Alpha = 1.0)
    ax.plot([cx2, cx2 + arm_len*0.45], [cy2, cy2 + arm_len], color=PATERNAL_COLOR, lw=w*2, solid_capstyle='round', zorder=3)
    ax.plot([cx2, cx2 + arm_len*0.45], [cy2, cy2 - arm_len], color=PATERNAL_COLOR, lw=w*2, solid_capstyle='round', zorder=3)
    
    # Paint Gene Region RED on BOTH Chromatid Arms (Left Faded, Right Solid)
    g2_tl_s = (cx2 - arm_len*0.18, cy2 + arm_len*0.40)
    g2_tl_e = (cx2 - arm_len*0.38, cy2 + arm_len*0.84)
    ax.plot([g2_tl_s[0], g2_tl_e[0]], [g2_tl_s[1], g2_tl_e[1]], color="#e74c3c", lw=w*2.5, solid_capstyle='round', alpha=0.35, zorder=4)

    g2_tr_s = (cx2 + arm_len*0.18, cy2 + arm_len*0.40)
    g2_tr_e = (cx2 + arm_len*0.38, cy2 + arm_len*0.84)
    ax.plot([g2_tr_s[0], g2_tr_e[0]], [g2_tr_s[1], g2_tr_e[1]], color="#e74c3c", lw=w*2.5, solid_capstyle='round', zorder=4)
    
    snp2 = None
    if show_variant_site:
        snp2 = ((g2_tr_s[0] + g2_tr_e[0])/2.0, (g2_tr_s[1] + g2_tr_e[1])/2.0)
        ax.plot(snp2[0], snp2[1], marker="o", markersize=9*scale, color=VARIANT_COLOR, markeredgecolor="#ffffff", markeredgewidth=1.2, zorder=6)

    ax.plot(cx2, cy2, marker="o", markersize=9*scale, color="#1b4f72", zorder=5)

    # Subplot Header & Labels
    ax.text(50.0, 92.0, title_lbl, fontsize=11, fontweight="bold", color="#1a252f", ha="center")
    ax.text(cx1, cy1 - arm_len - 6.0, "Maternal", fontsize=8.5, fontweight="bold", color=MATERNAL_COLOR, ha="center")
    ax.text(cx2, cy2 - arm_len - 6.0, "Paternal", fontsize=8.5, fontweight="bold", color=PATERNAL_COLOR, ha="center")

    return g1_tr_e, g2_tr_e, snp1, snp2

# Helper function to draw sequence track
def draw_custom_sequence_track_with_sister(ax, y_center, base_seq, is_maternal=True, gene_indices=[2, 3, 4], variant_idx=None):
    x_start = 28.0
    box_w = 5.6
    box_gap = 1.0
    box_h = 7.0   # Height of primary box
    faded_h = 3.2  # Height of faded sister box
    
    reg_color = MATERNAL_COLOR if is_maternal else PATERNAL_COLOR
    reg_bg    = "#fce4ec" if is_maternal else "#ebf5fb"
    reg_text  = "#c2185b" if is_maternal else "#1b4f72"
    
    # --- 1. FADED SISTER CHROMATID TRACK ---
    y_faded_bottom = y_center + 5.2
    
    for k in range(2, 0, -1):
        x_cont = x_start - k * (box_w + box_gap)
        b_box = patches.FancyBboxPatch((x_cont, y_faded_bottom), box_w, faded_h, 
                                       boxstyle="round,pad=0.02",
                                       facecolor=reg_bg, edgecolor=reg_color, lw=0.9, linestyle="--", alpha=0.25, zorder=2)
        ax.add_patch(b_box)
        
    for i, letter in enumerate(base_seq):
        x_pos = x_start + i * (box_w + box_gap)
        if variant_idx is not None and i == variant_idx:
            bg_c = VARIANT_COLOR if letter == "G" else "#fdebd0"
            border_c = "#5b2c6f" if letter == "G" else VARIANT_COLOR
            b_box = patches.FancyBboxPatch((x_pos, y_faded_bottom), box_w, faded_h, 
                                           boxstyle="round,pad=0.02",
                                           facecolor=bg_c, edgecolor=border_c, lw=1.2, alpha=0.35, zorder=2)
            ax.add_patch(b_box)
        elif i in gene_indices:
            b_box = patches.FancyBboxPatch((x_pos, y_faded_bottom), box_w, faded_h, 
                                           boxstyle="round,pad=0.02",
                                           facecolor="#fdebd0", edgecolor=GENE_COLOR, lw=1.0, alpha=0.35, zorder=2)
            ax.add_patch(b_box)
        else:
            b_box = patches.FancyBboxPatch((x_pos, y_faded_bottom), box_w, faded_h, 
                                           boxstyle="round,pad=0.02",
                                           facecolor=reg_bg, edgecolor=reg_color, lw=1.0, alpha=0.35, zorder=2)
            ax.add_patch(b_box)
            
    x_end_main = x_start + len(base_seq) * (box_w + box_gap)
    for k in range(2):
        x_cont = x_end_main + k * (box_w + box_gap)
        b_box = patches.FancyBboxPatch((x_cont, y_faded_bottom), box_w, faded_h, 
                                       boxstyle="round,pad=0.02",
                                       facecolor=reg_bg, edgecolor=reg_color, lw=0.9, linestyle="--", alpha=0.25, zorder=2)
        ax.add_patch(b_box)

    # --- 2. PRIMARY SOLID CHROMATID TRACK ---
    y_primary_bottom = y_center - 3.5
    
    for k in range(2, 0, -1):
        x_cont = x_start - k * (box_w + box_gap)
        b_box = patches.FancyBboxPatch((x_cont, y_primary_bottom), box_w, box_h, 
                                       boxstyle="round,pad=0.05",
                                       facecolor=reg_bg, edgecolor=reg_color, lw=1.3, linestyle="--", zorder=4)
        ax.add_patch(b_box)

    for i, letter in enumerate(base_seq):
        x_pos = x_start + i * (box_w + box_gap)
        
        if variant_idx is not None and i == variant_idx:
            if letter == "G":
                # ALT Variant Allele G: Solid PURPLE Background (#8e44ad) with Bold WHITE Text (#ffffff)!
                b_box = patches.FancyBboxPatch((x_pos, y_primary_bottom), box_w, box_h, 
                                               boxstyle="round,pad=0.05",
                                               facecolor=VARIANT_COLOR, edgecolor="#5b2c6f", lw=2.2, zorder=5)
                ax.add_patch(b_box)
                ax.text(x_pos + box_w/2.0, y_primary_bottom + box_h/2.0, letter, fontsize=10, fontweight="bold", color="#ffffff", ha="center", va="center", zorder=6)
            else:
                # REF Allele A: Orange Background (#fdebd0) with Purple Border (#8e44ad)
                b_box = patches.FancyBboxPatch((x_pos, y_primary_bottom), box_w, box_h, 
                                               boxstyle="round,pad=0.05",
                                               facecolor="#fdebd0", edgecolor=VARIANT_COLOR, lw=2.2, zorder=5)
                ax.add_patch(b_box)
                ax.text(x_pos + box_w/2.0, y_primary_bottom + box_h/2.0, letter, fontsize=9.5, fontweight="bold", color=GENE_COLOR, ha="center", va="center", zorder=6)
            
        elif i in gene_indices:
            b_box = patches.FancyBboxPatch((x_pos, y_primary_bottom), box_w, box_h, 
                                           boxstyle="round,pad=0.05",
                                           facecolor="#fdebd0", edgecolor=GENE_COLOR, lw=1.3, zorder=4)
            ax.add_patch(b_box)
            ax.text(x_pos + box_w/2.0, y_primary_bottom + box_h/2.0, letter, fontsize=9.5, fontweight="bold", color=GENE_COLOR, ha="center", va="center", zorder=5)
            
        else:
            b_box = patches.FancyBboxPatch((x_pos, y_primary_bottom), box_w, box_h, 
                                           boxstyle="round,pad=0.05",
                                           facecolor=reg_bg, edgecolor=reg_color, lw=1.3, zorder=4)
            ax.add_patch(b_box)
            ax.text(x_pos + box_w/2.0, y_primary_bottom + box_h/2.0, letter, fontsize=9.5, fontweight="bold", color=reg_text, ha="center", va="center", zorder=5)

    for k in range(2):
        x_cont = x_end_main + k * (box_w + box_gap)
        b_box = patches.FancyBboxPatch((x_cont, y_primary_bottom), box_w, box_h, 
                                       boxstyle="round,pad=0.05",
                                       facecolor=reg_bg, edgecolor=reg_color, lw=1.3, linestyle="--", zorder=4)
        ax.add_patch(b_box)

# ==============================================================================
# 1. TOP-LEFT PANEL: CASE 1 CHROMOSOMES (GENE LEVEL, NO VARIANT SITE)
# ==============================================================================
ax_tl.set_facecolor('#fafbfc')
ax_tl.set_xlim(0, 100)
ax_tl.set_ylim(0, 100)
ax_tl.axis('off')

b_tl = patches.FancyBboxPatch((0.5, 0.5), 99.0, 99.0, boxstyle="round,pad=0.2",
                              facecolor="#fafbfc", edgecolor="#d0d7de", lw=1.3)
ax_tl.add_patch(b_tl)

g1_mat_pt, g1_pat_pt, _, _ = draw_upright_chromosomes(ax_tl, center_x=50.0, center_y=45.0, 
                                                       title_lbl="Left Panel: Case 1 Chromosomes (Gene Level)", 
                                                       scale=0.85, show_variant_site=False)

# Gene Callout Header centered at X=50, pointing to BOTH Maternal & Paternal RED gene segments
ax_tl.text(50.0, 79.0, "Gene", fontsize=9.5, fontweight="bold", color="#c0392b", ha="center")
ax_tl.annotate("", xy=(g1_mat_pt[0], g1_mat_pt[1]), xytext=(45.0, 77.0),
               arrowprops=dict(arrowstyle="->", color="#e74c3c", lw=1.6))
ax_tl.annotate("", xy=(g1_pat_pt[0], g1_pat_pt[1]), xytext=(55.0, 77.0),
               arrowprops=dict(arrowstyle="->", color="#e74c3c", lw=1.6))

# ==============================================================================
# 2. BOTTOM-LEFT PANEL: CASE 2 CHROMOSOMES (VARIANT SITES ON BOTH CHROMOSOMES)
# ==============================================================================
ax_bl.set_facecolor('#fafbfc')
ax_bl.set_xlim(0, 100)
ax_bl.set_ylim(0, 100)
ax_bl.axis('off')

b_bl = patches.FancyBboxPatch((0.5, 0.5), 99.0, 99.0, boxstyle="round,pad=0.2",
                              facecolor="#fef9e7", edgecolor="#27ae60", lw=1.8)
ax_bl.add_patch(b_bl)

g2_mat_pt, g2_pat_pt, s2_maternal, s2_paternal = draw_upright_chromosomes(ax_bl, center_x=50.0, center_y=45.0, 
                                                                          title_lbl="Left Panel: Case 2 Chromosomes (Variant Site Level)", 
                                                                          scale=0.85, show_variant_site=True)

# Variant Site Callout centered at X=50, pointing to BOTH Maternal & Paternal PURPLE SNP dots
ax_bl.text(50.0, 79.0, "Variant Site", fontsize=9.5, fontweight="bold", color="#6c3483", ha="center")
ax_bl.annotate("", xy=(s2_maternal[0], s2_maternal[1]), xytext=(45.0, 77.0),
               arrowprops=dict(arrowstyle="->", color=VARIANT_COLOR, lw=1.6))
ax_bl.annotate("", xy=(s2_paternal[0], s2_paternal[1]), xytext=(55.0, 77.0),
               arrowprops=dict(arrowstyle="->", color=VARIANT_COLOR, lw=1.6))

# ==============================================================================
# 3. TOP-RIGHT PANEL: CASE 1 MAGNIFIED TRACK (HOMOZYGOUS A / A)
# ==============================================================================
ax_tr.set_facecolor('#ffffff')
ax_tr.set_xlim(0, 100)
ax_tr.set_ylim(0, 100)
ax_tr.axis('off')

b_tr = patches.FancyBboxPatch((0.5, 0.5), 99.0, 99.0, boxstyle="round,pad=0.2",
                              facecolor="#ffffff", edgecolor="#bdc3c7", lw=1.4)
ax_tr.add_patch(b_tr)

ax_tr.text(50.0, 88.0, "Right Panel: Case 1 — Homozygous (0/0 or 1/1)", 
           fontsize=10.5, fontweight="bold", color="#7f8c8d", ha="center")

# Track 1: Maternal Copy
ax_tr.text(4.0, 60.0, "Maternal:", fontsize=8.5, fontweight="bold", color=MATERNAL_COLOR)
draw_custom_sequence_track_with_sister(ax_tr, y_center=60.0, base_seq=["C", "T", "G", "A", "C", "A", "T"], is_maternal=True, gene_indices=[2, 3, 4], variant_idx=None)

# Track 2: Paternal Copy
ax_tr.text(4.0, 36.0, "Paternal:", fontsize=8.5, fontweight="bold", color=PATERNAL_COLOR)
draw_custom_sequence_track_with_sister(ax_tr, y_center=36.0, base_seq=["C", "T", "G", "A", "C", "A", "T"], is_maternal=False, gene_indices=[2, 3, 4], variant_idx=None)

# Result Box
res_tr = patches.FancyBboxPatch((5.0, 4.0), 90.0, 21.0, boxstyle="round,pad=0.1",
                                facecolor="#f4f6f7", edgecolor="#bdc3c7", lw=1.2)
ax_tr.add_patch(res_tr)
ax_tr.text(50.0, 16.0, "Uninformative Locus: Homozygous (0/0: A/A  or  1/1: G/G)", fontsize=9, fontweight="bold", color="#7f8c8d", ha="center")
ax_tr.text(50.0, 8.0, "Both copies carry identical sequence (e.g. A/A or G/G).\nCannot distinguish maternal vs. paternal RNA reads.", 
           fontsize=8, color="#566573", ha="center", va="center")

# ==============================================================================
# 4. BOTTOM-RIGHT PANEL: CASE 2 MAGNIFIED TRACK (HETEROZYGOUS A vs G TARGET!)
# ==============================================================================
ax_br.set_facecolor('#ffffff')
ax_br.set_xlim(0, 100)
ax_br.set_ylim(0, 100)
ax_br.axis('off')

# Prominent Green Highlight Outer Box
b_br = patches.FancyBboxPatch((0.5, 0.5), 99.0, 99.0, boxstyle="round,pad=0.2",
                              facecolor="#fef9e7", edgecolor="#27ae60", lw=2.2)
ax_br.add_patch(b_br)

# Target Header Badge
target_badge = patches.FancyBboxPatch((8.0, 84.0), 84.0, 12.0, boxstyle="round,pad=0.1",
                                        facecolor="#27ae60", edgecolor="none")
ax_br.add_patch(target_badge)
ax_br.text(50.0, 90.0, "Right Panel: Case 2 — Heterozygous (0/1) TARGET!", 
           fontsize=10.5, fontweight="bold", color="#ffffff", ha="center", va="center")

# Track 1: Maternal REF Allele A
ax_br.text(4.0, 60.0, "Maternal:", fontsize=8.5, fontweight="bold", color=MATERNAL_COLOR)
draw_custom_sequence_track_with_sister(ax_br, y_center=60.0, base_seq=["C", "T", "G", "A", "C", "A", "T"], is_maternal=True, gene_indices=[2, 3, 4], variant_idx=3)
ax_br.text(88.0, 60.0, "m6A Marked", fontsize=7.0, fontweight="bold", color="#1e8449", va="center")

# Track 2: Paternal ALT Allele G (SOLID PURPLE BACKGROUND #8e44ad)
ax_br.text(4.0, 36.0, "Paternal:", fontsize=8.5, fontweight="bold", color=PATERNAL_COLOR)
draw_custom_sequence_track_with_sister(ax_br, y_center=36.0, base_seq=["C", "T", "G", "G", "C", "A", "T"], is_maternal=False, gene_indices=[2, 3, 4], variant_idx=3)
ax_br.text(88.0, 36.0, "Unmethylated", fontsize=7.0, fontweight="bold", color="#1b4f72", va="center")

# Result Box - Highlighted Target
res_br = patches.FancyBboxPatch((5.0, 4.0), 90.0, 21.0, boxstyle="round,pad=0.1",
                                facecolor="#ffffff", edgecolor="#27ae60", lw=1.6)
ax_br.add_patch(res_br)
ax_br.text(50.0, 16.0, "Target Locus: Heterozygous Variant (0/1: A / G)", fontsize=9, fontweight="bold", color="#27ae60", ha="center")
ax_br.text(50.0, 8.0, "Heterozygous genotype (0/1: A vs G) enables allele-specific read separation.\nRequired to detect differential methylation (ASm6A) between maternal & paternal copies!", 
           fontsize=8, fontweight="bold", color="#1f4e79", ha="center", va="center")

# Save PNG and PDF
plt.savefig(out_png, dpi=300, bbox_inches='tight')
plt.savefig(out_pdf, bbox_inches='tight')
plt.close()

print(f"[SUCCESS] Gene vs Allele concept figure updated (Dual arrows pointing to BOTH Maternal and Paternal chromosomes; clean 'Variant Site' label) saved to:\n  PNG: {out_png}\n  PDF: {out_pdf}")
