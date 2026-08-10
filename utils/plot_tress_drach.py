#!/usr/bin/env python3
# [[utils/plot_tress_drach.py]]
# Goal: Visualize TRESS peak calling and DRACH motif sequence logo with vertically stretched letters on probability axis for Slide 7.

import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.transforms as mtransforms
from matplotlib.textpath import TextPath

# --- 1. Dynamic Pathing ---
def get_project_root():
    curr = os.getcwd()
    while curr != "/" and not os.path.exists(os.path.join(curr, "run")):
        curr = os.path.dirname(curr)
    return curr

project_root = get_project_root()
output_dir = os.path.join(project_root, "results", "utils", "figures")
os.makedirs(output_dir, exist_ok=True)

out_png = os.path.join(output_dir, "tress_drach_motif.png")
out_pdf = os.path.join(output_dir, "tress_drach_motif.pdf")

# Helper function: Draw vertically stretched nucleotide letter patch
def add_stretched_letter(ax, char, x_center, y_bottom, height, width, color, zorder=3):
    if height <= 0.001:
        return
    tp = TextPath((0, 0), char, size=1.0, prop=dict(weight='bold', family='sans-serif'))
    bbox = tp.get_extents()
    
    w_orig = bbox.width
    h_orig = bbox.height
    if w_orig == 0 or h_orig == 0:
        return
        
    transform = (
        mtransforms.Affine2D()
        .translate(-bbox.x0, -bbox.y0)
        .scale(width / w_orig, height / h_orig)
        .translate(x_center - width / 2.0, y_bottom)
        + ax.transData
    )
    
    patch = patches.PathPatch(tp, transform=transform, facecolor=color, edgecolor='none', zorder=zorder)
    ax.add_patch(patch)

# --- 2. Figure Setup ---
fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(13.5, 5.8), dpi=300, facecolor="#ffffff")
plt.rcParams['font.sans-serif'] = 'DejaVu Sans'

plt.subplots_adjust(left=0.06, right=0.96, bottom=0.15, top=0.88, wspace=0.26)

# ==============================================================================
# PANEL A: Peak Calling (IP vs Input)
# ==============================================================================
ax_a.set_facecolor('#ffffff')

x = np.linspace(0, 100, 500)

input_cov = 15 + 10 * np.exp(-((x - 75) / 15)**2)
ip_cov = 18 + 80 * np.exp(-((x - 40) / 10)**2) + 12 * np.exp(-((x - 75) / 15)**2)

# IP Track (Navy)
ax_a.plot(x, ip_cov, color='#1f4e79', lw=3.0, label='IP (m6A RNA)')
ax_a.fill_between(x, 0, ip_cov, color='#1f4e79', alpha=0.15)

# Input Track (Teal Dashed)
ax_a.plot(x, input_cov, color='#16a085', lw=2.5, linestyle='--', label='Input (Total RNA)')
ax_a.fill_between(x, 0, input_cov, color='#16a085', alpha=0.10)

# Highlight Peak Region
ax_a.axvspan(30, 50, color='#f39c12', alpha=0.22, zorder=0)

ax_a.text(40, 92, "IP-Enriched Peak\n(Normalized by Input)", 
          fontsize=10, fontweight='bold', color='#d35400', ha='center', va='bottom',
          bbox=dict(boxstyle="round,pad=0.4", fc="#fef9e7", ec="#f39c12", lw=1.2))

# Minimal Exon Line
ax_a.plot([5, 95], [-12, -12], color='#7f8c8d', lw=2)
ax_a.add_patch(patches.Rectangle((20, -17), 60, 10, facecolor="#34495e", edgecolor="none"))
ax_a.text(50, -28, "Exon / Transcript Region", fontsize=9, color="#7f8c8d", ha="center")

ax_a.set_ylim(-35, 115)
ax_a.set_xlim(0, 100)
ax_a.set_title("1. Peak Calling (IP / Input Normalization)", fontsize=12.5, fontweight='bold', color='#1a252f', pad=12)
ax_a.set_ylabel("Read Depth", fontsize=10.5, fontweight='bold', color='#2c3e50')
ax_a.legend(loc='upper right', frameon=True, facecolor='#ffffff', edgecolor='#d0d7de', fontsize=9.5)

ax_a.spines['top'].set_visible(False)
ax_a.spines['right'].set_visible(False)
ax_a.set_xticks([])

# ==============================================================================
# PANEL B: DRACH Motif with Stretched Letter Sequence Logo
# ==============================================================================
ax_b.set_facecolor('#ffffff')

# Nucleotide appearance probabilities at each position (Pos 1 to 5)
positions = [
    {"x": 1, "code": "D", "desc": "A / G / U\n(non-C)", "bases": [("A", 0.45, "#27ae60"), ("G", 0.35, "#f39c12"), ("U", 0.20, "#e74c3c")]},
    {"x": 2, "code": "R", "desc": "A / G\n(Purine)", "bases": [("A", 0.52, "#27ae60"), ("G", 0.48, "#f39c12")]},
    {"x": 3, "code": "A*", "desc": "m6A Site\n(Central A)", "bases": [("A", 1.00, "#d35400")]},
    {"x": 4, "code": "C", "desc": "Cytosine", "bases": [("C", 1.00, "#2980b9")]},
    {"x": 5, "code": "H", "desc": "A / C / U\n(non-G)", "bases": [("U", 0.40, "#e74c3c"), ("A", 0.35, "#27ae60"), ("C", 0.25, "#2980b9")]}
]

# Draw vertically stretched sequence logo
col_width = 0.58
for p in positions:
    x_pos = p["x"]
    y_bottom = 0.0
    
    # Highlight position 3 (central m6A Adenine)
    if p["code"] == "A*":
        bg = patches.Rectangle((x_pos - 0.36, 0.0), 0.72, 1.0, 
                               facecolor="#fef9e7", edgecolor="#f39c12", lw=1.5, zorder=1)
        ax_b.add_patch(bg)
        
        # -CH3 Stamp Badge above Y=1.0
        ch3_badge = patches.FancyBboxPatch((x_pos - 0.36, 1.04), 0.72, 0.12, boxstyle="round,pad=0.03",
                                           facecolor="#e67e22", edgecolor="#d35400", lw=1.2, zorder=4)
        ax_b.add_patch(ch3_badge)
        ax_b.text(x_pos, 1.10, "-CH3", fontsize=9.5, fontweight='bold', color="#ffffff", 
                  ha='center', va='center', zorder=5)

    for char, prob, col in p["bases"]:
        add_stretched_letter(ax_b, char, x_pos, y_bottom, prob, col_width, col, zorder=3)
        y_bottom += prob

ax_b.set_ylim(0, 1.28)
ax_b.set_xlim(0.3, 5.7)
ax_b.set_yticks([0.0, 0.25, 0.50, 0.75, 1.00])
ax_b.set_yticklabels(['0%', '25%', '50%', '75%', '100%'], fontsize=9.5, color='#2c3e50')
ax_b.set_ylabel("Appearance Chance (Probability)", fontsize=10.5, fontweight='bold', color='#2c3e50')

# Set X-axis tick labels
x_ticks = [p["x"] for p in positions]
x_labels = [f"{p['code']}\n({p['desc'].splitlines()[0]})" for p in positions]
ax_b.set_xticks(x_ticks)
ax_b.set_xticklabels(x_labels, fontsize=10, fontweight='bold', color='#2c3e50')

# Highlight position 3 label
ax_b.get_xticklabels()[2].set_color('#d35400')

ax_b.set_title("2. DRACH Consensus Motif (Sequence Logo)", fontsize=12.5, fontweight='bold', color='#1a252f', pad=12)

# Clean spines & grid
ax_b.spines['top'].set_visible(False)
ax_b.spines['right'].set_visible(False)
ax_b.spines['left'].set_color('#d0d7de')
ax_b.spines['bottom'].set_color('#d0d7de')
ax_b.grid(axis='y', linestyle='--', alpha=0.3, color='#bdc3c7')

# Bottom Callout Box with text wrapping
callout_text = (
    "Genome-wide Candidate m6A Sites\n"
    "Scanning TRESS peak islands for DRACH motif matches pinpoints target m6A loci."
)

ax_b.text(3.0, -0.32, callout_text, fontsize=9.5, color='#117864', ha='center', va='top',
          transform=ax_b.transData, clip_on=False,
          bbox=dict(boxstyle="round,pad=0.6", fc="#e8f8f5", ec="#17a589", lw=1.5))

# Save PNG and PDF
plt.savefig(out_png, dpi=300, bbox_inches='tight')
plt.savefig(out_pdf, bbox_inches='tight')
plt.close()

print(f"[SUCCESS] Figure saved to:\n  PNG: {out_png}\n  PDF: {out_pdf}")
