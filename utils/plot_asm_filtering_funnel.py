#!/usr/bin/env python3
# [[utils/plot_asm_filtering_funnel.py]]
# Goal: Visualize sequential ASM filtering with filter criteria on left of down arrows for Slide 9.

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

out_png = os.path.join(output_dir, "asm_filtering_funnel.png")
out_pdf = os.path.join(output_dir, "asm_filtering_funnel.pdf")
alt_png = os.path.join(project_root, "results", "03_variant_calling", "figures", "asm_filter_funnel.png")

# --- 2. Figure Setup ---
fig, ax = plt.subplots(figsize=(13.0, 6.4), dpi=300, facecolor="#ffffff")
plt.rcParams['font.sans-serif'] = 'DejaVu Sans'

ax.set_facecolor('#ffffff')
ax.set_xlim(0, 100)
ax.set_ylim(0, 100)
ax.axis('off')

# Concise Inverted Pyramid Tiers with Left Filter Criteria & Right Drop Counts
pyramid_tiers = [
    {
        "y": 76.0, "h": 18.0, "w": 90.0,
        "title": "1. Input Matrix",
        "metrics": "118,245 Rows  |  15,351 Loci  (39 Subjects)",
        "fc": "#ebf5fb", "ec": "#1f4e79", "tc": "#1b4f72",
        "is_highlight": False,
        "left_criterion": "IP >= 10, allele >= 2, q < 0.05, |AF-0.5| >= 0.2",
        "right_removed": "- 101,800 rows"
    },
    {
        "y": 52.0, "h": 16.0, "w": 70.0,
        "title": "2. Observation-Filtered",
        "metrics": "16,445 Rows Retained",
        "fc": "#fef9e7", "ec": "#f39c12", "tc": "#d35400",
        "is_highlight": False,
        "left_criterion": "Subjects >= 8 per locus (>= 2 MDD & >= 2 Ctrl)",
        "right_removed": "- 15,046 loci"
    },
    {
        "y": 28.0, "h": 16.0, "w": 52.0,
        "title": "3. Locus-Filtered",
        "metrics": "305 Loci Tested",
        "fc": "#e8f8f5", "ec": "#17a589", "tc": "#117864",
        "is_highlight": False,
        "left_criterion": "Sex-Adjusted Model (glmmTMB ~ Condition + Sex)",
        "right_removed": "FDR < 0.05"
    },
    {
        "y": 3.0, "h": 17.0, "w": 38.0,
        "title": "* 5 SIGNIFICANT DISCOVERIES *",
        "metrics": "MGAT5, BIN3 + IER2, TENM3, MTCH1",
        "fc": "#fef9e7", "ec": "#f1c40f", "tc": "#7d3c98",
        "is_highlight": True,
        "left_criterion": None,
        "right_removed": None
    }
]

# Render Vertical Pyramid Funnel
for idx, tier in enumerate(pyramid_tiers):
    x_left = 50.0 - tier["w"] / 2.0
    
    if tier["is_highlight"]:
        # Bold Highlight Container Box for Stage 4 Discoveries
        tier_box = patches.FancyBboxPatch((x_left, tier["y"]), tier["w"], tier["h"], boxstyle="round,pad=0.2",
                                          facecolor="#fef9e7", edgecolor="#f39c12", lw=2.8, zorder=3)
        ax.add_patch(tier_box)
        
        # Gold Banner Badge
        badge = patches.FancyBboxPatch((x_left + 2.0, tier["y"] + tier["h"] - 6.5), tier["w"] - 4.0, 5.5, boxstyle="round,pad=0.05",
                                       facecolor="#f39c12", edgecolor="none", zorder=4)
        ax.add_patch(badge)
        
        ax.text(50.0, tier["y"] + tier["h"] - 3.8, tier["title"], 
                fontsize=10.5, fontweight="bold", color="#ffffff", ha="center", va="center", zorder=5)
        
        ax.text(50.0, tier["y"] + 4.5, tier["metrics"], 
                fontsize=9.5, fontweight="bold", color="#6c3483", ha="center", va="center", zorder=5)
    else:
        # Standard Tier Box
        tier_box = patches.FancyBboxPatch((x_left, tier["y"]), tier["w"], tier["h"], boxstyle="round,pad=0.2",
                                          facecolor=tier["fc"], edgecolor=tier["ec"], lw=1.5, zorder=2)
        ax.add_patch(tier_box)
        
        ax.text(50.0, tier["y"] + tier["h"] - 5.0, tier["title"], 
                fontsize=10.5, fontweight="bold", color=tier["tc"], ha="center", zorder=3)
        
        ax.text(50.0, tier["y"] + 4.5, tier["metrics"], 
                fontsize=9.5, color="#2c3e50", ha="center", zorder=3)

    # Transition Arrow with Filter Criteria on LEFT and Drop Counts on RIGHT
    if tier["left_criterion"] is not None:
        next_y = pyramid_tiers[idx+1]["y"] + pyramid_tiers[idx+1]["h"]
        y_arrow_top = tier["y"]
        y_arrow_bot = next_y
        y_mid = (y_arrow_top + y_arrow_bot) / 2.0
        
        # Arrow down at X = 50.0
        ax.annotate("", xy=(50.0, y_arrow_bot + 0.4), xytext=(50.0, y_arrow_top - 0.4),
                    arrowprops=dict(arrowstyle="->", color="#c0392b" if "-" in tier["right_removed"] else "#34495e", lw=2.2), zorder=4)
        
        # Filter Criteria on LEFT of arrow (ha="right", x = 48.0)
        ax.text(48.0, y_mid, tier["left_criterion"], 
                fontsize=8.5, fontweight="bold", color="#c0392b" if "-" in tier["right_removed"] else "#2c3e50", 
                ha="right", va="center", zorder=5)
        
        # Drop Count / Output on RIGHT of arrow (ha="left", x = 52.0)
        ax.text(52.0, y_mid, tier["right_removed"], 
                fontsize=8.5, fontweight="bold", color="#c0392b" if "-" in tier["right_removed"] else "#16a085", 
                ha="left", va="center", zorder=5)

# Save PNG and PDF
plt.savefig(out_png, dpi=300, bbox_inches='tight')
plt.savefig(out_pdf, bbox_inches='tight')
if os.path.exists(os.path.dirname(alt_png)):
    plt.savefig(alt_png, dpi=300, bbox_inches='tight')
plt.close()

print(f"[SUCCESS] Symmetrical pyramid ASM filter funnel saved to:\n  PNG: {out_png}\n  PDF: {out_pdf}")
