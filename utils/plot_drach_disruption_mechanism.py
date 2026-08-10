#!/usr/bin/env python3
# [[utils/plot_drach_disruption_mechanism.py]]
# Goal: Visualize DRACH motif disruption mechanism (ASm6A).
# Layout: 2x2 grid
#   Left column  = RAW DRACH motif (sequence before any modification)
#   Right column = outcome after m6A methylation (CH3 marked on intact motif)
#   Top row = intact (maternal), Bottom row = disrupted (paternal)
#   Full-width bottom band = ASm6A consequence (skewed IP ratio)

import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches

plt.rcParams['mathtext.fontset'] = 'dejavusans'

def get_project_root():
    curr = os.getcwd()
    while curr != "/" and not os.path.exists(os.path.join(curr, "run")):
        curr = os.path.dirname(curr)
    return curr

project_root = get_project_root()
output_dir = os.path.join(project_root, "results", "utils", "figures")
os.makedirs(output_dir, exist_ok=True)

out_png = os.path.join(output_dir, "drach_disruption_mechanism.png")
out_pdf = os.path.join(output_dir, "drach_disruption_mechanism.pdf")

# Color Palette:
MATERNAL_COLOR = "#e91e63" # Pink
PATERNAL_COLOR = "#2980b9" # Blue
DRACH_COLOR    = "#27ae60" # Green for Intact DRACH motif
DISRUPT_COLOR  = "#e74c3c" # Red for Disrupted DRACH motif
ENZYME_COLOR   = "#8e44ad" # Purple for METTL3/14 m6A Writer Complex
METHYL_COLOR   = "#f39c12" # Gold for m6A methylation (CH3) stamp

def draw_seq_track(ax, y_center, base_seq, drach_indices, target_idx=None,
                   is_disrupted=False, is_maternal=True, show_ch3=False):
    """Draw a 7-base sequence track as rounded boxes. Optionally stamp CH3 (or 'no CH3') on the target base."""
    x_start = 24.0
    box_w, box_gap, box_h = 6.2, 1.0, 8.0

    reg_color = MATERNAL_COLOR if is_maternal else PATERNAL_COLOR
    reg_bg    = "#fce4ec" if is_maternal else "#ebf5fb"
    reg_text  = "#c2185b" if is_maternal else "#1b4f72"

    x_end_main = x_start + len(base_seq) * (box_w + box_gap)

    # Faded flanking boxes (neighboring sequence context)
    for k in range(2):
        xc = x_start - (k + 1) * (box_w + box_gap)
        ax.add_patch(patches.FancyBboxPatch((xc, y_center - box_h/2.0), box_w, box_h,
                     boxstyle="round,pad=0.05", facecolor=reg_bg, edgecolor=reg_color,
                     lw=0.9, linestyle="--", alpha=0.35, zorder=2))
        xc2 = x_end_main + k * (box_w + box_gap)
        ax.add_patch(patches.FancyBboxPatch((xc2, y_center - box_h/2.0), box_w, box_h,
                     boxstyle="round,pad=0.05", facecolor=reg_bg, edgecolor=reg_color,
                     lw=0.9, linestyle="--", alpha=0.35, zorder=2))

    base_boxes = {}
    for i, letter in enumerate(base_seq):
        x_pos = x_start + i * (box_w + box_gap)
        if is_disrupted and i == target_idx:
            bg, ec, tc, lw = DISRUPT_COLOR, "#922b21", "#ffffff", 2.2
        elif i in drach_indices:
            bg, ec, tc, lw = "#e8f8f5", DRACH_COLOR, "#1e8449", 1.6
        else:
            bg, ec, tc, lw = reg_bg, reg_color, reg_text, 1.3
        ax.add_patch(patches.FancyBboxPatch((x_pos, y_center - box_h/2.0), box_w, box_h,
                     boxstyle="round,pad=0.05", facecolor=bg, edgecolor=ec, lw=lw, zorder=4))
        ax.text(x_pos + box_w/2.0, y_center + box_h/2.0, letter, fontsize=10,
                fontweight="bold", color=tc, ha="center", va="center", zorder=5)
        base_boxes[i] = (x_pos, y_center - box_h/2.0, box_w, box_h)

    if show_ch3 and target_idx is not None:
        if is_disrupted:
            draw_no_ch3(ax, *base_boxes[target_idx])
        else:
            draw_ch3(ax, *base_boxes[target_idx])

    return x_start, box_w, box_gap

def draw_ch3(ax, x, y, w, h):
    """Stamp a CH3 (m6A) badge above the methylated base."""
    bx = x + w/2.0
    by = y + h + 1.0
    badge = patches.FancyBboxPatch((bx - 2.7, by), 5.4, 3.4, boxstyle="round,pad=0.15",
                                   facecolor="#fef9e7", edgecolor=METHYL_COLOR, lw=1.5, zorder=7)
    ax.add_patch(badge)
    ax.plot([bx, bx], [y + h - 0.5, by + 1.0], color=METHYL_COLOR, lw=1.3, zorder=6)
    ax.text(bx, by + 1.7, "CH3", fontsize=8, fontweight="bold", color="#b9770e",
            ha="center", va="center", zorder=8)
    ax.text(bx, by + 4.8, "m6A", fontsize=8, fontweight="bold", color="#b9770e",
            ha="center", va="center", zorder=8)

def draw_no_ch3(ax, x, y, w, h):
    """Greyed-out 'no CH3' badge above an unmodified base."""
    bx = x + w/2.0
    by = y + h + 1.0
    badge = patches.FancyBboxPatch((bx - 2.7, by), 5.4, 3.4, boxstyle="round,pad=0.15",
                                   facecolor="#eaeded", edgecolor="#7f8c8d", lw=1.2,
                                   linestyle="--", zorder=7)
    ax.add_patch(badge)
    ax.text(bx, by + 1.7, "no CH3", fontsize=6.5, fontweight="bold", color="#566573",
            ha="center", va="center", zorder=8)

def draw_bracket(ax, idx_list, x_start, box_w, box_gap, y_top, label, color, ls="-"):
    x0 = x_start + min(idx_list) * (box_w + box_gap)
    x1 = x_start + max(idx_list) * (box_w + box_gap) + box_w
    ax.plot([x0, x1], [y_top, y_top], color=color, lw=1.6, ls=ls, zorder=6)
    ax.plot([x0, x0], [y_top - 3.0, y_top], color=color, lw=1.6, zorder=6)
    ax.plot([x1, x1], [y_top - 3.0, y_top], color=color, lw=1.6, zorder=6)
    ax.text((x0 + x1)/2.0, y_top + 2.5, label, fontsize=7.5, fontweight="bold",
            color=color, ha="center", zorder=6)

def draw_enzyme(ax, label_text, blocked=False, y_bottom=20.0, cx=50.0):
    w, h = 20.0, 13.0
    if blocked:
        face, edge, color = "#fadbd8", DISRUPT_COLOR, "#922b21"
    else:
        face, edge, color = "#f5eef8", ENZYME_COLOR, ENZYME_COLOR
    ax.add_patch(patches.FancyBboxPatch((cx - w/2.0, y_bottom), w, h,
                 boxstyle="round,pad=0.2", facecolor=face, edgecolor=edge, lw=1.8, zorder=6))
    ax.text(cx, y_bottom + h/2.0, label_text, fontsize=8.5, fontweight="bold",
            color=color, ha="center", va="center", zorder=7)
    ax.annotate("", xy=(cx, y_bottom + h + 4.0), xytext=(cx, y_bottom + h),
                arrowprops=dict(arrowstyle="->", color=edge, lw=2.0,
                                linestyle="--" if blocked else "-"), zorder=6)

def setup_axes(ax):
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    ax.axis('off')
    ax.add_patch(patches.FancyBboxPatch((1, 1), 98, 98, boxstyle="round,pad=0.2",
                 facecolor="#ffffff", edgecolor="#bdc3c7", lw=1.3))

# --- Figure Setup (2x2 grid + bottom consequence band) ---
INTACT = dict(base=["C", "G", "G", "A", "C", "U", "A"], drach=[1, 2, 3, 4, 5], target=3)
DISRUPTED = dict(base=["C", "G", "G", "G", "C", "U", "A"], drach=[1, 2, 4, 5], target=3)

fig = plt.figure(figsize=(13.0, 10.0), dpi=300, facecolor="#ffffff")
gs = fig.add_gridspec(3, 2, height_ratios=[1.0, 1.0, 0.22], hspace=0.42, wspace=0.12,
                      left=0.06, right=0.94, bottom=0.06, top=0.88)

fig.suptitle("Molecular Mechanism: DRACH Motif Disruption & Allele-Specific m6A (ASm6A)",
             fontsize=13, fontweight="bold", color="#1a252f", y=0.99)

# Column headers
fig.text(0.28, 0.93, "1. Raw DRACH Motif", fontsize=10, fontweight="bold", color="#2c3e50", ha="center")
fig.text(0.72, 0.93, "2. After m6A Methylation", fontsize=10, fontweight="bold", color="#2c3e50", ha="center")

# ============================ ROW 1: INTACT (MATERNAL) ============================
ax00 = fig.add_subplot(gs[0, 0])
setup_axes(ax00)
ax00.text(50.0, 90.0, "Intact DRACH (Maternal copy)", fontsize=10, fontweight="bold", color=DRACH_COLOR, ha="center")
ax00.text(4.0, 42.0, "Maternal Copy:", fontsize=8.5, fontweight="bold", color=MATERNAL_COLOR)
x0, bw, bg = draw_seq_track(ax00, 52.0, INTACT["base"], INTACT["drach"], is_disrupted=False, is_maternal=True)
draw_bracket(ax00, INTACT["drach"], x0, bw, bg, 70.0, "Intact DRACH motif (GGACU)", DRACH_COLOR)

ax01 = fig.add_subplot(gs[0, 1])
setup_axes(ax01)
ax01.text(50.0, 90.0, "Methylated (Maternal copy)", fontsize=10, fontweight="bold", color=METHYL_COLOR, ha="center")
ax01.text(4.0, 42.0, "Maternal Copy:", fontsize=8.5, fontweight="bold", color=MATERNAL_COLOR)
x0, bw, bg = draw_seq_track(ax01, 52.0, INTACT["base"], INTACT["drach"], target_idx=INTACT["target"],
                            is_disrupted=False, is_maternal=True, show_ch3=True)
draw_bracket(ax01, INTACT["drach"], x0, bw, bg, 70.0, "Motif recognized", DRACH_COLOR)
draw_enzyme(ax01, "METTL3/14\nbinds", blocked=False, y_bottom=22.0)
ax01.add_patch(patches.FancyBboxPatch((10, 6), 80, 10, boxstyle="round,pad=0.1",
                 facecolor="#e8f8f5", edgecolor=DRACH_COLOR, lw=1.2))
ax01.text(50.0, 11.0, "✓ METHYLATED — CH3 deposited on the A", fontsize=8.5, fontweight="bold",
          color="#1e8449", ha="center", va="center")

# ========================= ROW 2: DISRUPTED (PATERNAL) ===========================
ax10 = fig.add_subplot(gs[1, 0])
setup_axes(ax10)
ax10.text(50.0, 90.0, "Disrupted DRACH (Paternal copy)", fontsize=10, fontweight="bold", color=DISRUPT_COLOR, ha="center")
ax10.text(4.0, 42.0, "Paternal Copy:", fontsize=8.5, fontweight="bold", color=PATERNAL_COLOR)
x0, bw, bg = draw_seq_track(ax10, 52.0, DISRUPTED["base"], DISRUPTED["drach"], target_idx=DISRUPTED["target"],
                            is_disrupted=True, is_maternal=False)
draw_bracket(ax10, DISRUPTED["drach"], x0, bw, bg, 70.0, "Disrupted motif (GGGCU)", DISRUPT_COLOR, ls="--")
ax10.text(50.0, 24.0, "SNP A → G breaks the methylation code", fontsize=8.5, fontweight="bold",
          color="#922b21", ha="center")

ax11 = fig.add_subplot(gs[1, 1])
setup_axes(ax11)
ax11.text(50.0, 90.0, "Unmethylated (Paternal copy)", fontsize=10, fontweight="bold", color="#7f8c8d", ha="center")
ax11.text(4.0, 42.0, "Paternal Copy:", fontsize=8.5, fontweight="bold", color=PATERNAL_COLOR)
x0, bw, bg = draw_seq_track(ax11, 52.0, DISRUPTED["base"], DISRUPTED["drach"], target_idx=DISRUPTED["target"],
                            is_disrupted=True, is_maternal=False, show_ch3=True)
draw_bracket(ax11, DISRUPTED["drach"], x0, bw, bg, 70.0, "Motif NOT recognized", DISRUPT_COLOR, ls="--")
draw_enzyme(ax11, "METTL3/14\nbinding BLOCKED ✗", blocked=True, y_bottom=22.0)
ax11.add_patch(patches.FancyBboxPatch((10, 6), 80, 10, boxstyle="round,pad=0.1",
                 facecolor="#f9ebea", edgecolor=DISRUPT_COLOR, lw=1.2))
ax11.text(50.0, 11.0, "✗ UNMETHYLATED — no CH3 deposited", fontsize=8.5, fontweight="bold",
          color="#922b21", ha="center", va="center")

# ========================= BOTTOM CONSEQUENCE BAND ================================
axb = fig.add_subplot(gs[2, :])
axb.set_xlim(0, 100)
axb.set_ylim(0, 100)
axb.axis('off')
axb.add_patch(patches.FancyBboxPatch((1, 8), 98, 84, boxstyle="round,pad=0.15",
             facecolor="#ffffff", edgecolor="#27ae60", lw=1.6))
axb.text(50.0, 62, "CONSEQUENCE: ALLELE-SPECIFIC METHYLATION (ASm6A)", fontsize=9,
         fontweight="bold", color="#27ae60", ha="center")
axb.text(50.0, 32, "Maternal allele methylated (CH3), paternal allele not → antibody IP reads SKEW toward the methylated allele\n(e.g., 85% A vs 15% G) → significant allelic imbalance in the methylated (IP) fraction.",
         fontsize=8, fontweight="bold", color="#1f4e79", ha="center", va="center")

# Save PNG and PDF
plt.savefig(out_png, dpi=300, bbox_inches='tight')
plt.savefig(out_pdf, bbox_inches='tight')
plt.close()

print(f"[SUCCESS] DRACH Motif Disruption mechanism visualization saved to:\n  PNG: {out_png}\n  PDF: {out_pdf}")
