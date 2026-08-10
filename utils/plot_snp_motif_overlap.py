#!/usr/bin/env python3
# [[utils/plot_snp_motif_overlap.py]]
# Goal: Visualize why one SNP site yields multiple SNP-DRACH overlaps
# (overlapping DRACH motif instances).

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

out_png = os.path.join(output_dir, "snp_motif_overlap.png")
out_pdf = os.path.join(output_dir, "snp_motif_overlap.pdf")

SNP_COLOR    = "#e74c3c"
MOTIF1_COLOR = "#27ae60"
MOTIF2_COLOR = "#1abc9c"

fig = plt.figure(figsize=(13.0, 8.0), dpi=300, facecolor="#ffffff")
gs = fig.add_gridspec(2, 1, height_ratios=[1.25, 0.75], hspace=0.32,
                      left=0.05, right=0.95, bottom=0.06, top=0.90)
fig.suptitle("Why 11,542 SNP sites yield 12,405 SNP\u2013DRACH overlaps",
             fontsize=13, fontweight="bold", color="#1a252f")

# ---------- TOP: sequence with overlapping motifs ----------
ax = fig.add_subplot(gs[0])
ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.axis('off')
ax.text(50, 92, "One SNP (red) falls inside TWO overlapping DRACH motifs",
        fontsize=10, fontweight="bold", color="#2c3e50", ha="center")

seq = list("AAGGACTGGAAC")
x0, bw, bg, bh = 8.0, 5.0, 1.0, 14.0
yseq = 52.0
motifA = [3, 4, 5, 6, 7]
motifB = [6, 7, 8, 9, 10]
snp_idx = 6

for i, letter in enumerate(seq):
    x = x0 + i * (bw + bg)
    if i == snp_idx:
        face, edge, tc = "#fdecea", SNP_COLOR, "#c0392b"
    elif i in motifA and i in motifB:
        face, edge, tc = "#d5f5e3", "#1e8449", "#145a32"
    elif i in motifA:
        face, edge, tc = "#e8f8f5", MOTIF1_COLOR, "#145a32"
    elif i in motifB:
        face, edge, tc = "#e0f2f1", MOTIF2_COLOR, "#0e6655"
    else:
        face, edge, tc = "#ecf0f1", "#bdc3c7", "#34495e"
    ax.add_patch(patches.FancyBboxPatch((x, yseq), bw, bh,
                 boxstyle="round,pad=0.05", facecolor=face, edgecolor=edge, lw=1.6, zorder=4))
    ax.text(x + bw/2, yseq + bh/2, letter, fontsize=9.5, fontweight="bold",
            color=tc, ha="center", va="center", zorder=5)

def xspan(idx_list):
    return (x0 + min(idx_list)*(bw+bg), x0 + max(idx_list)*(bw+bg) + bw)

xa1, xa2 = xspan(motifA)
xb1, xb2 = xspan(motifB)
xo1, xo2 = xspan([6, 7])

ax.add_patch(patches.FancyBboxPatch((xo1, 15), xo2 - xo1, 15.0,
             boxstyle="round,pad=0.05", facecolor="#fef9e7", edgecolor="#f39c12",
             lw=1.4, alpha=0.85, zorder=1))
ax.text((xo1 + xo2)/2, 36.5, "overlap", fontsize=8, fontweight="bold", color="#b9770e", ha="center")

ax.add_patch(patches.FancyBboxPatch((xa1, 22), xa2 - xa1, 5,
             boxstyle="round,pad=0.05", facecolor=MOTIF1_COLOR, edgecolor=MOTIF1_COLOR, lw=0, zorder=2))
ax.text((xa1 + xa2)/2, 28.5, "DRACH motif 1 (GGACT)", fontsize=8, fontweight="bold", color=MOTIF1_COLOR, ha="center")

ax.add_patch(patches.FancyBboxPatch((xb1, 15), xb2 - xb1, 5,
             boxstyle="round,pad=0.05", facecolor=MOTIF2_COLOR, edgecolor=MOTIF2_COLOR, lw=0, zorder=2))
ax.text((xb1 + xb2)/2, 21.5, "DRACH motif 2 (GGACT)", fontsize=8, fontweight="bold", color=MOTIF2_COLOR, ha="center")

xs = x0 + snp_idx*(bw + bg)
ax.annotate("SNP (A\u2192G)", xy=(xs + bw/2, yseq + bh + 3), xytext=(xs + bw/2 + 14, 82),
            arrowprops=dict(arrowstyle="->", color=SNP_COLOR, lw=1.8),
            color="#c0392b", fontsize=9, fontweight="bold", ha="center", zorder=6)

# ---------- BOTTOM: counting diagram ----------
axb = fig.add_subplot(gs[1])
axb.set_xlim(0, 100); axb.set_ylim(0, 100); axb.axis('off')

axb.add_patch(patches.FancyBboxPatch((8, 35), 26, 30, boxstyle="round,pad=0.2",
             facecolor="#ffffff", edgecolor=SNP_COLOR, lw=2.0))
axb.text(21, 56, "Distinct SNP site", fontsize=9, fontweight="bold", color="#2c3e50", ha="center")
axb.text(21, 44, "1", fontsize=16, fontweight="bold", color=SNP_COLOR, ha="center", va="center")

for yy in (55, 45):
    axb.annotate("", xy=(40, yy), xytext=(34, yy),
                 arrowprops=dict(arrowstyle="->", color="#7f8c8d", lw=1.8))

axb.add_patch(patches.FancyBboxPatch((40, 40), 24, 20, boxstyle="round,pad=0.2",
             facecolor="#e8f8f5", edgecolor=MOTIF1_COLOR, lw=1.8))
axb.text(52, 53, "overlap with motif 1", fontsize=8, fontweight="bold", color=MOTIF1_COLOR, ha="center")
axb.text(52, 44, "1", fontsize=12, fontweight="bold", color=MOTIF1_COLOR, ha="center")

axb.add_patch(patches.FancyBboxPatch((40, 18), 24, 20, boxstyle="round,pad=0.2",
             facecolor="#e0f2f1", edgecolor=MOTIF2_COLOR, lw=1.8))
axb.text(52, 31, "overlap with motif 2", fontsize=8, fontweight="bold", color=MOTIF2_COLOR, ha="center")
axb.text(52, 22, "1", fontsize=12, fontweight="bold", color=MOTIF2_COLOR, ha="center")

for yy in (50, 30):
    axb.annotate("", xy=(70, yy), xytext=(64, yy),
                 arrowprops=dict(arrowstyle="->", color="#7f8c8d", lw=1.8))

axb.add_patch(patches.FancyBboxPatch((70, 30), 24, 40, boxstyle="round,pad=0.2",
             facecolor="#fef9e7", edgecolor="#f39c12", lw=2.0))
axb.text(82, 62, "SNP\u2013motif overlaps", fontsize=9, fontweight="bold", color="#b9770e", ha="center")
axb.text(82, 46, "2", fontsize=16, fontweight="bold", color="#b9770e", ha="center")

plt.savefig(out_png, dpi=300, bbox_inches='tight')
plt.savefig(out_pdf, bbox_inches='tight')
plt.close()

print(f"[SUCCESS] SNP-motif overlap visualization saved to:\n  PNG: {out_png}\n  PDF: {out_pdf}")
