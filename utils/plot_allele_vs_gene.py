#!/usr/bin/env python3
# [[utils/plot_allele_vs_gene.py]]
# Goal: 2-Column Faceted Gene Plot showing both Allele Skew (Input vs IP) and Condition Shift (Control vs MDD).

import os
import sys
import numpy as np
import pandas as pd
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
matrix_file = os.path.join(project_root, "_data", "RData", "03_variant_calling", "asm_matrix.tsv")
sig_file = os.path.join(project_root, "results", "03_variant_calling", "tables", "12_asm_primary_significant.csv")

output_dir = os.path.join(project_root, "results", "utils", "figures")
os.makedirs(output_dir, exist_ok=True)

out_png = os.path.join(output_dir, "allele_vs_gene.png")
out_pdf = os.path.join(output_dir, "allele_vs_gene.pdf")

# --- 2. Load & Prepare Data ---
if not os.path.isfile(matrix_file):
    print(f"Error: Matrix file not found: {matrix_file}", file=sys.stderr)
    sys.exit(1)

df_matrix = pd.read_csv(matrix_file, sep="\t")

# Target 5 significant discovery genes
target_genes = ["MGAT5", "BIN3", "IER2", "TENM3", "MTCH1"]

# Filter matrix for target genes
df_sub = df_matrix[df_matrix["gene"].isin(target_genes)].copy()

# Compute ALT Allele Fractions: ALT = total - ref
df_sub["ip_alt"] = df_sub["ip_total"] - df_sub["ip_ref"]
df_sub["ip_af"]  = df_sub["ip_alt"] / df_sub["ip_total"]

df_sub["in_alt"] = df_sub["in_total"] - df_sub["in_ref"]
df_sub["in_af"]  = df_sub["in_alt"] / df_sub["in_total"]

# --- 3. Figure Layout (5 Rows x 2 Columns Faceted Grid) ---
fig, axes = plt.subplots(5, 2, figsize=(15.5, 9.5), dpi=300, facecolor="#ffffff", sharex=True)
plt.rcParams['font.sans-serif'] = 'DejaVu Sans'
plt.subplots_adjust(left=0.10, right=0.96, bottom=0.08, top=0.92, hspace=0.35, wspace=0.22)

# Color Palette:
col_ctrl = "#16a085" # Teal for Control
col_mdd  = "#e74c3c" # Coral Red for MDD
col_in   = "#7f8c8d" # Grey for Input baseline
col_ip   = "#8e44ad" # Purple for IP m6A pool

# Main Title
fig.suptitle("Allelic Skew (Input vs. IP) & Condition Shift (Control vs. MDD) Across 5 ASM Loci", 
             fontsize=13.5, fontweight="bold", color="#1a252f", y=0.96)

y_pos_bottom = 0.30
y_pos_top    = 0.60

for idx, gene in enumerate(target_genes):
    df_g = df_sub[df_sub["gene"] == gene]
    
    in_vals   = df_g["in_af"].dropna().values
    ip_vals   = df_g["ip_af"].dropna().values
    ctrl_vals = df_g[df_g["condition"] == "ctrl"]["ip_af"].dropna().values
    mdd_vals  = df_g[df_g["condition"] == "mdd"]["ip_af"].dropna().values
    
    # -------------------------------------------------------------------------
    # COLUMN 1 (LEFT): ALLELE SKEW (Input Baseline vs. IP Methylated Pool)
    # -------------------------------------------------------------------------
    ax_l = axes[idx, 0]
    ax_l.set_facecolor('#ffffff')
    ax_l.grid(True, linestyle="--", alpha=0.3, color="#bdc3c7", zorder=1)
    
    # Neutral reference line at 0.50
    ax_l.axvline(0.5, color="#7f8c8d", linestyle=":", lw=1.5, alpha=0.7, zorder=1)
    
    np.random.seed(100 + idx)
    j_in = np.random.uniform(-0.06, 0.06, len(in_vals))
    j_ip = np.random.uniform(-0.06, 0.06, len(ip_vals))
    
    ax_l.scatter(in_vals, y_pos_bottom + j_in, color=col_in, alpha=0.65, s=40, zorder=3, label="Input Baseline" if idx==0 else "")
    ax_l.scatter(ip_vals, y_pos_top + j_ip,    color=col_ip, alpha=0.75, s=45, zorder=3, label="IP Methylated Pool" if idx==0 else "")
    
    m_in = np.mean(in_vals) if len(in_vals) > 0 else 0.5
    m_ip = np.mean(ip_vals) if len(ip_vals) > 0 else 0.5
    delta_skew = m_ip - m_in
    
    ax_l.plot([m_in, m_in], [y_pos_bottom - 0.12, y_pos_bottom + 0.12], color=col_in, lw=3.0, zorder=4)
    ax_l.plot([m_ip, m_ip], [y_pos_top - 0.12,    y_pos_top + 0.12],    color=col_ip, lw=3.0, zorder=4)
    
    # Arrow connecting Input -> IP skew
    ax_l.annotate("", xy=(m_ip, 0.45), xytext=(m_in, 0.45),
                 arrowprops=dict(arrowstyle="->", color=col_ip, lw=1.8), zorder=2)
    
    sign_skew = "+" if delta_skew > 0 else ""
    badge_skew = f"IP Skew: {sign_skew}{delta_skew:.2f}"
    bx_l = max(m_ip, m_in) + 0.03
    if bx_l > 0.85:
        bx_l = min(m_ip, m_in) - 0.18
    ax_l.text(bx_l, 0.45, badge_skew, fontsize=8.0, fontweight="bold", color=col_ip, va="center",
              bbox=dict(boxstyle="round,pad=0.2", facecolor="#f5eef8", edgecolor="#8e44ad", lw=0.9), zorder=5)

    ax_l.set_ylim(0.0, 1.0)
    ax_l.set_yticks([y_pos_bottom, y_pos_top])
    ax_l.set_yticklabels(["Input", "IP Pool"], fontsize=8.5, fontweight="bold", color="#34495e")
    ax_l.text(-0.16, 0.50, f"Locus: {gene}", transform=ax_l.transAxes, fontsize=10.0, fontweight="bold", color="#1a252f", va="center", ha="right")

    # -------------------------------------------------------------------------
    # COLUMN 2 (RIGHT): CONDITION SHIFT (Control IP vs. MDD IP)
    # -------------------------------------------------------------------------
    ax_r = axes[idx, 1]
    ax_r.set_facecolor('#ffffff')
    ax_r.grid(True, linestyle="--", alpha=0.3, color="#bdc3c7", zorder=1)
    
    # Own Input Baseline reference line
    ax_r.axvline(m_in, color=col_in, linestyle="--", lw=1.5, alpha=0.7, zorder=1)
    
    np.random.seed(42 + idx)
    j_ctrl = np.random.uniform(-0.06, 0.06, len(ctrl_vals))
    j_mdd  = np.random.uniform(-0.06, 0.06, len(mdd_vals))
    
    ax_r.scatter(ctrl_vals, y_pos_bottom + j_ctrl, color=col_ctrl, alpha=0.75, s=45, zorder=3, label="Control IP" if idx==0 else "")
    ax_r.scatter(mdd_vals,  y_pos_top + j_mdd,     color=col_mdd,  alpha=0.75, s=45, zorder=3, label="MDD IP" if idx==0 else "")
    
    m_ctrl = np.mean(ctrl_vals) if len(ctrl_vals) > 0 else m_in
    m_mdd  = np.mean(mdd_vals)  if len(mdd_vals) > 0 else m_in
    delta_mdd = m_mdd - m_ctrl
    
    ax_r.plot([m_ctrl, m_ctrl], [y_pos_bottom - 0.12, y_pos_bottom + 0.12], color=col_ctrl, lw=3.0, zorder=4)
    ax_r.plot([m_mdd, m_mdd],   [y_pos_top - 0.12,    y_pos_top + 0.12],   color=col_mdd,  lw=3.0, zorder=4)
    
    # Arrow connecting Control -> MDD shift
    ax_r.annotate("", xy=(m_mdd, 0.45), xytext=(m_ctrl, 0.45),
                 arrowprops=dict(arrowstyle="->", color="#2c3e50", lw=1.8), zorder=2)
    
    sign_mdd = "+" if delta_mdd > 0 else ""
    badge_mdd = f"MDD Shift: {sign_mdd}{delta_mdd:.2f}"
    bx_r = max(m_mdd, m_ctrl) + 0.03
    if bx_r > 0.85:
        bx_r = min(m_mdd, m_ctrl) - 0.18
    ax_r.text(bx_r, 0.45, badge_mdd, fontsize=8.0, fontweight="bold", color="#2c3e50", va="center",
              bbox=dict(boxstyle="round,pad=0.2", facecolor="#f4f6f7", edgecolor="#bdc3c7", lw=0.9), zorder=5)

    ax_r.set_ylim(0.0, 1.0)
    ax_r.set_yticks([y_pos_bottom, y_pos_top])
    ax_r.set_yticklabels(["Control", "MDD"], fontsize=8.5, fontweight="bold", color="#34495e")

# Axis Labels & Column Titles
axes[-1, 0].set_xlim(0.0, 1.0)
axes[-1, 1].set_xlim(0.0, 1.0)

axes[-1, 0].set_xlabel(r"ALT Allele Fraction ($\mathrm{ALT}_{\mathrm{af}}$)", fontsize=10.5, fontweight="bold", color="#1a252f", labelpad=8)
axes[-1, 1].set_xlabel(r"ALT Allele Fraction ($\mathrm{ALT}_{\mathrm{af}}$)", fontsize=10.5, fontweight="bold", color="#1a252f", labelpad=8)

axes[0, 0].set_title("A. Allele Skew (Input Baseline vs. IP Pool)", fontsize=10.5, fontweight="bold", color="#1a252f", pad=10)
axes[0, 1].set_title("B. Condition Shift (Control IP vs. MDD IP)", fontsize=10.5, fontweight="bold", color="#1a252f", pad=10)

axes[0, 0].legend(loc="upper right", frameon=True, facecolor="#ffffff", edgecolor="#bdc3c7", fontsize=8.0)
axes[0, 1].legend(loc="upper right", frameon=True, facecolor="#ffffff", edgecolor="#bdc3c7", fontsize=8.0)

# Save PNG and PDF
plt.savefig(out_png, dpi=300, bbox_inches='tight')
plt.savefig(out_pdf, bbox_inches='tight')
plt.close()

print(f"[SUCCESS] 2-Column Faceted Allele vs Gene plot (showing BOTH Allele Skew and Condition Shift per gene) saved to:\n  PNG: {out_png}\n  PDF: {out_pdf}")
