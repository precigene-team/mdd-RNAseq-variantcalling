#!/usr/bin/env python3
# [[utils/make_asm_5gene_summary_table.py]]
# Goal: Generate summary statistics CSV table for 5 discovery loci with 4 separate N breakdown columns (MDD_M, MDD_F, Ctrl_M, Ctrl_F).

import os
import sys
import pandas as pd
import numpy as np

# --- 1. Dynamic Pathing Setup ---
def get_project_root():
    curr = os.getcwd()
    while curr != "/" and not os.path.exists(os.path.join(curr, "run")):
        curr = os.path.dirname(curr)
    return curr

project_root = get_project_root()
matrix_file = os.path.join(project_root, "_data", "RData", "03_variant_calling", "asm_matrix.tsv")
sig_file = os.path.join(project_root, "results", "03_variant_calling", "tables", "12_asm_primary_significant.csv")

output_table_dir = os.path.join(project_root, "results", "utils", "tables")
os.makedirs(output_table_dir, exist_ok=True)

out_csv = os.path.join(output_table_dir, "03_asm_5gene_summary_table.csv")

# --- 2. Load Input Files ---
if not os.path.isfile(matrix_file):
    print(f"[ERROR] Matrix file not found: {matrix_file}", file=sys.stderr)
    sys.exit(1)

if not os.path.isfile(sig_file):
    print(f"[ERROR] Significant loci file not found: {sig_file}", file=sys.stderr)
    sys.exit(1)

df_matrix = pd.read_csv(matrix_file, sep="\t")
df_sig    = pd.read_csv(sig_file)

# Compute ALT Allele Fractions
df_matrix["ip_alt"] = df_matrix["ip_total"] - df_matrix["ip_ref"]
df_matrix["ip_af"]  = df_matrix["ip_alt"] / df_matrix["ip_total"]

df_matrix["in_alt"] = df_matrix["in_total"] - df_matrix["in_ref"]
df_matrix["in_af"]  = df_matrix["in_alt"] / df_matrix["in_total"]

target_genes = ["MGAT5", "BIN3", "IER2", "TENM3", "MTCH1"]

summary_rows = []

for gene in target_genes:
    df_g = df_matrix[df_matrix["gene"] == gene]
    sig_row = df_sig[df_sig["gene"] == gene].iloc[0]
    
    # Subject-level grouping to get unique subject counts and means
    subj_df = df_g.groupby(["subject", "condition", "sex"])[["in_af", "ip_af"]].mean().reset_index()
    
    n_tot    = len(subj_df)
    n_mdd_m  = len(subj_df[(subj_df["condition"] == "mdd")  & (subj_df["sex"] == "male")])
    n_mdd_f  = len(subj_df[(subj_df["condition"] == "mdd")  & (subj_df["sex"] == "female")])
    n_ctrl_m = len(subj_df[(subj_df["condition"] == "ctrl") & (subj_df["sex"] == "male")])
    n_ctrl_f = len(subj_df[(subj_df["condition"] == "ctrl") & (subj_df["sex"] == "female")])
    
    in_mean   = subj_df["in_af"].dropna().mean() * 100.0
    ip_mean   = subj_df["ip_af"].dropna().mean() * 100.0
    skew_pct  = ip_mean - in_mean
    
    ctrl_mean = subj_df[subj_df["condition"] == "ctrl"]["ip_af"].dropna().mean() * 100.0
    mdd_mean  = subj_df[subj_df["condition"] == "mdd"]["ip_af"].dropna().mean() * 100.0
    mdd_pct   = mdd_mean - ctrl_mean
    
    locus_coord = f"chr{sig_row['chr']}:{sig_row['start']}-{sig_row['end']}"
    
    # 4 separate N columns: N_MDD_M (M*), N_MDD_F (F*), N_Ctrl_M (M), N_Ctrl_F (F)
    summary_rows.append({
        "Gene": gene,
        "Locus": locus_coord,
        "N_Total": n_tot,
        "N_MDD_M": n_mdd_m,
        "N_MDD_F": n_mdd_f,
        "N_Ctrl_M": n_ctrl_m,
        "N_Ctrl_F": n_ctrl_f,
        "Input_AF%": round(in_mean, 1),
        "IP_AF%": round(ip_mean, 1),
        "Skew_%": f"{skew_pct:+.1f}%",
        "Ctrl_AF%": round(ctrl_mean, 1),
        "MDD_AF%": round(mdd_mean, 1),
        "MDD_Shift_%": f"{mdd_pct:+.1f}%",
        "Coef_beta": round(sig_row["coef"], 3),
        "FDR": f"{sig_row['fdr']:.4f}",
        "Confidence": sig_row["confidence"]
    })

df_out = pd.DataFrame(summary_rows)

# Save CSV
df_out.to_csv(out_csv, index=False)

print("================================================================================================================")
print("                                 SUMMARY STATISTICS FOR 5 DISCOVERY ASm6A LOCI                                 ")
print("================================================================================================================")
print(df_out.to_string(index=False))
print("================================================================================================================")
print(f"[SUCCESS] Summary CSV Table with 4 separate N breakdown columns saved to:\n  {out_csv}")
