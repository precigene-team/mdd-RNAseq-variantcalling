# [[utils/qc02_cell_type_deconvolution.py]]
# Goal: Perform Marker Gene Profile (MGP) analysis to estimate relative cell-type proportions 
#       between MDD and Control subjects using VST-normalized expression counts.

import os
import pandas as pd
import numpy as np
from scipy import stats

project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
vst_counts_path = os.path.join(project_root, "results/02_dge_analysis/tables/02_RNAseq_vst_normalized_counts.csv")
sample_sheet_path = os.path.join(project_root, "_data/_metadata/sample_sheet.csv")

if not os.path.exists(vst_counts_path) or not os.path.exists(sample_sheet_path):
    print("Required data files not found. Ensure DGE has run.")
    exit(1)

# Define consensus brain cell-type marker genes (HGNC symbols)
cell_markers = {
    "Neurons": ["ENO2", "RBFOX3", "SNAP25", "SYT1", "MAP2"],
    "Astrocytes": ["GFAP", "ALDH1L1", "AQP4", "SOX9", "GJA1"],
    "Oligodendrocytes": ["MBP", "MOG", "MOBP", "OLIG2", "PLP1"],
    "Microglia": ["AIF1", "CD68", "TYROBP", "CSF1R", "C1QA"],
    "Endothelial": ["CLDN5", "PECAM1", "VWF", "CD34"]
}

# 1. Load VST counts
print("Loading VST-normalized expression matrix...")
df_vst = pd.read_csv(vst_counts_path)
df_vst.set_index("gene_name", inplace=True)

# 2. Load and clean metadata to match count columns
df_meta = pd.read_csv(sample_sheet_path)
# Keep only Input/RNA-seq samples from the metadata (matches "RNA-seq" or "input" case-insensitively)
df_meta = df_meta[df_meta["DataType"].str.lower().str.contains("rna-seq|input")].copy()
df_meta["Subject_id"] = df_meta["File1"].apply(lambda x: "_".join(x.split("_")[:3]) + "_input")

# Map Subject ID to SampleType
meta_dict = dict(zip(df_meta["Subject_id"], df_meta["SampleType"]))

# Filter out the 3 gated subjects from our metadata mapping
gated_subjects = ["female_mdd_2_input", "female_ctrl_6_input", "male_ctrl_1_input"]
sample_cols = [c for c in df_vst.columns if c in meta_dict and c not in gated_subjects]

print(f"Gating complete. Profiling {len(sample_cols)} subjects...")

# 3. Compute Marker Gene Profiles (MGP)
mgp_scores = {cell_type: [] for cell_type in cell_markers}
mgp_scores["Subject"] = sample_cols
mgp_scores["Group"] = [meta_dict[c] for c in sample_cols]
mgp_scores["Gender"] = ["female" if "female" in c else "male" for c in sample_cols]

for cell_type, markers in cell_markers.items():
    # Find markers present in our dataset
    valid_markers = [m for m in markers if m in df_vst.index]
    print(f"  {cell_type}: Found {len(valid_markers)} of {len(markers)} markers in dataset.")
    
    if not valid_markers:
        mgp_scores[cell_type] = [0.0] * len(sample_cols)
        continue
        
    # Extract marker expression
    expr = df_vst.loc[valid_markers, sample_cols]
    
    # Calculate MGP score per subject (first principal component or average expression)
    # Average Z-score is standard and highly robust
    expr_z = expr.apply(lambda x: (x - x.mean()) / (x.std() + 1e-8), axis=1)
    mean_scores = expr_z.mean(axis=0)
    
    mgp_scores[cell_type] = mean_scores.values

df_mgp = pd.DataFrame(mgp_scores)

# 4. Perform statistical tests (MDD vs Control)
print("\n" + "="*80)
print("🧠 CELL-TYPE HETEROGENEITY ANALYSIS (MDD VS CONTROL)")
print("="*80)

for cell_type in cell_markers.keys():
    ctrl_vals = df_mgp[df_mgp["Group"] == "ctrl"][cell_type]
    mdd_vals = df_mgp[df_mgp["Group"] == "mdd"][cell_type]
    
    # Overall t-test
    t_stat, p_val = stats.ttest_ind(ctrl_vals, mdd_vals, equal_var=False)
    
    # Stratified by sex
    males_ctrl = df_mgp[(df_mgp["Group"] == "ctrl") & (df_mgp["Gender"] == "male")][cell_type]
    males_mdd = df_mgp[(df_mgp["Group"] == "mdd") & (df_mgp["Gender"] == "male")][cell_type]
    t_stat_m, p_val_m = stats.ttest_ind(males_ctrl, males_mdd, equal_var=False)
    
    females_ctrl = df_mgp[(df_mgp["Group"] == "ctrl") & (df_mgp["Gender"] == "female")][cell_type]
    females_mdd = df_mgp[(df_mgp["Group"] == "mdd") & (df_mgp["Gender"] == "female")][cell_type]
    t_stat_f, p_val_f = stats.ttest_ind(females_ctrl, females_mdd, equal_var=False)
    
    print(f"\n📌 Cell Type: {cell_type}")
    print(f"  Overall Difference:  mean(Ctrl)={ctrl_vals.mean():.3f} vs mean(MDD)={mdd_vals.mean():.3f} (p-value = {p_val:.4f})")
    print(f"  Males Only:          mean(Ctrl)={males_ctrl.mean():.3f} vs mean(MDD)={males_mdd.mean():.3f} (p-value = {p_val_m:.4f})")
    print(f"  Females Only:        mean(Ctrl)={females_ctrl.mean():.3f} vs mean(MDD)={females_mdd.mean():.3f} (p-value = {p_val_f:.4f})")

print("="*80)

# Save scores to standardized results directory
out_path = os.path.join(project_root, "results/utils/tables/qc02_cell_type_mgp_scores.csv")
os.makedirs(os.path.dirname(out_path), exist_ok=True)
df_mgp.to_csv(out_path, index=False)
print(f"Saved MGP scores to: {out_path}\n")
