#!/usr/bin/env python3
# [[03_variant_calling/16_generate_pygenometracks.py]]
# Goal: Generate pyGenomeTracks .ini configuration files and execution shell script for all 5 FDR-significant ASm6A loci.

import os
import sys
import pandas as pd

def get_project_root():
    curr = os.path.abspath(os.getcwd())
    while curr != "/" and not os.path.exists(os.path.join(curr, "run")):
        curr = os.path.dirname(curr)
    return curr

project_root = get_project_root()
sig_table_file = os.path.join(project_root, "results/03_variant_calling/tables/12_asm_primary_significant.csv")
output_dir = os.path.join(project_root, "03_variant_calling")
figures_dir = os.path.join(project_root, "results/03_variant_calling/figures")

os.makedirs(figures_dir, exist_ok=True)

# Define 5 discovery loci specifications
loci_info = [
    {
        "gene": "BIN3",
        "chr": "chr8",
        "chr_num": "8",
        "start": 22620417,
        "end": 22621067,
        "view_start": 22619800,
        "view_end": 22621500,
        "subject": "male_mdd_10",
        "ip_bam": "_data/bam/male_mdd_10_m6A/male_mdd_10_m6A_Aligned.sortedByCoord.out.bam",
        "in_bam": "_data/bam/male_mdd_10_input/male_mdd_10_input_Aligned.sortedByCoord.out.bam",
        "ip_bam_alt": "_data/bam_sliced/male_mdd_10_m6A_BIN3.bam",
        "in_bam_alt": "_data/bam_sliced/male_mdd_10_input_BIN3.bam",
        "vcf": "_data/variant/raw_vcf/cohort_filtered.vcf.gz"
    },
    {
        "gene": "MGAT5",
        "chr": "chr2",
        "chr_num": "2",
        "start": 134448632,
        "end": 134449882,
        "view_start": 134447500,
        "view_end": 134451000,
        "subject": "male_mdd_10",
        "ip_bam": "_data/bam/male_mdd_10_m6A/male_mdd_10_m6A_Aligned.sortedByCoord.out.bam",
        "in_bam": "_data/bam/male_mdd_10_input/male_mdd_10_input_Aligned.sortedByCoord.out.bam",
        "ip_bam_alt": "_data/bam_sliced/male_mdd_10_m6A_BIN3.bam",
        "in_bam_alt": "_data/bam_sliced/male_mdd_10_input_BIN3.bam",
        "vcf": "_data/variant/raw_vcf/cohort_filtered.vcf.gz"
    },
    {
        "gene": "IER2",
        "chr": "chr19",
        "chr_num": "19",
        "start": 13152110,
        "end": 13152460,
        "view_start": 13151500,
        "view_end": 13153000,
        "subject": "male_mdd_10",
        "ip_bam": "_data/bam/male_mdd_10_m6A/male_mdd_10_m6A_Aligned.sortedByCoord.out.bam",
        "in_bam": "_data/bam/male_mdd_10_input/male_mdd_10_input_Aligned.sortedByCoord.out.bam",
        "ip_bam_alt": "_data/bam_sliced/male_mdd_10_m6A_BIN3.bam",
        "in_bam_alt": "_data/bam_sliced/male_mdd_10_input_BIN3.bam",
        "vcf": "_data/variant/raw_vcf/cohort_filtered.vcf.gz"
    },
    {
        "gene": "TENM3",
        "chr": "chr4",
        "chr_num": "4",
        "start": 182754567,
        "end": 182755167,
        "view_start": 182753800,
        "view_end": 182756000,
        "subject": "male_mdd_10",
        "ip_bam": "_data/bam/male_mdd_10_m6A/male_mdd_10_m6A_Aligned.sortedByCoord.out.bam",
        "in_bam": "_data/bam/male_mdd_10_input/male_mdd_10_input_Aligned.sortedByCoord.out.bam",
        "ip_bam_alt": "_data/bam_sliced/male_mdd_10_m6A_BIN3.bam",
        "in_bam_alt": "_data/bam_sliced/male_mdd_10_input_BIN3.bam",
        "vcf": "_data/variant/raw_vcf/cohort_filtered.vcf.gz"
    },
    {
        "gene": "MTCH1",
        "chr": "chr6",
        "chr_num": "6",
        "start": 36986056,
        "end": 36986206,
        "view_start": 36985200,
        "view_end": 36987000,
        "subject": "male_mdd_10",
        "ip_bam": "_data/bam/male_mdd_10_m6A/male_mdd_10_m6A_Aligned.sortedByCoord.out.bam",
        "in_bam": "_data/bam/male_mdd_10_input/male_mdd_10_input_Aligned.sortedByCoord.out.bam",
        "ip_bam_alt": "_data/bam_sliced/male_mdd_10_m6A_BIN3.bam",
        "in_bam_alt": "_data/bam_sliced/male_mdd_10_input_BIN3.bam",
        "vcf": "_data/variant/raw_vcf/cohort_filtered.vcf.gz"
    }
]

# Generate .ini files
sh_commands = ["#!/usr/bin/env bash", "# Execution script for pyGenomeTracks across all 5 discovery loci\n"]

for item in loci_info:
    gene = item["gene"]
    ini_filename = f"16_{gene.lower()}_tracks.ini"
    ini_path = os.path.join(output_dir, ini_filename)
    
    # Check BAM file availability (full original BAM vs sliced BAM fallback)
    ip_bam_path = item["ip_bam"] if os.path.exists(os.path.join(project_root, item["ip_bam"])) else item["ip_bam_alt"]
    in_bam_path = item["in_bam"] if os.path.exists(os.path.join(project_root, item["in_bam"])) else item["in_bam_alt"]
    vcf_path = item["vcf"] if os.path.exists(os.path.join(project_root, item["vcf"])) else "_data/bam_sliced/bin3_filtered.vcf.gz"

    ini_content = f"""[x-axis]
fontsize = 12
where = top

[spacer]
height = 0.5

[IP m6A Read Coverage]
file = {ip_bam_path}
title = MeRIP-seq m6A IP ({gene} Locus)
height = 3.5
color = #8856a7
number_of_bins = 300
min_value = 0
file_type = bam

[spacer]
height = 0.3

[Input Baseline Read Coverage]
file = {in_bam_path}
title = Total RNA Input Baseline
height = 2.5
color = #999999
number_of_bins = 300
min_value = 0
file_type = bam

[spacer]
height = 0.3

[Heterozygous Variants]
file = {vcf_path}
title = GATK Filtered SNPs
height = 2.0
color = #d9534f
file_type = vcf
border_color = black

[spacer]
height = 0.5
"""

    with open(ini_path, "w") as f:
        f.write(ini_content)
    
    print(f"[INFO] Generated track INI: {ini_path}")

    # Build pyGenomeTracks command for this locus
    out_png = f"results/03_variant_calling/figures/16_igv_track_{gene.lower()}.png"
    out_pdf = f"results/03_variant_calling/figures/16_igv_track_{gene.lower()}.pdf"
    region_str = f"{item['chr']}:{item['view_start']}-{item['view_end']}"
    
    sh_commands.append(f"echo '>>> Generating pyGenomeTracks figure for locus: {gene} ({region_str})...'")
    sh_commands.append(f"pyGenomeTracks --tracks 03_variant_calling/{ini_filename} --region {region_str} --outFileName {out_png} --dpi 300")
    sh_commands.append(f"pyGenomeTracks --tracks 03_variant_calling/{ini_filename} --region {region_str} --outFileName {out_pdf}\n")

# Write master shell wrapper
sh_script_path = os.path.join(output_dir, "16_plot_igv_tracks_all.sh")
with open(sh_script_path, "w") as f:
    f.write("\n".join(sh_commands))

os.chmod(sh_script_path, 0o755)
print(f"[INFO] Generated master shell script: {sh_script_path}")
print("[DONE] pyGenomeTracks generator setup complete.")
