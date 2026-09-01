#!/usr/bin/env Rscript
# [[utils/plot_gtex_gwas_integration.R]]
# Goal: 2-Facet Psychiatric GWAS Trait Matrix distinguishing Direct SNP hits (Red) from Gene-level eQTL hits (Navy) without No Match in legend.

# Steps:
#   1. Load prioritized candidates from results/03_variant_calling/alt/tables/04_prioritized_mdd_snps_summary.csv
#   2. Filter valid protein-coding gene symbols (excluding non-coding antisense transcripts and MESDC2)
#   3. Dynamically extract Top 5 GTEx Brain eQTL leads & Top 5 GWAS psychiatric/suicidality leads
#   4. Differentiate Direct SNP GWAS hits (marked in Red) from Gene-level eQTL hits (Navy)
#   5. Exclude 'No Match' from legend display
#   6. Construct 2-facet ggplot matrix and export publication-ready PNG/PDF

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

# --- 1. Dynamic Pathing Setup ---
get_project_root <- function() {
  curr <- getwd()
  while (curr != "/" && !file.exists(file.path(curr, "run"))) {
    curr <- dirname(curr)
  }
  return(curr)
}

project_root <- Sys.getenv("PROJECT_ROOT", unset = get_project_root())
output_dir <- file.path(project_root, "results", "03_variant_calling", "alt", "figures")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

utils_fig_dir <- file.path(project_root, "results", "utils", "figures")
dir.create(utils_fig_dir, showWarnings = FALSE, recursive = TRUE)

out_png <- file.path(output_dir, "04_gtex_gwas_multomic_integration.png")
out_pdf <- file.path(output_dir, "04_gtex_gwas_multomic_integration.pdf")
out_png_utils <- file.path(utils_fig_dir, "gtex_gwas_multomic_integration.png")

# --- 2. Load Data ---
csv_path <- file.path(project_root, "results", "03_variant_calling", "alt", "tables", "04_prioritized_mdd_snps_summary.csv")
if (!file.exists(csv_path)) {
  csv_path <- file.path(project_root, "results", "03_variant_calling", "alt", "tables", "03_gwas_catalog_neuro_matches.csv")
}
if (!file.exists(csv_path)) {
  stop(paste("ERROR: Input table not found at", csv_path))
}

df <- read.csv(csv_path, stringsAsFactors = FALSE)

# Filter for standard, annotated protein-coding gene symbols (excluding non-coding antisense transcripts and MESDC2)
is_valid_gene <- function(g) {
  !is.na(g) & g != "" & 
    !grepl("^(RP|AC|CTD|AL|AP00|BX|STAG3L5P)", g) & 
    !grepl("(-AS[0-9]+|-DT|-OT|\\.[0-9]+)$", g) &
    !g %in% c("MESDC2", "CDIPT-AS1")
}
df_valid <- df[is_valid_gene(df$target_gene_symbol), ]

# --- 3. Dynamic Candidate Extraction ---
# A. Top 5 GTEx Brain eQTL Leads
gtex_top5_genes <- head(unique(df_valid$target_gene_symbol), 5)

# B. Top 5 GWAS Leads (Suicidality priority, non-redundant)
all_trait_text <- paste(
  df_valid$direct_snp_gwas_psych_traits, 
  df_valid$target_gene_gwas_psych_traits, 
  df_valid$target_gene_gwas_all_traits,
  sep = "; "
)
is_suic <- grepl("suicid", tolower(all_trait_text))
is_psych <- grepl("depress|mood|affective|neuroticism|stress", tolower(all_trait_text))

df_gwas <- df_valid[is_suic | is_psych, ]
df_gwas$suic_priority <- ifelse(
  grepl("suicid", tolower(paste(df_gwas$direct_snp_gwas_psych_traits, df_gwas$target_gene_gwas_psych_traits, df_gwas$target_gene_gwas_all_traits))), 
  2, 1
)
df_gwas_sorted <- df_gwas[order(-df_gwas$suic_priority, -df_gwas$n_mdd, -df_gwas$n_subjects), ]

gwas_candidates <- unique(df_gwas_sorted$target_gene_symbol)
gwas_top5_genes <- head(setdiff(gwas_candidates, gtex_top5_genes), 5)

message("[INFO] Dynamically extracted Top 5 GTEx leads: ", paste(gtex_top5_genes, collapse = ", "))
message("[INFO] Dynamically extracted Top 5 GWAS leads: ", paste(gwas_top5_genes, collapse = ", "))

# Combine selected dynamic cohorts
all_target_genes <- c(gtex_top5_genes, gwas_top5_genes)
hits <- df_valid[df_valid$target_gene_symbol %in% all_target_genes, ]
hits <- hits[!duplicated(hits$target_gene_symbol), ]

# Trait Domains (for Y-axis)
domains <- c(
  "Mood & Behavior",
  "Suicidality & Crisis",
  "Stress & Neuroticism",
  "Cognitive Function",
  "Brain MRI / Treatment"
)

matrix_rows <- list()

for (i in seq_len(nrow(hits))) {
  gene <- hits$target_gene_symbol[i]
  
  # Assign Facet Group dynamically
  facet_group <- ifelse(gene %in% gtex_top5_genes, 
                        "Top 5 GTEx Brain eQTL Leads", 
                        "Top 5 Psychiatric GWAS Leads (Suicidality Focus)")
  
  direct_t_low <- tolower(paste(hits$direct_snp_gwas_psych_traits[i], hits$direct_snp_gwas_all_traits[i], sep = "; "))
  gene_t_low   <- tolower(paste(hits$target_gene_gwas_psych_traits[i], hits$target_gene_gwas_all_traits[i], sep = "; "))
  
  # Direct checks
  d_direct <- c(
    grepl("depress|mood|externalizing|behavior", direct_t_low),
    grepl("suicid", direct_t_low),
    grepl("neuroticism|personality|trauma|stress", direct_t_low),
    grepl("cognit|intellect|memory|educational|executive|pursuit", direct_t_low),
    grepl("brain|morphology|sulcal|antipsychotic|t2 brain|mri|cortical|corpus callosum", direct_t_low)
  )
  
  # Gene-level checks
  d_gene <- c(
    grepl("depress|mood|externalizing|behavior", gene_t_low),
    grepl("suicid", gene_t_low),
    grepl("neuroticism|personality|trauma|stress", gene_t_low),
    grepl("cognit|intellect|memory|educational|executive|pursuit", gene_t_low),
    grepl("brain|morphology|sulcal|antipsychotic|t2 brain|mri|cortical|corpus callosum", gene_t_low)
  )
  
  for (d_idx in seq_along(domains)) {
    hit_status <- "No Match"
    if (d_direct[d_idx]) {
      hit_status <- "Direct SNP Match"
    } else if (d_gene[d_idx]) {
      hit_status <- "Gene-Level eQTL Match"
    }
    
    matrix_rows[[length(matrix_rows) + 1]] <- data.frame(
      gene = gene,
      facet_group = facet_group,
      domain = domains[d_idx],
      hit_type = hit_status
    )
  }
}

mat_df <- do.call(rbind, matrix_rows)

# Order factors
mat_df$facet_group <- factor(
  mat_df$facet_group, 
  levels = c("Top 5 GTEx Brain eQTL Leads", "Top 5 Psychiatric GWAS Leads (Suicidality Focus)")
)

mat_df$gene <- factor(mat_df$gene, levels = c(gtex_top5_genes, gwas_top5_genes))
mat_df$domain <- factor(mat_df$domain, levels = rev(domains))
mat_df$hit_type <- factor(mat_df$hit_type, levels = c("No Match", "Gene-Level eQTL Match", "Direct SNP Match"))

# Palette
matrix_pal <- c(
  "No Match" = "#f1f5f9",
  "Gene-Level eQTL Match" = "#1e3a8a",
  "Direct SNP Match" = "#dc2626"
)

# Active legend breaks excluding No Match
legend_breaks <- c("Gene-Level eQTL Match", "Direct SNP Match")

# --- 4. Construct 2-Facet GWAS Matrix Plot ---
p_faceted <- ggplot(mat_df, aes(x = gene, y = domain)) +
  geom_tile(aes(fill = hit_type), color = "#ffffff", linewidth = 2.0) +
  geom_text(
    aes(label = ifelse(hit_type != "No Match", "✓", "")),
    color = "#ffffff",
    fontface = "bold",
    size = 7.5
  ) +
  facet_grid(. ~ facet_group, scales = "free_x", space = "free_x") +
  scale_fill_manual(
    values = matrix_pal,
    breaks = legend_breaks,
    name = NULL
  ) +
  scale_x_discrete(position = "top") +
  labs(
    title = "Psychiatric GWAS Colocalization: GTEx Brain vs Suicidality Leads",
    subtitle = "Phenotype colocalization distinguishing Direct SNP hits (Red) from Gene-level eQTLs (Navy)",
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_family = "sans") +
  theme(
    plot.title = element_text(size = 18.0, face = "bold", color = "#0f172a", margin = margin(b = 4)),
    plot.subtitle = element_text(size = 13.0, face = "italic", color = "#64748b", margin = margin(b = 14)),
    strip.background = element_rect(fill = "#f8fafc", color = "#cbd5e1", linewidth = 0.8),
    strip.text = element_text(size = 14.0, face = "bold", color = "#1e293b", margin = margin(t = 6, b = 6)),
    panel.grid = element_blank(),
    axis.text.x.top = element_text(size = 14.0, face = "bold", color = "#0f172a", margin = margin(b = 6)),
    axis.text.y = element_text(size = 14.0, face = "bold", color = "#1e293b"),
    legend.position = "bottom",
    legend.text = element_text(size = 13.5, face = "bold", color = "#334155"),
    legend.spacing.x = unit(0.4, "cm"),
    plot.margin = margin(t = 16, r = 20, b = 12, l = 15)
  )

# --- 5. Export 2-Facet Single Plot ---
ggsave(out_png, plot = p_faceted, width = 14.0, height = 6.2, dpi = 300, bg = "white")
ggsave(out_pdf, plot = p_faceted, width = 14.0, height = 6.2, dpi = 300, bg = "white")
ggsave(out_png_utils, plot = p_faceted, width = 14.0, height = 6.2, dpi = 300, bg = "white")

message("[SUCCESS] 2-Facet GWAS matrix generated (No Match excluded from legend):")
message("  PNG: ", out_png)
message("  PDF: ", out_pdf)
