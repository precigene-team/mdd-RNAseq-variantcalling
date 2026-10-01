#!/usr/bin/env Rscript
# [[utils/plot_vep_functional_consequences.R]]
# Goal: 2-Facet visualization of Ensembl VEP molecular consequences for the exact 10 prioritized candidate loci.

# Steps:
#   1. Load functional annotation candidate variants from results/tables/10_functional_annotation.csv
#   2. Extract the exact 10 candidate loci:
#        - Facet 1: Top 5 GTEx Leads (SPON2, CUL9, EPG5, FXYD2, CASQ1)
#        - Facet 2: Top 5 GWAS Suicidality Leads (PLS1, MCUR1, ZNF585B, SPPL3, RGPD8)
#   3. Map VEP consequences and impact categories
#   4. Construct 2-facet ggplot card/tile plot with 2x typography
#   5. Export high-res PNG and PDF

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
utils_fig_dir <- file.path(project_root, "results", "utils", "figures")
dir.create(utils_fig_dir, showWarnings = FALSE, recursive = TRUE)

out_png_utils <- file.path(utils_fig_dir, "functional_annotations.png")
out_pdf_utils <- file.path(utils_fig_dir, "functional_annotations.pdf")

# --- 2. Load Data ---
csv_path <- file.path(project_root, "results", "03_variant_calling", "tables", "10_functional_annotation.csv")
if (!file.exists(csv_path)) {
  stop("ERROR: Input table not found: ", csv_path, ". Run 10_functional_annotation.py first.")
}

df <- read.csv(csv_path, stringsAsFactors = FALSE)
if (!"target_gene_symbol" %in% names(df) && "gene" %in% names(df)) df$target_gene_symbol <- df$gene
if (!"vep_consequence" %in% names(df) && "consequence" %in% names(df)) df$vep_consequence <- df$consequence

# Exact 10 Candidates
gtex_top5_genes <- c("SPON2", "CUL9", "EPG5", "FXYD2", "CASQ1")
gwas_top5_genes <- c("PLS1", "MCUR1", "ZNF585B", "SPPL3", "RGPD8")
target_10 <- c(gtex_top5_genes, gwas_top5_genes)

hits <- df[df$target_gene_symbol %in% target_10, ]
hits <- hits[!duplicated(hits$target_gene_symbol), ]
hits <- hits[match(target_10, hits$target_gene_symbol), ]

# Format clean consequence names
clean_consequence <- function(cons) {
  if (grepl("3_prime_UTR", cons)) return("3' UTR Regulatory")
  if (grepl("5_prime_UTR", cons)) return("5' UTR Regulatory")
  if (grepl("missense", cons)) return("Missense Coding")
  if (grepl("synonymous", cons)) return("Synonymous Exonic")
  if (grepl("non_coding_transcript_exon", cons)) return("Non-Coding Exon")
  if (grepl("intron", cons)) return("Intronic Regulatory")
  return(cons)
}

hits$clean_cons <- sapply(hits$vep_consequence, clean_consequence)
hits$facet_group <- ifelse(hits$target_gene_symbol %in% gtex_top5_genes,
                           "Top 5 GTEx Brain eQTL Leads",
                           "Top 5 Psychiatric GWAS Leads (Suicidality Focus)")

hits$gene_label <- paste0(hits$target_gene_symbol, "\n(", hits$rsid, ")")

# Factor ordering
hits$facet_group <- factor(
  hits$facet_group,
  levels = c("Top 5 GTEx Brain eQTL Leads", "Top 5 Psychiatric GWAS Leads (Suicidality Focus)")
)
hits$gene_label <- factor(hits$gene_label, levels = rev(hits$gene_label))

# Color palette for VEP consequences
cons_pal <- c(
  "Missense Coding" = "#dc2626",         # Red
  "3' UTR Regulatory" = "#1e3a8a",       # Navy
  "5' UTR Regulatory" = "#2563eb",       # Blue
  "Synonymous Exonic" = "#64748b",       # Slate
  "Non-Coding Exon" = "#94a3b8",         # Light Slate
  "Intronic Regulatory" = "#cbd5e1"      # Light Gray
)

# --- 3. Construct 2-Facet VEP Plot ---
p_vep <- ggplot(hits, aes(x = n_mdd, y = gene_label, fill = clean_cons)) +
  geom_col(width = 0.65) +
  geom_text(
    aes(label = paste0(clean_cons, "  [", n_mdd, " MDD : ", n_ctrl, " CTRL]")),
    hjust = -0.05,
    size = 4.4,
    fontface = "bold",
    color = "#0f172a"
  ) +
  facet_grid(facet_group ~ ., scales = "free_y", space = "free_y") +
  scale_fill_manual(values = cons_pal, name = "VEP Consequence") +
  scale_x_continuous(limits = c(0, 18), breaks = seq(0, 14, 2)) +
  labs(
    title = "Ensembl VEP Functional Annotation for Top 10 Prioritized Loci",
    subtitle = "Molecular consequences and MDD carrier burden for GTEx brain leads and GWAS suicidality leads",
    x = "MDD Donor Carrier Count (N = 20)",
    y = NULL
  ) +
  theme_minimal(base_family = "sans") +
  theme(
    plot.title = element_text(size = 17.5, face = "bold", color = "#0f172a", margin = margin(b = 4)),
    plot.subtitle = element_text(size = 12.5, face = "italic", color = "#64748b", margin = margin(b = 14)),
    strip.background = element_rect(fill = "#f8fafc", color = "#cbd5e1", linewidth = 0.8),
    strip.text = element_text(size = 13.5, face = "bold", color = "#1e293b", margin = margin(t = 5, b = 5)),
    panel.grid.major.x = element_line(color = "#e2e8f0", linewidth = 0.5, linetype = "dashed"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line.x = element_line(color = "#cbd5e1", linewidth = 0.7),
    axis.title.x = element_text(size = 13.5, face = "bold", color = "#334155"),
    axis.text.x = element_text(size = 12.0, color = "#475569"),
    axis.text.y = element_text(size = 13.0, face = "bold", color = "#0f172a"),
    legend.position = "bottom",
    legend.text = element_text(size = 12.0, face = "bold", color = "#334155"),
    legend.title = element_text(size = 12.5, face = "bold", color = "#0f172a"),
    plot.margin = margin(t = 16, r = 25, b = 12, l = 15)
  )

# --- 4. Export 2-Facet Single Plot ---
ggsave(out_png_utils, plot = p_vep, width = 14.0, height = 7.0, dpi = 300, bg = "white")
ggsave(out_pdf_utils, plot = p_vep, width = 14.0, height = 7.0, dpi = 300, bg = "white")

message("[SUCCESS] 2-Facet functional annotation figure generated successfully:")
message("  PNG: ", out_png_utils)
message("  PDF: ", out_pdf_utils)
