#!/usr/bin/env Rscript
# [[03_variant_calling/13_plot_clinical_annotation_landscape.R]]
# Goal: Publication-grade visualization of functional consequences, clinical pathogenicity, and pharmacogenomic actionability.
#
# Steps:
#   1. Dynamic project root resolution and directory setup.
#   2. Ingest processed annotation tables (10_functional, 11_clinical, 12_pharmacogenomics).
#   3. Panel A: SnpEff Sequence Ontology consequences & coding impact distribution.
#   4. Panel B: ACMG/AMP clinical pathogenicity breakdown & ClinGen validity.
#   5. Panel C: Case-control carrier recurrence skew & brain PFC eQTL support.
#   6. Assemble multi-panel figure using patchwork and export to results/figures/13_clinical_annotation_landscape.

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(readr)
  library(scales)
})

# --- Step 1: Dynamic Pathing ---
initial_options <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("--file=", initial_options, value = TRUE)
script_dir <- if (length(file_arg) > 0) {
  dirname(normalizePath(sub("--file=", "", file_arg[1])))
} else {
  getwd()
}
project_root <- Sys.getenv("PROJECT_ROOT")
if (project_root == "") {
  curr <- normalizePath(script_dir)
  while (curr != "/" && !file.exists(file.path(curr, "run")) && !file.exists(file.path(curr, "CONTEXT.md"))) {
    curr <- dirname(curr)
  }
  project_root <- curr
}

tables_dir  <- file.path(project_root, "results", "tables")
mod_tables_dir <- file.path(project_root, "results", "03_variant_calling", "tables")
figures_dir <- file.path(project_root, "results", "figures")
mod_figures_dir <- file.path(project_root, "results", "03_variant_calling", "figures")

dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(mod_figures_dir, recursive = TRUE, showWarnings = FALSE)

find_table <- function(fname) {
  paths <- c(file.path(tables_dir, fname), file.path(mod_tables_dir, fname))
  p <- paths[file.exists(paths)][1]
  if (is.na(p)) stop(sprintf("Table %s not found in production results.", fname))
  return(p)
}

func_csv <- find_table("10_functional_annotation.csv")
clin_csv <- find_table("11_clinical_pathogenicity.csv")
pgx_csv  <- find_table("12_pharmacogenomics_actionability.csv")

# --- Step 2: Ingestion & Assertions ---
df_func <- read_csv(func_csv, show_col_types = FALSE)
df_clin <- read_csv(clin_csv, show_col_types = FALSE)
df_pgx  <- read_csv(pgx_csv, show_col_types = FALSE)

cat(sprintf("[INFO] Loaded %d annotated candidate variants across %d genes.\n", nrow(df_pgx), n_distinct(df_pgx$gene)))

# Theme setup
pub_theme <- theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12, color = "#1a252f", margin = margin(b = 6)),
    plot.subtitle = element_text(size = 9.5, color = "#566573", margin = margin(b = 8)),
    axis.title = element_text(face = "bold", size = 9.5, color = "#2c3e50"),
    axis.text = element_text(size = 8.5, color = "#2c3e50"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "#eaeded", linewidth = 0.5),
    legend.title = element_text(face = "bold", size = 9),
    legend.text = element_text(size = 8),
    legend.position = "bottom"
  )

# --- Step 3: Panel A - Sequence Ontology Consequences ---
so_summary <- df_func %>%
  count(consequence, impact) %>%
  arrange(desc(n))

p_a <- ggplot(so_summary, aes(x = reorder(consequence, n), y = n, fill = impact)) +
  geom_col(width = 0.7, color = "#2c3e50", linewidth = 0.2) +
  geom_text(aes(label = sprintf("n=%d", n)), hjust = -0.2, size = 3.0, color = "#2c3e50") +
  coord_flip(clip = "off") +
  scale_fill_manual(
    values = c("HIGH" = "#c0392b", "MODERATE" = "#e67e22", "LOW" = "#27ae60", "MODIFIER" = "#2980b9"),
    name = "SnpEff Impact:"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
  labs(
    title = "A. Sequence Ontology Functional Consequence",
    subtitle = "SnpEff GRCh38.115 molecular consequence & coding impact",
    x = NULL,
    y = "Variant Count"
  ) +
  pub_theme

# --- Step 4: Panel B - ACMG Pathogenicity & ClinGen Validity ---
acmg_summary <- df_clin %>%
  count(acmg_classification, clingen_validity) %>%
  mutate(
    acmg_classification = factor(acmg_classification, 
                                 levels = c("Likely_Pathogenic", "VUS", "Likely_Benign", "Benign"))
  )

p_b <- ggplot(acmg_summary, aes(x = acmg_classification, y = n, fill = clingen_validity)) +
  geom_col(position = "stack", width = 0.65, color = "#2c3e50", linewidth = 0.2) +
  scale_fill_manual(
    values = c(
      "Definitive" = "#8e44ad",
      "Moderate"   = "#2980b9",
      "Limited"    = "#e67e22",
      "Not_Curated_in_ClinGen" = "#95a5a6"
    ),
    name = "ClinGen Validity:"
  ) +
  labs(
    title = "B. Clinical Pathogenicity Classification",
    subtitle = "ACMG/AMP consensus categories cross-referenced with ClinGen",
    x = "ACMG/AMP Classification",
    y = "Variant Count"
  ) +
  pub_theme +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

# --- Step 5: Panel C - MDD Carrier Recurrence Skew & PFC eQTL ---
p_c <- ggplot(df_pgx, aes(x = mdd_fraction, y = n_subjects, color = is_gtex_ba9_eqtl, size = n_mdd)) +
  geom_point(alpha = 0.85, position = position_jitter(width = 0.02, height = 0.2, seed = 42)) +
  scale_color_manual(
    values = c("YES" = "#27ae60", "NO" = "#7f8c8d"),
    name = "GTEx Brain PFC eQTL:"
  ) +
  scale_size_continuous(range = c(2, 6), name = "MDD Carriers:") +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "#bdc3c7") +
  labs(
    title = "C. Case-Control Carrier Recurrence",
    subtitle = "MDD carrier fraction vs total cohort subjects with eQTL support",
    x = "MDD Carrier Fraction (Skew)",
    y = "Total Cohort Subjects (N)"
  ) +
  pub_theme

# --- Step 6: Assemble & Export ---
full_plot <- (p_a) / (p_b | p_c) +
  plot_annotation(
    title = "Transcriptome-Wide Variant Clinical Interpretation Landscape",
    subtitle = "Prioritized Postmortem PFC Input RNA-seq Heterozygous SNPs (SnpEff GRCh38.115 -> ACMG -> CPIC)",
    theme = theme(
      plot.title = element_text(face = "bold", size = 15, color = "#1a252f"),
      plot.subtitle = element_text(size = 11, color = "#566573")
    )
  )

pdf_out <- file.path(figures_dir, "13_clinical_annotation_landscape.pdf")
png_out <- file.path(figures_dir, "13_clinical_annotation_landscape.png")
mod_pdf_out <- file.path(mod_figures_dir, "13_clinical_annotation_landscape.pdf")
mod_png_out <- file.path(mod_figures_dir, "13_clinical_annotation_landscape.png")

ggsave(pdf_out, full_plot, width = 12.5, height = 10.0, device = cairo_pdf)
ggsave(png_out, full_plot, width = 12.5, height = 10.0, dpi = 300)
ggsave(mod_pdf_out, full_plot, width = 12.5, height = 10.0, device = cairo_pdf)
ggsave(mod_png_out, full_plot, width = 12.5, height = 10.0, dpi = 300)

cat(sprintf("[SUCCESS] Clinical landscape figures exported to:\n  - %s\n  - %s\n", pdf_out, png_out))
