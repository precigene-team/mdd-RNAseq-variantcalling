#!/usr/bin/env Rscript
# [[utils/plot_gtex_eqtl_summary.R]]
# Goal: Ultra-clean, minimal single-panel volcano plot of GTEx Brain BA9 eQTL matches with enrichment counts in legend and criteria in subtitle.

# Steps:
#   1. Load GTEx eQTL match results from results/03_variant_calling/alt/tables/02_gtex_brain_eqtl_matches.csv
#   2. Filter validated BA9 eQTL hits and calculate log10 p-values
#   3. Assign sub-threshold (< 1e-5) points to light gray
#   4. Calculate significant loci counts for MDD-Enriched (N=39) and CTRL-Enriched (N=13)
#   5. Construct volcano plot with enriched carrier criteria in subtitle and counts in legend
#   6. Export publication-ready PNG/PDF

suppressPackageStartupMessages({
  library(ggplot2)
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

out_png <- file.path(output_dir, "02_gtex_brain_eqtl_summary.png")
out_pdf <- file.path(output_dir, "02_gtex_brain_eqtl_summary.pdf")
out_png_utils <- file.path(utils_fig_dir, "gtex_brain_eqtl_summary.png")

# --- 2. Load and Prepare GTEx Match Data ---
csv_path <- file.path(project_root, "results", "03_variant_calling", "alt", "tables", "02_gtex_brain_eqtl_matches.csv")
if (!file.exists(csv_path)) {
  stop(paste("ERROR: Input table not found at", csv_path))
}

df <- read.csv(csv_path, stringsAsFactors = FALSE)
eqtls <- df[df$is_gtex_ba9_eqtl == "YES", ]

if (nrow(eqtls) == 0) {
  stop("ERROR: No validated BA9 eQTL hits found in table.")
}

eqtls$ba9_eqtl_pvalue <- as.numeric(eqtls$ba9_eqtl_pvalue)
eqtls$ba9_nes <- as.numeric(eqtls$ba9_nes)
eqtls$neg_log10p <- -log10(eqtls$ba9_eqtl_pvalue)

# Count significant loci (p < 1e-5)
n_mdd_sig  <- sum(eqtls$n_mdd >= 7 & eqtls$n_ctrl <= 2 & eqtls$neg_log10p >= 5.0)
n_ctrl_sig <- sum(eqtls$n_ctrl >= 7 & eqtls$n_mdd <= 2 & eqtls$neg_log10p >= 5.0)

mdd_label  <- paste0("MDD-Enriched (N = ", n_mdd_sig, ")")
ctrl_label <- paste0("CTRL-Enriched (N = ", n_ctrl_sig, ")")

# Categorize segregation group
eqtls$group <- "Sub-Threshold (p > 1e-5)"
eqtls$group[eqtls$n_mdd >= 7 & eqtls$n_ctrl <= 2 & eqtls$neg_log10p >= 5.0] <- mdd_label
eqtls$group[eqtls$n_ctrl >= 7 & eqtls$n_mdd <= 2 & eqtls$neg_log10p >= 5.0] <- ctrl_label

eqtls$group <- factor(
  eqtls$group, 
  levels = c(mdd_label, ctrl_label, "Sub-Threshold (p > 1e-5)")
)

# Colors matching cohesive slate/navy palette
pal <- c(
  setNames("#1e3a8a", mdd_label),
  setNames("#94a3b8", ctrl_label),
  "Sub-Threshold (p > 1e-5)" = "#cbd5e1"
)

# Active legend breaks excluding sub-threshold
legend_breaks <- c(mdd_label, ctrl_label)

# --- 3. Single-Panel Volcano Plot ---
p_volcano <- ggplot(eqtls, aes(x = ba9_nes, y = neg_log10p)) +
  # Confidence / Threshold Reference Lines
  geom_vline(xintercept = 0, color = "#94a3b8", linewidth = 0.7, linetype = "solid") +
  geom_vline(xintercept = c(-0.3, 0.3), color = "#cbd5e1", linewidth = 0.7, linetype = "dashed") +
  geom_hline(yintercept = 5.0, color = "#ea580c", linewidth = 0.9, linetype = "dashed") +
  
  # Scatter points (sub-threshold points are drawn with lighter opacity)
  geom_point(aes(color = group), size = 5.2, alpha = ifelse(eqtls$neg_log10p < 5.0, 0.45, 0.90)) +
  
  scale_color_manual(values = pal, breaks = legend_breaks) +
  scale_y_continuous(
    breaks = seq(0, 20, 5),
    limits = c(3.5, 23.5)
  ) +
  labs(
    title = "GTEx Frontal Cortex (BA9) eQTL Colocalization",
    subtitle = "Dashed line indicates threshold (p < 1e-5); Cohort enrichment defined as >= 7 carriers (<= 2 in opposing cohort)",
    x = "Normalized Effect Size (NES) in Frontal Cortex (BA9)",
    y = "-log10 (eQTL p-value)",
    color = NULL
  ) +
  theme_minimal(base_family = "sans") +
  theme(
    plot.title = element_text(size = 18.5, face = "bold", color = "#0f172a", margin = margin(b = 4)),
    plot.subtitle = element_text(size = 12.0, face = "italic", color = "#64748b", margin = margin(b = 10)),
    panel.background = element_rect(fill = "#fafbfc", color = "#cbd5e1", linewidth = 0.8),
    panel.grid.major = element_line(color = "#e2e8f0", linewidth = 0.5, linetype = "dotted"),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 17, face = "bold", color = "#334155"),
    axis.text = element_text(size = 14.5, color = "#475569"),
    legend.position = "top",
    legend.text = element_text(size = 14.0, face = "bold", color = "#334155"),
    legend.margin = margin(b = 6),
    plot.margin = margin(t = 14, r = 20, b = 12, l = 15)
  )

# --- 4. Export Single-Panel Plot ---
ggsave(out_png, plot = p_volcano, width = 9.2, height = 6.4, dpi = 300, bg = "white")
ggsave(out_pdf, plot = p_volcano, width = 9.2, height = 6.4, dpi = 300, bg = "white")
ggsave(out_png_utils, plot = p_volcano, width = 9.2, height = 6.4, dpi = 300, bg = "white")

message("[SUCCESS] GTEx volcano plot generated with dynamic sig counts in legend:")
message("  PNG: ", out_png)
message("  PDF: ", out_pdf)
