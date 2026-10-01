#!/usr/bin/env Rscript
# [[utils/plot_gtex_eqtl_summary.R]]
# Goal: Ultra-clean, minimal single-panel volcano plot of GTEx Brain BA9 eQTL matches with enrichment counts in legend and criteria in subtitle.

# Steps:
#   1. Load GTEx eQTL match results from results/tables/07_gtex_brain_eqtl_matches.csv
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
utils_fig_dir <- file.path(project_root, "results", "utils", "figures")
dir.create(utils_fig_dir, showWarnings = FALSE, recursive = TRUE)

out_png_utils <- file.path(utils_fig_dir, "gtex_brain_eqtl_summary.png")
out_pdf_utils <- file.path(utils_fig_dir, "gtex_brain_eqtl_summary.pdf")

# --- 2. Load and Prepare GTEx Match Data ---
csv_path <- file.path(project_root, "results", "03_variant_calling", "tables", "08_gtex_brain_eqtl_matches.csv")
if (!file.exists(csv_path)) {
  csv_path <- file.path(project_root, "results", "03_variant_calling", "tables", "07_gtex_brain_eqtl_matches.csv")
}
if (!file.exists(csv_path)) {
  stop("ERROR: Input GTEx table not found: ", csv_path)
}

df <- read.csv(csv_path, stringsAsFactors = FALSE)
# Deduplicate by variant_id to retain exact 29 candidate loci
df <- df[!duplicated(df$variant_id), ]

total_screened <- nrow(df)
hit_idx <- which(df$is_gtex_ba9_eqtl == "YES")
non_hit_idx <- which(df$is_gtex_ba9_eqtl != "YES")

n_hits <- length(hit_idx)
n_non_hits <- length(non_hit_idx)

# Hits data
hits_df <- df[hit_idx, ]
hits_df$ba9_eqtl_pvalue <- as.numeric(hits_df$ba9_eqtl_pvalue)
hits_df$ba9_nes <- as.numeric(hits_df$ba9_nes)
hits_df$neg_log10p <- -log10(hits_df$ba9_eqtl_pvalue)

# Order hits by effect size
hits_df <- hits_df[order(hits_df$ba9_nes), ]

# Format labels
hits_df$locus_label <- factor(
  paste0(hits_df$gene, " (", hits_df$n_mdd, ":", hits_df$n_ctrl, ")"),
  levels = paste0(hits_df$gene, " (", hits_df$n_mdd, ":", hits_df$n_ctrl, ")")
)

hits_df$dir <- ifelse(hits_df$ba9_nes < 0, "Repressive (NES < 0)", "Activating (NES > 0)")
hits_df$dir <- factor(hits_df$dir, levels = c("Repressive (NES < 0)", "Activating (NES > 0)"))
hits_df$stat_label <- sprintf("NES: %+.2f | p: %.1e", hits_df$ba9_nes, hits_df$ba9_eqtl_pvalue)

# --- 3. Panel A: Candidate Screening Funnel Bar ---
df_funnel <- data.frame(
  Category = factor(
    c("Validated BA9 eQTL", "Sub-Threshold in BA9"), 
    levels = c("Sub-Threshold in BA9", "Validated BA9 eQTL")
  ),
  Count = c(n_non_hits, n_hits),
  Label = c(
    paste0("Sub-Threshold\nN = ", n_non_hits, " (", sprintf("%.1f", 100 * n_non_hits / total_screened), "%)"),
    paste0("Validated Hits\nN = ", n_hits, " (", sprintf("%.1f", 100 * n_hits / total_screened), "%)")
  )
)

p_funnel <- ggplot(df_funnel, aes(x = "Candidates", y = Count, fill = Category)) +
  geom_bar(stat = "identity", width = 0.55, color = "#0f172a", linewidth = 0.6) +
  geom_text(
    aes(label = Label),
    position = position_stack(vjust = 0.5),
    color = c("#1e293b", "white"),
    fontface = "bold",
    size = 4.2
  ) +
  scale_fill_manual(values = c("Validated BA9 eQTL" = "#1e3a8a", "Sub-Threshold in BA9" = "#cbd5e1")) +
  scale_y_continuous(breaks = seq(0, 30, by = 5), limits = c(0, 30.5)) +
  labs(
    title = "A. Candidate Screening Yield",
    subtitle = paste0("Total Screened: N = ", total_screened, " Loci\n(n_MDD >= 8, n_CTRL <= 2)"),
    x = NULL,
    y = "Candidate Loci (N)"
  ) +
  theme_minimal(base_family = "sans") +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 14.5, face = "bold", color = "#0f172a", margin = margin(b = 4)),
    plot.subtitle = element_text(size = 11, face = "italic", color = "#64748b", margin = margin(b = 10)),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(size = 12.5, face = "bold", color = "#334155"),
    axis.text.y = element_text(size = 11.5, color = "#475569"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "#fafbfc", color = "#cbd5e1", linewidth = 0.8),
    plot.margin = margin(t = 12, r = 10, b = 12, l = 12)
  )

# --- 4. Panel B: Validated Regulatory Forest Plot ---
p_forest <- ggplot(hits_df, aes(x = ba9_nes, y = locus_label)) +
  geom_vline(xintercept = 0, color = "#94a3b8", linewidth = 0.8, linetype = "solid") +
  geom_vline(xintercept = c(-0.3, 0.3), color = "#cbd5e1", linewidth = 0.6, linetype = "dashed") +
  geom_segment(
    aes(x = 0, xend = ba9_nes, y = locus_label, yend = locus_label, color = dir),
    linewidth = 1.2, alpha = 0.85
  ) +
  geom_point(
    aes(color = dir, size = neg_log10p),
    alpha = 0.95
  ) +
  geom_text(
    aes(
      label = stat_label,
      hjust = ifelse(ba9_nes < 0, 1.15, -0.15)
    ),
    size = 3.6,
    fontface = "bold",
    color = "#1e293b"
  ) +
  scale_color_manual(values = c("Repressive (NES < 0)" = "#1e3a8a", "Activating (NES > 0)" = "#ea580c")) +
  scale_size_continuous(name = expression(-log[10](italic(p))), range = c(4.5, 9.0)) +
  scale_x_continuous(
    breaks = seq(-1.2, 1.3, by = 0.3),
    limits = c(-1.25, 1.35)
  ) +
  labs(
    title = "B. Validated Frontal Cortex (BA9) Regulatory Impact",
    subtitle = "Point size represents significance (-log10 p-value) | Stems indicate effect size (NES)",
    x = "Normalized Effect Size (NES) in Frontal Cortex (BA9)",
    y = "Candidate Locus (MDD : CTRL Ratio)",
    color = "Regulatory Direction"
  ) +
  theme_minimal(base_family = "sans") +
  theme(
    plot.title = element_text(size = 14.5, face = "bold", color = "#0f172a", margin = margin(b = 4)),
    plot.subtitle = element_text(size = 11, face = "italic", color = "#64748b", margin = margin(b = 10)),
    panel.background = element_rect(fill = "#fafbfc", color = "#cbd5e1", linewidth = 0.8),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(color = "#e2e8f0", linewidth = 0.5, linetype = "dotted"),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 12.5, face = "bold", color = "#334155"),
    axis.text.y = element_text(size = 11.5, face = "bold", color = "#1e293b"),
    axis.text.x = element_text(size = 11, color = "#475569"),
    legend.position = "right",
    legend.title = element_text(size = 10.5, face = "bold", color = "#334155"),
    legend.text = element_text(size = 10, color = "#475569"),
    plot.margin = margin(t = 12, r = 16, b = 12, l = 10)
  )

# --- 5. Export Dual-Panel Figure via Base Grid ---
suppressPackageStartupMessages(library(grid))

render_combined <- function() {
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(1, 2, widths = unit(c(3.6, 8.9), "null"))))
  print(p_funnel, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(p_forest, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
}

pdf(out_pdf_utils, width = 12.5, height = 6.2, bg = "white")
render_combined()
dev.off()

png(out_png_utils, width = 12.5, height = 6.2, units = "in", res = 300, bg = "white")
render_combined()
dev.off()

message("[SUCCESS] Concept 1 Dual-Panel figure generated:")
message("  PNG: ", out_png_utils)
message("  PDF: ", out_pdf_utils)
