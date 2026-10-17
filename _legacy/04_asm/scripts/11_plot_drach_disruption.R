#!/usr/bin/env Rscript
# [[04_asm/11_plot_drach_disruption.R]]
# Goal: Visualize DRACH motif disruption patterns using 100% stacked bar plots with clean D-R-A-C-H X-axis labels.

# Steps:
#   1. Resolve dynamic project root path.
#   2. Load annotated DRACH disruption table (05_drach_disrupting_snps.csv or legacy 09_drach_disrupting_snps.csv).
#   3. Deduplicate overlapping variants per (chr, pos, ref, alt, drach_name).
#   4. Compute position-wise disruption proportions across 5-mer motif positions (D-R-A-C-H).
#   5. Render publication-ready 100% stacked bar charts (Offset breakdown & Overall proportions).
#   6. Save figures as PDF and high-resolution PNG to results/04_asm/figures/.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
})

# --- 1. Dynamic Pathing ---
get_project_root <- function() {
  curr <- getwd()
  while (curr != "/" && !file.exists(file.path(curr, "run"))) {
    curr = dirname(curr)
  }
  return(curr)
}
project_root <- get_project_root()

# Input & Output Paths
input_csv <- file.path(project_root, "results/04_asm/tables/05_drach_disrupting_snps.csv")
if (!file.exists(input_csv)) {
  stop("Assertion Error: Input CSV not found at ", input_csv)
}

figures_dir <- file.path(project_root, "results/04_asm/figures")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

out_pdf <- file.path(figures_dir, "11_drach_disruption.pdf")
out_png <- file.path(figures_dir, "11_drach_disruption.png")

# --- 2. Load & Deduplicate Data ---
message("[INFO] Loading DRACH disruption table: ", input_csv)
df <- read_csv(input_csv, show_col_types = FALSE)

stopifnot("disruption_status" %in% colnames(df))
stopifnot("motif_offset" %in% colnames(df))

# Deduplicate cohort-level variant x motif overlaps
df_unique <- df %>%
  distinct(snp_chr, snp_pos, ref, alt, drach_name, .keep_all = TRUE)

message("[INFO] Deduplicated total overlaps: ", nrow(df_unique))

# Clean D - R - A - C - H IUPAC position labels (without Pos 0-4 prefixes)
offset_labels <- c(
  "0" = "D\n(A/G/U)",
  "1" = "R\n(A/G)",
  "2" = "A*\n(m6A Site)",
  "3" = "C\n(Cytosine)",
  "4" = "H\n(A/C/U)"
)

df_plot <- df_unique %>%
  filter(motif_offset %in% 0:4) %>%
  mutate(
    pos_label = factor(offset_labels[as.character(motif_offset)], levels = offset_labels),
    status_label = factor(ifelse(disruption_status == "DISRUPTS_DRACH", "Disrupting", "Preserving"), 
                          levels = c("Disrupting", "Preserving"))
  )

# --- 3. Panel A: 100% Stacked Bar Plot across D-R-A-C-H Positions ---
offset_summary <- df_plot %>%
  group_by(pos_label, status_label) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(pos_label) %>%
  mutate(
    total_pos = sum(count),
    pct = count / total_pos * 100
  )

color_palette <- c("Disrupting" = "#d9534f", "Preserving" = "#337ab7")

p1 <- ggplot(offset_summary, aes(x = pos_label, y = pct, fill = status_label)) +
  geom_bar(stat = "identity", position = position_stack(reverse = TRUE), width = 0.55, color = "white", linewidth = 0.4) +
  geom_text(
    aes(label = ifelse(pct >= 4, paste0(round(pct, 1), "%\n(", format(count, big.mark = ","), ")"), "")),
    position = position_stack(vjust = 0.5, reverse = TRUE),
    size = 3.3, fontface = "bold", color = "white"
  ) +
  scale_fill_manual(values = color_palette) +
  scale_y_continuous(labels = function(x) paste0(x, "%"), limits = c(0, 100.5), expand = c(0, 0)) +
  labs(
    title = "A. DRACH Motif Disruption Proportion by Motif Position",
    subtitle = "Position A* (central m6A Adenine) is 100% disrupted by any overlapping SNP",
    x = "DRACH Consensus Motif Positions",
    y = "Percentage of Overlapping SNPs (%)",
    fill = "Effect on Motif:"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(color = "grey30", size = 9.5, margin = margin(b = 8)),
    axis.text.x = element_text(angle = 0, hjust = 0.5, face = "bold", size = 10),
    legend.position = "top",
    legend.title = element_text(face = "bold", size = 10),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

# --- 4. Panel B: Overall Cohort Proportion (100% Stacked Bar) ---
overall_summary <- df_plot %>%
  group_by(status_label) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(
    pct = count / sum(count) * 100,
    x_group = "Cohort Total"
  )

p2 <- ggplot(overall_summary, aes(x = x_group, y = pct, fill = status_label)) +
  geom_bar(stat = "identity", position = position_stack(reverse = TRUE), width = 0.42, color = "white", linewidth = 0.4) +
  geom_text(
    aes(label = paste0(status_label, "\n", round(pct, 1), "%\n(", format(count, big.mark = ","), ")")),
    position = position_stack(vjust = 0.5, reverse = TRUE),
    size = 3.4, fontface = "bold", color = "white"
  ) +
  scale_fill_manual(values = color_palette) +
  scale_y_continuous(labels = function(x) paste0(x, "%"), limits = c(0, 100.5), expand = c(0, 0)) +
  labs(
    title = "B. Cohort Total",
    subtitle = paste0("N = ", format(nrow(df_unique), big.mark = ",")),
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(color = "grey30", size = 9.5, margin = margin(b = 8)),
    legend.position = "none",
    axis.text.x = element_text(face = "bold", size = 10.5),
    axis.text.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

# --- 5. Combine Panels & Save ---
if (requireNamespace("cowplot", quietly = TRUE)) {
  combined_plot <- cowplot::plot_grid(p1, p2, rel_widths = c(2.4, 1), align = "h", axis = "tb")
} else if (requireNamespace("patchwork", quietly = TRUE)) {
  combined_plot <- patchwork::wrap_plots(p1, p2, widths = c(2.4, 1))
} else {
  library(gridExtra)
  combined_plot <- gridExtra::marrangeGrob(list(p1, p2), ncol = 2, nrow = 1, widths = c(2.4, 1))
}

message("[INFO] Saving figure to: ", out_pdf)
ggsave(out_pdf, plot = combined_plot, width = 11, height = 5.5, dpi = 300)

message("[INFO] Saving figure to: ", out_png)
ggsave(out_png, plot = combined_plot, width = 11, height = 5.5, dpi = 300)

message("[DONE] DRACH disruption visualization with D-R-A-C-H labels complete.")
