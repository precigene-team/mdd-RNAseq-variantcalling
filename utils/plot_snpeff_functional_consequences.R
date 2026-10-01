#!/usr/bin/env Rscript
# [[utils/plot_snpeff_functional_consequences.R]]
# Goal: 2 Vertical bar plots for slides: (1) SnpEff Impact Tiers, (2) High-Impact Breakdown. Export to scratch/.

# Steps:
#   1. Dynamic project root resolution.
#   2. Ingest 07_functional_annotation.csv.
#   3. Plot 1: Vertical bars for impact tiers (MODIFIER, LOW, MODERATE, HIGH).
#   4. Plot 2: Vertical bars for high-impact consequences breakdown (N = 139).
#   5. Export both standalone figures directly to scratch/.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(scales)
})

# --- Step 1: Dynamic Pathing Setup ---
get_project_root <- function() {
  curr <- getwd()
  while (curr != "/" && !file.exists(file.path(curr, "run")) && !file.exists(file.path(curr, "CONTEXT.md"))) {
    curr <- dirname(curr)
  }
  return(curr)
}

project_root <- Sys.getenv("PROJECT_ROOT", unset = get_project_root())
tables_dir <- file.path(project_root, "results", "03_variant_calling", "tables")
figures_dir <- file.path(project_root, "results", "utils", "figures")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

in_csv <- file.path(tables_dir, "07_functional_annotation.csv")
if (!file.exists(in_csv)) {
  stop("[ERROR] Input annotation table not found: ", in_csv)
}

# --- Step 2: Ingest Data ---
cat("[INFO] Ingesting SnpEff annotations from:", in_csv, "\n")
df <- read_csv(in_csv, show_col_types = FALSE, progress = FALSE)
total_snps <- nrow(df)

# Base Slide Theme
theme_vertical_slide <- theme_minimal(base_family = "sans", base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 15, color = "#0f172a", margin = margin(b = 4)),
    plot.subtitle = element_text(size = 11.5, color = "#475569", margin = margin(b = 14)),
    axis.title.x = element_blank(),
    axis.title.y = element_text(face = "bold", size = 11.5, color = "#334155", margin = margin(r = 8)),
    axis.text.x = element_text(face = "bold", size = 12, color = "#0f172a"),
    axis.text.y = element_text(size = 10.5, color = "#475569"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "#e2e8f0", linewidth = 0.5, linetype = "dashed"),
    plot.margin = margin(t = 15, r = 20, b = 15, l = 20)
  )

# --- Step 3: Plot 1 - Impact Levels Only (Vertical Bars) ---
impact_summary <- df %>%
  count(impact) %>%
  mutate(
    impact = factor(impact, levels = c("MODIFIER", "LOW", "MODERATE", "HIGH")),
    pct = (n / total_snps) * 100,
    label_text = sprintf("%s\n(%0.1f%%)", comma(n), pct)
  )

tier_colors <- c(
  "MODIFIER" = "#2563eb", # Blue
  "LOW"      = "#16a34a", # Green
  "MODERATE" = "#ea580c", # Orange
  "HIGH"     = "#dc2626"  # Red
)

p1_impact <- ggplot(impact_summary, aes(x = impact, y = n, fill = impact)) +
  geom_col(width = 0.55, color = "#0f172a", linewidth = 0.3) +
  geom_text(
    aes(label = label_text),
    vjust = -0.3,
    size = 4.2,
    fontface = "bold",
    color = "#0f172a",
    lineheight = 0.95
  ) +
  scale_fill_manual(values = tier_colors, guide = "none") +
  scale_y_continuous(
    labels = comma,
    limits = c(0, 60000),
    breaks = seq(0, 50000, 10000),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    title = "SnpEff Functional Impact Distribution",
    subtitle = sprintf("Cohort-wide recurrent transcriptome SNPs (N = %s, >=8 subjects)", comma(total_snps)),
    y = "Variant Count"
  ) +
  theme_vertical_slide

# --- Step 4: Plot 2 - High-Impact Breakdown (Vertical Bars) ---
clean_high_term <- function(c) {
  c <- sub("&intron_variant", " (Intronic)", c)
  c <- sub("&splice_region_variant", " + Splice", c)
  c <- gsub("_variant", "", c)
  c <- gsub("_", " ", c)
  c <- sub("stop gained", "Stop Gained", c)
  c <- sub("splice donor", "Splice Donor", c)
  c <- sub("splice acceptor", "Splice Acceptor", c)
  c <- sub("start lost", "Start Lost", c)
  c <- sub("stop lost", "Stop Lost", c)
  return(c)
}

high_df <- df %>%
  filter(impact == "HIGH") %>%
  mutate(display_cons = clean_high_term(consequence))

total_high <- nrow(high_df)

high_summary <- high_df %>%
  count(display_cons) %>%
  arrange(desc(n)) %>%
  mutate(
    display_cons = factor(display_cons, levels = display_cons),
    pct = (n / total_high) * 100,
    label_text = sprintf("%d\n(%0.1f%%)", n, pct)
  )

p2_high <- ggplot(high_summary, aes(x = display_cons, y = n)) +
  geom_col(width = 0.55, fill = "#dc2626", color = "#0f172a", linewidth = 0.3) +
  geom_text(
    aes(label = label_text),
    vjust = -0.3,
    size = 4.0,
    fontface = "bold",
    color = "#0f172a",
    lineheight = 0.95
  ) +
  scale_y_continuous(
    limits = c(0, 78),
    breaks = seq(0, 70, 10),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    title = "High-Impact Variant Breakdown",
    subtitle = sprintf("Deleterious protein-truncating & essential splice variants (N = %d across %d genes)", 
                       total_high, n_distinct(high_df$gene)),
    y = "Variant Count"
  ) +
  theme_vertical_slide +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1, vjust = 1, face = "bold", size = 11, color = "#0f172a")
  )

# --- Step 5: Export to results/utils/figures/ ---
out_p1_png <- file.path(figures_dir, "snpeff_impact_levels.png")
out_p1_pdf <- file.path(figures_dir, "snpeff_impact_levels.pdf")
out_p2_png <- file.path(figures_dir, "snpeff_high_impact_breakdown.png")
out_p2_pdf <- file.path(figures_dir, "snpeff_high_impact_breakdown.pdf")

cat("[INFO] Exporting to results/utils/figures/ ...\n")
ggsave(out_p1_png, plot = p1_impact, width = 7.5, height = 5.8, dpi = 300, bg = "white")
ggsave(out_p1_pdf, plot = p1_impact, width = 7.5, height = 5.8, dpi = 300, bg = "white")

ggsave(out_p2_png, plot = p2_high, width = 8.5, height = 6.2, dpi = 300, bg = "white")
ggsave(out_p2_pdf, plot = p2_high, width = 8.5, height = 6.2, dpi = 300, bg = "white")

cat("[SUCCESS] Successfully generated separate vertical plots in results/utils/figures/:\n")
cat("  - Plot 1 (Impact Levels): ", out_p1_png, "\n")
cat("  - Plot 2 (High Impact):   ", out_p2_png, "\n")
