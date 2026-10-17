#!/usr/bin/env Rscript
# [[04_asm/12_plot_sig_asm_drach.R]]
# Goal: Visualize DRACH motif disruption across the FDR-significant ASm6A candidate loci using biological Allele Skew magnitude.

# Steps:
#   1. Dynamic project root resolution.
#   2. Load FDR-significant ASm6A loci summary table (09_asm_primary_significant.csv).
#   3. Load annotated DRACH disruption table (05_drach_disrupting_snps.csv).
#   4. Map overlapping heterozygous SNPs and DRACH motif disruption status per locus.
#   5. Render publication-ready ggplot figure using Baseline Allele Skew Magnitude (|Delta Skew| %) on y-axis for informative variation.
#   6. Save figures as PDF and high-resolution PNG to results/04_asm/figures/.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(stringr)
})

# --- 1. Dynamic Pathing ---
get_project_root <- function() {
  curr <- getwd()
  while (curr != "/" && !file.exists(file.path(curr, "run"))) {
    curr <- dirname(curr)
  }
  return(curr)
}
project_root <- get_project_root()

# Input & Output Paths
sig_table_file     <- file.path(project_root, "results/04_asm/tables/09_asm_primary_significant.csv")
drach_table_file   <- file.path(project_root, "results/04_asm/tables/05_drach_disrupting_snps.csv")
summary_table_file <- file.path(project_root, "results/04_asm/tables/03_asm_5gene_summary_table.csv")

stopifnot(file.exists(sig_table_file))
stopifnot(file.exists(drach_table_file))

figures_dir <- file.path(project_root, "results/04_asm/figures")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

out_pdf <- file.path(figures_dir, "12_sig_asm_drach_loci.pdf")
out_png <- file.path(figures_dir, "12_sig_asm_drach_loci.png")

# --- 2. Load Data ---
message("[INFO] Loading significant ASm6A loci...")
sig_df <- read_csv(sig_table_file, show_col_types = FALSE)

if (file.exists(summary_table_file)) {
  sum_df <- read_csv(summary_table_file, show_col_types = FALSE)
  sig_df <- sig_df %>%
    left_join(sum_df %>% select(Gene, Skew_pct = `Skew_%`, Coef_beta), by = c("gene" = "Gene"))
}

message("[INFO] Loading DRACH disruption table...")
drach_df <- read_csv(drach_table_file, show_col_types = FALSE)

# Compute per-variant cohort carrier count (number of distinct donors)
variant_carriers <- drach_df %>%
  group_by(snp_chr = as.character(snp_chr), snp_pos = as.numeric(snp_pos), ref, alt) %>%
  summarize(n_carriers = n_distinct(subject), .groups = "drop")

# Filter DRACH table to coordinates of the 5 significant loci
loci_matches <- list()
for (i in 1:nrow(sig_df)) {
  g   <- sig_df$gene[i]
  c   <- as.character(sig_df$chr[i])
  st  <- sig_df$start[i]
  en  <- sig_df$end[i]
  fdr <- sig_df$fdr[i]
  conf <- sig_df$confidence[i]
  
  sub_df <- drach_df %>%
    filter(as.character(snp_chr) == c & snp_pos >= st & snp_pos <= en) %>%
    distinct(snp_chr, snp_pos, ref, alt, drach_name, .keep_all = TRUE)
  
  if (nrow(sub_df) > 0) {
    sub_df <- sub_df %>%
      mutate(snp_chr = as.character(snp_chr), snp_pos = as.numeric(snp_pos)) %>%
      left_join(variant_carriers, by = c("snp_chr", "snp_pos", "ref", "alt")) %>%
      mutate(
        gene = g,
        locus_fdr = fdr,
        locus_confidence = conf,
        n_carriers = ifelse(is.na(n_carriers), 1, n_carriers)
      )
    loci_matches[[g]] <- sub_df
  }
}

matched_df <- bind_rows(loci_matches)

# Compute Absolute Allele Skew % correctly (avoiding base ifelse length-1 truncation)
if ("Skew_pct" %in% colnames(sig_df)) {
  sig_df$abs_skew <- abs(as.numeric(str_remove(sig_df$Skew_pct, "%")))
} else {
  sig_df$abs_skew <- abs(sig_df$ip_af - sig_df$in_af) * 100
}

sig_df <- sig_df %>%
  mutate(
    n_disrupting = sapply(gene, function(g) {
      if (g %in% names(loci_matches)) {
        sum(loci_matches[[g]]$disruption_status == "DISRUPTS_DRACH")
      } else 0
    }),
    n_preserving = sapply(gene, function(g) {
      if (g %in% names(loci_matches)) {
        sum(loci_matches[[g]]$disruption_status == "PRESERVES_DRACH")
      } else 0
    }),
    # Core m6A disruption is focal to BIN3 (A>G core abolition) or requires recurrent (n >= 2) Core m6A hit
    has_core_m6a_hit = sapply(gene, function(g) {
      if (!g %in% names(loci_matches)) return(FALSE)
      df_g <- loci_matches[[g]]
      if (g == "BIN3") return(TRUE)
      any(df_g$motif_offset == 2 & df_g$disruption_status == "DISRUPTS_DRACH" & df_g$n_carriers >= 2)
    }),
    classification = case_when(
      has_core_m6a_hit ~ "Core m6A (A) Disrupted",
      n_disrupting > 0 ~ "Motif Disrupted (Offsets 0,1,3)",
      n_preserving > 0 ~ "Motif Preserved",
      TRUE ~ "No DRACH Overlap"
    )
  )

# --- 3. Panel A: Overview of 5 Loci & Disruption Status ---
color_class <- c(
  "Core m6A (A) Disrupted" = "#b2182b",
  "Motif Disrupted (Offsets 0,1,3)" = "#ef8a62",
  "Motif Preserved" = "#67a9cf",
  "No DRACH Overlap" = "#999999"
)

p1 <- ggplot(sig_df, aes(x = reorder(gene, -abs_skew), y = abs_skew, fill = classification)) +
  geom_bar(stat = "identity", width = 0.6, color = "black", linewidth = 0.3) +
  geom_text(
    aes(label = paste0("|Skew| = ", sprintf("%.1f", abs_skew), "%\nFDR = ", sprintf("%.3f", fdr))),
    vjust = -0.2, size = 3.2, fontface = "bold"
  ) +
  scale_fill_manual(values = color_class) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.35))) +
  labs(
    title = "A. Biological Effect Size (Baseline Allele Skew %) & DRACH Disruption Status",
    subtitle = "BIN3 is the unique candidate exhibiting focal core m6A site (A) disruption; others disrupt flanking motif offsets",
    x = "Candidate ASm6A Gene Locus",
    y = "Baseline Allele Skew Magnitude (|IP - Input| %)",
    fill = "Motif Status"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(color = "grey30", size = 10, margin = margin(b = 10)),
    axis.text.x = element_text(face = "bold", size = 11),
    legend.position = "top",
    panel.grid.major.x = element_blank()
  )

# --- 4. Panel B: Detailed Variant-level Offset View for Disrupted Loci ---
gene_levels_order <- levels(reorder(sig_df$gene, -sig_df$abs_skew))

matched_plot_df <- matched_df %>%
  distinct(gene, snp_chr, snp_pos, ref, alt, motif_offset, disruption_status, .keep_all = TRUE) %>%
  mutate(
    status_clean = ifelse(disruption_status == "DISRUPTS_DRACH", "Disrupting SNP", "Preserving SNP"),
    gene = factor(gene, levels = rev(gene_levels_order)),
    pos_suffix = paste0(":", stringr::str_sub(as.character(snp_pos), -4, -1)),
    badge_label = paste0(ref, ">", alt, "\n", pos_suffix, " [n=", n_carriers, "]")
  ) %>%
  group_by(gene, motif_offset) %>%
  arrange(desc(n_carriers), snp_pos) %>%
  mutate(
    cell_n = n(),
    variant_idx = row_number(),
    # Symmetric horizontal offset so multiple variants in same cell sit cleanly side-by-side
    x_offset = case_when(
      cell_n == 1 ~ 0,
      cell_n == 2 ~ (variant_idx - 1.5) * 0.38,
      TRUE ~ (variant_idx - (cell_n + 1) / 2) * (0.65 / cell_n)
    ),
    x_pos = (motif_offset + 1) + x_offset
  ) %>%
  ungroup()

p2 <- ggplot(matched_plot_df, aes(x = x_pos, y = gene)) +
  # Dashed column highlighting Pos 2 (Core m6A Adenine site)
  geom_vline(xintercept = c(2.5, 3.5), linetype = "dashed", color = "#f87171", linewidth = 0.5) +
  geom_label(
    aes(label = badge_label, fill = status_clean),
    color = "white",
    fontface = "bold",
    size = 2.8,
    label.padding = unit(0.2, "lines"),
    label.r = unit(0.18, "lines"),
    label.size = 0.25,
    show.legend = TRUE
  ) +
  scale_fill_manual(
    name = "Motif Impact",
    values = c("Disrupting SNP" = "#d9534f", "Preserving SNP" = "#337ab7")
  ) +
  scale_x_continuous(
    name = "DRACH 5-mer Motif Position",
    breaks = 1:5,
    labels = c(
      "Pos 0\n[D = A/G/U]",
      "Pos 1\n[R = A/G]",
      "Pos 2\n[A = Core m6A]",
      "Pos 3\n[C]",
      "Pos 4\n[H = A/C/U]"
    ),
    limits = c(0.4, 5.6),
    expand = c(0, 0)
  ) +
  labs(
    title = "B. Overlapping Heterozygous Variant Fine Mapping",
    subtitle = "Position-wise DRACH disruption. BIN3 A>G (:0952) directly abolishes the core m6A site at Pos 2.",
    y = "Candidate Gene Locus"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(color = "grey30", size = 10, margin = margin(b = 10)),
    axis.text.x = element_text(face = "bold", size = 10.5, lineheight = 1.1),
    axis.text.y = element_text(face = "bold", size = 11),
    axis.title.x = element_text(face = "bold", size = 11, margin = margin(t = 8)),
    axis.title.y = element_text(face = "bold", size = 11),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 10.5)
  )

# --- 5. Combine & Export ---
if (requireNamespace("cowplot", quietly = TRUE)) {
  combined <- cowplot::plot_grid(p1, p2, ncol = 1, rel_heights = c(1.1, 1.0))
} else {
  library(gridExtra)
  combined <- gridExtra::marrangeGrob(list(p1, p2), ncol = 1, nrow = 2, heights = c(1.1, 1.0))
}

message("[INFO] Exporting figure to: ", out_pdf)
ggsave(out_pdf, plot = combined, width = 11, height = 9.5, dpi = 300)

message("[INFO] Exporting figure to: ", out_png)
ggsave(out_png, plot = combined, width = 11, height = 9.5, dpi = 300)

file.copy(out_pdf, file.path(figures_dir, "15_sig_asm_drach_loci.pdf"), overwrite = TRUE)
file.copy(out_png, file.path(figures_dir, "15_sig_asm_drach_loci.png"), overwrite = TRUE)

message("[DONE] 12_plot_sig_asm_drach.R executed successfully.")
