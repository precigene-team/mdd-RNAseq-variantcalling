#!/usr/bin/env Rscript
# [[03_variant_calling/15_plot_sig_asm_drach.R]]
# Goal: Visualize DRACH motif disruption across the 5 FDR-significant ASm6A candidate loci using biological Allele Skew magnitude.

# Steps:
#   1. Dynamic project root resolution.
#   2. Load 5 FDR-significant ASm6A loci summary table (03_asm_5gene_summary_table.csv / 12_asm_primary_significant.csv).
#   3. Load annotated DRACH disruption table (09_drach_disrupting_snps.csv).
#   4. Map overlapping heterozygous SNPs and DRACH motif disruption status per locus.
#   5. Render publication-ready ggplot figure using Baseline Allele Skew Magnitude (|Delta Skew| %) on y-axis for informative variation.
#   6. Save figures as PDF and high-resolution PNG to results/03_variant_calling/figures/.

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
summary_table_file <- file.path(project_root, "results/tables/03_asm_5gene_summary_table.csv")
sig_table_file     <- file.path(project_root, "results/03_variant_calling/tables/12_asm_primary_significant.csv")
drach_table_file   <- file.path(project_root, "results/03_variant_calling/tables/09_drach_disrupting_snps.csv")
if (!file.exists(drach_table_file)) {
  drach_table_file <- file.path(project_root, "results/tables/09_drach_disrupting_snps.csv")
}

stopifnot(file.exists(sig_table_file))
stopifnot(file.exists(drach_table_file))

figures_dir <- file.path(project_root, "results/03_variant_calling/figures")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

out_pdf <- file.path(figures_dir, "15_sig_asm_drach_loci.pdf")
out_png <- file.path(figures_dir, "15_sig_asm_drach_loci.png")

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
      mutate(
        gene = g,
        locus_fdr = fdr,
        locus_confidence = conf
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
    has_core_m6a_hit = sapply(gene, function(g) {
      if (g %in% names(loci_matches)) {
        any(loci_matches[[g]]$motif_offset == 2 & loci_matches[[g]]$disruption_status == "DISRUPTS_DRACH")
      } else FALSE
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
    subtitle = "MTCH1, BIN3, and MGAT5 exhibit strong allele-specific methylation enrichment (FDR < 0.02)",
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
matched_plot_df <- matched_df %>%
  mutate(
    variant_label = paste0("Chr", snp_chr, ":", snp_pos, " (", ref, ">", alt, ")"),
    offset_desc = case_when(
      motif_offset == 2 ~ "Pos 2 (Core m6A Site)",
      TRUE ~ paste0("Pos ", motif_offset)
    ),
    status_clean = ifelse(disruption_status == "DISRUPTS_DRACH", "Disrupting SNP", "Preserving SNP")
  )

p2 <- ggplot(matched_plot_df, aes(x = factor(snp_pos), y = gene, fill = status_clean)) +
  geom_tile(color = "white", linewidth = 0.8, height = 0.6) +
  geom_text(aes(label = paste0(ref, ">", alt, "\n[Off ", motif_offset, "]")), size = 3, fontface = "bold", color = "white") +
  scale_fill_manual(values = c("Disrupting SNP" = "#d9534f", "Preserving SNP" = "#4682b4")) +
  labs(
    title = "B. Overlapping Heterozygous Variant Fine Mapping",
    subtitle = "BIN3 carries a direct core m6A site disruption (Ref A > Alt G at offset 2)",
    x = "Genomic SNP Position",
    y = "Gene Locus",
    fill = "Variant Effect"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(color = "grey30", size = 10, margin = margin(b = 10)),
    axis.text.y = element_text(face = "bold", size = 11),
    legend.position = "right"
  )

# --- 5. Combine & Export ---
if (requireNamespace("cowplot", quietly = TRUE)) {
  combined <- cowplot::plot_grid(p1, p2, ncol = 1, rel_heights = c(1.2, 0.8))
} else {
  library(gridExtra)
  combined <- gridExtra::marrangeGrob(list(p1, p2), ncol = 1, nrow = 2, heights = c(1.2, 0.8))
}

message("[INFO] Exporting figure to: ", out_pdf)
ggsave(out_pdf, plot = combined, width = 10, height = 8.5, dpi = 300)

message("[INFO] Exporting figure to: ", out_png)
ggsave(out_png, plot = combined, width = 10, height = 8.5, dpi = 300)

message("[DONE] 15_plot_sig_asm_drach.R executed successfully.")
