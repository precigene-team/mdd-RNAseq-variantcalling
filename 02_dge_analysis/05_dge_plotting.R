# [[02_dge_analysis/05_dge_plotting.R]]
# Goal: Load multi-factor DGE and relative methylation results to generate Volcano and MA plots.
#
# Pipeline steps:
#   1. Load shrunken results RDS objects and raw results RDS objects for each contrast
#   2. Construct processed dataframes with pvalue, padj, and gene names
#   3. Process and format Volcano plot data (identify significant, top genes)
#   4. Render Volcano plot using ggplot2 and ggrepel, and save to results/02_dge_analysis/figures/ prefixed with 05_
#   5. Render MA plot using ggplot2 (log10 baseMean vs log2FC), and save to results/02_dge_analysis/figures/ prefixed with 05_

suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(ggrepel)
})

# Dynamically resolve project root
get_project_root <- function() {
  curr <- getwd()
  while (curr != "/" && !file.exists(file.path(curr, "run"))) {
    curr <- dirname(curr)
  }
  return(curr)
}
project_root <- get_project_root()

# Directories
rdata_dir <- file.path(project_root, "_data/RData/02_dge_analysis")
figures_dir <- file.path(project_root, "results/02_dge_analysis/figures")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

# Helper function to plot a single contrast
plot_contrast_results <- function(rds_prefix, volcano_name, ma_name, plot_title) {
  shrunk_file <- file.path(rdata_dir, paste0("04_", rds_prefix, "_res_shrunk.rds"))
  raw_file <- file.path(rdata_dir, paste0("04_", rds_prefix, "_res_raw.rds"))
  
  if (!file.exists(shrunk_file) || !file.exists(raw_file)) {
    message("  [SKIP] RDS files not found for contrast: ", rds_prefix)
    return(NULL)
  }
  
  message("  [PLOT] Rendering plots for: ", rds_prefix)
  
  res <- readRDS(shrunk_file)
  res_raw <- readRDS(raw_file)
  
  res_df <- as.data.frame(res)
  res_df$stat <- res_raw$stat
  res_df$pvalue <- res_raw$pvalue
  res_df$padj <- res_raw$padj
  res_df$gene_name <- rownames(res_df)
  
  # LFC threshold matching paper: FC > 1.2 => log2(1.2) ≈ 0.263
  LFC_THRESH <- log2(1.2)

  # 1. Volcano Plot
  volcano_df <- res_df %>%
    mutate(
      Regulation = case_when(
        !is.na(pvalue) & pvalue < 0.05 & log2FoldChange >= LFC_THRESH ~ "Upregulated",
        !is.na(pvalue) & pvalue < 0.05 & log2FoldChange <= -LFC_THRESH ~ "Downregulated",
        TRUE ~ "Non-significant"
      ),
      Regulation = factor(Regulation, levels = c("Downregulated", "Non-significant", "Upregulated")),
      negLog10P = -log10(pvalue)
    )
  
  # Select top 15 significant genes for labelling
  top_labels <- volcano_df %>%
    filter(Regulation != "Non-significant") %>%
    arrange(pvalue) %>%
    head(15)
  
  volcano_plot <- ggplot(volcano_df, aes(x = log2FoldChange, y = negLog10P)) +
    geom_point(aes(color = Regulation), alpha = 0.6, size = 1.5) +
    geom_vline(xintercept = c(-LFC_THRESH, LFC_THRESH), linetype = "dashed", color = "grey50") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
    geom_text_repel(
      data = top_labels,
      aes(label = gene_name),
      size = 3.5,
      max.overlaps = 20,
      fontface = "italic"
    ) +
    scale_color_manual(
      values = c("Downregulated" = "#1f78b4", "Non-significant" = "grey70", "Upregulated" = "#e31a1c"),
      labels = c("Downregulated" = "Downregulated", "Non-significant" = "Non-significant", "Upregulated" = "Upregulated")
    ) +
    theme_bw(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.title = element_blank(),
      legend.position = "bottom"
    ) +
    labs(
      title = paste0(plot_title, " - Volcano Plot"),
      x = "log2 Fold Change",
      y = "-log10 P-value"
    )
  
  volcano_file <- file.path(figures_dir, volcano_name)
  ggsave(
    filename = volcano_file,
    plot = volcano_plot,
    width = 8,
    height = 6.5,
    dpi = 300,
    bg = "white"
  )
  
  # 2. MA Plot (using ggplot2)
  ma_df <- res_df %>%
    mutate(
      Regulation = case_when(
        !is.na(pvalue) & pvalue < 0.05 & log2FoldChange >= LFC_THRESH ~ "Upregulated",
        !is.na(pvalue) & pvalue < 0.05 & log2FoldChange <= -LFC_THRESH ~ "Downregulated",
        TRUE ~ "Non-significant"
      ),
      Regulation = factor(Regulation, levels = c("Downregulated", "Non-significant", "Upregulated"))
    )
  
  ma_plot <- ggplot(ma_df, aes(x = baseMean, y = log2FoldChange)) +
    geom_point(aes(color = Regulation), alpha = 0.5, size = 1.2) +
    scale_x_log10(
      breaks = 10^(0:6),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    geom_hline(yintercept = 0, color = "#1f78b4", linetype = "dashed", linewidth = 0.8) +
    geom_hline(yintercept = c(-LFC_THRESH, LFC_THRESH), color = "grey50", linetype = "dotted") +
    scale_color_manual(
      values = c("Downregulated" = "#1f78b4", "Non-significant" = "grey70", "Upregulated" = "#e31a1c"),
      labels = c("Downregulated" = "Downregulated", "Non-significant" = "Non-significant", "Upregulated" = "Upregulated")
    ) +
    theme_bw(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.title = element_blank(),
      legend.position = "bottom"
    ) +
    labs(
      title = paste0(plot_title, " - MA Plot"),
      x = "Mean of Normalized Counts (log10 scale)",
      y = "log2 Fold Change"
    )
  
  ma_file <- file.path(figures_dir, ma_name)
  ggsave(
    filename = ma_file,
    plot = ma_plot,
    width = 8,
    height = 6.5,
    dpi = 300,
    bg = "white"
  )
}

# Define contrasts and output details
contrasts_list <- list(
  list(
    rds_prefix = "deg_male_input",
    volcano = "05_volcano_male_input.png",
    ma = "05_ma_male_input.png",
    title = "Male RNA-seq DGE (MDD vs Ctrl)"
  ),
  list(
    rds_prefix = "deg_female_input",
    volcano = "05_volcano_female_input.png",
    ma = "05_ma_female_input.png",
    title = "Female RNA-seq DGE (MDD vs Ctrl)"
  ),
  list(
    rds_prefix = "deg_male_ip",
    volcano = "05_volcano_male_ip.png",
    ma = "05_ma_male_ip.png",
    title = "Male IP DGE (MDD vs Ctrl)"
  ),
  list(
    rds_prefix = "deg_female_ip",
    volcano = "05_volcano_female_ip.png",
    ma = "05_ma_female_ip.png",
    title = "Female IP DGE (MDD vs Ctrl)"
  ),
  list(
    rds_prefix = "meth_male_interaction",
    volcano = "05_volcano_male_meth.png",
    ma = "05_ma_male_meth.png",
    title = "Male Relative Methylation (Interaction)"
  ),
  list(
    rds_prefix = "meth_female_interaction",
    volcano = "05_volcano_female_meth.png",
    ma = "05_ma_female_meth.png",
    title = "Female Relative Methylation (Interaction)"
  ),
  list(
    rds_prefix = "deg_sex_interaction_input",
    volcano = "05_volcano_sex_interaction_input.png",
    ma = "05_ma_sex_interaction_input.png",
    title = "Sex-Specific MDD Response\n(Gender × Condition Interaction)"
  ),
  list(
    rds_prefix = "meth_sex_interaction",
    volcano = "05_volcano_sex_interaction_meth.png",
    ma = "05_ma_sex_interaction_meth.png",
    title = "Sex-Specific m⁶A Methylation\n(Gender × Condition × Fraction)"
  )
)

# Render plots
message("\n--- Generating Volcano and MA plots for contrasts ---")
for (contrast in contrasts_list) {
  tryCatch({
    plot_contrast_results(
      rds_prefix = contrast$rds_prefix,
      volcano_name = contrast$volcano,
      ma_name = contrast$ma,
      plot_title = contrast$title
    )
  }, error = function(e) {
    message("  [ERROR] Plotting failed for ", contrast$rds_prefix, ": ", e$message)
  })
}

message("\n05_dge_plotting.R completed successfully.")
