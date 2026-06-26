# [[02_dge_analysis/05_dge_plotting.R]]
# Goal: Load shrunken DESeq2 results and generate Volcano and MA plots.
#
# Pipeline steps:
#   1. Load shrunken results RDS objects and raw results RDS objects
#   2. Load metadata to identify comparison groups
#   3. Process and format Volcano plot data (identify significant, top genes)
#   4. Render Volcano plot using ggplot2 and ggrepel, and save to results/figures/05_{DataType}_volcano.png
#   5. Render MA plot using base R plotMA, and save to results/figures/05_{DataType}_MA.png

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
figures_dir <- file.path(project_root, "results/figures")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

# Plot function
plot_dge_results <- function(datatype) {
  message("\n==============================")
  message("Plotting DGE for: ", datatype)
  message("==============================")

  shrunk_file <- file.path(rdata_dir, paste0("04_", datatype, "_res_shrunk.rds"))
  raw_file <- file.path(rdata_dir, paste0("04_", datatype, "_res_raw.rds"))
  dds_file <- file.path(rdata_dir, paste0("02_", datatype, "_dds.rds"))

  if (!file.exists(shrunk_file) || !file.exists(raw_file) || !file.exists(dds_file)) {
    stop("Required DGE results RDS files not found for ", datatype)
  }

  res <- readRDS(shrunk_file)
  res_raw <- readRDS(raw_file)
  dds <- readRDS(dds_file)

  # Check comparison groups
  groups <- levels(colData(dds)$SampleType)

  res_df <- as.data.frame(res)
  res_df$stat <- res_raw$stat
  res_df$pvalue <- res_raw$pvalue
  res_df$padj <- res_raw$padj
  res_df$gene_name <- rownames(res_df)

  # Volcano plot
  volcano_df <- res_df %>%
    mutate(
      Significant = !is.na(pvalue) & pvalue < 0.05 & abs(log2FoldChange) >= 1,
      negLog10P = -log10(pvalue)
    )

  top_labels <- volcano_df %>%
    filter(Significant) %>%
    arrange(pvalue) %>%
    head(20)

  volcano_plot <- ggplot(volcano_df, aes(x = log2FoldChange, y = negLog10P)) +
    geom_point(aes(color = Significant), alpha = 0.7, size = 1.5) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    geom_text_repel(data = top_labels, aes(label = gene_name), size = 3, max.overlaps = 20) +
    scale_color_manual(
      values = c("FALSE" = "grey70", "TRUE" = "red"),
      labels = c("FALSE" = "Non-significant", "TRUE" = "Significant")
    ) +
    theme_bw() +
    labs(
      title = paste0(datatype, " Volcano Plot (", groups[2], " vs ", groups[1], ")"),
      x = "log2 Fold Change",
      y = "-log10 P-value"
    )

  # Save volcano plot to results/figures/ with prefix 05_
  volcano_file <- file.path(figures_dir, paste0("05_", datatype, "_volcano.png"))
  ggsave(
    filename = volcano_file,
    plot = volcano_plot,
    width = 8,
    height = 6,
    dpi = 300
  )
  message("  [OUT] Saved Volcano plot to ", volcano_file)

  # MA plot to results/figures/ with prefix 05_
  ma_file <- file.path(figures_dir, paste0("05_", datatype, "_MA.png"))
  png(
    filename = ma_file,
    width = 2400,
    height = 1800,
    res = 300
  )

  plotMA(
    res,
    ylim = c(-5, 5),
    main = paste0(datatype, " MA Plot (", groups[2], " vs ", groups[1], ")")
  )

  dev.off()
  message("  [OUT] Saved MA plot to ", ma_file)

  message("Finished plotting for ", datatype)
}

datatypes <- c("RNAseq", "RIPseq")
for (dt in datatypes) {
  tryCatch({
    plot_dge_results(dt)
  }, error = function(e) {
    message("Error plotting ", dt, ": ", e$message)
  })
}

message("\n05_dge_plotting.R completed successfully.")
