# [[02_dge_analysis/08_enrichment_plotting.R]]
# Goal: Generate publication-grade GSEA NES Barplots for GO and Hallmark pathways.
#       Follows the aesthetic conventions of the neighboring RNAseq_pipeline.
#
# Steps:
#   1. Find GSEA results files results/02_dge_analysis/tables/07_*_gsea_*.csv
#   2. Load enrichment tables and format pathway names.
#   3. Select the top activated (NES > 0) and suppressed (NES < 0) pathways.
#   4. Render horizontal barplots with red/blue coloring and facet layout.
#   5. Save the plots as results/02_dge_analysis/figures/08_{contrast}_gsea_{type}_nes.png

suppressPackageStartupMessages({
  library(ggplot2)
  library(readr)
  library(dplyr)
  library(stringr)
  library(ggtext)
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
tables_dir <- file.path(project_root, "results/02_dge_analysis/tables")
figures_dir <- file.path(project_root, "results/02_dge_analysis/figures")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

# Helper: clean Hallmark/GO names for professional display
clean_pathway_name <- function(pathway_name) {
  # Clean MSigDB Hallmark prefix
  x <- sub("^HALLMARK_", "", pathway_name)
  # Replace underscores with spaces
  x <- gsub("_", " ", x)
  # Convert to Title Case
  x <- stringr::str_to_title(x)
  return(x)
}

# Find all GSEA files (both Hallmark and Functional GO)
gsea_files <- list.files(tables_dir, pattern = "^07_.*_gsea_.*\\.csv$", full.names = TRUE)

if (length(gsea_files) == 0) {
  stop("No GSEA results tables (07_*_gsea_*.csv) found in ", tables_dir)
}

message("\n--- Generating GSEA NES Barplots ---")

for (f in gsea_files) {
  filename <- basename(f)
  # Parse details from name
  contrast_name <- sub("^07_", "", sub("_gsea_.*\\.csv$", "", filename))
  gsea_type <- ifelse(grepl("hallmark", filename), "hallmark", "go")
  
  df <- read_csv(f, show_col_types = FALSE)
  
  if (nrow(df) == 0) {
    message("  [Skip] Empty file: ", filename)
    next
  }
  
  message("\n>>> Plotting: ", contrast_name, " (Type: ", gsea_type, ")")
  
  # Format pathways and split by direction
  plot_df <- df %>%
    mutate(
      Direction = ifelse(NES > 0, "Activated", "Suppressed"),
      Direction = factor(Direction, levels = c("Activated", "Suppressed")),
      Description = vapply(Description, clean_pathway_name, character(1))
    )
  
  # Pick top 10 activated and top 10 suppressed (ordered by absolute NES)
  top_active <- plot_df %>%
    filter(Direction == "Activated") %>%
    arrange(desc(abs(NES))) %>%
    head(10)
  
  top_suppressed <- plot_df %>%
    filter(Direction == "Suppressed") %>%
    arrange(desc(abs(NES))) %>%
    head(10)
  
  plot_subset <- bind_rows(top_active, top_suppressed)
  
  if (nrow(plot_subset) == 0) {
    message("    No significant terms to plot after filtering.")
    next
  }
  
  # Calculate dynamic plot height based on rows
  n_paths <- nrow(plot_subset)
  plot_height <- max(6, n_paths * 0.65)
  
  # Generate horizontal barplot
  p <- ggplot(plot_subset, aes(x = NES, y = reorder(Description, NES), fill = Direction)) +
    geom_col(width = 0.75, alpha = 0.85) +
    facet_grid(Direction ~ ., scales = "free_y", space = "free_y") +
    scale_fill_manual(
      values = c("Activated" = "#e31a1c", "Suppressed" = "#1f78b4"),
      name = "Enrichment Status"
    ) +
    scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 45)) +
    geom_vline(xintercept = 0, linewidth = 0.5, color = "grey30") +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.major.y = element_blank(),
      axis.text.y = element_text(size = 9, face = "bold", lineheight = 0.85),
      axis.title.y = element_blank(),
      plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5),
      strip.background = element_blank(),
      strip.text = element_blank(),
      legend.position = "bottom"
    ) +
    labs(
      title = paste0("GSEA NES - ", clean_pathway_name(contrast_name)),
      subtitle = paste0("Database: ", ifelse(gsea_type == "hallmark", "MSigDB Hallmark", "GO Biological Process")),
      x = "Normalized Enrichment Score (NES)"
    )
  
  # Save plot
  output_filename <- paste0("08_", contrast_name, "_gsea_", gsea_type, "_nes.png")
  output_path <- file.path(figures_dir, output_filename)
  
  ggsave(
    filename = output_path,
    plot = p,
    width = 8.5,
    height = plot_height,
    dpi = 300,
    bg = "white",
    limitsize = FALSE
  )
  
  message("    [SUCCESS] Saved plot to: ", basename(output_path))
}

message("\n08_enrichment_plotting.R completed successfully.")
