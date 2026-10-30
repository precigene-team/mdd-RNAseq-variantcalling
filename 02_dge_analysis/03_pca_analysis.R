# [[02_dge_analysis/03_pca_analysis.R]]
# Goal: Perform noise-filtered, sex-aware PCA and generate interactive 3D plot widgets.
#
# Steps:
#   1. Load normalized expression RDS objects (vst/rlog)
#   2. Retrieve raw sample metadata to map sex (Gender) covariates
#   3. Filter the expression matrix to the top 500 most variable genes
#   4. Compute PCA using prcomp() with center = TRUE and scale. = FALSE
#   5. Calculate percentage of variance explained by each PC
#   6. Save PCA coordinates and variance percent tables to results/02_dge_analysis/tables/
#   7. Save PCA prcomp R objects to _data/RData/02_dge_analysis/
#   8. Generate interactive 3D PCA plot widgets colored by sex-specific groups
#   9. Export plots as self-contained HTML files to results/02_dge_analysis/figures/
#
suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(readr)
  library(plotly)
  library(htmlwidgets)
  library(ggplot2)
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
tables_dir <- file.path(project_root, "results/02_dge_analysis/tables")
figures_dir <- file.path(project_root, "results/02_dge_analysis/figures")

dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

# Process PCA for each DataType
# Helper function to remove gender effect (limma removal or manual)
remove_gender_effect <- function(mat, genders) {
  if (requireNamespace("limma", quietly = TRUE)) {
    return(limma::removeBatchEffect(mat, batch = genders))
  } else {
    message("  [INFO] limma package not found. Using manual residualization for gender.")
    mat_corrected <- mat
    overall_mean <- rowMeans(mat)
    for (g in unique(genders)) {
      idx <- which(genders == g)
      if (length(idx) > 0) {
        gender_mean <- rowMeans(mat[, idx, drop = FALSE])
        mat_corrected[, idx] <- mat[, idx] - gender_mean + overall_mean
      }
    }
    return(mat_corrected)
  }
}

# Helper to run PCA, save tables, and save plot
run_and_save_pca <- function(mat, metadata_full, title_suffix, file_suffix, datatype, is_2d = FALSE) {
  # Filter to top 500 variable genes to reduce noise
  rv <- matrixStats::rowVars(mat)
  ntop <- 500
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]

  pca <- prcomp(
    t(mat[select, ]),
    center = TRUE,
    scale. = FALSE
  )
  
  var_exp <- 100 * (pca$sdev^2 / sum(pca$sdev^2))
  
  # PCA DataFrame
  pca_df <- data.frame(
    Sample = rownames(pca$x),
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    PC3 = pca$x[, 3],
    stringsAsFactors = FALSE
  ) %>%
    dplyr::left_join(metadata_full, by = "Sample") %>%
    dplyr::mutate(Group = paste(Gender, SampleType, sep = "_"))
  
  # Save PCA coordinates to results/02_dge_analysis/tables/
  write.csv(
    pca_df,
    file.path(tables_dir, paste0("03_", datatype, "_PCA_coordinates", file_suffix, ".csv")),
    row.names = FALSE,
    quote = FALSE
  )
  
  # Save PCA variance to results/02_dge_analysis/tables/
  write.csv(
    data.frame(
      PC = paste0("PC", seq_along(var_exp)),
      VariancePercent = var_exp
    ),
    file.path(tables_dir, paste0("03_", datatype, "_PCA_variance", file_suffix, ".csv")),
    row.names = FALSE,
    quote = FALSE
  )
  
  group_colors <- c(
    "female_ctrl" = "#fb9a99", # Light Red
    "female_mdd"  = "#b2182b", # Dark Red
    "male_ctrl"   = "#a6cee3", # Light Blue
    "male_mdd"    = "#1f78b4"  # Dark Blue
  )

  if (is_2d) {
    p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group)) +
      geom_point(size = 4, alpha = 0.8) +
      scale_color_manual(values = group_colors) +
      labs(
        title = paste("PCA -", datatype, title_suffix),
        x = sprintf("PC1 (%.1f%%)", var_exp[1]),
        y = sprintf("PC2 (%.1f%%)", var_exp[2])
      ) +
      theme_bw(base_size = 18) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.title = element_blank()
      )
    
    png_path <- file.path(figures_dir, paste0("03_", datatype, "_PCA_2D", file_suffix, ".png"))
    ggsave(
      filename = png_path,
      plot = p,
      width = 8,
      height = 6,
      dpi = 300,
      bg = "white"
    )
  } else {
    p <- plot_ly(
      data = pca_df,
      x = ~PC1,
      y = ~PC2,
      z = ~PC3,
      color = ~Group,
      colors = group_colors,
      text = ~Sample,
      type = "scatter3d",
      mode = "markers",
      marker = list(
        size = 6,
        opacity = 0.85
      ),
      hovertemplate = paste0(
        "<b>%{text}</b><br>",
        "Group: %{color}<br>",
        "PC1: %{x:.2f}<br>",
        "PC2: %{y:.2f}<br>",
        "PC3: %{z:.2f}<extra></extra>"
      )
    ) %>%
      layout(
        title = paste("PCA -", datatype, title_suffix),
        scene = list(
          xaxis = list(title = sprintf("PC1 (%.1f%%)", var_exp[1])),
          yaxis = list(title = sprintf("PC2 (%.1f%%)", var_exp[2])),
          zaxis = list(title = sprintf("PC3 (%.1f%%)", var_exp[3]))
        )
      )
    
    # Save Widget HTML
    saveWidget(
      p,
      file.path(figures_dir, paste0("03_", datatype, "_PCA_3D", file_suffix, ".html")),
      selfcontained = TRUE
    )
  }
  
  return(pca)
}

# Process PCA for each DataType
process_pca <- function(datatype) {
  message("\n==============================")
  message("PCA Analysis for: ", datatype)
  message("==============================")
  
  # Search for the normalized RDS file (vst or rlog normalized) prefixed with 02_
  norm_files <- list.files(
    rdata_dir, 
    pattern = paste0("^02_", datatype, "_.*_normalized\\.rds$"), 
    full.names = TRUE
  )
  
  if (length(norm_files) == 0) {
    stop("No normalized RDS file found for ", datatype)
  }
  
  norm_file <- norm_files[1]
  message("Loading normalized data: ", norm_file)
  rld <- readRDS(norm_file)
  norm_mat <- assay(rld)
  
  # Read metadata to retrieve Gender
  metadata_full <- readr::read_csv(
    file.path(project_root, "_data/_metadata/sample_sheet.csv"),
    show_col_types = FALSE
  ) %>%
    dplyr::mutate(Sample = sub("_1$", "", File1)) %>%
    dplyr::select(Sample, Gender, SampleType) %>%
    dplyr::distinct()

  # 1. Run and save original combined PCA
  message("Running original combined PCA...")
  pca_orig <- run_and_save_pca(norm_mat, metadata_full, "", "", datatype)
  saveRDS(pca_orig, file.path(rdata_dir, paste0("03_", datatype, "_pca.rds")))

  # Get genders in correct column order
  samples_ordered <- data.frame(Sample = colnames(norm_mat), stringsAsFactors = FALSE)
  metadata_aligned <- samples_ordered %>%
    dplyr::left_join(metadata_full, by = "Sample")
  genders <- tolower(metadata_aligned$Gender)
  
  # 2. Run and save sex-corrected combined PCA
  message("Running sex-corrected combined PCA...")
  norm_mat_corrected <- remove_gender_effect(norm_mat, genders)
  pca_corr <- run_and_save_pca(norm_mat_corrected, metadata_full, "(Sex-Corrected)", "_sex_corrected", datatype)
  saveRDS(pca_corr, file.path(rdata_dir, paste0("03_", datatype, "_pca_sex_corrected.rds")))
  
  # 3. Run and save Female-only PCA
  female_samples <- which(genders == "female")
  if (length(female_samples) > 0) {
    message("Running Female-only PCA...")
    norm_mat_female <- norm_mat[, female_samples, drop = FALSE]
    pca_female <- run_and_save_pca(norm_mat_female, metadata_full, "(Female-only)", "_female_only", datatype, is_2d = TRUE)
    saveRDS(pca_female, file.path(rdata_dir, paste0("03_", datatype, "_pca_female_only.rds")))
  } else {
    message("No female samples found for Female-only PCA.")
  }

  # 4. Run and save Male-only PCA
  male_samples <- which(genders == "male")
  if (length(male_samples) > 0) {
    message("Running Male-only PCA...")
    norm_mat_male <- norm_mat[, male_samples, drop = FALSE]
    pca_male <- run_and_save_pca(norm_mat_male, metadata_full, "(Male-only)", "_male_only", datatype, is_2d = TRUE)
    saveRDS(pca_male, file.path(rdata_dir, paste0("03_", datatype, "_pca_male_only.rds")))
  } else {
    message("No male samples found for Male-only PCA.")
  }
  
  message("Finished PCA variants for ", datatype)
}

datatypes <- c("RNAseq")
for (dt in datatypes) {
  tryCatch({
    process_pca(dt)
  }, error = function(e) {
    message("Error during PCA for ", dt, ": ", e$message)
  })
}

message("\n03_pca_analysis.R completed successfully.")
