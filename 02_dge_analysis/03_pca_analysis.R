# [[02_dge_analysis/03_pca_analysis.R]]
# Goal: Perform PCA and generate interactive 3D PCA plot widgets.
#
# Pipeline steps:
#   1. Load normalized expression RDS objects (vst/rlog) and metadata
#   2. Standardize sample name matching and align metadata
#   3. Transpose matrix and compute PCA using prcomp()
#   4. Calculate percentage of variance explained by each PC
#   5. Save PCA coordinates and variance percent tables to results/tables/
#   6. Save PCA prcomp R objects to _data/RData/02_dge_analysis/
#   7. Generate interactive 3D PCA plot widgets using plotly
#   8. Export plots as self-contained HTML files to results/figures/
#

suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(readr)
  library(plotly)
  library(htmlwidgets)
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
tables_dir <- file.path(project_root, "results/tables")
figures_dir <- file.path(project_root, "results/figures")

dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

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
  
  coldata <- colData(rld)
  
  # Run PCA
  pca <- prcomp(
    t(norm_mat),
    center = TRUE,
    scale. = TRUE
  )
  
  # Save PCA prcomp RDS object with prefix 03_
  saveRDS(pca, file.path(rdata_dir, paste0("03_", datatype, "_pca.rds")))
  
  # Compute variance explained
  var_exp <- 100 * (pca$sdev^2 / sum(pca$sdev^2))
  
  # PCA DataFrame
  pca_df <- data.frame(
    Sample = rownames(pca$x),
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    PC3 = pca$x[, 3],
    SampleType = coldata[rownames(pca$x), "SampleType"],
    stringsAsFactors = FALSE
  )
  
  # Save PCA coordinates to results/tables/ with prefix 03_
  write.csv(
    pca_df,
    file.path(tables_dir, paste0("03_", datatype, "_PCA_coordinates.csv")),
    row.names = FALSE,
    quote = FALSE
  )
  
  # Save PCA variance to results/tables/ with prefix 03_
  write.csv(
    data.frame(
      PC = paste0("PC", seq_along(var_exp)),
      VariancePercent = var_exp
    ),
    file.path(tables_dir, paste0("03_", datatype, "_PCA_variance.csv")),
    row.names = FALSE,
    quote = FALSE
  )
  
  # 3D PCA plot
  p <- plot_ly(
    data = pca_df,
    x = ~PC1,
    y = ~PC2,
    z = ~PC3,
    color = ~SampleType,
    colors = c("blue", "red"),
    text = ~Sample,
    type = "scatter3d",
    mode = "markers",
    marker = list(
      size = 6,
      opacity = 0.85
    ),
    hovertemplate = paste0(
      "<b>%{text}</b><br>",
      "SampleType: %{color}<br>",
      "PC1: %{x:.2f}<br>",
      "PC2: %{y:.2f}<br>",
      "PC3: %{z:.2f}<extra></extra>"
    )
  ) %>%
    layout(
      title = paste("PCA -", datatype),
      scene = list(
        xaxis = list(title = sprintf("PC1 (%.1f%%)", var_exp[1])),
        yaxis = list(title = sprintf("PC2 (%.1f%%)", var_exp[2])),
        zaxis = list(title = sprintf("PC3 (%.1f%%)", var_exp[3]))
      )
    )
  
  # Save Widget HTML to results/figures/ with prefix 03_
  saveWidget(
    p,
    file.path(figures_dir, paste0("03_", datatype, "_PCA_3D.html")),
    selfcontained = TRUE
  )
  
  message("Finished PCA for ", datatype)
}

datatypes <- c("RNAseq", "RIPseq")
for (dt in datatypes) {
  tryCatch({
    process_pca(dt)
  }, error = function(e) {
    message("Error during PCA for ", dt, ": ", e$message)
  })
}

message("\n03_pca_analysis.R completed successfully.")
