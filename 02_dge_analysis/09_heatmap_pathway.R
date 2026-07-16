# [[02_dge_analysis/09_heatmap_pathway.R]]
# Goal: Replicate pathway leading-edge ComplexHeatmap styling from the neighboring pipeline.

# Steps:
#   1. Load raw DESeq2 dds objects and GSEA Hallmark outputs.
#   2. Select the top activated and suppressed Hallmark pathways.
#   3. Extract leading edge genes for these pathways.
#   4. Compute VST expression, scale rows (z-score), and draw heatmaps using ComplexHeatmap.
#   5. Save pathway heatmaps to results/02_dge_analysis/figures/ prefixed with 09_.

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
  library(tidyverse)
  library(org.Hs.eg.db)
})

# --- Dynamic Pathing Setup ---
project_root <- Sys.getenv("PROJECT_ROOT")
if (project_root == "") {
  curr <- getwd()
  while (curr != "/" && !file.exists(file.path(curr, "run"))) { curr <- dirname(curr) }
  project_root = curr
}
data_dir     <- file.path(project_root, "_data")
results_root <- file.path(project_root, "results", "02_dge_analysis")
rdata_dir    <- file.path(data_dir, "RData", "02_dge_analysis")
tables_dir   <- file.path(results_root, "tables")
figures_dir  <- file.path(results_root, "figures")

dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

# Heatmap parameters matching neighbor pipeline
SINGLE_WIDTH <- 1200
col_fun <- colorRamp2(c(-1.5, 0, 1.5), c("#3B4CC0", "white", "#B40426"))

# Load gene mapper
mapper_file <- file.path(tables_dir, "01_gene_mapper.csv")
if (!file.exists(mapper_file)) {
  stop("Gene mapper file not found: ", mapper_file)
}
gene_mapper <- read_csv(mapper_file, show_col_types = FALSE)

# 1. Dynamically resolve the Ensembl ID column (ignoring peak_id)
valid_id_cols <- colnames(gene_mapper)[!grepl("peak", colnames(gene_mapper), ignore.case = TRUE)]

ensembl_col <- grep("^(gene_id|ensembl_gene_id|ensembl|id)$", valid_id_cols, ignore.case = TRUE, value = TRUE)[1]
if (is.na(ensembl_col)) {
  ensembl_col <- grep("id", valid_id_cols, ignore.case = TRUE, value = TRUE)[1]
}

# 2. Dynamically resolve the Gene Symbol column
symbol_col <- grep("^(gene_name|symbol|external_gene_name)$", colnames(gene_mapper), ignore.case = TRUE, value = TRUE)[1]
if (is.na(symbol_col)) {
  symbol_col <- grep("name|symbol", colnames(gene_mapper), ignore.case = TRUE, value = TRUE)[1]
}

# 3. Fallbacks if regex matching fails entirely
if (is.na(ensembl_col)) ensembl_col <- valid_id_cols[1]
if (is.na(symbol_col)) symbol_col <- colnames(gene_mapper)[min(2, ncol(gene_mapper))]

message("  [INFO] Mapping genes dynamically using columns: ", ensembl_col, " -> ", symbol_col)
gene_map_dict <- setNames(gene_mapper[[symbol_col]], gene_mapper[[ensembl_col]])

clean_pathway_name <- function(pathway_name) {
  x <- sub("^HALLMARK_", "", pathway_name)
  x <- gsub("_", " ", x)
  stringr::str_to_title(x)
}

# Find all GSEA Hallmark CSV tables
gsea_files <- list.files(tables_dir, pattern = "^07_.*_gsea_hallmark\\.csv$", full.names = TRUE)

if (length(gsea_files) == 0) {
  stop("No GSEA Hallmark tables (07_*_gsea_hallmark.csv) found in ", tables_dir)
}

message("\n--- Generating Pathway Leading-Edge Heatmaps (High-Fidelity) ---")

for (gsea_path in gsea_files) {
  # Extract contrast name from filename
  contrast <- sub("^07_", "", sub("_gsea_hallmark\\.csv$", "", basename(gsea_path)))
  
  # Load the corresponding dds object
  dds_file <- file.path(rdata_dir, paste0("04_", contrast, "_dds.rds"))
  if (!file.exists(dds_file)) {
    message("  [INFO] Skipping ", contrast, " (dds file not found).")
    next
  }
  
  df_gsea <- read_csv(gsea_path, show_col_types = FALSE)
  
  # Brain studies have subtle changes, so filter at nominal pvalue < 0.05
  df_sig <- df_gsea %>% filter(pvalue < 0.05)
  if (nrow(df_sig) == 0) {
    message("  [INFO] Skipping ", contrast, " (no significant Hallmark pathways at pvalue < 0.05).")
    next
  }
  
  # Select Top Pathways: up to 3 activated (NES > 0) and up to 3 suppressed (NES < 0)
  top_act <- df_sig %>% filter(NES > 0) %>% arrange(desc(NES)) %>% head(3)
  top_sup <- df_sig %>% filter(NES < 0) %>% arrange(NES) %>% head(3)
  sel_paths <- rbind(top_act, top_sup)
  
  dds <- readRDS(dds_file)
  vsd <- vst(dds, blind = FALSE)
  colData(vsd) <- droplevels(colData(vsd))
  
  gene_rows <- list()
  for (i in seq_len(nrow(sel_paths))) {
    path_name <- sel_paths$Description[i]
    path_clean <- clean_pathway_name(path_name)
    nes_val <- sel_paths$NES[i]
    
    # Extract leading edge genes (Ensembl IDs) and map to symbols
    genes_ens <- unique(unlist(strsplit(sel_paths$core_enrichment[i], "/")))
    genes_sym <- unname(gene_map_dict[genes_ens])
    genes_sym <- genes_sym[!is.na(genes_sym)]
    
    # Intersect with vsd row names (which are Gene Symbols)
    genes <- intersect(genes_sym, rownames(vsd))
    # Select top 10 genes per pathway for visualization clarity
    genes <- head(genes, 10)
    
    if (length(genes) > 0) {
      gene_rows[[length(gene_rows) + 1]] <- data.frame(
        Geneid = genes,
        gene_label = genes,
        pathway_display = stringr::str_wrap(paste0(path_clean, " (NES=", round(nes_val, 2), ")"), width = 25),
        stringsAsFactors = FALSE
      )
    }
  }
  
  if (length(gene_rows) == 0) next
  heatmap_df <- do.call(rbind, gene_rows)
  heatmap_df <- heatmap_df[!duplicated(heatmap_df$Geneid), ]
  
  # Extract expression matrix
  mat <- assay(vsd)[heatmap_df$Geneid, , drop=FALSE]
  rownames(mat) <- heatmap_df$gene_label
  
  # Clean up row names to prevent NA mapping errors
  valid_rows <- !is.na(rownames(mat))
  mat <- mat[valid_rows, , drop=FALSE]
  heatmap_df <- heatmap_df[valid_rows, ]

  # Format sample column names into shorthand immediately so pairing patterns work
  colnames(mat) <- colnames(mat) %>%
    str_replace_all("(?i)female", "F") %>%
    str_replace_all("(?i)male", "M") %>%
    str_replace_all("(?i)ctrl", "Ctrl") %>%
    str_replace_all("(?i)mdd", "MDD") %>%
    str_replace_all("(?i)input", "In") %>%
    str_replace_all("(?i)m6A_1|m6a_1|m6A|m6a", "IP") %>%
    str_replace_all("_+", "_") %>%
    str_remove("^_|_$")
  
  # Update the underlying vsd object names to match the new format
  colnames(vsd) <- colnames(mat)

  # Track if this contrast started with fraction design matrices
  had_fraction_originally <- "Fraction" %in% colnames(colData(vsd)) && length(unique(colData(vsd)$Fraction)) > 1
  is_paired_delta_mode <- FALSE

  # Scale fractions independently by checking the cleaned matrix column names directly
  input_cols <- which(grepl("_In$", colnames(mat)))
  ip_cols    <- which(grepl("_IP$", colnames(mat)))
  
  if (length(input_cols) > 0 && length(ip_cols) > 0) {
    message("  [INFO] Scaling Input and IP fractions independently...")
    mat_input <- mat[, input_cols, drop=FALSE]
    mat_ip    <- mat[, ip_cols, drop=FALSE]
    
    # Row-scale each fraction separately
    mat_input_scaled <- t(scale(t(mat_input)))
    mat_ip_scaled    <- t(scale(t(mat_ip)))
    
    mat_input_scaled[is.nan(mat_input_scaled)] <- 0
    mat_ip_scaled[is.nan(mat_ip_scaled)]       <- 0
    
    # Reassemble to tracking matrix
    mat_scaled <- mat
    mat_scaled[, input_cols] <- mat_input_scaled
    mat_scaled[, ip_cols]    <- mat_ip_scaled
    
    # --- Compute Delta Matrix (IP - Input) per paired Subject ---
    sample_names <- colnames(mat_scaled)
    subject_ids  <- unique(sub("_(In|IP)$", "", sample_names))
    
    delta_mat <- matrix(NA, nrow = nrow(mat_scaled), ncol = length(subject_ids),
                        dimnames = list(rownames(mat_scaled), subject_ids))
    
    for (sub_id in subject_ids) {
      in_idx  <- which(sample_names == paste0(sub_id, "_In"))
      ip_idx  <- which(sample_names == paste0(sub_id, "_IP"))
      
      if (length(in_idx) == 1 && length(ip_idx) == 1) {
        delta_mat[, sub_id] <- mat_scaled[, ip_idx] - mat_scaled[, in_idx]
      }
    }
    
    # Drop unmatched samples
    mat <- delta_mat[, colSums(is.na(delta_mat)) == 0, drop = FALSE]
    message("  [INFO] Delta matrix generated successfully for ", ncol(mat), " paired subjects.")
    
    # Update downstream annotation metadata to match the new delta columns
    vsd <- vsd[, which(colnames(vsd) %in% paste0(colnames(mat), "_In"))]
    colnames(vsd) <- sub("_In$", "", colnames(vsd))
    is_paired_delta_mode <- TRUE
    
  } else {
    # Standard scaling fallback for non-fraction baseline DGE designs
    mat <- t(scale(t(mat)))
  }
  
  # Cap values at +/- 1.5 for optimal visual dynamic range
  mat <- pmin(pmax(mat, -1.5), 1.5)
  
  # Define annotations dynamically depending on whether Fraction and Sex Interaction are active
  is_sex_interaction <- grepl("sex_interaction", contrast)
  
  if (is_sex_interaction) {
    if (is_paired_delta_mode) {
      # Columns are now collapsed to subjects, no distinct fraction columns remain in mat
      col_meta <- as.data.frame(colData(vsd)[, c("SampleType", "Gender"), drop=FALSE])
      colnames(col_meta) <- c("Condition", "Sex")
      col_meta$Condition <- tolower(as.character(col_meta$Condition))
      col_meta$Sex <- tolower(as.character(col_meta$Sex))
      
      col_order <- order(factor(col_meta$Sex, levels = c("male", "female")), col_meta$Condition)
      mat <- mat[, col_order, drop=FALSE]
      col_meta_ordered <- col_meta[col_order, , drop=FALSE]
      
      cond_colors <- c("ctrl" = "grey85", "mdd" = "#E41A1C")
      sex_colors <- c("male" = "#377EB8", "female" = "#F781BF")
      ha_list <- list(Condition = cond_colors, Sex = sex_colors)
      
      top_ha <- HeatmapAnnotation(
        df = col_meta_ordered,
        col = ha_list,
        show_legend = FALSE,
        annotation_name_gp = gpar(fontsize = 13, fontface = "bold")
      )
      
      column_split_factor <- factor(col_meta_ordered$Sex, levels = c("male", "female"))
    } else if (had_fraction_originally) {
      col_meta <- as.data.frame(colData(vsd)[, c("SampleType", "Fraction", "Gender"), drop=FALSE])
      colnames(col_meta) <- c("Condition", "Fraction", "Sex")
      col_meta$Condition <- tolower(as.character(col_meta$Condition))
      col_meta$Sex <- tolower(as.character(col_meta$Sex))
      col_meta$Fraction <- ifelse(tolower(as.character(col_meta$Fraction)) == "ip", "IP", "input")
      
      col_order <- order(
        factor(col_meta$Sex, levels = c("male", "female")), 
        factor(col_meta$Fraction, levels = c("input", "IP")), 
        col_meta$Condition
      )
      mat <- mat[, col_order, drop=FALSE]
      col_meta_ordered <- col_meta[col_order, , drop=FALSE]
      
      cond_colors <- c("ctrl" = "grey85", "mdd" = "#E41A1C")
      frac_colors <- c("input" = "#4DAF4A", "IP" = "#984EA3")
      sex_colors <- c("male" = "#377EB8", "female" = "#F781BF")
      ha_list <- list(Condition = cond_colors, Fraction = frac_colors, Sex = sex_colors)
      
      top_ha <- HeatmapAnnotation(
        df = col_meta_ordered,
        col = ha_list,
        show_legend = FALSE,
        annotation_name_gp = gpar(fontsize = 13, fontface = "bold")
      )
      
      column_split_factor <- list(
        factor(col_meta_ordered$Sex, levels = c("male", "female")),
        factor(col_meta_ordered$Fraction, levels = c("input", "IP"))
      )
    } else {
      col_meta <- as.data.frame(colData(vsd)[, c("SampleType", "Gender"), drop=FALSE])
      colnames(col_meta) <- c("Condition", "Sex")
      col_meta$Condition <- tolower(as.character(col_meta$Condition))
      col_meta$Sex <- tolower(as.character(col_meta$Sex))
      
      col_order <- order(factor(col_meta$Sex, levels = c("male", "female")), col_meta$Condition)
      mat <- mat[, col_order, drop=FALSE]
      col_meta_ordered <- col_meta[col_order, , drop=FALSE]
      
      cond_colors <- c("ctrl" = "grey85", "mdd" = "#E41A1C")
      sex_colors <- c("male" = "#377EB8", "female" = "#F781BF")
      ha_list <- list(Condition = cond_colors, Sex = sex_colors)
      
      top_ha <- HeatmapAnnotation(
        df = col_meta_ordered,
        col = ha_list,
        show_legend = FALSE,
        annotation_name_gp = gpar(fontsize = 13, fontface = "bold")
      )
      
      column_split_factor <- factor(col_meta_ordered$Sex, levels = c("male", "female"))
    }
  } else {
    # Non-sex interaction contrasts
    if (is_paired_delta_mode) {
      col_meta <- as.data.frame(colData(vsd)[, c("SampleType", "Gender"), drop=FALSE])
      colnames(col_meta) <- c("Condition", "Sex")
      col_meta$Condition <- tolower(as.character(col_meta$Condition))
      col_meta$Sex <- tolower(as.character(col_meta$Sex))
      
      col_order <- order(col_meta$Condition)
      mat <- mat[, col_order, drop=FALSE]
      col_meta_ordered <- col_meta[col_order, , drop=FALSE]
      
      cond_colors <- c("ctrl" = "grey85", "mdd" = "#E41A1C")
      sex_colors <- c("male" = "#377EB8", "female" = "#F781BF")
      ha_list <- list(Condition = cond_colors, Sex = sex_colors)
      
      top_ha <- HeatmapAnnotation(
        df = col_meta_ordered,
        col = ha_list,
        show_legend = FALSE,
        annotation_name_gp = gpar(fontsize = 13, fontface = "bold")
      )
      
      column_split_factor <- factor(col_meta_ordered$Condition, levels = c("ctrl", "mdd"))
    } else if (had_fraction_originally) {
      col_meta <- as.data.frame(colData(vsd)[, c("SampleType", "Fraction", "Gender"), drop=FALSE])
      colnames(col_meta) <- c("Condition", "Fraction", "Sex")
      col_meta$Condition <- tolower(as.character(col_meta$Condition))
      col_meta$Sex <- tolower(as.character(col_meta$Sex))
      col_meta$Fraction <- ifelse(tolower(as.character(col_meta$Fraction)) == "ip", "IP", "input")
      
      col_order <- order(factor(col_meta$Fraction, levels = c("input", "IP")), col_meta$Condition)
      mat <- mat[, col_order, drop=FALSE]
      col_meta_ordered <- col_meta[col_order, , drop=FALSE]
      
      cond_colors <- c("ctrl" = "grey85", "mdd" = "#E41A1C")
      frac_colors <- c("input" = "#4DAF4A", "IP" = "#984EA3")
      sex_colors <- c("male" = "#377EB8", "female" = "#F781BF")
      ha_list <- list(Condition = cond_colors, Fraction = frac_colors, Sex = sex_colors)
      
      top_ha <- HeatmapAnnotation(
        df = col_meta_ordered,
        col = ha_list,
        show_legend = FALSE,
        annotation_name_gp = gpar(fontsize = 13, fontface = "bold")
      )
      
      column_split_factor <- factor(col_meta_ordered$Fraction, levels = c("input", "IP"))
    } else {
      col_meta <- as.data.frame(colData(vsd)[, c("SampleType", "Gender"), drop=FALSE])
      colnames(col_meta) <- c("Condition", "Sex")
      col_meta$Condition <- tolower(as.character(col_meta$Condition))
      col_meta$Sex <- tolower(as.character(col_meta$Sex))
      
      col_order <- order(col_meta$Condition)
      mat <- mat[, col_order, drop=FALSE]
      col_meta_ordered <- col_meta[col_order, , drop=FALSE]
      
      cond_colors <- c("ctrl" = "grey85", "mdd" = "#E41A1C")
      sex_colors <- c("male" = "#377EB8", "female" = "#F781BF")
      ha_list <- list(Condition = cond_colors, Sex = sex_colors)
      
      top_ha <- HeatmapAnnotation(
        df = col_meta_ordered,
        col = ha_list,
        show_legend = FALSE,
        annotation_name_gp = gpar(fontsize = 13, fontface = "bold")
      )
      
      column_split_factor <- factor(col_meta_ordered$Condition, levels = c("ctrl", "mdd"))
    }
  }
  
  pathway_factor <- factor(heatmap_df$pathway_display, levels = unique(heatmap_df$pathway_display))
  
  n_pw <- nlevels(pathway_factor)
  pw_palette <- setNames(
    colorRampPalette(RColorBrewer::brewer.pal(max(3, min(n_pw, 8)), "Dark2"))(n_pw),
    levels(pathway_factor)
  )
  
  get_descriptive_title <- function(cntrst) {
    switch(cntrst,
      "deg_male_input" = "Male Baseline Gene Expression (RNA-seq)",
      "deg_female_input" = "Female Baseline Gene Expression (RNA-seq)",
      "deg_sex_interaction_input" = "Sex-Differential Gene Expression Interaction (RNA-seq)",
      "deg_male_ip" = "Male Absolute Methylated RNA Abundance (IP)",
      "deg_female_ip" = "Female Absolute Methylated RNA Abundance (IP)",
      "meth_male_interaction" = "Male Relative Methylation (IP vs Input Delta)",
      "meth_female_interaction" = "Female Relative Methylation (IP vs Input Delta)",
      "meth_sex_interaction" = "Sex-Differential Relative Methylation (3-Way Interaction Delta)",
      stringr::str_to_title(gsub("_", " ", cntrst))
    )
  }
  clean_title <- get_descriptive_title(contrast)
  wrapped_title <- stringr::str_wrap(clean_title, width = 60)
  
  ht <- Heatmap(
    mat,
    name = ifelse(is_paired_delta_mode, "Δ z-score", "z-score"),
    col = col_fun,
    row_split = pathway_factor,
    row_title_side = "right",
    row_title_rot = 0,
    row_title_gp = gpar(fontsize = 12, fontface = "bold"),
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    column_split = column_split_factor,
    column_gap = unit(2, "mm"),
    show_row_names = FALSE,
    right_annotation = rowAnnotation(
      id = anno_text(rownames(mat), gp = gpar(fontsize = 11)),
      Pathway = pathway_factor,
      col = list(Pathway = pw_palette),
      show_legend = FALSE,
      annotation_name_rot = 45,
      annotation_name_side = "bottom",
      annotation_name_gp = gpar(fontsize = 12, fontface = "bold"),
      show_annotation_name = c(id = FALSE, Pathway = TRUE)
    ),
    top_annotation = top_ha,
    column_title = wrapped_title,
    column_title_gp = gpar(fontsize = 16, fontface = "bold"),
    column_names_gp = gpar(fontsize = 8),
    column_names_rot = 45,
    show_heatmap_legend = FALSE
  )
  
  # Custom stacked legend list
  lgd_list <- list()
  lgd_list[[length(lgd_list) + 1]] <- Legend(title = "Condition", 
                                            at = c("ctrl", "mdd"), 
                                            labels = c("Control", "MDD"),
                                            legend_gp = gpar(fill = cond_colors))
  if (had_fraction_originally && !is_paired_delta_mode) {
    lgd_list[[length(lgd_list) + 1]] <- Legend(title = "Fraction", 
                                              at = c("input", "IP"), 
                                              labels = c("Input (RNA)", "IP (m6A)"),
                                              legend_gp = gpar(fill = frac_colors))
  }
  
  # Check if Sex variable exists inside ha_list mapping context
  if ("Sex" %in% names(ha_list)) {
    sex_colors <- ha_list[["Sex"]]
    lgd_list[[length(lgd_list) + 1]] <- Legend(title = "Sex", 
                                              at = c("male", "female"), 
                                              labels = c("Male", "Female"),
                                              legend_gp = gpar(fill = sex_colors))
  }
  
  lgd_list[[length(lgd_list) + 1]] <- Legend(title = ifelse(is_paired_delta_mode, "Δ z-score", "z-score"), 
                                            col_fun = col_fun, 
                                            at = c(-1.5, 0, 1.5))
  
  pd <- packLegend(list = lgd_list, direction = "vertical")
  
  # Save Heatmap plot with spacing parameters matching neighbor pipeline
  png_path <- file.path(figures_dir, paste0("09_heatmap_pathway_", contrast, ".png"))
  png(png_path, width = max(1200, 700 + ncol(mat)*45), height = 250 + nrow(mat)*25, res = 150)
  
  draw(ht, 
       annotation_legend_list = pd,
       heatmap_legend_side = "left", 
       annotation_legend_side = "left",
       padding = unit(c(5, 25, 5, 10), "mm"))
  
  dev.off()
  
  message("  [SUCCESS] Saved heatmap plot -> ", basename(png_path))
}

message("\n[OK] Pathway heatmaps completed successfully.\n")