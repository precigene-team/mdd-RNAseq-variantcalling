# [[02_dge_analysis/02_rna_processing.R]]
# Goal: Filter low-expression genes and normalize RNA-seq (Input) counts using DESeq2.
#
# Pipeline steps:
#   1. Load annotated counts and sample metadata sheet
#   2. Clean metadata: map sample names (strip _1) and align DataType to RNAseq
#   3. Intersect expression columns and metadata samples to align common samples
#   4. Split common metadata and counts by DataType
#   5. Filter out low-expressed genes (cumulative counts >= 10)
#   6. Save filtered count matrices to results/02_dge_analysis/tables/02_{DataType}_filtered_counts.csv
#   7. Load DESeqDataSet, estimate size factors, and normalize using vst() or rlog()
#   8. Save normalized matrices to results/02_dge_analysis/tables/02_{DataType}_{norm_method}_normalized_counts.csv
#   9. Save serialized dds and normalized RDS objects to _data/RData/02_dge_analysis/
#

suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(readr)
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

# Input files
expr_file <- file.path(project_root, "results/02_dge_analysis/tables/01_annotated_counts.csv")
meta_file <- file.path(project_root, "_data/_metadata/sample_sheet.csv")

# Output directories
rdata_dir <- file.path(project_root, "_data/RData/02_dge_analysis")
tables_dir <- file.path(project_root, "results/02_dge_analysis/tables")

dir.create(rdata_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

# Read data
expr_df <- read_csv(expr_file, show_col_types = FALSE)
metadata <- read_csv(meta_file, show_col_types = FALSE)

# Ensure first column is named gene_name and collapse any duplicates
colnames(expr_df)[1] <- "gene_name"

expr_df <- expr_df %>%
  mutate(gene_name = trimws(as.character(gene_name))) %>%
  filter(!is.na(gene_name), gene_name != "") %>%
  group_by(gene_name) %>%
  summarise(
    across(
      where(is.numeric),
      ~ sum(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  )

message("Expression matrix contains ", nrow(expr_df), " unique genes.")

#--------------------------------------------------
# Metadata Transformation & Alignment
#--------------------------------------------------
metadata <- metadata %>%
  mutate(
    Sample = sub("_1$", "", File1),
    # Map DataType to match downstream expectations (RNAseq)
    DataType = case_when(
      DataType == "RNA-seq" ~ "RNAseq",
      TRUE ~ DataType
    )
  )

required_cols <- c("Sample", "DataType", "SampleType")
missing_cols <- setdiff(required_cols, colnames(metadata))
if (length(missing_cols) > 0) {
  stop("Missing metadata columns: ", paste(missing_cols, collapse = ", "))
}

# Keep samples present in both tables
sample_cols <- colnames(expr_df)[-1]
common_samples <- intersect(sample_cols, metadata$Sample)
if (length(common_samples) == 0) {
  stop("No overlapping samples between expression matrix and metadata.")
}
message("Found ", length(common_samples), " common samples.")

expr_df <- expr_df %>% select(gene_name, all_of(common_samples))
metadata <- metadata %>% filter(Sample %in% common_samples)

#--------------------------------------------------
# Process each DataType independently
#--------------------------------------------------
sample_groups <- split(metadata, metadata$DataType)

process_datatype <- function(meta_subset, datatype, expr_df) {
  message("\nProcessing: ", datatype)
  
  # Remove samples with missing SampleType
  meta_subset <- meta_subset %>%
    filter(!is.na(SampleType), SampleType != "")
  
  if (nrow(meta_subset) == 0) {
    stop("No samples with valid SampleType for ", datatype)
  }
  
  sample_ids <- meta_subset$Sample
  
  # Build count matrix
  count_df <- expr_df %>%
    select(gene_name, all_of(sample_ids)) %>%
    as.data.frame()
  
  rownames(count_df) <- count_df$gene_name
  count_df$gene_name <- NULL
  
  count_mat <- as.matrix(count_df)
  storage.mode(count_mat) <- "integer"
  
  # Filter genes
  keep <- rowSums(count_mat) >= 10
  count_mat <- count_mat[keep, , drop = FALSE]
  message("Retained ", nrow(count_mat), " genes.")
  
  # Save filtered counts to results/02_dge_analysis/tables/ with prefix 02_
  write.csv(
    data.frame(gene_name = rownames(count_mat), count_mat, check.names = FALSE),
    file.path(tables_dir, paste0("02_", datatype, "_filtered_counts.csv")),
    row.names = FALSE,
    quote = FALSE
  )
  
  # Create colData
  coldata <- meta_subset %>%
    select(Sample, SampleType) %>%
    distinct() %>%
    as.data.frame()
  rownames(coldata) <- coldata$Sample
  coldata$Sample <- NULL
  coldata$SampleType <- factor(coldata$SampleType)
  
  # Reorder colData to match count matrix columns
  coldata <- coldata[colnames(count_mat), , drop = FALSE]
  
  # Design formula
  if (nlevels(coldata$SampleType) > 1) {
    design_formula <- ~ SampleType
  } else {
    warning("Only one SampleType detected for ", datatype, ". Using design = ~1")
    design_formula <- ~ 1
  }
  
  # DESeq2 object
  dds <- DESeqDataSetFromMatrix(
    countData = count_mat,
    colData = coldata,
    design = design_formula
  )
  
  dds <- estimateSizeFactors(dds)
  
  # Normalization
  # Set threshold to 20 samples to avoid extremely slow rlog execution for this cohort (39 samples)
  if (ncol(dds) > 20) {
    message("Using vst() normalization")
    rld <- vst(dds, blind = TRUE)
    norm_method <- "vst"
  } else {
    message("Using rlog() normalization")
    rld <- rlog(dds, blind = TRUE)
    norm_method <- "rlog"
  }
  
  norm_mat <- assay(rld)
  
  # Save normalized matrix to results/02_dge_analysis/tables/ with prefix 02_
  write.csv(
    data.frame(gene_name = rownames(norm_mat), norm_mat, check.names = FALSE),
    file.path(tables_dir, paste0("02_", datatype, "_", norm_method, "_normalized_counts.csv")),
    row.names = FALSE,
    quote = FALSE
  )
  
  # Save R objects with prefix 02_
  saveRDS(dds, file.path(rdata_dir, paste0("02_", datatype, "_dds.rds")))
  saveRDS(rld, file.path(rdata_dir, paste0("02_", datatype, "_", norm_method, "_normalized.rds")))
  
  message("Finished processing ", datatype)
}

# Run analysis
invisible(lapply(names(sample_groups), function(datatype) {
  process_datatype(sample_groups[[datatype]], datatype, expr_df)
}))

message("\n02_rna_processing.R completed successfully.")
