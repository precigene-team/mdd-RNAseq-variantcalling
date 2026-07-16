# [[02_dge_analysis/04_deg_analysis.R]]
# Goal: Run DESeq2 differential expression analysis (DGE) with stratified sex-specific contrasts
#       and relative methylation interaction modeling, saving results to tables/ and _data/.
#
# Steps:
#   1. Load dds RDS objects for RNA-seq (Input) and RIP-seq (IP)
#   2. Load metadata from sample_sheet.csv and annotate dds objects (Gender, Subject, Fraction)
#   3. Apply a relaxed, peak-aware QC library size gate (2M floor) to protect sample pairings
#   4. Run Baseline DGE (Input/RNA-seq only) for:
#      - Male DGE (stratified)
#      - Female DGE (stratified)
#      - Sex Interaction Model (sex-differential response)
#   5. Run Baseline DGE (IP/RIP-seq only) for:
#      - Male IP DGE (stratified)
#      - Female IP DGE (stratified)
#   6. Run Relative Methylation Interaction DGE on intersecting peaks for:
#      - Male Relative Methylation (paired interaction with nested subjects)
#      - Female Relative Methylation (paired interaction with nested subjects)
#      - Three-Way Interaction (sex-specific methylation change)
#   7. Shrink log fold changes globally using ashr to ensure statistical consistency across all phases
#   8. Save results tables (full and significant) to results/02_dge_analysis/tables/ prefixed with 04_
#   9. Save serialized dds and results objects to _data/RData/02_dge_analysis/ prefixed with 04_

suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(readr)
  library(ashr)
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
meta_file <- file.path(project_root, "_data/_metadata/sample_sheet.csv")

dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(meta_file)) {
  stop("Metadata file not found: ", meta_file)
}

# 1. Load and clean metadata sheet
message("[INFO] Loading metadata sheet...")
metadata <- read_csv(meta_file, show_col_types = FALSE) %>%
  mutate(
    Sample = sub("_1$", "", File1),
    Subject = sub("_(input|m6A)$", "", Sample),
    Fraction = ifelse(DataType == "RNA-seq", "Input", "IP"),
    Gender = factor(Gender),
    SampleType = factor(SampleType)
  ) %>%
  select(Sample, Subject, Gender, SampleType, Fraction) %>%
  distinct()

# Helper function to add full metadata annotations to a dds object
annotate_dds <- function(dds) {
  cd <- as.data.frame(colData(dds))
  cd$Sample <- rownames(cd)
  
  cd_annot <- cd %>%
    select(SampleType) %>%
    mutate(Sample = rownames(cd)) %>%
    left_join(metadata, by = "Sample")
  
  colData(dds)$Gender <- factor(cd_annot$Gender)
  colData(dds)$Subject <- factor(cd_annot$Subject)
  colData(dds)$Fraction <- factor(cd_annot$Fraction)
  
  # Set reference levels safely
  s_type <- factor(colData(dds)$SampleType)
  if ("ctrl" %in% levels(s_type)) {
    colData(dds)$SampleType <- relevel(s_type, ref = "ctrl")
  }
  
  gender_f <- factor(colData(dds)$Gender)
  if ("female" %in% levels(gender_f)) {
    colData(dds)$Gender <- relevel(gender_f, ref = "female")
  }
  
  fraction_f <- factor(colData(dds)$Fraction)
  if ("Input" %in% levels(fraction_f)) {
    colData(dds)$Fraction <- relevel(fraction_f, ref = "Input")
  }
  
  return(dds)
}

# Load dds RDS files
dds_rnaseq_file <- file.path(rdata_dir, "02_RNAseq_dds.rds")
dds_ripseq_file <- file.path(rdata_dir, "02_RIPseq_dds.rds")

if (!file.exists(dds_rnaseq_file) || !file.exists(dds_ripseq_file)) {
  stop("DDS RDS files not found in ", rdata_dir)
}

message("[INFO] Loading and annotating input DDS objects...")
dds_rnaseq <- annotate_dds(readRDS(dds_rnaseq_file))
dds_ripseq <- annotate_dds(readRDS(dds_ripseq_file))

# --- Dynamic QC Gating (Exclude extreme outlier subjects) ---
MIN_LIBRARY_SIZE <- 2e6 

# Calculate total peak-mapped counts for each sample
rna_lib_sizes <- colSums(counts(dds_rnaseq))
rip_lib_sizes <- colSums(counts(dds_ripseq))

# Identify subjects where either Input or IP library falls below the 2M floor
low_qual_rna_subs <- colData(dds_rnaseq)$Subject[rna_lib_sizes < MIN_LIBRARY_SIZE]
low_qual_rip_subs <- colData(dds_ripseq)$Subject[rip_lib_sizes < MIN_LIBRARY_SIZE]
excluded_subjects <- unique(c(as.character(low_qual_rna_subs), as.character(low_qual_rip_subs)))

if (length(excluded_subjects) > 0) {
  message("[QC GATE] Excluding subjects with at least one library < ", 
          format(MIN_LIBRARY_SIZE, scientific = FALSE), " counts: ", 
          paste(excluded_subjects, collapse = ", "))
  
  # Exclude from both datasets to maintain proper subject pairing for interaction models
  dds_rnaseq <- dds_rnaseq[, !colData(dds_rnaseq)$Subject %in% excluded_subjects]
  dds_ripseq <- dds_ripseq[, !colData(dds_ripseq)$Subject %in% excluded_subjects]
  
  # Drop unused subject levels to prevent design matrix rank deficiency
  colData(dds_rnaseq)$Subject <- droplevels(colData(dds_rnaseq)$Subject)
  colData(dds_ripseq)$Subject <- droplevels(colData(dds_ripseq)$Subject)
} else {
  message("[QC GATE] All subjects successfully passed the library size gate (", 
          format(MIN_LIBRARY_SIZE, scientific = FALSE), ").")
}

# Helper function to run DESeq2, shrink, and write tables
run_deseq_and_save <- function(dds_sub, design_formula, coef_regex, contrast_spec, output_name, use_model_matrix = FALSE) {
  if (use_model_matrix) {
    message("  [INFO] Building custom model matrix for formula: ", deparse(design_formula))
    mm <- model.matrix(design_formula, colData(dds_sub))
    # Drop columns that are all 0s (safeguard for nested designs with unequal groups)
    mm_clean <- mm[, colSums(abs(mm)) > 0, drop = FALSE]
    
    # Run DESeq2 with custom model matrix
    dds_sub <- DESeq(dds_sub, full = mm_clean)
  } else {
    design(dds_sub) <- design_formula
    dds_sub <- DESeq(dds_sub)
  }
  
  # Find matching coefficient using regex
  coef_names <- resultsNames(dds_sub)
  
  actual_coef <- NA
  if (!is.null(coef_regex) && !is.na(coef_regex)) {
    actual_coef <- grep(coef_regex, coef_names, value = TRUE)[1]
  }
  
  # Dynamically resolve contrast_spec regex patterns if it is a list
  resolved_contrast <- contrast_spec
  if (is.list(contrast_spec) && is.na(actual_coef)) {
    resolved_contrast <- lapply(contrast_spec, function(x) {
      if (length(x) > 1) {
        return(sapply(x, function(sub_regex) grep(sub_regex, coef_names, value = TRUE)[1]))
      } else {
        return(grep(x, coef_names, value = TRUE)[1])
      }
    })
  }
  
  # Get raw results (using name or contrast)
  if (!is.na(actual_coef)) {
    message("  [INFO] Extracting results using coefficient: ", actual_coef)
    res_raw <- results(dds_sub, name = actual_coef, independentFiltering = FALSE)
  } else {
    message("  [INFO] Extracting results using contrast specification")
    res_raw <- results(dds_sub, contrast = resolved_contrast, independentFiltering = FALSE)
  }
  
  # Shrink log fold changes globally using ashr for perfect statistical consistency
  res <- tryCatch({
    message("  [INFO] Shrinking log fold changes with ashr...")
    if (use_model_matrix || !is.na(actual_coef)) {
      lfcShrink(dds_sub, coef = actual_coef, type = "ashr")
    } else {
      lfcShrink(dds_sub, contrast = resolved_contrast, type = "ashr")
    }
  }, error = function(e) {
    message("  [WARNING] Shrinkage failed. Falling back to default results: ", e$message)
    res_raw
  })
  
  res_df <- as.data.frame(res)
  res_df$stat <- res_raw$stat
  res_df$pvalue <- res_raw$pvalue
  res_df$padj <- res_raw$padj
  res_df$gene_name <- rownames(res_df)
  
  # Reorder columns to make gene_name first
  res_df <- res_df %>% select(gene_name, everything())
  
  # Save full tables with prefix 04_
  write.csv(
    res_df,
    file.path(tables_dir, paste0("04_", output_name, ".csv")),
    row.names = FALSE,
    quote = FALSE
  )
  
  # Primary Significance Filter: pvalue < 0.05 and |log2FC| >= 0.2
  sig_df <- res_df %>%
    filter(
      !is.na(pvalue),
      pvalue < 0.05,
      abs(log2FoldChange) >= 0.2
    )
  
  # Save sig table with prefix 04_
  write.csv(
    sig_df,
    file.path(tables_dir, paste0("04_", output_name, "_sig.csv")),
    row.names = FALSE,
    quote = FALSE
  )
  
  message("  [SUCCESS] Written 04_", output_name, " tables. Significant genes: ", nrow(sig_df))
  
  # Save serialized objects
  saveRDS(res, file.path(rdata_dir, paste0("04_", output_name, "_res_shrunk.rds")))
  saveRDS(res_raw, file.path(rdata_dir, paste0("04_", output_name, "_res_raw.rds")))
  saveRDS(dds_sub, file.path(rdata_dir, paste0("04_", output_name, "_dds.rds")))
}

#--------------------------------------------------
# Phase A: Baseline DGE (Input/RNA-seq Only)
#--------------------------------------------------
message("\n--- Running Phase A: Baseline DGE (RNA-seq / Input) ---")

# A1. Male DGE
message("A1. Running Male Input DGE...")
dds_male_input <- dds_rnaseq[, colData(dds_rnaseq)$Gender == "male"]
colData(dds_male_input)$Subject <- droplevels(colData(dds_male_input)$Subject)
run_deseq_and_save(
  dds_sub = dds_male_input,
  design_formula = ~ SampleType,
  coef_regex = "SampleType_mdd_vs_ctrl",
  contrast_spec = c("SampleType", "mdd", "ctrl"),
  output_name = "deg_male_input"
)

# A2. Female DGE
message("A2. Running Female Input DGE...")
dds_female_input <- dds_rnaseq[, colData(dds_rnaseq)$Gender == "female"]
colData(dds_female_input)$Subject <- droplevels(colData(dds_female_input)$Subject)
run_deseq_and_save(
  dds_sub = dds_female_input,
  design_formula = ~ SampleType,
  coef_regex = "SampleType_mdd_vs_ctrl",
  contrast_spec = c("SampleType", "mdd", "ctrl"),
  output_name = "deg_female_input"
)

# A3. Sex Interaction Model (All Input samples)
message("A3. Running Sex-Differential Input Response (Interaction)...")
run_deseq_and_save(
  dds_sub = dds_rnaseq,
  design_formula = ~ Gender + SampleType + Gender:SampleType,
  coef_regex = "Gender.*male.*SampleType.*mdd",
  contrast_spec = NULL,
  output_name = "deg_sex_interaction_input"
)

#--------------------------------------------------
# Phase B: Baseline DGE (IP/RIP-seq Only)
#--------------------------------------------------
message("\n--- Running Phase B: Baseline DGE (RIP-seq / IP) ---")

# B1. Male IP DGE
message("B1. Running Male IP DGE...")
dds_male_ip <- dds_ripseq[, colData(dds_ripseq)$Gender == "male"]
colData(dds_male_ip)$Subject <- droplevels(colData(dds_male_ip)$Subject)
run_deseq_and_save(
  dds_sub = dds_male_ip,
  design_formula = ~ SampleType,
  coef_regex = "SampleType_mdd_vs_ctrl",
  contrast_spec = c("SampleType", "mdd", "ctrl"),
  output_name = "deg_male_ip"
)

# B2. Female IP DGE
message("B2. Running Female IP DGE...")
dds_female_ip <- dds_ripseq[, colData(dds_ripseq)$Gender == "female"]
colData(dds_female_ip)$Subject <- droplevels(colData(dds_female_ip)$Subject)
run_deseq_and_save(
  dds_sub = dds_female_ip,
  design_formula = ~ SampleType,
  coef_regex = "SampleType_mdd_vs_ctrl",
  contrast_spec = c("SampleType", "mdd", "ctrl"),
  output_name = "deg_female_ip"
)

#--------------------------------------------------
# Phase C: Relative Methylation (MeRIP-seq Interaction Model)
#--------------------------------------------------
message("\n--- Running Phase C: Relative Methylation Interactions ---")

# Combine counts on intersecting genes
common_genes <- intersect(rownames(dds_rnaseq), rownames(dds_ripseq))
combined_counts <- cbind(counts(dds_rnaseq)[common_genes, ], counts(dds_ripseq)[common_genes, ])
colData_combined <- rbind(as.data.frame(colData(dds_rnaseq)), as.data.frame(colData(dds_ripseq)))

dds_combined <- DESeqDataSetFromMatrix(
  countData = combined_counts,
  colData = colData_combined,
  design = ~ 1
)

# C1. Male Relative Methylation
message("C1. Running Male Relative Methylation (Interaction)...")
dds_male_combined <- dds_combined[, colData(dds_combined)$Gender == "male"]
colData(dds_male_combined)$Subject <- droplevels(colData(dds_male_combined)$Subject)

# Reconstruct nested subject factor for male subset
cd_male <- as.data.frame(colData(dds_male_combined))
cd_male <- cd_male %>%
  group_by(SampleType) %>%
  mutate(Subject.nested = factor(as.integer(factor(Subject)))) %>%
  ungroup()
colData(dds_male_combined)$Subject.nested <- factor(cd_male$Subject.nested)

run_deseq_and_save(
  dds_sub = dds_male_combined,
  design_formula = ~ SampleType + SampleType:Subject.nested + Fraction + Fraction:SampleType,
  coef_regex = "(Fraction.*IP.*SampleType.*mdd|SampleType.*mdd.*Fraction.*IP)",
  contrast_spec = NULL,
  output_name = "meth_male_interaction",
  use_model_matrix = TRUE
)

# C2. Female Relative Methylation
message("C2. Running Female Relative Methylation (Interaction)...")
dds_female_combined <- dds_combined[, colData(dds_combined)$Gender == "female"]
colData(dds_female_combined)$Subject <- droplevels(colData(dds_female_combined)$Subject)

# Reconstruct nested subject factor for female subset (unbalanced: 9 ctrl vs 10 mdd)
cd_female <- as.data.frame(colData(dds_female_combined))
cd_female <- cd_female %>%
  group_by(SampleType) %>%
  mutate(Subject.nested = factor(as.integer(factor(Subject)))) %>%
  ungroup()
colData(dds_female_combined)$Subject.nested <- factor(cd_female$Subject.nested)

run_deseq_and_save(
  dds_sub = dds_female_combined,
  design_formula = ~ SampleType + SampleType:Subject.nested + Fraction + Fraction:SampleType,
  coef_regex = "(Fraction.*IP.*SampleType.*mdd|SampleType.*mdd.*Fraction.*IP)",
  contrast_spec = NULL,
  output_name = "meth_female_interaction",
  use_model_matrix = TRUE
)

# C3. Three-Way Interaction (Sex-Specific Methylation Change)
message("C3. Running Three-Way Sex-Specific Methylation Change (Interaction)...")
colData(dds_combined)$Subject <- droplevels(colData(dds_combined)$Subject)

# Reconstruct nested subject factor for all samples
cd_all <- as.data.frame(colData(dds_combined))
cd_all <- cd_all %>%
  group_by(Gender, SampleType) %>%
  mutate(Subject.nested = factor(as.integer(factor(Subject)))) %>%
  ungroup()
colData(dds_combined)$Subject.nested <- factor(cd_all$Subject.nested)

run_deseq_and_save(
  dds_sub = dds_combined,
  design_formula = ~ Gender*SampleType + Gender:SampleType:Subject.nested + Fraction + Fraction:Gender + Fraction:SampleType + Fraction:Gender:SampleType,
  coef_regex = "(Fraction.*IP.*Gender.*male.*SampleType.*mdd|Gender.*male.*SampleType.*mdd.*Fraction.*IP)",
  contrast_spec = NULL,
  output_name = "meth_sex_interaction",
  use_model_matrix = TRUE
)

message("\n04_deg_analysis.R completed successfully.")
