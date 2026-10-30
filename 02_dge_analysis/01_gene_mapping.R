# [[02_dge_analysis/01_gene_mapping.R]]
# Goal: Map Ensembl Gene IDs (ENSG) in the gene-level count matrix to HGNC gene
#       symbols using the Ensembl GTF, restrict to the Input (RNA-seq) samples,
#       filter for valid biotypes, and collapse by gene.
#
# Pipeline steps:
#   1. Load raw gene-level count matrix (from featureCounts / Ensembl GTF, col 1 = ENSG)
#   2. Restrict the count matrix to the Input samples from the sample sheet
#   3. Build the ENSG -> gene symbol / biotype table from the Ensembl GTF
#   4. Join the annotation onto the expression count matrix
#   5. Filter the matrix strictly for valid biotypes (e.g. protein_coding)
#   6. Collapse duplicate HGNC symbols by sum (multiple ENSG per symbol)
#   7. Save annotated counts to results/02_dge_analysis/tables/01_annotated_counts.csv
#   8. Generate and save mapping table to results/02_dge_analysis/tables/01_gene_mapper.csv

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(yaml)
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

# Config
config_path <- file.path(project_root, "02_dge_analysis/config.yaml")
if (!file.exists(config_path)) {
  stop("Config file not found: ", config_path)
}
cfg            <- yaml::read_yaml(config_path)
valid_biotypes <- cfg$ensembl_map$valid_biotypes

# Paths
counts_path    <- file.path(project_root, "_data/featurecounts/raw_counts.csv")
meta_path      <- file.path(project_root, "_data/_metadata/sample_sheet.csv")
gtf_path       <- file.path(project_root, "_data/genome/Homo_sapiens.GRCh38.115.gtf")
out_dir        <- file.path(project_root, "results/02_dge_analysis/tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

rna_counts_path  <- file.path(out_dir, "01_annotated_counts.csv")
gene_mapper_path <- file.path(out_dir, "01_gene_mapper.csv")

# Assertions
if (!file.exists(counts_path)) stop("Raw counts file not found: ", counts_path)
if (!file.exists(meta_path))   stop("Sample sheet not found: ", meta_path)
if (!file.exists(gtf_path))    stop("Ensembl GTF not found: ", gtf_path)

# 1. Load raw gene-level count matrix
message("[LOAD] Reading raw gene counts from ", counts_path)
expr_raw <- read_csv(counts_path, show_col_types = FALSE)
colnames(expr_raw)[1] <- "ensg"
message("  Loaded matrix: ", nrow(expr_raw), " rows x ", ncol(expr_raw), " columns")

if (!all(vapply(colnames(expr_raw)[-1], function(x) is.numeric(expr_raw[[x]]), logical(1)))) {
  stop("Non-numeric count columns in ", counts_path)
}

# 2. Restrict the matrix to the Input (RNA-seq) libraries
metadata <- read_csv(meta_path, show_col_types = FALSE)
input_samples <- metadata %>%
  filter(DataType == "RNA-seq") %>%
  mutate(Sample = sub("_1$", "", File1)) %>%
  pull(Sample) %>%
  unique()

if (length(input_samples) == 0) {
  stop("No RNA-seq samples found in ", meta_path)
}
missing_samples <- setdiff(input_samples, colnames(expr_raw)[-1])
if (length(missing_samples) > 0) {
  stop("Input samples missing from the count matrix: ", paste(missing_samples, collapse = ", "))
}

expr_raw <- expr_raw %>% select(ensg, all_of(input_samples))
message("[INPUT] Keeping ", length(input_samples), " Input samples, dropping the IP libraries (",
        ncol(expr_raw) - 1, " count columns)")

# 3. Build the ENSG -> symbol / biotype table from the Ensembl GTF
# GTF attribute field: gene_id "ENSG..."; gene_name "SOME"; gene_biotype "protein_coding";
message("[ANNOT] Reading Ensembl GTF: ", gtf_path)
gtf_lines <- readLines(gtf_path)
gtf_lines <- gtf_lines[!startsWith(gtf_lines, "#") & trimws(gtf_lines) != ""]
attr <- sub("^[^\t]*\t[^\t]*\t[^\t]*\t[^\t]*\t[^\t]*\t", "", gtf_lines)

keep_lines <- grepl('gene_id "', attr)
gene_annot <- data.frame(
  ensg      = sub('.*gene_id "([^"]+)".*',    "\\1", attr[keep_lines]),
  gene_name = sub('.*gene_name "([^"]+)".*',  "\\1", attr[keep_lines]),
  gene_type = sub('.*gene_biotype "([^"]+)".*', "\\1", attr[keep_lines])
) %>%
  # Ensembl writes gene_id into gene_name when no symbol is assigned
  filter(gene_name != ensg) %>%
  distinct(ensg, .keep_all = TRUE)

message("  Genes with a symbol in the GTF: ", nrow(gene_annot))
stopifnot("No genes parsed from the GTF" = nrow(gene_annot) > 0)

# 4. Join ENSG -> symbol onto the count matrix
expr_ann <- expr_raw %>% left_join(gene_annot, by = "ensg")
n_unmapped <- sum(is.na(expr_ann$gene_name))
message("  ENSG mapped to a symbol: ", nrow(expr_ann) - n_unmapped, " / ", nrow(expr_ann))

# 5. Filter for valid biotypes
expr_ann <- expr_ann %>%
  filter(!is.na(gene_name), gene_name != "", gene_type %in% valid_biotypes)

message("[BIOTYPE] Rows remaining after biotype filter: ", nrow(expr_ann))
stopifnot("No genes match the valid biotypes" = nrow(expr_ann) > 0)

# 6. Collapse duplicate gene symbols by sum (multiple ENSG per symbol)
n_rows         <- nrow(expr_ann)
n_unique_genes <- n_distinct(expr_ann$gene_name)

if (n_rows == n_unique_genes) {
  message("  No duplicate gene symbols found, skipping collapse.")
  expr_collapsed <- expr_ann %>% select(-ensg, -gene_type) %>% relocate(gene_name)
} else {
  message("  Collapsing ", n_rows - n_unique_genes, " duplicate gene entries by sum...")
  expr_collapsed <- expr_ann %>%
    select(-ensg, -gene_type) %>%
    group_by(gene_name) %>%
    summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE)), .groups = "drop")
}

# Assertion: no duplicate gene names in output
stopifnot("Duplicate gene names remain after collapse" = n_distinct(expr_collapsed$gene_name) == nrow(expr_collapsed))

# 7. Save annotated counts
write.csv(expr_collapsed, rna_counts_path, row.names = FALSE, quote = FALSE)
message("  [OUT] Mapped counts saved to ", rna_counts_path)

# 8. Save gene mapper lookup table
gene_mapper <- expr_ann %>% select(ensg, gene_name, gene_type)
write.csv(gene_mapper, gene_mapper_path, row.names = FALSE, quote = FALSE)
message("  [OUT] Gene mapper saved to ", gene_mapper_path)

message("\n01_gene_mapping.R completed successfully.")
