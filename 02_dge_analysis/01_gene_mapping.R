# [[02_dge_analysis/01_gene_mapping.R]]
# Goal: Map Ensembl IDs to HGNC symbols and filter for valid biotypes using BioMart.
#
# Pipeline steps:
#   1. Load raw count matrix (from featureCounts output)
#   2. Clean Ensembl IDs by stripping version suffixes (e.g. .12) and collapse duplicates
#   3. Connect to Ensembl BioMart using Mirror = www
#   4. Query BioMart for HGNC symbol (gene_name) and Biotype (gene_type)
#   5. Deduplicate BioMart query results by keeping the first hit per Ensembl ID
#   6. Join annotations onto the expression count matrix
#   7. Filter the matrix strictly for valid biotypes (e.g. protein_coding)
#   8. Collapse duplicate HGNC symbols by sum
#   9. Save the annotated counts to results/tables/01_annotated_counts.csv
#  10. Generate and save mapping table to results/tables/01_gene_mapper.csv
#

suppressPackageStartupMessages({
  library(biomaRt)
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
cfg <- yaml::read_yaml(config_path)
valid_biotypes <- cfg$ensembl_map$valid_biotypes

# Paths
annotated_path <- file.path(project_root, "_data/featurecounts/raw_counts.csv")
out_dir <- file.path(project_root, "results/tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

rna_counts_path <- file.path(out_dir, "01_annotated_counts.csv")
gene_mapper_path <- file.path(out_dir, "01_gene_mapper.csv")

# Load raw counts
message("[LOAD] Reading raw counts from ", annotated_path)
expr_raw <- read_csv(annotated_path, show_col_types = FALSE)
message("  Loaded matrix: ", nrow(expr_raw), " rows x ", ncol(expr_raw), " columns")

# Strip version suffix and collapse duplicates if needed
expr_raw$gene_id <- expr_raw$Geneid
id_map <- expr_raw %>%
  select(gene_id) %>%
  mutate(gene_id_clean = sub("\\..*$", "", gene_id))

expr <- expr_raw %>%
  mutate(gene_id = sub("\\..*$", "", gene_id)) %>%
  select(-Geneid)
rm(expr_raw)

n_unique_stripped <- n_distinct(expr$gene_id)
n_duplicates <- nrow(expr) - n_unique_stripped

if (n_duplicates == 0) {
  message("  [STRIP] Unique gene_id = ", n_unique_stripped, " — no duplicates, skipping collapse.")
} else {
  message("  [STRIP] Unique gene_id = ", n_unique_stripped, " | duplicate rows to collapse = ", n_duplicates, " — collapsing by sum")
  expr <- expr %>%
    group_by(gene_id) %>%
    summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE)), .groups = "drop")
}

# Connect to Ensembl BioMart
message("[ENSEMBL] Connecting to hsapiens_gene_ensembl (mirror = www)...")
mart <- useEnsembl(
  biomart = "genes",
  dataset = "hsapiens_gene_ensembl",
  mirror = "www"
)

# Query BioMart
ensembl_ids <- unique(expr$gene_id)
message("[BIOMART] Querying Ensembl for HGNC symbols and biotypes...")
bm_raw <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol", "gene_biotype"),
  filters = "ensembl_gene_id",
  values = ensembl_ids,
  mart = mart,
  uniqueRows = TRUE
) %>%
  rename(
    gene_id = ensembl_gene_id,
    gene_name = hgnc_symbol,
    gene_type = gene_biotype
  )

n_bm_rows <- nrow(bm_raw)
n_bm_unique_id <- n_distinct(bm_raw$gene_id)
message("  BioMart rows returned: ", n_bm_rows, " | unique gene_ids: ", n_bm_unique_id)

# Retain first hit per gene_id
bm_first <- bm_raw %>% distinct(gene_id, .keep_all = TRUE)

# Join onto expression matrix
expr_ann <- expr %>% left_join(bm_first, by = "gene_id")
rm(expr)

# Filter by biotype
message("[BIOTYPE] Filtering for biotypes: ", paste(valid_biotypes, collapse = ", "))
expr_filt <- expr_ann %>%
  filter(
    !is.na(gene_name), gene_name != "",
    gene_type %in% valid_biotypes
  )
message("  Rows remaining after biotype filter: ", nrow(expr_filt))

# Collapse duplicate gene names
n_rows <- nrow(expr_filt)
n_unique_symbols <- n_distinct(expr_filt$gene_name)

if (n_rows == n_unique_symbols) {
  message("  No duplicate HGNC symbols found, skipping collapse.")
  expr_collapsed <- expr_filt
} else {
  message("  Collapsing duplicate HGNC symbols (", n_rows - n_unique_symbols, " duplicates) by sum...")
  expr_collapsed <- expr_filt %>%
    select(-gene_id, -gene_type) %>%
    group_by(gene_name) %>%
    summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE)), .groups = "drop")
}

# Save Outputs
write.csv(expr_collapsed, rna_counts_path, row.names = FALSE, quote = FALSE)
message("  [OUT] Mapped counts saved to ", rna_counts_path)

# Generate and save gene mapper mapping file
gene_mapper <- id_map %>%
  left_join(
    expr_filt %>% select(gene_id, gene_name, gene_type) %>% distinct(),
    by = c("gene_id_clean" = "gene_id")
  ) %>%
  filter(!is.na(gene_name), gene_name != "") %>%
  select(gene_id = gene_id_clean, gene_name, gene_type) %>%
  distinct(gene_name, .keep_all = TRUE)

write.csv(gene_mapper, gene_mapper_path, row.names = FALSE, quote = FALSE)
message("  [OUT] Gene mapper saved to ", gene_mapper_path)

message("\n01_gene_mapping.R completed successfully.")
