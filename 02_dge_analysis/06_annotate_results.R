# [[02_dge_analysis/06_annotate_results.R]]
# Goal: Map significant DEG gene names to Ensembl IDs, biotypes, genomic coordinates,
#       and export SRA run accessions linked to each Input (RNA-seq) cohort subgroup.
#
# Steps:
#   1. Load metadata from sample_sheet.csv and export SRA accession lists for the Input subgroups.
#   2. Load significant DGE results from results/02_dge_analysis/tables/04_*_sig.csv
#   3. Connect to Ensembl BioMart to fetch genomic coordinate and biotype annotations.
#   4. Join annotations to each significant DGE table.
#   5. Save the annotated significant tables as results/02_dge_analysis/tables/04_*_sig_annotated.csv

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

# Directories
tables_dir <- file.path(project_root, "results/02_dge_analysis/tables")
meta_file <- file.path(project_root, "_data/_metadata/sample_sheet.csv")

if (!file.exists(meta_file)) {
  stop("Metadata file not found: ", meta_file)
}

# 1. Export SRA Accessions per Input Cohort Subgroup
message("[INFO] Exporting SRA run accessions per Input subgroup...")
metadata <- read_csv(meta_file, show_col_types = FALSE) %>%
  filter(DataType == "RNA-seq")

# Generate list of accessions for each combination of Gender and SampleType
subgroups <- metadata %>%
  group_by(Gender, SampleType) %>%
  summarise(Accessions = paste(sort(unique(Accession)), collapse = ", "), .groups = "drop")

# Save accession lists to a text file
accessions_txt_file <- file.path(tables_dir, "06_subgroup_sra_accessions.txt")
write_tsv(subgroups, accessions_txt_file)
message("  [SUCCESS] SRA accession list written to ", accessions_txt_file)

# 2. Get list of significant tables to annotate
sig_files <- list.files(tables_dir, pattern = "^04_.*_sig\\.csv$", full.names = TRUE)

if (length(sig_files) == 0) {
  stop("No significant results tables (04_*_sig.csv) found in ", tables_dir)
}

# Extract all unique significant gene names across all tables to minimize BioMart query size
message("[INFO] Reading significant gene lists...")
all_sig_genes <- c()
for (f in sig_files) {
  df <- read_csv(f, show_col_types = FALSE)
  if (nrow(df) > 0) {
    all_sig_genes <- unique(c(all_sig_genes, df$gene_name))
  }
}

if (length(all_sig_genes) == 0) {
  message("[INFO] No significant genes detected across any tables. Skipping BioMart annotation.")
  q(save = "no")
}

message("  Total unique significant genes to annotate: ", length(all_sig_genes))

# 3. Connect to Ensembl BioMart
config_path <- file.path(project_root, "02_dge_analysis/config.yaml")
ensembl_version <- NULL
if (file.exists(config_path)) {
  cfg <- yaml::read_yaml(config_path)
  ensembl_version <- cfg$ensembl_map$ensembl_version
}

mart <- NULL
for (attempt in 1:5) {
  mart <- tryCatch({
    if (attempt == 1 && !is.null(ensembl_version)) {
      message("[ENSEMBL] Attempt 1: Connecting via useEnsembl with version = ", ensembl_version)
      useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", version = ensembl_version)
    } else if (attempt == 2) {
      message("[ENSEMBL] Attempt 2: Connecting via useEnsembl with mirror = 'useast'")
      useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", mirror = "useast")
    } else if (attempt == 3) {
      message("[ENSEMBL] Attempt 3: Connecting via useEnsembl with mirror = 'uswest'")
      useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", mirror = "uswest")
    } else if (attempt == 4) {
      message("[ENSEMBL] Attempt 4: Connecting via useEnsembl with mirror = 'www'")
      useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", mirror = "www")
    } else {
      message("[ENSEMBL] Attempt 5: Connecting via useMart with host = 'https://www.ensembl.org'")
      useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", host = "https://www.ensembl.org")
    }
  }, error = function(e) {
    message("  Connection failed: ", conditionMessage(e))
    NULL
  })
  if (!is.null(mart)) break
}

if (is.null(mart)) {
  stop("Failed to connect to Ensembl BioMart after multiple attempts.")
}

# Query genomic coordinates and biotypes
message("[BIOMART] Querying genomic coordinates and biotypes...")
bm_annot <- getBM(
  attributes = c("hgnc_symbol", "ensembl_gene_id", "gene_biotype", "chromosome_name", "start_position", "end_position", "strand", "description"),
  filters = "hgnc_symbol",
  values = all_sig_genes,
  mart = mart,
  uniqueRows = TRUE
) %>%
  rename(
    gene_name = hgnc_symbol,
    ensembl_id = ensembl_gene_id,
    gene_type = gene_biotype,
    chr = chromosome_name,
    start = start_position,
    end = end_position
  ) %>%
  mutate(
    # Format strand values as standard + / -
    strand = ifelse(strand == 1, "+", "-"),
    # Strip the Ensembl source brackets from the end of the description
    description = sub(" \\[Source:.*$", "", description)
  )

# Keep the first unique hit for each HGNC symbol to avoid duplicating rows in the join
bm_annot_uniq <- bm_annot %>% distinct(gene_name, .keep_all = TRUE)

# 4. Annotate and save tables
for (f in sig_files) {
  output_name <- sub("^04_", "06_", basename(f))
  df <- read_csv(f, show_col_types = FALSE)
  
  if (nrow(df) == 0) {
    message("  [INFO] Skipping empty table: ", basename(f))
    next
  }
  
  # Join annotations
  df_annotated <- df %>%
    left_join(bm_annot_uniq, by = "gene_name") %>%
    # Reorder columns to make annotations prominent
    select(gene_name, ensembl_id, gene_type, chr, start, end, strand, description, everything())
  
  out_file <- file.path(tables_dir, output_name)
  write_csv(df_annotated, out_file)
  message("  [SUCCESS] Written annotated table: ", basename(out_file))
}

message("\n06_annotate_results.R completed successfully.")
