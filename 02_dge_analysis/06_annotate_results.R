# [[02_dge_analysis/06_annotate_results.R]]
# Goal: Map significant DEG gene names to Ensembl IDs, biotypes, genomic coordinates,
#       and export SRA run accessions linked to each Input (RNA-seq) cohort subgroup.
#
# Steps:
#   1. Load metadata from sample_sheet.csv and export SRA accession lists for the Input subgroups.
#   2. Load significant DGE results from results/02_dge_analysis/tables/04_*_sig.csv
#   3. Fetch genomic coordinate and biotype annotations from the Ensembl REST gene lookup
#      (/lookup/symbol). Ensembl's legacy BioMart archive service is down, so biomaRt
#      cannot connect to any mirror; the REST endpoint answers.
#   4. Join annotations to each significant DGE table.
#   5. Save the annotated significant tables as results/02_dge_analysis/tables/04_*_sig_annotated.csv

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(jsonlite)
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

# Extract all unique significant gene names across all tables (one REST lookup per gene)
message("[INFO] Reading significant gene lists...")
all_sig_genes <- c()
for (f in sig_files) {
  df <- read_csv(f, show_col_types = FALSE)
  if (nrow(df) > 0) {
    all_sig_genes <- unique(c(all_sig_genes, df$gene_name))
  }
}

if (length(all_sig_genes) == 0) {
  message("[INFO] No significant genes detected across any tables. Skipping REST annotation.")
  q(save = "no")
}

message("  Total unique significant genes to annotate: ", length(all_sig_genes))

# 3. Fetch genomic coordinates and biotypes from Ensembl REST
#    One lookup per significant gene: /lookup/symbol/homo_sapiens/<symbol>
#    Ensembl REST rate limit is 15 requests/second, so requests are spaced.
rest_base <- "https://rest.ensembl.org/lookup/symbol/homo_sapiens/"
message("[REST] Querying Ensembl REST /lookup/symbol for genomic coordinates and biotypes...")

bm_annot <- data.frame()
for (g in all_sig_genes) {
  res <- NULL
  for (attempt in 1:3) {
    res <- tryCatch(
      jsonlite::fromJSON(paste0(rest_base, g)),
      error = function(e) {
        message("  [WARN] REST lookup for ", g, " (attempt ", attempt, "): ", conditionMessage(e))
        NULL
      }
    )
    if (!is.null(res) && !is.null(res$id)) break
    if (is.null(res$id)) {
      message("  [WARN] REST lookup for ", g, " (attempt ", attempt, "): no gene returned")
      res <- NULL
    }
    Sys.sleep(2)
  }
  if (is.null(res)) {
    message("  [WARN] No annotation returned for ", g)
    next
  }

  bm_annot <- rbind(bm_annot, data.frame(
    gene_name   = g,
    ensembl_id  = res$id,
    gene_type   = res$biotype,
    chr         = res$seq_region_name,
    start       = res$start,
    end         = res$end,
    # Format strand values as standard + / -
    strand      = ifelse(res$strand == 1, "+", "-"),
    # Strip the Ensembl source brackets from the end of the description
    description = if (is.null(res$description)) NA_character_ else
                    sub(" \\[Source:.*$", "", res$description),
    stringsAsFactors = FALSE
  ))
  Sys.sleep(0.1)
}

if (nrow(bm_annot) == 0) {
  stop("Failed to fetch annotations from Ensembl REST for ", length(all_sig_genes), " gene(s).")
}

message("  [SUCCESS] Annotated ", nrow(bm_annot), " of ", length(all_sig_genes), " significant genes")

# Keep the first unique hit for each gene name to avoid duplicating rows in the join
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
