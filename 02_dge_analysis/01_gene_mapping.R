# [[02_dge_analysis/01_gene_mapping.R]]
# Goal: Map TRESS consensus peak IDs to HGNC gene symbols using BioMart genomic
#       coordinate overlap, then filter for valid biotypes and collapse by gene.
#
# Pipeline steps:
#   1. Load raw peak-level count matrix (from featureCounts / consensus_peaks.saf)
#   2. Read peak genomic coordinates from the SAF annotation file
#   3. Connect to Ensembl BioMart
#   4. Query BioMart with chromosomal coordinates to retrieve overlapping genes
#   5. Deduplicate: for each peak, keep the first protein-coding gene hit
#   6. Join peak→gene mapping onto the expression count matrix
#   7. Filter the matrix strictly for valid biotypes (e.g. protein_coding)
#   8. Collapse duplicate HGNC symbols by sum (multiple peaks per gene)
#   9. Save annotated counts to results/02_dge_analysis/tables/01_annotated_counts.csv
#  10. Generate and save mapping table to results/02_dge_analysis/tables/01_gene_mapper.csv

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
cfg            <- yaml::read_yaml(config_path)
valid_biotypes <- cfg$ensembl_map$valid_biotypes

# Paths
counts_path    <- file.path(project_root, "_data/featurecounts/raw_counts.csv")
saf_path       <- file.path(project_root, "_data/peaks/consensus_peaks.saf")
out_dir        <- file.path(project_root, "results/02_dge_analysis/tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

rna_counts_path <- file.path(out_dir, "01_annotated_counts.csv")
gene_mapper_path <- file.path(out_dir, "01_gene_mapper.csv")

# Assertions
if (!file.exists(counts_path)) stop("Raw counts file not found: ", counts_path)
if (!file.exists(saf_path))    stop("Consensus peaks SAF not found: ", saf_path)

# 1. Load raw peak-level count matrix
message("[LOAD] Reading raw peak counts from ", counts_path)
expr_raw <- read_csv(counts_path, show_col_types = FALSE)
message("  Loaded matrix: ", nrow(expr_raw), " rows x ", ncol(expr_raw), " columns")
colnames(expr_raw)[1] <- "peak_id"

# 2. Load peak coordinates from SAF file
# SAF format: GeneID  Chr  Start  End  Strand
message("[LOAD] Reading peak coordinates from SAF: ", saf_path)
saf <- read.table(saf_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
colnames(saf) <- c("peak_id", "chr", "start", "end", "strand")

# Normalise chromosome names (remove 'chr' prefix if present to match Ensembl)
saf$chr <- sub("^chr", "", saf$chr)

# Keep only peaks present in the count matrix
saf <- saf %>% filter(peak_id %in% expr_raw$peak_id)
message("  SAF peaks retained after filter: ", nrow(saf))

# 3. Connect to Ensembl BioMart
ensembl_version <- cfg$ensembl_map$ensembl_version
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
      message("[ENSEMBL] Attempt 5: Connecting via useMart")
      useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", host = "https://www.ensembl.org")
    }
  }, error = function(e) {
    message("  Connection failed: ", conditionMessage(e))
    NULL
  })
  if (!is.null(mart)) break
}
if (is.null(mart)) stop("Failed to connect to Ensembl BioMart after multiple attempts.")

# 4. Query BioMart in batches using chromosomal coordinates
# BioMart chromosomal_region filter format: "chr:start:end:strand"
# where strand is 1 (forward) or -1 (reverse)
message("[BIOMART] Querying Ensembl by genomic coordinates (batched by chromosome)...")

saf$strand_bm <- ifelse(saf$strand == "+", 1, -1)
saf$region    <- paste(saf$chr, saf$start, saf$end, saf$strand_bm, sep = ":")

chroms <- unique(saf$chr)
bm_list <- vector("list", length(chroms))
names(bm_list) <- chroms

for (ch in chroms) {
  regions <- saf$region[saf$chr == ch]
  bm_chunk <- tryCatch({
    getBM(
      attributes = c("hgnc_symbol", "gene_biotype",
                     "chromosome_name", "start_position", "end_position"),
      filters    = "chromosomal_region",
      values     = regions,
      mart       = mart,
      uniqueRows = TRUE
    ) %>%
      mutate(chromosome_name = as.character(chromosome_name))
  }, error = function(e) {
    message("  [WARNING] BioMart query failed for chr", ch, ": ", e$message)
    NULL
  })
  if (!is.null(bm_chunk) && nrow(bm_chunk) > 0) {
    bm_list[[ch]] <- bm_chunk
  }
}

bm_raw <- bind_rows(bm_list)
message("  BioMart total gene hits returned: ", nrow(bm_raw))

if (nrow(bm_raw) == 0) {
  stop("BioMart returned 0 gene hits from coordinate queries. Check SAF chromosome names and Ensembl version.")
}

# 5. Map peaks to genes via coordinate overlap
# For each peak, find overlapping genes by matching chromosome and coordinate range
bm_raw <- bm_raw %>%
  rename(
    gene_name  = hgnc_symbol,
    gene_type  = gene_biotype,
    chr        = chromosome_name,
    gene_start = start_position,
    gene_end   = end_position
  ) %>%
  filter(gene_name != "", gene_type %in% valid_biotypes)

message("[OVERLAP] Mapping peaks to overlapping protein-coding genes...")
peak_gene_map <- saf %>%
  select(peak_id, chr, start, end) %>%
  left_join(bm_raw, by = "chr", relationship = "many-to-many") %>%
  filter(
    !is.na(gene_name),
    start <= gene_end,
    end   >= gene_start
  ) %>%
  distinct(peak_id, gene_name, gene_type) %>%
  # Keep one gene per peak (first alphabetically to be deterministic)
  group_by(peak_id) %>%
  slice(1) %>%
  ungroup()

n_mapped   <- n_distinct(peak_gene_map$peak_id)
n_total    <- nrow(saf)
n_genes    <- n_distinct(peak_gene_map$gene_name)
message("  Peaks mapped to a gene: ", n_mapped, " / ", n_total,
        " | Unique genes: ", n_genes)

stopifnot("No peaks could be mapped to genes" = n_mapped > 0)

# 6. Join peak→gene mapping onto count matrix
expr_ann <- expr_raw %>%
  left_join(peak_gene_map %>% select(peak_id, gene_name, gene_type), by = "peak_id") %>%
  filter(!is.na(gene_name), gene_name != "")

message("[BIOTYPE] Rows remaining after biotype filter: ", nrow(expr_ann))

# 7. Collapse duplicate gene names by sum (multiple peaks per gene)
n_rows         <- nrow(expr_ann)
n_unique_genes <- n_distinct(expr_ann$gene_name)

if (n_rows == n_unique_genes) {
  message("  No duplicate gene names found, skipping collapse.")
  expr_collapsed <- expr_ann %>% select(-peak_id, -gene_type) %>% relocate(gene_name)
} else {
  message("  Collapsing ", n_rows - n_unique_genes, " duplicate gene entries by sum...")
  expr_collapsed <- expr_ann %>%
    select(-peak_id, -gene_type) %>%
    group_by(gene_name) %>%
    summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE)), .groups = "drop")
}

# Assertion: no duplicate gene names in output
stopifnot("Duplicate gene names remain after collapse" = n_distinct(expr_collapsed$gene_name) == nrow(expr_collapsed))

# 8. Save annotated counts
write.csv(expr_collapsed, rna_counts_path, row.names = FALSE, quote = FALSE)
message("  [OUT] Mapped counts saved to ", rna_counts_path)

# 9. Save gene mapper lookup table
gene_mapper <- peak_gene_map %>%
  select(peak_id, gene_name, gene_type) %>%
  distinct(gene_name, .keep_all = TRUE)

write.csv(gene_mapper, gene_mapper_path, row.names = FALSE, quote = FALSE)
message("  [OUT] Gene mapper saved to ", gene_mapper_path)

message("\n01_gene_mapping.R completed successfully.")
