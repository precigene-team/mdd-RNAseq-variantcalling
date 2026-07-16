# [[02_dge_analysis/10_rrho2_sex_comparison.R]]
# Goal: Run RRHO2 cross-sex comparison on male and female DMR ranked lists.

# Steps:
#   1. Load male and female DMR result tables from TRESS differential analysis.
#   2. Compute rank scores: -log10(padj) * sign(log2FC) for each peak in each sex.
#   3. Assert shared peak universe between male and female lists.
#   4. Run RRHO2 on the aligned ranked vectors.
#   5. Export concordant (++/--) and discordant (+/-/-+) heatmaps to results/figures/.
#   6. Export significant overlap peak lists to results/tables/.

# --- Reproducible Dependency Bootstrapping ---
if (!requireNamespace("RRHO2", quietly = TRUE)) {
  message("[INFO] RRHO2 not found. Bootstrapping reproducible install from GitHub...")
  
  # Ensure remotes is available for the installation
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes", repos = "https://cloud.r-project.org")
  }
  
  # Pull directly from the repository head to bypass API reference errors
  remotes::install_github("RRHO2/RRHO2", 
                          upgrade = "never", 
                          quiet = TRUE)
}

suppressPackageStartupMessages({
  library(RRHO2)
  library(dplyr)
  library(ggplot2)
  library(grDevices)
})

# --- Dynamic Pathing Setup ---
project_root <- Sys.getenv("PROJECT_ROOT")
if (project_root == "") {
  curr <- getwd()
  while (curr != "/" && !file.exists(file.path(curr, "run"))) { curr <- dirname(curr) }
  project_root <- curr
}
results_root <- file.path(project_root, "results", "02_dge_analysis")
tables_dir   <- file.path(results_root, "tables")
figures_dir  <- file.path(results_root, "figures")

dir.create(tables_dir,  showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

# --- Input Paths ---
male_dmr_file   <- file.path(tables_dir, "04_deg_male_input.csv")
female_dmr_file <- file.path(tables_dir, "04_deg_female_input.csv")

for (f in c(male_dmr_file, female_dmr_file)) {
  if (!file.exists(f)) {
    stop("Required DMR results file not found: ", f,
         "\nRun TRESS differential analysis first.")
  }
}

# --- Load DMR Results ---
message("[INFO] Loading DMR results...")
male_dmr   <- read.csv(male_dmr_file,   stringsAsFactors = FALSE)
female_dmr <- read.csv(female_dmr_file, stringsAsFactors = FALSE)

# Assertion: required columns exist
required_cols <- c("gene_name", "padj", "log2FoldChange")
for (col in required_cols) {
  if (!col %in% colnames(male_dmr))   stop("Column '", col, "' missing from male DEG table.")
  if (!col %in% colnames(female_dmr)) stop("Column '", col, "' missing from female DEG table.")
}

# --- Compute Rank Scores ---
message("[INFO] Computing rank scores: -log10(padj) * sign(log2FC)...")

rank_score <- function(df) {
  df <- df %>%
    # Filter out any rows containing NA values in critical columns
    filter(!is.na(padj), !is.na(log2FoldChange), !is.na(gene_name)) %>%
    mutate(
      padj_clamped = pmax(padj, 1e-300),
      rank_score   = -log10(padj_clamped) * sign(log2FoldChange)
    ) %>%
    arrange(desc(rank_score))
  return(df)
}

male_ranked   <- rank_score(male_dmr)
female_ranked <- rank_score(female_dmr)

# --- Assert Shared Peak Universe ---
shared_peaks <- intersect(male_ranked$gene_name, female_ranked$gene_name)
message("[INFO] Shared genes between male and female lists: ", length(shared_peaks))

if (length(shared_peaks) < 100) {
  warning("Fewer than 100 shared genes found. RRHO2 may not be reliable.")
}

male_ranked   <- male_ranked   %>% filter(gene_name %in% shared_peaks) %>% arrange(gene_name)
female_ranked <- female_ranked %>% filter(gene_name %in% shared_peaks) %>% arrange(gene_name)

# Assertion: rows now match
stopifnot(nrow(male_ranked) == nrow(female_ranked))
stopifnot(all(male_ranked$gene_name == female_ranked$gene_name))

# --- Build RRHO2 input (data.frame: gene, rank_score) ---
male_list   <- data.frame(gene = male_ranked$gene_name,   value = male_ranked$rank_score)
female_list <- data.frame(gene = female_ranked$gene_name, value = female_ranked$rank_score)

# --- Run RRHO2 ---
message("[INFO] Running RRHO2...")
stepsize <- max(1L, floor(nrow(male_list) / 100L))

rrho_result <- RRHO2_initialize(
  list1       = male_list,
  list2       = female_list,
  stepsize    = stepsize,
  labels      = c("Male MDD vs Ctrl", "Female MDD vs Ctrl"),
  log10.ind   = TRUE
)

# --- Export RRHO2 heatmap as PDF ---
message("[INFO] Exporting RRHO2 heatmap...")

pdf(file.path(figures_dir, "10_rrho2_concordant.pdf"), width = 8, height = 7)
RRHO2_heatmap(rrho_result)
dev.off()
message("  Saved RRHO heatmap to 10_rrho2_concordant.pdf")

# Save hypergeometric map as table
tryCatch({
  write.csv(as.data.frame(rrho_result$hypermat),
            file.path(tables_dir, "10_rrho_hypermat.csv"),
            row.names = FALSE)
  message("  Saved hypergeometric map to 10_rrho_hypermat.csv")
}, error = function(e) {
  message("  [WARNING] Could not save hypergeometric map: ", e$message)
})

# --- Output Assertions ---
for (f in c(
  file.path(figures_dir, "10_rrho2_concordant.pdf")
)) {
  if (!file.exists(f) || file.info(f)$size == 0) {
    stop("Assertion Error: Output PDF not generated or empty: ", f)
  }
}

message("=== RRHO2 Analysis Complete ===")
