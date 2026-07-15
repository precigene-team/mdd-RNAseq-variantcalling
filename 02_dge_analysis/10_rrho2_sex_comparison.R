# [[02_dge_analysis/10_rrho2_sex_comparison.R]]
# Goal: Run RRHO2 cross-sex comparison on male and female DMR ranked lists.

# Steps:
#   1. Load male and female DMR result tables from TRESS differential analysis.
#   2. Compute rank scores: -log10(padj) * sign(log2FC) for each peak in each sex.
#   3. Assert shared peak universe between male and female lists.
#   4. Run RRHO2 on the aligned ranked vectors.
#   5. Export concordant (++/--) and discordant (+/-/-+) heatmaps to results/figures/.
#   6. Export significant overlap peak lists to results/tables/.

suppressPackageStartupMessages({
  library(RRHO)
  library(dplyr)
  library(ggplot2)
})

# --- Dynamic Pathing Setup ---
project_root <- Sys.getenv("PROJECT_ROOT")
if (project_root == "") {
  curr <- getwd()
  while (curr != "/" && !file.exists(file.path(curr, "run"))) { curr <- dirname(curr) }
  project_root <- curr
}
results_root <- file.path(project_root, "results")
tables_dir   <- file.path(results_root, "tables")
figures_dir  <- file.path(results_root, "figures")

dir.create(tables_dir,  showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

# --- Input Paths ---
male_dmr_file   <- file.path(tables_dir, "04_male_dmr_results.csv")
female_dmr_file <- file.path(tables_dir, "04_female_dmr_results.csv")

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
required_cols <- c("peak_id", "padj", "log2FC")
for (col in required_cols) {
  if (!col %in% colnames(male_dmr))   stop("Column '", col, "' missing from male DMR table.")
  if (!col %in% colnames(female_dmr)) stop("Column '", col, "' missing from female DMR table.")
}

# --- Compute Rank Scores ---
message("[INFO] Computing rank scores: -log10(padj) * sign(log2FC)...")

rank_score <- function(df) {
  df <- df %>%
    mutate(
      padj_clamped = pmax(padj, 1e-300),
      rank_score   = -log10(padj_clamped) * sign(log2FC)
    ) %>%
    arrange(desc(rank_score))
  return(df)
}

male_ranked   <- rank_score(male_dmr)
female_ranked <- rank_score(female_dmr)

# --- Assert Shared Peak Universe ---
shared_peaks <- intersect(male_ranked$peak_id, female_ranked$peak_id)
message("[INFO] Shared peaks between male and female lists: ", length(shared_peaks))

if (length(shared_peaks) < 100) {
  warning("Fewer than 100 shared peaks found. RRHO2 may not be reliable.")
}

male_ranked   <- male_ranked   %>% filter(peak_id %in% shared_peaks) %>% arrange(peak_id)
female_ranked <- female_ranked %>% filter(peak_id %in% shared_peaks) %>% arrange(peak_id)

# Assertion: rows now match
stopifnot(nrow(male_ranked) == nrow(female_ranked))
stopifnot(all(male_ranked$peak_id == female_ranked$peak_id))

# --- Build RRHO2 input (data.frame: gene, rank_score) ---
male_list   <- data.frame(gene = male_ranked$peak_id,   value = male_ranked$rank_score)
female_list <- data.frame(gene = female_ranked$peak_id, value = female_ranked$rank_score)

# --- Run RRHO2 ---
message("[INFO] Running RRHO2...")
rrho_result <- RRHO2_initialize(
  list1    = male_list,
  list2    = female_list,
  labels   = c("Male MDD vs Ctrl", "Female MDD vs Ctrl"),
  log10.ind = TRUE
)

# --- Export Heatmaps ---
message("[INFO] Exporting concordant heatmap...")
pdf(file.path(figures_dir, "10_rrho2_concordant.pdf"), width = 8, height = 7)
RRHO2_heatmap(rrho_result, type = "uu")  # upregulated in both
RRHO2_heatmap(rrho_result, type = "dd")  # downregulated in both
dev.off()

message("[INFO] Exporting discordant heatmap...")
pdf(file.path(figures_dir, "10_rrho2_discordant.pdf"), width = 8, height = 7)
RRHO2_heatmap(rrho_result, type = "ud")  # up male, down female
RRHO2_heatmap(rrho_result, type = "du")  # down male, up female
dev.off()

# --- Export Significant Overlap Peak Lists ---
# RRHO2 stores gene lists in the result object for each quadrant
for (quad in c("uu", "dd", "ud", "du")) {
  tryCatch({
    peak_list <- rrho_result[[quad]][["genelist"]]
    if (!is.null(peak_list) && length(peak_list) > 0) {
      out_file <- file.path(tables_dir, paste0("10_rrho2_", quad, "_peaks.csv"))
      write.csv(data.frame(peak_id = peak_list), out_file, row.names = FALSE)
      message("  Saved ", quad, " peaks (n=", length(peak_list), ") to ", out_file)
    }
  }, error = function(e) {
    message("  [WARNING] Could not extract ", quad, " gene list: ", e$message)
  })
}

# --- Output Assertions ---
for (f in c(
  file.path(figures_dir, "10_rrho2_concordant.pdf"),
  file.path(figures_dir, "10_rrho2_discordant.pdf")
)) {
  if (!file.exists(f) || file.info(f)$size == 0) {
    stop("Assertion Error: Output PDF not generated or empty: ", f)
  }
}

message("=== RRHO2 Analysis Complete ===")
