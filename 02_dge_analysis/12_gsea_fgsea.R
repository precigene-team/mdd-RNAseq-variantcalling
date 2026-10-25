# [[02_dge_analysis/12_gsea_fgsea.R]]
# Goal: Run preranked GSEA (fgsea) on the two input contrasts against MSigDB
#       Hallmark sets. Tables only; plotting lives in 13_gsea_figures.R.
#
# Steps:
#   1. Read the .rnk files written by 11_gsea_input.R (Wald stat ranking).
#   2. Load hallmark.gmt as the gene set collection.
#   3. Run fgseaMultilevel per contrast (minSize 15, maxSize 500, 1000 permutations).
#   4. Save 12_<contrast>_gsea_hallmark.csv to results/02_dge_analysis/tables/.
#   5. Print the top pathways per contrast so the null is visible immediately.

suppressPackageStartupMessages({
  library(fgsea)
  library(readr)
  library(dplyr)
})

set.seed(42)

get_project_root <- function() {
  curr <- getwd()
  while (curr != "/" && !file.exists(file.path(curr, "run"))) {
    curr <- dirname(curr)
  }
  return(curr)
}
project_root <- get_project_root()

gsea_dir   <- file.path(project_root, "results/02_dge_analysis/gsea")
tables_dir <- file.path(project_root, "results/02_dge_analysis/tables")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

contrasts <- c("deg_male_input", "deg_female_input")
gmt_path  <- file.path(gsea_dir, "hallmark.gmt")
gene_sets <- gmtPathways(gmt_path)

read_rnk <- function(path) {
  lines  <- readLines(path)
  body   <- lines[!startsWith(lines, "#")]
  parts  <- strsplit(body, "\t")
  stats  <- as.numeric(vapply(parts, `[`, character(1), 2))
  names(stats) <- vapply(parts, `[`, character(1), 1)
  stats[order(-stats)]
}

for (contrast in contrasts) {
  rnk_file <- file.path(gsea_dir, paste0(contrast, ".rnk"))
  if (!file.exists(rnk_file)) {
    message("  [Skip] no .rnk for: ", contrast)
    next
  }

  message("\n>>> ", contrast)
  stats <- read_rnk(rnk_file)
  message("  ", length(stats), " ranked genes, ", length(gene_sets), " gene sets")

  res <- fgseaMultilevel(pathways = gene_sets, stats = stats,
                         minSize = 15, maxSize = 500,
                         eps = 0.0, scoreType = "std", nproc = 1)

  res <- res %>%
    rename(pvalue = pval) %>%
    mutate(leadingEdge = sapply(leadingEdge, function(x) paste(x, collapse = ";"))) %>%
    arrange(pvalue)

  write_csv(res %>% select(pathway, size, ES, NES, pvalue, padj,
                           leadingEdge),
            file.path(tables_dir, paste0("12_", contrast, "_gsea_hallmark.csv")))

  top <- res %>%
    mutate(pathway = sub("^HALLMARK_", "", pathway)) %>%
    head(10)
  print(top %>% select(pathway, size, NES, pvalue, padj))

  message("  sets with padj < 0.25: ", sum(res$padj < 0.25, na.rm = TRUE),
          "  |  sets with padj < 0.05: ", sum(res$padj < 0.05, na.rm = TRUE),
          "  |  |NES| > 1.5: ", sum(abs(res$NES) > 1.5, na.rm = TRUE))
}

message("\n12_gsea_fgsea.R completed.")
