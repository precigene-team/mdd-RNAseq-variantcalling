# [[02_dge_analysis/11_gsea_input.R]]
# Goal: Write GSEA input files only. No enrichment run, no plots.
#
# Rationale: the single-gene tests are null (fold changes hug zero), so the
# next step is threshold-free set analysis on a *ranked* list. This script
# produces exactly what a preranked GSEA run consumes, for every contrast.
#
# Steps:
#   1. Load each 04_*_res_raw.rds and its dds (for phenotype labels).
#   2. Build two ranking metrics per gene:
#        stat      - DESeq2 Wald statistic (what 07_enrichment_analysis.R uses)
#        signed_p  - sign(log2FC) * -log10(pvalue) (what the lung_spatial
#                    preranked GSEA script uses; keeps effect direction)
#   3. Drop NA metrics, zero-count genes, and duplicate symbols.
#   4. Write, per contrast, into results/02_dge_analysis/gsea/:
#        <contrast>.rnk          GSEA desktop preranked, stat metric
#        <contrast>_signedp.rnk  same, signed_p metric
#        <contrast>_ranked.csv   full table (gene, both metrics, lfc, p, padj)
#        <contrast>.cls          phenotype file from dds colData
#   5. Write hallmark.gmt (MSigDB category H) for the desktop app.
#   6. Write gsea_manifest.csv: contrast, genes, samples, class counts.
#
# Run: ~/.conda/envs/mdd/bin/Rscript 02_dge_analysis/11_gsea_input.R

suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(readr)
  library(msigdbr)
})

get_project_root <- function() {
  curr <- getwd()
  while (curr != "/" && !file.exists(file.path(curr, "run"))) {
    curr <- dirname(curr)
  }
  return(curr)
}
project_root <- get_project_root()

rdata_dir <- file.path(project_root, "_data/RData/02_dge_analysis")
gsea_dir  <- file.path(project_root, "results/02_dge_analysis/gsea")
dir.create(gsea_dir, recursive = TRUE, showWarnings = FALSE)

# Input (RNA-seq) contrasts only: male and female MDD vs CTRL.
# No IP/RIP-seq, no sex interaction, no methylation contrasts.
contrasts <- c("deg_male_input", "deg_female_input")

# GSEA desktop .rnk: header line is 3 tab-separated fields, then GENE<TAB>score
write_rnk <- function(df, path, metric_col, header) {
  lines <- c(paste0("#", header, "\tGeneSet\t", metric_col),
             paste(df$gene, df[[metric_col]], sep = "\t"))
  writeLines(lines, path)
}

# GSEA .cls: #, n_samples n_classes, class names, then name<TAB>class
write_cls <- function(samples, classes, path) {
  cls <- factor(classes)
  lines <- c("#",
             paste(length(samples), nlevels(cls)),
             paste(levels(cls), collapse = " "),
             paste(samples, cls, sep = "\t"))
  writeLines(lines, path)
}

manifest <- data.frame()

for (contrast in contrasts) {
  res_file  <- file.path(rdata_dir, paste0("04_", contrast, "_res_raw.rds"))
  dds_file  <- file.path(rdata_dir, paste0("04_", contrast, "_dds.rds"))
  if (!file.exists(res_file)) {
    message("  [Skip] no res_raw for: ", contrast)
    next
  }

  message("\n>>> ", contrast)
  res_df <- as.data.frame(readRDS(res_file))
  res_df$gene <- rownames(res_df)

  # Two metrics. pvalue is capped so -log10 never returns Inf.
  res_df$signed_p <- sign(res_df$log2FoldChange) *
    (-log10(pmax(res_df$pvalue, 1e-300)))

  n_raw <- nrow(res_df)
  res_df <- res_df %>%
    filter(!is.na(stat), !is.na(signed_p),
           !is.na(gene), gene != "", baseMean > 0) %>%
    distinct(gene, .keep_all = TRUE)

  ranked <- res_df %>% arrange(desc(stat))

  write_rnk(ranked, file.path(gsea_dir, paste0(contrast, ".rnk")),
            "stat", paste0(contrast, " Wald statistic"))
  write_rnk(ranked %>% arrange(desc(signed_p)),
            file.path(gsea_dir, paste0(contrast, "_signedp.rnk")),
            "signed_p", paste0(contrast, " signed -log10 p"))

  write_csv(ranked %>% select(gene, stat, signed_p, log2FoldChange,
                              pvalue, padj, baseMean),
            file.path(gsea_dir, paste0(contrast, "_ranked.csv")))

  dds <- readRDS(dds_file)
  write_cls(colnames(dds), dds$SampleType,
            file.path(gsea_dir, paste0(contrast, ".cls")))

  tab <- table(dds$SampleType)
  manifest <- rbind(manifest, data.frame(
    contrast   = contrast,
    genes      = nrow(ranked),
    dropped    = n_raw - nrow(ranked),
    samples    = ncol(dds),
    classes    = paste0(names(tab), "=", tab, collapse = " "),
    max_abs_lfc = round(max(abs(ranked$log2FoldChange), na.rm = TRUE), 3),
    n_p_lt_05  = sum(ranked$pvalue < 0.05, na.rm = TRUE)
  ))
  message("  ", nrow(ranked), " genes, ", ncol(dds), " samples -> ",
          basename(file.path(gsea_dir, paste0(contrast, ".rnk"))))
}

# Hallmark gene sets as .gmt (symbol keys, matching the ranked lists)
h <- msigdbr(species = "Homo sapiens", collection = "H")
gmt <- h %>% group_by(gs_name) %>%
  summarise(genes = paste(gene_symbol, collapse = "\t")) %>%
  mutate(line = paste(gs_name, genes, sep = "\t"))
writeLines(gmt$line, file.path(gsea_dir, "hallmark.gmt"))
message("\nhallmark.gmt: ", nrow(gmt), " sets")

write_csv(manifest, file.path(gsea_dir, "gsea_manifest.csv"))
print(manifest)
message("\n11_gsea_input.R completed. Inputs in: ", gsea_dir)
