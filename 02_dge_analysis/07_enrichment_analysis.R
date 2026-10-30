# [[02_dge_analysis/07_enrichment_analysis.R]]
# Goal: Run functional enrichment analysis (GO Biological Process, KEGG, and MSigDB Hallmark)
#       using clusterProfiler, following the design from the neighboring RNAseq_pipeline.
#
# Steps:
#   1. Load raw DESeq2 results from _data/RData/02_dge_analysis/04_*_res_raw.rds
#   2. Convert Ensembl Gene IDs (or HGNC symbols) to Entrez IDs for KEGG ORA.
#   3. Run GO Over-Representation Analysis (ORA) on significant genes (pvalue < 0.05, |LFC| >= 0.2).
#   4. Run KEGG Over-Representation Analysis (ORA).
#   5. Run GSEA GO Biological Process using Wald-statistic (stat) ranking.
#   6. Run GSEA MSigDB Hallmark using Wald-statistic (stat) ranking.
#   7. Save all results tables to results/02_dge_analysis/tables/ prefixed with 07_.

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(msigdbr)
  library(readr)
  library(dplyr)
  library(stringr)
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
mapper_file <- file.path(tables_dir, "01_gene_mapper.csv")

if (!file.exists(mapper_file)) {
  stop("Gene mapper file not found: ", mapper_file)
}

# Load gene mapping (ensg → gene_name, gene_type)
message("[INFO] Loading gene mapper...")
gene_mapper <- read_csv(mapper_file, show_col_types = FALSE)

# Retrieve MSigDB Hallmark gene sets (gene symbol keys)
message("[INFO] Preparing MSigDB Hallmark gene sets...")
h_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, gene_symbol)

# List of contrasts to analyze
contrasts <- c(
  "deg_male_input",
  "deg_female_input",
  "deg_sex_interaction_input"
)

message("\n--- Running clusterProfiler Enrichment Analysis ---")

for (contrast in contrasts) {
  res_file <- file.path(rdata_dir, paste0("04_", contrast, "_res_raw.rds"))
  if (!file.exists(res_file)) {
    message("  [Skip] Raw results RDS not found for: ", contrast)
    next
  }
  
  message("\n>>> Processing contrast: ", contrast)
  res_obj <- readRDS(res_file)
  res_df <- as.data.frame(res_obj)
  res_df$gene_name <- rownames(res_df)
  
  # Join gene_name; mapper already keyed by gene_name so this is a no-op but
  # ensures consistent structure with future Ensembl-keyed mappers
  res_df <- res_df %>%
    filter(!is.na(gene_name), gene_name != "")
  
  # 1. Prepare Significant Genes for ORA (pvalue < 0.05 and |log2FoldChange| >= 0.2)
  sig_df <- res_df %>%
    filter(!is.na(pvalue), pvalue < 0.05, abs(log2FoldChange) >= 0.2)
  
  sig_symbols <- unique(sig_df$gene_name)
  
  if (length(sig_symbols) > 5) {
    message("  Running Over-Representation Analysis (ORA) on ", length(sig_symbols), " genes...")
    
    # Map SYMBOL to Entrez ID for KEGG
    gene_map <- tryCatch({
      bitr(sig_symbols, fromType = "SYMBOL", toType = c("ENTREZID", "ENSEMBL"), OrgDb = org.Hs.eg.db)
    }, error = function(e) return(NULL))
    
    # A. GO Biological Process ORA (SYMBOL key)
    ego <- tryCatch({
      enrichGO(
        gene          = sig_symbols,
        OrgDb         = org.Hs.eg.db,
        keyType       = "SYMBOL",
        ont           = "BP",
        pAdjustMethod = "BH",
        readable      = FALSE
      )
    }, error = function(e) {
      message("    [CAUTION] GO ORA failed: ", e$message)
      return(NULL)
    })
    
    if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
      write.csv(as.data.frame(ego), file = file.path(tables_dir, paste0("07_", contrast, "_go_ora.csv")), row.names = FALSE)
      message("    [SUCCESS] GO ORA results saved.")
    }
    
    # B. KEGG ORA (Entrez IDs)
    if (!is.null(gene_map)) {
      ekegg <- tryCatch({
        enrichKEGG(
          gene     = gene_map$ENTREZID,
          organism = "hsa"
        )
      }, error = function(e) {
        message("    [CAUTION] KEGG ORA failed: ", e$message)
        return(NULL)
      })
      
      if (!is.null(ekegg) && nrow(as.data.frame(ekegg)) > 0) {
        write.csv(as.data.frame(ekegg), file = file.path(tables_dir, paste0("07_", contrast, "_kegg_ora.csv")), row.names = FALSE)
        message("    [SUCCESS] KEGG ORA results saved.")
      }
    }
  } else {
    message("  [INFO] Too few significant genes (", length(sig_symbols), ") for ORA. Skipping.")
  }
  
  # 2. Prepare Ranked List for GSEA (Wald-statistic, named by gene symbol)
  res_clean <- res_df %>%
    filter(!is.na(stat), !is.na(gene_name)) %>%
    distinct(gene_name, .keep_all = TRUE)
  
  ranked_genes <- res_clean$stat
  names(ranked_genes) <- res_clean$gene_name
  ranked_genes <- sort(ranked_genes, decreasing = TRUE)
  
  if (length(ranked_genes) > 10) {
    message("  Running Gene Set Enrichment Analysis (GSEA)...")
    
    # A. GSEA GO Biological Process (SYMBOL key)
    ego_gsea <- tryCatch({
      gseGO(
        geneList     = ranked_genes,
        OrgDb        = org.Hs.eg.db,
        keyType      = "SYMBOL",
        ont          = "BP",
        minGSSize    = 10,
        maxGSSize    = 500,
        pvalueCutoff = 0.05,
        verbose      = FALSE
      )
    }, error = function(e) {
      message("    [CAUTION] GSEA GO failed: ", e$message)
      return(NULL)
    })
    
    if (!is.null(ego_gsea) && nrow(as.data.frame(ego_gsea)) > 0) {
      write.csv(as.data.frame(ego_gsea), file = file.path(tables_dir, paste0("07_", contrast, "_gsea_functional.csv")), row.names = FALSE)
      message("    [SUCCESS] GSEA GO results saved.")
    }
    
    # B. GSEA Hallmark
    egsea_hallmark <- tryCatch({
      GSEA(
        ranked_genes,
        TERM2GENE = h_t2g,
        minGSSize = 10,
        maxGSSize = 500,
        pvalueCutoff = 0.05,
        verbose = FALSE
      )
    }, error = function(e) {
      message("    [CAUTION] GSEA Hallmark failed: ", e$message)
      return(NULL)
    })
    
    if (!is.null(egsea_hallmark) && nrow(as.data.frame(egsea_hallmark)) > 0) {
      write.csv(as.data.frame(egsea_hallmark), file = file.path(tables_dir, paste0("07_", contrast, "_gsea_hallmark.csv")), row.names = FALSE)
      message("    [SUCCESS] GSEA Hallmark results saved.")
    }
  } else {
    message("  [INFO] Too few ranked genes for GSEA. Skipping.")
  }
}

message("\n07_enrichment_analysis.R completed successfully.")
