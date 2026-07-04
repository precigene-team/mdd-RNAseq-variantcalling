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
#   7. Save all results tables to results/tables/ prefixed with 07_.

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
tables_dir <- file.path(project_root, "results/tables")
mapper_file <- file.path(tables_dir, "01_gene_mapper.csv")

if (!file.exists(mapper_file)) {
  stop("Gene mapper file not found: ", mapper_file)
}

# Load gene mapping (HGNC symbol to Ensembl ID mapping)
message("[INFO] Loading gene mapper...")
gene_mapper <- read_csv(mapper_file, show_col_types = FALSE)

# Retrieve MSigDB Hallmark gene sets (ENSG keys)
message("[INFO] Preparing MSigDB Hallmark gene sets...")
h_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, ensembl_gene)

# List of contrasts to analyze
contrasts <- c(
  "deg_male_input",
  "deg_female_input",
  "deg_sex_interaction_input",
  "deg_male_ip",
  "deg_female_ip",
  "meth_male_interaction",
  "meth_female_interaction",
  "meth_sex_interaction"
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
  
  # Join Ensembl IDs from the gene mapper
  res_df <- res_df %>%
    left_join(gene_mapper, by = "gene_name") %>%
    filter(!is.na(gene_id))
  
  # 1. Prepare Significant Genes for ORA (pvalue < 0.05 and |log2FoldChange| >= 0.2)
  sig_df <- res_df %>%
    filter(!is.na(pvalue), pvalue < 0.05, abs(log2FoldChange) >= 0.2)
  
  sig_ensembl <- unique(sig_df$gene_id)
  
  if (length(sig_ensembl) > 5) {
    message("  Running Over-Representation Analysis (ORA) on ", length(sig_ensembl), " genes...")
    
    # Map Ensembl to Entrez ID for KEGG
    gene_map <- tryCatch({
      bitr(sig_ensembl, fromType = "ENSEMBL", toType = c("ENTREZID", "SYMBOL"), OrgDb = org.Hs.eg.db)
    }, error = function(e) return(NULL))
    
    # A. GO Biological Process ORA
    ego <- tryCatch({
      enrichGO(
        gene = sig_ensembl,
        OrgDb = org.Hs.eg.db,
        keyType = "ENSEMBL",
        ont = "BP",
        pAdjustMethod = "BH",
        readable = TRUE
      )
    }, error = function(e) {
      message("    [CAUTION] GO ORA failed: ", e$message)
      return(NULL)
    })
    
    if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
      write.csv(as.data.frame(ego), file = file.path(tables_dir, paste0("07_", contrast, "_go_ora.csv")), row.names = FALSE)
      message("    [SUCCESS] GO ORA results saved.")
    }
    
    # B. KEGG ORA
    if (!is.null(gene_map)) {
      ekegg <- tryCatch({
        enrichKEGG(
          gene = gene_map$ENTREZID,
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
    message("  [INFO] Too few significant genes (", length(sig_ensembl), ") for ORA. Skipping.")
  }
  
  # 2. Prepare Ranked List for GSEA (Using Wald-statistic 'stat')
  res_clean <- res_df %>%
    filter(!is.na(stat), !is.na(gene_id)) %>%
    distinct(gene_id, .keep_all = TRUE)
  
  ranked_genes <- res_clean$stat
  names(ranked_genes) <- res_clean$gene_id
  ranked_genes <- sort(ranked_genes, decreasing = TRUE)
  
  if (length(ranked_genes) > 10) {
    message("  Running Gene Set Enrichment Analysis (GSEA)...")
    
    # A. GSEA GO Biological Process
    ego_gsea <- tryCatch({
      gseGO(
        geneList = ranked_genes,
        OrgDb = org.Hs.eg.db,
        keyType = "ENSEMBL",
        ont = "BP",
        minGSSize = 10,
        maxGSSize = 500,
        pvalueCutoff = 0.05,
        verbose = FALSE
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
