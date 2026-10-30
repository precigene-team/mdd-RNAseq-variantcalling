# [[03_variant_calling/13_export_pathway_sets.R]]
# Goal: Export msigdbr gene sets to a flat TSV for the EXPLORATORY burden arm.
#       Network-free (msigdbr ships its data); runs once, cached by 14 and 16.
#
# Why collections are kept SEPARATE rather than lumped:
#   msigdbr C2 merges several databases into one "CURATED" bag and C5 merges the
#   three Gene Ontology branches into one "GO_BP" bag. Collapsing them would
#   (a) mislabel the source of any finding, and (b) corrupt the FDR denominator
#   -- 6,577 curated sets is a different test from 2,148 Reactome sets, and the
#   BH correction must be computed within the collection actually reported.
#   See docs/polygenic_burden_methodology.md section 5.3.
#
# Why this is EXPLORATORY and not the primary panel:
#   Pathway sets were designed for expression enrichment, not variant
#   collapsing, and they partition overlapping gene collections so their FDR is a
#   ranking device rather than a strict error rate. Reported in full, including
#   nulls; supports no claim on its own.
#
# Naming note: letter-suffix convention already present in this module
# (05_combine_gvcfs.sh, 06_genotype_gvcfs.sh) because msigdbr is an R package
# and 16 is Python. One language per script.
#
# Output:
#   13_pathway_sets.tsv   collection, set_name, set_size, genes (';'-joined)
#   Collections: HALLMARK, REACTOME, KEGG, WIKIPATHWAYS, BIOCARTA,
#                GO_BP, GO_CC, GO_MF
#
# Usage:
#   Rscript 03_variant_calling/13_export_pathway_sets.R

suppressPackageStartupMessages({
  library(msigdbr)
  library(data.table)
})

get_project_root <- function() {
  curr <- getwd()
  while (curr != "/" && !file.exists(file.path(curr, "run"))) curr <- dirname(curr)
  curr
}
project_root <- Sys.getenv("PROJECT_ROOT", get_project_root())

out_file <- file.path(project_root, "results/03_variant_calling/tables/13_pathway_sets.tsv")
dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)

# Map each source database to its canonical collection label. Keys are matched
# against msigdbr's gs_subcat for C2/C5 and gs_name prefix for Hallmark.
SUB2LABEL <- list(
  "CP:REACTOME"          = "REACTOME",
  "CP:KEGG_LEGACY"       = "KEGG",
  "CP:KEGG_MEDICUS"      = "KEGG",
  "CP:WIKIPATHWAYS"      = "WIKIPATHWAYS",
  "CP:BIOCARTA"          = "BIOCARTA",
  "GO:BP"                = "GO_BP",
  "GO:CC"                = "GO_CC",
  "GO:MF"                = "GO_MF"
)
# Anything not in the map above (HPO phenotype ontology, PID, etc.) is dropped:
# those are not gene sets suitable for variant collapsing, and mixing them into
# a burden collection would inflate the FDR denominator with unrelated entries.

collapse <- function(d, tag) {
  d <- as.data.table(d)
  d <- d[!is.na(gene_symbol) & nzchar(gene_symbol), ]
  if (!nrow(d)) return(NULL)
  key <- as.character(d$gs_subcat)
  # Blank / NA subcat means the source has no finer partition (Hallmark), so the
  # collection tag applies to the whole set.
  blank <- is.na(key) | !nzchar(trimws(key))
  lab <- rep(NA_character_, length(key))
  lab[blank] <- tag
  lookup <- !blank
  if (any(lookup)) {
    # match() against names(), NOT unlist() on the list: an unmatched key yields
    # a NULL element that unlist silently drops, which would shift every
    # subsequent label onto the wrong set.
    hit <- match(key[lookup], names(SUB2LABEL))
    lab[lookup] <- unlist(SUB2LABEL, use.names = FALSE)[hit]   # NA where unmapped
  }
  keep <- !is.na(lab)
  lab <- lab[keep]
  d <- d[keep]
  # Base assignment avoids data.table's scoping of a same-named local.
  d[["label"]] <- lab
  d[, .(genes = paste(sort(unique(gene_symbol)), collapse = ";"),
        set_size = as.integer(length(unique(gene_symbol)))),
    by = .(collection = label, set_name = gs_name)]
}

rows <- list()

message("[13] fetching Hallmark (H) ...")
h <- tryCatch(msigdbr(species = "Homo sapiens", category = "H"),
              error = function(e) { warning("H failed: ", conditionMessage(e)); NULL })
if (!is.null(h)) {
  x <- collapse(h, "HALLMARK")
  if (!is.null(x)) { rows[["HALLMARK"]] <- x; message("   HALLMARK: ", nrow(x), " sets") }
}

message("[13] fetching curated pathways (C2) ...")
c2 <- tryCatch(msigdbr(species = "Homo sapiens", category = "C2"),
               error = function(e) { warning("C2 failed: ", conditionMessage(e)); NULL })
if (!is.null(c2)) {
  x <- collapse(c2, "OTHER_CURATED")
  if (!is.null(x)) {
    for (lab in unique(x$collection)) {
      sub <- x[collection == lab]
      rows[[lab]] <- sub
      message("   ", lab, ": ", nrow(sub), " sets")
    }
  }
}

message("[13] fetching Gene Ontology (C5) ...")
c5 <- tryCatch(msigdbr(species = "Homo sapiens", category = "C5"),
               error = function(e) { warning("C5 failed: ", conditionMessage(e)); NULL })
if (!is.null(c5)) {
  x <- collapse(c5, "GO_OTHER")
  if (!is.null(x)) {
    for (lab in unique(x$collection)) {
      sub <- x[collection == lab]
      rows[[lab]] <- sub
      message("   ", lab, ": ", nrow(sub), " sets")
    }
  }
}

if (!length(rows)) stop("no msigdbr collection could be read; the exploratory analysis is skipped")

out <- rbindlist(rows, use.names = TRUE, fill = TRUE)
out[, set_name := gsub("[^A-Za-z0-9_]", "_", set_name)]

# Fixed order so downstream scripts and figure legends are stable.
ord <- c("HALLMARK", "REACTOME", "KEGG", "WIKIPATHWAYS", "BIOCARTA",
         "OTHER_CURATED", "GO_BP", "GO_CC", "GO_MF")
out[, collection := factor(collection, levels = c(ord, setdiff(levels(collection), ord)))]
setorder(out, collection, set_size, set_name)
out[, collection := as.character(collection)]

fwrite(out[, .(collection, set_name, set_size, genes)], out_file, sep = "\t")

message("[13] wrote ", nrow(out), " sets -> ", out_file)
for (lab in out[, .N, by = collection]$collection) {
  s <- out[collection == lab]
  message(sprintf("   %-14s %5d sets, median size %3d", lab, nrow(s), median(s$set_size)))
}