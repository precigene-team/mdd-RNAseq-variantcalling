# [[03_variant_calling/16_burden_control_vs_mdd_table.R]]
# Goal: Build the per-donor and per-set burden tables consumed by the 17_*
#       plotting scripts. TABLES ONLY -- draws nothing.
#
# Why the collections are kept SEPARATE:
#   17 computes BH-FDR WITHIN collection, so the denominator differs per tool by
#   two orders of magnitude (HALLMARK 50 sets vs GO_BP 5,831). Each 17_* script
#   draws one collection and states that collection's own denominator; pooling
#   them would present sets tested under different corrections as comparable.
#
# Reads:
#   11_burden_matrix.tsv.gz  donor x variant carriage with gene annotation
#   11_donor_totals.csv      donor, condition, sex, total variants (the offset)
#   13_pathway_sets.tsv     msigdbr sets, collection-labelled
#   14_pathway_burden.csv    statistics from 14 (optional; attached if present)
#
# Outputs:
#   16_burden_per_donor.csv  one row per donor per set, all collections
#   16_burden_summary.csv    group medians + p + FDR per set
#
# Usage:
#   Rscript 03_variant_calling/16_burden_control_vs_mdd_table.R

suppressPackageStartupMessages(library(data.table))

get_project_root <- function() {
  curr <- getwd()
  while (curr != "/" && !file.exists(file.path(curr, "run"))) curr <- dirname(curr)
  curr
}
project_root <- Sys.getenv("PROJECT_ROOT", get_project_root())

TAB <- file.path(project_root, "results/03_variant_calling/tables")
need <- function(f) {
  p <- file.path(TAB, f)
  if (!file.exists(p)) stop("missing input: ", p, "\nrun 15 and 16b first")
  fread(p)
}

mat  <- need("11_burden_matrix.tsv.gz")
tot  <- need("11_donor_totals.csv")
sets <- need("13_pathway_sets.tsv")
stats <- tryCatch(need("14_pathway_burden.csv"), error = function(e) NULL)

mat[, condition := tot$condition[match(donor, tot$donor)]]
genes_all <- unique(mat$gene)
message("[16] matrix: ", nrow(mat), " carrier rows, ", length(unique(mat$donor)), " donors")

# --- Set -> gene membership -------------------------------------------------
# One unnest for all 14k sets, then one join, instead of 14k merges and 14k
# re-splits of the ';'-joined gene strings.
memb <- sets[, .(set_idx = .I, gene = strsplit(genes, ";", fixed = TRUE)[[1]]),
            by = .(collection, set_name, set_size, genes)]
memb <- memb[gene %chin% genes_all]
message("[16] membership rows after cohort intersection: ", nrow(memb))

set_meta <- memb[, .(genes_in_cohort = uniqueN(gene)), by = .(set_idx, collection, set_name, set_size)]
set_meta <- set_meta[genes_in_cohort >= 3]
message("[16] sets with >=3 overlapping genes: ", nrow(set_meta), " of ", nrow(sets))

# --- Per-donor burden --------------------------------------------------------
# carriage is unique on (donor, gene), so summing over a set's genes gives the
# number of distinct variants that donor carries within the set.
carriage <- mat[, .(burden = .N), by = .(donor, gene)]

hits <- merge(memb[set_idx %in% set_meta$set_idx, .(set_idx, gene)],
              carriage, by = "gene", allow.cartesian = TRUE)
message("[16] membership x carriage hits: ", nrow(hits))

per_donor <- merge(hits[, .(burden = sum(burden)), by = .(set_idx, donor)],
                   set_meta[, .(set_idx, collection, set_name, genes_in_cohort)],
                   by = "set_idx")
# Every set must appear for every donor, including donors carrying no variant in
# it (burden 0). Build the full grid from set_meta so those rows keep their
# collection and name -- merging a bare donor x set_idx grid would leave them NA.
grid <- CJ(set_idx = set_meta$set_idx, donor = unique(mat$donor), unique = TRUE)
grid <- merge(grid, set_meta, by = "set_idx")
per_donor <- merge(grid, per_donor[, .(set_idx, donor, burden)],
                   by = c("set_idx", "donor"), all.x = TRUE)
per_donor[is.na(burden), burden := 0]
# sum() on an integer column stays integer for groups with no zero-fill, so the
# column is mixed-type across groups and data.table rejects the median().
per_donor[, burden := as.numeric(burden)]

pd <- merge(per_donor,
            tot[, .(donor, condition, sex, total_variants_detected)],
            by = "donor", all.x = TRUE)
pd[, grp := ifelse(condition == "mdd", "MDD (n=20)", "Control (n=19)")]
pd[, log2rr := log2((median(burden[grp == "MDD (n=20)"]) + 0.5) /
                    (median(burden[grp == "Control (n=19)"]) + 0.5)), by = set_idx]

if (!is.null(stats) && all(c("set_name", "p_glm_offset", "fdr_bh_within_collection") %in% names(stats))) {
  pd <- merge(pd, unique(stats[, .(set_name, p_glm_offset, fdr_bh_within_collection)]),
              by = "set_name", all.x = TRUE)
} else {
  warning("[16] 14_pathway_burden.csv not found or missing columns; p/FDR left NA")
  pd[, `:=`(p_glm_offset = NA_real_, fdr_bh_within_collection = NA_real_)]
}

fwrite(pd[, .(donor, condition, sex, collection, set_name, genes_in_cohort,
              burden, total_variants_detected, p_glm_offset, fdr_bh_within_collection)],
       file.path(TAB, "16_burden_per_donor.csv"))

summary_tab <- pd[, .(
  median_ctrl = median(burden[grp == "Control (n=19)"]),
  median_mdd  = median(burden[grp == "MDD (n=20)"]),
  log2rr      = log2rr[1],
  ratio       = 2^log2rr[1],
  p_glm       = p_glm_offset[1],
  fdr         = fdr_bh_within_collection[1],
  genes_in_cohort = genes_in_cohort[1]
), by = .(collection, set_name)]
setorder(summary_tab, collection, log2rr)
fwrite(summary_tab, file.path(TAB, "16_burden_summary.csv"))

message("[16] wrote 16_burden_per_donor.csv (", nrow(pd), " rows)")
message("[16] wrote 16_burden_summary.csv (", nrow(summary_tab), " sets)")
message("[16] --- per collection (n_sets is the BH-FDR denominator in 17) ---")
for (lab in summary_tab[, .N, by = collection]$collection) {
  s <- summary_tab[collection == lab]
  message(sprintf("   %-14s %5d sets | %d at FDR<0.05",
                  lab, nrow(s), sum(s$fdr < 0.05, na.rm = TRUE)))
}
message("[16] done")