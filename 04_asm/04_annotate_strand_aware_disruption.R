# [[04_asm/04_annotate_strand_aware_disruption.R]]
# Goal: Annotate whether each variant overlapping a DRACH motif abolishes METTL3/14 writer competence.
#
# Steps:
#   1. Load the bedtools -wa -wb variant x DRACH intersection BED (13 columns).
#   2. Map each variant onto its motif offset (0..4), strand-aware.
#   3. Complement alleles for minus-strand motifs via Biostrings.
#   4. Classify DISRUPTS_DRACH vs PRESERVES_DRACH and flag core m6A (offset 2).
#   5. Write the annotated overlap table and a cohort summary.
#
# Notes: IUPAC sets are transcript-oriented, so genomic alleles must be complemented
# first or every minus-strand verdict inverts (2,339 of 4,709 overlaps are minus
# strand). motif_id includes the instance start because motif names are not unique.
# See docs/04_asm_methods.md.

suppressPackageStartupMessages({
  library(Biostrings)
  library(data.table)
})

# --- Dynamic Pathing ---
file_arg   <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg[1])))
} else getwd()
project_root <- Sys.getenv("PROJECT_ROOT",
                           unset = normalizePath(file.path(script_dir, "..")))

intersect_file <- file.path(project_root, "_data", "logs", "drach", "tmp",
                            "snp_drach_intersect.bed")
out_dir   <- file.path(project_root, "results", "04_asm", "tables")
out_table <- file.path(out_dir, "04_drach_variant_overlaps.csv")
out_summary <- file.path(out_dir, "04_drach_overlap_summary.csv")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(intersect_file) || file.info(intersect_file)$size == 0) {
  stop("Assertion Error: intersection BED missing or empty: ", intersect_file,
       "\nRun 03_intersect_drach_genotypes.sh first.")
}

# --- IUPAC (TRANSCRIPT orientation, 5'->3') ---
# D = A/G/T (not C) | R = A/G | A = A (core m6A) | C = C | H = A/C/T (not G)
IUPAC_ALLOWED <- list("0"=c("A","G","T"), "1"=c("A","G"), "2"="A",
                      "3"="C", "4"=c("A","C","T"))

# --- Step 1: load intersection ---
# -a (variant BED): 0 chrom, 1 pos0, 2 pos1, 3 donor, 4 ref, 5 alt
# -b (motif BED)  : 6 chrom, 7 start0, 8 end0, 9 name, 10 score, 11 strand, 12 methyl_pos0
dt <- fread(intersect_file, header = FALSE, sep = "\t", quote = "",
            col.names = c("snp_chr", "snp_start0", "snp_pos", "donor", "ref", "alt",
                          "motif_chr", "motif_start0", "motif_end0", "motif_name",
                          "score", "strand", "methyl_pos0"))
if (nrow(dt) == 0L) stop("Assertion Error: intersection dataframe is empty.")

# --- Step 2: strand-aware motif offset ---
dt[, motif_offset := fifelse(strand == "+",
                             (snp_pos - 1L) - motif_start0,
                             motif_end0 - snp_pos)]
dt <- dt[motif_offset >= 0L & motif_offset <= 4L]

# --- Step 3: complement alleles for minus-strand motifs (Biostrings) ---
# Vectorized via DNAStringSet; NA/empty entries pass through untouched so a
# malformed row cannot abort the whole run.
comp_vec <- function(x) {
  out <- rep(NA_character_, length(x))
  ok <- !is.na(x) & nzchar(x)
  if (any(ok)) out[ok] <- as.character(complement(DNAStringSet(x[ok])))
  out
}
dt[, tx_ref := fifelse(strand == "-", comp_vec(ref), ref)]
dt[, tx_alt := fifelse(strand == "-", comp_vec(alt), alt)]
dt[, ref_allowed := mapply(function(o, b) b %in% IUPAC_ALLOWED[[as.character(o)]],
                           motif_offset, tx_ref, USE.NAMES = FALSE)]
dt[, alt_allowed := mapply(function(o, b) b %in% IUPAC_ALLOWED[[as.character(o)]],
                           motif_offset, tx_alt, USE.NAMES = FALSE)]

# --- Step 4: classify ---
dt[, disruption_status := fifelse(alt_allowed, "PRESERVES_DRACH", "DISRUPTS_DRACH")]
dt[, is_core_m6a := motif_offset == 2L]

# Invariant: any variant on the core adenine must abolish the motif
core_bad <- dt[is_core_m6a == TRUE & disruption_status != "DISRUPTS_DRACH"]
if (nrow(core_bad) > 0L) {
  stop("Invariant violated: core m6A variant preserves DRACH at ",
       nrow(core_bad), " row(s), e.g. ",
       core_bad$snp_chr[1], ":", core_bad$snp_pos[1])
}

# Motif sequence recorded by 01_find_drach_sites.R
dt[, motif_seq := toupper(sub("^.*:", "", motif_name))]
dt[, motif_base_at_offset := substr(motif_seq, motif_offset + 1L, motif_offset + 1L)]
# BED methyl_pos is 0-based; report 1-based to match VCF POS convention
dt[, methyl_pos_1based := methyl_pos0 + 1L]

dt[, snp_id := paste(snp_chr, snp_pos, ref, alt, sep = ":")]
# motif_name is NOT unique: a peak sequence can contain the same DRACH expansion
# at several offsets (128,619 duplicate names in drach_sites.bed). A single
# variant can therefore hit offset 0 of one motif instance (preserving) and
# offset 4 of another (disrupting) simultaneously. The motif instance start is
# appended so the join key identifies one specific motif instance.
dt[, motif_id := paste0(motif_name, ":", motif_start0)]
dt[, overlap_key := paste(snp_id, motif_id, sep = "__")]

# One status per (variant, motif instance) by construction
n_conflict <- dt[, .(n_status = uniqueN(disruption_status)), by = overlap_key][
  n_status > 1L, .N]
if (n_conflict > 0L) {
  stop("Invariant violated: ", n_conflict,
       " overlap_key(s) carry more than one disruption status")
}

out_cols <- c("snp_id", "snp_chr", "snp_pos", "donor", "ref", "alt",
              "tx_ref", "tx_alt", "strand", "motif_offset", "is_core_m6a",
              "motif_seq", "motif_base_at_offset", "ref_allowed", "alt_allowed",
              "disruption_status", "motif_start0", "motif_end0",
              "methyl_pos_1based", "motif_id", "overlap_key")
setorder(dt, snp_chr, snp_pos, motif_id, donor)
fwrite(dt[, ..out_cols], out_table)

# --- Step 5: summary (deduplicated to one row per variant x motif instance) ---
per <- unique(dt, by = "overlap_key")
n_donors <- uniqueN(dt$donor)
n_disrupt <- per[disruption_status == "DISRUPTS_DRACH", .N]
n_preserve <- per[disruption_status == "PRESERVES_DRACH", .N]
n_core <- per[is_core_m6a == TRUE & disruption_status == "DISRUPTS_DRACH", .N]
n_disrupt_fwd <- per[disruption_status == "DISRUPTS_DRACH" & strand == "+", .N]
n_disrupt_rev <- per[disruption_status == "DISRUPTS_DRACH" & strand == "-", .N]

summary_dt <- data.table(
  metric = c(
    "donors_with_overlapping_variant", "donor_x_variant_x_motif_rows",
    "distinct_variant_sites", "distinct_motif_instances",
    "distinct_variant_x_motif_overlaps", "disrupts_drach_overlaps",
    "preserves_drach_overlaps", "core_m6a_disruption_overlaps",
    "disrupts_forward_strand", "disrupts_reverse_strand"),
  value = c(
    n_donors, nrow(dt), uniqueN(dt$snp_id), uniqueN(dt$motif_id),
    nrow(per), n_disrupt, n_preserve, n_core, n_disrupt_fwd, n_disrupt_rev),
  note = c(
    "distinct donors carrying >=1 variant overlapping a DRACH motif",
    "per-donor variant x motif rows (het and hom-alt included)",
    "distinct (chr:pos:ref:alt) sites overlapping a DRACH motif",
    "distinct DRACH motif instances intersected",
    "deduplicated (variant + motif) overlaps",
    "overlaps where the variant abolishes DRACH writer competence",
    "overlaps where the variant leaves DRACH intact",
    "overlaps ablating the central methylated adenine",
    "disrupting overlaps on the plus strand",
    "disrupting overlaps on the minus strand (require allele complementation)")
)
fwrite(summary_dt, out_summary)

message("  Annotated overlaps -> ", out_table)
message("  Cohort summary     -> ", out_summary)
message(strrep("=", 72))
message("   Distinct variant x motif overlaps : ", format(nrow(per), big.mark = ","))
message("   DISRUPTS_DRACH                    : ", format(n_disrupt, big.mark = ","))
message("   PRESERVES_DRACH                   : ", format(n_preserve, big.mark = ","))
message("   Core m6A ablations                : ", format(n_core, big.mark = ","))
message("   Disrupting  + / -                 : ", n_disrupt_fwd, " / ", n_disrupt_rev)
message(strrep("=", 72))

if (!file.exists(out_table) || file.info(out_table)$size == 0) {
  stop("Assertion Error: output table empty or not written: ", out_table)
}