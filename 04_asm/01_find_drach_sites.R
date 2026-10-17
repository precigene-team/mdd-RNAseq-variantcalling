# [[04_asm/01_find_drach_sites.R]]
# Goal: Scan TRESS consensus peak sequences for all DRACH motif occurrences and convert each hit to genomic coordinates (BED).
#
# Steps:
#   1. Parse consensus_peaks.bed for peak coordinates.
#   2. Parse consensus_peaks.fa into a DNAStringSet.
#   3. Match the 18 concrete DRACH expansions on the forward strand.
#   4. Match them on the reverse complement and map offsets back to forward coordinates.
#   5. Write _data/peaks/drach_sites.bed and assert non-empty.
#
# Notes: matchPattern does not interpret IUPAC ambiguity codes, so DRACH is expanded
# to its 18 concrete sequences. Offset 0 is the leftmost base for '+' motifs and the
# rightmost for '-'; the methylated adenine sits at transcript offset 2 on both
# strands, so methyl_pos == genomic_start + 2. See docs/04_asm_methods.md.

suppressPackageStartupMessages(library(Biostrings))

# --- Dynamic Pathing ---
file_arg   <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg[1])))
} else getwd()
project_root <- Sys.getenv("PROJECT_ROOT",
                           unset = normalizePath(file.path(script_dir, "..")))

PEAKS_BED <- file.path(project_root, "_data", "peaks", "consensus_peaks.bed")
PEAKS_FA  <- file.path(project_root, "_data", "peaks", "consensus_peaks.fa")
OUT_BED   <- file.path(project_root, "_data", "peaks", "drach_sites.bed")

for (p in c(PEAKS_BED, PEAKS_FA)) {
  if (!file.exists(p)) stop("Assertion Error: missing input: ", p)
}

# --- DRACH definition: D=A/G/T, R=A/G, A, C, H=A/C/T  ->  3 * 2 * 3 = 18 ---
drach_patterns <- character(0)
for (d in c("A", "G", "T")) {
  for (r in c("A", "G")) {
    for (h in c("A", "C", "T")) {
      drach_patterns <- c(drach_patterns, paste0(d, r, "AC", h))
    }
  }
}
drach_patterns <- sort(drach_patterns)
stopifnot(length(drach_patterns) == 18L)
message("[INFO] DRACH consensus expanded to ", length(drach_patterns), " concrete sequences.")

MOTIF_LEN <- 5L

# --- Step 1: peak coordinates ---
message("[INFO] Loading consensus peaks BED...")
peaks <- read.delim(PEAKS_BED, header = FALSE, sep = "\t", quote = "",
                    col.names = c("chrom", "start", "end", "name", "score", "strand"))
message("[INFO] Loaded ", nrow(peaks), " peaks.")

# --- Step 2: peak sequences ---
read_fasta <- function(path) {
  lines <- readLines(path, warn = FALSE)
  hdr <- which(startsWith(lines, ">"))
  stopifnot(length(hdr) > 0L)
  body_len <- diff(c(hdr, length(lines) + 1L))
  nms  <- vapply(hdr, function(i) sub("^>", "", lines[i]), character(1))
  seqs <- vapply(seq_along(hdr),
                 function(k) toupper(paste0(lines[(hdr[k] + 1L):(hdr[k] + body_len[k] - 1L)],
                                           collapse = "")),
                 character(1))
  list(names = nms, seqs = seqs)
}

message("[INFO] Parsing consensus peak FASTA...")
fa <- read_fasta(PEAKS_FA)
peak_names <- sub("::.*$", "", fa$names)          # "peak_5603::5:81225659-81226459"
seqs_fwd   <- DNAStringSet(fa$seqs)
names(seqs_fwd) <- peak_names
message("[INFO] Loaded ", length(seqs_fwd), " peak sequences.")

# align BED rows to FASTA order
stopifnot(length(peak_names) == nrow(peaks))
idx_bed <- match(peak_names, peaks$name)
if (anyNA(idx_bed)) stop("Assertion Error: peak names differ between BED and FASTA")
peaks <- peaks[idx_bed, , drop = FALSE]
stopifnot(nchar(seqs_fwd) == peaks$end - peaks$start)

seqs_rev <- reverseComplement(seqs_fwd)

# --- Core matcher: returns data.frame(seq_index, offset, pattern) ---
# Biostrings::vmatchPattern returns a ByPos_MIndex; as.data.frame() expands it
# natively into one row per hit with `group` (1-based index into the subject set)
# and `start` (1-based). Using that avoids any manual index expansion.
scan_strand <- function(xset, strand_label) {
  out_seq <- integer(0); out_off <- integer(0); out_pat <- integer(0)
  for (pidx in seq_along(drach_patterns)) {
    m <- vmatchPattern(drach_patterns[pidx], xset)
    if (!length(m)) next
    df <- as.data.frame(m)
    if (nrow(df) == 0L) next
    n_hits <- nrow(df)
    out_seq <- c(out_seq, as.integer(df$group))
    out_off <- c(out_off, as.integer(df$start) - 1L)
    # NOTE: this R build cannot vectorize `times` in rep(), hence scalar times.
    out_pat <- c(out_pat, rep.int(pidx, n_hits))
    message("       ", strand_label, " ", drach_patterns[pidx], ": ",
            format(n_hits, big.mark = ","), " hits")
  }
  data.frame(seq_index = out_seq, offset = out_off, pat = out_pat)
}

# --- Step 3: forward strand ---
message("[STEP 3] Matching DRACH on the forward strand...")
fwd <- scan_strand(seqs_fwd, "fwd")

# --- Step 4: reverse strand ---
message("[STEP 4] Matching DRACH on the reverse complement...")
rev_hit <- scan_strand(seqs_rev, "rev")

# convert reverse-complement offsets back to forward coordinates
peak_len <- nchar(seqs_fwd)
rev <- data.frame(
  seq_index  = rev_hit$seq_index,
  offset     = peak_len[rev_hit$seq_index] - MOTIF_LEN - rev_hit$offset,
  pat        = rev_hit$pat,
  stringsAsFactors = FALSE
)

# --- Step 5: assemble records ---
make_records <- function(hits, strand) {
  if (nrow(hits) == 0L) {
    return(data.frame(chrom = character(0), gs = integer(0), ge = integer(0),
                      name = character(0), score = character(0),
                      strand = character(0), methyl = integer(0),
                      stringsAsFactors = FALSE))
  }
  i    <- hits$seq_index
  gs   <- as.integer(peaks$start[i]) + hits$offset
  data.frame(
    chrom  = peaks$chrom[i],
    gs     = gs,
    ge     = gs + MOTIF_LEN,
    name   = paste0(peak_names[i], ":DRACH_", strand, ":", drach_patterns[hits$pat]),
    score  = "0",
    strand = if (strand == "fwd") "+" else "-",
    methyl = gs + 2L,                 # transcript offset 2 on both strands
    stringsAsFactors = FALSE
  )
}

recs <- rbind(make_records(fwd, "fwd"), make_records(rev, "rev"))
message("[INFO] Found ", format(nrow(recs), big.mark = ","), " DRACH sites across all peaks.")

# sort by (chrom, genomic_start) using byte order to match the previous output
recs <- recs[order(recs$chrom, recs$gs, method = "radix"), , drop = FALSE]

# --- Step 6: write BED ---
message("[INFO] Writing DRACH sites to: ", OUT_BED)
con <- file(OUT_BED, "wt")
writeLines("# chrom\tstart\tend\tname\tscore\tstrand\tmethyl_genomic_pos", con)
writeLines(paste(recs$chrom, recs$gs, recs$ge, recs$name,
                 recs$score, recs$strand, recs$methyl, sep = "\t"), con)
close(con)

if (file.info(OUT_BED)$size == 0) stop("Assertion Error: output DRACH BED is empty.")
message("[DONE] ", format(nrow(recs), big.mark = ","),
        " DRACH sites written to ", OUT_BED)
message("       Next: run 02_verify_drach_reference.py, then 03_intersect_drach_genotypes.sh")