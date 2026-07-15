# [[01_core_alignment/09_merge_tress_peaks.R]]
# Goal: Merge TRESS peak coordinates across all samples into cohort consensus BED and SAF files.

# Steps:
#   1. Read subjects from sample sheet.
#   2. Load _peaks.xls files for all subjects and parse coordinates into GRanges.
#   3. Concatenate and merge (reduce) overlapping peaks across all subjects.
#   4. Export to BED (0-indexed start) and SAF (1-indexed start) formats.

suppressPackageStartupMessages({
  library(GenomicRanges)
})

# Dynamically resolve project root
project_root <- Sys.getenv("PROJECT_ROOT")
if (project_root == "") {
  get_project_root <- function() {
    curr <- getwd()
    while (curr != "/" && !file.exists(file.path(curr, "run"))) {
      curr <- dirname(curr)
    }
    return(curr)
  }
  project_root <- get_project_root()
}

# Paths
csv_file <- file.path(project_root, "_data/_metadata/sample_sheet.csv")
peaks_dir <- file.path(project_root, "_data/peaks")

metadata <- read.csv(csv_file, header = TRUE, stringsAsFactors = FALSE)
metadata$Sample <- sub("_1$", "", metadata$File1)
metadata$Subject <- sub("_(m6A|input)$", "", metadata$Sample)
subjects <- unique(metadata$Subject)

all_peaks <- list()
loaded_count <- 0

message("[INFO] Reading TRESS peaks for all subjects...")
for (subj in subjects) {
  xls_file <- file.path(peaks_dir, subj, paste0(subj, "_peaks.xls"))
  if (!file.exists(xls_file)) {
    message("  [WARNING] Peak file not found for subject: ", subj)
    next
  }
  
  # Read the tab-delimited xls file
  peaks_df <- tryCatch({
    read.table(xls_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  }, error = function(e) {
    message("  [ERROR] Failed to read peak file for subject ", subj, ": ", e$message)
    NULL
  })
  
  if (is.null(peaks_df) || nrow(peaks_df) == 0) {
    next
  }
  
  # Identify chromosome, start, end columns (TRESS columns are chr, start, end)
  chr_col <- if ("chr" %in% colnames(peaks_df)) "chr" else if ("seqnames" %in% colnames(peaks_df)) "seqnames" else NULL
  start_col <- if ("start" %in% colnames(peaks_df)) "start" else NULL
  end_col <- if ("end" %in% colnames(peaks_df)) "end" else NULL
  strand_col <- if ("strand" %in% colnames(peaks_df)) "strand" else NULL
  
  if (is.null(chr_col) || is.null(start_col) || is.null(end_col)) {
    message("  [WARNING] Columns chr/start/end not found in peak file for subject: ", subj)
    next
  }
  
  # Create GRanges
  strand_vals <- if (!is.null(strand_col)) peaks_df[[strand_col]] else "*"
  gr <- GRanges(
    seqnames = peaks_df[[chr_col]],
    ranges = IRanges(start = peaks_df[[start_col]], end = peaks_df[[end_col]]),
    strand = strand_vals
  )
  
  loaded_count <- loaded_count + 1
  all_peaks[[loaded_count]] <- gr
}

if (loaded_count == 0) {
  stop("No peaks could be loaded across subjects.")
}

message("[INFO] Loaded peaks from ", loaded_count, " subjects. Concatenating and merging...")

# Concatenate all GRanges
combined_gr <- do.call(c, all_peaks)

# Reduce to merge overlapping peaks
consensus_gr <- reduce(combined_gr)

message("[INFO] Merging complete. Consensus peaks: ", length(consensus_gr))

# Create output folder if not exist
dir.create(peaks_dir, recursive = TRUE, showWarnings = FALSE)

# 1. Export as BED format (0-indexed start, 1-indexed end)
bed_df <- data.frame(
  chr = as.character(seqnames(consensus_gr)),
  start = start(consensus_gr) - 1,
  end = end(consensus_gr),
  name = paste0("peak_", seq_along(consensus_gr)),
  score = 0,
  strand = as.character(strand(consensus_gr))
)
write.table(bed_df, file = file.path(peaks_dir, "consensus_peaks.bed"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
message("Saved consensus peaks to BED: ", file.path(peaks_dir, "consensus_peaks.bed"))

# 2. Export as SAF format (GeneID, Chr, Start, End, Strand) - NO HEADER
saf_df <- data.frame(
  GeneID = paste0("peak_", seq_along(consensus_gr)),
  Chr = as.character(seqnames(consensus_gr)),
  Start = start(consensus_gr),
  End = end(consensus_gr),
  Strand = as.character(strand(consensus_gr))
)
write.table(saf_df, file = file.path(peaks_dir, "consensus_peaks.saf"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
message("Saved consensus peaks to SAF: ", file.path(peaks_dir, "consensus_peaks.saf"))
