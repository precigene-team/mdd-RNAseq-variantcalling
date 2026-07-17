#!/usr/bin/env Rscript
# utils/annotate_consensus_peaks.R
# Goal: Annotate consensus_peaks.bed against GRCh38 GTF to count unique
#       transcripts and genes, replicating the paper's "30,279 peaks in
#       12,691 transcripts, avg 2.4 peaks/gene" metric.
#       Applies protein-coding biotype filter to match paper's gene universe.

suppressPackageStartupMessages({
  library(GenomicFeatures)
  library(GenomicRanges)
  library(AnnotationDbi)
  library(rtracklayer)
  library(dplyr)
})

# --- Paths ---
project_root <- getwd()
bed_file     <- file.path(project_root, "_data/peaks/consensus_peaks.bed")
gtf_file     <- file.path(project_root, "_data/genome/Homo_sapiens.GRCh38.115.gtf")
sqlite_file  <- file.path(project_root, "_data/genome/Homo_sapiens.GRCh38.115.sqlite")
out_file     <- file.path(project_root, "results/utils/peak_annotation_summary.txt")

stopifnot(file.exists(bed_file))
stopifnot(file.exists(gtf_file))

# --- 1. Load or build TxDb from GTF ---
if (file.exists(sqlite_file)) {
  message("[INFO] Loading TxDb from existing SQLite...")
  txdb <- loadDb(sqlite_file)
} else {
  message("[INFO] Building TxDb from GTF (this may take ~5 min)...")
  txdb <- makeTxDbFromGFF(gtf_file, format = "gtf")
  saveDb(txdb, file = sqlite_file)
  message("[INFO] Saved TxDb to: ", sqlite_file)
}

# --- 2. Load consensus peaks as GRanges ---
message("[INFO] Loading consensus peaks...")
bed <- read.table(bed_file, header = FALSE, sep = "\t",
                  col.names = c("chr", "start", "end", "name", "score", "strand"),
                  stringsAsFactors = FALSE)

peaks_gr <- GRanges(
  seqnames = bed$chr,
  ranges   = IRanges(start = bed$start + 1, end = bed$end),  # BED is 0-indexed
  strand   = bed$strand,
  peak_id  = bed$name
)
message("[INFO] Total consensus peaks: ", length(peaks_gr))

# --- 3. Get protein-coding gene IDs from GTF ---
message("[INFO] Reading gene biotypes from GTF (filtering to protein_coding)...")
gtf_gr <- import(gtf_file, format = "gtf", feature.type = "gene",
                 colnames = c("gene_id", "gene_biotype"))
protein_coding_genes <- unique(gtf_gr$gene_id[gtf_gr$gene_biotype == "protein_coding"])
message("[INFO] Total protein-coding genes in annotation: ", length(protein_coding_genes))

# --- 4. Get transcript GRanges from TxDb, filtered to protein-coding ---
message("[INFO] Extracting transcript annotations (protein-coding only)...")
tx_gr_all <- transcripts(txdb, columns = c("tx_id", "tx_name", "gene_id"))

# gene_id is stored as a CharacterList in TxDb — unlist to filter
gene_ids_unlisted <- unlist(mcols(tx_gr_all)$gene_id)
tx_is_pc <- vapply(mcols(tx_gr_all)$gene_id,
                   function(g) any(g %in% protein_coding_genes),
                   logical(1))
tx_gr <- tx_gr_all[tx_is_pc]
message("[INFO] Protein-coding transcripts retained: ", length(tx_gr), " / ", length(tx_gr_all))

# --- 5. Overlap peaks with protein-coding transcripts ---
message("[INFO] Overlapping peaks with transcripts...")
hits <- findOverlaps(peaks_gr, tx_gr, ignore.strand = FALSE)

# Peaks with at least one transcript overlap
peaks_in_tx <- length(unique(queryHits(hits)))
message("[INFO] Peaks overlapping at least one transcript: ", peaks_in_tx)

# --- 6. Count unique transcripts hit ---
tx_hits_df <- data.frame(
  peak_id   = peaks_gr$peak_id[queryHits(hits)],
  tx_name   = mcols(tx_gr)$tx_name[subjectHits(hits)],
  gene_id   = unlist(mcols(tx_gr)$gene_id[subjectHits(hits)])
)

n_unique_transcripts <- length(unique(tx_hits_df$tx_name))
n_unique_genes       <- length(unique(tx_hits_df$gene_id))
message("[INFO] Unique transcripts with peaks: ", n_unique_transcripts)
message("[INFO] Unique genes with peaks: ", n_unique_genes)

# --- 7. Compute peaks per gene ---
peaks_per_gene <- tx_hits_df %>%
  group_by(gene_id) %>%
  summarise(n_peaks = n_distinct(peak_id)) %>%
  ungroup()

avg_peaks_per_gene <- mean(peaks_per_gene$n_peaks)
message("[INFO] Avg peaks per gene: ", round(avg_peaks_per_gene, 2))

# --- 8. Write summary ---
summary_lines <- c(
  "=== Consensus Peak Annotation Summary (protein-coding only) ===",
  paste0("Total consensus peaks (BED):           ", length(peaks_gr)),
  paste0("Peaks overlapping PC transcripts:      ", peaks_in_tx),
  paste0("Unique PC transcripts with peaks:      ", n_unique_transcripts),
  paste0("Unique protein-coding genes with peaks:", n_unique_genes),
  paste0("Avg peaks per gene:                    ", round(avg_peaks_per_gene, 2)),
  "",
  "=== Paper benchmark (Mitsuhashi et al.) ===",
  "Total consensus peaks:  30,279",
  "Unique transcripts:     12,691",
  "Avg peaks per gene:     2.4"
)
writeLines(summary_lines, out_file)
cat(paste(summary_lines, collapse = "\n"), "\n")
message("[DONE] Summary written to: ", out_file)
