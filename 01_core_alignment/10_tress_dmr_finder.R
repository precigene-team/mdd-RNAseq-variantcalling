#!/usr/bin/env Rscript
# [[01_core_alignment/10_tress_dmr_finder.R]]
# Goal: Run TRESS DMR analysis for a sex-stratified cohort (male or female)
#       to identify differentially methylated regions (DMRs) between MDD and ctrl.

# Steps:
#   1. Parse sex argument (male / female) from command line.
#   2. Load sample sheet and filter to the appropriate sex.
#   3. Build IP.file and Input.file vectors for the cohort.
#   4. Assert BAM files exist on the HPC path.
#   5. Build variable data.frame (SampleType: ctrl/mdd).
#   6. Run TRESS_DMRfit + TRESS_DMRtest with paper-matched thresholds.
#   7. Save full results table and significant DMR table.

if (!requireNamespace("TRESS", quietly = TRUE)) {
  stop("The 'TRESS' package is not installed.\n",
       "Run: BiocManager::install('TRESS')")
}

suppressPackageStartupMessages({
  library(TRESS)
  library(dplyr)
})

# --- 1. Parse argument ---
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("ERROR: Sex argument required. Usage: Rscript 10_tress_dmr_finder.R <male|female>")
}
sex <- tolower(args[1])
if (!sex %in% c("male", "female")) {
  stop("ERROR: sex must be 'male' or 'female', got: ", sex)
}
message("[INFO] Running TRESS DMR analysis for sex: ", sex)

# --- 2. Resolve project root ---
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

# --- 3. Paths ---
csv_file    <- file.path(project_root, "_data/_metadata/sample_sheet.csv")
bam_dir     <- file.path(project_root, "_data/bam")
sqlite_file <- file.path(project_root, "_data/genome/Homo_sapiens.GRCh38.115.sqlite")
out_dir     <- file.path(project_root, "_data/dmr", sex)
results_dir <- file.path(project_root, "results/tables")
dir.create(out_dir,     recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

stopifnot(file.exists(csv_file))
stopifnot(file.exists(sqlite_file))

# --- 4. Load sample sheet & filter to sex + IP fraction ---
metadata <- read.csv(csv_file, header = TRUE, stringsAsFactors = FALSE)
# Derive subject and fraction from File1 (e.g. male_mdd_1_m6A_1 -> male_mdd_1, m6A)
metadata$Subject  <- sub("_(m6A|input)_[12]$", "", metadata$File1)
metadata$Fraction <- ifelse(grepl("_m6A_", metadata$File1), "IP", "Input")

# Filter to this sex, one row per subject per fraction
cohort <- metadata %>%
  filter(Gender == sex) %>%
  distinct(Subject, SampleType, Fraction) %>%
  arrange(Fraction, Subject)

n_ip    <- sum(cohort$Fraction == "IP")
n_input <- sum(cohort$Fraction == "Input")
message("[INFO] Cohort size — IP: ", n_ip, ", Input: ", n_input)
stopifnot(n_ip == n_input)

# --- 5. Build IP.file and Input.file vectors ---
# TRESS expects paths relative to InputDir (bam_dir)
ip_subjects    <- cohort$Subject[cohort$Fraction == "IP"]
input_subjects <- cohort$Subject[cohort$Fraction == "Input"]

# Subjects must be in the same order for IP / Input pairing
stopifnot(identical(sort(ip_subjects), sort(input_subjects)))
ip_subjects    <- sort(ip_subjects)
input_subjects <- sort(input_subjects)

ip_files    <- paste0(ip_subjects,    "_m6A/",   ip_subjects,    "_m6A_Aligned.sortedByCoord.out.bam")
input_files <- paste0(input_subjects, "_input/", input_subjects, "_input_Aligned.sortedByCoord.out.bam")

# --- 6. Assert BAM files exist ---
message("[INFO] Checking BAM file existence...")
for (f in c(ip_files, input_files)) {
  full_path <- file.path(bam_dir, f)
  if (!file.exists(full_path)) {
    stop("ERROR: BAM file not found: ", full_path)
  }
}
message("[INFO] All ", length(ip_files) + length(input_files), " BAM files confirmed.")

# --- 7. Build variable data.frame ---
# One row per subject (same order as IP/Input file vectors)
sample_info <- cohort %>%
  filter(Fraction == "IP") %>%
  arrange(Subject) %>%
  select(Subject, SampleType)

variable_df <- data.frame(
  SampleType = factor(sample_info$SampleType, levels = c("ctrl", "mdd")),
  row.names  = sample_info$Subject
)
message("[INFO] Sample variable table:")
print(variable_df)

# --- 8. Run TRESS_DMRfit (model fitting) ---
# Paper-matched thresholds: p < 0.05, FC > 1.2
model <- ~ SampleType
message("[INFO] Running TRESS_DMRfit (this may take 30–90 min on HPC)...")
DMR.fit <- TRESS_DMRfit(
  IP.file            = ip_files,
  Input.file         = input_files,
  Path_To_AnnoSqlite = sqlite_file,
  InputDir           = bam_dir,
  OutputDir          = out_dir,
  experimentName     = "DMRfit",
  variable           = variable_df,
  model              = model,
  binsize            = 50,
  WhichThreshold     = "pval_lfc",
  pval.cutoff        = 0.05,
  fdr.cutoff         = 0.05,
  lfc.cutoff         = log2(1.2),
  filetype           = "bam",
  IncludeIntron      = FALSE,
  filterRegion       = FALSE
)
message("[INFO] TRESS_DMRfit complete.")

# --- 9. Run TRESS_DMRtest (hypothesis testing, mdd vs ctrl) ---
message("[INFO] Running TRESS_DMRtest...")
DMR.test <- TRESS_DMRtest(DMR = DMR.fit, contrast = c(0, 1))
message("[INFO] TRESS_DMRtest complete.")

# --- 10. Save results ---
# Save full fit object (for potential downstream loading)
save(DMR.fit, file = file.path(out_dir, "DMRs.RData"))

# Save full test results as CSV
dmr_df <- as.data.frame(DMR.test)
out_full <- file.path(results_dir, paste0("10_dmr_", sex, "_full.csv"))
write.csv(dmr_df, out_full, row.names = FALSE, quote = FALSE)
message("[INFO] Full DMR results saved: ", out_full)

# Save significant DMRs (p < 0.05, |logOR| >= log2(1.2))
sig_dmr_df <- dmr_df %>%
  filter(!is.na(pvalue), pvalue < 0.05, abs(logOR) >= log2(1.2))
out_sig <- file.path(results_dir, paste0("10_dmr_", sex, "_sig.csv"))
write.csv(sig_dmr_df, out_sig, row.names = FALSE, quote = FALSE)
message("[INFO] Significant DMRs (", nrow(sig_dmr_df), "): ", out_sig)

message("[DONE] 10_tress_dmr_finder.R complete for sex: ", sex)
