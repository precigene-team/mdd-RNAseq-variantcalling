# [[01_core_alignment/08_tress_peak_calling.R]]
# Goal: Run TRESS peak calling on a single subject.

# Steps:
#   1. Load GenomicFeatures and build SQLite annotation from GTF if not present.
#   2. Load TRESS and verify package availability.
#   3. Parse subject name from command line argument.
#   4. Execute TRESS_peak for that subject.

if (!requireNamespace("TRESS", quietly = TRUE)) {
  stop("The 'TRESS' package is not installed. Please run the following in R within the active environment:\n",
       "  if (!require('BiocManager', quietly = TRUE)) install.packages('BiocManager')\n",
       "  BiocManager::install('TRESS')")
}

suppressPackageStartupMessages({
  library(TRESS)
  library(GenomicFeatures)
  library(AnnotationDbi)
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
bam_dir <- file.path(project_root, "_data/bam")
peaks_dir <- file.path(project_root, "_data/peaks")
gtf_file <- file.path(project_root, "_data/genome/Homo_sapiens.GRCh38.115.gtf")
sqlite_file <- file.path(project_root, "_data/genome/Homo_sapiens.GRCh38.115.sqlite")

# 1. Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("ERROR: Subject name must be provided as a command-line argument.")
}
subj <- args[1]

# 2. Build SQLite annotation database from GTF if missing
if (!file.exists(sqlite_file)) {
  message("[INFO] SQLite annotation database not found. Building from GTF...")
  txdb <- makeTxDbFromGFF(gtf_file, format = "gtf")
  saveDb(txdb, file = sqlite_file)
  message("[INFO] SQLite annotation database created successfully.")
}

# 3. Output directory for subject
subj_peaks_dir <- file.path(peaks_dir, subj)
dir.create(subj_peaks_dir, recursive = TRUE, showWarnings = FALSE)

# 4. Resolve BAM file paths
ip_bam <- file.path(bam_dir, paste0(subj, "_m6A"), paste0(subj, "_m6A_Aligned.sortedByCoord.out.bam"))
input_bam <- file.path(bam_dir, paste0(subj, "_input"), paste0(subj, "_input_Aligned.sortedByCoord.out.bam"))

if (!file.exists(ip_bam) || !file.exists(input_bam)) {
  stop("ERROR: BAM files for subject ", subj, " not found at paths:\n  ", ip_bam, "\n  ", input_bam)
}

# IP.file and Input.file paths relative to InputDir
ip_rel <- paste0(subj, "_m6A/", subj, "_m6A_Aligned.sortedByCoord.out.bam")
input_rel <- paste0(subj, "_input/", subj, "_input_Aligned.sortedByCoord.out.bam")

message(">>> Running TRESS_peak for subject: ", subj)
TRESS_peak(
  IP.file = ip_rel,
  Input.file = input_rel,
  Path_To_AnnoSqlite = sqlite_file,
  binsize = 50,
  WhichThreshold = "fdr_lfc",
  pval.cutoff0 = 1e-5,
  fdr.cutoff0 = 0.05,
  lfc.cutoff0 = log2(1.5),
  lowcount = 15,
  InputDir = bam_dir,
  OutputDir = subj_peaks_dir,
  experiment_name = subj,
  filetype = "bam",
  IncludeIntron = FALSE
)
message(">>> Completed TRESS_peak for subject: ", subj)
