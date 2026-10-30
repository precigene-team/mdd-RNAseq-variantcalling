# LEGACY — retired pipeline content

Scripts retired from `01_core_alignment/` and kept on disk under `_legacy/01_core_alignment/`.
`_legacy/` is gitignored (`.gitignore:69`), so these five files are **not part of the delivered
repository** — the deletion is the delivery, and `_legacy/` is the local backup.

## The retired MeRIP / TRESS arm

| Script | Role | Writes |
|---|---|---|
| `08_tress_peak_calling.sh` | Slurm wrapper, array `0-38%2`, 2 CPUs / 128G per task; skip guard on the subject's `peaks.xls` (`:66`); logs to `_data/logs/tress_peak/` | nothing directly — runs `Rscript 08_tress_peak_calling.R` at `:99` and `:102` |
| `08_tress_peak_calling.R` | `TRESS_peak()` on the IP fraction BAM only (`<subj>_m6A/<subj>_m6A_Aligned.sortedByCoord.out.bam`, `:62` and `:70`) | `_data/peaks/<subj>/` (`OutputDir`, `:58` and `:85`) |
| `09_merge_tress_peaks.R` | Merges and `reduce()`s every subject's peaks into one cohort peak set | `_data/peaks/consensus_peaks.bed` (`:108`) and `consensus_peaks.saf` (`:120`) — 33,123 rows |
| `10_tress_dmr_finder.sh` | Slurm submit per sex; skip guard `_data/dmr/<sex>/DMRs.RData` (`:66` and `:68`) | nothing directly |
| `10_tress_dmr_finder.R` | `TRESS_DMRfit` + `TRESS_DMRtest`, IP vs Input, MDD vs ctrl, paper-matched thresholds | `_data/dmr/<sex>/DMRs.RData` (`:146`), `results/01_core_alignment/tables/10_dmr_<sex>_full.csv` (`:150-151`), `10_dmr_<sex>_sig.csv` (`:157-158`) |

## Why retired

The MeRIP peak arm is not part of the delivered 3-module pipeline. `08_quantification.sh` counts
against `_data/genome/Homo_sapiens.GRCh38.115.gtf` (`:28`, `-F GTF` at `:92`), so the consensus
peak set is no longer a pipeline input, and `run`'s `core` section no longer invokes TRESS.

## Consumers of the retired artifacts, left untouched

These still read peak-derived files, and have no producer in the delivered tree:

- `02_dge_analysis/01_gene_mapping.R:6,44` — `_data/peaks/consensus_peaks.saf`
- `utils/ame_motif_enrichment.sh:25,31` — `consensus_peaks.bed` and `consensus_peaks.fa`
- `utils/annotate_consensus_peaks.R:18` — `consensus_peaks.bed`
- `utils/run_qc03_audit.sh:20` — `_data/peaks`
- `04_asm/01_find_drach_sites.R:26-27` — `consensus_peaks.bed` and `consensus_peaks.fa` (thesis arm)

The client-facing DGE tables in `results/02_dge_analysis/tables/` were built on the 33,123-peak
count matrix: `01_gene_mapper.csv` maps `peak_N` to the first protein-coding BioMart hit, and
`01_annotated_counts.csv` is keyed on those symbols. Gene-level counting is the replacement.
