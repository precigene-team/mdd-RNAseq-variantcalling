# LEGACY — retired pipeline content

Retired from `01_core_alignment/`: scripts kept on disk under `_legacy/01_core_alignment/`.
`_legacy/` is gitignored (`.gitignore:69`), so those five files are **not part of the delivered
repository** — the deletion is the delivery, and `_legacy/` is the local backup.

Retired from `02_dge_analysis/`: MeRIP blocks removed inside the scripts, **not archived** —
the previous versions live in the git history, so no `_legacy/` copy exists. Whole scripts
(`10`–`13`) are deleted from the tracked tree and archived under `_legacy/02_dge_analysis/`.

Every change that retires pipeline content is recorded here.

## The retired MeRIP / TRESS arm

| Script | Role | Writes |
|---|---|---|
| `08_tress_peak_calling.sh` | Slurm wrapper, array `0-38%2`, 2 CPUs / 128G per task; skip guard on the subject's `peaks.xls` (`:66`); logs to `_data/logs/tress_peak/` | nothing directly — runs `Rscript 08_tress_peak_calling.R` at `:99` and `:102` |
| `08_tress_peak_calling.R` | `TRESS_peak()` on the IP fraction BAM only (`<subj>_m6A/<subj>_m6A_Aligned.sortedByCoord.out.bam`, `:62` and `:70`) | `_data/peaks/<subj>/` (`OutputDir`, `:58` and `:85`) |
| `09_merge_tress_peaks.R` | Merges and `reduce()`s every subject's peaks into one cohort peak set | `_data/peaks/consensus_peaks.bed` (`:108`) and `consensus_peaks.saf` (`:120`) — 33,123 rows |
| `10_tress_dmr_finder.sh` | Slurm submit per sex; skip guard `_data/dmr/<sex>/DMRs.RData` (`:66` and `:68`) | nothing directly |
| `10_tress_dmr_finder.R` | `TRESS_DMRfit` + `TRESS_DMRtest`, IP vs Input, MDD vs ctrl, paper-matched thresholds | `_data/dmr/<sex>/DMRs.RData` (`:146`), `results/01_core_alignment/tables/10_dmr_<sex>_full.csv` (`:150-151`), `10_dmr_<sex>_sig.csv` (`:157-158`) |

### Why retired

The MeRIP peak arm is not part of the delivered 3-module pipeline. `08_quantification.sh` counts
against `_data/genome/Homo_sapiens.GRCh38.115.gtf` (`:28`, `-F GTF` at `:92`), so the consensus
peak set is no longer a pipeline input, and `run`'s `core` section no longer invokes TRESS.

## The retired MeRIP arm in `02_dge_analysis/`

| File | Removed | Role of the removed code | Outputs that no longer exist |
|---|---|---|---|
| `01_gene_mapping.R` | peak-era mapper: SAF load (`:44,53,61-72`), BioMart `chromosomal_region` query (`:103-141`), coordinate-overlap peak→gene with `slice(1)` one gene per peak (`:143-176`), `biomaRt` | mapped `peak_N` rows of the peak-level count matrix to HGNC symbols by genomic overlap (215 → 144 lines) | `01_annotated_counts.csv` and `01_gene_mapper.csv` are now gene-level: ENSG → HGNC from the Ensembl GTF, 39 Input samples, protein-coding only |
| `02_rna_processing.R` | `DataType == "m6A" ~ "RIPseq"` (`:73`), header lines `:2,6` | labelled the IP libraries as a second datatype so DESeq2 ran on them | `02_RIPseq_filtered_counts.csv`, `02_RIPseq_vst_normalized_counts.csv`, `02_RIPseq_dds.rds`, `02_RIPseq_vst_normalized.rds` |
| `03_pca_analysis.R` | `"RIPseq"` from `datatypes` (`:254`) | ran PCA on the IP datatype; hit `stop()` at `:195`, swallowed by `tryCatch` as a log line | `03_RIPseq_PCA_*.csv`, `03_RIPseq_pca*.rds`, `03_RIPseq_PCA_3D_*.html` |
| `04_deg_analysis.R` | Phase B (`:277-304`), Phase C (`:306-385`), `dds_ripseq` load + the `02_RIPseq_dds.rds` assertion (`:100-108`), IP half of the 2M library-size gate (`:115,119,129,133`), `Fraction` metadata + relevel (`:58,62,77,90-93`), header steps 5-6 (`:13-19`) and the `Phase A` label (`:225,227`) | MeRIP DE and relative-methylation contrasts on combined Input+IP counts with nested-subject designs (387 → 256 lines) | `04_deg_male_ip*`, `04_deg_female_ip*`, `04_meth_male_interaction*`, `04_meth_female_interaction*`, `04_meth_sex_interaction*` tables and RDS objects |
| `05_dge_plotting.R` | `deg_male_ip` (`:174-178`), `deg_female_ip` (`:180-184`), `meth_male_interaction` (`:186-190`), `meth_female_interaction` (`:192-196`), `meth_sex_interaction` (`:204-208`, title `Sex-Specific m⁶A Methylation\n(Gender × Condition × Fraction)`), header `:2` | volcano and MA plots for the five MeRIP contrasts; the missing-RDS gate at `:39-42` logged `[SKIP]` for each (226 → 196 lines, 8 → 3 contrasts) | `05_volcano_male_ip.png`, `05_ma_male_ip.png`, `05_volcano_female_ip.png`, `05_ma_female_ip.png`, `05_volcano_male_meth.png`, `05_ma_male_meth.png`, `05_volcano_female_meth.png`, `05_ma_female_meth.png`, `05_volcano_sex_interaction_meth.png`, `05_ma_sex_interaction_meth.png` |
| `07_enrichment_analysis.R` | `deg_male_ip`, `deg_female_ip`, `meth_male_interaction`, `meth_female_interaction`, `meth_sex_interaction` from `contrasts` (`:56-60`), stale `peak_id → gene_name` mapper comment (`:42`) | clusterProfiler GO/KEGG ORA and GSEA GO/Hallmark on the five MeRIP contrasts; the missing-RDS gate at `:67-70` logged `[Skip]` for each (195 → 190 lines, 8 → 3 contrasts) | `07_deg_male_ip_*.csv`, `07_deg_female_ip_*.csv`, `07_meth_male_interaction_*.csv`, `07_meth_female_interaction_*.csv`, `07_meth_sex_interaction_*.csv` |
| `09_heatmap_pathway.R` | Input/IP independent row-scaling and the Delta (IP − Input) matrix (`:165-218`), `had_fraction_originally` and `is_paired_delta_mode` (`:166-167`), the two `had_fraction_originally` annotation branches (`:250-280`, `:328-351`), the duplicated `is_paired_delta_mode` branches, the `m6A → IP` sample-name regex (`:158`), the 5 MeRIP titles in `get_descriptive_title` (`:390-394`), the Fraction legend (`:438-443`), the `Δ z-score` name expressions (`:403,454`) | per-fraction scaling, paired Subject delta heatmaps, Fraction column annotation and legend, MeRIP heatmap titles (474 → 308 lines; 48 insertions, 216 deletions) | `09_heatmap_pathway_deg_male_ip.png`, `09_heatmap_pathway_deg_female_ip.png`, `09_heatmap_pathway_meth_male_interaction.png`, `09_heatmap_pathway_meth_female_interaction.png`, `09_heatmap_pathway_meth_sex_interaction.png` |

`Subject` (`04_deg_analysis.R:57`) is kept: the library-size gate excludes by subject.
The three Input contrasts (`deg_male_input`, `deg_female_input`, `deg_sex_interaction_input`) are
unchanged and make up the whole delivered DGE analysis and its plots.

## The retired scripts `10`–`13` in `02_dge_analysis/`

Deleted from the tracked tree, archived in `_legacy/02_dge_analysis/`. The delivered DGE module
ends at `09_heatmap_pathway.R`.

| Script | Role | Outputs that no longer have a producer |
|---|---|---|
| `10_rrho2_sex_comparison.R` (147 lines) | RRHO2 cross-sex overlap on the male and female Input DEG tables (`tables_dir/04_deg_male_input.csv`, `04_deg_female_input.csv`, `:49-50`): rank score `-log10(padj) * sign(log2FC)` (`:74-84`), shared-gene universe assertion (`:89-102`), `RRHO2_initialize` + `RRHO2_heatmap` (`:112-125`). Invoked by `run:317` (line deleted) | `results/02_dge_analysis/figures/10_rrho2_concordant.pdf` (`:123`), `results/02_dge_analysis/tables/10_rrho_hypermat.csv` (`:130`) |
| `11_gsea_input.R` (130 lines) | Preranked GSEA input files for `deg_male_input` and `deg_female_input` only (`:47`); two ranking metrics per gene — DESeq2 Wald `stat` and `signed_p = sign(log2FC) * -log10(pvalue)` (`:81-82`); `.cls` built from `dds$SampleType` (`:103`) | `results/02_dge_analysis/gsea/<contrast>.rnk`, `<contrast>_signedp.rnk`, `<contrast>_ranked.csv`, `<contrast>.cls`, `hallmark.gmt`, `gsea_manifest.csv` (`:92-128`) |
| `12_gsea_fgsea.R` (80 lines) | `fgseaMultilevel` (`minSize = 15`, `maxSize = 500`, `eps = 0`, `scoreType = "std"`, `nproc = 1`) on the `.rnk` files from 11 against `hallmark.gmt`; the only `fgsea` consumer in the tree, and `fgsea` is absent from `environment.yml` | `tables/12_deg_male_input_gsea_hallmark.csv`, `12_deg_female_input_gsea_hallmark.csv` (`:66-68`) |
| `13_gsea_figures.R` (270 lines) | Publication figures over 12's tables: NES barplot per contrast, merged male/female dotplot with dumbbell segments joining each set's two NES values, `fgsea` enrichment plots for the top 3 male sets, male-vs-female NES scatter over all 50 sets, NES distribution | `figures/13_<contrast>_gsea_hallmark_nes.png` (`:99`), `13_gsea_hallmark_dotplot.png` (`:154`), `13_deg_male_input_gsea_hallmark_enrichment.png` (`:192`), `13_gsea_hallmark_nes_scatter.png` (`:228`), `13_gsea_hallmark_nes_distribution.png` (`:266`) |

`10`'s header (`:2,5,6,7,10`) and the `stop()` message (`:55`) were peak-era language ("DMR ranked
lists", "from TRESS differential analysis", "for each peak", "shared peak universe", "Run TRESS
differential analysis first"), and the variable names were `male_dmr_file`/`female_dmr_file`,
`male_dmr`/`female_dmr`, `shared_peaks`. The computation itself was gene-level correct, so the
whole file was retired rather than retargeted.

`11`–`13` duplicate work the delivered module already does: `07_enrichment_analysis.R` runs GSEA GO
Biological Process and MSigDB Hallmark with clusterProfiler, and `08_enrichment_plotting.R` renders
the NES barplots. `11`–`13` used `fgsea` on `.rnk` files with no `run` invocation
(`run` invoked only `10`, at `:317`).

## Consumers of the retired artifacts, left untouched

Peak-derived files, with no producer in the delivered tree:

- `utils/ame_motif_enrichment.sh:25,31` — `consensus_peaks.bed` and `consensus_peaks.fa`
- `utils/annotate_consensus_peaks.R:18` — `consensus_peaks.bed`
- `utils/run_qc03_audit.sh:20` — `_data/peaks`
- `04_asm/01_find_drach_sites.R:26-27` — `consensus_peaks.bed` and `consensus_peaks.fa` (thesis arm)

Retired-script consumers, with no producer in the delivered tree:

- `environment.yml:85` — `bioconductor-rrho`, the RRHO2 dependency `10_rrho2_sex_comparison.R` bootstrapped
- `docs/pipeline_problem_list.md:61-64` — Item 6 prescribes RRHO2 on "male and female DMR ranked lists"; `:115` lists it as pipeline step 10
- `README.md:61-63,136-139` — name three `03_variant_calling/` scripts that do not exist (`07_clinical_pathogenicity.py`, `08_pharmacogenomics_actionability.py`, `09_plot_clinical_annotation_landscape.R`); `README.md:53-76` omits 5 real ones (`06_build_transcriptome_snp_table.py`, `07_functional_annotation.py`, `08_map_gtex_brain_eqtls.py`, `10_prioritize_mdd_candidate_snps.py`, `14_slice_igv_bams.sh`); `README.md:86` claims retired ASM scripts are archived in `_legacy/04_asm/`, which does not exist; `README.md:88-90` lists 2 of the 25 files in `docs/`; `README.md:172` omits R libraries the delivered 04–09 use (`ashr`, `clusterProfiler`, `msigdbr`, `org.Hs.eg.db`, `ComplexHeatmap`, `circlize`, `RColorBrewer`, `tidyverse`, `ggtext`) — all reported, not fixed

`08_enrichment_plotting.R` has no MeRIP content: it discovers work by globbing `07_*_gsea_*.csv`
(`:47`) and parses the contrast name from the filename (`:58`), so the MeRIP plot names disappear
when 07 stops producing those tables.

The client-facing DGE tables in `results/02_dge_analysis/tables/` were built on the 33,123-peak
count matrix: `01_gene_mapper.csv` mapped `peak_N` to the first protein-coding BioMart hit, and
`01_annotated_counts.csv` was keyed on those symbols. Gene-level counting (78,899 features per
sample) is the replacement, so every 02 table must be regenerated.

## `02_dge_analysis/config.yaml` aligned to the delivered pipeline

The file is kept because `01_gene_mapping.R:33-35` requires it (`stop()` if absent) and reads
`cfg$ensembl_map$valid_biotypes` (`:38`), and `06_annotate_results.R:81` reads
`cfg$ensembl_map$ensembl_version`. Everything else was foreign to the delivered scope.
89 lines → 26 lines (16 insertions, 79 deletions).

Removed:
- header `:1-13` — the 10-step order naming `11_quantification.sh`, `12_tress_dmr_finder`, `15_pathway_analysis`, the MeRIP arm, and the TCGA cohort
- `cohorts:` `:19-22` — `mdd_input` / `mdd_m6a` / `tumor` / `normal` sample lists, consumed by no script
- `ensembl_map:input_file` and `output_dir` `:29-34` — the peak→gene mapper's paths, rewritten to
  `ensembl_version` and `valid_biotypes` (key name preserved; no script reads the input/output paths)
- `data_processing:` `:44-68` — `min_count`, `min_frac`, `mirbase_version`, `outliers` (`pca_threshold`
  and `filter_by_expr` are used, so they are kept); `mirbase_version` is the miRBase mapper's version
- `de:` `:74-89` — `logfc_cutoff`, `q_cutoff`, `het_p_cutoff`, `het_i2_cutoff`; `04_deg_analysis.R` hardcodes
  `p_adjust < 0.05` and `abs(log2FoldChange) > 1`, so the DE thresholds were not read by anything.
  `lfc_threshold` is kept: `07_enrichment_analysis.R:29` reads it (`lfc_cutoff` is kept too)

## Change log

| Date | Files | Change |
|---|---|---|
| 2026-10-30 | `01_core_alignment/08_tress_peak_calling.{R,sh}`, `09_merge_tress_peaks.R`, `10_tress_dmr_finder.{R,sh}` | MeRIP/TRESS arm deleted, archived in `_legacy/01_core_alignment/` |
| 2026-10-30 | `11_quantification.sh → 08_quantification.sh`, `12 → 09_counts_qc.sh`, `13 → 10_multiqc.sh` | gene-level featureCounts against the Ensembl GTF, numbering gap closed |
| 2026-10-30 | `README.md`, `CONTEXT.md`, `run`, `docs/pipeline_problem_list.md`, `docs/report.md` | scope reduced to the 39 Input samples stopping at variant calling; MeRIP objectives, contrast B/C, and the 04_asm CONTEXT item removed |
| 2026-10-30 | `02_dge_analysis/01_gene_mapping.R` | peak→gene coordinate mapper replaced by ENSG → HGNC from `Homo_sapiens.GRCh38.115.gtf`; matrix restricted to the 39 Input samples; protein-coding filter from `config.yaml:37-38`; 70 insertions, 141 deletions |
| 2026-10-30 | `02_dge_analysis/02_rna_processing.R` | `m6A → RIPseq` DataType mapping deleted; header retargeted to Input; 3 insertions, 4 deletions |
| 2026-10-30 | `02_dge_analysis/03_pca_analysis.R` | `"RIPseq"` deleted from `datatypes` (`:254`); 1 insertion, 1 deletion |
| 2026-10-30 | `02_dge_analysis/04_deg_analysis.R` | Phase B and Phase C deleted, IP dds load and IP QC-gate half deleted, `Fraction` removed, header retargeted to Input, `Phase A` label dropped, `A1./A2./A3.` renumbered `1./2./3.`; 27 insertions, 158 deletions (387 → 256 lines) |
| 2026-10-30 | `02_dge_analysis/05_dge_plotting.R` | 5 MeRIP contrasts (`deg_male_ip`, `deg_female_ip`, `meth_male_interaction`, `meth_female_interaction`, `meth_sex_interaction`) deleted from `contrasts_list`; header retargeted to Input; 1 insertion, 31 deletions (226 → 196 lines) |
| 2026-10-30 | `02_dge_analysis/06_annotate_results.R` | SRA accession export restricted to `DataType == "RNA-seq"` and `DataType` dropped from `group_by` (8 subgroups → 4 Input subgroups); header retargeted to Input; BioMart annotation source kept; 8 insertions, 7 deletions |
| 2026-10-30 | `02_dge_analysis/07_enrichment_analysis.R` | 5 MeRIP contrasts deleted from `contrasts`; `peak_id → gene_name` mapper comment corrected to `ensg → gene_name, gene_type` (matches what `01_gene_mapping.R` now writes); 2 insertions, 7 deletions (195 → 190 lines) |
| 2026-10-30 | `02_dge_analysis/09_heatmap_pathway.R` | per-fraction scaling + paired-Subject Delta matrix, `had_fraction_originally` / `is_paired_delta_mode`, the duplicated `is_paired_delta_mode` branches, the `m6A → IP` name regex, 5 MeRIP titles, the Fraction legend and the `Δ z-score` names deleted; standard row z-scoring is now the only path; 48 insertions, 216 deletions (474 → 308 lines) |
| 2026-10-31 | `10_rrho2_sex_comparison.R`, `11_gsea_input.R`, `12_gsea_fgsea.R`, `13_gsea_figures.R` (deleted), `run:317` (deleted) | retired the 4 scripts after `09_heatmap_pathway.R`; `run dge` ends at `09`; copies archived in `_legacy/02_dge_analysis/`; no reference to any of them remains outside `_legacy/` |
| 2026-10-31 | `02_dge_analysis/config.yaml` | 89 lines → 26 (16 insertions, 79 deletions): `ensembl_map:` kept as `ensembl_version` + `valid_biotypes` for `01:38` and `06:81`; header, `cohorts:`, `mirbase_version`, `outliers`, `min_count`/`min_frac`, `de:` cutoffs, and the `11`-`15` step numbering removed |
| 2026-10-31 | `README.md` | `02_dge_analysis` section aligned to the delivered 9 scripts (`:43-53` listing, section 2 command block extended to `09`), `:45` reworded to ENSG -> HGNC from the Ensembl GTF, Input samples only |
