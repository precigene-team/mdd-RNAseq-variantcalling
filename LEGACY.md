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

## The retired thesis arm `04_asm/`

Deleted from the tracked tree (`git rm -r 04_asm/`), archived in `_legacy/04_asm/` (sequence-level
m6A-SNP arm, 7 scripts, 1,676 lines) and `_legacy/04_asm/scripts/` (the read-ratio ASM generation,
14 files, 1,950 lines). The delivered pipeline is 3 modules, and the thesis arm was the
MeRIP-dependent fourth: every one of the 7 scripts takes TRESS peak artifacts as input.

| Script | Role |
|---|---|
| `01_find_drach_sites.R` (162 lines) | DRACH motif occurrences in `consensus_peaks.bed` + `consensus_peaks.fa` → BED |
| `02_verify_drach_reference.py` (287) | Hard QC gate on motif geometry read from peak sequence; exits non-zero on corruption |
| `03_intersect_drach_genotypes.sh` (102) | PASS biallelic genotypes × DRACH sites; delegates 02 then 04 |
| `04_annotate_strand_aware_disruption.R` (164) | Variant × motif instance → DISRUPTS/PRESERVES writer competence |
| `05_plot_m6asnp_funnel.py` (295) | Six-level discovery funnel + core-m6A site table consumed by 06 and 07 |
| `06_test_mdd_enrichment.py` (291) | MDD-vs-Ctrl enrichment, donor labels permuted within sex |
| `07_allele_specific_m6a.py` (375) | Per-donor IP-vs-Input ASm6A, concordant sites at >90% cross-donor agreement |

History, as the git log states it: `23ceab9` renamed 14 MeRIP/ASM files out of
`03_variant_calling/` into `04_asm/` (`R100`), so `run`'s `03_variant_calling/07_m6allele_asm.sh`
path was already stale before this pass. `eb75c6e` — "replace read-ratio ASM with sequence-level
m6A-SNP pipeline" — deleted that set and added the 7. `_legacy/04_asm/scripts/` is the `eb75c6e^`
state, and the path the docs already cite (`docs/asm_filter_flow.md:5-6`,
`docs/asm_analysis_methodology.md:15`, `docs/04_asm_methods.md:24`).

## Input-only `_data/_metadata/sample_sheet.csv`

`sample_sheet.csv` is now the header + 39 `DataType == "RNA-seq"` rows, in the original order:
10 male mdd, 10 male ctrl, 10 female mdd, 9 female ctrl. The 78-row cohort record (39 `m6A` + 39
`RNA-seq`) is copied verbatim to `_data/_metadata/sample_sheet_ip_included.csv`. 39 rows deleted.

Consumers reading `sample_sheet.csv` all get the Input rows: `01_core_alignment` `01:26`, `02:20`,
`03:21`, `04:30`, `05:27` (FastQC, trimming, rRNA, STAR — parse `File1`/`File2`),
`02_dge_analysis` `01:42`, `02:34`, `03:205`, `04:37`, `06:31`, `03_variant_calling` `04b:28`,
`05:30`, `15:53`. Effects: STAR aligns 39 samples, so `run`'s `RAW_BAM_LIST` no longer feeds IP
samples into `01_read_groups.sh`, `02_mark_duplicates.sh`, `03_split_n_cigar.sh` (the wasted
preprocessing the Input-only confirmation found), and `06:35-40`'s 39-VCF gate matches the
39 gVCFs `04` produces. No delivered script gates on a 78-sample count — the `78` literals in
`15:190,192,282,286`, `16:72,74`, `17:116-348` are rule widths.

Consumers left reading `sample_sheet.csv` that want the full cohort — **reported, not fixed**:
`_hpc/backup.sbatch:30`, `_hpc/generate_checksum_csv.sbatch:31`, `_hpc/migrate.sbatch:30`,
`_hpc/rename.sbatch:30`, `_hpc/verify.sbatch:30`, `utils/qc02_cell_type_deconvolution.py:12`. `docs/_archive/pipeline_problem_list.md:38,110` prescribe adding
`Age`/`RIN` columns to the sheet.

## Consumers of the retired artifacts, left untouched

Peak-derived files, with no producer in the delivered tree:

- `_legacy/04_asm/01_find_drach_sites.R:26-27` — `consensus_peaks.bed` and `consensus_peaks.fa` (archived thesis arm)

Retired-script consumers, with no producer in the delivered tree:

- `environment.yml:85` — `bioconductor-rrho`, the RRHO2 dependency `10_rrho2_sex_comparison.R` bootstrapped
- `docs/_archive/pipeline_problem_list.md:61-64` — Item 6 prescribes RRHO2 on "male and female DMR ranked lists"; `:115` lists it as pipeline step 10
- `README.md:63-66` — names 4 `03_variant_calling/` scripts that do not exist (`06_functional_annotation.py`, `07_clinical_pathogenicity.py`, `08_pharmacogenomics_actionability.py`, `09_plot_clinical_annotation_landscape.R`; the real script is `07_functional_annotation.py`); `README.md:55-79` omits 5 real ones (`06_build_transcriptome_snp_table.py`, `07_functional_annotation.py`, `08_map_gtex_brain_eqtls.py`, `10_prioritize_mdd_candidate_snps.py`, `14_slice_igv_bams.sh`); `README.md:82-83` lists 2 of the 10 files in `docs/` (the `_legacy/04_asm/` claim at the old `:86` is gone with the `04_asm` block removal); `README.md:148` omits R libraries the delivered 02/03 use (`ashr` in `04_deg_analysis.R`, `clusterProfiler`/`org.Hs.eg.db` in `07`/`09`, `msigdbr` in `07`/`16b`/`17`/`20`, `ComplexHeatmap`/`circlize`/`RColorBrewer`/`tidyverse` in `09`, `ggtext` in `08`) — all reported, not fixed

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

## References to the archived thesis arm, reported not fixed

- `README.md:83` links `docs/asm_analysis_methodology.md`, which is now `docs/_archive/asm_analysis_methodology.md`
- `docs/` holds no MeRIP-era files: the MeRIP/ASM-era docs are in `docs/_archive/` (11 files), the retrospective narrative docs are in gitignored `_legacy/docs/` (3 files)
- `environment.yml:71-72` — `bioconductor-tress` (served the retired `08/09/10` peak scripts) and `bioconductor-rrho` (served the retired `10_rrho2_sex_comparison.R`), no consumer in the delivered tree
- `environment.yml:83` cites `_legacy/04_asm/scripts/07_validate_external.R`, not present in the archive (the archive holds the 14 files from `eb75c6e^`)
- `03_variant_calling/15_build_burden_matrix.py:16` — comment cites `04_asm/07`
- `utils/case_control_null.py` — absent from the tree, never tracked (`git ls-tree -r HEAD` and `git log --all -- utils/case_control_null.py` find nothing), so `17_gene_burden_test.py:52`'s `import case_control_null as ccn` fails at runtime and the helper is not recoverable from git; `docs/polygenic_burden_methodology.md:238,363` describe its `permute_within` role
- the 5 `_hpc/` consumers of the 78-row sheet, and `utils/qc02_cell_type_deconvolution.py` reading `sample_sheet.csv` (Input-only is correct for it) — held by explicit choice; the 20 MeRIP-era `utils/` files are in gitignored `_legacy/utils/`

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
| 2026-10-31 | `run:194-232` (40 lines deleted) | removed the HPC variant-calling submissions of `06_ase_readcounter.sbatch`, `07_m6allele_asm.sh`, `08_find_drach_sites.py`, `09_intersect_drach_snps.sh` (all absent from the tree; paths stale since `23ceab9`) and the `_data/logs/ase`, `_data/logs/asm`, `_data/logs/drach` mkdirs, so the `run vc` submission path ends at Step 05 VariantFiltration |
| 2026-10-31 | `04_asm/` (7 scripts, 1,676 lines, deleted), `_legacy/04_asm/`, `_legacy/04_asm/scripts/` (14 files, 1,950 lines) | thesis arm moved out of the tracked tree into gitignored `_legacy/04_asm/`; the read-ratio ASM generation from `eb75c6e^` archived under `_legacy/04_asm/scripts/`, folding in the 11 copies previously staged in `_legacy/03_variant_calling/` |
| 2026-10-31 | `_data/_metadata/sample_sheet.csv` (39 rows deleted), `_data/_metadata/sample_sheet_ip_included.csv` (new, 78 rows, untracked) | delivered sheet reduced to the 39 Input rows in original order; the 78-row cohort record copied verbatim to `sample_sheet_ip_included.csv`; 13 delivered consumers read the reduced sheet, the 7 `_hpc/` and `utils/` consumers that need the full cohort left reading it |
| 2026-10-31 | `README.md:81-90`, `README.md:148-167` (30 lines deleted) | removed the `04_asm/` module listing and section 4 "`m^6A`-SNP Analysis" command block (6 script invocations + the outputs paragraph), so the README describes 3 modules and its section headings are now 1, 2, 3 + Environment |
| 2026-10-31 | `docs/04_asm_methods.md`, `docs/04_asm_postmortem.md`, `docs/asm_analysis_methodology.md`, `docs/asm_filter_flow.md`, `docs/asm_skew_and_shift_calculation.md` (5 files, `git mv`) | thesis-arm methodology, postmortem, filter-flow, and skew/shift write-ups archived in `docs/_archive/`; recorded as `R` renames, contents unchanged |
| 2026-10-31 | `environment.yml` (1 insertion, 14 deletions) | removed `r-glmmtmb`, `r-lme4` (served the retired `12_asm_compare.R`), the m6A-SNP comment block, `bioconductor-biostrings` (DRACH motif matching, thesis arm only), and `r-r.utils` (retired `09_plot_m6asnp_funnel.R`); kept `r-data.table` — its comment named `04_annotate_strand_aware_disruption.R` but it is used by `03_variant_calling/16b`, `18`, `20-26` — and kept `r-ashr` (`02_dge_analysis/04_deg_analysis.R`); YAML verified with `yaml.safe_load`, 46 dependencies |
| 2026-10-31 | `docs/04_asm_literature.md` -> `docs/_archive/asm_literature.md` | recorded as observed, not performed here: the MeRIP/ASM literature review (234 lines) sits in `docs/_archive/` as `asm_literature.md`, content byte-identical to HEAD, git state `D` unstaged + untracked; the "MeRIP-era files in `docs/`" bullet was corrected as the only edit, and now says `docs/` holds no MeRIP-era files |
| 2026-10-31 | `docs/cohort_limitations_and_resolutions.md` (2 insertions, 36 deletions), `docs/_archive/cohort_limitations_and_resolutions.md` (new, untracked, verbatim 71 lines) | the registry rows and sections QC-03 (MACS2 narrowPeak audit, `_data/peaks/`, `../04_peak_calling/02_audit_gene_vs_peak.py`), QC-04 (antibody cross-reactivity, spike-in absence, methylation stoichiometry) and QC-05 (TRESS DMR covariate inflation, `utils/filter_dmr_lowcount.sh`) removed — all three measure MeRIP-dependent artifacts; the intro names the Input-only RNA-seq cohort and points at the archived copy; QC-01 (individual RIN/PMI/pH covariates not public) and QC-02 (signature-based deconvolution on bulk RNA-seq) are the retained items, both Input-only |
| 2026-10-31 | `docs/_archive/` + `m6a_modifications.md` (88), `rna_vs_merip.md` (120), `asm_literature.md` (234), `pipeline_problem_list.md` (126), `asm_variant_calling_summaries.md` (169), `asm_methods.md` (335), `asm_postmortem.md` (288); `docs/report.md` -> `docs/upstream_report.md` (31); `_legacy/docs/` + `client_report.md` (234), `variant_calling_raw_review.md` (162), `thesis_defense_narrative_log.md` (208), and `_legacy/implementation_plan.md` (96, from `docs/implementation_plan.md`) | MeRIP-era docs moved out of `docs/` into `docs/_archive/`, and the retrospective narrative docs plus the implementation plan into gitignored `_legacy/`; the 3 `04_asm_*.md` were renamed `asm_*.md` and `report.md` became `upstream_report.md`; `asm_variant_calling_summaries.md` is the variant-calling-with-ASM doc, moved in its merged state; contents verbatim from HEAD, `client_report.md` having left the tracked tree in `b890d2e`; git state is `D`/`RD` deletions in `docs/` plus `??` untracked copies |
| 2026-10-31 | `docs/cohort_limitations_and_resolutions.md` (4 lines rewritten) | QC-01 item 4's 12M exonic-count gate replaced by the 2M gene-mapped floor in `04_deg_analysis.R:95` (`MIN_LIBRARY_SIZE <- 2e6`) and the 39 delivered Input rows; item 5 and the QC-01 status line follow the 2M wording; QC-02 item 1 repointed from `scratch/cell_type_deconvolution.py` to `utils/qc02_cell_type_deconvolution.py` |
| 2026-10-31 | `docs/cohort_limitations_and_resolutions.md:33` (1 line) | QC-02 item 2's "all 36 gated subjects" replaced by "the 39 Input subjects" — 36 was the 12M gate's surviving count, and the 2M floor in `04_deg_analysis.R:95` keeps all 39; the measured p-values are unchanged |
| 2026-10-31 | `docs/polygenic_burden_methodology.md:431,432` (2 lines) | the two `Related documents` cites repointed from `docs/04_asm_methods.md` and `docs/04_asm_postmortem.md` to `docs/_archive/asm_methods.md` and `docs/_archive/asm_postmortem.md`, the paths as they now exist |
