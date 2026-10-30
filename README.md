# MDD RNA-seq Analysis Pipeline

This repository hosts a modular, high-throughput bioinformatics pipeline for processing and analyzing paired-end human RNA-seq (Input) datasets. The study evaluates the transcriptome and germline variation of Major Depressive Disorder (MDD) vs. healthy controls (`ctrl`), stratified by sex (`male` and `female`). The analysis focuses on the 39 Input samples (20 MDD, 19 control; 20 male, 19 female) and stops at variant calling: `01_core_alignment` → `02_dge_analysis` → `03_variant_calling`. 

MeRIP-seq (m6A IP) analysis is out of scope for this repository.

---

## 🎯 Objectives

1. **Gene Expression Alterations (MDD vs. Control)**
   * **Differential Expression (RNA-seq DEGs)**: Identify differentially expressed genes (upregulated/downregulated) from the RNA-seq (Input) data and annotate gene names.
2. **Variant Identification (SNPs)**
   * Call SNPs on the RNA-seq (Input) data via GATK joint calling across the cohort, then hard filtration, functional annotation, and clinical actionability assessment.

---

## 🧬 Library Preparation Details

* **Organism**: *Homo sapiens* (postmortem brain tissue, prefrontal cortex)
* **Extracted Molecule**: Total RNA (fragmented to ~200nt)
* **Library Prep Kit**: SMARTer Stranded Total RNA-seq Kit v2 (Takara, #634412)
* **Strandedness**: Reverse-stranded (`-s 2` parameter in `featureCounts`), as Read 1 is antisense to the target RNA.
* **Sequencing Platform**: Illumina NovaSeq 6000 (Paired-end)

---

## 📂 Repository Structure

Below is the layout of the analysis workspace:

* **[`01_core_alignment/`](01_core_alignment)**: Core sequence processing pipeline.
  * [`01_raw_fastqc.sh`](01_core_alignment/01_raw_fastqc.sh): Quality control check on raw FastQ files.
  * [`02_trimming.sh`](01_core_alignment/02_trimming.sh): Trimmomatic adapter and quality trimming.
  * [`03_rrna_filter.sh`](01_core_alignment/03_rrna_filter.sh): Ribosomal RNA filtering.
  * [`04_trimmed_fastqc.sh`](01_core_alignment/04_trimmed_fastqc.sh): Post-trimming quality control.
  * [`05_star_align.sh`](01_core_alignment/05_star_align.sh): Map trimmed reads to the human reference.
  * [`06_qualimap_bamqc.sh`](01_core_alignment/06_qualimap_bamqc.sh) & [`07_qualimap_rnaseq.sh`](01_core_alignment/07_qualimap_rnaseq.sh): QC assessment on BAM alignments.
  * [`08_quantification.sh`](01_core_alignment/08_quantification.sh): Gene-level read quantification via featureCounts against the Ensembl GTF.
  * [`09_counts_qc.sh`](01_core_alignment/09_counts_qc.sh): QC assertions on counts matrix.
  * [`10_multiqc.sh`](01_core_alignment/10_multiqc.sh): Centralized report aggregator.

* **[`02_dge_analysis/`](02_dge_analysis)**: Expression statistical testing.
  * [`config.yaml`](02_dge_analysis/config.yaml): Ensembl annotation parameters (`ensembl_version`, `valid_biotypes`).
  * [`01_gene_mapping.R`](02_dge_analysis/01_gene_mapping.R): ENSG → HGNC symbol mapping from the Ensembl GTF, Input samples only.
  * [`02_rna_processing.R`](02_dge_analysis/02_rna_processing.R): Low-count filtering and normalization.
  * [`03_pca_analysis.R`](02_dge_analysis/03_pca_analysis.R): Cohort dimensionality reduction.
  * [`04_deg_analysis.R`](02_dge_analysis/04_deg_analysis.R): Differential expression testing.
  * [`05_dge_plotting.R`](02_dge_analysis/05_dge_plotting.R): Generate Volcano/MA plots.
  * [`06_annotate_results.R`](02_dge_analysis/06_annotate_results.R): Map genes to functional labels.
  * [`07_enrichment_analysis.R`](02_dge_analysis/07_enrichment_analysis.R): GO/KEGG over-representation and GSEA on significant genes.
  * [`08_enrichment_plotting.R`](02_dge_analysis/08_enrichment_plotting.R): GSEA NES barplots per contrast.
  * [`09_heatmap_pathway.R`](02_dge_analysis/09_heatmap_pathway.R): Leading-edge pathway heatmaps (row z-score).

* **[`03_variant_calling/`](03_variant_calling)**: GATK joint calling and clinical translation.
  * [`01_read_groups.sh`](03_variant_calling/01_read_groups.sh): Add read groups to aligned BAMs.
  * [`02_mark_duplicates.sh`](03_variant_calling/02_mark_duplicates.sh): Flag PCR duplicate alignments.
  * [`03_split_n_cigar.sh`](03_variant_calling/03_split_n_cigar.sh): Split reads spanning exon-exon junctions.
  * [`04_haplotype_caller.sh`](03_variant_calling/04_haplotype_caller.sh): Call raw variants using HaplotypeCaller in GVCF mode.
  * [`04b_combine_gvcfs.sh`](03_variant_calling/04b_combine_gvcfs.sh): Combine per-sample GVCFs into cohort GVCF.
  * [`04c_genotype_gvcfs.sh`](03_variant_calling/04c_genotype_gvcfs.sh): Joint-genotype cohort GVCF.
  * [`05_variant_filtration.sh`](03_variant_calling/05_variant_filtration.sh): Hard-filter SNPs and extract per-sample heterozygous calls.
  * [`06_functional_annotation.py`](03_variant_calling/06_functional_annotation.py): Sequence Ontology, CADD, SpliceAI, and PsychENCODE/PGC domain mapping.
  * [`07_clinical_pathogenicity.py`](03_variant_calling/07_clinical_pathogenicity.py): ClinVar cross-referencing and ACMG/AMP 5-tier classification.
  * [`08_pharmacogenomics_actionability.py`](03_variant_calling/08_pharmacogenomics_actionability.py): CPIC antidepressant guidelines and ClinPGx actionability.
  * [`09_plot_clinical_annotation_landscape.R`](03_variant_calling/09_plot_clinical_annotation_landscape.R): Multi-panel publication figure.
  * [`15_build_burden_matrix.py`](03_variant_calling/15_build_burden_matrix.py): Build the donor × variant carrier matrix for burden analysis. Data munging only.
  * [`16_fetch_opentargets_panel.py`](03_variant_calling/16_fetch_opentargets_panel.py): Fetch and cache the Open Targets MDD panel (genetic evidence only).
  * [`16b_export_pathway_sets.R`](03_variant_calling/16b_export_pathway_sets.R): Export msigdbr gene sets, collections kept separate, for the exploratory pathway arm.
  * [`17_gene_burden_test.py`](03_variant_calling/17_gene_burden_test.py): Four-arm burden test with detection-depth offset and within-sex permutation.
  * [`18_plot_burden_results.R`](03_variant_calling/18_plot_burden_results.R): Raw-vs-offset figure, Arm 1 volcano, Arm 2 panel.
  * [`20_burden_control_vs_mdd_table.R`](03_variant_calling/20_burden_control_vs_mdd_table.R): Per-donor and per-set burden tables. Draws nothing.
  * [`21_plot_burden_control_vs_mdd.R`](03_variant_calling/21_plot_burden_control_vs_mdd.R): Shared plot helper for the per-collection figures below. Sourced, not run.
  * [`22_plot_burden_reactome.R`](03_variant_calling/22_plot_burden_reactome.R): Control-vs-MDD burden, Reactome. Top-10 reference panel + Reactome synapse sets.
  * [`23_plot_burden_kegg.R`](03_variant_calling/23_plot_burden_kegg.R): Same layout, KEGG.
  * [`24_plot_burden_go_bp.R`](03_variant_calling/24_plot_burden_go_bp.R): Same layout, GO Biological Process.
  * [`25_plot_burden_go_cc.R`](03_variant_calling/25_plot_burden_go_cc.R): Same layout, GO Cellular Component.
  * [`26_plot_burden_go_mf.R`](03_variant_calling/26_plot_burden_go_mf.R): Same layout, GO Molecular Function.
  * Methodology: [`docs/polygenic_burden_methodology.md`](docs/polygenic_burden_methodology.md).

* **[`04_asm/`](04_asm)**: Sequence-level $m^6A$-SNP analysis — DRACH motif disruption in consensus m6A peaks, plus per-donor allele-specific m6A.
  * [`01_find_drach_sites.R`](04_asm/01_find_drach_sites.R): Scan consensus m6A peaks for the DRACH writer motif (Biostrings `vmatchPattern`, strand-aware).
  * [`02_verify_drach_reference.py`](04_asm/02_verify_drach_reference.py): Hard QC gate on the motif set; exits non-zero on corruption.
  * [`03_intersect_drach_genotypes.sh`](04_asm/03_intersect_drach_genotypes.sh): Join DRACH loci to cohort genotypes across all 39 donors; delegates to 02 then 04.
  * [`04_annotate_strand_aware_disruption.R`](04_asm/04_annotate_strand_aware_disruption.R): Classify each variant × motif instance as DISRUPTS/PRESERVES DRACH, with strand-aware complementation.
  * [`05_plot_m6asnp_funnel.py`](04_asm/05_plot_m6asnp_funnel.py): Single derivation point for the six-level funnel and the core $m^6A$ set.
  * [`06_test_mdd_enrichment.py`](04_asm/06_test_mdd_enrichment.py): Pooled MDD-vs-Ctrl test with donor-label permutation within sex.
  * [`07_allele_specific_m6a.py`](04_asm/07_allele_specific_m6a.py): Per-donor IP-vs-Input allele-specific m6A test; concordant sites called at >90% cross-donor agreement.
  * Retired read-ratio ASM scripts (M6Allele / ASEReadCounter / glmmTMB) are archived in [`_legacy/04_asm/`](_legacy/04_asm).
  * Methods: [`docs/04_asm_methods.md`](docs/04_asm_methods.md). Retrospective: [`docs/04_asm_postmortem.md`](docs/04_asm_postmortem.md).
* **[`docs/`](docs)**: Documentation and background theory.
  * [`cohort_limitations_and_resolutions.md`](docs/cohort_limitations_and_resolutions.md): Gating thresholds and deconvolution results.
  * [`asm_analysis_methodology.md`](docs/asm_analysis_methodology.md): ASM group comparison methodology and statistical models.
* **[`_data/`](_data)**: Reference annotation files, genome indices, and temporary directories.

---

## 🚀 Execution Workflow

### 1. Alignment & Quantification
Execute the pipeline shell scripts sequentially:
```bash
bash 01_core_alignment/01_raw_fastqc.sh
bash 01_core_alignment/02_trimming.sh
bash 01_core_alignment/03_rrna_filter.sh
bash 01_core_alignment/04_trimmed_fastqc.sh
bash 01_core_alignment/05_star_align.sh
bash 01_core_alignment/06_qualimap_bamqc.sh
bash 01_core_alignment/07_qualimap_rnaseq.sh
bash 01_core_alignment/08_quantification.sh
bash 01_core_alignment/09_counts_qc.sh
bash 01_core_alignment/10_multiqc.sh
```

### 2. Differential Expression (DGE) & Annotations
Map counts to HGNC symbols, filter low count noise, normalize, and run tests:
```bash
Rscript 02_dge_analysis/01_gene_mapping.R
Rscript 02_dge_analysis/02_rna_processing.R
Rscript 02_dge_analysis/03_pca_analysis.R
Rscript 02_dge_analysis/04_deg_analysis.R
Rscript 02_dge_analysis/05_dge_plotting.R
Rscript 02_dge_analysis/06_annotate_results.R
Rscript 02_dge_analysis/07_enrichment_analysis.R
Rscript 02_dge_analysis/08_enrichment_plotting.R
Rscript 02_dge_analysis/09_heatmap_pathway.R
```

### 3. Variant Calling & Clinical Translation
Joint genotyping on RNA-seq data followed by 3-stage clinical translation:
```bash
# GATK Joint Genotyping Core
bash 03_variant_calling/01_read_groups.sh
bash 03_variant_calling/02_mark_duplicates.sh
bash 03_variant_calling/03_split_n_cigar.sh
bash 03_variant_calling/04_haplotype_caller.sh
sbatch 03_variant_calling/04b_combine_gvcfs.sh
sbatch 03_variant_calling/04c_genotype_gvcfs.sh
sbatch 03_variant_calling/05_variant_filtration.sh

# Clinical Translation & Pharmacogenomics
python3 03_variant_calling/06_functional_annotation.py
python3 03_variant_calling/07_clinical_pathogenicity.py
python3 03_variant_calling/08_pharmacogenomics_actionability.py
Rscript 03_variant_calling/09_plot_clinical_annotation_landscape.R
```

### 4. $m^6A$-SNP Analysis & Allele-Specific m6A
Run in numeric order — script number == output number == execution order:
```bash
micromamba activate mdd

Rscript 04_asm/01_find_drach_sites.R            # ~8 s
python3 04_asm/02_verify_drach_reference.py      # ~40 s, hard gate
bash    04_asm/03_intersect_drach_genotypes.sh   # runs 02 + 04 internally
python3 04_asm/05_plot_m6asnp_funnel.py          # ~3 s
python3 04_asm/06_test_mdd_enrichment.py         # ~9 s (20k permutations)
python3 04_asm/07_allele_specific_m6a.py         # ~36 s

# tighter null for the enrichment test
N_PERM=100000 python3 04_asm/06_test_mdd_enrichment.py

# tighten the allele-specific concordance gate
CONCORDANCE=0.95 MIN_ALLELE_READS=3 python3 04_asm/07_allele_specific_m6a.py
```
Outputs land in `results/04_asm/tables/` (git-ignored; regenerate with the commands above).

---

## 🛠️ Environment & Prerequisites
Activate the custom Mamba environment configured for this analysis:
```bash
mamba activate mdd
```
Required software includes:
* **Java SDK** (v1.8+ or compatible runtime)
* **GATK Toolkit** (v4.x)
* **R Libraries**: `DESeq2`, `biomaRt`, `dplyr`, `yaml`, `ggplot2`, `ggrepel`, `plotly`, `htmlwidgets`
* **Python Libraries**: `pandas`, `pathlib`, `concurrent.futures`, `subprocess`
