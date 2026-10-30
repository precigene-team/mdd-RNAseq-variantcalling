# MDD RNA-seq Analysis Pipeline

This repository hosts a modular, high-throughput bioinformatics pipeline for processing and analyzing paired-end human RNA-seq (Input) datasets. The study evaluates the transcriptome and germline variation of Major Depressive Disorder (MDD) vs. healthy controls (`ctrl`), stratified by sex (`male` and `female`). The analysis focuses on the 39 Input samples (20 MDD, 19 control; 20 male, 19 female) and stops at variant calling: `01_core_alignment` → `02_dge_analysis` → `03_variant_calling`. 

MeRIP-seq (m6A IP) analysis is out of scope for this repository.

---

## 🎯 Objectives

1. **Gene Expression Alterations (MDD vs. Control)**
   * **Differential Expression (RNA-seq DEGs)**: Identify differentially expressed genes (upregulated/downregulated) from the RNA-seq (Input) data and annotate gene names.
2. **Variant Identification (SNPs)**
   * Call SNPs on the RNA-seq (Input) data via GATK joint calling across the cohort, then hard filtration, functional annotation, the MDD candidate review table, and the polygenic burden test.

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

* **[`03_variant_calling/`](03_variant_calling)**: GATK joint calling, functional annotation, and polygenic burden testing.
  * [`01_read_groups.sh`](03_variant_calling/01_read_groups.sh): Add read groups to aligned BAMs.
  * [`02_mark_duplicates.sh`](03_variant_calling/02_mark_duplicates.sh): Flag PCR duplicate alignments.
  * [`03_split_n_cigar.sh`](03_variant_calling/03_split_n_cigar.sh): Split reads spanning exon-exon junctions.
  * [`04_haplotype_caller.sh`](03_variant_calling/04_haplotype_caller.sh): Call raw variants using HaplotypeCaller in GVCF mode.
  * [`05_combine_gvcfs.sh`](03_variant_calling/05_combine_gvcfs.sh): Combine per-sample GVCFs into cohort GVCF.
  * [`06_genotype_gvcfs.sh`](03_variant_calling/06_genotype_gvcfs.sh): Joint-genotype cohort GVCF.
  * [`07_variant_filtration.sh`](03_variant_calling/07_variant_filtration.sh): Hard-filter SNPs and extract per-sample heterozygous calls.
  * [`08_build_transcriptome_snp_table.py`](03_variant_calling/08_build_transcriptome_snp_table.py): Aggregate the per-sample filtered VCFs into the cohort transcriptome SNP master table (PASS, 0/1, DP>=10, GQ>=20) with carrier counts by condition and sex.
  * [`09_functional_annotation.py`](03_variant_calling/09_functional_annotation.py): SnpEff GRCh38.115 annotation of recurrent SNPs (>=8 subjects): Sequence Ontology consequence, impact tier, HGVS c/p.
  * [`10_prioritize_mdd_candidate_snps.py`](03_variant_calling/10_prioritize_mdd_candidate_snps.py): MDD-enriched protein-altering candidates (`n_ctrl <= 2` & `n_mdd >= 8`) with the external evidence URLs. Terminal review table.
  * [`11_build_burden_matrix.py`](03_variant_calling/11_build_burden_matrix.py): Build the donor × variant carrier matrix for burden analysis. Data munging only.
  * [`12_fetch_opentargets_panel.py`](03_variant_calling/12_fetch_opentargets_panel.py): Fetch and cache the Open Targets MDD panel (genetic evidence only).
  * [`13_export_pathway_sets.R`](03_variant_calling/13_export_pathway_sets.R): Export msigdbr gene sets, collections kept separate, for the exploratory pathway arm.
  * [`14_gene_burden_test.py`](03_variant_calling/14_gene_burden_test.py): Four-arm burden test with detection-depth offset and within-sex permutation.
  * [`15_plot_burden_results.R`](03_variant_calling/15_plot_burden_results.R): Raw-vs-offset figure, Arm 1 volcano, Arm 2 panel.
  * [`16_burden_control_vs_mdd_table.R`](03_variant_calling/16_burden_control_vs_mdd_table.R): Per-donor and per-set burden tables. Draws nothing.
  * [`17_plot_burden_control_vs_mdd.R`](03_variant_calling/17_plot_burden_control_vs_mdd.R): Shared plot helper for the per-collection figures below. Sourced, not run.
  * [`18_plot_burden_reactome.R`](03_variant_calling/18_plot_burden_reactome.R): Control-vs-MDD burden, Reactome. Top-10 reference panel + Reactome synapse sets.
  * [`19_plot_burden_kegg.R`](03_variant_calling/19_plot_burden_kegg.R): Same layout, KEGG.
  * [`20_plot_burden_go_bp.R`](03_variant_calling/20_plot_burden_go_bp.R): Same layout, GO Biological Process.
  * [`21_plot_burden_go_cc.R`](03_variant_calling/21_plot_burden_go_cc.R): Same layout, GO Cellular Component.
  * [`22_plot_burden_go_mf.R`](03_variant_calling/22_plot_burden_go_mf.R): Same layout, GO Molecular Function.
  * Methodology: [`docs/polygenic_burden_methodology.md`](docs/polygenic_burden_methodology.md).

* **[`docs/`](docs)**: Documentation and background theory.
  * [`aspera_setup.md`](docs/aspera_setup.md): Aspera CLI setup and configuration for fast downloads on the HPC.
  * [`cohort_limitations_and_resolutions.md`](docs/cohort_limitations_and_resolutions.md): Gating thresholds and deconvolution results.
  * [`environment.md`](docs/environment.md): Full `mdd` package inventory, which script loads each package, and the rationale.
  * [`mapping_qc_qa.md`](docs/mapping_qc_qa.md): RNA-seq mapping and quality-control Q&A summary.
  * [`pipeline_integrity_checklist.md`](docs/pipeline_integrity_checklist.md): Data-integrity checks for hybrid local-HPC processing.
  * [`polygenic_burden_methodology.md`](docs/polygenic_burden_methodology.md): Design rationale for the `11`–`22` burden test, permutation offset, and stated limits.
  * [`split_n_cigar_qa.md`](docs/split_n_cigar_qa.md): SplitNCigar and HaplotypeCaller execution progress Q&A.
  * [`upstream_report.md`](docs/upstream_report.md): Upstream alignment and QC report.
  * **[`docs/_archive/`](docs/_archive)**: 14 archived MeRIP/TRESS-era documents, including `multiqc_summary.md` and the two QC verdict reports (156 FastQC reports across 78 libraries; the IP arm is out of scope, so those numbers cannot be reproduced from the 39 Input libraries).
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

### 3. Variant Calling, Annotation & Burden Testing
Joint genotyping on RNA-seq data, then the transcriptome SNP table, functional annotation, the MDD candidate review table, and the polygenic burden test:
```bash
# GATK Joint Genotyping Core
bash 03_variant_calling/01_read_groups.sh
bash 03_variant_calling/02_mark_duplicates.sh
bash 03_variant_calling/03_split_n_cigar.sh
bash 03_variant_calling/04_haplotype_caller.sh
sbatch 03_variant_calling/05_combine_gvcfs.sh
sbatch 03_variant_calling/06_genotype_gvcfs.sh
sbatch 03_variant_calling/07_variant_filtration.sh

# Transcriptome SNP table, SnpEff annotation, candidate URL table
python3 03_variant_calling/08_build_transcriptome_snp_table.py
python3 03_variant_calling/09_functional_annotation.py
python3 03_variant_calling/10_prioritize_mdd_candidate_snps.py

# Polygenic Burden Test
python3 03_variant_calling/11_build_burden_matrix.py
python3 03_variant_calling/12_fetch_opentargets_panel.py
Rscript 03_variant_calling/13_export_pathway_sets.R
python3 03_variant_calling/14_gene_burden_test.py
Rscript 03_variant_calling/15_plot_burden_results.R
Rscript 03_variant_calling/16_burden_control_vs_mdd_table.R
Rscript 03_variant_calling/18_plot_burden_reactome.R
Rscript 03_variant_calling/19_plot_burden_kegg.R
Rscript 03_variant_calling/20_plot_burden_go_bp.R
Rscript 03_variant_calling/21_plot_burden_go_cc.R
Rscript 03_variant_calling/22_plot_burden_go_mf.R
# 17_plot_burden_control_vs_mdd.R is the shared plot helper sourced by 18-22; not run directly.
# The same 08-22 chain in one command: ./run vc-post runs it locally; ./run vc submits
# it as the vc_annotation_burden Slurm job after the filtration job.
```

---

## 🛠️ Environment & Prerequisites
Activate the custom Mamba environment configured for this analysis:
```bash
mamba activate mdd
```
Full package inventory, which script loads each one, and the rationale: [`docs/environment.md`](docs/environment.md).
Required software includes:
* **Java SDK** (v1.8+ or compatible runtime)
* **GATK Toolkit** (v4.x)
* **R Libraries**: `tidyverse` (ggplot2, dplyr, readr, stringr), `DESeq2`, `biomaRt`, `limma`, `clusterProfiler`, `org.Hs.eg.db`, `ashr`, `msigdbr`, `data.table`, `yaml`, `ggtext`, `ggrepel`, `patchwork`, `ComplexHeatmap`, `circlize`, `RColorBrewer`, `scales`, `plotly`, `htmlwidgets`
* **Python Libraries**: `pandas`, `pathlib`, `concurrent.futures`, `subprocess`
