# MDD RNA-seq & MeRIP-seq Analysis Pipeline

This repository hosts a modular, high-throughput bioinformatics pipeline for processing and analyzing paired-end human RNA-seq and MeRIP-seq (m6A IP) datasets. The study evaluates transcriptomic and epitranscriptomic profiles of Major Depressive Disorder (MDD) vs. healthy controls (`ctrl`), stratified by sex (`male` and `female`).

---

## 🎯 Objectives

1. **Gene Expression & Methylation Alterations (MDD vs. Control)**
   * **Baseline Expression (RNA-seq DEGs)**: Identify differentially expressed genes (upregulated/downregulated) from the standard RNA-seq (Input) data and annotate gene names.
   * **Methylation Evaluation (RIP-seq vs. RNA-seq)**: Detect differential methylation by comparing RIP-seq DEGs normalized by RNA-seq DEGs ($[\text{DEGs RIP}] / [\text{DEGs RNA}]$) between MDD and Control groups to identify hypo/hyper-methylation.
2. **Variant Identification (SNPs)**
   * Identify SNPs primarily within the RIP-seq (m6A) fraction. Because RIP-seq library features may limit coverage, additional SNP variant calling will be performed on the RNA-seq (Input) dataset for validation and comparative screening.

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
  * [`08_tress_peak_calling.R`](01_core_alignment/08_tress_peak_calling.R): Peak calling via TRESS in R.
  * [`09_merge_tress_peaks.R`](01_core_alignment/09_merge_tress_peaks.R): Merge peaks into cohort consensus coordinates.
  * [`10_quantification.sh`](01_core_alignment/10_quantification.sh): Peak-level read quantification via featureCounts.
  * [`11_counts_qc.sh`](01_core_alignment/11_counts_qc.sh): QC assertions on counts matrix.
  * [`12_multiqc.sh`](01_core_alignment/12_multiqc.sh): Centralized report aggregator.
* **[`02_dge_analysis/`](02_dge_analysis)**: Expression & methylation statistical testing.
  * [`config.yaml`](02_dge_analysis/config.yaml): Parameter configurations.
  * [`01_gene_mapping.R`](02_dge_analysis/01_gene_mapping.R): Ensembl database mapping via BioMart.
  * [`02_rna_processing.R`](02_dge_analysis/02_rna_processing.R): Low-count filtering and normalization.
  * [`03_pca_analysis.R`](02_dge_analysis/03_pca_analysis.R): Cohort dimensionality reduction.
  * [`04_deg_analysis.R`](02_dge_analysis/04_deg_analysis.R): Differential expression testing.
  * [`05_dge_plotting.R`](02_dge_analysis/05_dge_plotting.R): Generate Volcano/MA plots.
  * [`06_annotate_results.R`](02_dge_analysis/06_annotate_results.R): Map genes to functional labels.
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
* **[`04_asm/`](04_asm)**: Epitranscriptomic allele-specific methylation (ASM) and DRACH motif disruption.
  * [`01_ase_readcounter.sbatch`](04_asm/01_ase_readcounter.sbatch): GATK ASEReadCounter on IP BAMs.
  * [`02_input_ase_readcounter.sbatch`](04_asm/02_input_ase_readcounter.sbatch): ASEReadCounter on Input (baseline) BAMs.
  * [`03_m6allele_asm.sh`](04_asm/03_m6allele_asm.sh): M6Allele statistical ASM detection per subject.
  * [`04_find_drach_sites.py`](04_asm/04_find_drach_sites.py): Identify reference DRACH motifs in consensus peaks.
  * [`05_annotate_drach_disruption.py`](04_asm/05_annotate_drach_disruption.py) & [`06_intersect_drach_snps.sh`](04_asm/06_intersect_drach_snps.sh): Annotate SNPs disrupting DRACH motifs.
  * [`07_drach_analysis.sbatch`](04_asm/07_drach_analysis.sbatch): Slurm orchestrator for DRACH analysis.
  * [`08_asm_build_matrix.py`](04_asm/08_asm_build_matrix.py): Build tidy per-(locus, subject) allele count matrix.
  * [`09_asm_compare.R`](04_asm/09_asm_compare.R): Beta-binomial GLMM group comparison (MDD vs. Ctrl, sex-stratified).
  * [`10_asm_compare_viz.R`](04_asm/10_asm_compare_viz.R): Visualizations and diagnostic plots for ASM.
  * [`11_plot_drach_disruption.R`](04_asm/11_plot_drach_disruption.R): Cohort-wide DRACH motif disruption overview.
  * [`12_plot_sig_asm_drach.R`](04_asm/12_plot_sig_asm_drach.R): Fine-mapping figure of significant ASM DRACH loci.
  * [`13_bin3_tracks.ini`](04_asm/13_bin3_tracks.ini) & [`13_plot_igv_track_bin3.sh`](04_asm/13_plot_igv_track_bin3.sh): pyGenomeTracks and IGV report for BIN3 locus.
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
Rscript 01_core_alignment/08_tress_peak_calling.R
Rscript 01_core_alignment/09_merge_tress_peaks.R
bash 01_core_alignment/10_quantification.sh
bash 01_core_alignment/11_counts_qc.sh
bash 01_core_alignment/12_multiqc.sh
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

### 4. Allele-Specific Methylation (ASM) & DRACH Fine-Mapping
Run ASM quantification, beta-binomial modeling, and track generation:
```bash
# Allelic counting
sbatch 04_asm/01_ase_readcounter.sbatch
sbatch 04_asm/02_input_ase_readcounter.sbatch

# DRACH disruption analysis
sbatch 04_asm/07_drach_analysis.sbatch

# Matrix assembly & statistical modeling
python3 04_asm/08_asm_build_matrix.py
Rscript 04_asm/09_asm_compare.R
Rscript 04_asm/10_asm_compare_viz.R
Rscript 04_asm/11_plot_drach_disruption.R
Rscript 04_asm/12_plot_sig_asm_drach.R
bash 04_asm/13_plot_igv_track_bin3.sh
```

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
