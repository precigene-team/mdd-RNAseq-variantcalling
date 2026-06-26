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
<!-- 
Below is the layout of the analysis workspace:

* **[`01_core_alignment/`](file:///mnt/22T/phongdinh/RNAseq_pipeline_depression/01_core_alignment)**: Core sequence processing pipeline.
  * [`fastqc/`](file:///mnt/22T/phongdinh/RNAseq_pipeline_depression/01_core_alignment/fastqc): Sequence quality check.
    * [`initial_qc.py`](file:///mnt/22T/phongdinh/RNAseq_pipeline_depression/01_core_alignment/fastqc/initial_qc.py): Parallelized pre-trim QC.
    * [`after_trim_QC.py`](file:///mnt/22T/phongdinh/RNAseq_pipeline_depression/01_core_alignment/fastqc/after_trim_QC.py): Parallelized post-trim QC.
  * [`trim/`](file:///mnt/22T/phongdinh/RNAseq_pipeline_depression/01_core_alignment/trim): Adapter and quality trimming.
    * [`trimmomatic.py`](file:///mnt/22T/phongdinh/RNAseq_pipeline_depression/01_core_alignment/trim/trimmomatic.py): Multi-threaded adapter trimming via Trimmomatic.
  * [`star/`](file:///mnt/22T/phongdinh/RNAseq_pipeline_depression/01_core_alignment/star): Genome alignment.
    * [`star.py`](file:///mnt/22T/phongdinh/RNAseq_pipeline_depression/01_core_alignment/star/star.py): Alignment of trimmed paired-end reads to the human reference.
  * [`qualimap/`](file:///mnt/22T/phongdinh/RNAseq_pipeline_depression/01_core_alignment/qualimap): Alignment assessment.
    * [`qualimap_bamqc.py`](file:///mnt/22T/phongdinh/RNAseq_pipeline_depression/01_core_alignment/qualimap/qualimap_bamqc.py) & [`qualimap_rnaseq.py`](file:///mnt/22T/phongdinh/RNAseq_pipeline_depression/01_core_alignment/qualimap/qualimap_rnaseq.py): BAM checks.
  * [`count/`](file:///mnt/22T/phongdinh/RNAseq_pipeline_depression/01_core_alignment/count): Gene abundance quantification.
    * [`featurecounts.py`](file:///mnt/22T/phongdinh/RNAseq_pipeline_depression/01_core_alignment/count/featurecounts.py): Reads mapping counts matrix creation.
* **[`02_dge_analysis/`](file:///mnt/22T/phongdinh/RNAseq_pipeline_depression/02_dge_analysis)**: Expression & methylation statistical testing.
  * [`config.yaml`](file:///mnt/22T/phongdinh/RNAseq_pipeline_depression/02_dge_analysis/config.yaml): Parameter configurations.
  * [`gene_map.R`](file:///mnt/22T/phongdinh/RNAseq_pipeline_depression/02_dge_analysis/gene_map.R): Ensembl database mapping via BioMart.
  * [`rna_processing.R`](file:///mnt/22T/phongdinh/RNAseq_pipeline_depression/02_dge_analysis/rna_processing.R): Low-count filtering, variance stabilization (VST/rlog), and 3D PCA plot generation.
  * [`deg.R`](file:///mnt/22T/phongdinh/RNAseq_pipeline_depression/02_dge_analysis/deg.R): Differential expression testing via DESeq2.
* **[`03_variant_calling/`](file:///mnt/22T/phongdinh/RNAseq_pipeline_depression/03_variant_calling)**: GATK variant calling on RNA-seq data.
  * [`ReadGroups.txt`](file:///mnt/22T/phongdinh/RNAseq_pipeline_depression/03_variant_calling/ReadGroups.txt): Read-group addition.
  * [`MarkDuplicates.txt`](file:///mnt/22T/phongdinh/RNAseq_pipeline_depression/03_variant_calling/MarkDuplicates.txt): Duplicate flagging.
* **[`docs/`](file:///mnt/22T/phongdinh/RNAseq_pipeline_depression/docs)**: Documentation and background theory.
  * [`rna_vs_merip.md`](file:///mnt/22T/phongdinh/RNAseq_pipeline_depression/docs/rna_vs_merip.md): Detailed comparison of RNA-seq and MeRIP-seq.
* **[`_data/`](file:///mnt/22T/phongdinh/RNAseq_pipeline_depression/_data)**: Reference annotation files, genome indices, and temporary directories.
* **[`utils/`](file:///mnt/22T/phongdinh/RNAseq_pipeline_depression/utils)**: Common scripts for parsing log metrics and resolving mapping anomalies. -->

---

## 🚀 Execution Workflow

<!-- ### 1. Data Retrieval and Migration
Run the batch migration script to transfer files from the local raw directory:
```bash
sbatch migrate.sbatch
```

### 2. Alignment & Quantification
Execute the Python automation scripts in [`01_core_alignment/`](file:///mnt/22T/phongdinh/RNAseq_pipeline_depression/01_core_alignment) sequentially:
```bash
python 01_core_alignment/fastqc/initial_qc.py
python 01_core_alignment/trim/trimmomatic.py
python 01_core_alignment/fastqc/after_trim_QC.py
python 01_core_alignment/star/star.py
python 01_core_alignment/qualimap/qualimap_rnaseq.py
python 01_core_alignment/count/featurecounts.py
```

### 3. Differential Expression (DGE) & Annotations
Map counts to HGNC symbols, filter low count noise, normalize, and run tests:
```bash
Rscript 02_dge_analysis/gene_map.R
Rscript 02_dge_analysis/rna_processing.R
Rscript 02_dge_analysis/deg.R
```

### 4. Variant Calling
Configure BAM file structure and flag PCR duplication artifacts:
```bash
bash 03_variant_calling/ReadGroups.txt
bash 03_variant_calling/MarkDuplicates.txt
``` -->

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
