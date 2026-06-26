# Project Analysis Roadmap & TODO list

Below are the outstanding tasks required to complete the MDD RNA-seq and MeRIP-seq epitranscriptomic analysis.

---

## 📋 Outstanding Tasks

### 📈 1. Implement Extended DGE & Interaction Contrasts
* [ ] Modify [`02_dge_analysis/deg.R`](file:///mnt/22T/phongdinh/RNAseq_pipeline_depression/02_dge_analysis/deg.R) or create a new script to perform sex-specific and interaction-based contrasts:
  * Male MDD vs Control (Input)
  * Female MDD vs Control (Input)
  * Male MDD vs Control (m6A IP)
  * Female MDD vs Control (m6A IP)
  * Interaction (Relative Methylation change: `Condition:Fraction`) for Males, Females, and overall cohort.
* [ ] Implement threshold filtering prioritizing $|\text{log2FC}| \ge 2$ and $p < 0.05$.
* [ ] Save the resulting differential tables and generate sex-stratified Volcano and MA plots.

### 🏷️ 2. Integrate Downstream Functional Annotations
* [ ] Implement/extend [`02_dge_analysis/gene_map.R`](file:///mnt/22T/phongdinh/RNAseq_pipeline_depression/02_dge_analysis/gene_map.R) to map DEG gene IDs to extra functional descriptors:
  * HGNC symbol & Ensembl Gene ID.
  * Gene Biotype (confirm protein coding or explore non-coding RNAs).
  * Genomic coordinates.
  * SRA run accession identifiers linked to each cohort subgroup.

### 🧬 3. Variant Calling & BQSR Calibration Verification
* [ ] Check local availability of GATK known variant databases (e.g., dbSNP, 1000 Genomes VCFs) on the server.
  * If reference datasets are missing, implement variant calling skipping the BQSR step, or document the download commands.
* [ ] Set up the GATK variant calling scripts under [`03_variant_calling/`](file:///mnt/22T/phongdinh/RNAseq_pipeline_depression/03_variant_calling):
  * Write `SplitNCigarReads.sh` to handle splice junction splitting.
  * Write `HaplotypeCaller.sh` to run parallel variant calling on processed BAM files.
  * Write `VariantFiltration.sh` to filter raw SNPs according to RNA-seq filtering criteria.
* [ ] Run SNP variant calling on both:
  * Input RNA-seq samples.
  * MeRIP-seq m6A IP samples.

### 📊 4. Centralized MultiQC & BAM QC
* [ ] Move previous individual task MultiQC results (from `_data/qc_raw/multiqc/` and `_data/qc_trim/multiqc/`) to `_legacy/`.
* [ ] Keep the `qc_raw`, `qc_trim` FastQC results, and counts.
* [ ] Create a centralized `multiqc_config.yaml` specifying section order chronologically (FastQC -> STAR -> Picard/Qualimap -> featureCounts -> Biotypes).
* [ ] Aggregate all reports (raw FastQC, trimmed FastQC, STAR logs, Qualimap QC, featureCounts) into a single consolidated MultiQC report.


