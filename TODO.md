# Project Analysis Roadmap & TODO List

This document lists the completed and remaining tasks for the MDD RNA-seq and MeRIP-seq epitranscriptomic analysis pipeline.

---

## 📋 Project Status Checklist

### 📈 1. Extended DGE & Interaction Contrasts
* [x] Modify/create R scripts to perform sex-specific and interaction-based contrasts (relative methylation change: `Condition:Fraction`).
* [x] Implement filter thresholds ($|\text{log2FC}| \ge 2$, $p < 0.05$).
* [x] Generate sex-stratified Volcano and MA plots.
* **Status**: **Completed**. Results saved in `results/02_dge_analysis/tables/` and `results/02_dge_analysis/figures/`.

### 🏷️ 2. Downstream Functional Annotations
* [x] Map DEG/DMG gene IDs to HGNC symbols, biotypes, and genomic coordinates.
* [x] Resolve and map subgroup SRA accession identifiers to cohort groups.
* **Status**: **Completed**. Mapping tables saved to `results/02_dge_analysis/tables/01_gene_mapper.csv` and `results/02_dge_analysis/tables/06_subgroup_sra_accessions.txt`.

### 🛡️ 3. Quality Control Registry (QC Resolutions)
* [x] **QC-01: Postmortem Confounders**: Calibrate and apply the 12M exonic counts gate (gated 36 subjects, filtered out `male_ctrl_1`).
* [x] **QC-02: Cell-Type Heterogeneity**: Run cell-type deconvolution (marker gene profiles) to verify glia/astrocyte loss (confirmed astrocyte reduction in MDD, $p = 0.0389$).
* [x] **QC-03: Gene-level vs. Peak-level Resolution**: Run MACS2 peak-calling on the 36 gated subjects and perform a cohort-wide overlap audit to validate gene body counts concordance.
* [x] **QC-04: Antibody occupancy limitations**: Document and interpret methylation shifts as relative interaction ratios only.
* **Status**: **Completed**. Documented in [cohort_limitations_and_resolutions.md](docs/cohort_limitations_and_resolutions.md).

### 🧬 4. Variant Calling (GATK Best Practices)
* [x] Write GATK processing scripts under `03_variant_calling/` (ReadGroups, MarkDuplicates, SplitNCigarReads, HaplotypeCaller, VariantFiltration).
* [x] Execute local worker pipeline script (`scratch/local_worker_pipeline.sh`) to run CPU-heavy SplitNCigarReads on aligned BAMs.
  * **Current State**: Active. Two local workers are currently processing SplitNCigarReads for the remaining cohort samples.
* [x] Execute HaplotypeCaller (`04_haplotype_caller.sh`) on the split BAMs to identify raw SNPs. 
* [ ] Apply VariantFiltration (`05_variant_filtration.sh`) to identify high-confidence m6A-SNPs and allele-specific methylation (ASM) patterns.
* **Status**: **In Progress** (HaplotypeCalling processing).

### 📊 5. Centralized MultiQC Aggregation
* [x] Create the centralized configuration [multiqc_config.yaml](multiqc_config.yaml) to organize reports chronologically.
* [x] Run the final MultiQC execution command to consolidate raw/trimmed FastQC, STAR alignment, Qualimap, and featureCounts metrics.
* **Status**: **Completed**. Final report generated at `_data/multiqc/rna_seq_pipeline_summary.html`.

- Update Peak Calling in README. 
- Format TODO
- Reaudit DGE analysis