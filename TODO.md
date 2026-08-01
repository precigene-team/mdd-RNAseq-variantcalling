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
* [x] Execute GATK pipeline (Steps 01-05) for all 78 BAMs on HPC via Slurm array.
* [x] Step 06 ASE ReadCounter completed for all 39 subjects.
* **Status**: **Completed**.

### 📊 5. Centralized MultiQC Aggregation
* [x] Create the centralized configuration [multiqc_config.yaml](multiqc_config.yaml) to organize reports chronologically.
* [x] Run the final MultiQC execution command to consolidate raw/trimmed FastQC, STAR alignment, Qualimap, and featureCounts metrics.
* **Status**: **Completed**. Final report generated at `_data/multiqc/rna_seq_pipeline_summary.html`.

### 🛠️ 6. Pipeline Refinements & Reorganization
* [x] Move cell-type deconvolution script from `scratch/cell_type_deconvolution.py` to `utils/cell_type_deconvolution.py` and redirect output scores to `results/utils/tables/cell_type_mgp_scores.csv`.
* [FAILED/DELEGATED_TO_DMR] Reaudit DGE analysis (specifically inspect `05_ma_male_meth.png` for potential cloud tilting or normalization biases).
* [ ] Refactor `02_dge_analysis/08_enrichment_plotting.R` to show divided facet headers for Activated/Suppressed pathways.
* [x] Format and clean up the TODO list.

### 📖 7. Learning & Literature Review
* **Bulk RNA-seq Mentorship (Active)**:
  * [x] Module 1: File Formats (FASTQ vs. FASTA vs. SAM/BAM) (Completed)
  * [x] Module 2: Normalization Math (TPM, FPKM, and DESeq2) (Completed)
  * [x] Module 3: Confounders & Cell-Type Deconvolution (Completed)
* **RIP-seq/m6A Epitranscriptomics Mentorship (Future)**:
  * [ ] Read RIP-seq/m6A review paper
    * [ ] Document the paper's DOI/URL link for tracking
    * [ ] Focus on the "Bulk RNA-seq vs. RIP-seq" transition:
      * Identify how RIP-seq normalization handles IP (immunoprecipitated) samples vs. Input/IgG controls, compared to standard Bulk RNA-seq normalization
      * Note why peak calling (finding local sequence enrichment) is required rather than gene-level expression quantification
    * [ ] Detail the metagene profile differences (3' UTR/stop codon enrichment in MeRIP-seq vs. uniform coverage in RNA-seq)
    * [ ] Establish benchmarks: QC (FRiP, SNR), Peak Caller logic (MACS2 vs. RIPSeeker/Piranha), and Biological Validation (DRACH motif)
    * [ ] Review m6A biology fundamentals (Writers, Erasers, Readers)
    * [ ] Extract and log key methodologies/comparison tables verbatim

---

## 🧬 Active Workstream (July 2026)

### Step 07 — M6Allele ASM
- [x] Fix: bcftools filter to biallelic het SNPs (crashed on multi-allelic)
- [x] Fix: graceful `if` wrap for Java crashes (IndexOutOfBounds, NumberFormatException)
- [x] Fix: auto-scale `-Xmx` from `SLURM_MEM_PER_NODE` (not hardcoded)
- [x] Fix: `--mem=16G` → `32G` for java heap headroom
- [ ] **IN PROGRESS** — fix committed, re-submission pending on HPC
- [ ] Count final success/failure and re-run any problematic subjects individually

### Steps 08-09 — DRACH Annotation
- [x] DRACH sites identified: 602,657 sites from 33,123 peaks
- [x] SNP×DRACH overlap: 211,284 total, 122,959 motif-disrupting
- [x] DRACH CSV pulled locally: `results/03_variant_calling/tables/09_drach_disrupting_snps.csv` (22 MB, 211K rows)

### TRESS DMR Analysis
- [x] Male DMR: 8,508 significant DMRs (completed)
- [x] Female DMR: 3,669 significant DMRs (completed)
- [x] `lowcount=15` post-hoc filter script: `utils/filter_dmr_lowcount.sh`
- [x] Filtered outputs: `results/utils/tables/10_dmr_{male,female}_filtered.csv` (8,495 / 3,590 DMRs)
- [x] QC-05 documented: missing RIN/age covariates inflate DMR counts vs paper (239/696)

### Downstream Integration (for presentation)
- [ ] Overlap m6A-SNPs (DRACH CSV) with DMR regions (male/female)
- [ ] Write integration script / summary table
- [ ] Prepare figures for mid-internship presentation

---

## 📋 Paper Replication Tracker (Mitsuhashi et al. 2024, iScience)

### DMR Parameter Alignment
| Parameter         | Paper                 | Ours                          | Match?                                  |
| ----------------- | --------------------- | ----------------------------- | --------------------------------------- |
| Tool              | TRESS                 | TRESS                         | ✅                                       |
| Binsize           | 50                    | 50                            | ✅                                       |
| lowcount          | 15                    | default (~10)                 | ⚠️ post-hoc filtered (negligible effect) |
| Covariates        | RIN + age             | None                          | ❌ not publicly available                |
| Threshold         | p<0.05, FC>1.2        | p<0.05, FC>1.2                | ✅                                       |
| Within-sex male   | 239 DMRs              | 8,495 (our unfiltered: 8,508) | ❌ covariate-driven gap                  |
| Within-sex female | 696 DMRs              | 3,590 (our unfiltered: 3,669) | ❌ covariate-driven gap                  |
| Interaction       | 629 DMRs (152M, 477F) | —                             | ❓ not run yet                           |

### Replication Steps
- [ ] Understand paper's DEG pipeline (DESeq2 + Subread, see STAR Methods)
- [ ] Run paper's sex-disease interaction DMR contrast
- [ ] Compare our DMR gene lists to paper's Supplementary Tables S3/S5
- [ ] GO enrichment on DMR genes (clusterProfiler, ViSEAGO)
- [ ] PPI network of shared DMR-DEG clusters (STRING + Cytoscape)
- [ ] RRHO analysis (male vs female methylation rank concordance)


---
# SNP Downstream

* [x] Re-check metadata, ensure each subject has correct Input/m6A-IP pair.
* [x] Use STAR-aligned BAMs from both Input and IP. (Steps 01-02 on HPC)
* [x] Run MarkDuplicates. (Step 02)
* [x] Run SplitNCigarReads separately for each BAM: (Step 03)
    * [x] All Input samples.
    * [x] All m6A-IP samples.
* [x] Verify split BAMs: no corruption, .bai index exists, correct read group and sample names.
* [x] Use only Input RNA-seq for SNP calling: (Steps 04-05)
    * [x] Run HaplotypeCaller for each Input sample. (per-sample GVCF mode — Step 04)
    * [x] Joint genotyping across all Input samples. (CombineGVCFs 04b + GenotypeGVCFs 04c)
    * [x] Keep only biallelic SNPs. (SelectVariants --restrict-alleles-to BIALLELIC — Step 05)
    * [x] Filter low-quality SNPs. (VariantFiltration FS > 30.0, QD < 2.0, -f PASS — Step 05)
    * [x] Extract heterozygous 0/1 SNPs per individual. (bcftools view -g het — Step 05)
    * [x] Prioritize SNPs with DP ≥ 10, GQ ≥ 20, and sufficient reads supporting both REF and ALT. (bcftools view -i 'FMT/DP >= 10 && FMT/GQ >= 20' — Step 05, 2.25M het sites across 39 subjects)
* [x] Call m6A peaks by comparing m6A-IP with Input:
    * [x] Each IP sample must be paired with the correct Input from the same individual.
    * [x] Used exomePeak2 per-subject.
    * [x] Export peak list in BED format.
    * [x] Create consensus peaks (`consensus_peaks.bed`, 33,123 peaks).
* [x] Intersect Input heterozygous SNP list with m6A peak regions: (Steps 08-09)
    * [x] Keep only het SNPs located within or near m6A-enriched regions.
    * [x] Do not call SNPs independently from IP data. (note: also annotated DRACH motifs — 211K SNP×motif overlaps, `09_drach_disrupting_snps.csv`)
* [x] At selected SNPs, count alleles in both BAMs from the same individual: (Step 06 — ASEReadCounter)
    * [x] REF count and ALT count in Input.
    * [x] REF count and ALT count in IP.
    * [x] Used GATK ASEReadCounter.
* [x] Calculate per-SNP metrics:
    * [x] REF fraction = REF / (REF + ALT).
    * [x] ALT fraction can be calculated similarly.
    * [x] Compare REF/ALT ratios between IP and Input for the same individual.
* [~] Only consider as allele-specific m6A candidate when: (Step 07 — M6Allele, in progress)
    * [ ] SNP is a reliable heterozygous call in Input.
    * [ ] Sufficient coverage in both Input and IP.
    * [ ] Allele ratio in IP is significantly skewed relative to Input.
    * [ ] SNP lies within an m6A peak region.
    * [ ] Not dominated by mapping bias, strand bias, or a handful of outlier reads.
* [ ] After obtaining per-subject ASm6A:
    * [ ] Compare allele imbalance levels between MDD and control groups.
    * [ ] Statistical model must adjust for at least sex.
    * [ ] Use ~ sex + condition; to test sex-specific differences, use ~ sex * condition.
* [ ] Interpretation notes:
    * [ ] A SNP lying in an m6A peak does NOT mean that SNP is the methylated nucleotide.
    * [ ] Only conclude that the allele carrying the SNP has a different enrichment level in the m6A region compared to the other allele.
    * [ ] Do not take a SNP that appears only in IP and conclude it is an m6A-associated SNP.