---
title: "Allele-Specific $m^6A$ Epitranscriptomic Signal Identification in Major Depressive Disorder (MDD) Patients"
date: "2026-09-26"
---

# Intro

## Context

- **MDD Burden**: ~300M people globally; 30–40% treatment-resistant; high relapse rate.
- **Target Brain Region**: Prefrontal cortex (vmPFC / Brodmann Area 9) — key executive hub regulating emotion and mood.
- **Transcriptomics vs Epigenetics**:
  - Bulk RNA-seq GWAS finds risk loci, but >90% are non-coding regulatory variants.
  - Steady-state mRNA expression misses post-transcriptional dynamics and protein turnover.
- **Why $m^6A$?**:
  - Neurons are highly polarized (long axons/dendrites) -> need local mRNA storage rather than immediate somatic translation.
  - $m^6A$ tag acts as a molecular switch for mRNA localization, stability, and activity-dependent local translation at synaptic spines.
  - Controlled by **Writers** (`METTL3/14`, `WTAP`) that deposit tags, **Erasers** (`FTO`, `ALKBH5`) that remove them, and **Readers** (`YTHDF1-3`) that trigger translation or decay.
  - **The SNP Disruption Mechanism**: A mutation inside the `DRACH` motif prevents Writers from adding the tag -> Readers never bind -> mRNA fails to be properly translated or stabilized at the synapse.
  - Stress reprograms $m^6A$ marks on neuroplasticity genes.

## Problem

- **MDD Molecular Signals Hard to Pinpoint**:
  - Multifactorial etiology; postmortem human brain tissue is scarce ($N=39$).
  - Bulk transcriptomic signals are fuzzy/weak and confounded by cell-type shifts.
- **Standard MeRIP-seq Limitations (Fact-Checked)**:
  - Standard pipelines (exomePeak, TRESS) stop at **Differential Methylation Regions (DMRs)** between pooled Cases vs Controls.
  - Bulk DMR convolved: can't distinguish true methylation shifts from transcript abundance changes.
  - Ignores individual genetic variation: treats cohorts in aggregate, blind to personal heterozygous alleles.

## Solution / Goal

- **Allele-Specific Methylation (ASM) Pipeline**:
  - Combine paired Input RNA-seq + MeRIP-seq IP at single-nucleotide resolution.
  - Call personal heterozygous SNPs (`0/1`) from unbiased Input reads.
  - Identify $m^6A$-SNPs that physically break canonical **`DRACH` motifs** ($D=A/G/U, R=A/G, A=m^6A, C=C, H=A/C/U$).
  - Measure allelic skew in IP vs Input to detect preferential methylation.
- **Target Prioritization & Downstream Testing**:
  - Discover high-confidence disease-associated and sex-dimorphic ASM genes.
  - Pinpointing broken `DRACH` sites directly justifies downstream functional assays to test target rescue and translation.
  - In silico disruption modeling (AlphaGenome / AlphaFold).

<!-- pagebreak -->

# Method

- **Cohort / Data**:
  - Public human postmortem vmPFC MeRIP-seq ([GSE275676](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE275676), Mitsuhashi et al. 2024).
  - 39 donors, 78 paired libraries (39 Input, 39 IP):
    - Male MDD: $n = 10$ | Male Control: $n = 10$
    - Female MDD: $n = 10$ | Female Control: $n = 9$
- **Tech Stack**:
  - **Bash / SLURM**: Pipeline orchestration, HPC job management.
  - **Alignment & QC**: `fastp`, `STAR` (rRNA pre-mapping depletion & GRCh38 alignment), `featureCounts`.
  - **Variant Calling**: GATK 4.x (`SplitNCigarReads`, `HaplotypeCaller` in GVCF mode, `GenotypeGVCFs`).
  - **ASM Modeling**: `ASEReadCounter`, R (`glmmTMB`, `DESeq2`), Python (`pandas`, `pyGenomeTracks`).
- **Workflow Pipeline**:
  1. *QC & Depletion*: Adapter trimming (`fastp`) -> `STAR` rRNA pre-mapping depletion.
  2. *Peak Calling*: TRESS across cohort -> 33,123 consensus $m^6A$ peaks.
  3. *Unbiased Variant Discovery*: Run GATK **on Input libraries only** -> prevents antibody IP capture bias.
  4. *Motif Overlap*: Scan peaks for 5-nt `DRACH` -> classify SNPs into **Disrupting** vs **Preserving**.
  5. *ASM Quantification*: GATK `ASEReadCounter` on Input and IP BAMs at all het sites.
  6. *Statistical Model*: Beta-Binomial GLMM (`glmmTMB`):
     $$\text{cbind}(\text{IP}_{\text{REF}}, \text{IP}_{\text{ALT}}) \sim \text{Condition} + \text{Sex} + \text{offset}(\text{logit}(\text{Input}_{\text{AF}}))$$
     - Per-subject Input AF offset absorbs reference-mapping bias.
     - Tests condition effect while adjusting for sex and overdispersion.

# Results (Preliminary)

## Pipeline & Epitranscriptomic Filtering

- **Variant Funnel**: Identified millions of expressed heterozygous calls $\to$ thousands of peak-overlapping loci $\to$ high-frequency direct disruption of the central `DRACH` recognition core.
- **Statistical Gating**: Beta-Binomial GLMM with personal Input offset successfully isolated a tight set of high-confidence ASM loci ($\text{FDR} < 0.05$) from noisy transcriptomic background.

## Candidate Qualifying & QC Metrics

- **Genotype & Calling Quality**:
  - GATK hard-filtered heterozygous calls (`0/1`): $\text{FS} \le 30.0$, $\text{QD} \ge 2.0$, $\text{GQ} \ge 20$, $\text{DP} \ge 10$.
  - Minimum of $\ge 2$ supporting reads per allele in unbiased Input baseline.
- **Cross-Subject Recurrence**:
  - Minimum coverage across $\ge 8$ individuals with balanced case vs. control representation.
- **Allelic Read Depth in Target Fraction**:
  - Robust IP read support (hundreds to thousands of allelic counts per locus) to prevent low-depth binomial noise.
- **Reference-Mapping Bias Control**:
  - Personal Input allelic fraction offset (`logit(Input AF)`) applied to absorb reference alignment bias.
- **Statistical Significance & Threshold Stability**:
  - Stringent multiple testing threshold: $\text{FDR} < 0.05$ (Benjamini-Hochberg).
  - High effect size filter ($|\Delta\text{AF}| \ge 0.20$).


<!-- pagebreak -->

# Conclusion & Outlook

- **Takeaway**:
  - Demonstrated that standard DMRs miss allele-level regulatory events.
  - Proof-of-concept pipeline successfully maps personal genetic variants to sex-dimorphic epitranscriptomic dysregulation in MDD vmPFC.
- **Proposed Next Steps**:
  - More samples for statistical integrity
  - Prioritize top teased candidates for wet-lab validation (targeted MeRIP-qPCR / Nanopore direct RNA-seq).
  - In silico disruption modeling (AlphaGenome / RNA secondary structure folding).
  - Reporter minigene assays in human astrocyte/neuronal cell lines.
