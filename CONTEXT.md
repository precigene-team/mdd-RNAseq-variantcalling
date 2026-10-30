# Project Context: Transcriptomics of Major Depressive Disorder (MDD)

---

## 🧠 1. Biological and Clinical Background

Major Depressive Disorder (MDD) is a severe, multi-faceted neuropsychiatric condition with high lifetime prevalence, and there is a large unmet need for objective biological markers. This study looks for measurable signals that separate MDD patients from healthy controls: transcriptional differences measured by RNA-seq, and germline variant differences measured by RNA-seq-based SNP calling.

The analyses are run on postmortem **prefrontal cortex (PFC)** tissue and stratified by sex, because MDD prevalence, symptomatology, and pharmacological response differ markedly between males and females (females are diagnosed at roughly twice the rate).

---

## 📊 2. Cohort Metadata & Design

The study is conducted on a cohort of **39 human postmortem brain subjects** (PFC tissue). The delivered analysis uses the paired-end **RNA-seq (Input)** libraries — 39 libraries, one per subject:
1. **Input Fraction (`RNA-seq`)**: Total sheared RNA, representing baseline transcriptional abundance.

### Subject Stratification Table
| Cohort Group | Number of Subjects | Input Libraries | Accession Range |
|---|---|---|---|
| **Male MDD** | 10 | 10 | [`SRR30401480`](_data/_metadata/sample_sheet.csv) - [`SRR30401481`](_data/_metadata/sample_sheet.csv) & [`SRR30401472`](_data/_metadata/sample_sheet.csv) - [`SRR30401480`](_data/_metadata/sample_sheet.csv) |
| **Male Control** | 10 | 10 | [`SRR30401482`](_data/_metadata/sample_sheet.csv) - [`SRR30401491`](_data/_metadata/sample_sheet.csv) |
| **Female MDD** | 10 | 10 | [`SRR30401505`](_data/_metadata/sample_sheet.csv) - [`SRR30401513`](_data/_metadata/sample_sheet.csv) & [`SRR30401492`](_data/_metadata/sample_sheet.csv) |
| **Female Control** | 9 | 9 | [`SRR30401514`](_data/_metadata/sample_sheet.csv) - [`SRR30401522`](_data/_metadata/sample_sheet.csv) |

---

## 📈 3. Differential Expression Contrast Design

To resolve the transcriptomic mechanisms of MDD, the analysis implements sex-stratified contrasts on the Input libraries:

### A. Total RNA-seq Expression Changes (Input)
Isolates alterations in steady-state transcript levels:
$$\text{Contrast}_{\text{Male}} = \text{Male MDD}_{\text{Input}} \text{ vs } \text{Male Ctrl}_{\text{Input}}$$
$$\text{Contrast}_{\text{Female}} = \text{Female MDD}_{\text{Input}} \text{ vs } \text{Female Ctrl}_{\text{Input}}$$
* **Fold Change Prioritization**: Target genes with $|\text{log2FC}| \ge 2$ to capture high-effect transcriptional alterations.

---

## 🧬 4. Variant Calling Architecture

The pipeline executes RNA-seq variant discovery and clinical interpretation in one dedicated module:
1. **Variant Calling & Clinical Interpretation (`03_variant_calling/`)**: GATK joint calling on RNA-seq inputs followed by a 3-stage clinical translation pipeline:
   - **Stage 1 (Functional)**: In silico deleteriousness (CADD, SpliceAI) and psychiatric brain domain localization (PsychENCODE, PGC MDD GWAS).
   - **Stage 2 (Pathogenicity)**: ACMG/AMP 5-tier classification and ClinVar cross-referencing.
   - **Stage 3 (Pharmacogenomics)**: CPIC clinical practice guidelines and ClinPGx antidepressant actionability (metabolism & target sensitivity).
