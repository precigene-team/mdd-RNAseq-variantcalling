# Project Context: Epitranscriptomics & Transcriptomics of Major Depressive Disorder (MDD)

---

## 🧠 1. Biological and Clinical Background

Major Depressive Disorder (MDD) is a severe, multi-faceted neuropsychiatric condition with high lifetime prevalence. While extensive transcriptomic analyses have identified gene expression changes in postmortem human brains, the post-transcriptional mechanisms regulating gene expression—specifically **RNA modifications (epitranscriptomics)**—remain largely uncharacterized in psychiatric disorders.

### The Role of N6-Methyladenosine ($m^6A$)
* **Definition**: $m^6A$ is the most prevalent, dynamic internal modification on eukaryotic messenger RNA (mRNA).
* **Function**: It governs RNA splicing, stability, localization, and translation efficiency. In neurons, $m^6A$ regulates synaptic plasticity, local translation in dendrites, and behavioral responses to stress.
* **Regulation**: The methylation landscape is balanced by:
  * **Writers** (methyltransferases): `METTL3`, `METTL14`, `WTAP`.
  * **Erasers** (demethylases): `FTO`, `ALKBH5`.
  * **Readers** (RNA-binding proteins): `YTHDF1-3`, `YTHDC1-2`.

### Sex-Specific Dysregulation
MDD shows prominent sex differences in prevalence (females are twice as likely to be diagnosed), symptomatology, and pharmacological response. Recent evidence shows that transcriptomic and $m^6A$ epitranscriptomic landscapes in the human **Prefrontal Cortex (PFC)** are regulated in a highly sex-dependent manner, with distinct molecular networks showing opposite alterations in males vs. females.

---

## 📊 2. Cohort Metadata & Design

The study is conducted on a cohort of **39 human postmortem brain subjects** (PFC tissue). The experimental design employs paired-end **MeRIP-seq (m6A-seq)**, generating two libraries per subject:
1. **Input Fraction (`RNA-seq`)**: Total sheared RNA, representing baseline transcriptional abundance.
2. **IP Fraction (`m6A`)**: Methylated RNA captured using an anti-$m^6A$ antibody, representing the abundance of methylated transcripts.

### Subject Stratification Table
| Cohort Group | Number of Subjects | Total Libraries (Input + IP) | Accession Range |
|---|---|---|---|
| **Male MDD** | 10 | 20 | [`SRR30401480`](file:///mnt/22T/phongdinh/RNAseq_pipeline_depression/_data/_metadata/sample_sheet.csv) - [`SRR30401481`](file:///mnt/22T/phongdinh/RNAseq_pipeline_depression/_data/_metadata/sample_sheet.csv) & [`SRR30401472`](file:///mnt/22T/phongdinh/RNAseq_pipeline_depression/_data/_metadata/sample_sheet.csv) - [`SRR30401480`](file:///mnt/22T/phongdinh/RNAseq_pipeline_depression/_data/_metadata/sample_sheet.csv) |
| **Male Control** | 10 | 20 | [`SRR30401482`](file:///mnt/22T/phongdinh/RNAseq_pipeline_depression/_data/_metadata/sample_sheet.csv) - [`SRR30401491`](file:///mnt/22T/phongdinh/RNAseq_pipeline_depression/_data/_metadata/sample_sheet.csv) |
| **Female MDD** | 10 | 20 | [`SRR30401505`](file:///mnt/22T/phongdinh/RNAseq_pipeline_depression/_data/_metadata/sample_sheet.csv) - [`SRR30401513`](file:///mnt/22T/phongdinh/RNAseq_pipeline_depression/_data/_metadata/sample_sheet.csv) & [`SRR30401492`](file:///mnt/22T/phongdinh/RNAseq_pipeline_depression/_data/_metadata/sample_sheet.csv) |
| **Female Control** | 9 | 18 | [`SRR30401514`](file:///mnt/22T/phongdinh/RNAseq_pipeline_depression/_data/_metadata/sample_sheet.csv) - [`SRR30401522`](file:///mnt/22T/phongdinh/RNAseq_pipeline_depression/_data/_metadata/sample_sheet.csv) |

---

## 📈 3. Differential Expression & Methylation Contrast Design

To resolve the transcriptomic and epitranscriptomic mechanisms of MDD, the analysis implements three types of contrast models:

### A. Total RNA-seq Expression Changes (Input)
Isolates alterations in steady-state transcript levels:
$$\text{Contrast}_{\text{Male}} = \text{Male MDD}_{\text{Input}} \text{ vs } \text{Male Ctrl}_{\text{Input}}$$
$$\text{Contrast}_{\text{Female}} = \text{Female MDD}_{\text{Input}} \text{ vs } \text{Female Ctrl}_{\text{Input}}$$
* **Fold Change Prioritization**: Target genes with $|\text{log2FC}| \ge 2$ to capture high-effect transcriptional alterations.

### B. m6A IP Expression (IP only)
Measures the absolute abundance of methylated transcript copies:
$$\text{Contrast}_{\text{Male IP}} = \text{Male MDD}_{\text{IP}} \text{ vs } \text{Male Ctrl}_{\text{IP}}$$
$$\text{Contrast}_{\text{Female IP}} = \text{Female MDD}_{\text{IP}} \text{ vs } \text{Female Ctrl}_{\text{IP}}$$

### C. Differential Methylation (IP vs Input Interaction)
A gene's IP abundance can change because of:
1. Changes in the relative methylation level (the fraction of molecules methylated).
2. Changes in the total transcript abundance (underlying gene transcription).

To isolate true relative methylation changes (independent of transcription levels), we analyze the **interaction contrast**:
$$\text{Design} = \sim \text{Condition} + \text{Fraction} + \text{Condition} \times \text{Fraction}$$
The interaction term evaluates:
$$\Delta m^6A = \left(\frac{\text{IP}}{\text{Input}}\right)_{\text{MDD}} \div \left(\frac{\text{IP}}{\text{Input}}\right)_{\text{Ctrl}}$$
This isolates actual post-transcriptional $m^6A$ methylation differences.

---

## 🧬 4. Variant Calling Goals

By executing RNA-seq based variant calling (using GATK Best Practices) across both `RNA-seq` (Input) and `m6A` (IP) fractions, the project aims to identify:
1. **m6A-SNPs (m6A-modifying variants)**: Genomic single nucleotide polymorphisms (SNPs) located within or near consensus $m^6A$ motifs (`DRACH`) that disrupt or create methylation sites.
2. **Allele-Specific Methylation (ASM)**: Overlap of variants in MeRIP-seq IP fractions showing biased allele ratios compared to Input fractions, suggesting that one allele is preferentially methylated.
