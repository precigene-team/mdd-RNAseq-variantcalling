# Understanding RNA-seq vs. MeRIP-seq (m6A-seq)

This document provides the theoretical background explaining the differences between standard **Bulk RNA-seq** (Input fraction) and **MeRIP-seq** (Immunoprecipitation fraction), and how they are used together to analyze RNA methylation dynamics in Major Depressive Disorder (MDD).

---

## 🔬 1. Conceptual Comparison

| Feature | Bulk RNA-seq (Input Control) | MeRIP-seq / m6A-seq (IP Fraction) |
| :--- | :--- | :--- |
| **Primary Target** | All transcript molecules in the sample. | Only transcript fragments containing $m^6A$ modifications. |
| **Selection Method** | Poly(A) enrichment or ribosomal RNA depletion. | Immunoprecipitation using an $m^6A$-specific antibody. |
| **Output Information** | Total gene/transcript expression levels. | Absolute abundance of methylated transcripts. |
| **Biological Insight** | Transcriptional activation and RNA degradation rate. | Post-transcriptional chemical modification of RNA. |
| **Typical Analysis** | Differential Gene Expression (DGE). | Peak Calling, Differential Methylation (interaction analysis). |

---

## ⚙️ 2. How MeRIP-seq Works

MeRIP-seq (Methylated RNA Immunoprecipitation Sequencing) is an antibody-based technique used to map $m^6A$ modifications across the transcriptome:

```mermaid
graph TD
    A[Total RNA Extraction] --> B[RNA Fragmentation ~100-200 nt]
    B --> C[Divide into two fractions]
    C --> D[Input Fraction: Standard Library Prep]
    C --> E[IP Fraction: Anti-m6A Immunoprecipitation]
    E --> F[Elute Methylated RNA fragments]
    F --> G[IP Library Prep]
    D --> H[Sequencing: Total Expression Baseline]
    G --> I[Sequencing: Methylated RNA Sites]
```

1. **Fragmentation**: RNA is fragmented into short segments (~100–200 nucleotides) before selection. This is critical because the antibody pulls down the entire fragment containing the modification. Fragmenting allows high-resolution localization of the modification site (the "peak").
2. **Immunoprecipitation (IP)**: An antibody specific to $m^6A$ is incubated with the fragments. The antibody-bound fragments are precipitated, and unmethylated fragments are washed away.
3. **Sequencing**: Both the **Input** (unprecipitated total RNA) and **IP** (precipitated methylated RNA) libraries are sequenced.

---

## 🧮 3. Why We Need Both: The Interaction Logic

Analyzing the IP fraction alone is **insufficient** to determine if RNA methylation has changed. Consider the following scenarios:

### Scenario A: Transcriptional Upregulation (No Methylation Change)
* **Control**: 100 transcripts of Gene X. 10 are methylated (10% methylation rate).
  * *Input count* $\propto$ 100, *IP count* $\propto$ 10.
* **MDD**: 200 transcripts of Gene X due to gene activation. 20 are methylated (10% methylation rate).
  * *Input count* $\propto$ 200, *IP count* $\propto$ 20.
* **Result**: The IP signal doubled in MDD. However, the methylation rate remained exactly 10%. This is transcriptional upregulation, not differential methylation.

### Scenario B: Hyper-methylation (True Differential Methylation)
* **Control**: 100 transcripts of Gene X. 10 are methylated (10% methylation rate).
  * *Input count* $\propto$ 100, *IP count* $\propto$ 10.
* **MDD**: 100 transcripts of Gene X. 40 are methylated (40% methylation rate).
  * *Input count* $\propto$ 100, *IP count* $\propto$ 40.
* **Result**: The IP signal quadrupled in MDD while the Input signal remained unchanged. This represents true hyper-methylation of Gene X in MDD.

### Mathematical Normalization via DESeq2
To isolate true methylation changes from changes in gene expression, we test the **interaction term** in a generalized linear model:
$$\text{Design} = \sim \text{Condition} + \text{Fraction} + \text{Condition:Fraction}$$
Where:
* `Condition`: `mdd` vs. `ctrl`
* `Fraction`: `IP` vs. `Input`
* `Condition:Fraction`: The interaction term.

A statistically significant interaction term indicates that the ratio of $\frac{\text{IP}}{\text{Input}}$ differs between MDD and Control, confirming **differential methylation**.

---

## 🧬 4. Variant Calling in Input vs. IP

Performing SNP variant calling on standard RNA-seq vs. MeRIP-seq IP libraries serves two distinct purposes:
1. **Input RNA-seq**: Identifies genomic variants (SNPs) expressed in the subject's brain tissue.
2. **IP Fraction**: Identifies variants present in the methylated pool of transcripts. 
   * By comparing allele ratios in the **IP** vs. **Input** libraries, researchers can detect **Allele-Specific Methylation (ASM)**—where one allele of a heterozygous SNP is preferentially methylated over the other.
