# Allele Skew and Condition Shift Calculation in ASm6A Analysis

> **SUPERSEDED — read-ratio ASM, retained as historical record.**
> The allele-skew and condition-shift definitions below belong to the retired
> window-level pipeline. The current module (`04_asm/01`–`07`) uses a per-donor
> IP-vs-Input design in which each donor is self-normalising, so the cross-group
> condition-shift contrast is no longer the primary statistic. See
> `docs/04_asm_methods.md` and `docs/04_asm_postmortem.md`.

This document outlines the mathematical definitions, biological significance, and step-by-step calculations for **Allele Skew** and **Condition Shift** in the Allele-Specific m6A Methylation (ASm6A) pipeline.

---

## 1. Overview of Key Metrics

| Metric | Comparison | Question Answered | Formula (Plain Text) |
| :--- | :--- | :--- | :--- |
| **Allele Fraction (AF)** | Allele Reads vs. Total | What proportion of reads carries the reference/alternate allele? | `AF = Allele_Reads / Total_Reads` |
| **Single-Sample Skew** | IP Pool vs. 50% Null | Does methylation favor one allele over another in an individual? | `Skew = |IP_AF - 0.50|` |
| **Baseline Allele Skew** | IP Pool vs. Input Baseline | How much does the methylated pool deviate from total expression? | `Delta_Skew = IP_AF - Input_AF` |
| **Empirical Condition Shift**| MDD IP vs. Control IP | How much does allelic methylation fraction shift in MDD? | `MDD_Shift_% = (Mean_IP_AF_MDD - Mean_IP_AF_Ctrl) * 100` |
| **Model Shift (Log-Odds)** | MDD vs. Control (GLMM) | What is the statistical effect size of disease on allele odds? | `Beta_MDD from glmmTMB(cbind(Ref, Alt) ~ Condition + Sex + offset(logit(Input_AF)))` |

---

## 2. Mathematical Formulations

### A. Allele Fraction (AF)
For any given locus and sequencing pool (Input or IP):
```text
AF = Allele_Reads / Total_Reads
```
* In this pipeline:
  * `in_af = in_ref / in_total` (or smoothed with pseudo-counts: `(in_ref + 0.5) / (in_total + 1.0)`)
  * `ip_af = ip_ref / ip_total`

---

### B. Allele Skew (Input vs. IP Methylated Pool)

1. **Deviation from Neutral 50/50 (`12_asm_compare.R` filtering)**:
   ```text
   Skew = |ip_af - 0.50|
   ```
   * **Filter Threshold (`SKEW_CUT = 0.20`)**: Requires at least a **70:30** allelic preference (`|0.70 - 0.50| = 0.20`) to eliminate minor stochastic noise.

2. **Baseline-Corrected Allele Skew (`15_plot_sig_asm_drach.R` / Summary Tables)**:
   ```text
   Delta_Skew = ip_af - in_af
   Abs_Skew_% = |ip_af - in_af| * 100%
   ```
   * Controlling against `in_af` removes reference-mapping bias (e.g., aligner preference for REF).

---

### C. Condition Shift (Control IP vs. MDD IP)

1. **Empirical Percentage Shift**:
   ```text
   MDD_Shift_% = (Mean_IP_AF_MDD - Mean_IP_AF_Ctrl) * 100%
   ```
   * Positive value: MDD increases preference for the target allele.
   * Negative value: MDD decreases preference for the target allele.

2. **Statistical Effect Size via Beta-Binomial GLMM (`glmmTMB`)**:
   ```text
   logit(P(IP_Allele)) = Beta_0 + Beta_MDD * Condition + Beta_Sex * Sex + offset(logit(in_af))
   ```
   * **Shift Odds Ratio (OR)**:
     ```text
     OR = exp(Beta_MDD)
     ```
   * `Beta_MDD`: The log-odds shift in allele preference in MDD patients relative to controls, adjusted for sex and baseline Input RNA balance.

---

## 3. Concrete Example: MGAT5 Discovery Locus

Coordinates: `chr2:134448632-134449882` (from [12_asm_primary_significant.csv](../results/03_variant_calling/tables/12_asm_primary_significant.csv))

### Raw Cohort Observations
* **Total subjects covered**: N = 11 (6 MDD, 5 Controls)
* **Input Baseline (`in_af`)**: `50.0%` (0.500)
* **Control IP Average (`ip_af_ctrl`)**: `85.9%` (0.859)
* **MDD IP Average (`ip_af_mdd`)**: `92.7%` (0.927)

---

### Step 1: Calculate Allele Skew
* **Control Baseline Skew**:
  ```text
  Skew_Ctrl = 85.9% - 50.0% = +35.9%
  ```
* **Biological Meaning**: In healthy controls, m6A methylation is strongly skewed toward the reference allele (35.9% higher than neutral expectation).

---

### Step 2: Calculate Condition Shift
* **Empirical MDD Shift**:
  ```text
  MDD_Shift_% = 92.7% - 85.9% = +6.8%
  ```
* **Statistical GLMM Parameters**:
  * `Beta_MDD` = `+0.966`
  * `Standard Error (SE)` = `0.242`
  * `p-value` = `6.52e-05`
  * `BH-FDR` = `0.0176` (FDR < 0.05)
  * `Odds Ratio` = `exp(0.966) = 2.63`

---

### Summary Interpretation
1. **Allele Skew (+35.9% in Controls)** confirms that *MGAT5* exhibits strong **Allele-Specific Methylation (ASm6A)** where the REF allele is preferentially methylated over the ALT allele.
2. **Condition Shift (+6.8% shift, OR = 2.63 in MDD)** demonstrates a statistically significant **hyper-skew in Major Depressive Disorder**, meaning MDD patients have 2.63-fold higher odds of depositing m6A on the preferred allele compared to matched controls.
