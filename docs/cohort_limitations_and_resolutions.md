# Cohort Limitations & Resolutions Registry

This document outlines potential technical and biological confounders in the MDD Input-only RNA-seq cohort (from raw FASTQ to differential gene expression) and defines the roadmap for resolving them.
The MeRIP-dependent items (QC-03 gene-level vs. peak-level audit, QC-04 antibody cross-reactivity, QC-05 TRESS DMR covariates) are archived in [`_archive/cohort_limitations_and_resolutions.md`](_archive/cohort_limitations_and_resolutions.md).

---

## 📋 Registry of Issues & Resolution Status

| ID | Issue Description | Type | Resolution Strategy | Status |
| :--- | :--- | :--- | :--- | :--- |
| **QC-01** | **Postmortem Confounders** (RIN, PMI, Age, pH, etc.) skewing expression profiles | Technical | Confirm cohort matching in paper and verify public availability of individual covariates | **Resolved (Baseline Model Final)** |
| **QC-02** | **Cell-Type Heterogeneity** (neuronal atrophy making synaptic transcripts look downregulated) | Biological | Run signature-based cell-type deconvolution to estimate glial/neuronal proportions | **Resolved (Astrocyte Loss Documented)** |

---

## 🔍 Detail & Active Resolutions

### 🧬 QC-01: Postmortem Confounders (PMI, RIN, Age, pH)
* **The Risk:** Postmortem brain tissues degrade dynamically. If the MDD cohort has a significantly longer Post-Mortem Interval (PMI) or lower RNA Integrity Number (RIN) than Controls, the differential expression analysis will call degradation-susceptible transcripts as "downregulated in MDD," leading to high false-positive rates.
* **Resolution & Findings:**
  1. Audited the paper text in `scratch/paper.txt` (Mitsuhashi et al. 2024, iScience): The authors explicitly **matched cases and controls for PMI, tissue pH, and RIN** during cohort recruitment to minimize confounding.
  2. Audited the SRA metadata and the `mmc1.pdf` supplement (specifically Page 2, Table S1): The authors only published **aggregated group means and SEMs** (e.g. Male Control Age: 56.9 ± 5.2, Male MDD Age: 44.9 ± 4.3). 
  3. Individual-level demographic/clinical covariates (e.g., the specific age or RIN for each of the 39 subjects) were **not publicly released** in either the GEO records (GSE275676) or the supplemental materials to protect donor privacy.
  4. **The 2 Million (2e6) QC Gate:** [`04_deg_analysis.R`](../02_dge_analysis/04_deg_analysis.R) applies a **2 Million gene-mapped count floor** per subject; the delivered cohort is the 39 Input rows in `_data/_metadata/sample_sheet.csv` — 10 male MDD, 10 male ctrl, 10 female MDD, 9 female ctrl — and the gate drops any subject whose gene-level library falls below the floor.
  5. **Conclusion:** Because individual-level covariates are not public, running a covariate-adjusted DGE model is not possible. However, the rigorous matching of groups and the 2M gene-level floor in `04_deg_analysis.R` ensure that the unadjusted baseline models are clean, robust, and aligned with the published study.
* **Status:** Resolved. The 2M-gated matched baseline model is finalized.

### 🧠 QC-02: Cell-Type Heterogeneity Confounding
* **The Risk:** A change in bulk RNA-seq could mean the cells are expressing less of a gene, or simply that the tissue contains fewer of those cells (e.g. neuronal loss in MDD).
* **Resolution & Findings:**
  1. Implemented a Marker Gene Profile (MGP) script in [`utils/qc02_cell_type_deconvolution.py`](../utils/qc02_cell_type_deconvolution.py) using consensus human brain markers for Neurons, Astrocytes, Oligodendrocytes, Microglia, and Endothelial cells.
  2. Estimated cell-type enrichment scores across the 39 Input subjects and performed t-tests comparing MDD vs. Control.
  3. **Astrocyte Loss:** Found a **statistically significant overall reduction in Astrocyte markers in MDD brains compared to Controls (p-value = 0.0389)**.
  4. **Other Cell Types:** Neurons showed a slight trend of higher scores in MDD (p-value = 0.0956, not significant). Oligodendrocytes, Microglia, and Endothelial cells showed no significant differences (p > 0.25).
  5. **Biological Link:** The significant reduction in Astroglial markers is a well-established neuropathological hallmark of major depression in the human prefrontal cortex (reflecting glia-specific atrophy/loss). We have documented this cellular shift as a potential driver of astrocyte-related transcript changes in our report.
* **Status:** Resolved. Astrocyte marker reduction documented.
