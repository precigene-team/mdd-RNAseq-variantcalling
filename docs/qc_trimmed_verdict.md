# Trimmed Reads Bioinformatics QC Verdict Report

**Project:** MDD RNA-seq & MeRIP-seq Analysis  
**Cohort Size:** 39 subjects (78 libraries: 39 Input fraction, 39 IP fraction, paired-end)  
**Date:** June 26, 2026  
**Status:** **PROCEEDABLE TO ALIGNMENT**

---

## 📋 1. Executive Summary & Verdict

We have completed the trimmed quality control assessment of the reads in the `_data/trimmed_fastq` directory using MultiQC (results located in `scratch/multitrim/`).

*   **Overall Sequence Quality:** **Outstanding**. The quality filtering trimmed away poor bases and adapters while preserving the integrity and depth of the libraries. Average Q30 rate is **93.66%**.
*   **Sequencing Depth:** **Excellent**. Average read depth is **133.43M** reads per library (R1 + R2). The minimum depth sample (`female_ctrl_6_input`) contains **45.11M** reads, which is more than double the standard 20M QC threshold.
*   **Read Survival Rate:** **Excellent**. An average of **98.36%** of raw reads survived filtering. No libraries suffered from excessive read loss (minimum survival is **96.22%**).
*   **Adapter Removal:** **Highly Effective**. An average of **12.07%** of reads contained adapter sequences that were successfully trimmed. 
*   **Verdict:** **Proceed directly to Step 04 (STAR Alignment)**. The trimmed libraries are exceptionally clean and fully optimized for alignment.
*   **Analysis Reference**: Compiled via the analysis script [summarize_trim_qc.py](../scratch/summarize_trim_qc.py). (ON HPC)

---

## 📊 2. Detailed QC Metrics Summary

### Section A: Sequencing Depth & Filtering Results (in Millions of Reads)
*   **Total Libraries:** 78 libraries
*   **Mean Depth:** 133.43M reads per library (combined R1/R2)
*   **Standard Deviation:** 48.93M reads
*   **Minimum Depth:** 45.11M reads (`female_ctrl_6_input`)
*   **Maximum Depth:** 249.32M reads (`female_ctrl_2_m6A`)
*   **Verdict:** **Pass**. No libraries suffer from poor sequencing depth or dropout after filtering.

### Section B: Quality and Base Composition
*   **Q30 Rate (%):** Mean = **93.66%** (Min = 91.60% `female_ctrl_3_input` / Max = 94.68% `male_mdd_4_m6A`)
*   **Read Survival Rate (%):** Mean = **98.36%** (Min = 96.22% `male_mdd_4_input` / Max = 99.48% `female_mdd_4_m6A`)
*   **Adapter Content (%):** Mean = **12.07%** (Min = 0.88% `female_mdd_2_m6A` / Max = 35.54% `male_mdd_4_input`)
*   **Duplication Rate (%):** Mean = **28.80%** (Min = 7.28% `male_mdd_9_input` / Max = 49.62% `female_mdd_8_m6A`)
*   **GC Content (%):** Mean = **55.17%** (Min = 50.16% `male_mdd_9_input` / Max = 60.32% `male_mdd_4_input`)

---

## 🚨 3. Quality Verification Checks

*   **Read Depth Check:** ✅ All libraries have $\ge$ 20M passed reads (minimum is 45.11M).
*   **Survival Rate Check:** ✅ All libraries have $\ge$ 90% read survival rate (minimum is 96.22%).
*   **GC Distribution:** ✅ Stable GC composition (50–60%) across the cohort, matching human transcriptome standards.