# Raw Reads Bioinformatics QC Verdict Report

**Project:** MDD RNA-seq & MeRIP-seq Analysis  
**Cohort Size:** 39 subjects (78 libraries: 39 Input fraction, 39 IP fraction, paired-end, total of 156 FastQC files)  
**Date:** June 26, 2026  
**Status:** **PROCEEDABLE WITH TRIMMING**

---

## 📋 1. Executive Summary & Verdict

We have completed the raw quality control assessment of the sequencing reads in the `_data/qc_raw` directory using MultiQC (per-library FastQC reports in `_data/qc_raw/`, aggregated by [`10_multiqc.sh`](../../01_core_alignment/10_multiqc.sh) into `_data/multiqc/`).

*   **Overall Sequence Quality:** **Excellent**. All libraries show high-confidence base calling with perfect Phred score profiles across the entire read length (101 bp).
*   **Sequencing Depth:** **Excellent**. Average read depth is **67.92M** reads per library. All 78 libraries significantly exceed the typical minimum 10M read depth threshold (minimum sample has **23.23M** reads).
*   **Adapter Contamination:** **High (expected)**. Out of 156 files, **52 failed** and **26 warned** on adapter content checks. 
*   **Verdict:** **Proceed directly to Step 02 (Trimming)**. Raw data is highly robust, but adapter trimming is strictly necessary to clean contamination before mapping.
*   **Analysis Reference:** Detailed stats compiled from the FastQC reports in `_data/qc_raw/`, produced on the HPC by [`01_raw_fastqc.sh`](../../01_core_alignment/01_raw_fastqc.sh).

---

## 📊 2. Detailed QC Metrics Summary

### Section A: Sequencing Depth (in Millions of Reads)
*   **Total FastQC Reports:** 156 files
*   **Mean Depth:** 67.92M reads per file
*   **Standard Deviation:** 25.06M reads
*   **Minimum Depth:** 23.23M reads (`female_ctrl_6_input_1`)
*   **Maximum Depth:** 125.63M reads (`female_ctrl_2_m6A_1`)
*   **Verdict:** **Pass**. No libraries suffer from poor sequencing depth or dropout.

### Section B: Quality and Base Composition
*   **Per Base Sequence Quality:** **100% Pass (156/156)**. Median quality score is consistently $>30$ (Phred) across all bases.
*   **Per Sequence Quality Scores:** **100% Pass (156/156)**. High-quality average Phred distribution.
*   **Basic Statistics:** **100% Pass (156/156)**. Clear documentation of conventional Sanger/Illumina 1.9 encoding.

### Section C: Expected Biological/Technical Warning and Failures
Standard RNA-seq and MeRIP-seq libraries typically trigger warnings or failures in several FastQC categories. These are not quality issues:

1.  **Sequence Duplication Levels (141 Fails, 14 Warns, 1 Pass):**
    *   *Cause:* PCR duplication of low-input cDNA in the IP (m6A) fractions and extremely high expression of specific transcripts in postmortem brain tissue.
2.  **Per Base Sequence Content (153 Fails, 3 Warns):**
    *   *Cause:* Caused by sequence composition bias during random hexamer priming (first 10–15 bp) and biological enrichment of m6A consensus motifs (`DRACH`) in the IP fraction.
3.  **Per Sequence GC Content (35 Fails, 113 Warns, 8 Passes):**
    *   *Cause:* Typical transcriptome-wide GC skew and peak distribution deviations due to highly expressed housekeeping or mitochondrial genes.
4.  **Overrepresented Sequences (6 Fails, 141 Warns, 9 Passes):**
    *   *Cause:* Highly expressed specific transcripts or adapter sequences.

---

## 🚨 3. Key Issues & Remediation Action Plan

### Issue: Adapter Contamination
*   **Metric:** 52 files failed and 26 files warned.
*   **Observation:** Significant adapter read-through exists in 50% of the files (particularly reverse reads).
*   **Remediation:** Trimming must be executed using Trimmomatic with adapter databases specified in the pipeline configurations to clean reads.

### Verification of Outliers:
*   No outliers detected below 10M total raw reads. All samples show superb starting read depths.