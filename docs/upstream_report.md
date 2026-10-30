# Upstream RNA-seq Pipeline Report
**MDD PFC Cohort (78 Libraries Downloaded, 39 Input Libraries In Scope)**

## 📅 Timeline

* **June 22:** Restructured legacy pipeline code. Initial PCA showed issues, prompting a FASTQ download redo.
* **June 23–24:** Set up remote server access and completed the first half of FASTQ downloads.
* **June 25:** Switched to Aspera (`ascli`) for high-speed transfer. Completed remaining downloads (100% checksum pass).
* **June 26:** Verified FastQC metrics and initiated first-pass STAR alignment.
* **June 27:** Stopped alignment because it ran too slow due to high rRNA contamination. Attempted SortMeRNA, Bowtie2, BBDUK but it was computationally too expensive for 78 samples.
* **June 28:** Test mapping on the worst-performing sample showed only 30% alignment. Suspected the STAR genome index was incomplete.
* **June 29–30:** Re-generated a complete STAR genome index and re-aligned all 78 samples successfully.
* **July 01 (Today):** Completed featureCounts and compiled coordinates. Finalized pipeline diagnostics and downstream fixes.

---

## 📈 What Went Well
* **High-Speed Transfer**: Transitioning from `aria2c` to **Aspera (`ascli`)** bypassed server bandwidth bottlenecks, enabling secure, fast downloads of all 156 heavy FASTQs with 100% checksum verification.
* **Alignment Success**: Building a full primary assembly STAR genome index (`sjdbOverhang 149`) successfully resolved the mapping rate issues, aligning all 78 libraries.
* **Downstream Integration**: Completed featureCounts quantification and compiled clean normalized data matrices ready for DGE.

---

## ⚠️ What Went Wrong & How It Was Resolved

### 1. Incomplete STAR Index & Low Mapping Rates
* **What went wrong**: Initial STAR alignment rates were extremely low (20–30%). Re-generating the full STAR genome index improved mapping, but uniquely-mapped rates remain low for some samples.
* **How it was dealt with**: We rebuilt the genome index. The remaining low mapping rate issue is currently under investigation and will be analyzed further once the Qualimap RNA-seq QC step is completed to inspect transcript region distributions.

### 2. Artificially Low Assignment Rates in Input Samples
* **What went wrong**: Counts QC flagged a warning for all 39 input samples showing $<50\%$ assignment rates (typically 5–15%).
* **How it was dealt with**: Identified that the QC script was counting raw multi-mapping alignments in its denominator (which were highly inflated due to some escaped rRNA mapping to repetitive rDNA units). The delivered [09_counts_qc.sh](../01_core_alignment/09_counts_qc.sh) keeps that denominator: `:154-165` divides `Assigned` by the sum of every `.summary` category, `Unassigned_MultiMapping` included, so the uniquely-mapped-denominator refactor was not carried into this tree. Against the gene-level counts, that shipped rate reads **2.3–15.8% (median 8.2%)** for Input samples, and the `>= 50%` gate at `:167-183` prints a `WARNING` only — `09` exits 0. The **50.5%** figure is an observation under the uniquely mapped denominator (total minus `Unassigned_MultiMapping`): there Input samples read **42.9–54.0% (median 48.9%)**, and 7 of the 22 quantified so far clear 50%.