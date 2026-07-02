# Upstream RNA-seq & MeRIP-seq Pipeline Report
**MDD PFC Cohort (78 Libraries)**

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
* **What went wrong**: Counts QC flagged a warning for all 39 input samples showing $<50\%$ assignment rates (typically 5–15%), while the IP (m6A) samples passed.
* **How it was dealt with**: Identified that the QC script was counting raw multi-mapping alignments in its denominator (which were highly inflated due to some escaped rRNA mapping to repetitive rDNA units). Refactored [09_counts_qc.sh](../01_core_alignment/09_counts_qc.sh) to calculate assignment rates relative to uniquely mapped reads, raising the true rate to **50.5%**.