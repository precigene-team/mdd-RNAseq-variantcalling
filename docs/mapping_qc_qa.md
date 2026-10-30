# RNA-seq Mapping & Quality Control Q&A Summary

This document summarizes the technical discussion and decisions made regarding upstream mapping issues, library size gating, and statistical power in the MDD Input RNA-seq cohort.

---

### Q1: Does low count mapping upstream affect the downstream processes?
**Answer:** Yes, absolutely. Low mapping rates or low read counts upstream have a cascading negative effect:
1. **DGE (DESeq2):** Low counts mean high random sampling variation. DESeq2 estimates high dispersions for these genes, reducing statistical power and failing to call genes as significant. Furthermore, DESeq2's independent filtering will discard low-expression genes before testing.
2. **Variant Calling (GATK):** `HaplotypeCaller` requires sufficient read coverage at a specific position to call a variant. Low counts lead to high false-negative rates.

---

### Q2: What are the worst samples in our cohort?
**Answer:** The lowest count samples are Input samples:
* **`female_mdd_2_input`** — **5.8 Million** counts (extremely low)
* **`female_ctrl_6_input`** — **9.9 Million** counts
* **`male_ctrl_1_input`** — **11.6 Million** counts
* **`male_mdd_2_input`** — **12.1 Million** counts
* **`female_mdd_8_input`** — **12.2 Million** counts

---

### Q3: If STAR has a low mapping rate, does it affect counting? What are the standard thresholds for bulk RNA-seq and DEGs?
**Answer:** Yes. `featureCounts` cannot count what STAR fails to align. If STAR has a low unique mapping rate, the exonic counts will be low.

Standard quality thresholds are:
* **Total Read Depth:** 20M–30M+ reads (Good), < 10M (Red Flag).
* **STAR Unique Alignment:** > 75%–85% (Good), < 50% (Red Flag).
* **featureCounts Exonic Assignment:** > 60% (Good), < 30% (Red Flag).
* **DEG Filtering:** Nominal $p < 0.05$ (standard for noisy psychiatric cohorts) and $|\log_2\text{FC}| \ge 0.2$ (equivalent to a $\approx 15\%$ expression shift in complex brain tissues).

---

### Q4: Why do Input samples have such high multi-mapping rates (up to 96%)?
**Answer:** This is a signature of **incomplete Ribosomal RNA (rRNA) depletion** during the Input library preparation: rRNA makes up ~90% of total RNA. If the lab's rRNA depletion step is incomplete, the library is dominated by rRNA. Since rRNA genes are highly repetitive, STAR maps them to many locations, and `featureCounts` discards them as multimappers.

For example, `female_mdd_2_input` had **360 Million total reads**, but **96.5% (347.8M)** were discarded as multimappers, leaving only **1.7%** assigned to exons.

---

### Q5: Can we rerun STAR with a stricter filtration method to increase unique assignment rates?
**Answer:** No. 
1. **rRNA repeats** map to multiple locations with 100% identity. Stricter thresholds will not make them unique.
2. **Unique reads with mismatches** (from SNPs or sequencing errors) will be rejected under stricter filters, turning uniquely mapped reads into unmapped ones and **further decreasing** your unique exonic counts.
3. Changing the multimap limit (e.g., `--outFilterMultimapNmax 1`) simply reclassifies multi-mappers as unmapped, leaving the unique counts unchanged.

Unique mapping is determined by the **physical sequence of the reads**, not aligner thresholds.

---

### Q6: How should we filter out these low-quality samples?
**Answer:** Rather than hardcoding a blacklist of sample IDs, we implement an **objective QC Library Size Gate** in `02_dge_analysis/04_deg_analysis.R`:
* **Gate:** `MIN_LIBRARY_SIZE <- 2e6` (2 Million gene-mapped counts), applied to the Input samples at `02_dge_analysis/04_deg_analysis.R:95`.
* **Runtime exclusion:** `colSums(counts(dds_rnaseq))` is compared with the floor, failing subjects are dropped, their `Subject` levels are `droplevels()`d to keep the design matrix full rank, and the names are printed by the `[QC GATE]` message (`:101-112`). No sample list is hardcoded.
* **Retired pairing rule:** the 10 Million floor and the Input+IP pairing-preservation rule belonged to the 78-library two-fraction analysis. The MeRIP/IP arm is out of scope (`README.md:5`), so the gate constrains Input samples only; `female_mdd_2`, `female_ctrl_6`, `male_ctrl_1` is the 10M-era exclusion set, not the shipped 2M result.

---

### Q7: Does the 10 Million threshold exist as a standard in literature?
**Answer:** Yes:
* **ENCODE Consortium Guidelines:** Recommend >20M–25M uniquely mapped reads for bulk RNA-seq but set the warning flag/quality gate at **10 Million uniquely mapped reads**.
* **Conesa et al. 2016 (Genome Biology):** Highlights 20M–30M as standard for mammalian DGE, but establishes 10M as the lower bound below which statistical power degrades heavily for mid-to-low expression genes.
* **Tarazona et al. 2011 (Genome Research):** Shows that below 10M mapped reads, Poisson sampling noise inflates the False Discovery Rate (FDR).

---

### Q8: Why shouldn't we raise the gate threshold to 20 Million?
**Answer:** Raising the gate to 20M would exclude **7 subjects** (18% of your cohort), reducing your sample size from 39 to 32 (only 8 Ctrl vs 8 MDD per sex group).

In statistical genomics, **biological replication (number of subjects) is a much stronger driver of statistical power than read depth**, provided the libraries are above the 10M stabilization threshold. Throwing away 7 biological replicates to eliminate slight depth variation would severely underpower the study.

---

### Q9: Do we have biological replicates in this study?
**Answer:** Yes. Every subject (e.g., `male_ctrl_1`, `male_ctrl_2`, etc.) represents a separate, independent human postmortem donor. Thus, you have 39 independent biological replicates, not technical replicates.
