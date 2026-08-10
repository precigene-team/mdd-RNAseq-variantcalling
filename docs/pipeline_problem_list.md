# MeRIP-seq Pipeline — Problem List & Remediation Plan

*Identified via codebase audit against MeRIPseqPipe (Bao et al., 2022) and Mitsuhashi et al. (2024).*

---

## 🟡 Important — Affects Statistical Validity

### 1. Asymmetric rRNA Filtering
- **Module**: `01_core_alignment/03_rrna_filter.sh` → `01_core_alignment/05_star_align.sh`
- **Problem**: `03_rrna_filter.sh` explicitly skips all `m6A` (IP) samples with a `data_type != "RNA-seq"` gate. Only Input libraries are rRNA-depleted before STAR alignment. Consequently, in `05_star_align.sh`, RNA-seq samples map from `_data/rrna_filtered/` while m6A samples map from raw `_data/trimmed_fastq/`.
- **Empirical Audit Result** *(scratch/rrna_audit/rrna_contamination_report.tsv)*: A coordinate-based intersection of all 36 existing MACS2 narrowPeak files against 552 GRCh38 rRNA gene loci showed **0/36 samples exceed 0.02% rRNA peak overlap**. The theoretical prediction of "enormous false-positive rRNA peaks" was **not confirmed**. MACS2's local background modeling (`--slocal`/`--llocal`) partially rescues this by normalizing against local IP signal rather than absolute Input depth, attenuating the expected artifact at rRNA loci.
- **Actual Consequence**: The rRNA asymmetry does not produce visible false-positive peaks in the called set. The real damage is more subtle and operates through three mechanisms:
  1. **Library size normalization bias**: rRNA reads retained in the IP BAM inflate the total IP library size. DESeq2 size factor estimation then systematically under-normalizes mRNA genes in the IP fraction relative to Input, biasing the `Fraction:SampleType` interaction term.
  2. **featureCounts assignment distortion**: rRNA reads map to rRNA gene IDs in the GTF and are assigned count mass to those genes in the IP matrix, inflating their apparent IP expression.
  3. **Background noise floor inflation**: Elevated IP background genome-wide raises MACS2's pseudo-replicate noise estimate slightly, potentially suppressing low-enrichment real peaks at borderline significance.
- **Severity Reclassification**: Downgraded from 🔴 Critical to 🟡 Important. The asymmetry is a real methodological inconsistency that should be corrected for rigorous publication, but it is not invalidating the existing peak calls.
- **Fix**: Remove the `RNA-seq`-only gate in `03_rrna_filter.sh` to run rRNA depletion on **both** fractions. Re-run rRNA filter + STAR alignment on the **36 IP samples only**. The 36 Input BAMs are untouched.

### 2. Gene-Level vs. Peak-Level Quantification
- **Module**: `01_core_alignment/10_quantification.sh` → `02_dge_analysis/04_deg_analysis.R`
- **Problem**: `10_quantification.sh` uses `featureCounts` with the GTF file as annotation, counting reads across entire exon bodies of annotated genes. `04_deg_analysis.R` then runs the DESeq2 `Fraction:SampleType` interaction model on these whole-gene counts.
- **Consequence**: m6A methylation is physically localized (~150bp near stop codons). Counting across thousands of bases of non-methylated exons dilutes the IP enrichment signal, drastically reducing the statistical power of the interaction term and producing systematic false negatives.
- **Fix**:
  1. Run TRESS to generate high-confidence peak coordinates (replaces MACS2 — see Item 4).
  2. Merge TRESS peak coordinates across all 36 samples into a cohort-wide consensus BED file.
  3. Run `featureCounts` with `-a consensus_peaks.bed -F SAF` to generate peak-level counts.
  4. Run the DESeq2 interaction model on the peak-count matrix instead of the gene-count matrix.

---

## 🟡 Important — Affects Statistical Validity

### 3. Missing Age and RIN Covariates
- **Module**: `02_dge_analysis/04_deg_analysis.R`
- **Problem**: Mitsuhashi et al. (2024), using the same GSE275676 postmortem vmPFC cohort, explicitly include `Age` and `RIN` (RNA Integrity Number) as covariates in all DESeq2 and TRESS models. Your current design formula omits them.
- **Consequence**: RIN directly measures RNA degradation, which correlates strongly with the proportion of reads mapping to exons vs. introns. Age is a known confounder in postmortem brain transcriptomics. Omitting both introduces systematic bias in differential expression estimates.
- **Fix**: Retrieve `Age` and `RIN` values from GEO accession `GSE275676`. Add as columns to `_data/_metadata/sample_sheet.csv`. Update design formulas to `~ Age + RIN + SampleType` (and corresponding interaction terms in Phase C).

### 4. Peak Calling Tool: Switch to TRESS, Drop MACS2
- **Module**: `04_peak_calling/01_peak_calling.sbatch`
- **Problem**: The pipeline uses only `MACS2`, which was designed for ChIP-seq and does not model RNA expression-level variation as background.
- **Empirical Evidence**: FRiP audit revealed 60,000–105,000 peaks per sample with FRiP scores of 0.77–0.86 across all 36 samples. Published MeRIP-seq studies report 5,000–20,000 peaks at FRiP 0.2–0.5. MACS2 is over-calling by 4–16× because it calls local enrichment wherever IP/Input ratio exceeds the local background — which at RNA-seq sequencing depths occurs across the entire expressed transcriptome.
- **Should MACS2 be kept alongside TRESS?** **No.** An earlier recommendation to "use both" was incorrect and is retracted here. The reasoning was that MACS2 produces coordinate BED files needed for ChIPseeker, AME, and Guitar. This is wrong on two counts:
  1. TRESS produces genomic coordinates (chr, start, end) as part of its output — those coordinates are directly usable by ChIPseeker, AME, and Guitar.
  2. Feeding 80,000 MACS2-over-called coordinates into downstream annotation or motif tools would propagate the same false-positive problem into the biological interpretation.
- **Consequence of keeping MACS2**: Motif enrichment on 80,000 spurious peaks will dilute the true GGACH signal with background sequence. ChIPseeker annotation will show false 3' UTR enrichment. FRiP metrics become meaningless as a QC measure.
- **Fix**: Replace MACS2 entirely with `TRESS` in R. Use TRESS-derived coordinates for all downstream steps.
  - Parameters: `lowcount = 15`, `Binsize = 50`, `FDR < 0.05`, `FC > 1.5`
  - TRESS handles both peak detection and differential methylation in one framework using a hierarchical negative binomial model with replicate-aware variance estimation.

---

## 🟢 Gaps — Required for Biological Validation

### 5. No Motif Verification (GGACH / AME)
- **Module**: Missing — no motif enrichment step exists
- **Problem**: No computational check confirms that called peaks actually contain the m6A consensus `GGACH` motif. Without this, there is no proof that the antibody specifically captured m6A-modified transcripts rather than random fragments.
- **Fix**: Run `AME` from MEME Suite on peak sequences with a `GGACH` motif matrix. Calculate accumulative proportion per sample as in Mitsuhashi et al. (2024).

### 6. No Cross-Sex RRHO Comparison
- **Module**: Missing — no RRHO step exists in `02_dge_analysis/`
- **Problem**: The paper uses Rank-Rank Hypergeometric Overlap (`RRHO2`) to compare male vs. female methylation profiles genome-wide, identifying sex-specific convergent and divergent dysregulation patterns.
- **Fix**: Run `RRHO2` R package on your male and female DMR ranked lists (ranked by `-log10(p) * sign(log2FC)`).

---

## 🔬 Novel Contributions — Analyses Beyond the Published Paper

Mitsuhashi et al. (2024) does not perform any variant-level or genetic analysis. Your GATK pipeline (`03_variant_calling/`) already produces the foundational data to enable two genuine novel contributions.

### 7. Allele-Specific Methylation (ASM)
- **What it is**: At each heterozygous SNP position within an m6A peak, test whether the maternal and paternal alleles carry different m6A methylation levels. A SNP that disrupts a `DRACH` motif on one allele will selectively eliminate METTL3 recognition on that allele, producing a measurable asymmetric IP/Input allele ratio.
- **Why it is novel**: Mitsuhashi et al. report only *average* methylation across subjects. ASM operates at the allele level within individuals, distinguishing *genetically encoded* methylation changes from *environmentally driven* ones. No postmortem MDD study has published allele-specific m6A analysis.
- **What you already have**: GATK HaplotypeCaller VCFs on all 36 Input samples from `03_variant_calling/`. The remaining step is running `GATK ASEReadCounter` on the IP BAMs at heterozygous SNP positions to extract allelic IP read counts.
- **Tool**: M6Allele (Zhang, Z., Tang, K., et al., 2025, GigaScience 14:giaf040, DOI: [10.1093/gigascience/giaf040](https://doi.org/10.1093/gigascience/giaf040)) provides a Bayesian model specifically for this, accepting GATK VCF + IP/Input BAM pairs as input.
- **Caveats**: RNA-seq-derived SNPs have ~5–10% higher FDR than matched WGS. Filter strictly: PASS filter, biallelic, heterozygous only. m6A-QTL (genotype → peak intensity) requires 100+ samples — **do not attempt with n=36.**

### 8. m6A-SNP Motif Disruption Analysis
- **What it is**: Intersect your cohort's called SNPs (from GATK VCFs) with the `DRACH` motif positions inside called m6A peaks. A variant that converts a `DRACH` site to a non-DRACH sequence is a direct candidate causal variant for loss-of-methylation at that site.
- **Why it is novel**: Directly links the genetic architecture of the postmortem MDD cohort to the observed epitranscriptomic dysregulation, providing mechanistic evidence rather than purely statistical association.
- **What you already have**: GATK VCFs on Input samples + TRESS peak BED coordinates. This is a pure coordinate intersection — no new sequencing or alignment required.

---

## ✅ What Does NOT Need to Change

| Module                                     | Status              | Reason                                                                                |
| :----------------------------------------- | :------------------ | :------------------------------------------------------------------------------------ |
| 36 Input (RNA-seq) BAMs                    | ✅ Keep as-is        | Correctly trimmed, rRNA-filtered, and STAR-aligned                                    |
| GATK variant calling on Input samples      | ✅ Keep as-is        | HaplotypeCaller on Input is the correct and only required approach for ASM            |
| GATK on IP samples                         | ✅ Not needed        | ASM uses Input-derived SNP coordinates + `ASEReadCounter` on IP BAMs                  |
| DESeq2 statistical engine                  | ✅ Keep              | Correct choice; used by the paper for DEG                                             |
| Sex-stratified design                      | ✅ Keep              | More rigorous than the paper (paper simply stratifies; you model an interaction term) |
| `04_deg_analysis.R` interaction model math | ✅ Keep architecture | Architecturally sound; only needs covariate addition                                  |
| `01_gene_mapping.R` BioMart annotation     | ✅ Keep              | Fully valid biotype filtering and HGNC symbol mapping                                 |

---

## 📋 Execution Priority Order

```
── REMEDIATION (fix existing pipeline) ─────────────────────────────────────
 1. Fix 03_rrna_filter.sh (remove IP skip gate)              [Important]
 2. Re-run STAR alignment on 36 IP samples only              [Important]
 3. Run TRESS peak calling on corrected IP BAMs              [Important]
 4. Build cohort consensus BED from TRESS peaks              [Important]
 5. Re-run featureCounts at peak-level (SAF mode)            [Critical]
 6. Retrieve Age + RIN from GEO GSE275676                    [Important]
 7. Add covariates to sample_sheet.csv                       [Important]
 8. Run TRESS differential methylation analysis              [Critical]

── VALIDATION (prove biological specificity) ────────────────────────────────
 9. Run AME motif enrichment on TRESS peak sequences (GGACH) [Validation]
10. Run RRHO2 cross-sex comparison on DMR ranked lists       [Validation]

── NOVEL CONTRIBUTIONS (beyond the published paper) ─────────────────────────
11. Run GATK ASEReadCounter on IP BAMs at VCF SNP positions  [Novel]
12. Run M6Allele ASM analysis on IP/Input BAM pairs          [Novel]
13. Intersect GATK VCFs with DRACH motifs in TRESS peaks     [Novel]
```

---

*Reference: Mitsuhashi et al. (2024) iScience 27(12) — DOI: [10.1016/j.isci.2024.111316](https://doi.org/10.1016/j.isci.2024.111316)*
*Reference: Bao et al. (2022) MeRIPseqPipe — DOI: [10.1093/bioinformatics/btac025](https://doi.org/10.1093/bioinformatics/btac025)*
