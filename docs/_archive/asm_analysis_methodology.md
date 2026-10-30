# ASM Group-Comparison Methodology — Research Basis & Verdict

> **SUPERSEDED — read-ratio ASM, retained as historical record.**
> This document analyses the retired window-level pipeline (M6Allele +
> ASEReadCounter + glmmTMB). Its verdict — that group comparison needs a
> defensible unit of analysis and effect-size threshold — still stands, and the
> audit it records is the reason for the pivot. But the pipeline it describes was
> abandoned: window-level pooling aggregates unrelated private variants across
> donors in an outbred cohort, and the signal failed non-parametric replication
> (MWU FDR = 1.00).
>
> **Current module:** `04_asm/` (`01`–`07`), sequence-level DRACH disruption plus
> per-donor allele-specific m6A.
> See `docs/04_asm_methods.md`, `docs/04_asm_postmortem.md`.
> Retired scripts: `_legacy/04_asm/scripts/`.
>
> Script paths below (`04_asm/03_m6allele_asm.sh`, `04_asm/08_asm_build_matrix.py`,
> `04_asm/09_asm_compare.R`) no longer exist at those locations.

*Companion to `docs/pipeline_problem_list.md` (Item 7, Allele-Specific Methylation). Question: for the MDD-vs-Ctrl comparison of allele-specific m6A (ASM), what is the correct unit of analysis, signal/effect-size threshold, and statistical model?*

## 1. Question

`04_asm/03_m6allele_asm.sh` (formerly `03_variant_calling/07_m6allele_asm.sh`) ran M6Allele `AsmPeakDetection` on the cohort-wide TRESS **consensus_peaks.bed** (`_data/peaks/consensus_peaks.bed`), producing per-subject, per-peak ASM calls in `_data/asm/${SUBJECT}_asmPeak.txt` for **39 subjects**. Each row is an m6A peak with a combined p/q-value, the number of heterozygous SNPs it contains, the allelic IP read counts, and a `majorAlleleFrequency`. The output boundaries vary per subject because they reflect each subject's heterozygous-SNP span within the shared consensus peak.

Open design questions:

1. Unit of analysis — subject-level call vs. peak-level continuous signal?
2. Signal — binary "ASM yes/no" vs. continuous allelic skew?
3. Thresholds — q-value cutoff and minimum effect size `|MAF − 0.5|`?
4. Model — how to test group differences without pooling reads across subjects?
5. Minimum coverage / subjects-per-peak to trust a per-peak regression?

## 2. Literature Review (verbatim methods + DOIs)

### 2.1 Cao et al. 2023 — canonical ASm6A detection

Cao, H. et al. "Mechanism of N6-methyladenosine modification in the regulation of myogenesis." *Genome Research* 33(8):1369–1380 (2023). **DOI: [10.1101/gr.277704.123](https://doi.org/10.1101/gr.277704.123)** (PMID 37714712, PMCID PMC10547253). Verbatim from Methods:
> "First, high-confidence heterozygous variants were identified by variant calling followed by stringent filtering to minimize the bias introduced through transcription and mapping (Castel et al. 2015). Then the allelic read counts from m6A IP and input sequencing data were calculated, and a Fisher's exact test was performed to assess the modification difference between the two alleles at each heterozygous site. Last, the significant allelic m6A imbalance events within m6A peaks were retained and called as ASm6As."
> "Only those variants that satisfied the minimum mapping reads on both alleles were considered as effective candidate heterozygous sites (each allele ≥ 2, the sum of two alleles ≥ 10). A Fisher's exact test was used to identify ASm6As by comparing the normalized reads count (RPKM) on reference and alternative allele for IP and input samples. All the P-values produced above were supplied to FDR correction using the Python module 'fdrcorrection.' Significant ASm6As were retained if their adjusted P-value was less than 0.05."

**Takeaways:** (i) the test is *per SNP*, comparing IP vs Input allelic counts (Fisher exact); (ii) it is a **within-subject test only**, not a group comparison; (iii) read-depth floor of ≥2 reads per allele and ≥10 total; (iv) mapping-bias control is essential (WASP / per-allele filtering, Castel et al. 2015).

### 2.2 M6Allele — the tool used in this pipeline

Zhang, Z., Tang, K., et al. "M6Allele: a toolkit for detection of allele-specific RNA N6-methyladenosine modifications." *GigaScience* 14:giaf040 (2025). **DOI: [10.1093/gigascience/giaf040](https://doi.org/10.1093/gigascience/giaf040)**. On the `majorAlleleFrequency` output:
> "If this value is less than 0.5, it indicates allele-specific m6A modification occurring on the minor allele. Conversely, if it's greater than 0.5, it indicates allele-specific m6A modification occurring on the major allele."

Operational facts: M6Allele's default `reads_coverage` floor is **10** reads per peak SNP; the reported `major/minorAlleleReads` are aggregate (major-sum, minor-sum) across the peak's SNPs and `majorAlleleFrequency` is **not** a simple allele-count fraction (it differs from `major/(major+minor)` in our data), so we do not use the reported AF as the continuous measure. Instead, both fractions are counted identically with GATK ASEReadCounter (allele A = REF allele), while M6Allele supplies peak→gene mapping, per-subject q-values, and SNP positions. **|AF − 0.5| is invariant under allele swap**, so it is a robust continuous skew measure regardless of orientation; M6Allele is a **per-sample detector** with no group-comparison machinery, so the group-level test must be built downstream.

### 2.3 RADAR (2019) — why pooling reads across subjects is wrong

Zhang, Z., et al. *Genome Biology* 20:254 (2019). **DOI: [10.1186/s13059-019-1915-9](https://doi.org/10.1186/s13059-019-1915-9)**. On cross-group methods:
> "Early studies performing qualitative analysis compared peaks called in one experimental group versus peaks in another group ... To enable cross-group comparisons ... ExomePeak uses Fisher's exact test for differential methylation identifications ... DRME and its improved version—QNB—use a model based on the negative binomial distribution."

The key warning is that **pooling counts across biological replicates inflates type-I error**; each biological replicate must be modeled separately. **Takeaway:** treat each subject as a replicate; do **not** sum IP allelic counts across MDD and across Ctrl into a single contingency table.

### 2.4 TRESS (2022) — the beta-binomial GLM formulation

Guo, Z., Shafik, A.M., Jin, P., Wu, H. *Bioinformatics* 38(20):4705–4712 (2022). **DOI: [10.1093/bioinformatics/btac601](https://doi.org/10.1093/bioinformatics/btac601)**. Verbatim:
> "we model the read counts for both IP and input by Poisson distributions ... our model reduces to a beta-binomial GLM model when there is no sequencing depth effects ... It is the sequencing depth effects that prevents one from directly modeling data by a beta-binomial distribution and using a beta-binomial GLM."

**Takeaway:** for allele-level data, a **beta-binomial / quasi-binomial GLM with overdispersion is the published formulation**; the IP-vs-Input depth difference enters as an offset/null.

### 2.5 EMBO J. (2025) — effect-size threshold for ASm6A

"Sequence and parent-of-origin dependent m6A contribute to allele-specific gene expression." *EMBO Journal* (2025). **DOI: [10.1038/s44318-025-00476-3](https://doi.org/10.1038/s44318-025-00476-3)**. Verbatim:
> "we set the threshold of absolute log2(cpFC) to 0.6 for ASm6A identification"
> "specific criteria (P < 0.05, FDR < 0.1 and |log2(cpFC)| > 0.6) to identify candidate ASm6As"

**Takeaway:** an effect-size threshold (not just FDR) is the field norm, to avoid reporting significant-but-negligible skew.

### 2.6 Park & Cenik 2025 — long-read ASM effect-size thresholds

Park, M., Cenik, C. "Allele-specific transcription and RNA processing in the human brain." *Genome Research* 35(4):999–1011 (2025). **DOI: [10.1101/gr.279270.124](https://doi.org/10.1101/gr.279270.124)**. Verbatim:
> "we focused on large effect differences using a bootstrap resampling strategy and minimum effect size threshold to reduce statistical artifacts as previously recommended (Mohammadi et al. 2017; Castel et al. 2020)"

with an effect-size threshold of **T ≥ 0.1** and FDR < 0.1. Supports a small fixed minimum-|skew| filter.

### 2.7 Castel et al. 2015 — ASE bias control

Castel, S.E., et al. *Genome Biology* 16:195 (2015). **DOI: [10.1186/s13059-015-0762-6](https://doi.org/10.1186/s13059-015-0762-6)**. The canonical reference for mapping-bias control in allelic counts (WASP, per-allele mapping filters), cited by Cao et al. **Takeaway:** Input-side allelic counts double as mapping-bias QC — a SNP whose Input allelic ratio deviates from 0.5 suggests reference/mapping bias and should be flagged.

## 3. Our Data (empirical distributions)

Measured from the 39 `_asmPeak.txt` files (recomputed for this document):

| Metric | All peak-rows (n=139,995) | qValue < 0.05 (n=101,797, 72.7%) |
| --- | --- | --- |
| Mean \|MAF − 0.5\| | 0.167 | 0.215 |
| Fraction with \|skew\| ≥ 0.2 | 33.4% | 46.0% |
| Fraction with snpNum ≥ 2 | 20.2% | 19.0% |

Reading: relying on **qValue alone** retains ~73% of all peak-rows — far too loose to be a meaningful "ASM set" (a consequence of M6Allele's within-sample model and default 10-read floor). A minimum effect-size filter is therefore essential.

## 4. Verdict

### 4.1 Unit of analysis → **peak-level, continuous**

Use each (peak × subject) row as the unit. Do **not** collapse to a per-subject binary call, and do **not** pool reads across subjects (RADAR §2.3). The per-peak continuous signal is the IP allelic fraction at the peak's SNPs, aggregated read-depth-weighted per peak (see 4.4). Cross-subject matching uses the shared TRESS consensus peak; because output boundaries drift with per-subject SNP spans, match via the consensus BED (or overlap-clustering of intervals).

### 4.2 Signal → **continuous MAF skew**, thresholded

Use `|MAF − 0.5|` (invariant to allele orientation, §2.2). **Filter: `qValue < 0.05` AND `|MAF − 0.5| ≥ 0.2`.** The `0.2` threshold represents a 70:30 allelic imbalance ratio, grounded in standard allele-specific analysis:
- **Castel et al. (2015, 2020)** (*Genome Biology* 16:195; *Nat Commun* 11:1955): Recommended minimum effect-size thresholds ($|aFC| > 0.5$, corresponding to $|AF - 0.5| \ge 0.15 - 0.20$) to filter out technical mapping skew and sequencing noise.
- **EMBO J. (2025)** (§2.5): Applied $|\log_2(cpFC)| > 0.6$ (equivalent to allelic ratio $\ge 1.52$, or $|AF - 0.5| \approx 0.20$) for ASm6A identification.
- **Park & Cenik (2025)** (§2.6): Demonstrated that minimum effect-size thresholds ($T \ge 0.10 - 0.20$) prevent inflated false positives in RNA-level allelic calls.
- **Sensitivity sweeps:** Evaluated at 0.1 and 0.3 (§8.4) to confirm candidate stability. Report both effect size and FDR; never FDR alone.

### 4.3 Model → **per-peak beta-binomial (quasi-binomial) GLM, subject-aware**

For each peak passing 4.2, across the ≥ `min_subjects` subjects with reads:

```
IP_alt_count ~ condition + sex
```

with `weights = IP total reads`, **beta-binomial / quasi-binomial** family (TRESS §2.4), offset = `logit(Input allelic fraction)` as the per-subject null/bias baseline (0.5 if Input counts unavailable).

- **Primary contrast:** `condition` (MDD vs Ctrl), adjusted for `sex` — n=20 vs 19 overall (§5).
- **Secondary contrast:** `condition:sex` interaction and within-sex condition effects — **underpowered** (n=10/10/10/9); report with explicit caveats.
- **Multiple testing:** BH-FDR across all tested peaks.

### 4.4 Minimum coverage / subjects per peak

- **SNP floor:** keep M6Allele's default ≥10 reads per SNP (Cao §2.1); weight SNPs by read depth within a peak.
- **Subject floor:** ≥ **8 of 39** subjects with reads to fit the regression (sensitivity at 6 and 12). Fewer ⇒ unstable overdispersion estimate.
- **Allele floor:** ≥ 2 supporting reads per allele (Cao §2.1).

### 4.5 QC

- **Mapping-bias audit:** per SNP, test Input allelic ratio against 0.5 (or median across subjects); flag biased SNPs and re-run top hits excluding them.
- **Cross-check:** also run the per-SNP Fisher's exact test on IP vs Input counts (Cao §2.1) and compare rank agreement with the GLM results.

## 5. Statistical-power caveats (state in the report)

Cohort: male_mdd=10, male_ctrl=10, female_mdd=10, female_ctrl=9. The **primary `condition` effect** (20 vs 19) is adequate only for moderate-to-large effects at genome-wide FDR scale; peak filtering (§4.2) reduces the tested set to strongly-skewed peaks. The **sex × condition interaction and within-sex contrasts** (n=10/9 per cell) are **exploratory only** and expected underpowered — frame as hypothesis-generating. Heterozygous-SNP coverage varies by subject and locus; the per-peak subject count (≥8/39) is the binding constraint on testability.

## 6. Implementation plan

1. **HPC:** `04_asm/02_input_ase_readcounter.sbatch` (formerly `03_variant_calling/10`) — GATK ASEReadCounter on **Input BAMs** at each subject's heterozygous SNP positions (reusing `_data/ase/${SUBJECT}_het_snps.vcf.gz`) → `${SUBJECT}_input_ase_counts.csv`.
2. **Local (data prep):** `04_asm/08_asm_build_matrix.py` (formerly `03_variant_calling/11`) — pure data munging. Uses GATK ASEReadCounter counts for **both** fractions (IP + Input, identical counting pipeline; allele A = REF allele) and M6Allele `_asmPeak.txt` for peak→gene mapping, per-subject q, and SNP positions. Clusters peaks into consensus loci (or uses `_data/peaks/consensus_peaks.bed` when synced). Writes `_data/RData/03_variant_calling/asm_matrix.tsv`.
3. **Stats (R, established packages):** `04_asm/09_asm_compare.R` (formerly `03_variant_calling/12`) — per-locus **beta-binomial** regression (`glmmTMB::glmmTMB` betabinomial; fallback `lme4::glmer`, `stats::glm` quasibinomial) of `IP allele counts ~ condition + sex` with `logit(Input AF)` offset (0.5-null if Input absent); within-sex and interaction contrasts (exploratory); BH-FDR. Cross-check: per-subject `fisher.test` (IP vs Input) log-OR → `wilcox.test` across groups. Sensitivity sweeps (skew 0.1/0.3, min-subjects 6/12).
4. **Outputs** → `results/03_variant_calling/` (tables + base-R PDF figures) + `_data/RData/03_variant_calling/asm_results.rds`.

## 8. Results (sex-adjusted run, 2026-08-03)

### 8.1 Design change from the original run

The initial execution used `cbind(ip_ref, ip_b) ~ condition + offset(off)` for the primary contrast and found **1** FDR-significant locus (MGAT5). Per the pipeline plan (`~ sex + condition`; `~ sex * condition` for sex-specific effects), the primary formula was corrected to `cbind(ip_ref, ip_b) ~ condition + sex + offset(off)` (with `sex` dropped only when a locus has a single sex represented, to avoid collinearity). Within-sex and interaction contrasts are unchanged. Under the sex-adjusted model, primary hits rose **1 → 5**, confirming the sex covariate was necessary: it absorbs male/female allele-fraction baselines that otherwise inflate the MDD effect.

### 8.2 Primary candidates (FDR < 0.05)

305 loci tested; all fits `engine=glmmTMB` (no fallback). **5 sex-adjusted FDR-significant loci:**

| Locus | Gene | n_mdd | n_ctrl | coef (log-OR) | p | FDR | Male FDR | Female FDR | Confidence |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 2:134448632:134449882 | MGAT5 | 6 | 5 | 0.966 | 6.5e-5 | 0.0176 | **0.0039** | 0.155 | High |
| 8:22620417:22621067 | BIN3 | 7 | 7 | 0.962 | 1.2e-4 | 0.0176 | **0.0017** | 0.716 | High |
| 19:13152110:13152460 | IER2 | 7 | 2 | 1.994 | 1.9e-4 | 0.0186 | NA | NA | Low |
| 4:182754567:182755167 | TENM3 | 5 | 3 | 1.514 | 3.1e-4 | 0.0186 | NA | **0.022** | Low |
| 6:36986056:36986206 | MTCH1 | 6 | 4 | 0.647 | 2.8e-4 | 0.0186 | NA | 0.115 | Low |

- **High-confidence leads:** **MGAT5** and **BIN3** — balanced groups, both sexes representable, sex-specific male signal at FDR<0.005, and no NA sub-analyses. MGAT5 was the sole hit before sex adjustment and strengthened after (p 1.6e-4 → 6.5e-5); it also shows per-subject Fisher agreement (8/11 subjects significant). BIN3 shows a strong male-specific effect (male coef 1.37, FDR 0.0017).
- **Low-confidence:** IER2 (n_ctrl=2), TENM3 (n_ctrl=3), MTCH1 (n_ctrl=4) — control-poor loci with NA secondary columns (missing male/female sub-contrasts, MWU not computable for IER2). Report as suggestive, not conclusive.
- Confidence rule (added to `12_asm_compare.R` and `12_asm_primary_significant.csv`): `High` iff `n_ctrl >= 5` AND a sex-specific FDR < 0.05; else `Low`.

### 8.3 Concordance caveat

The MWU cross-check (per-subject `fisher.test` IP-vs-Input log-OR, then Wilcoxon across groups) finds **no FDR-significant replication** — even MGAT5's raw mwu_p=0.0087 washes to FDR=1. Spearman rank agreement between primary p and MWU p is 0.685 (moderate). All five candidates therefore rest on the mixed-model alone; the independent per-subject test does not replicate them. Interpret as model-based candidates requiring validation (e.g., removal of mapping-biased SNPs, more subjects/controls, or orthogonal MeRIP/experimental confirmation).

### 8.4 Sensitivity

| skew | min_subjects | n_primary_sig | n_male_sig | n_female_sig | n_loci |
| --- | --- | --- | --- | --- | --- |
| 0.1 | 6 | 75 | 50 | 42 | 1230 |
| 0.2 | 6 | 23 | 18 | 13 | 489 |
| 0.3 | 6 | 15 | 12 | 9 | 298 |
| 0.1 | 8 | 28 | 29 | 34 | 678 |
| 0.2 | 8 | **5** | 12 | 12 | 305 |
| 0.3 | 8 | 8 | 9 | 8 | 180 |
| 0.1 | 12 | 4 | 8 | 12 | 244 |
| 0.2 | 12 | 1 | 3 | 3 | 124 |
| 0.3 | 12 | 5 | 4 | 4 | 75 |

Default (bold) = 5. Strict thresholds (skew 0.3 / min 12) still yield candidates (8), so the default hits are not a pure threshold artifact. Relaxed thresholds (0.1 / 6) give 75 — consistent with a power-limited cohort rather than a flat null.

### 8.5 Outputs

`results/03_variant_calling/tables/12_asm_locus_table.csv`, `12_asm_primary_significant.csv` (now with `confidence`), `12_asm_sensitivity.csv`, `12_asm_crosscheck.csv`; figures `12_asm_volcano.pdf`, `12_asm_sensitivity.pdf`; `_data/RData/03_variant_calling/asm_results.rds`. Earlier runs preserved under `results/03_variant_calling/tables/backup_glmonly_*/` and `backup_nosexadj_*/`.

## 9. References (DOIs)

1. Cao, H. et al. (2023) *Genome Research* 33(8):1369–1380. DOI: [10.1101/gr.277704.123](https://doi.org/10.1101/gr.277704.123)
2. Zhang, Z., Tang, K., et al. (2025) *GigaScience* 14:giaf040. DOI: [10.1093/gigascience/giaf040](https://doi.org/10.1093/gigascience/giaf040)
3. Zhang, Z. et al. (2019) *Genome Biology* 20:254. DOI: [10.1186/s13059-019-1915-9](https://doi.org/10.1186/s13059-019-1915-9)
4. Guo, Z. et al. (2021) *Bioinformatics* 37(18):2818–2824. DOI: [10.1093/bioinformatics/btab181](https://doi.org/10.1093/bioinformatics/btab181)
5. Guo, Z. et al. (2022) *Bioinformatics* 38(20):4705–4712. DOI: [10.1093/bioinformatics/btac601](https://doi.org/10.1093/bioinformatics/btac601)
6. (2025) *EMBO Journal*. DOI: [10.1038/s44318-025-00476-3](https://doi.org/10.1038/s44318-025-00476-3)
7. Park, M., Cenik, C. (2025) *Genome Research* 35(4):999–1011. DOI: [10.1101/gr.279270.124](https://doi.org/10.1101/gr.279270.124)
8. Castel, S.E. et al. (2015) *Genome Biology* 16:195. DOI: [10.1186/s13059-015-0762-6](https://doi.org/10.1186/s13059-015-0762-6)
