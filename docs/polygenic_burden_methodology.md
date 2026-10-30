# Polygenic Burden Analysis — Methodology and Design Rationale

**Date:** 17 October 2026
**Module:** `03_variant_calling` (steps 11–15)
**Status:** Methodology locked. Implementation pending.

---

## 1. Question and scope

The client requires a genetic-risk deliverable for the MDD cohort. The agreed scope (`docs/thesis_defense_narrative_log.md` §2.1) is **polygenic pathway burden analysis**: shift the unit of analysis from individual SNPs to cumulative gene-set burden, because single-SNP association testing at N = 39 yields nothing after multiple-testing correction.

This document records *why the analysis is shaped the way it is*. Every quantitative claim below was measured on this cohort's data; every methodological claim is sourced. DOIs were resolved against Crossref or the issuing repository and are listed in §14. Citations that could not be verified are quarantined in §15 rather than guessed.

---

## 2. A polygenic score is not computable in this cohort

This is a hard blocker, not a matter of effort, and it must be stated plainly to the client before any expectation is set.

### 2.1 What a PGS requires

The authoritative protocol defines the required inputs as exactly two datasets:

> "PRS analyses can be characterized by the two input data sets that they require: (i) base (GWAS) data: summary statistics ... genome-wide, and (ii) target data: genotypes and phenotype(s) in individuals of the target sample."
> — Choi et al., *Nature Protocols* 2020;15, doi:10.1038/s41596-020-0353-1

On sample size, the same source is explicit:

> "We recommend performing PRS analyses that involve association testing on target sample sizes of at least 100 individuals (or effective sample sizes over 100 for case/control data) and caution against analyses that utilise base data with low h²snp and small target sample size."

It further specifies GWAS-grade target QC: "genotyping rate > 0.99, sample missingness < 0.02, heterozygosity P > 1x10-6, minor allele frequency (MAF) > 1%, imputation 'info score' > 0.8".

### 2.2 What this cohort actually has

Verified against NCBI E-utilities for GSE275676 (the source series, Mitsuhashi et al., *iScience* 2024, doi:10.1016/j.isci.2024.111316):

| Field | Value |
|---|---|
| Title | m6A-seq in the human ventromedial prefrontal cortex in major depressive disorder |
| Type | Expression profiling by high throughput sequencing; Other |
| Samples | 78 (39 donors × IP + Input) |
| WGS / genotyping array / DNA-derived data | **none** |

Two independent failures therefore apply:

1. **No genome-wide genotypes exist for these donors.** The variants available are RNA-seq calls — restricted to expressed sequence, hard-called, never imputed. A PGS scored on these would assign missing values across most of the genome and would be biased toward highly-expressed genes.
2. **N = 39 is below the recommended floor of 100** for PRS association testing, before any cohort-specific consideration.

### 2.3 What is already solved

The *base* data — historically the hard part — is free and unrestricted. The Psychiatric Genomics Consortium MDD2025 release is **CC BY 4.0**, DOI 10.6084/m9.figshare.27061255.v5, covering roughly 500,000 participants across multiple phenotype definitions (EHR, clinical/interview, questionnaire), including variants excluding UK Biobank to permit independent validation.

**Implication for the client:** a genuine PGS becomes achievable for this cohort the moment DNA from the postmortem tissue is genotyped or sequenced and imputed. Nothing about the analysis design is blocked — only the tissue. This should be raised as a scope question, not attempted with substitutes.

---

## 3. Rare-variant burden methods are structurally inapplicable

The obvious substitute — rare-variant gene-set burden — also fails, for a reason specific to how these variants were called.

### 3.1 The measured allele frequency spectrum

Derived from the cohort variant table `results/03_variant_calling/tables/09_functional_annotation.csv` (61,685 variants), using AF = (n_het + 2·n_hom_alt) / (2·n_subjects):

| MAF band | All variants | Coding (HIGH/MODERATE/LOW) |
|---|---|---|
| no alt alleles observed | 4,678 (7.6%) | 621 (5.4%) |
| **< 1%** | **0 (0.0%)** | **0 (0.0%)** |
| 1–5% | 1,649 (2.7%) | 315 (2.7%) |
| 5–20% | 9,309 (15.1%) | 1,431 (12.5%) |
| **> 20%** | **46,049 (74.7%)** | **9,126 (79.4%)** |

**There is not one variant below 1% MAF in this cohort.**

### 3.2 Why

This is calling-method bias compounded by an upstream filter. RNA-seq variant calling recovers only variants visible in expressed sequence at sufficient depth; a genuinely private variant in a gene covered at 5× will not pass calling thresholds. Direct measurement of this effect: RNA sequencing recovers approximately 40% of the variants identifiable by DNA sequencing, at greater than 80% precision (Leonard et al., *BMC Genomics* 2024, doi:10.1186/s12864-024-10645-z).

A second, more specific cause applies to the table used here: `09_functional_annotation.csv` contains only variants seen in ≥8 of 39 donors, so low-recurrence variants were removed before this analysis ever ran. See §4.2.

The practical consequence is that this dataset contains **common** germline variants that RNA-seq happens to observe, not the rare variants that burden testing is designed for.

### 3.3 The methods that would be inapplicable

Every dominant gene-set burden framework is gated on minor allele count at or below 1%:

- **SKAT** (Wu et al., *Am J Hum Genet* 2011;89, doi:10.1016/j.ajhg.2011.05.029) and **SKAT-O** (Lee et al., *Am J Hum Genet* 2012;91, doi:10.1016/j.ajhg.2012.06.007)
- **SAIGE-GENE** (doi:10.1038/s41588-020-0621-6) and **SAIGE-GENE+** (*Nat Genet* 2022;54, doi:10.1038/s41588-022-01178-w), which uses multiple MAF cutoffs and saddlepoint approximation for MAF ≤ 1%
- **STAAR** (*Nat Genet* 2020;52, doi:10.1038/s41588-020-0676-4), which aggregates burden across functional annotations

Applied to a dataset containing zero variants below 1% MAF, these would filter out the entire cohort. Running them would produce an empty result that could be misread as evidence of absence.

A second, independent objection to the classic burden statistic: it assumes that variants within a set act in the same direction with similar magnitude (Lee et al., *Biostatistics* 2012;13, doi:10.1093/biostatistics/kxs014). That assumption is wrong for MDD, a highly polygenic trait whose architecture is dominated by many small common effects rather than a few rare large ones — a point also made for m6A-related heritability by Zhang et al. (*Nat Genet* 2020;52, doi:10.1038/s41588-020-0644-z).

### 3.4 Reframe

The defensible framing is **common-variant genetic load**: how much genotyped variation does each donor carry in a given gene or gene set, and does that differ by case status after controlling for detection probability. This is honest about what the data contains, and it sits closer to what a PGS does than to what rare-variant burden testing does.

---

## 4. The variant table is pre-filtered to ≥8 carriers

### 4.1 An earlier claim, and why it was wrong

An earlier draft of this document reported that `09_functional_annotation.csv` contained **1,942 rows with byte-identical carrier counts (`n_subjects=8`, `n_mdd=4`, `n_ctrl=4`)**, described them as placeholder padding, and excluded them. That finding was incorrect and the filter has been removed. It is recorded here because the same misreading produced an incorrect attribution in §3.

The error was generalising from a partial view: the tail of the file was printed, where low-carrier rows had clustered, and a run of identical values was read as fabrication without checking the distribution.

**Three measurements refute the placeholder claim:**

| Check | Result |
|---|---|
| `n_mdd + n_ctrl == n_subjects` | True in **all 61,685 rows** — `n_subjects` is exactly the carrier count |
| Distribution of `n_subjects` | min **8**, median 15, mean 17.2, max 39. 55,224 rows exceed 8 |
| Frequency at each carrier count | 8:6,461 · 9:5,506 · 10:4,816 · 11:4,052 · 12:3,577 · 13:3,253 — smooth decay, **no spike at 8** |

Placeholder padding produces a spike. There is none. 8 is the **minimum**, not a cap.

The group split is likewise unremarkable. `(4,4)` is the modal split, but it sits among near neighbours: `(5,3)` 1,804 · `(5,4)` 1,710 · `(6,4)` 1,328 · `(5,5)` 1,284. A roughly balanced 20/19 cohort produces exactly this.

So 1,942 rows at `(8,4,4)` is the product of the **modal carrier count** and the **modal group split**. Nothing anomalous.

**Cost of the erroneous filter:** 254 of 11,493 coding variants (2.2%) and 8 of 433 panel variants were deleted. No conclusion changed — every result in §9b is null with or without them — but legitimate measurements were being discarded on a premise that did not hold.

### 4.2 The real constraint

`08_build_transcriptome_snp_table.py` bins loci by recurrence:

```
Private (1 Subject)
Low Recurrence (2-4 Subjects)
Moderate Recurrence (5-7 Subjects)
High Recurrence (>=8 Subjects)   <-- the only tier present in this table
```

`09_functional_annotation.csv` contains **only** the ≥8-carrier tier. Its minimum is exactly 8 and no row falls below.

**Consequence:** private and low-recurrence variants are absent by construction. The burden analysis never saw a private variant — not because the assay could not detect one, but because this table deliberately drops them. This is a stronger and more specific statement than the RNA-seq depth argument in §3.2, and it is the operative reason no rare-variant method applies here.

**Handling:** no carrier-count filter. Unplaced contigs (`GL*`, `KI*`, `JH*`, `chrUn*`) are still excluded — 116 rows — because those are genuinely not assignable to a locus.

---

## 5. Gene-set selection methodology

This is the design decision most likely to draw reviewer scrutiny, and the first proposal was rejected for good reason.

### 5.1 Design principle: no hand-picked gene sets

Gene sets must not be assembled by picking entries whose names sound biologically relevant. That practice fails in three ways:

1. **Selection bias by construction.** A set chosen for thematic fit cannot be distinguished, by its construction, from a set chosen to produce a result.
2. **Arbitrary boundary.** The boundary depends on whoever assembled it, so it is neither reproducible nor auditable — and whatever was left out is invisible.
3. **Non-independence.** Pathway databases partition overlapping gene collections (Reactome alone offers 7,670 sets), so false-discovery rates computed across them are not meaningful.

Two of these three are design problems rather than outright errors, which is why curated pathway resources are retained below as a clearly-labelled exploratory analysis. What is prohibited is their use as the primary basis for inference.

One practical note recorded during design: the originally requested set `HALLMARK_SYNAPTIC_PLASTICITY` does not exist in current MSigDB. Verified against msigdbr 26.1.1, which returns 50 Hallmark sets, none of them that one. The set has been retired upstream, so the requested gene set was not buildable from Hallmark as literally written.

### 5.2 Primary analysis: Open Targets, MDD genetic evidence

**Source: Ochoa et al., *Nucleic Acids Research* 2023;51, doi:10.1093/nar/gkac1046. Disease `MONDO_0002009`, restricted to the `genetic_association` evidence datatype.**

Three properties make this the primary basis for inference:

- **External.** Derived from large-scale GWAS, entirely independent of these 39 donors.
- **Non-circular.** Filtering to *genetic evidence only* deliberately excludes literature and clinical evidence, which would otherwise echo the same publications this project cites.
- **Quantitative, not thematic.** Membership is determined by a published evidence score, and the threshold is declared before the data are examined.

Measured panel sizes:

| Genetic evidence score | Panel genes | Present in cohort | Qualifying variants |
|---|---|---|---|
| ≥ 0.50 | 204 | 103 | 222 |
| **≥ 0.35 (primary threshold)** | **384** | **181** | **441** |
| ≥ 0.20 | 559 | 239 | 548 |
| ≥ 0.10 | 673 | 288 | 640 |

Highest-ranked genes are established MDD loci — `NEGR1` (0.918), `CTTNBP2` (0.898), `LRFN5` (0.896), `SORCS3` (0.872), `TCF4` (0.867), `STK24` (0.862), `NCAM1` (0.858), `DRD2` (0.854). Genes in the panel that carry variants in this cohort include `ANK3`, `CACNA1C`, `ASTN2`, `MAGI2`, `PCLO`, `CSMD1`, `NAV2`, `TENM2`, `HTT` — a useful sanity signal that the external panel and the cohort data describe the same biology.

The panel is cached to disk with a pinned retrieval date, so the analysis does not depend on a live API call. This matters practically: during this project m6AVar, RMVar and m6A-TSHub all refused connection from this host, so network-dependent steps have already proven fragile here.

### 5.3 Secondary analysis: Hallmark and Reactome (exploratory only)

Curated pathway resources are retained as an **exploratory, post-hoc analysis**, reported alongside but never load-bearing for inference:

- Source: msigdbr 26.1.1, exported with the collections kept **separate**: HALLMARK 50 sets, REACTOME 1,839, KEGG 844, WIKIPATHWAYS 925, BIOCARTA 292, GO_BP 7,538, GO_CC 1,080, GO_MF 1,872 (14,440 exported, 11,234 with ≥3 genes overlapping the cohort).
- The three Gene Ontology branches are tested as three collections, not one. GO:BP alone is 5,831 tests and GO:CC 858, so pooling them would let BP sets set the correction threshold for CC and MF sets that carry different gene definitions.
- Every set is tested and reported, including nulls. Nothing is filtered after seeing results.
- Multiple-testing correction is applied across tested sets within each collection, and the correction denominator is reported per collection.
- Results from this analysis are labelled exploratory in every output table and figure, and carry no claim that survives without corroboration from the gene-level or panel analyses.

This analysis exists because pathway-level burden is a biologically interpretable unit that the Open Targets panel cannot supply — the panel is disease-specific and will never contain a pathway such as synaptic adhesion that is simply not implicated in MDD GWAS. It is worth computing. It is not worth believing on its own.

Two independent cautions apply, both already recorded in §5.1: the collections are non-independent, so the corrected FDR should be read as a ranking device rather than a strict error rate; and sets were designed for expression enrichment, not variant collapsing, so they carry no validation as burden units.

### 5.4 Gene Ontology as an alternative neuronal definition

Gene Ontology is already in the analysis as three separate collections (GO:BP, GO:CC, GO:MF), tested independently so that the 5,831 BP tests do not set the correction threshold for the 858 CC tests. The GO sets are matched to synapse and neurotransmission vocabulary by name, which is the pragmatic reading of this section's request but is a name match, not a subtree expansion.

If a formally complete neuronal definition is wanted, the defensible construction is a **complete subtree** of the Gene Ontology rooted at one canonical term, expanded mechanically — not a list of neuronal-sounding keys. Candidate roots: `GO:0045202` synapse (1,659 directly annotated genes, verified), `GO:0044459` neuron part, `GO:0043005` neuron projection. Taking the whole subtree means synaptic plasticity, glutamatergic and GABAergic terms are included by construction and nothing is left out by oversight.

Not implemented in the current pipeline. Recorded here as the recommended route should a neuron-restricted analysis be required later.

---

## 6. Normalization

This is the part that explains the previously reported inconsistent results, and it is the substantive methodological contribution of this analysis.

### 6.1 The measured problem

Counting variants per donor directly from `_data/variant/filtered_vcf/`:

| Group | n donors | Mean variants per donor |
|---|---|---|
| MDD | 20 | 62,472 |
| Control | 19 | 49,825 |
| **Ratio** | | **1.254×** |

Per-donor range across the cohort: **7,611 to 131,865 — a 17-fold spread.**

Burden is a count. Unnormalized, MDD donors read as more burdened by 25% before any biology is involved. Normalize it, and the set of genes clearing any filter shifts, so the headline result changes with the normalization chosen. That is not noise; it is the predictable behaviour of an unnormalized count comparison. The identical 1.254× ratio was independently measured in the `04_asm` case-control test, so this is one confound presenting in two modules.

### 6.2 Layer 1 — count model with library-size offset (primary)

Burden is modelled as a count response with total detected variation per donor as a known offset. This is the same normalization principle used for RNA-seq read counts, where library size enters as `offset = log(total)` rather than as a fitted covariate — see DESeq2 (Love, Huber & Anders, *Genome Biology* 2014;15, doi:10.1186/s13059-014-0550-8) and edgeR (Robinson, McCarthy & Smyth, *Bioinformatics* 2010;26, doi:10.1093/bioinformatics/btp616).

```
log E[burden_{d,G}] = log(total_variants_d) + β0 + β1·condition + β2·sex + β3·(condition × sex)
```

Negative binomial, not Poisson. Overdispersion is guaranteed given the 17-fold spread; Poisson would understate the variance and produce anti-conservative p-values.

### 6.3 Layer 2 — permutation null preserving call count

Condition labels are permuted **within sex**, so every donor retains their own call count. This is the mechanism that prevents the depth imbalance from masquerading as signal, and it is already implemented in `04_asm/06_test_mdd_enrichment.py`. It is implemented inside `14_gene_burden_test.py`, so the script has no import dependency; the archived `04_asm/06_test_mdd_enrichment.py` used the same logic.

Duplicating the permutation would repeat exactly the failure mode found and fixed in `07_allele_specific_m6a.py`, where two implementations of one statistic disagreed between a results CSV and the prose — the CSV reported annotation counts over a different denominator than the concordance-filtered headline.

### 6.4 Layer 3 — the confound the offset does *not* fix

An offset on total variants corrects global depth. It does **not** correct **gene-specific detection probability**. If a donor expresses a gene less, fewer variants are detectable there regardless of genotype. RNA-seq variant calling inherits this expression bias directly, and this cohort already shows differential expression between groups (astrocyte loss, p = 0.0389, per QC-02).

Required mitigations:

- Report burden as a **rate** over variants for which the donor was informative in that set, alongside the raw count.
- Run a sensitivity analysis restricted to gene sets with uniform coverage across donors.
- Treat any gene-set result as descriptive of *detected* variation, never of genetic risk.

### 6.5 Outlier handling

`male_ctrl_1` has 13,970 variants, near the bottom of the 17-fold range, and was already excluded at the RNA-seq analysis under QC-01. All primary tests are run with and without it, and both are reported.

---

## 7. Statistical model and testing structure

### 7.1 Testing structure

| Analysis | Gene set | Role | Multiple-testing correction |
|---|---|---|---|
| **Gene-level burden (unbiased)** | every gene with ≥1 qualifying variant (5,831 tested) | **primary inference** | Benjamini–Hochberg across genes |
| **Panel burden (Open Targets)** | MDD genetic evidence ≥ 0.35 (384 panel genes) | interpretable secondary, declared in advance | panel-level |
| **Pathway burden (curated)** | msigdbr Hallmark / C2 / C5 (§5.3) | exploratory only | BH within each collection; denominator reported |
| **Convergence** | gene-level burden significant ∩ panel member | credible set | — |

Gene-level burden carries the inference because no gene set was chosen by hand. Panel burden is interpretable because it is external and evidence-scored. Together they make the conclusion defensible in either direction: a gene-level hit absent from the panel is exploratory, while a panel gene with no gene-level support is evidence against a burden effect in this cohort.

Pathway burden is reported for completeness and biological interpretability. It is labelled exploratory in every output and supports no claim on its own.

### 7.2 Qualifying variant set

SnpEff impact classes HIGH, MODERATE and LOW (11,486 qualifying variants cohort-wide; panel count in §9b). Annotation follows Cingolani et al., *Fly* 2012;6, doi:10.4161/fly.19695.

LOW impact (synonymous) variants are included deliberately. Synonymous variants are not inert: they alter mRNA secondary structure, codon optimality, translation kinetics, and RNA-binding/m6A motifs. Supporting precedent from this project: the `HAX1` core-ablation locus (chr1:154273441) is annotated `synonymous_variant`, LOW impact, so excluding LOW would discard a locus that the module's own allele-specific analysis recovered independently.

---

## 8. Sensitivity analyses

Every one of these is reported, not cherry-picked:

1. Count versus rate normalization (§6.4).
2. With and without `male_ctrl_1`.
3. Impact strata reported separately: HIGH, MODERATE, LOW, and combined.
4. All four Open Targets thresholds (0.50 / 0.35 / 0.20 / 0.10).
5. MAF strata: 5–20% and >20%.
6. Sex-stratified estimates, and the condition × sex interaction.
7. Permutation count sensitivity (20k versus 100k).

---

## 9. Expected outcome

**A null result is the expected outcome**, and this should be stated before the analysis runs.

A preliminary unnormalized calculation on the panel gives qualifying variants with a mean burden of ~247 carrier-events per donor, MDD and controls differing by ~1%. This apparent agreement is itself misleading: because MDD donors carry 1.254× more variants overall, a depth-proportional MDD burden would be ~25% higher still. The raw numbers resemble a null only because panel burden does not scale with depth the way total load does. Applying the offset tightens rather than reveals.

At N = 39 with a 1.25× coverage imbalance, residual variance after offsetting will be large. A controlled null is nonetheless a real deliverable: it discharges the agreed analysis scope, it quantifies the confound, and "no gene set shows differential burden once detection depth is controlled" is a defensible client answer. The before/after-offset comparison figure is informative on its own, because it shows how much of any apparent burden difference is sequencing depth.

---

## 9b. Result — NULL across all four analyses

Run on 17 October 2026, `N_PERM=20000`, seed 20261017. Every number below is read from a file written by `14_gene_burden_test.py`.

| Analysis | Tests | Result |
|---|---|---|
| **Gene-level burden (unbiased)** | 5,918 genes | 50 nominal p<0.05 (**0.9%**); **0 at BH-FDR<0.05** |
| **Panel burden (Open Targets)** | 4 thresholds | offset RR 1.031–1.075, all p > 0.69; within-sex permutation p 0.70–0.85 |
| **Pathway burden (curated)** | 11,234 sets | **0 at FDR<0.05** in any of the 8 collections — REACTOME 1,497 · GO_BP 5,831 · GO_MF 1,345 · GO_CC 858 · WIKIPATHWAYS 788 · KEGG 643 · BIOCARTA 222 · HALLMARK 50 |
| **Convergence** | — | 0 (gene-level burden produced no FDR-significant gene) |

Sensitivity, all null:

| Analysis | raw ratio | offset RR | p (offset) | RR no offset | perm p |
|---|---|---|---|---|---|
| all genes | 1.014 | 1.071 | 0.713 | 0.944 | 0.752 |
| all genes, exclude `male_ctrl_1` | — | 1.071 | 0.711 | — | — |
| impact HIGH only | 1.006 | 1.028 | 0.879 | 0.951 | 0.558 |
| impact MODERATE only | 1.014 | 1.060 | 0.752 | 0.939 | 0.738 |
| impact LOW only | 1.014 | 1.079 | 0.687 | 0.945 | 0.775 |

**No gene or gene set shows differential burden once detection depth is controlled. No MDD enrichment claim is made.**

### Two features of the result that need stating plainly

**1. The offset changes the sign of the point estimate.** Across every sensitivity row, the unadjusted rate ratio sits at ~0.94–0.95 while the depth-corrected ratio sits at ~1.03–1.08. Both are null, but they are not the same number, and the unadjusted one is the wrong one. This is the quantitative demonstration that the 1.254× depth ratio was distorting earlier comparisons — and it is why "the results were inconsistent" between analyses.

**2. Gene-level burden is under-powered by construction, and its nominal rate is below chance.** Only 0.9% of genes reach nominal p<0.05 where ~5% would be expected by chance. This is not a defect in the correction; it reflects the sparsity of per-gene burden. Median carriage is 29 carriers spread across 39 donors — **under one variant per donor** — so a per-gene case-control difference has almost nothing to detect. Estimated dispersion `alpha` has median 0, i.e. most genes fit Poisson rather than being overdispersed; the conservativeness comes from tiny counts, not from over-correction.

**Consequence for interpretation:** the gene-level analysis cannot distinguish a true null from an underpowered one. The aggregate analyses (panel and pathway) pool many genes and are better powered, and they are null as well — but they too rest on 39 donors. The honest summary is that this cohort cannot resolve burden differences of the size that common-variant MDD architecture predicts, and that conclusion is a statement about cohort size and assay, not evidence that burden is equal.

### What would change the answer

1. Genome-wide genotypes (see §2.3), enabling a real PGS powered by PGC MDD2025 rather than by 39 donors.
2. A larger cohort. Nothing else substitutes for n.
3. Single-base methylation measurement, so "DRACH-core adenine" can be tested for actual methylation rather than motif presence.

---

## 10. Implementation

One responsibility per script, following the module convention that script number equals output number equals execution order. `16b` is R because msigdbr is an R package and `16` is Python; one language per script.

| Script | Sole responsibility | Output |
|---|---|---|
| `11_build_burden_matrix.py` | Parse GT from 39 per-sample VCFs into a donor × variant carrier matrix; join annotation; apply the impact filter | `11_burden_matrix.tsv.gz`, `11_donor_totals.csv` |
| `12_fetch_opentargets_panel.py` | Single API call; cache panel with pinned retrieval date | `12_opentargets_mdd_panel.tsv` |
| `13_export_pathway_sets.R` | Export msigdbr sets, collections kept separate | `13_pathway_sets.tsv` |
| `14_gene_burden_test.py` | Collapse to gene burden; Analyses 1–3; negative binomial with offset; within-sex permutation; BH-FDR within collection; sensitivity sweeps | `14_*.csv` |
| `15_plot_burden_results.R` | Before/after-offset comparison; gene-level burden FDR volcano; panel composition | `15_*.pdf/.png` |
| `16_burden_control_vs_mdd_table.R` | Per-donor and per-set burden tables. Draws nothing | `16_burden_per_donor.csv`, `16_burden_summary.csv` |
| `17_plot_burden_control_vs_mdd.R` | Shared plot helper, sourced by 18–22. Not run directly | — |
| `18_plot_burden_reactome.R` | Reactome: top-10 reference panel + Reactome synapse sets | `18_burden_reactome.pdf/.png` |
| `19_plot_burden_kegg.R` | KEGG | `19_burden_kegg.pdf/.png` |
| `20_plot_burden_go_bp.R` | GO:BP | `20_burden_go_bp.pdf/.png` |
| `21_plot_burden_go_cc.R` | GO:CC | `21_burden_go_cc.pdf/.png` |
| `22_plot_burden_go_mf.R` | GO:MF | `22_burden_go_mf.pdf/.png` |

Plus the case-control null statistics, implemented inside `14_gene_burden_test.py`: `permute_within_sex()`, `empirical_pvalue()`, `offset_count_model()`, `bh_fdr()`.

Figure 18–22 share one layout: the LEFT panel is the ten largest effects **within that collection** (restricted to sets overlapping ≥10 cohort genes so the ranking is not filled by three-gene quotients), the RIGHT panel is that collection's synapse / neurotransmission sets. Both panels are scoped to one collection so a single database's BH-FDR denominator governs the whole figure, and so no half of a figure is a copy shared with the other four. HALLMARK and KEGG have too few sets on the synapse vocabulary for a meaningful strip, so their right panel falls back to ranks 11+ of the same collection — skipping ranks 1–10, which panel A already shows — and the subtitle states the substitution.

Sets are matched to the synapse vocabulary by name, per collection, in `17_plot_burden_control_vs_mdd.R`. The GO:MF vocabulary is deliberately restricted to named neurotransmitter systems (glutamate, GABA, NMDA, kainate, adrenergic, dopaminergic, cholinergic, serotonergic) rather than the bare word "receptor", which in GO:MF is also satisfied by taste and interleukin receptor families and would fill the panel with sets unrelated to synaptic transmission.

```bash
micromamba activate mdd
python3 03_variant_calling/11_build_burden_matrix.py    # ~40 s
python3 03_variant_calling/12_fetch_opentargets_panel.py # seconds, then offline
Rscript  03_variant_calling/13_export_pathway_sets.R
python3 03_variant_calling/14_gene_burden_test.py        # permutation-bound
Rscript  03_variant_calling/15_plot_burden_results.R
Rscript  03_variant_calling/16_burden_control_vs_mdd_table.R
Rscript  03_variant_calling/18_plot_burden_reactome.R    # and 19–22
```

---

## 11. Limitations to state before a reviewer does

1. **Not a PGS**, and not presented as one (§2).
2. **Not rare-variant burden.** Zero variants below 1% MAF (§3.1).
3. **RNA-seq calling bias.** Only common, well-expressed variants are detectable; absence of a variant is not evidence of reference genotype (§3.2).
4. **Expression confounding.** Burden reflects detected variation, which tracks expression, which differs by group (§6.4).
5. **Power.** N = 39, with a female-control stratum of n = 9. Within-sex permutation in the female stratum rests on 10 versus 9.
6. **Postmortem confounders unmodelled.** Sex is stratified; age, PMI and RIN are not yet checked against the sample sheet.
7. **HPC dependency.** Per-sample VCFs and the cohort joint VCF live on `hpc-01` at `/mnt/22T/phongdinh/mdd_RNAseq_pipeline/_data/`. Re-calling is deliberately avoided, so this analysis reads existing files only.
8. **Degenerate rows** in the variant table must be excluded and the exclusion reported (§4).

---

## 12. References

All DOIs below were resolved against Crossref or the issuing repository on 17 October 2026.

1. Choi SW et al. Tutorial: a guide to performing polygenic risk score analyses. *Nature Protocols* 2020;15. doi:10.1038/s41596-020-0353-1
2. Mitsuhashi H et al. Altered m6A RNA methylation profiles in depression implicate the dysregulation of discrete cellular functions in males and females. *iScience* 2024. doi:10.1016/j.isci.2024.111316
3. Leonard A et al. RNA-DNA differences in variant calls from cattle tissues. *BMC Genomics* 2024. doi:10.1186/s12864-024-10645-z
4. Wu Q et al. Rare-variant association testing for sequencing data with the sequence kernel association test (SKAT). *Am J Hum Genet* 2011;89. doi:10.1016/j.ajhg.2011.05.029
5. Lee S et al. Optimal unified approach for rare-variant association testing with application to small-sample case-control studies (SKAT-O). *Am J Hum Genet* 2012;91. doi:10.1016/j.ajhg.2012.06.007
6. Lee S et al. Optimal tests for rare variant effects in sequencing association studies. *Biostatistics* 2012;13. doi:10.1093/biostatistics/kxs014
7. Scalable generalized linear mixed model for region-based association tests in large biobanks (SAIGE-GENE). *Nat Genet* 2020;52. doi:10.1038/s41588-020-0621-6
8. SAIGE-GENE+ improves the efficiency and accuracy of set-based rare variant association tests. *Nat Genet* 2022;54. doi:10.1038/s41588-022-01178-w
9. Dynamic incorporation of multiple in silico functional annotations empowers rare variant association analysis (STAAR). *Nat Genet* 2020;52. doi:10.1038/s41588-020-0676-4
10. Ochoa D et al. The next-generation Open Targets Platform: reimagined, redesigned, rebuilt. *Nucleic Acids Research* 2023;51. doi:10.1093/nar/gkac1046
11. Zhang Z et al. Genetic analyses support the contribution of mRNA N6-methyladenosine (m6A) modification to human disease heritability. *Nat Genet* 2020;52. doi:10.1038/s41588-020-0644-z
12. Love MI et al. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology* 2014;15. doi:10.1186/s13059-014-0550-8
13. Robinson MD et al. edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. *Bioinformatics* 2010;26. doi:10.1093/bioinformatics/btp616
14. Cingolani P et al. A program for annotating and predicting the effects of single nucleotide polymorphisms. *Fly* 2012;6. doi:10.4161/fly.19695
15. Liberzon A et al. The Molecular Signatures Database Hallmark Gene Set Collection. *Cell Systems* 2015;1:417-425. doi:10.1016/j.cels.2015.12.004
16. Psychiatric Genomics Consortium. GWAS summary statistics for major depression (PGC MDD2025), v5. FigShare. doi:10.6084/m9.figshare.27061255.v5. License CC BY 4.0. *(DOI issued by DataCite; not indexed in Crossref. Verified via the FigShare API record.)*

---

## 13. Citations that could NOT be verified

Recorded rather than guessed. Both are referenced elsewhere in this project and should be resolved before any external submission.

| Reference | Used for | Status |
|---|---|---|
| Cao et al. 2023 | >90% cross-donor concordance criterion for allele-specific m6A (`04_asm/07_allele_specific_m6a.py`) | **DOI unverified.** Cited in `TODO.md`; searches returned no matching record. Confirm before citing externally. |
| PRS-CS (Choi, Achour & Speed) | the LD-aware PGS construction tool that would be used *if* tissue genotypes become available | **DOI unverified.** Referenced descriptively in §2.3 only; no DOI asserted. |

---

## 14. Related documents

- `docs/_archive/asm_methods.md` — m6A-SNP module methods; source of the shared permutation requirement
- `docs/_archive/asm_postmortem.md` — the coverage-confound finding this methodology addresses
- `docs/thesis_defense_narrative_log.md` §2.1 — the agreed analysis scope this implements, and §3.3 the retracted per-site carrier claims
- `results/03_variant_calling/tables/09_functional_annotation.csv` — cohort variant table (61,685 rows)