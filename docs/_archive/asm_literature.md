# Methodological literature for m6A–genotype analysis

Why this document exists: the m6A-SNP funnel in module 04 was built largely from
assumption and biological intuition. This file records what the published
literature actually does, so each design choice can be justified or corrected
against a precedent rather than against intuition.

Read this before changing any threshold or proposing a new analysis.

---

## 1. The three claims that shaped module 04

| Claim in the pipeline | Literature status | Source |
|---|---|---|
| A SNV ablating the DRACH offset-2 adenine removes the m6A site | **Supported in principle, but motif presence ≠ methylation.** Only ~5% of available DRACH motifs are methylated in vivo. | §2.1 |
| m6A state is genetically variable and trait-relevant | **Supported.** m6A-QTL contribute to disease heritability at levels comparable to splicing QTLs. | §3.1 |
| Single-base m6A changes are detectable in MeRIP data | **Not supported at this resolution.** Detectable effects are orders of magnitude below those often reported; MeRIP cannot quantify methylation fraction. | §3.3 |

---

## 2. Motif context: DRACH is necessary but not sufficient

### 2.1 The ~5% figure

Single-molecule direct RNA sequencing across species reports that m6A
modifications are found in only **approximately 5% of available DRACH motifs**.
Sequence specificity alone does not explain m6A distribution.

> "the presence of the consensus sequence DRACH … Nevertheless, the sequence
> specificity alone cannot explain the observed selectivity in m6A distribution,
> as m6A modifications are only found within a relatively small fraction of the
> available DRACH motifs, accounting for approximately 5%."
> — *Single-Molecule Direct RNA Sequencing Reveals the Shaping of Epitranscriptome
> Across Multiple Species* (bioRxiv 2023.11.16.567334)

**Consequence for module 04.** We scan 602,657 DRACH motifs within IP peaks.
Calling the resulting set "m6A sites" overstates what is known. The defensible
label is **"adenosines in a DRACH context within an IP-enriched region."** A peak
raises the prior above 5% but does not establish it, and we hold no
single-base-resolution methylation measurement to confirm any individual site.

### 2.2 Non-DRACH m6A sites are real and numerous

DeepRM maps m6A across human cell lines and finds **10.1% (HEK293T) / 10.2%
(HeLa) of m6A sites lie in non-DRACH contexts**, consistent with miCLIP2 (6.33%)
and GLORI (14.7%). Almost all (96%) deviate from DRACH by a single nucleotide,
with GGACG and GGAUU most frequent.

Critically, non-DRACH sites differ in kind, not just in sequence:

> "The m6A modification stoichiometry distribution of the non-DRACH sites was
> unimodal and centered at low modification stoichiometry, in contrast to the
> clear bimodal distribution observed for the DRACH sites."

> — *Comprehensive discovery of m6A sites in the human transcriptome using DeepRM*,
> Nature Communications (2025) [10.1038/s41467-025-67417-w]

**Consequence.** A DRACH-focused analysis is not obviously wrong — DRACH sites
are the high-stoichiometry, bimodal ones, so they are the more tractable target.
But the funnel excludes ~10% of genuine m6A sites by construction. This is a
defensible scoping decision, not a complete survey, and should be stated as such.

### 2.3 Methylation is co-transcriptional and promoter-biased

m6A deposition occurs both co-transcriptionally and post-transcriptionally, and
the distribution shifts by transcript region, with strong enrichment in CDS and
3′UTR and near stop codons.

> — *Comprehensive analysis of m6A-seq data reveals distinct features of conserved
> and unique m6A sites in mammals*, RNA (2025)

---

## 3. How genotype → m6A is actually tested

Two independent groups have solved this. Both converge on the same principle:
**effect must be consistent across independent samples.**

### 3.1 m6A-QTL: peak-level, permutation-calibrated

Zhang et al. mapped m6A QTL in 60 Yoruba LCLs.

> "We found that m6A QTLs are largely independent of expression and splicing QTLs
> and are enriched with binding sites of RNA-binding proteins, RNA
> structure-changing variants and transcriptional features. … Through integration
> with data from genome-wide association studies, we show that m6A QTLs contribute
> to the heritability of various immune and blood-related traits at levels
> comparable to splicing QTLs and roughly half of expression QTLs."
> — *Genetic analyses support the contribution of mRNA N6-methyladenosine (m6A)
> modification to human disease heritability*, Nature Genetics 52:939–949 (2020)
> [10.1038/s41588-020-0644-z]

Method, precisely:

- **Phenotype is peak-level**, log2(IP/Input) per joint peak — *not* a single base.
- Joint peak calling across samples via **MeRIPtools**, then filtered for
  non-zero reads in every sample.
- Phenotype normalised for **GC content and IP efficiency**; standardized and
  quantile-normalised.
- **Linear model, 15 PCs** as covariates (PC count chosen to maximise power).
- **cis window 100 kb** (most associations fell within this).
- **FastQTL**, 1,000 permutations, **beta approximation** → empirical p per peak.
- **Storey q-value** across peaks for FDR.

The same design was applied to human tissue, including brain, by Xiong et al.
(129 high-quality m6A maps, 91 individuals, 53 brain), with promoter regions
added to the search space because of co-transcriptional deposition.

> — *Genetic drivers of m6A methylation in human brain, lung, heart and muscle*,
> Nature (2021) [PMC9112289]

**Why this matters here.** Peak-level. Our variants sit inside ~300 bp median
consensus peaks. A 1 bp ablation cannot be resolved as a phenotype in this
framework — the closest available phenotype is the whole peak's IP/Input, which
is what step 07 (pooled enrichment) approximates.

### 3.2 ASm6A: allele-specific, concordance-gated

Cao et al. built allele-specific m6A directly from MeRIP-seq — the design that
*does* resolve single nucleotides.

Pipeline: high-confidence heterozygous variants → per-allele read counts in IP
and input → **Fisher's exact test** per site → retain allelic imbalance inside
m6A peaks → require the same bias direction across samples.

> "Genetic variants with high concordance (>90%) were considered as putative
> functional variants. A total of 1184 variants were identified across samples."

Validation was reporter-based: 26 variant pairs plus 8 negative controls in
HEK293T. All >90%-concordance variants reproduced the tissue bias; **only 1 of 8
(12.5%)** negative controls produced any allele-specific modification. Tissue vs
cell-line allelic ratio correlated at **R² = 0.74**.

> — *Allele-specific RNA N6-methyladenosine modifications reveal functional genetic
> variants in human tissues* (PMC10547253)

Two hard requirements we should note:

- **≥2 reads on each allele** ("Loci with reads on any allele below two were discarded")
- **Cross-sample concordance** as the confidence gate, not a p-value threshold

**This is what the retired read-ratio ASM path was reaching for.** It failed for
an unrelated reason — ~600 bp windows pooled unrelated private variants from
different donors into one signal, manufacturing group differences from
singletons. The concept was sound; the window was wrong.

### 3.3 Resolution limits — the hard constraint

> "Although MeRIP-seq can reveal approximate sites of m6A, it cannot be used to
> quantitatively measure the fraction of transcript copies that are methylated."
> "Our results suggest that the scale of statistically detectable m6A changes in
> response to various stimuli is orders of magnitude lower than the scale of
> changes reported in many studies. However, we also found that statistical
> detection could miss the majority of changed sites when using only 2–3
> replicates."
> — *Limits in the detection of m6A changes using MeRIP/m6A-seq*, Scientific
> Reports (2020) [10.1038/s41598-020-63355-3]

Long-read work independently confirms the dilution problem and quantifies the
hit rate:

> "In total, 41 m6A sites had SNPs in D, R, or H positions, with 6 of these (14%)
> classified as ASM."

and of 23 ASM sites detected by long reads, 19 appeared in MeRIP-seq but **only 4
had sufficient coverage** for allele-bias measurement. The paper attributes
inflated site counts in the MeRIP literature partly to antibody specificity
artifacts.

> — *Long-read RNA sequencing reveals allele-specific N6-methyladenosine
> modifications*, Genome Research 35:999 (2025)

**Consequence.** Our measurement geometry is: 1 bp perturbation, 50 bp counting
bins, ~200 bp MeRIP fragments, 300 bp median peak. Against a background where
detectable m6A change is already at the assay's edge, a single-base ablation
inside a 300 bp peak is not a viable per-site phenotype. This is the strongest
argument for a **pooled/aggregate** test rather than per-site IP loss.

---

## 4. Prediction tools and databases (currently unreachable)

Useful for annotation once network access is restored, or as a comparison set.

| Resource | What it provides | Note |
|---|---|---|
| **m6ASNP** (GigaScience 2018) | Random forest over 55,548 human m6A sites from miCLIP; predicts gain/loss. 133,394 gain + 214,884 loss variants from dbSNP. Functional **loss** variants are more conserved and more deleterious than gain variants. | [10.1093/gigascience/giy035] |
| **m6AVar** (NAR 2017) | Curated functional variants affecting m6A. | Web-only; host refused connection |
| **RMVar** | RNA modification–associated variants. | Web-only; host refused connection |
| **m6A-TSHub** | Tissue-specific m6A hub datasets. | Web-only; host refused connection |
| **GWAS Catalog REST** | Disease associations. | HTTP 410/404 |

Note that m6ASNP's conservation finding — that m6A **loss** variants are under
stronger selective pressure than gain variants — is a mild, independent
support for focusing on disruption rather than creation. Module 04 currently
analyses loss only, which aligns.

---

## 5. Design decisions this literature supports or contradicts

| Module 04 choice | Verdict |
|---|---|
| Core adenine = DRACH offset 2 | **Supported.** Transcript offset 2 on both strands. |
| Minus-strand allele complementation | **Required.** Half the disrupting set (2,198/4,413) is minus-strand. |
| Focus on motif **loss**, not gain | **Supported** by m6ASNP conservation analysis. |
| **No recurrence gate** in the funnel | **Supported.** Filtering on donor count before mechanistic steps conditions the m6A set on an arbitrary cutoff. |
| Label level 6 "core m6A ablation" | **Overstated.** Rename to DRACH-context core adenine in an IP peak (§2.1). |
| Per-site `n_mdd > 6, n_ctrl <= 2` selection | **Weakest link.** `n_mdd` *is* the case-control skew; filtering on it and calling the result "MDD-exclusive" approaches defining the answer. Replace with the pooled test. |
| Per-site IP loss as the phenotype | **Not supported** at this resolution (§3.3). |
| Allele-specific design with concordance gating | **The validated path** (§3.2), if pursued. Requires per-allele read counting from BAMs. |

---

## 6. What would make the analysis defensible

Ordered by strength of precedent.

1. **Pooled enrichment with permutation calibration** — closest to §3.1's
   calibration philosophy, applied to the question we can actually answer with
   peak-resolution data. Implemented as step 07.
2. **Peak-level IP/Input association with genotype**, m6A-QTL style — linear
   model, GC + IP-efficiency + PC covariates, FastQTL-style empirical p, Storey
   q. Needs peak-level phenotypes across all 39 donors, which the TRESS bin
   counts can supply.
3. **Allele-specific m6A from BAMs**, ASm6A style — per-allele IP/Input counts
   restricted to heterozygotes, Fisher exact, >90% concordance gate, ≥2 reads
   per allele. This is the only design in the literature that resolves single
   nucleotides, and it is the correct target for a dissertation claim about
   m6A-SNPs.

Option 3 is a project, not a script. It requires per-base alignments from the
HPC (read-only) and a window design that does not pool unrelated variants.
