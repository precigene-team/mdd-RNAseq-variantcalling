# Module 04 — m6A-SNP Methods, Trials and Interpretation

Companion to the `04_asm/` pipeline. Rationale, design decisions and negative
results live here; the scripts carry only a short `Goal`/`Steps` header.

---

## 1. Why the read-ratio ASM path was retired

The original approach quantifies allele-specific methylation from MeRIP-seq read
counts: GATK `ASEReadCounter` → M6Allele peak quantification → beta-binomial
`glmmTMB` over ~600 bp peak windows. A data-integrity audit found three problems.

1. **Window-pooling lumping.** In outbred humans a ~600 bp window aggregates many
   unrelated private SNPs from different donors into one "locus". This manufactures
   apparent group differences from singletons — the `BIN3` core variant is carried by
   **one donor**, not the 14 the window implied.
2. **Replication failure.** The independent Mann-Whitney U cross-check returns
   **FDR = 1.00**. None of the five FDR-significant loci replicate; rank agreement
   between the mixed model and MWU is 0.685.
3. **Fragment dilution.** The ~200 bp antibody pull-down inherently dilutes
   single-base resolution, so per-base IP/Input allele ratios are not interpretable.

Scripts and results are retained in `_legacy/04_asm/` as the evidence base for this
critique. Nothing there should be re-run.

---

## 2. The current definition

An **$m^6A$-SNP** is a variant that, inside an experimentally defined m6A peak,
abolishes METTL3/14 writer competence by destroying part of the catalytic `DRACH`
consensus. This is sequence-level and needs no MeRIP read-counting.

| Genotype | Motif state | Outcome |
|---|---|---|
| `0/0` | intact on both alleles | no disruption |
| `0/1` | one intact, one abolished | allele-specific loss |
| `1/1` | abolished on both alleles | complete biallelic ablation |

The het-only restriction was removed: `0/1` and `1/1` are both valid $m^6A$-SNPs,
and removing it unifies Module 03 and Module 04 on the same cohort VCFs.

### The strand bug this replaced

The retired annotator compared **genomic** alleles against **transcript** IUPAC sets
without complementing minus-strand motifs, inverting every minus-strand verdict — a
real `T>A` at a minus-strand core adenine read as *preserving* when it is a core
ablation. **2,339 of the 4,709 disrupting overlaps are minus-strand**, so this was not
a corner case. `04_annotate_strand_aware_disruption.R` complements with Biostrings
before the IUPAC test.

---

## 3. Pipeline

| Step | Script | Output |
|---|---|---|
| 01 | `01_find_drach_sites.R` | `_data/peaks/drach_sites.bed` |
| 02 | `02_verify_drach_reference.py` | `02_drach_reference_qc.csv` (+ `_failures.csv` on failure) |
| 03 | `03_intersect_drach_genotypes.sh` | delegates to 02 then 04 |
| 04 | `04_annotate_strand_aware_disruption.R` | `04_drach_variant_overlaps.csv`, `04_drach_overlap_summary.csv` |
| 05 | `05_plot_m6asnp_funnel.py` | `05_funnel_levels.csv`, `05_core_m6a_sites.csv`, `05_core_pool_recurrence.csv`, `05_m6asnp_funnel.pdf/.png` |
| 06 | `06_test_mdd_enrichment.py` | `06_mdd_enrichment_summary.csv`, `06_mdd_enrichment_strata.csv`, `06_call_rate_confound.csv` |
| 07 | `07_allele_specific_m6a.py` | `07_asm_sites.csv`, `07_asm_donor_site.csv`, `07_asm_summary.csv` |

Numbering convention: script number = output number = execution order. `05` is the
single derivation point for the funnel; `06` reads `05_core_m6a_sites.csv` rather
than re-deriving, so figure and test can never disagree.

See `docs/04_asm_literature.md` for the published precedent behind each design
choice and for the claims that the literature does **not** support.

```bash
mamba activate mdd
Rscript 04_asm/01_find_drach_sites.R
python3 04_asm/02_verify_drach_reference.py
bash    04_asm/03_intersect_drach_genotypes.sh
python3 04_asm/05_plot_m6asnp_funnel.py
python3 04_asm/06_test_mdd_enrichment.py
python3 04_asm/07_allele_specific_m6a.py       # ~36 s

# stricter confidence gates
CONCORDANCE=0.95 MIN_ALLELE_READS=3 python3 04_asm/07_allele_specific_m6a.py
```

---

## 4. Biostrings motif scanning

`01` uses `vmatchPattern`. Two constraints discovered empirically:

- `matchPattern` **does not** interpret IUPAC ambiguity codes. Passing `"DRACH"`
  directly matches nothing, so the consensus is expanded to its 18 concrete
  sequences (`D`=A/G/T × `R`=A/G × `H`=A/C/T) and matched explicitly.
- `matchPattern` takes no `strand` argument in Biostrings 2.70.1. Minus-strand
  matching is done on `reverseComplement()` with offsets converted back. Motif
  offset 0 is the leftmost base for `+` motifs, the rightmost for `-`; the
  methylated adenine sits at transcript offset 2 on **both** strands, so
  `methyl_pos == genomic_start + 2` always.

The Biostrings implementation reproduces the validated output set-identically
(602,657 motifs; sorted files byte-identical to the original scan).

### `02` is deliberately independent

The QC gate re-derives every motif from the peak FASTA using plain Python string
handling — **not** Biostrings. If both the producer and the verifier used the same
library, the gate could no longer detect a library-level bug. It re-derives the 5-mer
and its reverse complement, checks IUPAC validity, and confirms `methyl_pos`.

Verified to fail correctly: injecting three corrupted motif records produced
`EXIT=1` and precise QC rows.

---

## 5. Two correctness traps found and fixed

### Non-unique motif identifiers

A peak sequence can contain the same `DRACH` expansion at several offsets —
**128,619 duplicate motif names** in `drach_sites.bed`. A single variant can therefore
disrupt one motif instance while preserving another (e.g. `5:67166241` sits at offset 0
of one `AGACA` and offset 4 of another). `motif_id` therefore includes the motif
instance start. The previous keying silently lost 2 real motif instances.

### A vacuous Fisher test

The RNA-seq VCFs record **only non-reference calls**, so there is no explicit `0/0`:

| genotype | count |
|---|---|
| `0/1` | 14,439 |
| `1/1` | 7,867 |
| `no_call` | 161,345 |
| `0/0` | **0** |

Restricting the Fisher denominator to genotyped donors collapses every contingency
table to "carriers vs zero non-carriers", returning p = 1 unconditionally and a fake
"zero nominally significant". The non-carrier group is now the **full cohort
complement** (20 MDD / 19 Ctrl), with `n_mdd_no_call` / `n_ctrl_no_call` reporting how
much of that complement is unobserved rather than confirmed reference.

---

## 6. Results

### Discovery funnel (derived in `05`)

Every level is recomputed from source data by `05_plot_m6asnp_funnel.py`. The script
asserts that each level is a strict subset of the one above it, so a logic error
fails loudly rather than silently inflating a level. **No recurrence threshold is
applied**: filtering on donor count before the mechanistic steps conditions the
m6A set on an arbitrary cutoff, so recurrence is reported as a distribution and
deferred (§6.2).

| # | Level | Count | Criterion applied |
|---|---|---|---|
| 1 | Donor × variant observations | 2,196,117 | PASS biallelic non-reference calls, 39 donors |
| 2 | Distinct expressed loci | 855,260 | deduplicated `chr:pos:ref:alt` |
| 3 | Within consensus m6A peaks | 50,161 | overlaps a TRESS IP-vs-Input region |
| 4 | Overlapping a DRACH motif | 7,243 | writer recognition sequence |
| 5 | Disrupts DRACH | 4,314 | abolishes writer competence at ≥1 instance |
| 6 | **Core m6A ablation** | **1,284** | abolishes the offset-2 catalytic adenine |

#### Two traps in this derivation

- **Peak containment is half-open.** `consensus_peaks.bed` is 0-based `[start,end)`.
  A 1-based SNP position `p` occupies 0-based `p-1`, so the test is
  `start < p <= end`. Treating the BED `end` as inclusive admits variants one base
  past every peak and inflates level 3 by ~150 loci (7,869 → 8,029 with a
  recurrence gate in place; 50,161 vs the correct value without it).
- **The overlap table is per donor.** `04_drach_variant_overlaps.csv` has one row
  per (donor, variant, motif), so a variant carried by several donors contributes
  several rows. All levels deduplicate to one row per (variant, motif instance)
  before counting.

Level 1 → 2 is a **unit change** (observations → loci), not a biological filter.
The −1,340,857 drop is deduplication and should not be read as variants removed.

### 6.2 Recurrence is deferred, not applied

64.8% of the 1,284 core sites are seen in exactly one donor; 9.0% are seen in 21–39.
The distribution is in `05_core_pool_recurrence.csv`. Any future threshold should be
chosen against that distribution rather than by default.

A threshold of ≥8 was previously inherited from the retired ASM path, where it was
the sample-size half of a paired criterion ("≥8 subjects, ≥2 MDD & ≥2 Ctrl") needed
to estimate a case/control effect. That rationale does not transfer here: the
m6A-SNP question is whether ablating a site reduces IP, which needs a handful of
`1/1` donors, not eight total carriers. Requiring eight carriers also discards the
rare variants most likely to have a large effect.

### 6.3 MDD enrichment — pooled test, null result (`06`)

Per-site selection on `n_mdd` was **rejected as circular**: `n_mdd` *is* the
case-control skew, so filtering sites by it and calling the survivors "MDD-exclusive"
approaches defining the answer. That script is archived. `06` instead asks whether
the core set as a whole is enriched, calibrated by permuting donor labels.

| Metric | Value |
|---|---|
| Core sites / carriers | 1,284 / 6,375 |
| Core pooled MDD fraction | 0.5272 |
| Donor base rate (20/39) | 0.5128 |
| **MDD share of all variant CALLS** | **0.5689** |
| **MDD/Ctrl per-donor call-rate ratio** | **1.254** (male 1.424, female 1.087) |
| Permutation statistic / z / p | +91.8 / +0.83 / **0.4174** |
| Background (non-peak) MDD fraction | 0.5691 |
| **Core − background z** | **−10.01** |

**No MDD enrichment in the core set** (permutation p = 0.42).

#### The coverage confound dominates the cohort

MDD donors are 51.3% of the cohort but contribute **56.9% of all variant calls**,
averaging 62,472 calls per donor against 49,825 for controls. That 1.25× ratio is
sequencing depth, not biology — and it is what produces an apparent 0.5691 MDD
fraction in the *background*. Singletons make this vivid: a variant seen in one
donor should be MDD with probability 20/39 = 0.513, yet the 642,830 singleton loci
average 0.635, falling monotonically to 0.513 at 21–39 carriers.

Consequences:

- Any case-control carrier analysis in this cohort is confounded by library depth
  until call counts are normalised (offset, coverage, or a matched call threshold).
- The core set is significantly **depleted** relative to background (z = −10.0),
  not enriched. Earlier in this project's history a z of +2.30 was read as a
  modest enrichment; that was the confound, and it does not survive permutation.

Sex strata (compare against both nulls):

| Stratum | Carriers | MDD fraction | Donor base rate | Call-rate base rate |
|---|---|---|---|---|
| male | 3,329 | 0.5305 | 0.5000 | 0.5875 |
| female | 3,046 | 0.5236 | 0.5263 | 0.5469 |
| both | 6,375 | 0.5272 | 0.5128 | 0.5689 |

Donor-label permutation is the correct null precisely because it preserves each
donor's total call count, so the coverage imbalance is carried into the null
distribution and cannot masquerade as signal.

### 6.4 Statistical ceiling

With 20 MDD / 19 Ctrl, **the best attainable per-site Fisher p is ≈ 0.020**, and no
configuration does better. `n_mdd` and `n_ctrl` rise together because common
variants are called in both arms, so no site reaches the 20 MDD / 0 Ctrl
configuration (p ≈ 1e-11) needed for significance after correction. Across all
1,284 core sites, **min BH-FDR = 1.0**.

Only 260 of the 1,284 core sites have any homozygous-alt carrier, and 200 have ≥2.
For the remaining 1,044 no dosage contrast exists at all.

**Level 6 is not "m6A sites."** Motif presence does not establish methylation —
direct-RNA work finds only ~5% of available DRACH motifs are methylated in vivo.
The set is adenosines in a DRACH context inside an IP-enriched region. See
`docs/04_asm_literature.md` §2.1.

**Resolution limits the achievable design.** A 1 bp ablation measured in 50 bp
bins inside ~200 bp MeRIP fragments at 300 bp median peak width is at the edge of
what the assay resolves; MeRIP cannot quantify methylation fraction at all. A
per-site IP-loss phenotype is not viable here — see `docs/04_asm_literature.md` §3.3.

---

## 7. Population frequency is a falsification control, **not** a selection criterion

MDD is polygenic: essentially all genuine risk loci are **common** — that is why
polygenic scores use common SNPs. Selecting by rarity therefore reduces the analysis
to single-variant rare-variant association testing, which in MDD also requires very
large N and solves nothing. It would also discard candidates for being lncRNAs rather
than for any biological reason. An earlier version of `07` ranked variants as
"COMMON (>10%) — not interpretable as risk"; that was wrong and was removed.

MAF is used only to supply the null expectation

$$\mathbb{E}[\text{carriers}] = n \times \big(1 - (1 - \text{MAF})^2\big)$$

and to test whether observed carriers exceed it. Among 28 T1/T2 candidates with
external data on 15:

| Verdict | n | Detail |
|---|---|---|
| enriched vs expectation | 4 | `DPP6`, `POU6F1`, `ASB8;PFKM`, `LINC01550` |
| as expected | 0 | — |
| depleted vs expectation | 11 | all MAF ≥ 0.12, incl. `RASGRF2` 0.16 |
| no population data | 13 | incl. `FAM13A` (gnomAD-only, no global freq) |

**obs/exp falls monotonically with MAF** (0.33 at MAF 0.27 → 0.89 at MAF 0.12), so the
recurrence ladder is largely a population-frequency artifact.

Observed counts sit *below* expectation because RNA-seq calling under-genotypes — 88%
of donor × locus cells are `no_call`. In-house carrier counts are therefore not a
substitute for external population frequency.

---

## 8. Environment quirks in this build

- **`rep()` / `rep.int()` cannot vectorize `times`.** `rep(1L, times=c(4L,1L))` fails
  with *"invalid 'times' argument"*. Use `Biostrings` `as.data.frame()` expansion,
  `data.table::fwrite`, or scalar-`times` `rep.int()`.
- **`IRanges` uses `start`/`width` slots**, not `start`/`end`, and `seqnames()` does not
  dispatch reliably. Access `obj@seqnames` and `obj@ranges@start` directly.
- **`data.table` `j` cannot see a column computed in the same call.** Aggregate, then
  merge in a second pass.
- **`DT[order(-col)]` evaluates `order()` outside the table scope** — bare column
  names must be qualified with `DT$col`.
- **`setorder()` rejects expressions** such as `is.na(x)`; add an explicit rank column.
- **`fread()` needs `R.utils`** to read `.gz` directly (added to `environment.yml`).
- **`msigdbr()` returns a tibble**, whose `unique()` dispatch loses `data.table`'s `.`
  NSE. Coerce with `as.data.table()` first.
- **`Content-Type` is not a valid R argument name** (hyphen) — use `Accept` only.

---

## 9. The core open question

Every analysis above asks a question N=39 cannot answer. The untested prediction is
**mechanistic and does not need disease labels**:

> At core-$m^6$A-ablating sites, IP signal should follow
> `IP(0/0) > IP(0/1) > IP(1/1)`.

A `1/1` carrier has no catalytic adenine on either allele, so the antibody should
recruit less RNA there. This is answerable with existing per-donor data
(`_data/peaks/<donor>/<donor>_BinsAndCounts.rda`, 4,309,757 50 bp bins with IP and Input
counts) and 261 testable sites.

Expected outcome is a **dose-response, not absence**: the antibody pulls down ~200 bp
fragments containing other m6A sites, and DRACH is the consensus rather than the only
methylated context.

Known limitations: 50 bp bins against a 1 bp event; consensus peaks are cohort-derived,
so a `1/1` donor contributed nothing to peak calling — a genuine selection effect and
the strongest objection to anticipate.