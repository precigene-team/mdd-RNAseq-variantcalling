# 04_asm Postmortem — From Read-Ratio ASM to Sequence-Level m6A-SNPs

**Date:** 17 October 2026
**Scope:** Module `04_asm` (epitranscriptomic allele-specific methylation), Oct 2026 pivot.
**Status:** Work paused and landed. Pipeline is complete and reproducible; the
MDD case-control question returned a **null**, and that null is the finding.

Two audiences, two halves. **Part 1** is an engineering retrospective for the
mentor — what was tried, why it was abandoned, what the traps were. **Part 2** is a
defensibility brief for the thesis committee — what a hostile reviewer will attack,
and the honest answer to each.

---

## Part 1 — Engineering retrospective

### 1.1 The arc

| Phase | Approach | Outcome |
|---|---|---|
| Initial | M6Allele + GATK ASEReadCounter + glmmTMB beta-binomial, ~600 bp windows | **Failed.** Null / non-replicating |
| Audit | Why did it fail? | Window-pooling lumping artifact; MWU replication FDR = 1.00 |
| Pivot | Sequence-level DRACH motif destruction ($m^6A$-SNPs), exact base-pair | Mechanistically clean, assay-limited |
| Test | Pooled MDD-vs-Ctrl, donor-label permutation within sex | **Null.** p = 0.42 |
| ASM v2 | Per-donor IP-vs-Input allele imbalance, concordance > 90% | **Positive but thin.** 8,711 sites, median 2-3 donors |

The pivot is the substance of this module. The read-ratio pipeline is not a
preliminary attempt that got improved — it was **structurally incapable** of
answering the question in this cohort, and finding that out is what the module
actually produced.

### 1.2 Why window-level ASM could not work here

M6Allele and peers assign variants to ~600 bp peak windows and test group skew in
allele ratios within each window. In an inbred animal model a window's variants
belong to one animal, so pooling is harmless. In 39 outbred humans a window
collects *unrelated private variants from different donors*, each present in one
person. Pooling them manufactures a group-level skew out of singletons.

Concretely: the `BIN3` core mutation was carried by **1 donor**. Non-parametric
Mann-Whitney U replication gave **FDR = 1.00**. The signal did not survive a
change of test.

This was not a tuning failure. No threshold on a pooled statistic fixes an
artifact that is in the pooling itself.

### 1.3 Two analyses rejected on principle, not on results

Both are archived in `_legacy/04_asm/` and both are worth recording, because the
temptation to reinstate them will recur.

**Per-site case-control selection** (`06_select_mdd_exclusive_core_m6a.py`).
Selecting sites because their MDD carrier count is high, then reporting that
count as evidence of MDD association, is circular — `n_mdd` *is* the skew under
test. There is no held-out split, so the ranking cannot fail; it only restates the
selection rule. This is worse than an underpowered test, because it looks
significant.

**Recurrence threshold** (`06_m6asnp_recurrence.R`). The inherited "≥8 carriers"
cutoff came from the retired ASM path, where it was half of a paired
"≥8 subjects, ≥2 MDD & ≥2 Ctrl" criterion. It does not transfer: this question
needs a handful of `1/1` donors, not eight total carriers, and the cutoff discards
the rare variants most likely to have a large effect. 64.8% of the 1,284 core
sites are singletons — the cutoff would have deleted most of the signal it was
meant to find. Recurrence is now *reported as a distribution*
(`05_core_pool_recurrence.csv`), not applied as a filter.

### 1.4 The data hazard that silently corrupted joins

`contig` dtype differs between the IP and Input ASE files — string on one side,
mixed int/str on the other. A naive `merge` on `(contig, position)` silently
drops **5-72% of shared sites per donor**. Worst case `female_ctrl_1`: 58,835
shared sites became 16,632 after the join — a 72% loss that raises no error.

Fix: cast `contig` to string on both sides before joining.

> **Standing caution:** any earlier analysis touching these ASE files is suspect
> until re-checked against this. Silent row loss on an inner join is the failure
> mode that survives review, because every downstream number still looks
> plausible.

### 1.5 The coverage confound (why the case-control test is null)

MDD donors are **20/39 (51.3%)** of the cohort but **56.9% of all variant calls**,
because MDD libraries yield more calls (mean 62,472 vs 49,825 per donor).

That 1.25x ratio alone produces a background MDD fraction of **0.5691** —
indistinguishable from any "enriched" candidate list. Measured against that
background, the core set is significantly **depleted** (z = **−10.0**).

The pooled permutation test (p = 0.42) is null because it is the *correct* test:
it permutes donor labels **within sex**, preserving each donor's call count, so
the depth imbalance cannot masquerade as signal. A test that did not preserve call
counts would have reported a confident positive.

**Consequence:** no case-control carrier analysis in this cohort is interpretable
until library depth is normalised. This is not fixable by a covariate after the
fact in a 39-donor design.

### 1.6 Two correctness traps worth remembering

- **Peak containment is half-open.** BED is 0-based `[start, end)`; a 1-based SNP
  at position `p` sits at 0-based `p-1`, so the test is `start < p <= end`.
  Treating BED `end` as inclusive admits a variant one base past every peak.
- **The overlap table is per donor.** All funnel levels deduplicate to one row per
  (variant, motif instance) *before* counting, or level counts inflate.

And one bug found and fixed during this cleanup, recorded because it nearly
shipped:

- **`07_asm_summary.csv` mixed two denominators.** The concordance-filtered
  headline (8,711 sites) sat directly above annotation counts (5,056 / 455 / 151)
  computed over *all* sites with ≥2 significant donors, with nothing labelling the
  difference. It read as a funnel; it was not. `TODO.md` quoted the
  concordance-filtered set (3,355 / 544 / 317 / 109) — the correct numbers. Fixed
  by emitting both denominators under explicit `_conc_passed` and
  `_all_ge2_donor_sites` suffixes. Two documents disagreed about the same
  quantity; a machine-readable file would have caught it immediately.

### 1.7 Cohort composition — verified

From `_data/_metadata/sample_sheet.csv` (77 library rows → 39 donors):

| | MDD | Ctrl | total |
|---|---|---|---|
| male | 10 | 10 | 20 |
| female | 10 | **9** | 19 |
| total | **20** | **19** | 39 |

Note the **female control stratum is short one donor (n=9)**. Since `06` stratifies
by sex, the female permutation runs on 10 vs 9 — thin enough that the
within-sex null should be treated as indicative, not exact. This was previously
stated the wrong way round (19 MDD / 20 Ctrl) in the defense narrative; corrected.

---

## Part 2 — Defensibility brief

State these before a reviewer does.

### 2.1 What can be claimed

- **Allele-specific methylation is detectable and peaks are enriched.** The
  strongest claim in the module. Depth-stratified in-peak fraction rises
  monotonically with depth — 1.00x at ≥2 reads/allele to **1.45x at ≥20**
  (background 0.191, p ≈ 6e-274). Sampling noise would shrink under a stricter
  depth floor; growth is the signature of real signal.
- **DRACH motif disruption is a sequence-level fact.** Abolishing the motif
  removes METTL3/14 writer competence on that allele. This needs no
  case-control comparison and is not in question.
- **The pipeline is correct and reproducible.** `01` output is set-identical to
  the validated baseline (602,657 motifs). `02` is a hard gate with a verified
  negative control. `07` re-runs to identical counts.

### 2.2 What cannot be claimed

- **No MDD enrichment.** p = 0.42. The core set is depleted vs background
  (z = −10.0). Do not present a case-control carrier list.
- **No causal claim.** N = 39 cannot support it. The mechanism is plausible and
  sequence-verified; the disease link is untested.
- **These are not "m6A sites."** They are adenosines in a DRACH context inside an
  IP-enriched region. Direct-RNA work finds only ~5% of available DRACH motifs are
  methylated in vivo. The level-6 label is "DRACH-core adenine" for this reason.

### 2.3 The four attacks, and the answers

**"You tested 157,693 sites with no FDR correction."**
Correct, and it is why the per-site p-value is *not* the significance argument.
The significance argument is the depth-stratified peak enrichment (§2.1), which is
a single test with a monotone dose-response. Nominal per-site p-values are
reported for ranking only. This is a defensible position provided it is stated
plainly rather than discovered by the committee.

**"Median support is 2-3 donors."**
The weakest point. Only 1,853 sites have ≥4 donors, 752 have ≥6; in the core set
29 of 109 have ≥4. Concordance on 2 donors is weak evidence. A `MIN_DONORS` gate is
the obvious next tightening — `MIN_DONORS` is already a documented knob, it is
simply not yet set.

**"You selected inside peaks, then tested peak enrichment."**
Peak ascertainment circularity, and it is real. Peaks were called by TRESS using
IP-vs-Input, so testing IP enrichment inside them re-uses the same signal. The
depth-stratification mitigates but does not eliminate this. m6A-QTL precedent uses
input-matched peak definitions for exactly this reason.

**"Your assay cannot resolve single bases."**
It cannot, and no computational fix exists. ~200 bp MeRIP fragments, 50 bp
counting bins, ~300 bp median peaks, 1 bp perturbation. Detectable m6A change is
already at the assay's edge (doi:10.1038/s41598-020-63355-3). A per-site
methylation-fraction phenotype is not a viable target.

### 2.4 Open items carried forward

Deliberately **not** fixed in this cleanup:

1. **`.gitignore` excludes `_legacy/`** while 11 files were force-added into it.
   They are tracked-but-ignored, so any *new* archive file silently never appears
   in `git status`. Five later-archived scripts (`05_m6asnp_dosage_matrix.py`,
   `06_m6asnp_recurrence.R`, `06_select_mdd_exclusive_core_m6a.py`,
   `07_validate_external.R`, `08_rank_by_mechanism.R`) plus
   `_legacy/04_asm/README.md` exist on disk but are **not in the repository**. A
   clone gets the first batch only.
2. **No result is reproducible from a clone.** `results/` is git-ignored, so every
   headline number here lives in prose only. The summary tables total < 10 KB and
   would be worth tracking; the heavy ones (`07_asm_donor_site.csv` 81 MB,
   `04_drach_variant_overlaps.csv` 7.3 MB) should stay ignored.
3. **`r-r.utils` is an unused dependency.** It was required by the retired
   `09_plot_m6asnp_funnel.R`; that step is now `05_plot_m6asnp_funnel.py`.
   Comment in `environment.yml` marks it safe to drop.
4. **Two known loci were not recovered** by the ASM design — `RASGRF2`
   chr5:81226098 and `POU6F1` chr12:51189658. Both were `1/1`-only, i.e. no
   heterozygotes, so a per-donor design **structurally cannot** test them.
   Expected, not a miss. Say so before it is raised.
5. **Postmortem confounders unmodelled.** Sex is stratified in `06`; age, PMI and
   RIN are not yet checked against the sample sheet. The female-ctrl n=9 imbalance
   (§1.7) compounds this.
6. **No external replication.** m6AVar, RMVar and m6A-TSHub hosts refuse
   connection from this machine; GWAS Catalog returns 410/404. m6ASNP (GigaScience
   2018) predicts 133,394 gain + 214,884 loss variants and would give an
   independent comparison set — **not yet attempted**.

---

## Appendix B. A retracted finding: the "degenerate rows"

Worth recording because it is the same failure mode as §1.4, and because it briefly
propagated into code and documentation.

**The claim.** `07_functional_annotation.csv` was reported to contain 1,942 rows
(3.1%) with byte-identical carrier counts (`n_subjects=8`, `n_mdd=4`, `n_ctrl=4`),
described as placeholder padding written when per-gene listings were truncated.
The rows concentrate in the largest brain-expressed genes, which is what made the
fabrication story plausible. A filter was added to exclude them.

**Why it was wrong.** The pattern was inferred from a partial view — the tail of
the file, where low-carrier rows had clustered — without checking the
distribution. Three checks refuted it:

- `n_mdd + n_ctrl == n_subjects` in all 61,685 rows, so `n_subjects` is exactly
  the carrier count.
- `n_subjects` has min **8**, median 15, max 39. 8 is a floor, not a cap, and
  55,224 rows exceed it.
- Carrier counts decay smoothly (8:6,461 · 9:5,506 · 10:4,816 · 11:4,052) with no
  spike at 8. Fabrication would produce a spike.
- `(4,4)` is the modal group split, among near neighbours: `(5,3)` 1,804,
  `(5,4)` 1,710, `(6,4)` 1,328.

**Cost.** 254 of 11,493 coding variants (2.2%) and 8 of 433 panel variants were
deleted. 87 genes dropped out of the analysis entirely, having been reduced to
zero retained carriers. No conclusion changed — the burden results are null either
way — but legitimate measurements were discarded on a premise that did not hold.
Removed 18 October 2026.

**The real constraint.** `06_build_transcriptome_snp_table.py` bins loci by
recurrence, and this table contains only the `High Recurrence (>=8 Subjects)`
tier. Private and low-recurrence variants are therefore absent **by construction**,
which is a firmer explanation for the missing rare variants than RNA-seq depth
alone.

**The lesson, same as §1.4.** Both failures are the same shape: a pattern noticed
in a slice of the data was generalised to the whole dataset without checking the
marginal distribution. Row-count error (§1.4) and frequency-floor error (here) are
both caught by looking at the distribution of a single column before acting on it.
The difference from §1.4 is that one was a code bug and this one was an
interpretation error wearing a code bug as a disguise.

Full detail: `docs/polygenic_burden_methodology.md` §4.

## Reproducing

```bash
micromamba activate mdd
Rscript   04_asm/01_find_drach_sites.R            # ~8 s
python3   04_asm/02_verify_drach_reference.py      # ~40 s, hard gate
bash      04_asm/03_intersect_drach_genotypes.sh   # delegates 02 + 04
python3   04_asm/05_plot_m6asnp_funnel.py          # ~3 s
python3   04_asm/06_test_mdd_enrichment.py         # ~9 s
python3   04_asm/07_allele_specific_m6a.py         # ~36 s
```

Outputs → `results/04_asm/tables/` and `results/04_asm/figures/` (git-ignored).

## References

- `docs/04_asm_methods.md` — module methods, correctness traps, full results
- `docs/04_asm_literature.md` — what each design choice is justified against
- `TODO.md` — locked decisions and per-phase status
- Zhang 2020 *Nat Genet* — m6A-QTL permutation calibration precedent
- Cao 2023 — >90% cross-donor concordance criterion