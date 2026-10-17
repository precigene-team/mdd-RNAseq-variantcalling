# Project Analysis Roadmap & TODO List

This document lists the completed and remaining tasks for the MDD RNA-seq and MeRIP-seq epitranscriptomic analysis pipeline.

---

## 📋 Project Status Checklist

### 📈 1. Extended DGE & Interaction Contrasts

* [x] Modify/create R scripts to perform sex-specific and interaction-based contrasts (relative methylation change: `Condition:Fraction`).
* [x] Implement filter thresholds ($|\text{log2FC}| \ge 2$, $p < 0.05$).
* [x] Generate sex-stratified Volcano and MA plots.
* **Status**: **Completed**. Results saved in `results/02_dge_analysis/tables/` and `results/02_dge_analysis/figures/`.

### 🏷️ 2. Downstream Functional Annotations

* [x] Map DEG/DMG gene IDs to HGNC symbols, biotypes, and genomic coordinates.
* [x] Resolve and map subgroup SRA accession identifiers to cohort groups.
* **Status**: **Completed**. Mapping tables saved to `results/02_dge_analysis/tables/01_gene_mapper.csv` and `results/02_dge_analysis/tables/06_subgroup_sra_accessions.txt`.

### 🛡️ 3. Quality Control Registry (QC Resolutions)

* [x] **QC-01: Postmortem Confounders**: Calibrate and apply the 12M exonic counts gate (gated 36 subjects, filtered out `male_ctrl_1`).
* [x] **QC-02: Cell-Type Heterogeneity**: Run cell-type deconvolution (marker gene profiles) to verify glia/astrocyte loss (confirmed astrocyte reduction in MDD, $p = 0.0389$).
* [x] **QC-03: Gene-level vs. Peak-level Resolution**: Run MACS2 peak-calling on the 36 gated subjects and perform a cohort-wide overlap audit to validate gene body counts concordance.
* [x] **QC-04: Antibody occupancy limitations**: Document and interpret methylation shifts as relative interaction ratios only.
* **Status**: **Completed**. Documented in [cohort_limitations_and_resolutions.md](docs/cohort_limitations_and_resolutions.md).

### 🧬 4. Variant Calling (GATK Best Practices)

* [x] Write GATK processing scripts under `03_variant_calling/` (ReadGroups, MarkDuplicates, SplitNCigarReads, HaplotypeCaller, VariantFiltration).
* [x] Execute GATK pipeline (Steps 01-05) for all 78 BAMs on HPC via Slurm array.
* [x] Step 06 ASE ReadCounter completed for all 39 subjects.
* **Status**: **Completed**.

### 📊 5. Centralized MultiQC Aggregation

* [x] Create the centralized configuration [multiqc_config.yaml](multiqc_config.yaml) to organize reports chronologically.
* [x] Run the final MultiQC execution command to consolidate raw/trimmed FastQC, STAR alignment, Qualimap, and featureCounts metrics.
* **Status**: **Completed**. Final report generated at `_data/multiqc/rna_seq_pipeline_summary.html`.

### 🛠️ 6. Pipeline Refinements & Reorganization

* [x] Move cell-type deconvolution script from `scratch/cell_type_deconvolution.py` to `utils/cell_type_deconvolution.py` and redirect output scores to `results/utils/tables/cell_type_mgp_scores.csv`.
* [FAILED/DELEGATED_TO_DMR] Reaudit DGE analysis (specifically inspect `05_ma_male_meth.png` for potential cloud tilting or normalization biases).
* [ ] Refactor `02_dge_analysis/08_enrichment_plotting.R` to show divided facet headers for Activated/Suppressed pathways.
* [x] Format and clean up the TODO list.

### 📖 7. Learning & Literature Review

* **Bulk RNA-seq Mentorship (Active)**:
  * [x] Module 1: File Formats (FASTQ vs. FASTA vs. SAM/BAM) (Completed)
  * [x] Module 2: Normalization Math (TPM, FPKM, and DESeq2) (Completed)
  * [x] Module 3: Confounders & Cell-Type Deconvolution (Completed)
* **RIP-seq/m6A Epitranscriptomics Mentorship (Future)**:
  * [ ] Read RIP-seq/m6A review paper
    * [ ] Document the paper's DOI/URL link for tracking
    * [ ] Focus on the "Bulk RNA-seq vs. RIP-seq" transition:
      * Identify how RIP-seq normalization handles IP (immunoprecipitated) samples vs. Input/IgG controls, compared to standard Bulk RNA-seq normalization
      * Note why peak calling (finding local sequence enrichment) is required rather than gene-level expression quantification
    * [ ] Detail the metagene profile differences (3' UTR/stop codon enrichment in MeRIP-seq vs. uniform coverage in RNA-seq)
    * [ ] Establish benchmarks: QC (FRiP, SNR), Peak Caller logic (MACS2 vs. RIPSeeker/Piranha), and Biological Validation (DRACH motif)
    * [ ] Review m6A biology fundamentals (Writers, Erasers, Readers)
    * [ ] Extract and log key methodologies/comparison tables verbatim

---

## 🧬 ACTIVE PIVOT — $m^6A$-SNP Sequence-Level Analysis (Oct 2026)

**Status:** Rebuilding `04_asm/` around the locked pivot. Read-ratio ASM is **deprecated**.
**Execution order = file number** (`01` → `07`), matching `01_core_alignment` / `02_dge_analysis`.

### Locked Decisions

| Decision | Value |
|---|---|
| Analysis unit | Consensus m6A peak + exact DRACH motif, **strand-aware** |
| Read-ratio ASM (M6Allele / ASEReadCounter / glmmTMB) | **Dropped** — window-pooling lumping artifact in outbred humans |
| Genotype classes | `0/0` intact · `0/1` allele-specific loss · `1/1` biallelic ablation |
| Het-only restriction | **Removed** — `0/1 ∪ 1/1` are both valid $m^6A$-SNPs |
| Recurrence threshold in the funnel | **Removed entirely.** Filtering on donor count before the mechanistic steps conditions the m6A set on an arbitrary cutoff and discards the rare variants most likely to have a large effect. Reported as a distribution instead. |
| Per-site case-control selection (`n_mdd > 6`, `n_ctrl <= 2`) | **Rejected as circular** — `n_mdd` *is* the skew. Archived. Replaced by the pooled test. |
| Case-control test design | **Pooled + donor-label permutation within sex.** Preserves each donor's call count, so the 1.25× MDD/Ctrl coverage ratio cannot masquerade as signal. |
| Null calibration | **Permutation, not nominal p** — m6A-QTL precedent (Zhang 2020 *Nat Genet*) |
| Confidence criterion (future ASM) | **Cross-sample concordance >90%**, not a p-value threshold (Cao 2023 precedent) |
| Level 6 label | **"DRACH-core adenine"**, not "m6A site" — motif presence ≠ methylation (~5% of DRACH motifs are methylated in vivo) |
| Cohort composition | **20 MDD / 19 Ctrl** (39 donors; 20 male, 19 female) |
| Motif matching | **Biostrings `vmatchPattern`** — DRACH expanded to 18 concrete sequences |
| Strand conversion | **Biostrings `reverseComplement` / `complement`** |
| Validation roadmap | PGC GWAS coloc → PsychENCODE replication → clinical sub-phenotype |

### Pipeline (run in numeric order)

| Step | Script | Output(s) |
|---|---|---|
| 01 | `01_find_drach_sites.R` | `_data/peaks/drach_sites.bed` *(602,657 motifs)* |
| 02 | `02_verify_drach_reference.py` | `02_drach_reference_qc.csv` + `02_drach_reference_failures.csv` *(only on failure)* |
| 03 | `03_intersect_drach_genotypes.sh` | delegates to 02 then 04 |
| 04 | `04_annotate_strand_aware_disruption.R` | `04_drach_variant_overlaps.csv`, `04_drach_overlap_summary.csv` |
| 05 | `05_plot_m6asnp_funnel.py` | `05_funnel_levels.csv`, `05_core_m6a_sites.csv`, `05_core_pool_recurrence.csv`, `05_m6asnp_funnel.pdf/.png` |
| 06 | `06_test_mdd_enrichment.py` | `06_mdd_enrichment_summary.csv`, `06_mdd_enrichment_strata.csv`, `06_call_rate_confound.csv` |
| 07 | `07_allele_specific_m6a.py` | `07_asm_sites.csv`, `07_asm_donor_site.csv`, `07_asm_summary.csv` |

**Convention:** script number == output number == execution order. `05` is the
single derivation point for the funnel; `06` reads `05_core_m6a_sites.csv` rather
than re-deriving, so figure and test cannot disagree.

**Archived to `_legacy/04_asm/`** (retired, superseded, circular, or blocked):
`05_m6asnp_dosage_matrix.py`, `06_m6asnp_recurrence.R` (recurrence gate rejected as
cherry-picking), `07_validate_external.R` (hosts unreachable), `08_rank_by_mechanism.R`,
`06_select_mdd_exclusive_core_m6a.py` (per-site selection on the outcome statistic).

> **Tracking gap:** `.gitignore` excludes `_legacy/`, and the 11 files moved out of
> `04_asm/` were force-added into it. These 5 later-archived scripts, plus
> `_legacy/04_asm/README.md`, exist on disk but are **not in the repository** — a
> clone gets the earlier batch only. Recorded in `docs/04_asm_postmortem.md` §2.4.

### Run commands

```bash
cd /home/mashxp/Documents/Bioinformatics/mdd_RNAseq_pipeline
micromamba activate mdd

Rscript 04_asm/01_find_drach_sites.R            # ~8 s
python3 04_asm/02_verify_drach_reference.py      # ~40 s, hard gate
bash    04_asm/03_intersect_drach_genotypes.sh   # runs 02 + 04 internally
python3 04_asm/05_plot_m6asnp_funnel.py          # ~3 s
python3 04_asm/06_test_mdd_enrichment.py         # ~9 s (20k permutations)
python3 04_asm/07_allele_specific_m6a.py         # ~36 s

N_PERM=100000 python3 04_asm/06_test_mdd_enrichment.py   # tighter null
CONCORDANCE=0.95 MIN_ALLELE_READS=3 python3 04_asm/07_allele_specific_m6a.py
```

### Phase 0–2 — COMPLETE

* [x] Pull `consensus_peaks.fa` from hpc-01 (18,831,795 B, 33,123 seqs); FASTA↔BED 0 mismatches
* [x] 3 GB genome FASTA confirmed **unnecessary** — all motifs lie inside peaks
* [x] `01` rewritten in R/Biostrings; output **set-identical** to the validated Python + HPC baseline
      (602,657 motifs, 49 chromosomes, 302,562 fwd / 300,095 rev, sorted files byte-identical)
* [x] `02` QC gate passes 602,657/602,657 motifs; **negative control verified** (3 injected
      corruptions → `EXIT=1`, precise QC rows)
* [x] `03` no longer uses `-g het` → 2,196,117 donor×variant rows across all 39 donors
* [x] `04` strand-aware complementation via Biostrings `complement()`
* [x] Fixed **`motif_name` non-uniqueness** (128,619 duplicates): one variant can disrupt one
      motif instance and preserve another. `motif_id` now includes the instance start
* [x] Legacy read-ratio scripts + results moved to `_legacy/04_asm/` (see its local-only
      `README.md` — not tracked, see the tracking-gap note above)
* [x] Renumbered to contiguous `01`–`07`; the SLURM orchestrator retired

### Current Results

| Metric | Value |
|---|---|
| DRACH motifs in consensus peaks | 602,657 |
| Distinct variants overlapping a motif | 7,392 |
| Distinct variant × motif-instance overlaps | 7,968 |
| **DISRUPTS_DRACH** | **4,709** (2,370 `+` / 2,339 `−`) |
| PRESERVES_DRACH | 3,259 |
| DRACH-core ablations (offset 2), all loci | 1,404 |

> **2,339 minus-strand disruptions** were previously inverted by the old annotator, which
> compared genomic alleles against transcript IUPAC sets without complementing.

### Phase 3 — Funnel (COMPLETE)

Six levels, each a strict subset of the one above, all derived in `05` with an
assertion that fails loudly rather than silently inflating a level.

| # | Level | Count |
|---|---|---|
| 1 | Donor × variant observations | 2,196,117 |
| 2 | Distinct expressed loci | 855,260 |
| 3 | Within consensus m6A peaks | 50,161 |
| 4 | Overlapping a DRACH motif | 7,243 |
| 5 | Disrupts DRACH | 4,314 |
| 6 | **DRACH-core adenine (offset 2)** | **1,284** |

**Level 6 is NOT "m6A sites."** Motif presence does not establish methylation:
direct-RNA work finds only ~5% of available DRACH motifs are methylated in vivo.
These are adenosines in a DRACH context inside an IP-enriched region. Renamed
accordingly — see `docs/04_asm_literature.md` §2.1.

Two correctness traps, both documented in `docs/04_asm_methods.md` §6:
- **Peak containment is half-open.** BED is 0-based `[start,end)`; a 1-based SNP
  position `p` sits at 0-based `p-1`, so the test is `start < p <= end`. Treating
  the BED `end` as inclusive admits variants one base past every peak.
- **The overlap table is per donor.** All levels deduplicate to one row per
  (variant, motif instance) before counting.

Recurrence is **deferred, not applied** — 64.8% of the 1,284 core sites are
singletons and 9.0% appear in 21–39 donors (`05_core_pool_recurrence.csv`).
A threshold of ≥8 was inherited from the retired ASM path, where it was the
sample-size half of a paired "≥8 subjects, ≥2 MDD & ≥2 Ctrl" criterion. That
rationale does not transfer: this question needs a handful of `1/1` donors, not
eight total carriers, and requiring eight discards the rare variants most likely
to have a large effect.

### Phase 4 — MDD enrichment test (COMPLETE, NULL RESULT)

`06_test_mdd_enrichment.py`. Pooled, not per-site, because per-site selection on
`n_mdd` is circular — `n_mdd` *is* the case-control skew.

| Metric | Value |
|---|---|
| Core sites / carriers | 1,284 / 6,375 |
| Core pooled MDD fraction | 0.5272 |
| Donor base rate (20/39) | 0.5128 |
| **MDD share of variant CALLS** | **0.5689** |
| **MDD/Ctrl per-donor call-rate ratio** | **1.254** (male 1.424, female 1.087) |
| Permutation statistic / z / p | +91.8 / +0.83 / **0.4174** |
| Background (non-peak) MDD fraction | 0.5691 |
| **Core − background z** | **−10.01** |

**No MDD enrichment.** Once donor labels are permuted within sex — which
preserves each donor's call count and therefore the coverage imbalance — the core
set shows no case-control skew (p = 0.42).

**The coverage confound dominates everything.** MDD donors are 51.3% of the
cohort but 56.9% of all variant calls, because MDD libraries yield more calls
(mean 62,472 vs 49,825). That 1.25× ratio — not biology — produces an apparent
enrichment of 0.5691 in the *background*, and the core set is significantly
*depleted* relative to it (z = −10.0). Any case-control carrier analysis in this
cohort is confounded by library depth until this is normalised.

### Phase 5 — Allele-specific m6A, corrected design (COMPLETE, POSITIVE)

`07_allele_specific_m6a.py`. The legacy read-ratio ASM approach was retired for a
window-pooling artifact (~600 bp windows pooling unrelated private SNPs from
different donors). The fix is structural: compare the two alleles **within each
heterozygous donor**, IP vs that donor's own Input. This is self-normalising —
degradation and depth affect both alleles equally, and allele-specific mapping
bias cancels because Input carries the same penalty.

Data was already local: `_data/ase/` has per-variant `refCount`/`altCount` for
39 IP + 39 Input files. No new HPC work.

| Metric | Value |
|---|---|
| Donors | 39 (20 MDD / 19 Ctrl) |
| Donor-site cells tested (≥2 reads/allele both fractions) | 571,937 |
| Distinct sites | 157,693 |
| Cells surviving \|Δ\| ≥ 0.1 and exact Fisher | 92,866 |
| Cells at Fisher p < 0.05 | 75,767 |
| Sites with ≥2 significant donors | 14,940 |
| **Sites at concordance > 90%** | **8,711** |
| …in an m6A peak | 3,355 |
| …overlapping a DRACH motif | 544 |
| …**disrupting DRACH** | **317** |
| …**core m6A ablating (offset 2)** | **109** |

Median \|Δ\| in the core set = 0.237 (Cao 2023 reported ~0.48 at their top
subset, same order of magnitude).

#### Evidence the signal is biological, not artifact

Allele imbalance could easily be a mapping or depth artifact, so it was tested
against the null expectation of peak localisation:

| min reads/allele | cells | in-peak fraction | enrichment | Fisher p |
|---|---|---|---|---|
| 2 | 75,767 | 0.346 | 1.00× | 1.0 |
| 5 | 44,040 | 0.393 | 1.13× | 5e-58 |
| 10 | 26,872 | 0.449 | 1.30× | 2e-195 |
| 20 | 14,589 | 0.503 | 1.45× | 6e-274 |

**Enrichment strengthens with depth.** Low-depth sampling noise would shrink or
vanish under a stricter depth floor; this grows to 1.45× (background in-peak
fraction is 0.191). That is the signature of real signal.

Also note: 33,632 cells favour the alt allele vs 42,135 favour ref — roughly
symmetric, as expected if the direction is driven by which allele carries the
methylation rather than by a systematic bias toward one.

#### Independent re-discovery

Two loci previously verified by a separate sequence-level audit re-appear at
100% concordance:

- **`HAX1` chr1:154273441 T>C** — verified core ablation earlier (offset 2, minus strand)
- **chr19:16815256 A>T** — was in the earlier MDD-exclusive shortlist

Neither was given to the ASM pipeline, so this is an independent recovery.

#### Honest caveats

- **Median support is thin.** Median donors per concordant site is 2–3. Only
  1,853 sites have ≥4 donors and 752 have ≥6. In the core set, 29 of 109 have ≥4
  donors. Concordance on 2 donors is weak evidence; a `MIN_DONORS` gate is the
  obvious next tightening.
- **No FDR correction was applied across 157,693 sites.** The depth-stratified
  peak enrichment is the honest significance argument, not the per-site p-value.
- **Two known loci were not recovered** (`RASGRF2` chr5:81226098, `POU6F1`
  chr12:51189658) — both were `1/1`-only sites in the earlier list, i.e. no
  heterozygotes, so the per-donor design structurally cannot test them. That is
  expected, not a failure.
- **Not yet validated experimentally.** Cao 2023 confirmed top hits with reporter
  assays in HEK293T. That is the natural follow-up and it is what makes this a
  discovery rather than an association.

#### Data hazard found and fixed

`contig` dtype differs between the IP and Input ASE files (string vs mixed
int/str). A naive merge silently drops **5–72% of shared sites per donor**
(`female_ctrl_1`: 58,835 shared → 16,632 joined). Contig is cast to string on
both sides before joining. Any earlier analysis touching these files is suspect.

#### Runtime

36 s. An earlier version used a pure-Python hypergeometric loop; on deep cells
(`totalCount` up to 29,301) the support range is tens of thousands of terms, so
that version timed out. Replaced with scipy's C implementation plus Cao's own
effect-size threshold (|Δ| ≥ 0.1) applied before the exact test — faster and
closer to the reference method.

### Known Limitations (state these before a reviewer does)

1. **Resolution.** 1 bp perturbation, 50 bp counting bins, ~200 bp MeRIP
   fragments, 300 bp median peak. MeRIP cannot quantify methylation fraction, and
   detectable m6A change is already at the assay's edge
   ([10.1038/s41598-020-63355-3]). A per-site IP-loss phenotype is not viable.
2. **No single-base methylation measurement.** No miCLIP/PA-m6A-seq/SELECT here, so
   "is this adenine actually methylated?" is unanswerable. ~5% of DRACH motifs are.
3. **Peak ascertainment circularity.** Peaks were called by TRESS using IP-vs-Input,
   so selecting inside them and then testing IP enrichment re-uses the same signal.
4. **Dosage power.** Only 260 of 1,284 core sites have any homozygous-alt carrier;
   200 have ≥2. For 1,044 there is no dosage contrast at all.
5. **Coverage confound.** See Phase 4. Not yet normalised.
6. **No external replication.** GWAS Catalog returns 410/404; m6AVar, RMVar and
   m6A-TSHub hosts refuse connection. m6ASNP (GigaScience 2018) predicts 133,394
   gain + 214,884 loss variants and would give an independent comparison set.
7. **Postmortem confounders.** Sex, age, PMI, RIN not modelled. Sex is stratified
   in `06` as a start; the others need checking against the sample sheet.

### Environment Note

`bioconductor-biostrings`, `r-data.table`, `r-r.utils` declared in `environment.yml`.
Pipeline runs in the `mdd` env: `micromamba run -n mdd python 04_asm/05_plot_m6asnp_funnel.py`

### Reference Docs

* `docs/04_asm_methods.md` — module 04 methods, traps, results
* `docs/04_asm_literature.md` — **methodological literature**, what each design choice
  is justified against, and which claims are unsupported
* `docs/04_asm_postmortem.md` — **retrospective**: why read-ratio ASM was abandoned, the
  rejected analyses and why, the `contig` dtype hazard, the coverage confound, plus a
  defensibility brief and the open items carried forward
* `docs/thesis_defense_narrative_log.md` — defense narrative; §3.3 records the
  **retracted** per-site 8-vs-2 carrier claims

## 🧬 Active Workstream (July 2026)

### Step 07 — M6Allele ASM (SUPERSEDED)

The M6Allele Java step never completed, but it proved unnecessary: the ASE counts
it would have produced already exist per-variant in `_data/ase/`, and the
per-donor analysis that replaced it does not need a Java tool. See Phase 5.

* [x] Fix: bcftools filter to biallelic het SNPs (crashed on multi-allelic)
* [x] Fix: graceful `if` wrap for Java crashes (IndexOutOfBounds, NumberFormatException)
* [x] Fix: auto-scale `-Xmx` from `SLURM_MEM_PER_NODE` (not hardcoded)
* [x] Fix: `--mem=16G` → `32G` for java heap headroom
* [~] abandoned in favour of `07_allele_specific_m6a.py`

### Steps 08-09 — DRACH Annotation

* [x] DRACH sites identified: 602,657 sites from 33,123 peaks
* [x] SNP×DRACH overlap: 211,284 total, 122,959 motif-disrupting
* [x] DRACH CSV pulled locally: `results/03_variant_calling/tables/09_drach_disrupting_snps.csv` (22 MB, 211K rows)

### TRESS DMR Analysis

* [x] Male DMR: 8,508 significant DMRs (completed)
* [x] Female DMR: 3,669 significant DMRs (completed)
* [x] `lowcount=15` post-hoc filter script: `utils/filter_dmr_lowcount.sh`
* [x] Filtered outputs: `results/utils/tables/10_dmr_{male,female}_filtered.csv` (8,495 / 3,590 DMRs)
* [x] QC-05 documented: missing RIN/age covariates inflate DMR counts vs paper (239/696)

### Downstream Integration (for presentation)

* [ ] Overlap m6A-SNPs (DRACH CSV) with DMR regions (male/female)
* [ ] Write integration script / summary table
* [ ] Prepare figures for mid-internship presentation

---

## 📋 Paper Replication Tracker (Mitsuhashi et al. 2024, iScience)

### DMR Parameter Alignment

| Parameter         | Paper                 | Ours                          | Match?                                  |
| ----------------- | --------------------- | ----------------------------- | --------------------------------------- |
| Tool              | TRESS                 | TRESS                         | ✅                                       |
| Binsize           | 50                    | 50                            | ✅                                       |
| lowcount          | 15                    | default (~10)                 | ⚠️ post-hoc filtered (negligible effect) |
| Covariates        | RIN + age             | None                          | ❌ not publicly available                |
| Threshold         | p<0.05, FC>1.2        | p<0.05, FC>1.2                | ✅                                       |
| Within-sex male   | 239 DMRs              | 8,495 (our unfiltered: 8,508) | ❌ covariate-driven gap                  |
| Within-sex female | 696 DMRs              | 3,590 (our unfiltered: 3,669) | ❌ covariate-driven gap                  |
| Interaction       | 629 DMRs (152M, 477F) | —                             | ❓ not run yet                           |

### Replication Steps

* [ ] Understand paper's DEG pipeline (DESeq2 + Subread, see STAR Methods)
* [ ] Run paper's sex-disease interaction DMR contrast
* [ ] Compare our DMR gene lists to paper's Supplementary Tables S3/S5
* [ ] GO enrichment on DMR genes (clusterProfiler, ViSEAGO)
* [ ] PPI network of shared DMR-DEG clusters (STRING + Cytoscape)
* [ ] RRHO analysis (male vs female methylation rank concordance)

---

# SNP Downstream

* [x] Re-check metadata, ensure each subject has correct Input/m6A-IP pair.
* [x] Use STAR-aligned BAMs from both Input and IP. (Steps 01-02 on HPC)
* [x] Run MarkDuplicates. (Step 02)
* [x] Run SplitNCigarReads separately for each BAM: (Step 03)
  * [x] All Input samples.
  * [x] All m6A-IP samples.
* [x] Verify split BAMs: no corruption, .bai index exists, correct read group and sample names.
* [x] Use only Input RNA-seq for SNP calling: (Steps 04-05)
  * [x] Run HaplotypeCaller for each Input sample. (per-sample GVCF mode — Step 04)
  * [x] Joint genotyping across all Input samples. (CombineGVCFs 04b + GenotypeGVCFs 04c)
  * [x] Keep only biallelic SNPs. (SelectVariants --restrict-alleles-to BIALLELIC — Step 05)
  * [x] Filter low-quality SNPs. (VariantFiltration FS > 30.0, QD < 2.0, -f PASS — Step 05)
  * [x] Extract heterozygous 0/1 SNPs per individual. (bcftools view -g het — Step 05)
  * [x] Prioritize SNPs with DP ≥ 10, GQ ≥ 20, and sufficient reads supporting both REF and ALT. (bcftools view -i 'FMT/DP >= 10 && FMT/GQ >= 20' — Step 05, 2.25M het sites across 39 subjects)
* [x] Call m6A peaks by comparing m6A-IP with Input:
  * [x] Each IP sample must be paired with the correct Input from the same individual.
  * [x] Used exomePeak2 per-subject.
  * [x] Export peak list in BED format.
  * [x] Create consensus peaks (`consensus_peaks.bed`, 33,123 peaks).
* [x] Intersect Input heterozygous SNP list with m6A peak regions: (Steps 08-09)
  * [x] Keep only het SNPs located within or near m6A-enriched regions.
  * [x] Do not call SNPs independently from IP data. (note: also annotated DRACH motifs — 211K SNP×motif overlaps, `09_drach_disrupting_snps.csv`)
* [x] At selected SNPs, count alleles in both BAMs from the same individual: (Step 06 — ASEReadCounter)
  * [x] REF count and ALT count in Input.
  * [x] REF count and ALT count in IP.
  * [x] Used GATK ASEReadCounter.
* [x] Calculate per-SNP metrics:
  * [x] REF fraction = REF / (REF + ALT).
  * [x] ALT fraction can be calculated similarly.
  * [x] Compare REF/ALT ratios between IP and Input for the same individual.
* [~] Only consider as allele-specific m6A candidate when: (Step 07 — M6Allele, in progress)
  * [ ] SNP is a reliable heterozygous call in Input.
  * [ ] Sufficient coverage in both Input and IP.
  * [ ] Allele ratio in IP is significantly skewed relative to Input.
  * [ ] SNP lies within an m6A peak region.
  * [ ] Not dominated by mapping bias, strand bias, or a handful of outlier reads.
* [ ] After obtaining per-subject ASm6A:
  * [ ] Compare allele imbalance levels between MDD and control groups.
  * [ ] Statistical model must adjust for at least sex.
  * [ ] Use ~ sex + condition; to test sex-specific differences, use ~ sex * condition.
* [ ] Interpretation notes:
  * [ ] A SNP lying in an m6A peak does NOT mean that SNP is the methylated nucleotide.
  * [ ] Only conclude that the allele carrying the SNP has a different enrichment level in the m6A region compared to the other allele.
  * [ ] Do not take a SNP that appears only in IP and conclude it is an m6A-associated SNP.

---

## 📝 External Internship Report Tracker (`_report/`)

### Front matter

* [x] Cover page
  * [x] University/Faculty name (VNU-HCM — Ho Chi Minh City University of Technology — Faculty of Chemical Engineering)
  * [x] Report title: EXTERNAL INTERNSHIP — Precision Gene Joint Stock Company
  * [x] Student: Dinh Chan Phong — Student ID: 2352902
  * [x] Company supervisor: Dr. Luu Phuc Loi
  * [x] Academic supervisor: Dr. Phan Thi Huyen
  * [x] HCMC, July 2026
* [ ] Inner title page
* [x] Acknowledgements
* [ ] Evaluation by the internship unit (company + supervisor, signature, stamp)
* [ ] Table of contents
* [ ] List of abbreviations
* [ ] List of tables
* [ ] List of figures
* [ ] Report content
* [ ] References
* [ ] Appendix (if any)

### Chapter 1. Overview of the internship unit

* [ ] History of establishment and development, location, functions and responsibilities
  * [/] Development history/milestones — asking the company (startup, likely none)
* [ ] Personnel management and operations (organizational chart, hiring process, number/qualifications of staff)
  * [x] Hiring process — via direct contact with Dr. Luu Phuc Loi
  * [/] Employment system, staff count/qualifications — asking the company
* [x] Raw materials and products/services (properties, specifications, substitutability, supply availability; products/services, quality requirements, market demand)
  * [x] Raw materials/inputs (raw sequencing data) — §1.6
  * [x] Products/services — §1.4
  * [x] Quality requirements and market demand of products/services — §1.4.4–1.4.5
* [ ] Site layout diagram — user to provide layout later

### Chapter 2. Production (experiment) process and equipment

* [ ] Block diagram of the production/experiment process
  * [ ] Receiving, inspection, and storage of raw materials
  * [ ] Finished product inspection
  * [ ] Each step of the process: purpose and requirements
* [ ] Main equipment: basic construction, operating principle, function, operation
* [ ] Troubleshooting / incidents and remedies
* [ ] Actual experimental results obtained (if any)

### Chapter 3. Auxiliary facilities

* [ ] Power supply, water supply, boiler, compressed air systems (if any)
* [ ] Environmental treatment systems: drainage, water treatment, waste treatment, exhaust gas treatment

### Chapter 4. Quality control

* [ ] Quality management tools at the internship site (ISO, GMP, HACCP, ...)

### Chapter 5. Engineer duties

* [ ] Understand the duties/responsibilities of an engineer at the internship unit
* [ ] Relationship of the engineer with other fields at the internship unit

### Chapter 6. Conclusion

* [ ] Summary of internship results
* [ ] Lessons learned / recommendations (if any)

### Back matter

* [ ] References (ordered by first citation; book/journal names in italics)
* [ ] Appendix (if any)

### Report figures/tables to insert (Sections 2.5–2.8)

* [ ] **Fig — Pipeline schematic** for §2.5 (variant-oriented analysis flow: ReadGroups → MarkDuplicates → SplitNCigarReads → HaplotypeCaller GVCF → GenotypeGVCFs → filter). No artifact exists — need to create/draw.
* [x] **Table 2.4 — QC thresholds** for §2.5.9 (FS ≤ 30, QD ≥ 2.0, DP ≥ 10, GQ ≥ 20) — text already contains it; matches `05_variant_filtration.sh`.
* [x] **§2.6.2 TRESS peak summary table** — built via `03_variant_calling/13_make_tress_peak_summary.py` from `asm_matrix.tsv`: `results/03_variant_calling/tables/tress_peak_summary.csv` (per-cohort) + `tress_peak_per_subject.csv` (per-subject). 15,351 loci / 8,793 genes / 39 subjects.
* [x] **Fig. 2.x — ASM filtering funnel** for §2.6.8 — render the mermaid in `docs/asm_filter_flow.md` (118,245 → 85,453 → 16,445 → 312 → 305).
* [x] **Table — DRACH-disrupting SNPs** for §2.6.4 — **RERUN + FIXED** (`09_intersect_drach_snps.sh` now `-g het`): per-subject het counts are now distinct (507–1,349 sites/subject), and the reproducible cohort summary is auto-written to `results/03_variant_calling/tables/09_drach_summary.csv` (11,542 distinct het sites / 11,578 motif instances / 12,405 distinct overlaps / 6,975 disrupting / 5,432 preserving / 2 both). Table inserted as **Table 2.5** in `_report/C2.5.md` §2.6.4. ✅ Old "211,284 / 122,959" figures superseded.
* [x] **Table 2.6 — 5 significant ASM loci** for §2.7.4 — **regenerated from `asm_results.rds`** (was stale/empty): `results/03_variant_calling/tables/12_asm_primary_significant.csv` now has 5 rows (MGAT5, BIN3, IER2, TENM3, MTCH1, fdr<0.05).
* [x] **Fig. — Volcano plot** for §2.7.4 — `results/03_variant_calling/figures/12_asm_volcano.pdf`.
* [ ] **Fig. — Sensitivity barplot** for §2.7.4 (or Supplementary) — `results/03_variant_calling/figures/12_asm_sensitivity.pdf`.
* [x] (optional) **Table — Sensitivity grid** — **regenerated** from `asm_results.rds`: `results/03_variant_calling/tables/12_asm_sensitivity.csv` (9 rows, non-zero; primary sig at 0.2/8 = 5).
* [ ] Decide primary vs. Supplementary placement (recommend: volcano + significant-loci table → main; sensitivity → Supplementary).

### Formatting / Layout

* [ ] Font: Times New Roman, size 13, line spacing 1.5
* [ ] Margins: top 2 cm, bottom 2 cm, right 2 cm, left 2.5 cm
* [ ] Page numbers: bottom-right corner
  * [ ] Acknowledgements → end of list of figures: numbered i, ii, iii, iv, v, ...
  * [ ] Chapter 1 → end of references: numbered 1, 2, 3, ...
  * [ ] Appendix: no page numbers
* [ ] Subsections numbered by chapter (max 4 levels: 1 → 1.1 → 1.1.1 → 1.1.1.1)
* [ ] Figures and tables numbered by chapter (e.g., Figure 3.1, Table 3.2)

---

## 🧬 Candidate Loci Functional Annotation (Mentor Task)

* [ ] Multi-layer functional annotation for 5 significant ASM candidate loci (**MGAT5**, **BIN3**, **IER2**, **TENM3**, **MTCH1**):
  * [x] **Layer 1 (Ensembl VEP)**: Map coordinates, variant consequences (coding/intronic/UTR), and affected transcripts
  * [x] **Layer 2 (dbSNP)**: Retrieve rsIDs, minor allele frequencies (MAF), and variant info
  * [ ] **Layer 3 (GTEx Portal)**: Audit tissue-specific eQTL and sQTL associations (prioritize Frontal Cortex BA9 / Brain)
  * [ ] **Layer 4 (GWAS Catalog)**: Cross-reference associations with MDD, psychiatric traits, and neurological phenotypes
  * [ ] **Layer 5 (Literature Curation)**: Query PubMed/Google Scholar (`[Gene] + MDD/depression/brain/stress/m6A`)
  * [ ] **Deliverable**: Compile summary table (`Gene | SNP/rsID | Variant Consequence | eQTL/sQTL | Brain/MDD Evidence | m6A Evidence | Reference`)
