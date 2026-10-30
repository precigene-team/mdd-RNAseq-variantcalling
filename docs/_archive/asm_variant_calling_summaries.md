# Outline — Variant Calling: Purpose & Plan

*Course goal: understand what the `03_variant_calling/` step does and why, and be able to judge whether the ASM group-comparison plan is appropriate.*

**Lessons in order (we follow these; no skipping, no straying):**

1. **The biology setup** — DNA → RNA → protein; what m6A is; MeRIP-seq (IP vs Input); SNPs and heterozygous positions; why allele-level analysis is powerful; the two novel contributions.
2. **What variant calling does** — reads → alignment → GATK steps → VCF; why RNA-seq variant calling is special.
3. **From variants to allele counts** — ASEReadCounter: counting reads per allele in IP and Input fractions.
4. **Detecting ASM in one person** — M6Allele: is one allele more methylated than the other within a peak?
5. **The new group-comparison plan** — `11_asm_build_matrix.py` + `12_asm_compare.R`: the statistics and every judgment call.
6. **How to judge whether a plan is sound** — a reviewer's checklist.

**Method:** teacher asks one question at a time; student answers from intuition; every Q&A is recorded in `variant_calling_raw_review.md`; after each lesson a summary is appended to `variant_calling_summaries.md`.

# Lesson Summaries — Variant Calling: Purpose & Plan

*Summaries appended here after each completed lesson.*

---

## Lesson 1 — The biology setup

**Big picture:** We study RNA, not DNA, because DNA is the fixed blueprint (same in every cell, never changes), while RNA shows what a cell is *actually doing* at that moment — which genes are switched on, which are suppressed.

**Key ideas (built up via Q&A):**

1. **Two copies of every gene.** One from each parent; usually identical, sometimes differing by one letter. If the two copies differ, that spot is a **SNP**; seeing both letters in the reads means the person is **heterozygous**. The comparison that matters is *inside one person* — tagging each RNA read as "copy A" or "copy B" and counting them — not comparing against other people or population frequency.
2. **Two-tube experiment (MeRIP-seq).** From the same brain RNA we make an **Input** tube (all RNA) and an **IP** tube (only methylated RNA, fished out by an antibody against the m⁶A "stamp"). Counting reads per allele in each tube lets us see whether the methylated fraction is biased toward one copy.
3. **ASM in one line.** If Input is balanced (50/50) but IP is skewed (80/20), one allele is **preferentially** methylated — that's Allele-Specific Methylation (ASM).
4. **Why the Input is the ruler.** The Input ratio is the *baseline*: if one copy is simply transcribed more, Input is already skewed — so a skewed IP would be boring. The interesting signal is the IP skew *beyond* what Input predicts. Deviation is measured relative to Input.

**Misconceptions corrected:**

- "Methylation is *only* on one allele" → it can be *partial* (80/20, 70/30); the word is **preferential**.
- "Is the variant unique in the population?" → that's allele-frequency-in-a-population, a separate question; here the ruler is the person's own two copies.

**Terms introduced:** DNA vs RNA roles; SNP; heterozygous / homozygous; Input vs IP fractions; m⁶A stamp; preferential vs exclusive; ASM.

**So what (course goal):** The whole variant-calling step exists to tag reads by allele and compare IP vs Input ratios — and our group analysis will do this per-peak across the 39 subjects.

---

## Lesson 2 — What variant calling does

**Big picture:** From millions of short RNA reads, variant calling finds the positions where a person's letters differ from the reference — by mapping reads to the reference and trusting only positions backed by enough reads.

**Key ideas (built up via Q&A):**

1. **Mapping.** Each read is placed at its position on a single "standard text" (the reference genome). Mismatches are allowed; that's how differences are found. The reference is a patchwork of a few anonymous donors — not a SNP database, and not one ethnicity's genome.
2. **Depth beats noise.** One odd letter is plausibly a sequencing error; a deviation seen in many reads is a real variant. That's the whole art of calling.
3. **RNA coverage is patchy — two consequences.** (a) A silent gene produces zero reads → a **hole** in the output → that position can never be tested for ASM in that person. (b) Reads that span splice junctions contain an "N" in their CIGAR string; `SplitNCigarReads` cuts them so base-counting happens inside the correct exon. Same root cause: reads come from transcripts, not the whole genome.
4. **Call in Input, not IP.** Beyond IP being sparser, IP is *biased*: it is enriched for methylated molecules, so its allele ratio is already tipped. Discovering variants in IP would make the discovered positions depend on methylation — circular, and it would inflate ASM. Input gives an unbiased list of the person's letters; the same positions are then counted in both tubes.

**Misconceptions corrected:**

- "The reference genome contains the common SNP variants" → it does not; common variants live in separate databases (dbSNP/gnomAD) used only for annotation afterward.
- Ethnicity/gender/lifestyle cannot *mask* letters — the letters live in the reads; those factors affect read depth (and sex enters the model as a covariate later), never the alleles.
- "IP is avoided only because it has less data" → the disqualifying reason is bias, not just sparsity.
- "The 'hole' is related to SplitNCigarReads" → they are separate consequences of the same patchy-RNA cause.

**Terms introduced:** reference genome; alignment/mapping; BAM; CIGAR string, "N" (splice gap); SplitNCigarReads; HaplotypeCaller; VCF; depth as error-defense; exon; discovery vs counting.

**So what (course goal):** The pipeline's variant calling (Steps 1–4, Input only) produces the unbiased list of heterozygous positions per person. Next: counting reads per allele at those positions — Lesson 3.

---

## Lesson 3 — From variants to allele counts

**Big picture:** The VCF is the "catalog" of a person's heterozygous positions, but the actual per-allele *counts* must come from counting reads again, in both Input and IP, with the same rules.

**Key ideas (built up via Q&A):**

1. **VCF gives "which positions"; BAM gives "how many."** The VCF was built from Input only, so it holds no IP counts; ASM needs counts in both fractions. Even its Input depth is the caller's filtered view, not a neutral count.
2. **Reference-allele (read-mapping) bias.** At a heterozygous A/G SNP, the allele matching the reference (A) aligns more cleanly; the differing allele's reads carry a mismatch and can fail or mis-map. So a raw A:G ratio is tilted toward the reference allele even with no biological skew.
3. **Why Input is the null.** Because of that bias, a fixed 0.5 is not the right baseline. The person's own Input allelic fraction is the natural ruler: ASM is IP skew *relative to* Input, and Chromosome-relative correction absorbs the reference-mapping artifact (Castel et al. 2015).

**Misconceptions corrected:**

- "We count 'to be sure the variant is real'" → the VCF already called it; we count for the ratio, and because IP data only exists in the reads.
- The raw A:G ratio read off a BAM is NOT automatically clean 50/50 — it carries a systematic reference-allele tilt.

**Terms introduced:** per-allele count; reference allele vs alternate allele; major-mapping (reference) bias; Input allelic fraction (Input-as-null); ASEReadCounter.

**So what (course goal):** This is the precise, bias-aware read count that the per-peak ASM model consumes. Next: combining individual noisy SNPs within an m6A peak into a per-person "is one allele methylated more" call — Lesson 4.

---

## Lesson 4 — Detecting ASM in one person

**Big picture:** For one person, decide whether one allele is preferentially methylated within an m6A peak — by pooling all the peak's SNP counts and comparing IP vs Input.

**Key ideas (built up via Q&A):**

1. **Pool, don't pick.** A single SNP has ~10–15 reads (a coin flip). Summing counts across all heterozygous SNPs inside the peak gives one combined IP ratio with real weight (15+20+25 → 60 reads).
2. **Compare the peak's IP ratio to the same person's Input ratio** — not to a fixed 0.5 (Lesson 3's reference-bias point). Deviation = preferential methylation.
3. **Gap ≠ significance.** The ratio gap is the *estimate*. Whether it's "significant" is decided by a hypothesis test (e.g., Fisher's exact, beta-binomial) that weighs the gap against the read depth behind it: 0.85-vs-0.52 with ~10 reads is noise; with 10,000 reads it's overwhelming.
4. **One person ≠ the project's answer.** A striking single-person skew is evidence of ASM *in that person only*. The MDD-vs-Control question needs many people — Lesson 5.

**Misconceptions corrected:**

- FDR and machine learning are NOT involved inside one person's one peak; BH-FDR is the *group* multiple-testing step (many loci × subjects), and ML is not used anywhere.

**Terms introduced:** pooling/summing counts per peak; IP vs Input per-person ratio; estimate vs significance (test); per-person ASM "call" vs group question.

**So what (course goal):** We now know how to get a per-person per-peak ASM signal. Next: testing whether that signal biologically differs between MDD and Control across the cohort — Lesson 5 (the actual analysis plan).

---

## Lesson 5 — The new group-comparison plan

**Big picture:** Turn every subject × peak into one comparable number, then test whether MDD's number is higher than Control's — while absorbing the two things that could fake a difference.

**Key ideas (built up via Q&A):**

1. **One number per person: their skew = IP fraction relative to their own Input fraction** — not the raw IP. Each person's Input is their personal ruler (Lesson 3), so it folds out the reference/mapping baseline.
2. **Two artifacts the test must absorb:**
   - Reference-mapping bias → placed as an offset (`logit(Input AF)`), so the null is per-person, not 0.5.
   - Sex imbalance → a covariate (`condition + sex`) so the MDD-vs-Ctrl estimate is sex-adjusted; within-sex splits kept exploratory (the cohort is 10M/10F MDD vs 10M/9F Ctrl — unbalanced).
3. **Reliability lives in the read pile.** A skew of 0.30 from 2,000 reads vastly out-weights the same 0.30 from 15. So the plan fits a **beta-binomial regression on the raw read counts** (each person weighted by their evidence) instead of a plain test on the skew fractions treated as equally certain.
4. Per-peak test → correct across peaks with **BH-FDR**; checkpoint a sensitivity sweep over thresholds (the "judgment calls" to review).

**Misconceptions corrected:**

- Compare skew-relative-to-own-Input, not raw IP fractions, across people.
- "Sex imbalance" is a confounder that must be a covariable, not ignored ("sample-type imbalance").
- The big MDD-vs-Ctrl number is not the raw IP gap but the sex-adjusted, read-weighted difference in per-person skew.

**Terms introduced:** per-person skew; confounding; covariate/sex-adjustment; offset (per-person Input as null); read-weighted / beta-binomial regression; BH-FDR across peaks; sensitivity sweep.

**So what (course goal):** The student can now read the plan (`11_asm_build_matrix.py` → `12_asm_compare.R`) and see *why* each design choice exists. Next: the reviewer's checklist — is this plan sound? — Lesson 6.

---

## Lesson 6 — How to judge whether a plan is sound

**Big picture:** Flip from builder to reviewer and inspect the judgment calls and the model's honesty — not the biology, but the choices that could make results fragile or misleading.

**Key ideas (built up via Q&A):**

1. **Cutoff choices are bias/variance knobs.** `MIN_IP_READS`/`SKEW_CUT`/`MIN_SUBJECTS`: strict → fewer, more robust, but dropped sparse/real signals (false negatives); lax → more peaks, but one noisy person can manufacture a hit (false positives). There's no free "best" number.
2. **Distinguish two cutoffs students conflate:** `SKEW_CUT` = signal *magnitude* within one peak/person; `MIN_SUBJECTS` = cross-person *reproducibility*. "Should recur in all cohorts" is robustness, which belongs to min-subjects/FDR — not the skew cut.
3. **The sweep is the hedge.** `12_asm_compare.R` sweeps `MIN_SUBJECTS ∈ {6,8,12}` and `SKEW ∈ {0.1,0.2,0.3}`. A finding that survives the sweep is trustworthy; one that flips with the cutoff is fragile. That's what a good plan does instead of betting on one number.
4. **Model honesty / engine fallback.** glmmTMB and lme4 give each SUBJECT its own baseline (per-person ruler); the bare `glm` averages everyone onto ONE shared baseline, so genuine per-person ruler variation gets misread as signal. The report must say which engine ran — the output's `engine` column does exactly that, so a `glm` row is read with caution.

**Misconceptions corrected:**

- "Coherent/sparse" vs really a `glmmTMB-vs-glm` (per-subject ruler) distinction.
- The fancier engine gives per-*subject*, not per-*cohort*, baselines.
- A single hard-coded cutoff is a review liability; the versatility of a sweep is the defense.

**Terms introduced:** bias/variance knob; sensitivity sweep; model `engine` field; per-subject random effect; fallback vs primary estimator.

**So what (course goal):** The student can now *judge* the ASM group-comparison plan in `12_asm_compare.R` and `docs/asm_analysis_methodology.md`: it's sound because discovery is Input-only (blind), reference-mapping bias is absorbed as per-person Input offset, sex is a covariate, reliability is count-weighted (beta-binomial), and cutoffs / engine are reported and swept — so no single reviewer's objection (min-subjects harsh, glm fallback, unbalanced sex) can silently decide the answer.

---

## Course complete (Lessons 1–6)

The student progressed from biology to a defensible reviewer of the entire `03_variant_calling/` ASM plan:

1. Biology setup (RNA vs DNA; m6A; two copies/SNPs; IP vs Input; ASM).
2. Variant calling (mapping, depth-vs-error, patchy RNA, Input-only discovery).
3. Allele counts (VCF=which positions, BAM=how many; reference-mapping bias).
4. Per-person ASM (pool SNPs per peak; Input-as-null; estimate vs test).
5. Group plan (per-person skew; offset + sex covariate; read-weighted beta-binomial; BH-FDR; sweep).
6. Judgment (cutoffs as a two-sided knife; sweep hedge; engine `engine` honesty).

**Skill tailoring applied throughout:** concrete examples, never pre-empt answers, visualize before asking, simplest metaphor on "don't get it" — all logged in the skill's Tailoring Log.
