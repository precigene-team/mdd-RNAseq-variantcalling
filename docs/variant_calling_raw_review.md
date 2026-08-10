# Raw Review — Variant Calling: Purpose & Plan

*Course started: 2026-08-02. Method: teacher asks, student answers from intuition, every exchange recorded verbatim. Summaries live in `variant_calling_summaries.md`.*

---

## Lesson 1 — The biology setup

### Q1 (Lesson 1)
You're studying RNA from 39 donated brains — never DNA. If you wanted to know what a cell was "doing" at the moment of death, which would be a better clue: reading the person's DNA (the fixed blueprint every cell shares) or reading the RNA (the working copies currently being made)? Why?

> Answer (student): "To know what cells are doing. RNA is the key since RNA displays the current state of that person abundance in transcriptomic products activated genes and suppressed genes. While DNA can't tell whether certain genes are activated or suppressed"

Teacher's next move: Confirmed the intuition (DNA = static blueprint / potential; RNA = current activity). Introduced that every cell holds two copies of each gene (one from each parent) and that the two copies can differ by one letter. Q2 builds on this toward SNPs.

### Q2 (Lesson 1)
Every cell carries two copies of each gene — one inherited from your mother, one from your father. Most of the time the two copies match letter-for-letter, but occasionally they differ at a single letter. Now imagine a person whose two copies of a gene differ at one letter — and both copies are active in the same cell. What could you compare? And what would that comparison be able to tell you that looking at the gene "on average" could not?

> Answer (student): "We would infer the difference in 1, in the person's replicated reads itself? Or 2 if the variant is unique among the population compare to the other samples? And we could determine whether the SNP is heterozygous or homozygous?"

Teacher's next move: Affirmed the correct parts (we can tag reads by which copy they came from; we can tell heterozygous vs homozygous from the reads). Redirected the "population uniqueness" tangent — that's about how common a variant is in a population, not what this per-person analysis needs. Q3 introduces the two-library (Input vs IP) setup to lead to the allele-skew question.

### Q3 (Lesson 1)
For each person we make two tubes from the same brain RNA: Tube 1 "Input" = all the RNA; Tube 2 "IP" = only the methylated RNA, fished out with an antibody that grabs the m⁶A "stamps." Suppose we count, in each tube, how many reads came from copy A vs copy B. Input = 50/50 balanced, but IP = 80/20 one-sided. What would you conclude about the two copies and their m⁶A stamps? And why is the Input tube the important "ruler" to compare against?

> Answer (student): "Tube 2 is the most important one, if the methylated section got pulled down in different ratio this suggested that the section is only on one of the alleles. And to determine the 'deviation' we must use the Input to get the baseline amount of reads?"

Teacher's next move: Affirmed strongly — Input as baseline is exactly the right intuition. Refined "only on one allele" → "preferentially on one allele" (skew can be partial). Added the crucial nuance: the Input ruler rules out the boring explanation (if copy A is simply transcribed more, Input is already skewed, so IP skew beyond Input skew is what matters — deviation measured relative to Input). **Lesson 1 complete.**

---

## Lesson 2 — What variant calling does

### Q1 (Lesson 2)
Picture the raw data: millions of short RNA reads per person, just strings of letters, and you want to count copy A vs copy B at one specific letter-position in the genome. Before any counting is possible, what's the first thing you have to figure out for each read — where it came from, or how to spell-check it? And if a read's letters were changed by noise during sequencing, how would you even know which position it belongs to?

> Answer (student): "We need the aid of RNA seq to read out the transcripts then map it to a consensus genome. The consensus genome database also include the common SNP variants. And while mapping (BAM) since the amount of reads is abundant the algorithmn can detect unique changes: If the deviation amount is abundant enough to make statistical difference, then it might be a variant, else it's just a sequencing error"

Teacher's next move: Affirmed the two correct pillars (1) reads are mapped/aligned to a reference genome; (2) depth is the defense against sequencing error — enough supporting reads = variant, too few = error. Corrected one misconception: the reference genome is a single representative sequence and does NOT contain the common SNP variants as alternatives; common variants live in separate annotation databases (dbSNP/gnomAD) used afterward. Q2 introduces why RNA-seq variant calling is special (reads come only from expressed genes).

> Student question (clarifying): "So for the 'consensus' genome, I'm not sure if the sample ethnicity, gender or lifestyle, age or many factors could potentially mask SNP or not?"

> Teacher's answer: The reference is a single "standard text" — a composite of a small number of anonymous donors, so it is no single ethnicity's genome. A person's real letters are always carried in their own reads and show up as differences; population only changes which differences are common, and our analysis never uses population frequency. Sex does not mask letters — it changes copy count (males 1 X, females 2 X), and sex enters our model later as a covariate, not at mapping. Age and lifestyle change RNA amounts, not DNA letters, so they cannot mask a SNP. Mild truth: if the reference is very different from a person in some region, reads there can map less well — but our ASM ruler is the person's own two copies, so both alleles still get read and compared regardless of which one matches the reference. Noted as a genuinely good question; rewarded and reassured.

(Continuing after the student's interjection, back to Q2: why RNA-seq variant calling is special — silent genes produce no reads.)

> Student request: "Can you give me concrete examples? Please update my skill?" → Skill `socratic-teaching` updated (Tailoring Log 2026-08-03: always give concrete worked examples after any abstract concept, without waiting to be asked).

> Student correction: "Wait i haven't answered Q2 at all. Can you ask again?" — Teacher had overstepped by giving worked examples that pre-empted the student's answer. Logged in skill Tailoring Log (2026-08-03): never answer the question you just asked; worked examples may illustrate *context*, not supply the answer. Q2 re-posed, answer pending.

> Answer (student) to Q2: "It would look like a hole. And that position wouldn't be tested for ASM."

> Student follow-up question: "Is this related to SplitNCigar Bams?"

Teacher's next move: Confirmed the hole answer. Explained that the "hole" and SplitNCigarReads are two separate facts about RNA-seq variant calling: (a) silent gene → zero reads → hole (expression), unrelated to read splitting; (b) SplitNCigarReads handles reads that DO exist but span an exon-exon junction (CIGAR "N"), splitting them so base-counting is done within the right exon. Both facts share one root cause: RNA coverage is patchy (only expressed exons). Noted that HaplotypeCaller in this pipeline runs on Input only.

### Q3 (Lesson 2)
Your pipeline calls variants on the Input tube only — HaplotypeCaller is never run on the IP tube. The IP tube is the same person's RNA, just filtered for the m⁶A-stamped molecules. Why do you think we would discover variants in Input but never in IP?

> Answer (student): "The IP has even less data than the Input it only consist of the m6a regions."

Teacher's next move: Affirmed the sparsity point (IP has fewer reads → noisier). Added the deeper, disqualifying reason — BIAS, not just quantity: IP is, by construction, enriched for the m⁶A-stamped molecules, and if one allele is preferentially methylated the IP ratio is already tipped (Lesson 1's 80/20). Calling variants in IP would make the *set of discovered positions* depend on methylation skew — circular, and it would inflate ASM. Discovery in Input = an unbiased list of the person's letters; the same positions are then counted in BOTH tubes. **Lesson 2 complete.**

> Student request: "Can you give me some examples for L2 Q3? Then wait" → Teacher supplied worked examples (see chat: sparsity; homozygous-collapse at extreme ASM; inflation of the discoverable set). Waiting for the student's own answer to Q3 before advancing.

> Answer (student, after examples): "The input tubes presents a more coherent read mappings more than the IP ones which are only enriched in m6a regions and discarded the non m6a regions"

Teacher's next move: The student again framed it as data volume / mapping coherence (IP "discards non-m6A"), which is the sparsity side — still not the disqualifying BIAS/circularity point. Re-framed: IP is a FILTER; calling in IP makes the discoverable catalog depend on the very methylation signal we measure. Re-asked to elicit the bias point in the student's own words before closing Lesson 2.

### Q4 (Lesson 2)
A heterozygous position where the A copy is heavily stamped (IP runs ~99/1) — if we built the catalog from IP, what would the variant caller most likely decide about that one SNP? And why is that the worst possible way to build a catalog?

> Answer (student): "It would discard that as an SNP as a sequencing error?"

Teacher's next move: Confirmed — at ~99/1 the caller drops the SNP (looks not-heterozygous / like an error), and that position is the strongest ASM signal. The student articulated the bias/circularity point, which closes Lesson 2. The takeaway: discovery must be blind to the signal; Input is that blind ruler. **Lesson 2 complete (final).**

---

## Lesson 3 — From variants to allele counts

### Q1 (Lesson 3)
After HaplotypeCaller we hold a VCF — a list of heterozygous positions, each with a genotype like A/G, plus per-site fields (read depth). We also still hold the original BAM with every raw read. To measure ASM we need per-allele counts (how many A-reads, how many G-reads) at each position, in both Input and IP. Should those counts come from the numbers baked into the VCF, or from going back to the raw reads in the BAM? What's risky about trusting the VCF's own numbers?

> Answer (student): "We need both? To ensure the variant in the VCF is surely a variant"

Teacher's next move: Affirmed we go back to the reads (don't trust the VCF's arithmetic). Refined the *reason*: it's not to re-prove the variant (the VCF already called it). (a) The VCF came from Input only — it holds NO IP allele counts; ASM needs both fractions, so the IP numbers can only come from the BAM. (b) Even for Input, the VCF's depth is the caller's filtered output with known reference-allele bias. So: VCF = "which positions"; BAM = "how many of each allele."

### Q2 (Lesson 3)
A heterozygous SNP where A matches the reference genome and G differs from it, at the same position. When counting aligned reads, which allele is counted slightly more reliably / aligns back more easily — the reference-matching A or the differing G? If that bias is real, what does it do to the raw A:G ratio you read off the BAM?

> Answer (student): "The A is more skewed?"

Teacher's next move: Confirmed — the reference-matching allele A is favored (its reads align cleanly; the differing allele's reads carry a mismatch and can fail or mis-map). So the raw A:G ratio is tilted toward the reference, even with no biological skew. This reference/mapping bias is why ASEReadCounter is used and why the Input allelic fraction (not a fixed 0.5) is the per-subject null in the model (Castel et al. 2015). **Lesson 3 complete.**

---

## Lesson 4 — Detecting ASM in one person

### Q1 (Lesson 4)
We have, for one person, per-allele counts at single SNP positions inside an m6A peak. A single SNP may hold only ~15 reads — a coin flip. A peak usually has several heterozygous SNPs. If you had to judge whether this person's peak is "skewed toward one allele," how would you combine the noisy SNP counts — trust one loud SNP, or lump them? And what would you compare the combined ratio against?

> Answer (student): "It would account for the FDR and certain machine learning algorithmns right?"

Teacher's next move: Re-framed — FDR and ML belong at the *group* level (many loci × subjects, Lessons 5–6), not inside one person's single peak; and no ML is used at all (principled counts + established statistics). For one person, one peak the honest move is simple POOLING: sum the read counts across all heterozygous SNPs in the peak into one combined IP ratio, compared against that person's own combined Input ratio. Pooling ~10-read SNPs (e.g., 15+20+25 = 60 reads) turns a coin flip into a readable ratio. This is exactly M6Allele's per-peak aggregation.

### Q2 (Lesson 4)
One person, one peak — pooled IP fraction 0.85 toward allele A, pooled Input fraction 0.52 toward A. In plain words, what does this tell us about that person's m6A at this peak — and is this already "ASM significant," or does it just *point* at something still needing a test?

> Answer (student): "I would guess it's ASM significant?"

Teacher's next move: Confirmed it strongly *points to* preferential methylation (0.85 vs 0.52 = one allele clearly stamped more). Corrected the leap to "significant": the gap is the *estimate*, but significance is decided by a test that weighs the gap against how many reads back it (0.85-vs-0.52 with ~10 reads is noise; with 10,000 reads it's overwhelming). Fisher's exact / beta-binomial turn estimate + depth into a p-value. Also distinguished: even a "suspicious" single-person call is only evidence of ASM *in that person*; the project's MDD-vs-control question needs aggregating people — Lesson 5. **Lesson 4 complete.**

---

## Lesson 5 — The new group-comparison plan

### Q1 (Lesson 5) — re-posed in concrete table form
[Table: five rows, one peak in GABRB1; columns label / IP fraA / Input fraA.]
(a) If we compared everyone's raw IP fraction directly, what is the problem — what does each person carry in that isn't about MDD? (b) To make MDD rows comparable to Ctrl rows, what one number per person do you reduce each row to? (Fill in: "How far is this person's ______ from their own ______.")

> Answer (student): "what caught my eye is male mdd has more skewed IP to Input ratio so the contrast between male cntl vs male mdd is better addressed I would say"

Teacher's next move: Affirmed the student spotted the actual signal — the male MDD rows show IP−Input gaps near 0.30 (0.84→0.54, 0.81→0.53) while controls sit near 0.01 (0.50→0.49, 0.52/0.53→0.5x); that skew is exactly what we hunt for. Also affirmed the sex-stratified instinct (male MDD vs male Ctrl), which foreshadows sex control. Completed the fill-in: the one number per person is their **skew = IP fraction relative to their own Input fraction** (their Input is that person's ruler) — not the raw IP. Both "person-to-person reference baseline" and "sex" are artifacts the group test must absorb.

> Student feedback: "For some reason I found your questions difficult to visualize and answer, you might want to adjust that." → Skill `socratic-teaching` Tailoring Log (2026-08-03): pose questions with an immediately imageable frame (named gene, tiny filled-in table with real numbers, fill-in-the-blank). Lesson 5 Q1 re-posed in a concrete table form.

### Q2 (Lesson 5)
Two things can slip in and fake an MDD-vs-Control difference. (1) The measurement artifact every brain drags in (Lesson 3's ruler). (2) The thing about WHO sits in each column (note: 3 of 5 ctrl rows are female, 0 of 2 mdd rows) — if that thing changes skew and is unbalanced across columns, it masquerades as MDD.

> Answer (student): "1. RMI? 2. imbalance in sample types?"

Teacher's next move: Affirmed both. (1) = reference-mapping bias (accepted "RMI" as reference-mapping; each person's Input AF is the baseline/ruler). (2) = sex imbalance (a participant-composition imbalance: MDD and Ctrl have different sex counts, and sex can shift skew — cohort is 10M/10F MDD vs 10M/9F Ctrl). Explained how the real plan absorbs each: reference-mapping bias → `offset(logit(Input AF))` so the null is the person's own fraction, not 0.5; sex → covariate (`condition + sex`) so the condition estimate is sex-adjusted, with within-sex contrasts kept exploratory.

### Q3 (Lesson 5)
Same peak, two MDD people, skew ~0.30 each: male_mdd_7 backed by 2,000 IP reads; male_mdd_3 backed by 15. (Table shown.) If we dropped all skews into a plain MDD-vs-Ctrl test ignoring the read pile, what are we throwing away? Why does fitting raw read counts (beta-binomial) beat comparing skews as if equally certain?

> Answer (student): "The mdd7 is more trustworthy"

Teacher's next move: Confirmed — trust/reliability lives in the read pile, so a 0.30 from 2,000 reads vastly outweighs the same 0.30 from 15. A framework that fits the raw counts (beta-binomial regression) weights each person by their evidence automatically; a plain test on skew fractions treats everyone's number as equally certain and wastes the reliability information. **Lesson 5 complete.**

---

## Lesson 6 — How to judge whether a plan is sound

### Q1 (Lesson 6)
The plan hard-codes cutoffs (MIN_IP_READS=10, SKEW_CUT=0.20, MIN_SUBJECTS=8). Pick one; give the two-sided cost — what strictness gains vs loses; what being too lax vs too strict invites. Walk it like a reviewer, not a defender.

> Answer (student): "For the MINIPREADS, if there's not enough reads in the IP we can't possibly sure whether the section has enough confidence, skew cut means the ratio is enough to shift the balance and such event should occur in ALL cohorts hence min subj = 8. I think min subj = 8 is a bit too strict. Just a feeling, and skew cut could lax out a little. But I think the paper has these addressed already"

Teacher's next move: Affirmed the reviewer instinct and the two-sided instinct. Corrected a conflation: SKEW_CUT (how far a single peak must lean off 50/50 — signal magnitude, per person) and MIN_SUBJECTS (how many people must pass at a peak — cross-subject sample size) are distinct; "should occur in all cohorts" is reproducibility/robustness, which belongs to min-subjects and to FDR, not to the skew cut. Accepted the "bit too strict on min_subj/skew" as a testable feeling, not a verdict. Key point: the plan does NOT bet on single cutoffs — `12_asm_compare.R` runs sweeps (MIN_SUBJECTS ∈ {6,8,12}, SKEW ∈ {0.1,0.2,0.3}), the hedge is that conclusions are checked across cutoffs.

### Q2 (Lesson 6)
The plan tries three tools in order: the careful one (glmmTMB), the middle one (lme4), then the bare one (glm). The bare `glm` is the only one that skips giving each person their own baseline — it treats all 39 people as interchangeable with one shared baseline. As reviewer: what's the real danger if a peak ends up tested with the bare tool? What should the report say so you're not silently reading a bare-tool result as if it were the careful one?

> Answer (student): "I don't get this question" → Skill Tailoring Log (2026-08-03): re-ask with the simplest concrete metaphor, one idea at a time, no jargon. Re-posed in plain form in chat.

> Answer (student, plain re-ask): "glmmTMB is more dynamic and it correctly assigns the correct input bases for each cohort? While glm regress everything into 1 single averaged out base?"

Teacher's next move: Confirmed correct in substance — glmmTMB/lme4 give each SUBJECT its own baseline (personal ruler); the bare glm averages everyone onto ONE shared baseline, so genuine per-person ruler differences get misread as signal. Refined the terms slightly: per-*subject* (not per-cohort), and the thing given is a baseline/ruler (their own Input), not "input bases." The output's `engine` column tells the reader exactly which tool ran, so a `glm` row is read with caution. **Lesson 6 complete.**

### Course wrap-up
Full course done: lessons 1–6 flowed from biology → variant calling → allele counts → per-person ASM → group analysis → judging the plan. The student can now read the pipeline scripts and methodology and defend/scrutinize each design choice (Input-only discovery, Input-as-null offset, sex covariate, count-weighted beta-binomial, BH-FDR, sensitivity sweeps, engine reporting).
