# Course Outline — Variant Calling: Purpose & Plan

*Course goal: understand what the `03_variant_calling/` step does and why, and be able to judge whether the ASM group-comparison plan is appropriate.*

**Lessons in order (we follow these; no skipping, no straying):**

1. **The biology setup** — DNA → RNA → protein; what m6A is; MeRIP-seq (IP vs Input); SNPs and heterozygous positions; why allele-level analysis is powerful; the two novel contributions.
2. **What variant calling does** — reads → alignment → GATK steps → VCF; why RNA-seq variant calling is special.
3. **From variants to allele counts** — ASEReadCounter: counting reads per allele in IP and Input fractions.
4. **Detecting ASM in one person** — M6Allele: is one allele more methylated than the other within a peak?
5. **The new group-comparison plan** — `11_asm_build_matrix.py` + `12_asm_compare.R`: the statistics and every judgment call.
6. **How to judge whether a plan is sound** — a reviewer's checklist.

**Method:** teacher asks one question at a time; student answers from intuition; every Q&A is recorded in `variant_calling_raw_review.md`; after each lesson a summary is appended to `variant_calling_summaries.md`.
