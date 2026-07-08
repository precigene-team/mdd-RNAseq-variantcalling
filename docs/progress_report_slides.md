---
marp: true
theme: uncover
class: invert
paginate: true
style: |
  @import url('https://fonts.googleapis.com/css2?family=Inter:wght@300;400;600;700&display=swap');

  :root {
    --bg: #0f111a;
    --surface: #1a1d2e;
    --border: #2a2d3e;
    --accent: #60a5fa;
    --text: #e2e8f0;
    --text-dim: #94a3b8;
  }

  section {
    background: var(--bg);
    color: var(--text);
    font-family: 'Inter', system-ui, sans-serif;
    font-size: 22px;
    padding: 50px 60px;
    line-height: 1.6;
    justify-content: flex-start;
  }

  h1 {
    color: var(--accent);
    font-size: 1.5em;
    font-weight: 700;
    border-left: 4px solid var(--accent);
    padding-left: 16px;
    margin: 0 0 20px 0;
  }
  h2 {
    color: var(--accent);
    font-size: 0.95em;
    font-weight: 600;
    margin: 0 0 10px 0;
  }
  h3 {
    color: var(--text-dim);
    font-weight: 400;
    font-size: 1.05em;
    margin: 0 0 14px 0;
  }

  ul { margin: 4px 0 10px 0; padding-left: 20px; }
  li { margin-bottom: 5px; font-size: 0.88em; line-height: 1.45; }
  strong { color: #f1f5f9; }
  code { background: var(--surface); border: 1px solid var(--border); padding: 1px 5px; border-radius: 3px; font-size: 0.82em; color: var(--accent); }
  em { color: var(--text-dim); font-style: normal; }
  a { color: var(--accent); text-decoration: none; }

  table { font-size: 0.78em; border-collapse: collapse; width: 100%; margin-top: 6px; }
  th { background: var(--surface); color: var(--accent); border: 1px solid var(--border); padding: 6px 8px; font-weight: 600; text-align: left; }
  td { border: 1px solid var(--border); padding: 6px 8px; }

  img { border-radius: 6px; border: 1px solid var(--border); object-fit: contain; }

  .half { display: grid; grid-template-columns: 1fr 1fr; gap: 30px; align-items: start; }
  .text-fig { display: grid; grid-template-columns: 1fr 1fr; gap: 24px; align-items: center; }
  .text-fig img { max-height: 420px; width: 100%; }
  .fig-only { display: flex; justify-content: center; align-items: center; gap: 20px; }
  .fig-only img { max-height: 450px; }
  .fig-pair { display: grid; grid-template-columns: 1fr 1fr; gap: 16px; align-items: center; }
  .fig-pair img { max-height: 400px; width: 100%; }

  .badge { display: inline-block; background: var(--surface); border: 1px solid var(--border); padding: 2px 6px; border-radius: 4px; margin: 2px; font-size: 0.68em; color: var(--accent); font-family: monospace; }

  /* ── Title slide ── */
  section.title {
    justify-content: center;
    text-align: center;
    padding: 60px 80px;
  }
  section.title h1 {
    border-left: none;
    padding-left: 0;
    font-size: 1.9em;
    line-height: 1.3;
    margin-bottom: 16px;
  }
  section.title h3 { font-size: 1.1em; color: var(--text-dim); margin-bottom: 24px; }
  section.title p { font-size: 0.9em; color: var(--text-dim); margin: 3px 0; }

  /* ── Filler / question slide ── */
  section.filler {
    justify-content: center;
    text-align: center;
    padding: 80px 100px;
    background: linear-gradient(135deg, #0f111a 0%, #1a1d2e 100%);
  }
  section.filler h1 {
    border-left: none;
    padding-left: 0;
    font-size: 1.55em;
    line-height: 1.5;
    font-weight: 600;
  }
---

<!-- _class: title -->

# Prefrontal Cortex Transcriptomics & m⁶A Epitranscriptomics in MDD

### Pipeline Progress Report

**Scope:** Paired MeRIP-seq (Input vs. IP)
**Model:** Sex-stratified baseline & interaction contrasts
**Key Targets:** DGE · DMG · MACS2 Validation · GATK Variant Calling

---

# Biological & Clinical Background

<div class="half">
<div>

## m⁶A Epitranscriptomics
* Most prevalent internal mRNA modification — regulates splicing, stability, translation
* **Writers:** `METTL3/14`, `WTAP`
* **Erasers:** `FTO`, `ALKBH5`
* **Readers:** `YTHDF1-3`, `YTHDC1-2`
* Implicated in synaptic plasticity & stress response

</div>
<div>

## Sex-Specific Dysregulation
* Females: 2× MDD prevalence vs. males
* PFC transcriptomic & methylation networks are highly sex-dependent
* Males vs. females show opposing molecular alterations
* Justifies **sex-stratified** analysis design

</div>
</div>

---

# Cohort Stratification & Library Prep

<div class="half">
<div>

## MeRIP-Seq Design
* Paired-end from postmortem PFC tissue (**39 subjects**)
* **Input fraction:** steady-state total transcriptome
* **IP fraction:** anti-m⁶A antibody capture
* **Strandedness:** reverse-stranded (antisense R1)

</div>
<div>

## Subject Gating

| Group | Initial | Final | Libraries |
| :--- | :---: | :---: | :---: |
| Male MDD | 10 | 10 | 20 |
| Male Ctrl | 10 | 9 | 18 |
| Female MDD | 10 | 10 | 20 |
| Female Ctrl | 9 | 9 | 18 |
| **Total** | **39** | **38** | **76** |

`male_ctrl_1` removed (11.6M < 12M gate)

</div>
</div>

---

# Upstream: Alignment & MultiQC Diagnostics

<div class="half">
<div>

## STAR Alignment
* Initial mapping rates low (20–30%) for Input fraction
* Genome index rebuilt with primary assembly (`sjdbOverhang 149`)
* Rates improved slightly but remain suboptimal
* **Cause:** incomplete rRNA depletion during Input library prep — rRNA fragments multi-map to hundreds of repetitive rDNA loci

</div>
<div>

## MultiQC Summary Statistics

| Metric | IP (m⁶A) | Input |
| :--- | :---: | :---: |
| **Unique Align %** | 93.25% | 52.60% |
| **Multi-map %** | 2.49% | 29.03% |
| **Duplication %** | 64.52% | 57.55% |
| **Assign Rate %** | 68.24% | 8.30% |

*IP self-depletes rRNA via antibody specificity*

</div>
</div>

---

# Upstream: featureCounts Assignment Refactoring

<div class="half">
<div>

## Problem
* Reported featureCounts assignment rate was **<50%** for Input samples
* Denominator inflated by multi-mapping rDNA alignments that featureCounts discards

## Resolution
* Refactored QC to compute assignment rate relative to **uniquely mapped reads only**
* **True assignment rate: 50.5%**

</div>
<div>

## Worst Samples (Library Size)

| Sample | Exonic Counts |
| :--- | :---: |
| `female_mdd_2_input` | 5.8M ⚠️ |
| `female_ctrl_6_input` | 9.9M |
| `male_ctrl_1_input` | 11.6M |
| `male_mdd_2_input` | 12.1M |
| `female_mdd_8_input` | 12.2M |

**Gate:** 12M exonic counts threshold applied

</div>
</div>

---
<!-- _class: filler -->

# How do we ensure exonic counts represent valid human genes rather than rRNA contamination noise?

---


# QC-01 · Postmortem Confounders (PMI, RIN, Age, pH)

<div class="half">
<div>

## The Risk
* Postmortem tissues degrade dynamically — longer PMI or lower RIN in MDD cohort would call **degradation artifacts as DEGs**
* Individual-level covariates were **not publicly released** (GSE275676) — donor privacy protection
* Without per-subject RIN/PMI values, covariate adjustment in the model is not possible

</div>
<div>

## Resolution & Outcome
* Authors explicitly matched cases and controls for **PMI, tissue pH, and RIN** during cohort recruitment (Mitsuhashi et al. 2024, *iScience*)
* **12M exonic counts gate** calibrated by matching published cohort sizes — correctly excludes `male_ctrl_1` (11.6M), aligning with paper's reported 9 Male Controls
* **Conclusion:** Unadjusted baseline models are robust — group-level matching compensates for the absence of individual covariates

</div>
</div>

---

# QC-02 · Cell-Type Heterogeneity Confounding

<div class="half">
<div>

## The Risk
* Bulk RNA-seq cannot distinguish: *is a gene downregulated*, or *are fewer of those cells present*?
* Neuronal atrophy and glial loss in MDD would make cell-type-specific transcripts appear differentially expressed

## Method
* Marker Gene Profile (MGP) deconvolution across 36 gated subjects
* Consensus human brain markers for: Neurons, Astrocytes, Oligodendrocytes, Microglia, Endothelial cells
* MDD vs. Control t-tests per cell-type

</div>
<div>

## Outcome

| Cell Type | MDD vs. Ctrl | *p*-value |
| :--- | :--- | :---: |
| **Astrocytes** | ↓ Significant | **0.0389** |
| Neurons | ↑ Trend | 0.0956 |
| Oligodendrocytes | No change | > 0.25 |
| Microglia | No change | > 0.25 |
| Endothelial | No change | > 0.25 |

**Astrocyte marker loss** is a well-established neuropathological hallmark of MDD in human PFC — documented as a biological driver, not a technical artifact

</div>
</div>

---

# QC-03 · Gene-Level vs. Peak-Level Resolution

<div class="half">
<div>

## The Risk
* Gene-body read summation can **cancel out site-specific shifts** — a gain at one m⁶A site may be masked by a loss at another within the same gene body
* Pure gene-level DGE may miss isoform- or site-specific methylation changes

## Method
* MACS2 narrowPeak calling on all 36 gated subjects (HPC)
* Exonic overlap audit via `02_audit_gene_vs_peak.py` — maps peak coordinates to exonic intervals across the cohort

</div>
<div>

## Outcome
* **Exonic mapping rate:** uniform **~61–63%** across all 36 subjects → confirms technical pipeline consistency
* **Female DMG concordance:** **81%** — 20 of 25 interaction DMGs are physically supported by narrowPeaks → high biological confidence
* **Male DMG concordance:** **33%** — 1.3 of 4 DMGs supported → male cohort shows higher individual peak heterogeneity
* **Conclusion:** Gene-level model is well-validated by independent peak evidence, particularly for the female interaction results

</div>
</div>

---

# QC-04 · Antibody Cross-Reactivity & Spike-in Absence

<div class="half">
<div>

## The Risk
* Anti-m⁶A antibodies are semi-quantitative and can capture **non-m⁶A** modified sites (e.g. m⁶Am, m¹A)
* Without synthetic RNA spike-in controls sequenced in the same run, it is impossible to calculate **absolute methylation stoichiometry** (fraction of methylated copies per transcript)

## Audit Finding
* Paper text confirmed: cDNA libraries were generated **without synthetic spike-in controls** (Mitsuhashi et al. 2024)
* Spike-in sequences were absent from the raw FASTQ files

</div>
<div>

## Resolution & Scope of Results
* All interaction DGE results express **relative abundance shifts** — the ratio of IP/Input enrichment in MDD relative to Control
$$\Delta m^6A = \frac{(\text{IP}/\text{Input})_{\text{MDD}}}{(\text{IP}/\text{Input})_{\text{Ctrl}}}$$
* This is the **standard methodology** for clinical human postmortem MeRIP-seq cohorts
* Absolute stoichiometry is beyond scope — results represent **differential relative methylation**, which is interpretable and biologically meaningful

</div>
</div>

---
<!-- _class: filler -->

# Differential Gene Expression Analysis

---

# DGE Step 1 · Gene Mapping & Pre-processing

<div class="half">
<div>

## `01_gene_mapping.R`
* Connects to BioMart → maps Ensembl IDs to HGNC symbols and chromosomal coordinates
* Filters out duplicate gene entries
* Output: `01_annotated_counts.csv`, `01_gene_mapper.csv`

</div>
<div>

## `02_rna_processing.R`
* Low-count filter: count ≥ 10 in ≥ 50% of samples
* **60,600 raw features → ~18,000 high-confidence transcripts**
* DESeq2 VST normalization
* Outputs: `02_RIPseq_filtered_counts.csv`, `02_RNAseq_filtered_counts.csv`

</div>
</div>

---
<!-- _class: filler -->

# Does global variance reflect biological factors like sex or disease, or is it dominated by technical batch effects?

---

# DGE Step 2 · PCA: RNAseq Female

<div class="text-fig">
<div>

## Female RNAseq (Input) PCA
* MDD & Control samples **do not** separate distinctly
* High inter-subject heterogeneity typical of postmortem brain tissue
* Subtle disease signal requires statistical contrasts to extract
* *Method:* VST-normalized PCA via `03_pca_analysis.R`

</div>
<div>

![](../results/02_dge_analysis/figures/03_RNAseq_PCA_2D_female_only.png)

</div>
</div>

---

# DGE Step 2 · PCA: RNAseq Male

<div class="text-fig">
<div>

## Male RNAseq (Input) PCA
* Mirrors female cohort — interleaved MDD/Control
* Confirms within-subject paired models are required
* *Method:* VST-normalized PCA via `03_pca_analysis.R`

</div>
<div>

![](../results/02_dge_analysis/figures/03_RNAseq_PCA_2D_male_only.png)

</div>
</div>

---

# DGE Step 2 · PCA: RIPseq Female

<div class="text-fig">
<div>

## Female RIPseq (IP) PCA
* Similar pattern to RNAseq — no clear MDD vs. Control clustering
* The IP fraction captures m⁶A-enriched transcripts, showing the methylation landscape also varies by individual, not condition alone
* *Method:* VST-normalized PCA via `03_pca_analysis.R`

</div>
<div>

![](../results/02_dge_analysis/figures/03_RIPseq_PCA_2D_female_only.png)

</div>
</div>

---

# DGE Step 2 · PCA: RIPseq Male

<div class="text-fig">
<div>

## Male RIPseq (IP) PCA
* Comparable to female IP — no distinct condition clusters
* Reinforces the need for paired interaction models to detect relative methylation shifts
* *Method:* VST-normalized PCA via `03_pca_analysis.R`

</div>
<div>

![](../results/02_dge_analysis/figures/03_RIPseq_PCA_2D_male_only.png)

</div>
</div>

---

# DGE Step 2 · PCA: 3D Interactive Views

* **3D PCA visualizations** generated as interactive HTML widgets for both RNAseq and RIPseq fractions:
  * `03_RNAseq_PCA_3D.html` — full cohort, uncorrected
  * `03_RNAseq_PCA_3D_sex_corrected.html` — sex-effect removed
  * `03_RIPseq_PCA_3D.html` — full cohort, uncorrected
  * `03_RIPseq_PCA_3D_sex_corrected.html` — sex-effect removed
* **Key observation:** PCA reveals a **distinct split between male and female** samples on PC1/PC2, but MDD vs. Control is **not clearly visualized** — confirming the sex-stratified model design
* 3D views available in `results/02_dge_analysis/figures/`

---
<!-- _class: filler -->

# Are steady-state mRNA levels altered in PFC of MDD subjects in a sex-specific manner?

---

# DGE Step 3 · DESeq2 Contrast Design

<div class="half">
<div>

## Baseline Expression Contrast
$$\text{MDD}_{\text{Input}} \;\text{vs}\; \text{Ctrl}_{\text{Input}}$$
* Tests steady-state transcription changes
* Run separately for **Female** and **Male** cohorts
* Also: **Sex-Interaction** contrast (opposing effects across sexes)

## Thresholds
* |log₂FC| ≥ 0.2, *p* < 0.05
* Relaxed from 0.5 — retains subtle PFC brain signatures

</div>
<div>

## IP Fraction Contrast
$$\text{MDD}_{\text{IP}} \;\text{vs}\; \text{Ctrl}_{\text{IP}}$$
* Tests m⁶A-enriched transcript abundance
* Run separately for **Female** and **Male** cohorts

## Interaction Contrast (Methylation)
$$\Delta m^6A = \frac{(\text{IP}/\text{Input})_{\text{MDD}}}{(\text{IP}/\text{Input})_{\text{Ctrl}}}$$
* Isolates **relative methylation** changes from transcription changes

</div>
</div>

---

# DGE Step 3A · Female Baseline DEGs — Volcano

<div class="text-fig">
<div>

## Female Input: 38 DEGs
* ↑ `CARTPT` (1.13), `PLVAP` (0.74), `IRAK3` (0.40), `NECTIN3` (0.28)
* ↓ `DTHD1` (−2.05), `DRC1` (−0.40), `ARSJ` (−0.36)
* *Source:* `06_deg_female_input_sig.csv`

</div>
<div>

![](../results/02_dge_analysis/figures/05_volcano_female_input.png)

</div>
</div>

---

# DGE Step 3A · Female Baseline DEGs — MA Plot

<div class="text-fig">
<div>

## Female Input: MA Distribution
* MA plot shows log₂FC vs. mean normalized expression
* Most genes cluster near zero — 38 significant genes span a broad expression range
* Confirms no systematic bias in fold-change estimation

</div>
<div>

![](../results/02_dge_analysis/figures/05_ma_female_input.png)

</div>
</div>

---

# DGE Step 3A · Male Baseline DEGs — Volcano

<div class="text-fig">
<div>

## Male Input: 66 DEGs
* ↑ `AP1M2` (1.48), `CACNA1S` (0.68), `CELSR1` (0.69)
* ↓ `DGKK` (−2.64), `CD14` (−1.23), `CASP1` (−0.98), `BCL2` (−0.25)
* *Source:* `06_deg_male_input_sig.csv`

</div>
<div>

![](../results/02_dge_analysis/figures/05_volcano_male_input.png)

</div>
</div>

---

# DGE Step 3A · Male Baseline DEGs — MA Plot

<div class="text-fig">
<div>

## Male Input: MA Distribution
* MA plot confirms fold-change symmetry across expression levels
* More DEGs than female (66 vs. 38) — male cohort shows broader immune/apoptotic signature

</div>
<div>

![](../results/02_dge_analysis/figures/05_ma_male_input.png)

</div>
</div>

---

# DGE Step 3A · IP Fraction DEGs — Female

<div class="fig-pair">

![](../results/02_dge_analysis/figures/05_volcano_female_ip.png)

![](../results/02_dge_analysis/figures/05_ma_female_ip.png)

</div>

*Source:* `06_deg_female_ip_sig.csv`

---

# DGE Step 3A · IP Fraction DEGs — Male

<div class="fig-pair">

![](../results/02_dge_analysis/figures/05_volcano_male_ip.png)

![](../results/02_dge_analysis/figures/05_ma_male_ip.png)

</div>

*Source:* `06_deg_male_ip_sig.csv`

---
<!-- _class: filler -->

# Do baseline transcriptomic changes converge on shared pathway signatures?

---

# DGE Step 3B · Female Baseline Pathways (Gene-Level → Pathway-Level)

<div class="text-fig">
<div>

## Enriched GO (ORA)
* **Cilium movement** — *p* = 2.01 × 10⁻¹⁵, 29 genes (e.g. `CFAP73`, `DRC1`, `DNAI2`)
* **Axoneme assembly** — *p* = 1.80 × 10⁻¹³
* **Microtubule bundle formation** — *p* = 2.44 × 10⁻¹¹

</div>
<div>

![](../results/02_dge_analysis/figures/08_deg_female_input_gsea_go_nes.png)

</div>
</div>

---

# DGE Step 3B · Female Baseline — Hallmark GSEA & Heatmap

<div class="fig-pair">

![](../results/02_dge_analysis/figures/08_deg_female_input_gsea_hallmark_nes.png)

![](../results/02_dge_analysis/figures/09_heatmap_pathway_deg_female_input.png)

</div>

---

# DGE Step 3B · Male Baseline Pathways (Gene-Level → Pathway-Level)

<div class="text-fig">
<div>

## Enriched GO (ORA)
* **Extrinsic apoptotic signaling** — *p* = 8.49 × 10⁻¹⁰, 37 genes (e.g. `BCL2`, `CASP1`, `TNF`)
* **IL-6 production** — *p* = 7.53 × 10⁻⁹
* **Type II interferon response** — *p* = 1.68 × 10⁻⁸

</div>
<div>

![](../results/02_dge_analysis/figures/08_deg_male_input_gsea_go_nes.png)

</div>
</div>

---

# DGE Step 3B · Male Baseline — Hallmark GSEA & Heatmap

<div class="fig-pair">

![](../results/02_dge_analysis/figures/08_deg_male_input_gsea_hallmark_nes.png)

![](../results/02_dge_analysis/figures/09_heatmap_pathway_deg_male_input.png)

</div>

---

# DGE Step 3B · IP Fraction Pathways — Female

<div class="fig-pair">

![](../results/02_dge_analysis/figures/08_deg_female_ip_gsea_go_nes.png)

![](../results/02_dge_analysis/figures/09_heatmap_pathway_deg_female_ip.png)

</div>

*Hallmark:* `08_deg_female_ip_gsea_hallmark_nes.png`

---

# DGE Step 3B · IP Fraction Pathways — Male

<div class="fig-pair">

![](../results/02_dge_analysis/figures/08_deg_male_ip_gsea_go_nes.png)

![](../results/02_dge_analysis/figures/09_heatmap_pathway_deg_male_ip.png)

</div>

*Hallmark:* `08_deg_male_ip_gsea_hallmark_nes.png`

---

# DGE Step 3B · Sex-Interaction Baseline Pathways

<div class="text-fig">
<div>

## Sex-Interaction Input
* Genes with **opposing** fold-change directions between male and female
* *Source:* `06_deg_sex_interaction_input_sig.csv`

</div>
<div>

![](../results/02_dge_analysis/figures/08_deg_sex_interaction_input_gsea_go_nes.png)

</div>
</div>

---

# DGE Step 3B · Sex-Interaction — Hallmark & Heatmap

<div class="fig-pair">

![](../results/02_dge_analysis/figures/08_deg_sex_interaction_input_gsea_hallmark_nes.png)

![](../results/02_dge_analysis/figures/09_heatmap_pathway_deg_sex_interaction_input.png)

</div>

---
<!-- _class: filler -->

# How do we isolate post-transcriptional m⁶A shifts from changes in total transcription?

---

# DGE Step 4A · Female Interaction DMGs — Volcano

<div class="text-fig">
<div>

## Interaction Contrast
$$\Delta m^6A = \frac{(\text{IP}/\text{Input})_{\text{MDD}}}{(\text{IP}/\text{Input})_{\text{Ctrl}}}$$

**25 Female DMGs:**

<span class="badge">C16orf86</span><span class="badge">C19orf53</span><span class="badge">CCBE1</span><span class="badge">CCDC169</span><span class="badge">CDC123</span><span class="badge">CREB3L4</span><span class="badge">FXYD5</span><span class="badge">GGN</span><span class="badge">GNE</span><span class="badge">H1-4</span><span class="badge">KCNK10</span><span class="badge">KPNA2</span><span class="badge">LAMTOR5</span><span class="badge">MAN1A2</span><span class="badge">NECTIN2</span><span class="badge">OR2L5</span><span class="badge">PINK1</span><span class="badge">RAB9B</span><span class="badge">RPL36AL</span><span class="badge">RPL8</span><span class="badge">RPS13</span><span class="badge">RTCB</span><span class="badge">TBX2</span><span class="badge">TFPT</span><span class="badge">USE1</span>

</div>
<div>

![](../results/02_dge_analysis/figures/05_volcano_female_meth.png)

</div>
</div>

---

# DGE Step 4A · Female Interaction DMGs — MA Plot

<div class="text-fig">
<div>

## Female Interaction: MA Distribution
* Shows the spread of interaction log₂FC across mean expression
* 25 DMGs validated at |log₂FC| ≥ 0.2, *p* < 0.05

</div>
<div>

![](../results/02_dge_analysis/figures/05_ma_female_meth.png)

</div>
</div>

---

# DGE Step 4A · Male Interaction DMGs — Volcano

<div class="text-fig">
<div>

## Interaction Contrast
$$\Delta m^6A = \frac{(\text{IP}/\text{Input})_{\text{MDD}}}{(\text{IP}/\text{Input})_{\text{Ctrl}}}$$

**4 Male DMGs:**
<span class="badge">AK6</span><span class="badge">C19orf81</span><span class="badge">C1QL4</span><span class="badge">GRIN3B</span>

**2 Sex-Interaction DMGs** *(opposing shifts):*
<span class="badge">HMGCL</span><span class="badge">TLL2</span>

</div>
<div>

![](../results/02_dge_analysis/figures/05_volcano_male_meth.png)

</div>
</div>

---

# DGE Step 4A · Male Interaction DMGs — MA Plot

<div class="text-fig">
<div>

## Male Interaction: MA Distribution
* Only 4 DMGs — reflects lower signal-to-noise in the male methylation cohort
* `GRIN3B` is the dominant hit (synaptic NMDA receptor subunit)

</div>
<div>

![](../results/02_dge_analysis/figures/05_ma_male_meth.png)

</div>
</div>

---
<!-- _class: filler -->

# Do differentially methylated transcripts represent coordinated biological processes?

---

# DGE Step 4B · Female Interaction Pathways

<div class="text-fig">
<div>

## Female Hallmark GSEA
* `PROTEIN_SECRETION` — NES = −2.22, *p* = 2.07 × 10⁻⁷
* `MYC_TARGETS_V1` — NES = −1.86
* `MTORC1_SIGNALING` — NES = −1.84
* `OXIDATIVE_PHOSPHORYLATION` — NES = −1.77

</div>
<div>

![](../results/02_dge_analysis/figures/08_meth_female_interaction_gsea_go_nes.png)

</div>
</div>

---

# DGE Step 4B · Female Interaction — Hallmark & Heatmap

<div class="fig-pair">

![](../results/02_dge_analysis/figures/08_meth_female_interaction_gsea_hallmark_nes.png)

![](../results/02_dge_analysis/figures/09_heatmap_pathway_meth_female_interaction.png)

</div>

---

# DGE Step 4B · Male Interaction Pathways

<div class="text-fig">
<div>

## Male Synaptic Regulation
* Hyper-methylation of **`GRIN3B`** — glutamate NMDA receptor subunit 3B
* Suggests post-transcriptional gating of synaptic receptors in male PFC
* Limited pathway enrichment due to low DMG count (n = 4)

</div>
<div>

![](../results/02_dge_analysis/figures/08_meth_male_interaction_gsea_go_nes.png)

</div>
</div>

---

# DGE Step 4B · Male Interaction — Hallmark & Heatmap

<div class="fig-pair">

![](../results/02_dge_analysis/figures/08_meth_male_interaction_gsea_hallmark_nes.png)

![](../results/02_dge_analysis/figures/09_heatmap_pathway_meth_male_interaction.png)

</div>

---

# DGE Step 4B · Methylation Sex-Interaction Pathways

<div class="text-fig">
<div>

## Sex-Interaction Methylation
* Genes with **opposing** methylation shifts between male and female
* 2 DMGs: `HMGCL`, `TLL2`
* *Source:* `06_meth_sex_interaction_sig.csv`

</div>
<div>

![](../results/02_dge_analysis/figures/08_meth_sex_interaction_gsea_go_nes.png)

</div>
</div>

---
<!-- _class: filler -->

# Do gene-level methylation changes correspond to physical m⁶A peaks called by MACS2?

---

# DGE Step 5 · Peak-Calling Validation (QC-03)

<div class="half">
<div>

## Overlap Audit Strategy
* MACS2 narrowPeak calling on 36 gated subjects
* Exonic overlap via `02_audit_gene_vs_peak.py`
* **Uniform exonic mapping: ~61–63%** across all subjects

</div>
<div>

## Concordance Results
* **Female DMGs:** **81%** (20/25 supported by narrowPeaks)
* **Male DMGs:** **33%** (1.3/4 supported)
* Male results show higher subject-level peak heterogeneity
* Female interaction model is well-validated by independent peak evidence

</div>
</div>

---

# Future Work · Variant Calling Pipeline

<div class="half">
<div>

## GATK Pipeline Progress
- [x] Read Groups (`01_read_groups.sh`)
- [x] PCR duplicate flagging (`02_mark_duplicates.sh`)
- [x] Exon-exon split reads (`03_split_n_cigar.sh`)
- [ ] **HaplotypeCaller** — *in progress*, local CPU workers

</div>
<div>

## Next Steps
* Complete HC traversal across all chromosomes
* Apply GATK hard filtration → high-confidence SNPs
* Map SNPs to **DRACH** m⁶A consensus motif
* Identify allele-specific methylation (ASM) via IP/Input allele ratio deviation

</div>
</div>

---

# Future Work · Epitranscriptomic Elaboration

<div class="half">
<div>

## m⁶A-SNP Mapping
* Map high-confidence SNPs to DRACH consensus: [A/G/U][A/G]**A**C[A/C/U]
* Predict mutations that create or destroy m⁶A sites (motif disruption analysis)

</div>
<div>

## Allele-Specific Methylation (ASM)
* Compare allele ratios in paired IP vs. Input fractions
* Significant deviation → allele-specific methylation (parental allele bias)
* **Integration:** Link DMGs/m⁶A-SNPs to clinical postmortem covariates

</div>
</div>
