# Conda Environment (`mdd`)

`environment.yml` is the package list only: `name`, `channels`, three grouping
headers, and the dependencies. The rationale for the list lives here.

## How to build it

```bash
micromamba env create -f environment.yml        # creates env "mdd"
micromamba activate mdd
```

To add packages to an environment that already exists, rather than rebuilding it:

```bash
micromamba install -n mdd -c conda-forge bioconductor-biostrings r-data.table
```

`environment.yml` pins `python=3.9` and `r-base=4.3`, and uses `conda-forge`

The env built on this machine is not at those pins: `python` 3.9.23 (matches), `r-base` **4.4.1** (the pin asks 4.3), `r-data.table` 1.18.6.1, `bioconductor-biostrings` 2.74.0, `bioconductor-biomart` 2.62.0. A clean rebuild from `environment.yml` therefore does not reproduce the environment the scripts were developed against.
then `bioconda`. `defaults` is deliberately **not** enabled: bioconda packages
resolve against `conda-forge`, and mixing `defaults` in breaks the pins.

The JVM the GATK and Qualimap calls run on comes in transitively with `gatk4`
and `qualimap`; `openjdk` is not declared separately, which is why README can say the
environment needs a Java SDK without listing a Java package.

## Upstream Core Tools

| Package | Where it is used |
| --- | --- |
| `star` | `01_core_alignment/03_rrna_filter.sh`, `01_core_alignment/05_star_align.sh` |
| `fastp` | `01_core_alignment/02_trimming.sh` |
| `fastqc` | `01_core_alignment/01_raw_fastqc.sh`, `01_core_alignment/04_trimmed_fastqc.sh` |
| `qualimap` | `01_core_alignment/06_qualimap_bamqc.sh`, `01_core_alignment/07_qualimap_rnaseq.sh` |
| `subread` | `featureCounts` in `01_core_alignment/08_quantification.sh`, `01_core_alignment/09_counts_qc.sh` |
| `multiqc` | `01_core_alignment/10_multiqc.sh`, `multiqc_config.yaml` |
| `picard` | `multiqc_config.yaml` (BAM header fields MultiQC reports) |
| `samtools` | `03_variant_calling/01`–`04`, `utils/slice_igv_bams.sh` |
| `gatk4` | `03_variant_calling/01`–`07` |
| `bcftools` | `03_variant_calling/07_variant_filtration.sh` |
| `aspera-cli` | `_hpc/download.sbatch` — see `docs/aspera_setup.md` |

Declared for armoured capability, no consumer in the delivered scripts:
`bedtools` (interval arithmetic), `trimmomatic` (adapter trimming; `fastp` is
what 01 actually runs), `meme` + `maven` (the motif-enrichment arm archived in
`_legacy/utils/ame_motif_enrichment.sh`).

## Python Stack

| Package | Where it is used |
| --- | --- |
| `pandas` | `03_variant_calling/09`, `10`, `11`, `14`, `utils/qc02_cell_type_deconvolution.py` |
| `numpy` | `03_variant_calling/09`, `10`, `14`, the four `utils/plot_gatk*.py` / `plot_alt*.py`, `utils/qc02_cell_type_deconvolution.py` |
| `scipy` | `utils/qc02_cell_type_deconvolution.py`, and the backend `statsmodels` solves against |
| `statsmodels` | `03_variant_calling/14_gene_burden_test.py` — negative-binomial GLM with a `log(total_variants)` offset, plus the Benjamini-Hochberg FDR |
| `matplotlib` | the `utils/plot_*.py` figures |

`numpy` and `scipy` used to be pulled in transitively by `pandas`. They are
declared explicitly because `14` depends on them directly, and the offset is
the mechanism that cancels the measured 1.254x MDD/Control coverage imbalance
(`docs/polygenic_burden_methodology.md` section 6) — a transitive pin is not a
contract.

Declared, no consumer in the delivered scripts: `tqdm`, `requests`, `pyfaidx`,
`pygenometracks`, `igv-reports`, `lxml`. Their consumers (reference FASTA
builder, IGV session slicer, curation CSV export, Ensembl REST fetcher) were
removed from `03_variant_calling/` in the Input-only cleanup. `pip` is the
anchor for the `pip:` sub-list some builds append.

## R Stack

The delivered scripts load these libraries:

| Library | Package | Where |
| --- | --- | --- |
| `tidyverse` | `r-tidyverse` | `02_dge_analysis/09_heatmap_pathway.R` |
| `dplyr` / `readr` / `stringr` / `yaml` | `r-dplyr` `r-readr` `r-stringr` `r-yaml` | `02_dge_analysis/01`, `06`, `07`, `08` |
| `DESeq2` | `r-deseq2` | `02_dge_analysis/02`, `03`, `04`, `05`, `09` |
| `ashr` | `bioconductor-ashr` | `02_dge_analysis/04_deg_analysis.R:166,168` — the shrink estimator (`type = "ashr"`) |
| `ggplot2` | `r-ggplot2` | `02_dge_analysis/03`, `05`, `08`, `03_variant_calling/15`, `17`, `20`, `utils/plot_snpeff_functional_consequences.R` |
| `ggrepel` | `r-ggrepel` | `02_dge_analysis/05_dge_plotting.R` |
| `ggtext` | `r-ggtext` | `02_dge_analysis/08_enrichment_plotting.R` |
| `scales` | `r-scales` | `utils/plot_snpeff_functional_consequences.R:128` |
| `plotly` + `htmlwidgets` | `r-plotly` `r-htmlwidgets` | `02_dge_analysis/03_pca_analysis.R` |
| `ComplexHeatmap` + `circlize` + `RColorBrewer` | `r-complexheatmap` `r-circlize` `r-rcolorbrewer` | `02_dge_analysis/09_heatmap_pathway.R:15,16` |
| `patchwork` | `r-patchwork` | `03_variant_calling/17_plot_burden_control_vs_mdd.R:26` |
| `biomaRt` | `bioconductor-biomart` | `02_dge_analysis/06_annotate_results.R` |
| `clusterProfiler` | `r-clusterprofiler` | `02_dge_analysis/07_enrichment_analysis.R` |
| `org.Hs.eg.db` | `r-org.hs.eg.db` | `02_dge_analysis/07`, `09` |
| `msigdbr` | `r-msigdbr` | `02_dge_analysis/07`, `03_variant_calling/13_export_pathway_sets.R` |
| `data.table` | `r-data.table` | `03_variant_calling/13`–`22` |

Declared, not loaded by any delivered script: `r-limma` (the linear-model arm
of a previous DGE design), `r-biocmanager` (package-install helper),
`r-jsonlite` (transitive), `bioconductor-go.db` (the Gene Ontology annotation
package `clusterProfiler` resolves GO term names from; declared explicitly so
`02_dge_analysis/07_enrichment_analysis.R` runs offline).

### Packages removed from the environment

`bioconductor-tress`, `bioconductor-rrho` (the DGE plots, removed from 02),
`bioconductor-apeglm` (the previous DE shrink estimator, replaced by `ashr`).

### Packages added to the environment

`r-patchwork`, `r-circlize`, `r-rcolorbrewer`, `r-scales` — each is loaded by
a delivered script and was previously installed only as a transitive dependency.
A direct `library()` call deserves a pin, so they are declared. `bioconductor-biomart` is restored here: it is declared at HEAD, installed in env `mdd` (2.62.0), and loaded by `02_dge_analysis/06_annotate_results.R`, so a package list that omits it cannot build an environment the script can run in.

## m6A-SNP resources

The m6A-SNP cross-validation arms (m6AVar, RMVar, m6A-TSHub) have **no conda
package**: they are web resources, not installable libraries. This environment
cannot reach them, so cross-validation is done by hand through the URL columns
`03_variant_calling/10_prioritize_mdd_candidate_snps.py` attaches to every
candidate row (OMIM, DECIPHER, PubMed, gnomAD gene and variant, GENCC,
EBI GWAS Catalog).
