# RNA-seq Reanalysis of a SARS-CoV-2 Infection Dataset: Differential Expression and Functional Enrichment

This project demonstrates a reproducible RNA-seq analysis workflow using a publicly available SARS-CoV-2 infection dataset.

Differential expression
via **DESeq2**, functional enrichment via **GSEA** (GO Biological Process /
Cellular Component / Molecular Function, and KEGG pathways).


📊 **[RNA-seq Analysis Report](https://claude.ai/code/artifact/821e1224-4e65-40cc-b8d3-4f124b78354c)**

## Contents

- [Data](#data)
- [Computational Methods](#computational-methods)
- [Environment Setup](#environment-setup)
- [Running the Analysis Workflow](#running-the-analysis-workflow)
- [Directory Structure](#directory-structure)
- [Output Files](#output-files)
- [Notes](#notes)
- [References](#references)


## Data

**Source:** GEO accession [GSE164073](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE164073) <br>
**File:** `data/GSE164073_Eye_count_matrix.csv` — count matrix, 27,946 genes × 18 samples <br>
**Design:** 3 tissues (cornea, limbus, sclera) × 2 conditions (mock, SARS-CoV-2-infected) × 3 replicates <br>


## Computational Methods

### Principal Component Analysis



### Differential expression

Differential expression (DE) analysis was performed between SARS-CoV-2-infected and mock samples in each ocular tissue region (cornea, limbus, and sclera) using DESeq2 (R/Bioconductor package). 
Genes with zero counts across all samples were removed prior to analysis. DESeq2 estimates size factors using the median-of-ratios method for normalization and models read counts 
using a negative binomial distribution to account for biological variability and overdispersion.  Statistical significance was assessed using the Wald test, with p-values adjusted 
using the Benjamini–Hochberg false discovery rate (FDR) correction.. Genes with adjusted p-value (padj) < 0.05 and |log2FoldChange| > 1.5 were considered significantly differentially expressed.

**Normalization:** DESeq2 median-of-ratios (size factors) for the DE test
  itself; VST (variance-stabilizing transformation, `blind = TRUE`) is used
  separately, only for visualization (PCA, sample-distance heatmap) — never
  for the statistical test.

### GSEA Enrichment

GSEA was performed to identify enriched functional pathways in.
Uses **GSEA**, not over-representation analysis (ORA). GSEA ranks *every*
  gene tested in a tissue by `log2FoldChange` and asks whether a gene set is
  enriched toward the top (`NES > 0`, "Upregulated") or bottom (`NES < 0`,
  "Downregulated") of that ranking — no significance cutoff or separate
  up/down gene list is needed. This matters here specifically because an
  ORA approach (hard cutoff + separate up/down lists) left cornea and sclera
  with only 1 significant "down" gene each — too few to test — while GSEA,
  using the full ranked list, still recovered hundreds of suppressed terms
  for both.
**Gene sets:** GO (BP/CC/MF) via `org.Hs.eg.db`; KEGG pathways via live
  query (see [Known caveats](#known-caveats)).
**Redundant GO term removal:** `simplify()` (GOSemSim, Wang semantic
  similarity, cutoff = 0.7 — clusterProfiler's own default) collapses
  near-duplicate terms (e.g. "DNA replication" / "DNA-templated DNA
  replication" / "chromosome segregation" describing the same underlying
  process) down to one representative term per cluster, keeping the most
  significant.
**Reproducibility:** `set.seed(42)` before all `gseGO()`/`gseKEGG()` calls
  — these use permutation-based p-values (`fgsea`) which are otherwise
  non-deterministic between runs.

### Data Analysis
- Significant genes per tissue (padj < 0.05, |log2FC| > 1.5), CoV2 vs. mock.


## Environment Setup

Dependencies are pinned in `renv.lock` (145 packages: R 4.5.3, Bioconductor 3.22).  In the project root, open R:

```r
install.packages("renv")   # if not already installed
renv::restore()
```

renv::restore() recreates the project's package environment using the versions recorded in `renv.lock` without modifying global R package library.

Key package versions:

| Package | Version | Source |
|---|---|---|
| DESeq2 | 1.50.2 | Bioconductor |
| clusterProfiler | 4.18.4 | Bioconductor |
| org.Hs.eg.db | 3.22.0 | Bioconductor |
| enrichplot | 1.30.4 | Bioconductor |
| GOSemSim | 2.36.0 | Bioconductor |
| EnhancedVolcano | 1.28.2 | Bioconductor |
| pheatmap | 1.0.13 | CRAN |
| ggplot2 | 4.0.3 | CRAN |


(Full list: all 145 entries in `renv.lock`, which also includes every dependency of the above.)


## Running the Analysis Workflow

After completing setup using the R envirnoment specified in `renv.lock`, the analysis workflow can be executed in two ways:


### A. Manual Execution

Run `deseq2.R` first, followed by `gsea.R`.

```bash
cd RNA-seq/ # run from the project root

Rscript src/deseq2.R
Rscript src/gsea.R
```

Both scripts accept `--key=value` CLI flags.  Any flag left unset keeps its default value:


`deseq2.R`:

```bash
Rscript scripts/deseq2.R \
  --counts_file=GSE164073_Eye_count_matrix.csv \
  --padj_cutoff=0.05 \
  --lfc_cutoff=1.5 \
  --min_count_sum=0 \
  --tissues=cornea,limbus,sclera
```

`gsea.R`:

```bash
Rscript scripts/gsea.R \
  --tissues=cornea,limbus,sclera \
  --use_simplify=TRUE \
  --simplify_cutoff=0.7
```


### B. Snakemake

Snakemake requires a separate environment from the R project environment managed by `renv`.  To create the environment:

```bash
mamba create -n snakemake_env -c bioconda -c conda-forge \
  --no-channel-priority snakemake-minimal=9.23.1
```

A `Snakefile` at the project root automates execution of the DESeq2 and GSEA steps.  Run from the project root with the `snakemake_env` environment active:

```bash
conda activate snakemake_env

snakemake -n           # dry run
snakemake --cores 3 --config rscript=/full/path/to/Rscript
```


## Directory Structure

```
RNA-seq/
├── data/
│   └── GSE164073_Eye_count_matrix.csv   # input file, count matrix
├── src/
│   ├── deseq2.R                         # differential expression analysis
│   └── gsea.R                           # GSEA enrichment analysis
├── Snakefile                            # optional Snakemake workflow
├── renv.lock                            # pinned R package versions
├── .Rprofile                            # activates renv for this project
├── renv/                                # renv infrastructure
├── results/
│   ├── <tissue>/
│   │   ├── de_tables/                   # DE result tables (.tsv)
│   │   ├── plots/                       # PCA, volcano, MA, heatmaps (.png)
│   │   ├── qc/                          # QC logs (.txt)
│   │   └── gsea/                        # GSEA tables (.tsv) and dotplots (.png)
│   ├── summary_all_tissues.tsv          # DE gene counts across tissues
│   └── gsea_summary.tsv                 # GSEA term counts across tissues
└── report.html                          # interactive report
```


## Output files

Per tissue, in `results/<tissue>/`:

- **`de_tables/`** — `<tissue>_all_genes_CoV2_vs_mock.tsv` (all genes), `<tissue>_significant_CoV2_vs_mock.tsv` (significant genes), `<tissue>_normalized_counts.tsv` (DESeq2-normalized counts)
- **`plots/`** — `pca.png`, `sample_distance_heatmap.png`, `volcano.png`, `ma_plot.png`, `top_de_genes_heatmap.png`
- **`qc/qc_summary.txt`** — filtering stats, size factors, DE gene count
- **`gsea/`** — `<tissue>_GSEA_GO_{BP,CC,MF}.tsv`, `<tissue>_GSEA_KEGG.tsv` result tables, plus a corresponding `_dotplot.png` for each category

`report.html` — interactive summary report containing QC, differential expression, and enrichment analysis plots across tissues.

## Notes


## References

Eriksen AZ, Møller R, Makovoz B, Uhl SA, tenOever BR, Blenkinsop TA. SARS-CoV-2 infects human adult donor eyes and hESC-derived ocular epithelium. Cell Stem Cell. 2021 Jul 1;28(7):1205-1220.e7. doi: 10.1016/j.stem.2021.04.028. Epub 2021 May 17. PMID: 34022129; PMCID: PMC8126605.

