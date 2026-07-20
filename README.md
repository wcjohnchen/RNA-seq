# SARS-CoV-2 Ocular Tissue Response — RNA-seq Reanalysis

Independent reanalysis of a public RNA-seq dataset examining the transcriptional
response of human ocular tissue to SARS-CoV-2 infection. Differential expression
via **DESeq2**, functional enrichment via **GSEA** (GO Biological Process /
Cellular Component / Molecular Function, and KEGG pathways).

This is a personal reanalysis for learning/portfolio purposes, independent of
the original authors. See [Citation](#citation--data-provenance) below.

📊 **[RNA-seq report](https://claude.ai/code/artifact/821e1224-4e65-40cc-b8d3-4f124b78354c)** — every plot below, in one self-contained page with navigation and click-to-enlarge.

## Contents

- [Data](#data)
- [Methods](#methods)
- [Results summary](#results-summary)
- [Setup](#setup)
- [Running the pipeline](#running-the-pipeline)
- [Directory structure](#directory-structure)
- [Output files](#output-files)
- [Known caveats](#known-caveats)
- [Citation / data provenance](#citation--data-provenance)

---

## Data

- **Source:** GEO accession [GSE164073](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE164073)
- **File:** `GSE164073_Eye_count_matrix.csv` — raw gene-level count matrix,
  27,946 genes × 18 samples
- **Design:** 3 tissues (cornea, limbus, sclera) × 2 conditions (mock,
  SARS-CoV-2-infected) × 3 replicates = 18 samples. Sample columns are named
  `MW<n>_<tissue>_<condition>_<replicate>`, e.g. `MW1_cornea_mock_1`.

If the CSV isn't present in this directory, download it from the GEO
accession above (Supplementary file) and place it at the project root.

## Methods

### Differential expression (`scripts/deseq2.R`)

- **Each tissue is modeled independently** — its own size factors, dispersion
  trend, and contrast — rather than pooling all 18 samples into one
  `~tissue + condition` design. Cornea, limbus, and sclera are different
  enough tissue types that forcing a shared dispersion trend across them
  would be a worse fit than modeling each separately.
- **Gene filtering:** genes with zero counts across *all* samples within a
  given tissue's own 6-sample subset are dropped before fitting (drops
  ~7,900–8,500 of 27,946 genes per tissue). This is a per-tissue filter,
  not a global one — a gene expressed in one tissue but silent in another is
  only dropped from the tissue where it's actually silent.
- **Normalization:** DESeq2 median-of-ratios (size factors) for the DE test
  itself; VST (variance-stabilizing transformation, `blind = TRUE`) is used
  separately, only for visualization (PCA, sample-distance heatmap) — never
  for the statistical test.
- **Test:** Wald test via `DESeq()`, contrast = CoV2 vs. mock.
- **Fold-change shrinkage:** `lfcShrink(..., type = "normal")`.
- **Significance cutoff:** `padj < 0.05` and `|log2FoldChange| > 1.5`
  (CLI-configurable, see below).

### Enrichment — GSEA (`scripts/gsea.R`)

- Uses **GSEA**, not over-representation analysis (ORA). GSEA ranks *every*
  gene tested in a tissue by `log2FoldChange` and asks whether a gene set is
  enriched toward the top (`NES > 0`, "Upregulated") or bottom (`NES < 0`,
  "Downregulated") of that ranking — no significance cutoff or separate
  up/down gene list is needed. This matters here specifically because an
  ORA approach (hard cutoff + separate up/down lists) left cornea and sclera
  with only 1 significant "down" gene each — too few to test — while GSEA,
  using the full ranked list, still recovered hundreds of suppressed terms
  for both.
- **Gene sets:** GO (BP/CC/MF) via `org.Hs.eg.db`; KEGG pathways via live
  query (see [Known caveats](#known-caveats)).
- **Redundant GO term removal:** `simplify()` (GOSemSim, Wang semantic
  similarity, cutoff = 0.7 — clusterProfiler's own default) collapses
  near-duplicate terms (e.g. "DNA replication" / "DNA-templated DNA
  replication" / "chromosome segregation" describing the same underlying
  process) down to one representative term per cluster, keeping the most
  significant.
- **Reproducibility:** `set.seed(42)` before all `gseGO()`/`gseKEGG()` calls
  — these use permutation-based p-values (`fgsea`) which are otherwise
  non-deterministic between runs.

## Results summary

Significant genes per tissue (padj < 0.05, |log2FC| > 1.5), CoV2 vs. mock:

| Tissue | Significant | Up | Down |
|---|---|---|---|
| Cornea | 9 | 8 | 1 |
| Limbus | 16 | 4 | 12 |
| Sclera | 14 | 13 | 1 |

GSEA enriched terms per tissue/category (after `simplify()` for GO):

| Tissue | GO BP | GO CC | GO MF | KEGG |
|---|---|---|---|---|
| Cornea | 334 (211↑ / 123↓) | 84 | 100 | 95 |
| Limbus | 188 (107↑ / 81↓) | 42 | 38 | 71 |
| Sclera | 238 (91↑ / 147↓) | 43 | 45 | 104 |

Biological read: all three tissues show an NF-κB/chemokine/complement
inflammatory signature under infection, consistent with the original study's
reported conclusion — but each tissue leans toward a different dominant
program (cornea: chemokine/inflammatory; limbus: epithelial/transport genes,
possibly tissue-identity rather than infection response; sclera:
complement/oxidative stress — `CFB`, `C3`, `SOD2`). See `report.html` for the
full figure set.

## Setup

Dependencies are pinned in `renv.lock` (145 packages: R 4.5.3, Bioconductor
3.22). From R, in the project root:

```r
install.packages("renv")   # if not already installed
renv::restore()
```

This installs every package at the exact recorded version into an isolated,
project-local library — it does not touch your system/global R library.

**First-time disk cost:** ~2.6 GB (dominated by `org.Hs.eg.db`, the human
gene-annotation database, at ~405 MB). This is a **one-time, per-machine**
cost — `renv` shares a single package cache across all `renv`-managed
projects on the same machine, so a second project needing the same packages
costs nothing extra.

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
| RColorBrewer | 1.1-3 | CRAN |
| gtable | 0.3.6 | CRAN |

(Full list: all 145 entries in `renv.lock`, which also includes every
transitive dependency of the above.)

## Running the pipeline

After [Setup](#setup) (`renv::restore()` — always required, either way),
pick one of two ways to actually run the scripts. Both call the exact same
`deseq2.R`/`gsea.R` underneath and produce identical results — Snakemake is
purely a convenience layer for automatic run-ordering and only-rerunning
what changed, not a different computation path.

- **Manual** — run the two scripts directly. No extra install needed.
- **Snakemake** — install Snakemake separately, then let it handle
  run-order and incremental re-runs automatically.

### Manual

Run in this order — `gsea.R` reads `deseq2.R`'s output tables.

```bash
cd GSE164073_Eye_count_matrix/     # must run from the project root (relative paths)

Rscript scripts/deseq2.R
Rscript scripts/gsea.R
```

Both scripts accept `--key=value` CLI flags; any flag left unset keeps its
default (shown below).

**`deseq2.R`:**

```bash
Rscript scripts/deseq2.R \
  --counts_file=GSE164073_Eye_count_matrix.csv \
  --padj_cutoff=0.05 \
  --lfc_cutoff=1.5 \
  --min_count_sum=0 \
  --tissues=cornea,limbus,sclera
```

**`gsea.R`:**

```bash
Rscript scripts/gsea.R \
  --tissues=cornea,limbus,sclera \
  --use_simplify=TRUE \
  --simplify_cutoff=0.7
```

Example — loosen the fold-change cutoff and run only one tissue:

```bash
Rscript scripts/deseq2.R --lfc_cutoff=1 --tissues=cornea
```

### Snakemake

A `Snakefile` at the project root wraps `deseq2.R` and `gsea.R` **unmodified**
— it calls them through the `--tissues=<tissue>` flag they already support,
so each tissue runs as an independent job. Snakemake then handles run-order
automatically (instead of you remembering to run `deseq2.R` before `gsea.R`)
and only re-runs the specific tissue(s) whose inputs actually changed.

Snakemake needs its own environment, separate from the one `renv.lock`
manages — it's a workflow tool, not an R package, and a recent version's
dependencies conflicted with packages already pinned for the R/Bioconductor
side. `--no-channel-priority` is required here because this machine's
`.condarc` sets `channel_priority: strict`, which otherwise makes the
solver fail with spurious dependency conflicts:

```bash
mamba create -n snakemake_env -c bioconda -c conda-forge \
  --no-channel-priority snakemake-minimal=9.23.1
```

Run from the project root, with `snakemake_env` active:

```bash
conda activate snakemake_env

snakemake -n           # dry run — show what would execute
snakemake --cores 3    # run for real, up to 3 tissues in parallel
```

`RSCRIPT` in the `Snakefile` defaults to plain `Rscript` — this works as
long as it's reachable on `PATH` from wherever `snakemake` runs, relying on
the project's `.Rprofile` to auto-activate `renv` (and therefore the
correct package versions) whenever R starts here. If Snakemake lives in a
separate environment with no R installed in it at all (as `snakemake_env`
deliberately does here, to avoid the dependency conflict above), override
the path at the command line instead of editing the `Snakefile`:

```bash
snakemake --cores 3 --config rscript=/full/path/to/Rscript
```

The two cross-tissue summary files (`results/summary_all_tissues.tsv`,
`results/gsea_summary.tsv`) aren't produced by the R scripts in this mode —
each script writes its summary by overwriting one shared file, which works
when it processes all tissues in a single call but would clobber itself if
Snakemake invoked it three separate times. A dedicated
`aggregate_summaries` rule in the `Snakefile` rebuilds both summaries in
plain Python directly from each tissue's already-correct output files
instead, after all three tissues finish — no R code changes needed for
this either.

**Caveat:** GSEA term counts can show minor (single-term, boundary-case)
variation between runs even with `set.seed(42)` — most likely `fgsea`'s
internal parallelism consuming random numbers in a run-order-dependent way
for genes sitting right at the significance threshold. BP, CC, and KEGG
have reproduced exactly across every run so far; GO MF has differed by
exactly one term in isolated cases. This is a `fgsea`/`clusterProfiler`
limitation, not something the pipeline or Snakemake wrapper introduces.

## Directory structure

```
GSE164073_Eye_count_matrix/
├── GSE164073_Eye_count_matrix.csv   # input count matrix (from GEO)
├── renv.lock                        # pinned package versions
├── .Rprofile                        # auto-activates renv for this project
├── renv/                            # renv infrastructure (library/ is gitignored)
├── Snakefile                        # optional Snakemake orchestration (see below)
├── scripts/
│   ├── deseq2.R                     # DE analysis (run first)
│   └── gsea.R                       # GSEA enrichment (run second)
├── results/
│   ├── <tissue>/de_tables/          # DE result tables (.tsv)
│   ├── <tissue>/plots/              # PCA, volcano, MA, heatmaps (.png)
│   ├── <tissue>/qc/                 # QC log per tissue
│   ├── <tissue>/gsea/               # GSEA tables + dotplots per ontology
│   ├── summary_all_tissues.tsv      # DE gene counts, all tissues
│   └── gsea_summary.tsv             # GSEA term counts, all tissues
├── report.html                      # standalone interactive report (all plots embedded)
└── report_link.txt                  # hosted copy of report.html, if published
```

## Output files

Per tissue, in `results/<tissue>/`:

- **`de_tables/`** — `<tissue>_all_genes_CoV2_vs_mock.tsv` (every tested
  gene), `<tissue>_significant_CoV2_vs_mock.tsv` (filtered to the cutoff
  above), `<tissue>_normalized_counts.tsv` (DESeq2-normalized counts)
- **`plots/`** — `pca.png`, `sample_distance_heatmap.png`, `volcano.png`,
  `ma_plot.png`, `top_de_genes_heatmap.png`
- **`qc/qc_summary.txt`** — filtering stats, size factors, DE gene count
- **`gsea/`** — `<tissue>_GSEA_GO_{BP,CC,MF}.tsv` /
  `<tissue>_GSEA_KEGG.tsv` result tables, plus a matching `_dotplot.png`
  for each category

`report.html` bundles every plot above into one self-contained,
navigable page (sidebar per tissue/section, click-to-enlarge, stat-tile
overview) — open it directly in a browser, no server needed.

## Known caveats

- **KEGG requires live internet access.** `gseKEGG()` queries KEGG's REST
  API on every run rather than using a bundled offline copy — Bioconductor's
  old offline `KEGG.db` package was deprecated after KEGG restricted bulk
  redistribution of their database around 2011; free access remains
  available via their API for this kind of per-run academic use, just not
  as a redistributable static file. If the machine running this has no
  internet access, GO results still complete normally (they use the local
  `org.Hs.eg.db` package), but KEGG results for that run will be empty
  (the `gseKEGG()` call is wrapped in `tryCatch()`, so this doesn't crash
  the script). Because of this, exact KEGG term counts aren't pinned the
  way GO/DESeq2 results are — they reflect KEGG's database as of whenever
  the script was last run, not a fixed snapshot.
- **GSEA term counts can still differ from the numbers in this README** if
  re-run without `renv::restore()` first (i.e. on different package
  versions than pinned) — `fgsea`'s permutation test is seeded (`42`), so
  results are reproducible *given identical package versions*, but the
  seed doesn't guarantee identical output across different `clusterProfiler`
  or `fgsea` releases.
- **`deseq2.R` must be run before `gsea.R`** — the latter reads
  the former's `de_tables/*.tsv` output directly.

## Citation / data provenance

Original data and study:

> Eriksen AZ, Møller R, Makovoz B, Uhl SA, tenOever BR, Blenkinsop TA.
> **SARS-CoV-2 infects human adult donor eyes and hESC-derived ocular
> epithelium.** *Cell Stem Cell.* 2021 Jul 1;28(7):1205-1220.e7.
> doi: [10.1016/j.stem.2021.04.028](https://doi.org/10.1016/j.stem.2021.04.028).
> PMID: [34022129](https://pubmed.ncbi.nlm.nih.gov/34022129/).

Raw data: GEO [GSE164073](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE164073).

All analysis in this repository (DESeq2 pipeline, GSEA, plots, this
README) is an independent reanalysis and is **not** the original authors'
code or figures. The original study's own enrichment analysis used Enrichr
with a top-100-differentially-expressed-gene list; this repository instead
uses DESeq2 with an explicit significance cutoff for DE calling and GSEA
(full ranked list, no cutoff) for enrichment — a different methodology, so
exact term-level results are not expected to match the original paper
one-to-one, though the broad biological conclusion (NF-κB-mediated
chemokine/inflammatory response, interferon-stimulated genes) does
reproduce independently.
