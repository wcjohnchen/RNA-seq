#!/usr/bin/env Rscript
# GSEA (Gene Set Enrichment Analysis): GO (BP/CC/MF) and KEGG.
# GSEA ranks all genes tested in tissues by log2FoldChange. 
# NES (normalized enriched score) > 0 indicates enrichment among upregulated genes. 
# NES < 0 indicates enrichment among downregulated genes.
# GO terms are simplified using simplify() from clusterProfiler,
# which applies GOSemSim semantic similarity analysis (Wang method) to identify and remove redundant GO terms.
# Highly similar GO terms are grouped, and the most significant term is used as the representative term.
# KEGG pathways are not simplified because they are already curated with less redundancy.

suppressMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(ggplot2)
  library(GOSemSim)
})

set.seed(42)

# CLI overrides, e.g.:
# Rscript gsea.R --tissues=cornea,sclera --use_simplify=FALSE --simplify_cutoff=0.8
# any flag left unset keeps its default below
parse_cli_args <- function(defaults) {
  for (arg in commandArgs(trailingOnly = TRUE)) {
    m <- regmatches(arg, regexec("^--([a-zA-Z_]+)=(.*)$", arg))[[1]]
    if (length(m) == 3 && m[2] %in% names(defaults)) defaults[[m[2]]] <- m[3]
  }
  defaults
}

cli <- parse_cli_args(list(
  tissues         = "cornea,limbus,sclera",
  use_simplify    = "TRUE",
  simplify_cutoff = "0.7"
))

tissues <- strsplit(cli$tissues, ",")[[1]]
ontologies <- c("BP", "CC", "MF")
ont_names <- c(BP = "Biological Process", CC = "Cellular Component", MF = "Molecular Function")
use_simplify <- as.logical(cli$use_simplify)
simplify_cutoff <- as.numeric(cli$simplify_cutoff)  # clusterProfiler's own default

if (use_simplify) {
  cat("Precomputing GO semantic similarity data (once per ontology, reused across tissues)...\n")
  sem_data <- lapply(ontologies, function(o) godata("org.Hs.eg.db", ont = o))
  names(sem_data) <- ontologies
}

summary_rows <- list()

sign_labeller <- as_labeller(c(activated = "Upregulated", suppressed = "Downregulated"))

save_gsea_plot <- function(gse, out_path, title) {
  n_show <- min(15, nrow(as.data.frame(gse)))
  p <- dotplot(gse, showCategory = n_show, split = ".sign", label_format = 55) +
    facet_grid(. ~ .sign, labeller = sign_labeller) +
    ggtitle(title) +
    theme(axis.text.y = element_text(size = 16, lineheight = 0.9))
  png(out_path, width = 2600, height = 650 + n_show * 165, res = 150)
  print(p)
  dev.off()
}

capitalize <- function(x) paste0(toupper(substring(x, 1, 1)), substring(x, 2))

for (t in tissues) {
  t_label <- capitalize(t)
  tab <- read.delim(sprintf("results/%s/de_tables/%s_all_genes_CoV2_vs_mock.tsv", t, t))
  tab <- tab[!is.na(tab$log2FoldChange), ]

  gene_list <- tab$log2FoldChange
  names(gene_list) <- tab$gene_id
  gene_list <- sort(gene_list, decreasing = TRUE)

  out_dir <- sprintf("results/%s/gsea", t)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  cat(sprintf("\n==== %s: %d genes ranked by log2FoldChange ====\n", t, length(gene_list)))

  ##### GSEA GO
  for (ont in ontologies) {
    gse <- tryCatch(
      gseGO(geneList = gene_list, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = ont,
            pAdjustMethod = "BH", pvalueCutoff = 0.05, minGSSize = 5, eps = 0, verbose = FALSE),
      error = function(e) { cat("  GO", ont, "error:", conditionMessage(e), "\n"); NULL }
    )
    n_before <- if (is.null(gse)) 0 else nrow(as.data.frame(gse))
    if (use_simplify && !is.null(gse) && n_before > 0) {
      gse <- tryCatch(
        simplify(gse, cutoff = simplify_cutoff, by = "p.adjust", select_fun = min,
                 measure = "Wang", semData = sem_data[[ont]]),
        error = function(e) { cat("  simplify() failed, using unsimplified result:", conditionMessage(e), "\n"); gse }
      )
    }
    df <- if (is.null(gse)) data.frame() else as.data.frame(gse)
    n_up <- sum(df$NES > 0); n_down <- sum(df$NES < 0)
    cat(sprintf("  GO %s (%s): %d terms before simplify -> %d after (%d up / NES>0, %d down / NES<0)\n",
                ont, ont_names[ont], n_before, nrow(df), n_up, n_down))
    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      tissue = t, category = sprintf("GO_%s", ont), n_terms = nrow(df), n_up = n_up, n_down = n_down)

    if (nrow(df) == 0) next
    write.table(df, sprintf("%s/%s_GSEA_GO_%s.tsv", out_dir, t, ont), sep = "\t", quote = FALSE, row.names = FALSE)
    save_gsea_plot(gse, sprintf("%s/%s_GSEA_GO_%s_dotplot.png", out_dir, t, ont),
                    sprintf("%s: GSEA GO %s", t_label, ont_names[ont]))
  }

  ##### GSEA KEGG
  gene_map <- suppressMessages(bitr(names(gene_list), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db))
  gene_map <- gene_map[!duplicated(gene_map$SYMBOL), ]
  merged <- merge(data.frame(SYMBOL = names(gene_list), lfc = gene_list), gene_map, by = "SYMBOL")
  kegg_list <- merged$lfc
  names(kegg_list) <- merged$ENTREZID
  kegg_list <- sort(kegg_list, decreasing = TRUE)

  gsekegg <- tryCatch(
    gseKEGG(geneList = kegg_list, organism = "hsa", pAdjustMethod = "BH",
            pvalueCutoff = 0.05, minGSSize = 5, eps = 0, verbose = FALSE),
    error = function(e) { cat("  KEGG error:", conditionMessage(e), "\n"); NULL }
  )
  dfk <- if (is.null(gsekegg)) data.frame() else as.data.frame(gsekegg)
  n_up <- sum(dfk$NES > 0); n_down <- sum(dfk$NES < 0)
  cat(sprintf("  KEGG: %d terms (%d up / NES>0, %d down / NES<0)\n", nrow(dfk), n_up, n_down))
  summary_rows[[length(summary_rows) + 1]] <- data.frame(
    tissue = t, category = "KEGG", n_terms = nrow(dfk), n_up = n_up, n_down = n_down)

  if (nrow(dfk) > 0) {
    gsekegg_readable <- tryCatch(setReadable(gsekegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID"), error = function(e) gsekegg)
    write.table(as.data.frame(gsekegg_readable), sprintf("%s/%s_GSEA_KEGG.tsv", out_dir, t),
                sep = "\t", quote = FALSE, row.names = FALSE)
    save_gsea_plot(gsekegg, sprintf("%s/%s_GSEA_KEGG_dotplot.png", out_dir, t),
                    sprintf("%s: GSEA KEGG pathways", t_label))
  }
}

summary_df <- do.call(rbind, summary_rows)
write.table(summary_df, "results/gsea_summary.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
cat("\n==== GSEA summary ====\n")
print(summary_df)
cat("\nDone. Per-tissue outputs in results/<tissue>/gsea/\n")
