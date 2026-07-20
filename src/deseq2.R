#!/usr/bin/env Rscript
# DE analysis using DESeq2: CoV2 vs mock in cornea, limbus, and sclera.

suppressMessages({
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
  library(EnhancedVolcano)
  library(RColorBrewer)
  library(gtable)
})

draw_pheatmap_matrix_border <- function(ph, border_color = "black", border_lwd = 3) {
  idx <- which(ph$gtable$layout$name == "matrix")
  pos <- ph$gtable$layout[idx, ]
  border_grob <- grid::rectGrob(gp = grid::gpar(col = border_color, fill = NA, lwd = border_lwd))
  ph$gtable <- gtable::gtable_add_grob(ph$gtable, border_grob, t = pos$t, l = pos$l, b = pos$b, r = pos$r,
                                        name = "matrix_border", z = Inf)
  grid::grid.newpage()
  grid::grid.draw(ph$gtable)
}

# CLI overrides, e.g.:
# Rscript deseq2.R --padj_cutoff=0.01 --lfc_cutoff=1 --tissues=cornea,sclera
# any flag left unset keeps its default below
parse_cli_args <- function(defaults) {
  for (arg in commandArgs(trailingOnly = TRUE)) {
    m <- regmatches(arg, regexec("^--([a-zA-Z_]+)=(.*)$", arg))[[1]]
    if (length(m) == 3 && m[2] %in% names(defaults)) defaults[[m[2]]] <- m[3]
  }
  defaults
}

cli <- parse_cli_args(list(
  counts_file   = "GSE164073_Eye_count_matrix.csv",
  padj_cutoff   = "0.05",
  lfc_cutoff    = "1.5",
  min_count_sum = "0",
  tissues       = "cornea,limbus,sclera"
))

counts_file   <- cli$counts_file
padj_cutoff   <- as.numeric(cli$padj_cutoff)
lfc_cutoff    <- as.numeric(cli$lfc_cutoff)  # log2FC cutoff: > +lfc_cutoff up, < -lfc_cutoff down
min_count_sum <- as.numeric(cli$min_count_sum)
tissues       <- strsplit(cli$tissues, ",")[[1]]

custom_ma_plot <- function(res_df, alpha, col_sig = "navyblue",
                            col_nonsig = "gray60", col_line = "red", main = "") {
  keep <- res_df$baseMean != 0 & !is.na(res_df$log2FoldChange)
  mean_vals <- res_df$baseMean[keep]
  lfc_vals <- res_df$log2FoldChange[keep]
  is_sig <- !is.na(res_df$padj[keep]) & res_df$padj[keep] < alpha

  log2_mean <- log2(mean_vals)
  ylim <- c(-1, 1) * quantile(abs(lfc_vals[is.finite(lfc_vals)]), probs = 0.99, na.rm = TRUE) * 1.1
  pch_vals <- ifelse(lfc_vals < ylim[1], 6, ifelse(lfc_vals > ylim[2], 2, 1)) # triangle capped points are used for genes whose log2FC are outside the plotting range.
  plot(log2_mean, pmax(ylim[1], pmin(ylim[2], lfc_vals)),
       pch = pch_vals, cex = 0.45,
       col = ifelse(is_sig, col_sig, col_nonsig),
       xlab = "log2(mean of normalized counts)", ylab = "log fold change",
       ylim = ylim, main = main)
  abline(h = 0, lwd = 1, col = col_line)
}

##### load count matrix
counts_all <- read.csv(counts_file, row.names = 1, check.names = FALSE)
counts_all <- as.matrix(counts_all)
storage.mode(counts_all) <- "integer"

sample_names <- colnames(counts_all)
parts <- strsplit(sample_names, "_")
tissue_vec    <- vapply(parts, `[`, character(1), 2)
condition_vec <- vapply(parts, `[`, character(1), 3)

metadata_all <- data.frame(sample_id = sample_names,
                            tissue = tissue_vec,
                            condition = condition_vec,
                            row.names = sample_names)

cat(sprintf("Loaded %d genes x %d samples\n", nrow(counts_all), ncol(counts_all)))
print(table(metadata_all$tissue, metadata_all$condition))

capitalize <- function(x) paste0(toupper(substring(x, 1, 1)), substring(x, 2))

analyze_tissue <- function(tissue_name) {
  tissue_label <- capitalize(tissue_name)
  cat(sprintf("\n==== %s: mock vs CoV2 ====\n", tissue_name))

  out_de   <- file.path("results", tissue_name, "de_tables")
  out_plot <- file.path("results", tissue_name, "plots")
  out_qc   <- file.path("results", tissue_name, "qc")

  qc_log <- character()
  log_msg <- function(...) {
    msg <- sprintf(...)
    cat(msg, "\n")
    qc_log <<- c(qc_log, msg)
  }

  samples <- rownames(metadata_all)[metadata_all$tissue == tissue_name]
  counts_sub <- counts_all[, samples]
  meta_sub <- metadata_all[samples, ]
  meta_sub$condition <- factor(meta_sub$condition, levels = c("mock", "CoV2"))

  keep <- rowSums(counts_sub) > min_count_sum
  counts_filt <- counts_sub[keep, ]
  log_msg("Kept %d / %d genes after filtering (rowSum > %d, i.e. drop all-zero genes only)",
          sum(keep), length(keep), min_count_sum)

  dds <- DESeqDataSetFromMatrix(countData = counts_filt,
                                 colData = meta_sub,
                                 design = ~condition)
  dds <- DESeq(dds)
  log_msg("Size factors: %s",
          paste(sprintf("%s=%.3f", colnames(dds), sizeFactors(dds)), collapse = ", "))

  res <- results(dds, contrast = c("condition", "CoV2", "mock"), alpha = padj_cutoff)
  res_shrunk <- lfcShrink(dds, contrast = c("condition", "CoV2", "mock"),
                           res = res, type = "normal")

  res_df <- as.data.frame(res_shrunk)
  res_df$gene_id <- rownames(res_df)
  res_df <- res_df[, c("gene_id", setdiff(colnames(res_df), "gene_id"))]
  res_df <- res_df[order(res_df$padj), ]

  sig_df <- subset(res_df, !is.na(padj) & padj < padj_cutoff & abs(log2FoldChange) > lfc_cutoff)

  write.table(res_df, file.path(out_de, sprintf("%s_all_genes_CoV2_vs_mock.tsv", tissue_name)),
              sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(sig_df, file.path(out_de, sprintf("%s_significant_CoV2_vs_mock.tsv", tissue_name)),
              sep = "\t", quote = FALSE, row.names = FALSE)

  norm_counts <- counts(dds, normalized = TRUE)
  write.table(data.frame(gene_id = rownames(norm_counts), norm_counts, check.names = FALSE),
              file.path(out_de, sprintf("%s_normalized_counts.tsv", tissue_name)),
              sep = "\t", quote = FALSE, row.names = FALSE)

  log_msg("DE testing complete: %d significant genes (padj < %.2f, |log2FC| > %.1f)",
          nrow(sig_df), padj_cutoff, lfc_cutoff)

  ##### plot
  vsd <- vst(dds, blind = TRUE)

  pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
  pct_var <- round(100 * attr(pca_data, "percentVar"))
  p_pca <- ggplot(pca_data, aes(PC1, PC2, color = condition)) +
    geom_point(size = 3, shape = 1, stroke = 1.2) +
    xlab(sprintf("PC1: %d%% variance", pct_var[1])) +
    ylab(sprintf("PC2: %d%% variance", pct_var[2])) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    ggtitle(sprintf("%s: PCA of VST-normalized counts", tissue_label))
  png(file.path(out_plot, "pca.png"), width = 900, height = 750, res = 150)
  print(p_pca)
  dev.off()

  sample_dists <- dist(t(assay(vsd)))
  dist_mat <- as.matrix(sample_dists)
  ph_dist <- pheatmap(dist_mat,
                       annotation_col = meta_sub["condition"],
                       color = colorRampPalette(rev(brewer.pal(9, "Blues")))(255),
                       border_color = NA,
                       cluster_rows = FALSE,
                       cluster_cols = FALSE,
                       silent = TRUE)
  png(file.path(out_plot, "sample_distance_heatmap.png"), width = 6, height = 5, units = "in", res = 150)
  draw_pheatmap_matrix_border(ph_dist)
  dev.off()

  png(file.path(out_plot, "ma_plot.png"), width = 1200, height = 900, res = 150)
  custom_ma_plot(as.data.frame(res_shrunk), alpha = padj_cutoff, col_sig = "navyblue",
                  main = sprintf("%s: MA plot (CoV2 vs mock)", tissue_label))
  dev.off()

  #  -log10(padj) values are capped at 100 to prevent extremely small values from stretching the y-axis.
  y_axis_cap <- 100
  res_df$padj_capped <- pmax(res_df$padj, 10^(-y_axis_cap))

  p_volcano <- EnhancedVolcano(res_df,
                                lab = res_df$gene_id,
                                x = "log2FoldChange",
                                y = "padj_capped",
                                ylab = bquote(~-Log[10]~italic(P[adj])),
                                ylim = c(0, y_axis_cap),
                                pCutoff = padj_cutoff,
                                FCcutoff = lfc_cutoff,
                                title = tissue_label,
                                subtitle = "CoV2 vs mock",
                                drawConnectors = TRUE,
                                widthConnectors = 0.4,
                                min.segment.length = 0,
                                max.overlaps = Inf,
                                gridlines.major = FALSE,
                                gridlines.minor = FALSE)

  png(file.path(out_plot, "volcano.png"), width = 1400, height = 1400, res = 150)
  print(p_volcano)
  dev.off()

  top_genes <- head(sig_df$gene_id[order(sig_df$padj)], 30)
  if (length(top_genes) >= 2) {
    mat <- assay(vsd)[top_genes, ]
    mat <- mat - rowMeans(mat)
    ph_top <- pheatmap(mat,
                        annotation_col = meta_sub["condition"],
                        show_rownames = length(top_genes) <= 40,
                        border_color = NA,
                        cluster_rows = TRUE,
                        cluster_cols = FALSE,
                        silent = TRUE)
    png(file.path(out_plot, "top_de_genes_heatmap.png"), width = 7, height = 8, units = "in", res = 150)
    draw_pheatmap_matrix_border(ph_top)
    dev.off()
  } else {
    log_msg("Skipped top-gene heatmap: fewer than 2 significant genes")
  }

  writeLines(qc_log, file.path(out_qc, "qc_summary.txt"))
  data.frame(tissue = tissue_name, n_significant = nrow(sig_df),
             n_up = sum(sig_df$log2FoldChange > 0), n_down = sum(sig_df$log2FoldChange < 0))
}

summary_rows <- lapply(tissues, analyze_tissue)
summary_df <- do.call(rbind, summary_rows)
write.table(summary_df, "results/summary_all_tissues.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
cat("\n==== Summary across tissues ====\n")
print(summary_df)
cat("\nDone. Per-tissue outputs in results/<tissue>/{de_tables,plots,qc}\n")
