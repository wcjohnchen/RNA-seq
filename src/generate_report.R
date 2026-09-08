#!/usr/bin/env Rscript
# Regenerates report.html from src/report_template.html, using each tissue's actual results/ files

suppressMessages({
  library(base64enc)
})

TISSUES <- c("cornea", "limbus", "sclera")
template_path <- "src/report_template.html"
output_path   <- "report.html"

log_msg <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), sprintf(...), "\n", sep = "")
}

fmt_int <- function(x) format(round(x), big.mark = ",", scientific = FALSE, trim = TRUE)

##### Step 1/3: read aggregate summaries

log_msg("Step 1/3: Reading results/ files...")

de_summary   <- read.delim("results/summary_all_tissues.tsv")
gsea_summary <- read.delim("results/gsea_summary.tsv")

tokens <- list()

for (t in TISSUES) {
  de_row <- de_summary[de_summary$tissue == t, ]
  stopifnot(nrow(de_row) == 1)
  tokens[[paste0("SIG_", t)]]  <- fmt_int(de_row$n_significant)
  tokens[[paste0("UP_", t)]]   <- fmt_int(de_row$n_up)
  tokens[[paste0("DOWN_", t)]] <- fmt_int(de_row$n_down)

  gsea_cat <- function(category) {
    row <- gsea_summary[gsea_summary$tissue == t & gsea_summary$category == category, ]
    stopifnot(nrow(row) == 1)
    row
  }
  bp <- gsea_cat("GO_BP")
  tokens[[paste0("GOBP_", t)]]      <- fmt_int(bp$n_terms)
  tokens[[paste0("GOBP_UP_", t)]]   <- fmt_int(bp$n_up)
  tokens[[paste0("GOBP_DOWN_", t)]] <- fmt_int(bp$n_down)
  tokens[[paste0("GOCC_", t)]] <- fmt_int(gsea_cat("GO_CC")$n_terms)
  tokens[[paste0("GOMF_", t)]] <- fmt_int(gsea_cat("GO_MF")$n_terms)
  tokens[[paste0("KEGG_", t)]] <- fmt_int(gsea_cat("KEGG")$n_terms)


  norm_counts_path <- file.path("results", t, "de_tables", sprintf("%s_normalized_counts.tsv", t))
  header <- colnames(read.delim(norm_counts_path, nrows = 0, check.names = FALSE))
  sample_cols <- setdiff(header, "gene_id")
  tokens[[paste0("N_MOCK_", t)]] <- fmt_int(sum(grepl("_mock_", sample_cols)))
  tokens[[paste0("N_COV2_", t)]] <- fmt_int(sum(grepl("_CoV2_", sample_cols)))


  qc_path <- file.path("results", t, "qc", "qc_summary.txt")
  qc_lines <- readLines(qc_path, warn = FALSE)
  pca_line <- grep("^PCA variance:", qc_lines, value = TRUE)
  if (length(pca_line) == 0) {
    stop(sprintf(
      "No 'PCA variance:' line in %s -- rerun src/deseq2.R for this tissue first (older qc_summary.txt files predate that log line).",
      qc_path
    ))
  }
  pc1 <- sub(".*PC1=([0-9]+)%.*", "\\1", pca_line[1])
  tokens[[paste0("PC1_VAR_", t)]] <- pc1

  log_msg("  %s: %s significant (%s up, %s down) | GO_BP %s terms | PC1=%s%%",
          t, tokens[[paste0("SIG_", t)]], tokens[[paste0("UP_", t)]],
          tokens[[paste0("DOWN_", t)]], tokens[[paste0("GOBP_", t)]], pc1)
}

##### Step 2/3: encode figures

log_msg("Step 2/3: Encoding figures...")

image_map <- list()
for (t in TISSUES) {
  image_map[[paste0("IMG_", t, "_pca")]]                       <- file.path("results", t, "plots", "pca.png")
  image_map[[paste0("IMG_", t, "_sample_distance_heatmap")]]   <- file.path("results", t, "plots", "sample_distance_heatmap.png")
  image_map[[paste0("IMG_", t, "_volcano")]]                   <- file.path("results", t, "plots", "volcano.png")
  image_map[[paste0("IMG_", t, "_ma_plot")]]                   <- file.path("results", t, "plots", "ma_plot.png")
  image_map[[paste0("IMG_", t, "_top_de_genes_heatmap")]]      <- file.path("results", t, "plots", "top_de_genes_heatmap.png")
  image_map[[paste0("IMG_", t, "_GSEA_GO_BP")]]  <- file.path("results", t, "gsea", sprintf("%s_GSEA_GO_BP_dotplot.png", t))
  image_map[[paste0("IMG_", t, "_GSEA_GO_CC")]]  <- file.path("results", t, "gsea", sprintf("%s_GSEA_GO_CC_dotplot.png", t))
  image_map[[paste0("IMG_", t, "_GSEA_GO_MF")]]  <- file.path("results", t, "gsea", sprintf("%s_GSEA_GO_MF_dotplot.png", t))
  image_map[[paste0("IMG_", t, "_GSEA_KEGG")]]   <- file.path("results", t, "gsea", sprintf("%s_GSEA_KEGG_dotplot.png", t))
}

for (name in names(image_map)) {
  path <- image_map[[name]]
  if (!file.exists(path)) stop(sprintf("Missing figure: %s", path))
  tokens[[name]] <- base64enc::base64encode(path)
}
log_msg("Encoded %d figures", length(image_map))

##### Step 3/3: fill the template and write report.html

log_msg("Step 3/3: Filling template and writing report.html...")

html <- paste(readLines(template_path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")

for (name in names(tokens)) {
  html <- gsub(paste0("{{", name, "}}"), tokens[[name]], html, fixed = TRUE)
}

remaining <- regmatches(html, gregexpr("\\{\\{[A-Za-z0-9_]+\\}\\}", html))[[1]]
if (length(remaining) > 0) {
  stop("Unfilled template token(s) remain: ", paste(unique(remaining), collapse = ", "))
}

writeLines(html, output_path, useBytes = TRUE)
log_msg("Saved report.html")
