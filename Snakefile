# Snakemake workflow for the RNA-seq DE + GSEA pipeline.
# Wraps src/deseq2.R, src/gsea.R, and src/generate_report.R
#
# Uses the same R environment as manual execution (README section 3):
# Snakemake invokes plain `Rscript`, which auto-activates renv via
# .Rprofile since every rule runs with cwd = project root -- no separate
# conda env needed for R packages.
#
# onstart below runs `renv::restore()` twice before any rule executes, so
# a fresh clone's R packages get installed automatically. Twice, not once,
# because Bioconductor's deep dependency chains (e.g. DESeq2 -> ... ->
# DelayedArray, enrichplot -> RSQLite) can hit a parallel-install ordering
# race on the first pass, where a package finishes building just before
# one of its own dependencies does and fails; the second pass installs
# only whatever failed, and by then its dependencies already exist.
# Confirmed live: first pass failed on DESeq2 + enrichplot, second pass
# installed both cleanly with nothing else to redo (everything already
# correct was skipped, not reinstalled).
#
# Usage:
#   snakemake -n                 # dry run
#   snakemake --cores 3          # real run
#   snakemake --cores 3 --config rscript=/full/path/to/Rscript  # override which Rscript to use

RSCRIPT = config.get("rscript", "Rscript")

onstart:
    shell(RSCRIPT + " -e 'renv::restore(prompt=FALSE)' && " + RSCRIPT + " -e 'renv::restore(prompt=FALSE)'")

TISSUES = ["cornea", "limbus", "sclera"]
ONTOLOGIES = ["BP", "CC", "MF"]

DESEQ2_PLOTS = ["pca.png", "sample_distance_heatmap.png", "volcano.png",
                "ma_plot.png", "top_de_genes_heatmap.png"]
GSEA_DOTPLOTS = ["GSEA_GO_BP_dotplot.png", "GSEA_GO_CC_dotplot.png",
                 "GSEA_GO_MF_dotplot.png", "GSEA_KEGG_dotplot.png"]

rule all:
    input:
        expand("results/{t}/de_tables/{t}_all_genes_CoV2_vs_mock.tsv", t=TISSUES),
        expand("results/{t}/de_tables/{t}_significant_CoV2_vs_mock.tsv", t=TISSUES),
        expand("results/{t}/de_tables/{t}_normalized_counts.tsv", t=TISSUES),
        expand("results/{t}/plots/{plot}", t=TISSUES, plot=DESEQ2_PLOTS),
        expand("results/{t}/qc/qc_summary.txt", t=TISSUES),
        expand("results/{t}/gsea/{t}_GSEA_GO_{ont}.tsv", t=TISSUES, ont=ONTOLOGIES),
        expand("results/{t}/gsea/{t}_GSEA_KEGG.tsv", t=TISSUES),
        expand("results/{t}/gsea/{t}_{plot}", t=TISSUES, plot=GSEA_DOTPLOTS),
        "results/summary_all_tissues.tsv",
        "results/gsea_summary.tsv",
        "report.html"

rule deseq2:
    input:
        counts = "data/GSE164073_Eye_count_matrix.csv",
        script = "src/deseq2.R"
    output:
        "results/{tissue}/de_tables/{tissue}_all_genes_CoV2_vs_mock.tsv",
        "results/{tissue}/de_tables/{tissue}_significant_CoV2_vs_mock.tsv",
        "results/{tissue}/de_tables/{tissue}_normalized_counts.tsv",
        expand("results/{{tissue}}/plots/{plot}", plot=DESEQ2_PLOTS),
        "results/{tissue}/qc/qc_summary.txt"
    shell:
        RSCRIPT + " {input.script} --tissues={wildcards.tissue} --counts_file={input.counts}"

rule gsea:
    input:
        de_table = "results/{tissue}/de_tables/{tissue}_all_genes_CoV2_vs_mock.tsv",
        script = "src/gsea.R"
    output:
        expand("results/{{tissue}}/gsea/{{tissue}}_GSEA_GO_{ont}.tsv", ont=ONTOLOGIES),
        "results/{tissue}/gsea/{tissue}_GSEA_KEGG.tsv",
        expand("results/{{tissue}}/gsea/{{tissue}}_{plot}", plot=GSEA_DOTPLOTS)
    shell:
        RSCRIPT + " {input.script} --tissues={wildcards.tissue}"

rule aggregate_summaries:
    input:
        sig_tables = expand("results/{t}/de_tables/{t}_significant_CoV2_vs_mock.tsv", t=TISSUES),
        gsea_go = expand("results/{t}/gsea/{t}_GSEA_GO_{ont}.tsv", t=TISSUES, ont=ONTOLOGIES),
        gsea_kegg = expand("results/{t}/gsea/{t}_GSEA_KEGG.tsv", t=TISSUES)
    output:
        de_summary = "results/summary_all_tissues.tsv",
        gsea_summary = "results/gsea_summary.tsv"
    run:
        import csv

        def read_tsv(path):
            with open(path) as f:
                return list(csv.DictReader(f, delimiter="\t"))

        # ---- DE summary: one row per tissue ----
        with open(output.de_summary, "w", newline="") as f:
            w = csv.writer(f, delimiter="\t")
            w.writerow(["tissue", "n_significant", "n_up", "n_down"])
            for t, path in zip(TISSUES, input.sig_tables):
                rows = read_tsv(path)
                n_sig = len(rows)
                n_up = sum(1 for r in rows if float(r["log2FoldChange"]) > 0)
                n_down = sum(1 for r in rows if float(r["log2FoldChange"]) < 0)
                w.writerow([t, n_sig, n_up, n_down])

        # ---- GSEA summary: one row per tissue x category ----
        go_paths = {(t, ont): p for t in TISSUES for ont in ONTOLOGIES
                    for p in input.gsea_go if f"/{t}/" in p and f"_GO_{ont}.tsv" in p}
        kegg_paths = {t: p for t in TISSUES for p in input.gsea_kegg if f"/{t}/" in p}

        with open(output.gsea_summary, "w", newline="") as f:
            w = csv.writer(f, delimiter="\t")
            w.writerow(["tissue", "category", "n_terms", "n_up", "n_down"])
            for t in TISSUES:
                for ont in ONTOLOGIES:
                    rows = read_tsv(go_paths[(t, ont)])
                    n_up = sum(1 for r in rows if float(r["NES"]) > 0)
                    n_down = sum(1 for r in rows if float(r["NES"]) < 0)
                    w.writerow([t, f"GO_{ont}", len(rows), n_up, n_down])
                rows = read_tsv(kegg_paths[t])
                n_up = sum(1 for r in rows if float(r["NES"]) > 0)
                n_down = sum(1 for r in rows if float(r["NES"]) < 0)
                w.writerow([t, "KEGG", len(rows), n_up, n_down])

rule generate_report:
    input:
        template = "src/report_template.html",
        script = "src/generate_report.R",
        de_summary = "results/summary_all_tissues.tsv",
        gsea_summary = "results/gsea_summary.tsv",
        norm_counts = expand("results/{t}/de_tables/{t}_normalized_counts.tsv", t=TISSUES),
        qc = expand("results/{t}/qc/qc_summary.txt", t=TISSUES),
        deseq2_plots = expand("results/{t}/plots/{plot}", t=TISSUES, plot=DESEQ2_PLOTS),
        gsea_dotplots = expand("results/{t}/gsea/{t}_{plot}", t=TISSUES, plot=GSEA_DOTPLOTS)
    output:
        "report.html"
    shell:
        RSCRIPT + " {input.script}"
