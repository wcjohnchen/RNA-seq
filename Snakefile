# Snakemake workflow for the RNA-seq DE + GSEA pipeline.
# Wraps src/deseq2.R and src/gsea.R
# Usage:
#   snakemake -n                 # dry run
#   snakemake --cores 3          # real run
#   snakemake --config rscript=/full/path/to/Rscript --cores 3

RSCRIPT = config.get("rscript", "Rscript")

TISSUES = ["cornea", "limbus", "sclera"]
ONTOLOGIES = ["BP", "CC", "MF"]

DESEQ2_PLOTS = ["pca.png", "sample_distance_heatmap.png", "volcano.png",
                "ma_plot.png", "top_de_genes_heatmap.png"]

rule all:
    input:
        expand("results/{t}/de_tables/{t}_all_genes_CoV2_vs_mock.tsv", t=TISSUES),
        expand("results/{t}/de_tables/{t}_significant_CoV2_vs_mock.tsv", t=TISSUES),
        expand("results/{t}/de_tables/{t}_normalized_counts.tsv", t=TISSUES),
        expand("results/{t}/plots/{plot}", t=TISSUES, plot=DESEQ2_PLOTS),
        expand("results/{t}/qc/qc_summary.txt", t=TISSUES),
        expand("results/{t}/gsea/{t}_GSEA_GO_{ont}.tsv", t=TISSUES, ont=ONTOLOGIES),
        expand("results/{t}/gsea/{t}_GSEA_KEGG.tsv", t=TISSUES),
        "results/summary_all_tissues.tsv",
        "results/gsea_summary.tsv"

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
        "results/{tissue}/gsea/{tissue}_GSEA_KEGG.tsv"
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
