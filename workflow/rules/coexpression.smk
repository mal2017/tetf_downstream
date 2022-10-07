rule filter_models:
    """
    merged: all models -> mean of replicates
    filtered: merged -> all reps significant and no sign switches
    extreme: filtered -> filter at some extreme quantile of abs(coefs) (see script for exact choice)
    """
    input:
        expand("data/linear_models/{m}/{r}/lm.tidy.corrected.tsv.gz",m=["male_model_01","female_model_01"],r=[0,1,2])
    output:
        merged_tsv =  "results/analysis/coexpression/merged_models.tsv.gz",
        filtered_tsv =  "results/analysis/coexpression/filtered_models.tsv.gz",
        extreme_tsv =  "results/analysis/coexpression/extreme_models.tsv.gz"
    script:
        "../scripts/coexpression/consensus_lms.R"

rule plot_intro_heatmap:
    input:
        lkup = rules.make_gene_symbol_lookup.output.tsv,
        filtered_mods = rules.filter_models.output.filtered_tsv,
        extreme_mods = rules.filter_models.output.extreme_tsv,
        merged_mods = rules.filter_models.output.merged_tsv,
    output:
        heats = "results/plots/intro_heatmaps.gg-list.rds",
        gene_hcs = "results/analysis/coexpression/coex-clust.genes.hc-list.rds",
        te_hcs = "results/analysis/coexpression/coex-clust.tes.hc-list.rds",
    script:
        "../scripts/coexpression/plot_intro_heatmap.R"

rule plot_intro_upsetplot:
    input:
        mods = rules.filter_models.output.filtered_tsv,
    output:
        rds="results/plots/plot_intro_upsetplot.rds"
    script:
        "../scripts/coexpression/plot_intro_upsetplot.R"

rule plot_intro_histograms:
    input:
        mods = rules.filter_models.output.filtered_tsv,
    output:
        rds="results/plots/plot_intro_histograms.rds"
    script:
        "../scripts/coexpression/plot_intro_histograms.R"

rule plot_intro_ncoex_scatter:
    input:
        mods = rules.filter_models.output.filtered_tsv,
    output:
        rds="results/plots/plot_intro_ncoex_scatter.rds"
    script:
        "../scripts/coexpression/plot_intro_ncoex_scatter.R"

rule plot_piRNA_genes_in_lms:
    input:
        mods = rules.filter_models.output.filtered_tsv,
        pirna = rules.make_pirna_gene_list.output.tsv,
    output:
        volcano = "results/plots/plot_pirna_genes_in_lms.volc.rds",
        easy_bar = "results/plots/plot_pirna_genes_in_lms.easy_bar.rds",
        mean_coef_ecdf = "results/plots/plot_pirna_genes_in_lms.mean_coef_ecdf.rds",
        ntotal_ecdf = "results/plots/plot_pirna_genes_in_lms.ntotal_ecdf.rds",
    script:
        "../scripts/coexpression/plot_piRNA_genes_in_lms.R"

rule plot_coef_waterfall:
    input:
        mods = rules.filter_models.output.filtered_tsv
    output:
        rds = "results/plots/plot_coef_waterfall.rds"
    script:
        "../scripts/coexpression/plot_coef_waterfall.R"


rule plot_tx_related_in_lms:
    input:
        mods = rules.filter_models.output.filtered_tsv,
        tfs = "data/Drosophila_melanogaster_TF.txt",
        cofacs = "data/Drosophila_melanogaster_TF_cofactors.txt",
        lkup = rules.make_gene_symbol_lookup.output.tsv,
    output:
        rds = "results/plots/plot_tx_related_in_lms.rds"
    script:
        "../scripts/coexpression/plot_tx-related_in_lms.R"



