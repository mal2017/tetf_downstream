rule plot_n_models:
    input:
        mods = config.get("MERGED_MODELS")
    output:
        png = "results/plots/n_models.png",
        rds = "results/plots/n_models.rds"
    script:
        "../scripts/lm_descriptive_viz/plot_n_models.R"

rule plot_sig_coefs_alluvial:
    input:
        mods = config.get("MERGED_MODELS")
    output:
        png = "results/plots/sig_coefs_alluvial.png",
        rds = "results/plots/sig_coefs_alluvial.rds"
    script:
        "../scripts/lm_descriptive_viz/plot_sig_coefs_alluvial.R"

rule plot_intro_heatmap:
    input:
        mods = config.get("MERGED_MODELS")
    output:
        rds = "results/plots/intro_heatmaps.gg-list.rds",
        pdf = "results/plots/intro_heatmaps.pdf",
    script:
        "../scripts/lm_descriptive_viz/plot_intro_heatmap.R"

rule plot_binary_heatmap:
    input:
        mods = config.get("MERGED_MODELS"),
        classes = config.get("TE_CLASSES")
    output:
        rds = "results/plots/binary_heatmap.rds",
        pdf = "results/plots/binary_heatmap.pdf",
    script:
        "../scripts/lm_descriptive_viz/binary_map.R"

rule plot_intro_ncoex_scatter:
    input:
        mods = config.get("MERGED_MODELS")
    output:
        rds="results/plots/intro_ncoex_scatter.rds",
        png="results/plots/intro_ncoex_scatter.png"
    script:
        "../scripts/lm_descriptive_viz/plot_intro_ncoex_scatter.R"

rule plot_male_vs_female_signal:
    input:
        mods = config.get("MERGED_MODELS")
    output:
        rds="results/plots/male_vs_female_signal.rds",
        png="results/plots/male_vs_female_signal.png"
    script:
        "../scripts/lm_descriptive_viz/plot_male_vs_female_signal.R"

rule plot_variance_explained_overview_boxplot:
    input:
        mods = config.get("MERGED_MODELS")
    output:
        rds="results/plots/variance_explained_overview_boxplot.rds",
        png="results/plots/variance_explained_overview_boxplot.png"
    script:
        "../scripts/lm_descriptive_viz/plot_variance_explained_overview_boxplots.R"

rule plot_piRNA_genes_in_lms:
    input:
        mods = config.get("MERGED_MODELS"),
        pirna = rules.make_pirna_gene_list.output.tsv,
    output:
        ntotal_ecdf_rds = "results/plots/pirna_genes_in_lms.ntotal_ecdf.rds",
        ntotal_ecdf_png = "results/plots/pirna_genes_in_lms.ntotal_ecdf.png",
    script:
        "../scripts/coexpression/plot_piRNA_genes_in_lms.R"


# rule plot_intro_histograms:
#     input:
#         mods = rules.filter_models.output.filtered_tsv,
#     output:
#         rds="results/plots/plot_intro_histograms.rds"
#     script:
#         "../scripts/coexpression/plot_intro_histograms.R"




# rule plot_coef_waterfall:
#     input:
#         mods = rules.filter_models.output.filtered_tsv
#     output:
#         rds = "results/plots/plot_coef_waterfall.rds"
#     script:
#         "../scripts/coexpression/plot_coef_waterfall.R"


# rule plot_tx_related_in_lms:
#     input:
#         mods = rules.filter_models.output.filtered_tsv,
#         tfs = "data/Drosophila_melanogaster_TF.txt",
#         cofacs = "data/Drosophila_melanogaster_TF_cofactors.txt",
#         lkup = rules.make_gene_symbol_lookup.output.tsv,
#     output:
#         rds = "results/plots/plot_tx_related_in_lms.rds"
#     script:
#         "../scripts/coexpression/plot_tx-related_in_lms.R"



