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



