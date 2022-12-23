rule plot_pan_coex_te_expression_slaidina:
    input:
        sce = config.get("SLAIDINA_SCE"),
        gsea = rules.ourKD_gsea.output.rds,
    output:
        png = "results/plots/pan_coex_te_expression_slaidina.png",
        rds = "results/plots/pan_coex_te_expression_slaidina.rds"
    script:
        "../scripts/single-cell-rna/plot_pan_coex_te_expression.R"

rule plot_pirna_gene_expression_slaidina:
    input:
        sce = config.get("SLAIDINA_SCE"),
        pirna = rules.make_pirna_gene_list.output.tsv,
    output:
        png = "results/plots/pirna_gene_expression_slaidina.png",
        rds = "results/plots/pirna_gene_expression_slaidina.rds"
    script:
        "../scripts/single-cell-rna/plot_pirna_gene_expression.R"

# rule plot_cell_labels_slaidina:
#     input:
#         sce = config.get("SLAIDINA_SCE"),
#         gsea = rules.ourKD_gsea.output.rds,
#     output:
#         png = "results/plots/cell_labels_slaidina.png",
#         rds = "results/plots/cell_labels_slaidina.rds"
#     script:
#         "../scripts/single-cell-rna/plot_cell_labels_slaidina.R"