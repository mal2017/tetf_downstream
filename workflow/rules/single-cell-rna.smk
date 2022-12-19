rule plot_pan_coex_te_expression_slaidina:
    input:
        sce = config.get("SLAIDINA_SCE"),
        gsea = rules.ourKD_gsea.output.rds,
    output:
        png = "results/plots/pan_coex_te_expression_slaidina.png",
        rds = "results/plots/pan_coex_te_expression_slaidina.rds"
    script:
        "../scripts/single-cell-rna/plot_pan_coex_te_expression.R"
