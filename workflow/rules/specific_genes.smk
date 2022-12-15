rule plot_pirna_genes_in_lms:
    input:
        mods = config.get("MERGED_MODELS"),
        pirna = rules.make_pirna_gene_list.output.tsv
    output:
        rds = "results/plots/pirna_genes_in_lms.rds",
        #stats_rds = "results/analysis/specific_genes/pirna_genes_in_lms.rds",
        png = "results/plots/pirna_genes_in_lms.png",
    script:
        "../scripts/specific_genes/plot_piRNA_genes_in_lms.R"