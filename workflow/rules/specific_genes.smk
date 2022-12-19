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

rule plot_pirna_genes_in_our_kd:
    input:
        res = rules.this_study_kd_deseq2.output.grs,
        pirna = rules.make_pirna_gene_list.output.tsv
    output:
        rds = "results/plots/pirna_genes_in_our_kd.rds",
        png = "results/plots/pirna_genes_in_our_kd.png",
    script:
        "../scripts/specific_genes/plot_piRNA_genes_in_our_kd.R"