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

rule plot_pirna_score_hist:
    input:
        mods = config.get("MERGED_MODELS"),
        pirna = rules.make_pirna_gene_list.output.tsv
    output:
        rds = "results/plots/pirna_score_hist.rds",
        png = "results/plots/pirna_score_hist.png",
    script:
        "../scripts/specific_genes/plot_piRNA_score_hist.R"

rule plot_pirna_genes_in_our_kd:
    input:
        res = rules.this_study_kd_deseq2.output.grs,
        pirna = rules.make_pirna_gene_list.output.tsv
    output:
        rds = "results/plots/pirna_genes_in_our_kd.rds",
        png = "results/plots/pirna_genes_in_our_kd.png",
    script:
        "../scripts/specific_genes/plot_piRNA_genes_in_our_kd.R"


rule plot_pirna_genes_in_our_kd_all:
    input:
        res = rules.this_study_kd_deseq2.output.grs,
        pirna = rules.make_pirna_gene_list.output.tsv
    output:
        rds = "results/plots/pirna_genes_in_our_kd_all.rds",
        png = "results/plots/pirna_genes_in_our_kd_all.png",
    script:
        "../scripts/specific_genes/plot_piRNA_genes_in_our_kd_all.R"


rule plot_kd_gene_coefs_boxplot:
    input:
        mods = config.get("MERGED_MODELS"),
    output:
        rds = "results/plots/kd_gene_coefs_boxplot.rds",
        png = "results/plots/kd_gene_coefs_boxplot.png",
    script:
        "../scripts/specific_genes/plot_kd_gene_coefs_boxplot.R"