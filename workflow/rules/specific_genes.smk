rule plot_pirna_genes_in_lms:
    input:
        mods = config.get("MERGED_MODELS"),
        pirna = rules.make_pirna_gene_list.output.tsv
    output:
        rds = "results/plots/pirna_genes_in_lms.rds",
        #stats_rds = "results/analysis/specific_genes/pirna_genes_in_lms.rds",
        png = "results/plots/pirna_genes_in_lms.png",
        json = "results/stats/pirna_genes_in_lms.json",
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

rule remap_peaks_near_pirna_genes_contingency:
    input:
        pirna = rules.make_pirna_gene_list.output.tsv,
        remap = rules.annotate_fixed_insertions.output.remap,
        txdb = rules.make_txdb.output.txdb,
        deseq_gr = rules.this_study_kd_deseq2.output.grs,
    output:
        res_rds = "results/analysis/specific_genes/remap_peaks_near_pirna_genes_contingency.rds",
        rds = "results/plots/remap_peaks_near_pirna_genes_contingency.pirna-more-regulated.rds",
        png = "results/plots/remap_peaks_near_pirna_genes_contingency.pirna-more-regulated.png",
        kd_chip_intersect_rds = "results/analysis/specific_genes/remap_peaks_near_pirna_genes_contingency.kd-chip-intersect.rds",
        json = "results/stats/remap_peaks_near_pirna_genes_contingency.json",
    threads:
        4
    script:
        "../scripts/specific_genes/remap_peaks_near_pirna_genes_contingency.R"


# rule remap_peaks_near_pirna_genes_dist:
#     """
#     orthologous approach that shows that the peaks are very nearby piRNA genes
#     """
#     input:
#         pirna = rules.make_pirna_gene_list.output.tsv,
#         remap = rules.annotate_fixed_insertions.output.remap,
#         txdb = rules.make_txdb.output.txdb,
#     output:
#         rds = "results/analysis/specific_genes/remap_peaks_near_pirna_genes.rds",
#         tsv = "results/analysis/specific_genes/remap_peaks_near_pirna_genes.tsv",
#     script:
#         "../scripts/specific_genes/remap_peaks_near_pirna_genes.R"
