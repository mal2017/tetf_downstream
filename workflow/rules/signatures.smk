rule gene_group_gsea:
    input:
        lms = config.get("MERGED_MODELS"),
        tfs = config.get("TFS"),
        cofacs = config.get("COFACS"),
        pirna = rules.make_pirna_gene_list.output.tsv,
    output:
        rds = "results/analysis/signatures/gene_group_gsea.tbl.rds",
    script:
        "../scripts/signatures/gene_group_gsea.R"

rule plot_gene_group_gsea_volc:
    input:
        gene_group_gsea = rules.gene_group_gsea.output.rds
    output:
        png = "results/plots/gene_group_gsea.volc.png",
        rds = "results/plots/gene_group_gsea.volc.rds",
    script:
        "../scripts/signatures/plot_gene_group_gsea_volc.R"

rule plot_gene_group_gsea_table:
    input:
        gene_group_gsea = rules.gene_group_gsea.output.rds
    output:
        png = "results/plots/gene_group_gsea.table.png",
        rda = "results/plots/gene_group_gsea.table.rda",
    script:
        "../scripts/signatures/plot_gene_group_gsea_table.R"

rule ourKD_gsea:
    input:
        mods = config.get("MERGED_MODELS"),
        res = rules.this_study_kd_deseq2.output.grs
    output:
        rds = "results/analysis/signatures/ourKD_gsea.tbl.rds",
    script:
        "../scripts/signatures/ourKD_gsea.R"


rule plot_ourKD_gsea:
    input:
        gsea_tbl = rules.ourKD_gsea.output.rds,
    output:
        rds = "results/plots/ourKD_gsea.rds",
        png = "results/plots/ourKD_gsea.png"
    script:
        "../scripts/signatures/plot-ourKD-gsea.R"

rule per_te_topgo:
    input:
        mods = config.get("MERGED_MODELS"),
    output:
        tsv = "results/analysis/signatures/per_te_topgo.tsv.gz",
    script:
        "../scripts/signatures/topgo_by_te.R"

rule plot_per_te_topgo_heatmaps:
    input:
        tsv = rules.per_te_topgo.output.tsv,
    output:
        rds = "results/plots/per_te_topgo_heatmaps.rds",
        pdf = "results/plots/plot_per_te_topgo_heatmaps.pdf"
    script:
        "../scripts/signatures/plot_per_te_topgo_heatmaps.R"

rule plot_per_te_topgo_overview:
    input:
        tsv = rules.per_te_topgo.output.tsv,
    output:
        rds = "results/plots/per_te_topgo_overview.rds",
        png = "results/plots/plot_per_te_topgo_overview.png"
    script:
        "../scripts/signatures/plot_per_te_topgo_overview.R"



rule s2rplus_coex_te_gsea_by_de:
    """
    All is a tibble with rowns containing nested tbls and/or orther S4 objects.
    pairs is a flat tibble that includes all matches in which the tested signature==KD.
    """
    input:
        deg = rules.s2rplus_limma.output.tsv,
        coex = config.get("MERGED_MODELS")
    output:
        all = "results/analysis/signatures/s2rplus_te_gsea.all.tbl.rds",
        pairs = "results/analysis/signatures/s2rplus_te_gsea.pairs.tbl.rds",
    script:
        "../scripts/signatures/tfrnai_gsea_de.R"

rule plot_tfrnai_gsea:
    """
    Overview of all GSEA results comparing predicted TE signatures to KDs.
    """
    input:
        limma = rules.s2rplus_limma.output.tsv,
        lkup = rules.make_gene_symbol_lookup.output.tsv,
        mods = config.get("MERGED_MODELS"),
        gsea_pairs = rules.s2rplus_coex_te_gsea_by_de.output.pairs,
    output:
        rds = "results/plots/plot_tfrnai_gsea.plot_list.rds",
        png_bar = "results/plots/plot_tfrnai_gsea.bar.png",
        png_volc = "results/plots/plot_tfrnai_gsea.volc.png",
    script:
        "../scripts/signatures/plot-tfrnai-gsea.R"


rule plot_oi_s2rplus_te_gsea:
    """
    Random walk plots of significant GSEA results.
    """
    input:
        all = rules.s2rplus_coex_te_gsea_by_de.output.all,
        pairs = rules.s2rplus_coex_te_gsea_by_de.output.pairs,
    output:
        rds = "results/plots/plot_oi_s2rplus_te_gsea.rds",
        png  = "results/plots/plot_oi_s2rplus_te_gsea.png",
    script:
        "../scripts/signatures/plot-oi-s2rplus-gsea.R"


