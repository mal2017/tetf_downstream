rule s2rplus_coex_te_gsea:
    """
    All is a tibble with rowns containing nested tbls and/or orther S4 objects.
    pairs is a flat tibble that includes all matches in which the tested signature==KD.
    """
    input:
        deg = rules.s2rplus_limma.output.tsv,
        coex = rules.filter_models.output.filtered_tsv,
    output:
        all = "results/analysis/signatures/s2rplus_te_gsea.all.tbl.rds",
        pairs = "results/analysis/signatures/s2rplus_te_gsea.pairs.tbl.rds",
    threads:
        6
    script:
        "../scripts/signatures/tfrnai_gsea.R"

rule plot_oi_s2rplus_te_gsea:
    input:
        all = rules.s2rplus_coex_te_gsea.output.all,
        pairs = rules.s2rplus_coex_te_gsea.output.pairs,
    output:
        o = "results/plots/plot_oi_s2rplus_te_gsea.rds"
    script:
        "../scripts/signatures/plot-oi-s2rplus-gsea.R"


rule gene_group_gsea:
    input:
        lms = rules.filter_models.output.filtered_tsv,
        tfs = "data/Drosophila_melanogaster_TF.txt",
        cofacs = "data/Drosophila_melanogaster_TF_cofactors.txt",
        pirna = rules.make_pirna_gene_list.output.tsv,
    output:
        rds = "results/analysis/signatures/gene_group_gsea.tbl.rds",
    script:
        "../scripts/signatures/gene_group_gsea.R"

rule plot_gene_group_gsea:
    input:
        gene_group_gsea = rules.gene_group_gsea.output.rds
    output:
        top10 = "results/plots/plot_gene_group_gsea.top10.rds",
        random_walk = "results/plots/plot_gene_group_gsea.random_walk.rds",
    script:
        "../scripts/signatures/plot_gene_group_gsea.R"


rule plot_tfrnai_gsea:
    input:
        limma = rules.s2rplus_limma.output.tsv,
        lkup = rules.make_gene_symbol_lookup.output.tsv,
        filtered_mods = rules.filter_models.output.filtered_tsv,
        gsea_pairs = rules.s2rplus_coex_te_gsea.output.pairs,
    output:
        rds = "results/plots/plot_tfrnai_gsea.plot_list.rds"
    script:
        "../scripts/signatures/plot-tfrnai-gsea.R"

rule ourKD_gsea:
    input:
        filtered_mods = rules.filter_models.output.filtered_tsv,
        res = rules.this_study_kd_deseq2.output.grs
    output:
        rds = "results/analysis/signatures/ourKD_gsea.tbl.rds",
    script:
        "../scripts/signatures/ourKD_gsea.R"


rule plot_ourKD_gsea:
    input:
        lkup = rules.make_gene_symbol_lookup.output.tsv,
        gsea_tbl = rules.ourKD_gsea.output.rds,
    output:
        rds = "results/plots/plot_ourKD_gsea.plot_list.rds"
    script:
        "../scripts/signatures/plot-ourKD-gsea.R"