rule generate_rankings:
    """
    genes ranked by several metrics from each dataset for leading edge analysis
    """
    input:
        mods = config.get("MERGED_MODELS"),
        reps = config.get("INDEPENDENT_DATASET"),
        replicated = rules.independent_dataset_correlation.output.rds # for easily getting list of pairs that are sig in each dataset
    output:
        rds = "results/analysis/enrichment/rankings.rds"
    script:
        "../scripts/enrichment/generate_rankings.R"

rule fb_gg_gsea:
    input:
        rds =  rules.generate_rankings.output.rds,
        pirna = rules.make_pirna_gene_list.output.tsv,
        zad = rules.get_zad_genes.output.tsv,
    output:
        rds =  "results/analysis/enrichment/fb_gg_gsea.rds",
    script:
        "../scripts/enrichment/fb_gg_gsea.R"

rule plot_fb_gg_gsea_dotplots:
    input:
        rds = rules.fb_gg_gsea.output.rds,
    output:
        rds = "results/plots/fb_gg_gsea_dotplots.rds",
        png = "results/plots/fb_gg_gsea_dotplots.png",
    script:
        "../scripts/enrichment/fb_gg_gsea_dotplots.R"

rule plot_fb_gg_gsea_randomwalk:
    input:
        rds = rules.fb_gg_gsea.output.rds,
    output:
        rds = "results/plots/fb_gg_gsea_randomwalk.rds",
        png = "results/plots/fb_gg_gsea_randomwalk.png",
    script:
        "../scripts/enrichment/fb_gg_gsea_randomwalk.R"