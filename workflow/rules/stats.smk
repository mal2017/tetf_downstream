rule stats_descriptive_lms:
    input:
        mods = config.get("MERGED_MODELS"),
    output:
        json = "results/stats/descriptive_lms.json"
    script:
        "../scripts/stats/descriptive_lms.R"

rule stats_pan_nfi_cg16779_vs_all:
    input:
        mods = config.get("MERGED_MODELS"),
    output:
        json = "results/stats/pan_nfi_cg16779_vs_all.json"
    script:
        "../scripts/stats/pan_nfi_cg16779_vs_all.R"

rule s2rplus_stats:
    input:
        mods = config.get("MERGED_MODELS"),
        limma = rules.s2rplus_limma.output.tsv,
        gsea_pairs = rules.s2rplus_coex_te_gsea_by_de.output.pairs,
        pirna = rules.make_pirna_gene_list.output.tsv,
    output:
        json = "results/stats/s2rplus_stats.json"
    script:
        "../scripts/stats/s2rplus_stats.R"

rule our_kd_stats:
    input:
        gg_pirna_in_kds = rules.plot_pirna_genes_in_our_kd_all.output.rds,
    output:
        json = "results/stats/our_kd_stats.json"
    script:
        "../scripts/stats/our_kd_stats.R"

rule our_kd_signatures_stats:
    input:
        gsea_tbl = rules.ourKD_gsea.output.rds,
    output:
        json = "results/stats/our_kd_signatures_stats.json"
    script:
        "../scripts/stats/our_kd_signatures_stats.R"


rule collect_stats:
    """
    collects stats with the expected structure: model/stat_group/data=statistic/value
    """
    input:
        rules.stats_descriptive_lms.output.json,
        rules.plot_pirna_genes_in_lms.output.json,
        rules.s2rplus_stats.output.json,
        rules.our_kd_stats.output.json,
        rules.our_kd_signatures_stats.output.json,
        rules.remap_peaks_near_pirna_genes_contingency.output.json,
        rules.stats_pan_nfi_cg16779_vs_all.output.json,
    output:
        json = touch("results/stats/collected_stats.json")
    conda:
        "../envs/jq.yaml"
    shell:
        """
        jq -s . {input} > {output}
        """