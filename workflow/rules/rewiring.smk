rule get_candidate_rewiring_events:
    input:
        mods = config.get("MERGED_MODELS"),
        gtf = config.get("GTF"),
        tfs = config.get("TFS"),
        se = config.get("THIS_STUDY_DGRP_RNA"),
        unique_insertions = rules.annotate_insertion_penetrance.output.penetrance,
        all_insertions = rules.annotate_insertion_penetrance.output.all_ins,
    output:
        rds = "results/analysis/rewiring/candidate_rewiring_events.rds"
    params:
        tss_dist = config["REWIRING_ANALYSIS"].get("TSS_DIST"),
        max_tes_per_gene = config["REWIRING_ANALYSIS"].get("MAX_TES_PER_GENE"),
        coef_quantile_cutoff = config["REWIRING_ANALYSIS"].get("COEF_QUANTILE_CUTOFF"),
        min_pct_fixation = config["REWIRING_ANALYSIS"].get("MIN_PCT_FIXATION"),
        max_pct_fixation = config["REWIRING_ANALYSIS"].get("MAX_PCT_FIXATION"),
    script:
        "../scripts/rewiring/regulatory_network_rewiring_01.R"


rule plot_candidate_rewiring_events_scatter:
    input:
        rds = rules.get_candidate_rewiring_events.output.rds
    output:
        rds="results/figs/candidate_rewiring_events_scatter.rds",
        png="results/figs/candidate_rewiring_events_scatter.png"
    script:
        "../scripts/rewiring/plot_candidate_rewiring_events_scatter.R"

#rule plot_candidate_rewiring_events_exemplary: