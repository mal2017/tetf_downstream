rule independent_dataset_correlation:
    input:
        mods = config.get("MERGED_MODELS"),
        indep = config.get("INDEPENDENT_DATASET")
    output:
        rds = "results/analysis/independent_dataset/replicate_dataset_correlation.rds",
        json = "results/stats/independent_dataset_correlation.json"
    script:
        "../scripts/independent_dataset/replicate_dataset_correlation.R"

rule plot_independent_dataset_correlation:
    input:
        rds = rules.independent_dataset_correlation.output.rds
    output:
        rds = "results/figures/independent_dataset/independent_dataset_correlation.rds"
    script:
        "../scripts/independent_dataset/plot_independent_dataset_correlation.R"