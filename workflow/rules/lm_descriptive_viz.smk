rule plot_n_models:
    input:
        mods = config.get("MERGED_MODELS")
    output:
        png = "results/plots/n_models.png",
        rds = "results/plots/n_models.rds"
    script:
        "../scripts/lm_descriptive_viz/plot_n_models.R"

rule plot_sig_coefs_alluvial:
    input:
        mods = config.get("MERGED_MODELS")
    output:
        png = "results/plots/sig_coefs_alluvial.png",
        rds = "results/plots/sig_coefs_alluvial.rds"
    script:
        "../scripts/lm_descriptive_viz/plot_sig_coefs_alluvial.R"

rule plot_intro_heatmap:
    input:
        mods = config.get("MERGED_MODELS")
    output:
        rds = "results/plots/intro_heatmaps.gg-list.rds",
        pdf = "results/plots/intro_heatmaps.pdf",
    script:
        "../scripts/lm_descriptive_viz/plot_intro_heatmap.R"

rule plot_binary_heatmap:
    input:
        mods = config.get("MERGED_MODELS"),
        classes = config.get("TE_CLASSES")
    output:
        rds = "results/plots/binary_heatmap.rds",
        pdf = "results/plots/binary_heatmap.pdf",
    script:
        "../scripts/lm_descriptive_viz/binary_map.R"

rule plot_intro_ncoex_scatter:
    input:
        mods = config.get("MERGED_MODELS")
    output:
        rds="results/plots/intro_ncoex_scatter.rds",
        png="results/plots/intro_ncoex_scatter.png"
    script:
        "../scripts/lm_descriptive_viz/plot_intro_ncoex_scatter.R"

rule plot_male_vs_female_signal:
    input:
        mods = config.get("MERGED_MODELS")
    output:
        rds="results/plots/male_vs_female_signal.rds",
        png="results/plots/male_vs_female_signal.png"
    script:
        "../scripts/lm_descriptive_viz/plot_male_vs_female_signal.R"

rule plot_variance_explained_overview_boxplot:
    input:
        mods = config.get("MERGED_MODELS")
    output:
        rds="results/plots/variance_explained_overview_boxplot.rds",
        png="results/plots/variance_explained_overview_boxplot.png"
    script:
        "../scripts/lm_descriptive_viz/plot_variance_explained_overview_boxplots.R"

rule plot_intro_histograms:
    input:
        mods = config.get("MERGED_MODELS")
    output:
        rds="results/plots/plot_intro_histograms.rds",
        png = "results/plots/plot_intro_histograms.png",
    script:
        "../scripts/lm_descriptive_viz/plot_intro_histograms.R"


rule plot_sex_specific_barplot:
    input:
        mods = config.get("MERGED_MODELS")
    output:
        rds = "results/plots/sex_specific_barplot.rds",
        png = "results/plots/sex_specific_barplot.png",
    script:
        "../scripts/lm_descriptive_viz/plot_sex_specific_barplot.R"

rule plot_consistency:
    """
    Evaluate consistency of conclusions for terms that are shared by all models for a given TE.
    e.g. copy number and wolbachia
    """
    input:
        mods = config.get("MERGED_MODELS")
    output:
        rds = "results/plots/consistency.rds",
        png = "results/plots/consistency.png",
    script:
        "../scripts/lm_descriptive_viz/plot_consistency.R"

rule plot_exemplary_scatters:
    input:
        mods = config.get("MERGED_MODELS")
    output:
        rds = "results/plots/exemplary_scatters.rds", # no png because this is many plots
    script:
        "../scripts/lm_descriptive_viz/plot_exemplary_scatters.R"