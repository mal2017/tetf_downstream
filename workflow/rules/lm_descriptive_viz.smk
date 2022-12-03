rule plot_n_models:
    input:
        mods = config.get("MERGED_MODELS")
    output:
        png = "results/figs/n_models.png",
        rds = "results/figs/n_models.rds"
    script:
        "../scripts/lm_descriptive_viz/plot_n_models.R"

rule plot_sig_coefs_alluvial:
    input:
        mods = config.get("MERGED_MODELS")
    output:
        png = "results/figs/sig_coefs_alluvial.png",
        rds = "results/figs/sig_coefs_alluvial.rds"
    script:
        "../scripts/lm_descriptive_viz/plot_sig_coefs_alluvial.R"

rule plot_intro_heatmap:
    input:
        mods = config.get("MERGED_MODELS")
    output:
        rds = "results/figs/intro_heatmaps.gg-list.rds",
        pdf = "results/figs/intro_heatmaps.pdf",
    script:
        "../scripts/lm_descriptive_viz/plot_intro_heatmap.R"

rule plot_binary_heatmap:
    input:
        mods = config.get("MERGED_MODELS"),
        classes = config.get("TE_CLASSES")
    output:
        rds = "results/figs/binary_heatmap.rds",
        pdf = "results/figs/binary_heatmap.pdf",
    script:
        "../scripts/lm_descriptive_viz/binary_map.R"

rule plot_intro_ncoex_scatter:
    input:
        mods = config.get("MERGED_MODELS")
    output:
        rds="results/figs/intro_ncoex_scatter.rds",
        png="results/figs/intro_ncoex_scatter.png"
    script:
        "../scripts/lm_descriptive_viz/plot_intro_ncoex_scatter.R"

rule plot_male_vs_female_signal:
    input:
        mods = config.get("MERGED_MODELS")
    output:
        rds="results/figs/male_vs_female_signal.rds",
        png="results/figs/male_vs_female_signal.png"
    script:
        "../scripts/lm_descriptive_viz/plot_male_vs_female_signal.R"

rule plot_variance_explained_overview_boxplot:
    input:
        mods = config.get("MERGED_MODELS")
    output:
        rds="results/figs/variance_explained_overview_boxplot.rds",
        png="results/figs/variance_explained_overview_boxplot.png"
    script:
        "../scripts/lm_descriptive_viz/plot_variance_explained_overview_boxplots.R"