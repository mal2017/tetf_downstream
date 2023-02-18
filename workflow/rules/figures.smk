rule figure1:
  input:
    cartoon=config.get("LM_CARTOON"),
    heatmaps=rules.plot_intro_heatmap.output.rds,
    ncoex_scatter=rules.plot_intro_ncoex_scatter.output.rds
  output:
    pdf = "results/figures/figure1.pdf"
  params:
    plot_title= "Figure 1"
  script:
    "../scripts/figures/figure1.R"

rule figure1_supp1:
  input:
    var_exp_box = rules.plot_variance_explained_overview_boxplot.output.rds,
    alluvial = rules.plot_sig_coefs_alluvial.output.rds,
    mf_corr = rules.plot_male_vs_female_signal.output.rds,
    sex_specific_bar = rules.plot_sex_specific_barplot.output.rds,
    filtering_barplot = rules.plot_n_models.output.rds,
    coex_te_hist = rules.plot_intro_histograms.output.rds,
    consistency = rules.plot_consistency.output.rds,
  output:
    pdf = "results/figures/figure1_supp1.pdf"
  params:
    plot_title= "Supplementary Figure 1 (related to Figure 1)"
  script:
    "../scripts/figures/figure1_supp1.R"


rule figure2:
  input:
    pirna_box = rules.plot_pirna_genes_in_lms.output.rds,
    gg_gsea_volc = rules.plot_gene_group_gsea_volc.output.rds,
  output:
    pdf = "results/figures/figure2.pdf"
  params:
    plot_title= "Figure 2"
  script:
    "../scripts/figures/figure2.R"

rule figure2_supp1:
  input:
    #pirna_hist = rules.plot_pirna_score_hist.output.rds,
    gg_gsea_table = rules.plot_gene_group_gsea_table.output.png,
  output:
    pdf = "results/figures/figure2_supp1.pdf"
  params:
    plot_title= "Supplementary Figure 2 (related to Figure 2)"
  script:
    "../scripts/figures/figure2_supp1.R"

rule figure3:
  input:
    our_kds = rules.plot_ourKD_gsea.output.rds,
    #s2rnai = rules.plot_tfrnai_gsea.output.rds,
  output:
    pdf = "results/figures/figure3.pdf"
  params:
    plot_title= "Figure 3"
  script:
    "../scripts/figures/figure3.R"

rule figure3_supp1:
  """
  pan scatter
  """
  input:
    exemplary_scatters = rules.plot_exemplary_scatters.output.rds,
  output:
    pdf = "results/figures/figure3_supp1.pdf"
  params:
    plot_title= "Supplementary Figure 3 (related to Figure 3)"
  script:
    "../scripts/figures/figure3_supp1.R"

rule figure3_supptable:
  input:
    kd_info = rules.plot_kd_info.output.rds,
  output:
    pdf = "results/figures/figure3_supptable.pdf"
  params:
    plot_title= "Supplementary Table 1 (related to Figure 3)"
  script:
    "../scripts/figures/figure3_supptable.R"

rule figure3_supp2:
  input:
    kd_gene_coefs_box = rules.plot_kd_gene_coefs_boxplot.output.rds,
    #waterfall = rules.plot_s2rplus_lfc_waterfall.output.rds,
    volcs = rules.plot_this_study_kd_deseq2.output.rds,
  output:
    pdf = "results/figures/figure3_supp2.pdf"
  params:
    plot_title= "Supplementary Figure 4 (related to Figure 3)"
  script:
    "../scripts/figures/figure3_supp2.R"

rule figure4:
  input:
    our_kds = rules.plot_pirna_genes_in_our_kd.output.rds,
    motif_gt = rules.plot_combined_motif_analysis_table.output.png,
    slaidina = rules.plot_pan_coex_te_expression_slaidina.output.rds,
  output:
    pdf = "results/figures/figure4.pdf"
  params:
    plot_title= "Figure 4"
  script:
    "../scripts/figures/figure4.R"

rule figure4_supptable:
  input:
    proximity = rules.remap_peaks_near_pirna_genes_contingency.output.kd_chip_intersect_rds,
  output:
    pdf = "results/figures/figure4_supptable.pdf"
  params:
    plot_title= "Supplementary Table 2 (related to Figure 4)"
  script:
    "../scripts/figures/figure4_supptable.R"

rule figure4_supp1:
  input:
    our_kds_all = rules.plot_pirna_genes_in_our_kd_all.output.rds,
    pirna_slaidina = rules.plot_pirna_gene_expression_slaidina.output.rds,
    motifs = rules.plot_combined_motif_analysis_table.output.rda,
  output:
    pdf = "results/figures/figure4_supp1.pdf"
  params:
    plot_title= "Supplementary Figure 5 (related to Figure 4)"
  script:
    "../scripts/figures/figure4_supp1.R"


rule local_r_info:
  output:
    package_info = "results/figures/local_r_pkgs.tsv",
    platform_info = "results/figures/local_r_platform.tsv",
    external_session_info = "results/figures/local_r_external.tsv"
  script:
    "../scripts/figures/local_r_pkgs.R"

rule figures:
  input:
    rules.figure1.output.pdf,
    rules.figure2.output.pdf,
    rules.figure3.output.pdf,
    rules.figure4.output.pdf,

    rules.figure1_supp1.output.pdf,
    rules.figure2_supp1.output.pdf,
    rules.figure3_supp1.output.pdf,
    rules.figure3_supp2.output.pdf,
    rules.figure4_supp1.output.pdf,    
    
    rules.figure3_supptable.output.pdf,
    rules.figure4_supptable.output.pdf,
  output:
    "results/figures/all_figures.pdf"
  shell:
    """
    gs -dNOPAUSE -sDEVICE=pdfwrite -sOUTPUTFILE={output} -dBATCH {input}
    """