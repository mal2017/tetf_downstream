rule s2rplus_limma:
  input:
    se = config.get("S2RPLUS_RNAI_SE"),
    runselect = config.get("S2RPLUS_RNAI_RUNSELECTOR"),
    batch = config.get("S2RPLUS_RNAI_BATCH"),
    mods = config.get("MERGED_MODELS"),
  output:
    tsv = "results/analysis/deg/s2rplus.res.tsv.gz",
  script:
    "../scripts/differential_expression/basic_full_limma.R"


rule this_study_kd_deseq2:
  input:
    se = config.get("THIS_STUDY_RNAI"),
  output:
    grs = "results/analysis/deg/ourKD.de.grs.rds",
    dds = "results/analysis/deg/ourKD.dds.list.rds"
  script:
    "../scripts/differential_expression/ourKD_deseq2.R"


rule plot_this_study_kd_deseq2:
  input:
    deseq_gr = rules.this_study_kd_deseq2.output.grs,
    mods = config.get("MERGED_MODELS"),
  output:
    rds = "results/plots/this_study_kd_deseq2.rds",
    png = "results/plots/this_study_kd_deseq2.png"
  script:
    "../scripts/differential_expression/plot_ourKD_deseq2.R"

rule plot_kd_info:
  input:
    dds =  rules.this_study_kd_deseq2.output.dds,
  output:
    rds = "results/plots/kd_info.rds",
    #png = "results/plots/kd_info.png"
  script:
    "../scripts/differential_expression/plot_kd_info.R"

rule plot_s2rplus_lfc_waterfall:
  input:
    deg = rules.s2rplus_limma.output.tsv,
    #lkup = rules.make_gene_symbol_lookup.output.tsv,
    #mods = config.get("MERGED_MODELS"),
  output:
    rds = "results/plots/s2rplus_lfc_waterfall.rds",
    png = "results/plots/s2rplus_lfc_waterfall.png"
  script:
    "../scripts/differential_expression/plot_s2rplus_lfc_waterfall.R"