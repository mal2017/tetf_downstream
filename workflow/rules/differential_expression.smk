# rule s2rplus_limma:
#   input:
#     se = "data/tfrnai.se.gene.0.rds",
#     runselect = "data/full_tfrnai_srarunselector.txt",
#     batch = "data/batch_data.tsv.gz",
#   output:
#     tsv = "results/analysis/deg/s2rplus.res.tsv.gz",
#   script:
#     "../scripts/differential_expression/basic_full_limma.R"


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
    rds = "results/figs/this_study_kd_deseq2.rds",
    png = "results/figs/this_study_kd_deseq2.png"
  script:
    "../scripts/differential_expression/plot_ourKD_deseq2.R"