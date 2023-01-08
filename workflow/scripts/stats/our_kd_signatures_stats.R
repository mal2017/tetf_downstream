library(tidyverse)
library(clusterProfiler)

#gsea_tbl_path <- "results/analysis/signatures/ourKD_gsea.tbl.rds"
gsea_tbl_path <- snakemake@input[["gsea_tbl"]]

gsea_tbl <- read_rds(gsea_tbl_path) 


gsea_tbl %>%
  dplyr::select(comparison, NES, p.adjust) %>%
  mutate(model="all",stat_group = "te_signatures_gsea") %>%
  nest(-model,-stat_group) %>%
  jsonlite::write_json(snakemake@output$json, prettify=T)
