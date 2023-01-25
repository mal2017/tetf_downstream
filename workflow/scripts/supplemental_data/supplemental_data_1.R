library(tidyverse)

x <- read_tsv("upstream/final-models.collected-info.tsv.gz")

x %>%
  filter(significant_model & significant_x) %>%
  dplyr::select(sex=model, feature.x, gene_symbol, feature.y, score = estimate.qnorm, p.value_anova_x, adj_p.value_anova_x) %>%
  write_tsv("~/Downloads/SupplementaryData1.tsv.gz")
