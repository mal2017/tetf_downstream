library(tidyverse)
library(jsonlite)

# ------------------------------------------------------------------------------
# Read in data 
# ------------------------------------------------------------------------------
dat <- ifelse(exists("snakemake"), 
               snakemake@input$gg_pirna_in_kds, 
               "results/plots/pirna_genes_in_our_kd.rds") %>%
  read_rds()


to_write_json <- dat$data %>%
  filter(padj < 0.1) %>%
  count(RNAi) %>%
  pivot_wider(names_from = "RNAi", values_from = "n") %>%
  tibble(model=NA, stat_group="sig_pirna_our_kd", .) %>%
  nest(-model,-stat_group)


write_json(to_write_json, snakemake@output$json, pretty=T)

