library(tidyverse)

#mods_path <- "results/analysis/coexpression/filtered_models.tsv.gz"
mods_path <- snakemake@input[["mods"]]

mods <- read_tsv(mods_path)


g_genes <- mods %>%
  count(model,feature.x) %>%
  ggplot(aes(n,fill=model)) +
  geom_histogram(position="dodge")


g_tes <- mods %>%
  count(model,feature.y) %>%
  ggplot(aes(n,fill=model)) +
  geom_histogram(position="dodge")



write_rds(list(TEs=g_tes, genes=g_genes),snakemake@output[["rds"]])
