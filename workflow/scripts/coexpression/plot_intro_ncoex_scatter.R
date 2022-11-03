library(tidyverse)
library(ggdensity)

#mods_path <- "results/resources/extreme_models.tsv.gz"
mods_path <- snakemake@input[["mods"]]

mods <- read_tsv(mods_path)

probs <- c(0.3,0.6,0.9)

genes_df <- mods %>%
  count(model,feature=feature.x,relationship) %>%
  pivot_wider(names_from = relationship, values_from = n,values_fill = 0)

te_df <- mods %>%
  count(model,feature=feature.y,relationship) %>%
  pivot_wider(names_from = relationship, values_from = n,values_fill = 0)


g_tes <- ggplot(te_df,aes(pos,neg)) +
  geom_hdr_points(probs=probs,position = "jitter") +
  xlab("positively coexpressed genes") +
  ylab("negatively coexpressed genes") +
  ggtitle("TEs")

g_genes <- ggplot(genes_df,aes(pos,neg)) +
  geom_hdr_points(probs=probs,position = "jitter") +
  xlab("positively coexpressed TEs") +
  ylab("negatively coexpressed TEs") +
  ggtitle("genes")
  





write_rds(list(tes = g_tes,genes=g_genes),snakemake@output[["rds"]])