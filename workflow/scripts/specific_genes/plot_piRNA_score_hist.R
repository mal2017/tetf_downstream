library(tidyverse)
library(tidytext)
library(patchwork)

mods <- ifelse(exists("snakemake"), snakemake@input[["mods"]], "upstream/final-models.collected-info.tsv.gz") %>%
  read_tsv() %>%
  filter(significant_x)

pirna_tbl <- ifelse(exists("snakemake"), snakemake@input[["pirna"]], "results/resources/pirna_pathway.tsv") %>%
  read_tsv()

g <- mods %>%
  mutate(Czech2013 = gene_symbol %in% filter(pirna_tbl,in.Czech13)$gene_symbol) %>%
  mutate(Handler2013 = gene_symbol %in% filter(pirna_tbl,in.Handler13)$gene_symbol) %>%
  dplyr::select(model, gene_symbol, feature.y, Czech2013, Handler2013, estimate.qnorm) %>%
  pivot_longer(contains("2013"),names_to = "study") %>%
  filter(value) %>%
  ggplot(aes(estimate.qnorm)) +
  geom_histogram() +
  facet_grid(model~study)
  
ggsave(snakemake@output[["png"]],g)
saveRDS(g,snakemake@output[["rds"]])