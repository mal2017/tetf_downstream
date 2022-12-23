library(tidyverse)

#mods_path <- "upstream/final-models.collected-info.tsv.gz"
mods_path <- snakemake@input[["mods"]]
mods <- read_tsv(mods_path) %>% filter(significant_x)


mods <- mods %>% filter(gene_symbol  %in% c("NfI","pan","vvl","awd","ct","Unr","mamo","CG16779"))

g <- mods %>%
  mutate(gene_symbol = fct_reorder(gene_symbol, estimate.qnorm)) %>%
  ggplot(aes(gene_symbol,estimate.qnorm,color=model, fill=model)) +
  geom_boxplot(fill="white") +
  geom_point(position=position_jitterdodge())


ggsave(snakemake@output[["png"]],g)
saveRDS(g,snakemake@output[["rds"]])