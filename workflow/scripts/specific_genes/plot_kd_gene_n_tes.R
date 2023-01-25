library(tidyverse)

#mods_path <- "upstream/final-models.collected-info.tsv.gz"
mods_path <- snakemake@input[["mods"]]
mods <- read_tsv(mods_path) %>% filter(significant_x)

mods <- mods %>% filter(gene_symbol  %in% c("NfI","pan","CG16779"))

g <- mods %>%
  count(model, gene_symbol) %>%
  mutate(gene_symbol = fct_reorder(gene_symbol, n)) %>%
  ggplot(aes(gene_symbol,n,color=model, fill=model)) +
  geom_col(position="dodge") +
  xlab("gene") +
  ylab("coexpressed TEs") 


g + theme_classic() + scale_fill_grey() + 
  scale_color_grey()

ggsave(snakemake@output[["png"]],g)
saveRDS(g,snakemake@output[["rds"]])