library(tidyverse)
library(tidytext)

#topgo_fl <- "results/analysis/signatures/per_te_topgo.tsv.gz"
topgo_fl <- snakemake@input[["tsv"]]

df <- read_tsv(topgo_fl) %>% 
  #filter(Annotated < 250) %>%
  filter(pval < 0.05)

theme_set(theme(text = element_text(size=5), axis.text.x = element_text(angle=45,hjust=1)))


df2 <- df %>%
  group_by(GO.ID,dir) %>% add_tally(name = "n_tes_show_enrichment") %>% # add
  group_by(name,,ont,dir) %>% add_tally(name = "n_terms_enriched") %>% # add
  arrange(ont,dir,-n_tes_show_enrichment,-n_terms_enriched) %>%
  mutate(name = reorder_within(name,n_terms_enriched,within = list(ont,dir))) %>%
  mutate(GO.ID = reorder_within(GO.ID,n_tes_show_enrichment,within = list(ont,dir))) %>%
  mutate(Term = reorder_within(Term,n_tes_show_enrichment,within = list(ont,dir))) %>%
  nest(-ont,-dir,-GO.ID,-Term,-n_tes_show_enrichment) %>%
  group_by(ont,dir) %>%
  slice_max(n_tes_show_enrichment,n=10) %>%
  arrange(-n_tes_show_enrichment)


g <- df2 %>% 
  unnest(data) %>%
  #filter(ont == "BP" & dir == "both") %>%
  ggplot(aes(name,Term,size=-log10(pval))) +
  scale_y_reordered() +
  scale_x_reordered() +
  geom_point() +
  facet_wrap(dir~ont,scales = "free") +
  scale_fill_manual(values = c("TRUE"="black","FALSE"="white")) +
  theme_linedraw() +
  theme(axis.text.x = element_text(angle=45,hjust=1), text=element_text(size=5))


write_rds(df2,snakemake@output[["rds"]])

ggsave(snakemake@output[["png"]],g,width = 10, height = 10)