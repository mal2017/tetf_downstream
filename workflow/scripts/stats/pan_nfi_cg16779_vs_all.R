library(tidyverse)
library(jsonlite)

# Read in data
dat <- ifelse(exists("snakemake"), 
              snakemake@input$mods, 
              "upstream/final-models.collected-info.tsv.gz") %>%
  read_tsv()

sig_pairs <- dat %>% filter(significant_x)

n_sig_pairs <- sig_pairs %>% dplyr::select(feature.x, feature.y) %>% distinct() %>% nrow()

res <- sig_pairs %>%
  group_by(model, gene_symbol) %>%
  mutate(estimate.qnorm = abs(estimate.qnorm)) %>%
  summarise(across("estimate.qnorm",.fns = c(mean=mean, median=median))) %>%
  group_by(model) %>%
  mutate(across(contains("estimate.qnorm"), .f=list(cume_dist=cume_dist), .names = "{.col}_{.fn}"),.keep = "all",) %>%
  ungroup()

res_long <- res %>% pivot_longer(-c(model, gene_symbol), names_to = "metric", values_to = "value")

res_with_cutoff <- res  %>% 
  filter(gene_symbol %in% c("pan","NfI","CG16779")) %>%
  dplyr::select(model, gene_symbol, contains("cume_dist")) %>%
  pivot_longer(-c(model, gene_symbol), names_to = "metric", values_to = "cutoff") %>%
  left_join(res_long, by = c("model","metric"),suffix=c(".cutoff",".target"))

res_n_better <- res_with_cutoff %>% 
  filter(value > cutoff) %>%
  dplyr::select(-model) %>%
  distinct() %>%
  count(gene_symbol.cutoff, metric) %>%
  mutate(stat_group = "n_better_than") %>%
  mutate(model = "combined") %>%
  nest(-model, -stat_group)

res_pctile <- res  %>% 
  filter(gene_symbol %in% c("pan","NfI","CG16779")) %>%
  arrange(gene_symbol) %>%
  mutate(stat_group = "relative_rank_genes_of_interest") %>%
  nest(-model, -stat_group) 


write_json(bind_rows(res_n_better, res_pctile), snakemake@output$json, prettify=T)
