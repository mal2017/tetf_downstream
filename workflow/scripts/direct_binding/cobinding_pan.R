library(tidyverse)
library(plyranges)
library(rtracklayer)

remap <- read_rds("results/resources/remap.gr.rds")

filter_by_overlaps(remap$pan, remap$NfI)

remap$pan
remap$NfI

remap2 <- remap %>% unlist() %>% mutate(ChIP = names)

pan.anno <- remap$pan %>% join_overlap_left(remap2) %>% 
  as_tibble()

pan.anno <- pan.anno %>% unite("locus",-ChIP, sep = "_") %>% distinct()

n_sites <-  length(unique(pan.anno$locus))

pan.anno <- pan.anno %>% 
  count(ChIP) %>% 
  mutate(pct.cobound = n/n_sites)

pan.anno %>% arrange(pct.cobound) %>% print(n=50)


lms <- read_tsv("upstream/final-models.collected-info.tsv.gz")

lms %>% 
  filter(gene_symbol %in% c("sv","pan")) %>% 
  dplyr::select(model, gene_symbol,feature.y, estimate.qnorm) %>%
  pivot_wider(names_from = gene_symbol, values_from = estimate.qnorm) %>%
  ggplot(aes(sv, pan)) +
  geom_point() +
  facet_wrap(~model) +
  ggpubr::stat_cor()

# now for each ChIP, find the percentage of peaks also bound by pan
others.anno <- remap2 %>% join_overlap_left(mutate(remap$pan,ChIP2="pan")) %>%
  as_tibble()

others.anno <- others.anno %>% unite("locus",-c(ChIP,ChIP2), sep = "_") %>% distinct()

others.anno <- others.anno %>% group_by(ChIP) %>% summarize(pct.pan.bound = sum(ChIP2=="pan",na.rm = T)/dplyr::n())

others.anno %>% arrange(pct.pan.bound) %>% print(n=50)
