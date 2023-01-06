library(magrittr)
library(tidyverse)

tsvs <- Sys.glob("results/analysis/motifs/xstreme_per_tf/pan/xstreme.tsv")
tsvs <- snakemake@input[["memes"]] %>% paste0("/xstreme.tsv")
tsvs <- tsvs %>% set_names(.,str_extract(.,"(?<=tf\\/).+(?=\\/xstreme.tsv)"))
res <- tsvs %>% 
  map_df(read_tsv, comment="#",.id = "te_group")

shuf_tsvs <- Sys.glob("results/analysis/motifs/xstreme_per_tf_shuffled/*/*/xstreme.tsv")
shuf_tsvs <- snakemake@input[["shuf"]] %>% paste0("/xstreme.tsv")
shuf_tsvs <- shuf_tsvs %>% set_names(.,str_extract(.,"(?<=tf_shuffled\\/).+(?=\\/xstreme.tsv)"))

#https://stats.stackexchange.com/questions/299449/is-empirical-fdr-denominated-by-measurements-or-by-experiments
shuf <- shuf_tsvs %>%
  map_df(read_tsv, comment = "#",.id="te_group") %>%
  #group_by(te_group) %>%
  #slice_min(EVALUE) %>%
  #ungroup() %>%
  separate(te_group, into =c("te_group","rep"), sep = "\\/")

res2 <- bind_rows(result = res, shuffled = shuf, .id = "type") %>%
  group_by(te_group) %>%
  mutate(estimated_fdr = map_dbl(EVALUE, ~{sum(shuf$EVALUE <=.x)/nrow(shuf)})) %>%
  ungroup()

res2 <- res2 %>% filter(type == "result")

write_tsv(res2, snakemake@output[["tsv"]])
