library(tidyverse)

library(corrr)

mods <- read_tsv("upstream/final-models.collected-info.tsv.gz")

#tfs <-read_tsv("resources/Drosophila_melanogaster_TF.txt")

#mods <- mods %>% filter(feature.x %in% tfs$Ensembl)

cdf <- mods %>% 
  #filter(model == "female") %>%
  split(.,.$model) %>% 
  map_df(
    .f = ~{.x %>% 
        dplyr::select(estimate,feature.y,gene_symbol) %>%
        #dplyr::group_by(estimate.qnorm, gene_symbol) %>%
        #dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
        #dplyr::filter(n > 1L)
        pivot_wider(names_from = gene_symbol, values_from = estimate,values_fill = 0) %>%
        corrr::correlate() %>%
        corrr::stretch() %>%
        filter(x=="pan")},
    .id="model")


cdf %>%
  filter(y=="Jarid2")
  group_by(y) %>%
  summarise(r=mean(r)) %>% # if this goes to production do the fisher's z
  arrange(-abs(r)) %>%
  pull(y) %>% head(50) %>% walk(message)
  print(n=50)
  filter(str_detect(y,"pygo"))
