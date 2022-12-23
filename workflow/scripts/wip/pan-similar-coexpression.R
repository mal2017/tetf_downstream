library(tidyverse)
library(corrr)

mods <- read_tsv("upstream/final-models.collected-info.tsv.gz")

tfs <-read_tsv("resources/Drosophila_melanogaster_TF.txt")

mods <- mods %>% filter(feature.x %in% tfs$Ensembl)

cdf <- mods %>% 
  filter(model == "female") %>%
  dplyr::select(estimate,feature.y,gene_symbol) %>%
  #dplyr::group_by(estimate.qnorm, gene_symbol) %>%
  #dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  #dplyr::filter(n > 1L)
  pivot_wider(names_from = gene_symbol, values_from = estimate,values_fill = 0) %>%
  corrr::correlate()


cdf %>% corrr::network_plot()


cdf_long <- cdf %>% corrr::stretch() %>%
  filter()


cdf_long %>%
  filter(x == "pan" & y == "sd") %>%
  arrange(-abs(r))
