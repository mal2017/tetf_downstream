library(tidyverse)
library(tidyverse)
library(circlize)

#mods_path <- "upstream/final-models.collected-info.tsv.gz"
mods_path <- snakemake@input[["mods"]]

# considering only valid, significant models
mods <- vroom::vroom(mods_path) %>% filter(significant_model)

cutoff <- 0.1


# cool plot in theory but absurdly slow.
mods %>%
  mutate(across(contains("sumsq"),~{.x/total_variance})) %>%
  dplyr::select(model,contains("sumsq")) %>%
  mutate(across(contains("sumsq"),str_remove,"sumsq_anova_")) %>%
  mutate(id = row_number()) %>%
  pivot_longer(-c(model,id)) %>%
  ggplot(aes(name,value,group=id)) +
  geom_line(alpha=0.001) +
  facet_wrap(~model)
