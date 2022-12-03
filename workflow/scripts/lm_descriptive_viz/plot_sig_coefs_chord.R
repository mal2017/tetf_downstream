library(tidyverse)
library(circlize)

#mods_path <- "upstream/final-models.collected-info.tsv.gz"
mods_path <- snakemake@input[["mods"]]

# considering only valid, significant models
mods <- vroom::vroom(mods_path) %>% filter(significant_model)

cutoff <- 0.1

# need a matrix with pairwise counts to make a chord diagram
get_pair_counts <- function(a,b) {
  mods %>%
    filter(if_all(c(a,b),~{.<cutoff})) %>%
    count(model)
}

df <- colnames(mods)[grepl("adj_p.value_anova",colnames(mods))] %>%
  expand_grid(a=.,b=.) %>%
  mutate(n = map2(a,b,get_pair_counts)) %>%
  unnest(n)

# prettier names
dfs <- df %>% mutate(across(c(a,b),~str_remove(.x,"adj_p.value_anova_"))) %>%
  split(.,.$model) %>%
  map(dplyr::select,-model)

# this isn't great beacuse it overrepresents the actual number of significant models
# and doesn't describe the incidence of models with > 2 sigif coefs
#chordDiagram(dfs$female,symmetric = F, link.target.prop = F)


mods %>% count(model,adj_p.value_anova_scaled.copies.y < 0.1)
df

