library(tidyverse)
library(ggheatmap)
library(uwot)
library(matrixStats)
library(BiocNeighbors)
library(igraph)
library(paletteer)
library(tidygraph)
library(ggraph)
library(rstatix)

#lkup_path <- "results/resources/gene_symbol_lookup.tsv.gz"
lkup_path <- snakemake@input[["lkup"]]

#filtered_mods_path <- "results/resources/extreme_models.tsv.gz"
filtered_mods_path <- snakemake@input[["filtered_mods"]]

#merged_mods_path <- "results/resources/merged_models.tsv.gz"
merged_mods_path <- snakemake@input[["merged_mods"]]


lkup <- read_tsv(lkup_path)

filtered_mods <- read_tsv(filtered_mods_path)

mods <- read_tsv(merged_mods_path) %>% filter(feature.y %in% filtered_mods$feature.y)

mats <- mods %>%
  split(.,.$model) %>%
  map(~dplyr::select(.,feature.x,feature.y,mean_estimate.qnorm)) %>%
  map(distinct) %>%
  map(pivot_wider,names_from="feature.y", values_from = mean_estimate.qnorm, values_fill=0) %>%
  map(column_to_rownames,"feature.x") %>%
  map(as.matrix)

mats <- mats %>% 
  map(~{.x[!rowAlls(.x,value=0),]})

# the abs is important - we don't actually care if a given gene is positively or
# negatively correlated, we know there are multiple plausible ways a negative regulator could 
# be positively correlated and vice versa - we're just getting at groups of TEs
# that are responsive to the same pathways
mats <- mats %>% 
  map(abs) %>% 
  map(t)

tg <- mats %>%
  map(t) %>%
  map(as_tibble,rownames="teA") %>%
  map(column_to_rownames,"teA") %>%
  map(cor_pmat) %>%
  #map(as.matrix) %>%
  map_df(as_tibble,rownames = "teA",.id="model") %>%
  pivot_longer(c(-model,-teA),names_to = "teB",values_to = "corr") %>%
  filter(teA!=teB) %>%
  #filter(corr > 0.6) %>%
  mutate(d=1-corr) %>%
  group_by(model,teA) %>%
  slice_min(d,n = 1) %>%
  nest(-model) %>%
  mutate(graph_tbl = map(data,.f = ~{ tidygraph::tbl_graph(edges = .,directed = F)}))


tg$graph_tbl[[1]] %>%
  ggraph(layout = "fr") +
  geom_edge_link() +
  geom_node_label(aes(label = name), vjust = 0.4) +
  scale_fill_paletteer_d("ggsci::default_igv")


algs <- list(leiden05 = . %>% cluster_leiden(.,resolution_parameter = 0.5, objective_function = "modularity"),
             leiden1 = . %>% cluster_leiden(.,resolution_parameter = 1, objective_function = "modularity"),
  leiden15 = . %>% cluster_leiden(.,resolution_parameter = 1.5, objective_function = "modularity"),
              leiden20 = . %>% cluster_leiden(.,resolution_parameter = 2,objective_function = "modularity"),
              leiden25 = . %>% cluster_leiden(.,resolution_parameter = 2.5,objective_function = "modularity"))


get_comms <- function(g,fx) {
  #message("getting communities",as.character(substitute(func)))
  set.seed(2022)
  g %>% 
    fx %>%
    membership() %>%
    enframe("feature","community") %>%
    mutate(community = as.character(community))
}


comms <- expand_grid(tg,enframe(algs,"algorithm","func")) %>%
  mutate(community = map2(graph_tbl,func,.f=~get_comms(.x,.y)))

comms <- comms %>%
  dplyr::select(model,algorithm,community) %>%
  unnest(community)


Gs <- tg %>% dplyr::select(-data) %>% deframe() %>% map(as.igraph)

write_rds(Gs,snakemake@output[["igraph"]])
write_rds(comms,snakemake@output[["comms"]])














