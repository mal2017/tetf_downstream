library(tidyverse)
library(ggheatmap)
library(uwot)
library(matrixStats)
library(BiocNeighbors)
library(igraph)
library(paletteer)

#lkup_path <- "results/resources/gene_symbol_lookup.tsv.gz"
lkup_path <- snakemake@input[["lkup"]]

#merged_mods_path <- "results/analysis/coexpression/merged_models.tsv.gz"
merged_mods_path <- snakemake@input[["merged_mods"]]

#filtered_mods_path <- "results/analysis/coexpression/filtered_models.tsv.gz"
filtered_mods_path <- snakemake@input[["filtered_mods"]]

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

k_nei <- 1#round(sqrt(nrow(mats$female_model_01)))

knns <- map(mats,findKNN,k = k_nei) %>% # 12 was good by eye which makes sense given sqrt(134) ~ 11.5
  map(`[[`,"index") %>% 
  map_df(as_tibble,rownames="idx0",.id="model") %>%
  pivot_longer(-c(idx0,model),names_to = "kthNN",values_to = "idx") %>%
  mutate(idx=as.character(idx))

knn_lkup <- mats %>%
  map(rownames) %>%
  map_df(enframe,name="idx",value="feature",.id="model") %>%
  mutate(idx=as.character(idx))
  
edges <- left_join(knns,knn_lkup,by=c("idx0"="idx","model"="model")) %>%
  dplyr::select(model,idx0,feature0=feature,idx) %>%
  left_join(knn_lkup, by=c("idx","model")) %>%dplyr::select(model,a=feature0,b=feature) %>%
  filter(a!=b) %>%
  distinct()

Gs <- split(edges,edges$model) %>%
  map(dplyr::select, -model) %>%
  map(as.matrix) %>%
  map(graph_from_edgelist,directed=F) %>%
  map(simplify)

get_comms <- function(func,x=Gs) {
  #message("getting communities",as.character(substitute(func)))
  set.seed(2022)
  map(x,func) %>%
    map(membership) %>%
    map_df(enframe, "feature","community",.id="model") %>%
    mutate(community = as.character(community)) 
}


# try many different comm detection approaches and compare/evaluate later
comms <- list(
              #leiden01 = ~cluster_leiden(.x,resolution_parameter = 0.1,objective_function = "modularity"),
              #leiden025 = ~cluster_leiden(.x,resolution_parameter = 0.25,objective_function = "modularity"),
              #leiden075 = ~cluster_leiden(.x,resolution_parameter = 0.75,objective_function = "modularity"),
              #leiden05 = ~cluster_leiden(.x,resolution_parameter = 0.5,objective_function = "modularity"), # nice
              #leiden10 = ~cluster_leiden(.x,resolution_parameter = 1,objective_function = "modularity"),
              leiden15 = ~cluster_leiden(.x,resolution_parameter = 1.5,objective_function = "modularity"),
              leiden20 = ~cluster_leiden(.x,resolution_parameter = 2,objective_function = "modularity"),
              leiden25 = ~cluster_leiden(.x,resolution_parameter = 2.5,objective_function = "modularity")
              #lead.eigen = ~cluster_leading_eigen(.x), # nice
              #fast.greedy = ~cluster_fast_greedy(.x),
              #label.prop = ~cluster_label_prop(.x),
              #edge.betwe = ~cluster_edge_betweenness(.x),
              #walktrap = ~cluster_infomap(.x),
              #infomap = ~cluster_infomap(.x),
              #spinglass = ~cluster_spinglass(.x)
              ) %>% # nice
  map_df(get_comms,.id="algorithm") %>%
  dplyr::select(model,algorithm,feature,community) 


write_rds(Gs,snakemake@output[["igraph"]])
write_rds(comms,snakemake@output[["comms"]])














