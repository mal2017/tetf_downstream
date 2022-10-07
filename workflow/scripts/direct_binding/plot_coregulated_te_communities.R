library(ggraph)
library(tidyverse)
library(graphlayouts)
library(paletteer)

#Gs_path <- "results/analysis/direct_binding/coregulated_te_communities.igraph_list.rds"
Gs_path <- snakemake@input[["igraph"]]


#comms_path <- "results/analysis/direct_binding/coregulated_te_communities.communities_tbl.rds"
comms_path <- snakemake@input[["comms"]]

Gs <- read_rds(Gs_path)
comms <-read_rds(comms_path)


graphs_tbl <- Gs %>% 
  enframe(name = "model",value = "igraph") %>%
  left_join(comms, by="model") %>%
  nest(data=c(-model,-igraph,-algorithm))



get_plottable_graph <- function(gr,df) {
  tidygraph::as_tbl_graph(gr) %>%
    left_join(filter(df, feature %in% names(igraph::V(gr))),by=c(name="feature"))
}

plot_graph <- . %>%
  {x <- .;
  set.seed(2022);
  ggraph(graph=x,layout = 'stress') +
  geom_edge_link(alpha=0.1,check_overlap = T) +
  #geom_node_point(size = 8,aes(color=community)) +
  geom_node_label(aes(label = name,fill=community),colour = 'black', vjust = 0.4) +
  scale_fill_paletteer_d("ggsci::default_igv") +
  #guides(color="none",fill="none") +
  theme_graph()
  }


graphs_tbl <- graphs_tbl %>%
  mutate(tidygr = map2(igraph,data,get_plottable_graph)) %>%
  mutate(gg = map(tidygr,plot_graph))

# break glass in case we want titles
#graphs_tbl <- graphs_tbl %>%
#  unite(title,model,algorithm,sep = "~") %>%
#  dplyr::select(title,gg) %>%
#  mutate(gg = map2(gg,title,.f=~{.x + ggtitle(.y)}))


#graphs_tbl %>%
#  filter(str_detect(title,"female")) %>%
#  pull(gg)


write_rds(graphs_tbl,snakemake@output[["rds"]])




# set.seed(2022)
# umaps <- map(mats, ~umap(.,n_neighbors = 12,spread = 4))
# 
# umaps2 <- umaps %>% 
#   map_df(as_tibble, rownames="feature",.id="model") %>%
#   dplyr::rename(UMAP1=V1,UMAP2=V2)
# 
# umaps2 %>%
#   left_join(comms, by = c("model", "feature")) %>%
#   mutate(c2 = paste(model,community,sep=".")) %>%
#   filter(algorithm == "leiden05") %>%
#   filter(model=="male_model_01") %>%
#   ggplot(aes(UMAP1,UMAP2,fill=c2)) +
#   geom_label(aes(label=feature)) +
#   #geom_point() +
#   facet_wrap(~model+algorithm) +
#   scale_fill_paletteer_d("ggsci::default_igv") +
#   guides(color="none",fill="none")



#write_rds(umaps2,snakemake@output[["umap_coords"]])