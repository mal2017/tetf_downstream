library(tidyverse)
library(igraph)
library(tidygraph)
library(ggraph)
library(paletteer)

#filtered_mods_path <- "results/analysis/coexpression/filtered_models.tsv.gz"
filtered_mods_path <- snakemake@input[["filtered_mods"]]

filtered_mods <- read_tsv(filtered_mods_path)

#kmer_dist_path <- "results/analysis/direct_binding/te_mash.txt"
kmer_dist_path <- snakemake@input[["kmer_dist"]]
kmer_dists <- read_tsv(kmer_dist_path,col_names = c("seqA","seqB","D","p","shared_hashes"))

kmer_dists <- kmer_dists %>%
  filter(seqA!=seqB) %>%
  filter(seqA %in% filtered_mods$feature.y & seqB %in% filtered_mods$feature.y)

k_nei <- 1#sqrt(nrow(as.matrix(kmer_dists)),)

tg <- kmer_dists %>%
  filter(p<1e-10 & D!=1) %>%
  group_by(seqA) %>%
  slice_min(D,n=k_nei) %>%
  mutate(weight = 1-D) %>%
  tidygraph::tbl_graph(#nodes = distinct(dplyr::select(kmer_dists,name=seqA)),
                       edges = .,directed = F)

tg <- tg %>%
  mutate(louv = group_louvain(weights = weight))


#tg %>%
#  ggraph(layout = "nicely") +
#  geom_edge_link() +
#  geom_node_label(aes(label = name,fill=as.character(louv)),colour = 'white', vjust = 0.4) +
#  scale_fill_paletteer_d("ggsci::default_igv")
  

write_rds(tg,snakemake@output[["tg"]])
