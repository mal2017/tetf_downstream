# ------------------------------------------------------------------------------
library("treeio")
library(tidytree)
library(caret)
library(dendextend)
library(ggtree)
library(paletteer)
library(seriation)
library(fossil)
library(ape)
library(phytools)
library(TreeTools)
library(tidyverse)
library(phangorn)
library(aricode)
library(dendsort)
library(tidygraph)
library(ggraph)
library(vegan)

comms <- read_rds("results/analysis/direct_binding/coregulated_te_communities.communities_tbl.rds") #%>%
  #filter(algorithm=="leiden05")

tg <- read_rds("results/analysis/direct_binding/sequence_similarity_te_communities.tbl_graph.rds")



labs1 <- comms %>% filter(model == "female_model_01" & algorithm == "leiden25") %>% dplyr::select(3,4) %>% deframe() %>% 
  map_dbl(as.numeric)

labs2 <- as_tibble(tg) %>% deframe()


# only consider nodes calssified in both
labs1 <- labs1[intersect(names(labs1),names(labs2))]
labs2 <- labs2[intersect(names(labs1),names(labs2))]
  
#cutree(kmer_hc, k=max(as.numeric(labs1)))[names(labs1)]

# in case anyone says cant use AMI with hierarchies
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7439807/

# https://scikit-learn.org/stable/modules/generated/sklearn.metrics.adjusted_mutual_info_score.html

aricode::clustComp(labs1,labs2)


tg %>%
  left_join(dplyr::filter(comms,model == "male_model_01" & algorithm == "leiden15"),by=c(name="feature")) %>%
  ggraph(layout = "kk") +
  geom_edge_link(alpha=0.1,check_overlap = T) +
  geom_node_label(aes(label=name,fill=community)) +
  scale_fill_paletteer_d("ggsci::default_igv")






#kmer_dists <- read_rds("results/analysis/direct_binding/te_kmer_dist_mat.rds")
#kmer_dists <- as.matrix(kmer_dists) %>% {.[rownames(.) %in% comms$feature,colnames(.) %in% comms$feature]} %>% as.dist()
# https://www.r-bloggers.com/2021/07/three-ways-to-check-and-fix-ultrametric-phylogenies/
# http://phytools.org/mexico2018/ex/2/Intro-to-phylogenies.html
#kmer_phylo<- phangorn::upgma(kmer_dists) # centroid????
# conversion  works fine if we use upgma, not nj-based methods
#kmer_hc <- as.hclust(kmer_phylo)
#kmer_phylo <- dendsort(kmer_hc) %>% as.phylo()


# F1 via caret? https://www.statology.org/f1-score-in-r/
#https://cran.r-project.org/web/packages/dendextend/vignettes/Cluster_Analysis.html#similaritydifference-between-various-clustering-algorithms
#dendlist <- dendlist(kmer = as.dendrogram(hclust(as.dist(as.matrix(kmer_dists)[to_plot_tes,to_plot_tes]),method = "ward.D")))



get_tree_plot <- . %>% mutate(label=feature) %>%
  tidytree::left_join(kmer_phylo,., by="label") %>% 
  {ggtree(.,layout = "ape") + 
  #geom_tree()  +
  #layout_
  geom_tippoint(aes(color=community),size=5) +
  geom_tiplab2(aes(label=paste0(label," (",community,")")),size=3) +
  theme_tree() +
  scale_color_paletteer_d("ggsci::default_igv") #+ # ggsci::default_igv
  #guides(color="none")
    }


gg_tbl <- comms %>%
  nest(data=c(-model,-algorithm)) %>%
  mutate(gg=map(data,get_tree_plot))


gg_tbl %>%
  filter(model=="female_model_01") %>%
  pull(gg)



filter(comms,model == "female_model_01" & algorithm == "leiden15") %>%
  mutate(label=feature) %>%
  tidytree::left_join(kmer_phylo,., by="label") %>% 
  ggtree(layout = "ape") + 
  #geom_tree()  +
  #layout_
  geom_tippoint(aes(color=community),size=5) +
  #geom_tiplab2(aes(label=paste0(label," (",community,")")),size=3) +
  theme_tree() +
  scale_color_paletteer_d("ggsci::default_igv") +
  guides(color="none")












cor.dendlist(dendlist,method = "common_nodes")
cor.dendlist(dendlist,method = "FM_index", k=10)

Bk_plot(dendlist[[3]], dendlist[[2]],
        rejection_line_asymptotic = T,add_E = T)

Bk_plot(dendlist[[1]], dendlist[[2]],
        rejection_line_asymptotic = T,add_E = T)

bk_res <- Bk_plot(dendlist[[1]], dendlist[[3]],
                  rejection_line_asymptotic = T,add_E = T)



dendlist %>% dendlist(which = c(1,3)) %>% 
  ladderize %>% 
  #set("branches_k_color", k=20) %>% 
  untangle(method = "step2side", k_seq = 3:40) %>%
  tanglegram(faster = TRUE, common_subtrees_color_branches = TRUE)

xx <- filter(comms,model == "female_model_01") %>% 
  mutate(label = feature) %>%
  #tidytree::left_join(as.phylo(te_hcs$female_model_01),., by="label") 
  tidytree::left_join(kmer_hc,., by="label")

xy <- filter(comms,model == "male_model_01") %>% 
  mutate(label = feature) %>%
  #tidytree::left_join(as.phylo(te_hcs$male_model_01),., by="label")
  tidytree::left_join(kmer_hc,., by="label")

trs <- list(male = xy, female = xx)
class(trs) <- 'treedataList'




#https://yulab-smu.top/treedata-book/chapter4.html
ggplot(trs, aes(x, y)) + 
  geom_tree()  +
  layout_circular() +
  geom_tippoint(aes(color=community)) +
  geom_tiplab2(aes(label=paste0(label," (",community,")")),hjust=-0.1) +
  theme_tree() +
  scale_color_paletteer_d("ggsci::default_igv") +
  guides(color="none") +
  facet_wrap(~.id)




