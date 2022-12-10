library(tidyverse)
library(ggrepel)

#tfs_path <- "resources/Drosophila_melanogaster_TF.txt"
#cofacs_path <- "resources/Drosophila_melanogaster_TF_cofactors.txt"

tfs_path <- snakemake@input[["tfs"]]
cofacs_path <- snakemake@input[["cofacs"]]

#kmer_dist_mat_path <- "results/analysis/direct_binding/te_kmer_dist_mat.rds"
kmer_dist_mat_path <- snakemake@input[["kmer_dist"]]
kmer_dist_mat <- read_rds(kmer_dist_mat_path)

#boots_path <- "results/analysis/direct_binding/each_gene_te_kmer_dist.tbl.rds"
boots_path <- snakemake@input[["boots"]]

#null_path <- "results/analysis/direct_binding/null_model.rds"
null_path <- snakemake@input[["null_model"]]

tfs <- read_tsv(tfs_path)
cofacs <- read_tsv(cofacs_path)


null_model <- read_rds(null_path)
boots <- read_rds(boots_path)

null_model_df0 <- null_model %>% enframe(name = "n_TEs") %>% unnest(value) %>% filter(n_TEs!="1") %>% 
  dplyr::rename(dist = "value") %>%
  mutate(n_TEs = as.numeric(n_TEs)) %>%
  mutate(gene_symbol= "null.model") %>%
  mutate(group = "expected")

null_model_df <- null_model %>% enframe(name = "n_TEs") %>% unnest(value) %>% filter(n_TEs!="1") %>% group_by(n_TEs) %>% summarise(expected = median(value)) %>%
  mutate(n_TEs = as.numeric(n_TEs))
# ------------------------------------------------------------------------------
# adding some info

boots <- boots %>%
  mutate(similarity_class = case_when(padj < 0.1 ~ "sig.",
                               T~"n.s.")) %>%
  mutate(type = case_when(feature.x %in% cofacs$Ensembl~"Tx.related",
                          feature.x %in% tfs$Ensembl~"Tx.related",
                          T~"other"))


# ------------------------------------------------------------------------------
# Q2.0 WHich and how many pairs make the cut?: Obs expected

g_dist_diff_by_nte <- boots %>%
  dplyr::select(gene_symbol,n_TEs,dist = in.dist) %>%
  mutate(group =  "observed") %>%
  bind_rows(filter(null_model_df0,n_TEs %in% boots$n_TEs)) %>%
  ggplot(aes(as.factor(n_TEs),dist,fill=group)) +
  geom_boxplot(outlier.shape = NA)

g_dist_diff_all <- boots %>%
  dplyr::select(gene_symbol,n_TEs,dist = in.dist) %>%
  mutate(group =  "observed") %>%
  bind_rows(filter(null_model_df0,n_TEs %in% boots$n_TEs)) %>%
  ggplot(aes(group,dist)) +
  geom_violin(draw_quantiles = c(0.25,0.5,0.75),trim = T) +
  ggpubr::stat_compare_means()


g_tf_groups_not_more_similar <- boots %>%
  filter(similarity_class == "sig.") %>%
  left_join(null_model_df) %>%
  ggplot(aes(type,(in.dist))) +
  #geom_violin(draw_quantiles = c(0.25,0.5,0.75),trim = T) +
  geom_boxplot() +
  ggpubr::stat_compare_means()

g_tfs_have_shifted_ranking <- boots %>%
  filter(similarity_class == "sig.") %>%
  left_join(null_model_df) %>%
  dplyr::select(similarity_class,type,gene_symbol,in.dist,expected) %>%
  group_by(type) %>%
  arrange(log2(in.dist/expected),in.dist) %>%
  mutate(rnk = row_number()) %>%
  filter(rnk <= 100) %>%
  ggplot(aes(rnk,log2(in.dist/expected),color=type)) +
    geom_point()


o <- list(obsexp_by_n_tes = g_dist_diff_by_nte,
          obsexp_box = g_dist_diff_all,
          tfs_have_shifted_ranking = g_tfs_have_shifted_ranking,
     tf_te_groups_not_more_similar = g_tf_groups_not_more_similar)

write_rds(o,snakemake@output[["rds"]])


pdf(snakemake@output[["pdf"]],onefile = T)
o
dev.off()