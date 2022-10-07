library(tidyverse)

#tfs_path <- "data/Drosophila_melanogaster_TF.txt"
#cofacs_path <- "data/Drosophila_melanogaster_TF_cofactors.txt"

tfs_path <- snakemake@input[["tfs"]]
cofacs_path <- snakemake@input[["cofacs"]]

#kmer_dist_mat_path <- "results/analysis/direct_binding/te_kmer_dist_mat.rds"
kmer_dist_mat_path <- snakemake@input[["kmer_dist"]]
kmer_dist_mat <- read_rds(kmer_dist_mat_path)

#boots_path <- "results/analysis/direct_binding/each_gene_te_kmer_dist.tbl.rds"
boots_path <- snakemake@input[["boots"]]

boots <- read_rds(boots_path)

tfs <- read_tsv(tfs_path)
cofacs <- read_tsv(cofacs_path)

#extreme_mods_path <- "results/analysis/coexpression/extreme_models.tsv.gz"
#extreme_mods_path <- snakemake@input[["extreme"]]

#extreme <- read_tsv(extreme_mods_path)

# ------------------------------------------------------------------------------
# adding some info

boots <- boots %>%
  mutate(similarity_class = case_when(padj < 0.05 & in.dist < matched.dist.mean ~ "sig.",
                               T~"n.s.")) %>%
  mutate(type = case_when(feature.x %in% cofacs$Ensembl~"Tx.related",
                          feature.x %in% tfs$Ensembl~"Tx.related",
                          T~"other"))


# ------------------------------------------------------------------------------
# Sanity checking
exemplary_boot <- boots %>%
  filter(in.dist < matched.dist.mean) %>%
  filter(n_TEs > 5) %>%
  filter(padj < 0.1) %>%
  arrange(in.dist - matched.dist.mean) %>%
  .[1,]

# sanity - this does equal the pval I calculated previously
#wilcox.test(exemplary_boot$matched.dists[[1]], mu = exemplary_boot$in.dist)

tidy_dist_mat <- function(tes,mat = kmer_dist_mat) {
  submat <- as.matrix(mat)[tes,tes]
  as_tibble(submat,rownames = "te1") %>%
    pivot_longer(-te1,names_to = "te2", values_to = "D")
} 

exemplary_to_plot <- bind_rows(coexpressed_group =exemplary_boot$TEs %>% enframe(name="rep","TEs"),
     bootstrap=exemplary_boot$matched_boots[[1]] %>% enframe(name = "rep","TEs"),.id="group_type") %>%
  mutate(data = map(TEs,tidy_dist_mat)) %>%
  dplyr::select(-TEs) %>%
  unnest(data)
  
g_exemplary_dist_mat <- exemplary_to_plot %>%
   ggplot(aes(te1,te2,fill=D)) + 
  geom_tile() +
  scale_fill_distiller() +
  theme(axis.text.x = element_text(angle=90,hjust=1)) +
  facet_wrap(~group_type + rep,scales="free", ncol=2)

g_exemplary_dist_jitter <- exemplary_to_plot %>%
  filter(te1!=te2) %>%
  unite(col = "grp",group_type,rep,sep = ".") %>%
  mutate(grp = fct_reorder(grp,D)) %>%
  mutate(grp = fct_relevel(grp,"coexpressed_group.1")) %>%
  ggplot(aes(grp,D)) +
  geom_jitter() +
  theme(axis.text.x = element_text(angle=45,hjust=1))


# ------------------------------------------------------------------------------
# Q1: Are the groups we identified significantly more similar than randomly selected groups?

g_are_te_groups_more_similar_than_expected <- boots %>%
  dplyr::select(model,gene_symbol,in.dist,matched.dist.mean) %>%
  pivot_longer(-c(gene_symbol,model),names_to = "group",values_to = "D") %>%
  ggplot(aes(group,D)) +
  ggpubr::stat_compare_means() +
  geom_boxplot() +
  facet_wrap(~model)

# ------------------------------------------------------------------------------
# Q2.0 WHich and how many pairs make the cut?: Obs expected

g_tfs_obs_exp <- boots %>%
  filter(feature.x %in% c(tfs$Ensembl,cofacs$Ensembl)) %>%
  dplyr::select(model,gene_symbol,in.dist,matched.dist.mean,similarity_class) %>%
  ggplot(aes(matched.dist.mean,in.dist)) +
  geom_point(data = . %>% filter(similarity_class != "sig."),color="gray") +
  geom_point(data = . %>% filter(similarity_class == "sig."),color="red") +
  geom_abline(intercept = 0, slope=1) +
  xlab("expected") +
  ylab("observed")


# ------------------------------------------------------------------------------
# Q2.2: WHich and how many pairs make the cut?: nTE assoc
g_related_pairs_scatter <- boots %>%
  dplyr::select(model,n_TEs,gene_symbol,in.dist,matched.dist.mean,padj,similarity_class) %>%
  arrange(similarity_class) %>%
  ggplot(aes(x=n_TEs,in.dist,color=similarity_class)) +
  #geom_point(aes(y=matched.dist.mean),color="gray") +
  geom_point() +
  facet_wrap(~model) +
  theme_classic() +
  ylab("mean intra-group dist") +
  coord_cartesian(xlim=c(2,NA)) +
  scale_x_continuous(expand = expansion(add = 1))

# ----------------------------------------------------------------------------
# obs: most of the low-dist groups are small
g_lowdist_groups_are_small <- boots %>%
  dplyr::select(model,n_TEs,gene_symbol,in.dist,matched.dist.mean) %>%
  pivot_longer(-c(model,gene_symbol,n_TEs),names_to = "group",values_to = "D") %>%
  mutate(bin = cut_interval(as.integer(n_TEs),length = 10)) %>%
  ggplot(aes(bin,D,fill=group)) +
  geom_violin(draw_quantiles = 0.5,scale = "width", trim = T) +
  facet_wrap(~model,ncol=1) +
  theme_classic() +
  theme(axis.text.x = element_text(angle=45,hjust=1)) +
  scale_fill_grey(start = 0.4, end = 0.8)


# ------------------------------------------------------------------------------
# Q:are TFs more likely to have small in-group dists, suggesting direct binding/coregulation?
# strangely, no


g_counts_in_bin <- boots %>%
  filter(similarity_class == "sig.") %>%
  ungroup() %>%
  mutate(bin = cut_interval(as.integer(n_TEs),length = 3)) %>%
  mutate(similarity_class = fct_relevel(similarity_class,"sig.")) %>%
  mutate(type = fct_relevel(type,c("TF","cofactor","other"))) %>%
  #group_by(model,type,similarity_class,bin) %>%
  #tally() %>%
  ggplot(aes(bin)) +
  facet_wrap(~model) +
  geom_bar(aes(fill=type),position="fill") +
  geom_text(data = . %>% count(model,bin), aes(y=0.5,label=n)) +
  coord_flip() +
  scale_y_continuous(labels=scales::percent) +
  ylab("") + xlab("")

g_tf_groups_not_more_similar <- boots %>%
  filter(similarity_class == "sig.") %>%
  mutate(type = fct_relevel(type,c("TF","cofactor","other"))) %>%
  ggplot(aes(type,(in.dist-matched.dist.mean))) +
  facet_wrap(~model) +
  geom_boxplot()

g_tfs_have_shifted_ranking <- boots %>%
  filter(similarity_class == "sig.") %>%
  dplyr::select(model,similarity_class,type,gene_symbol,in.dist,matched.dist.mean) %>%
  group_by(model,type) %>%
  arrange(log2(in.dist/matched.dist.mean),in.dist) %>%
  mutate(rnk = row_number()) %>%
  filter(rnk <= 100) %>%
  ggplot(aes(rnk,log2(in.dist/matched.dist.mean),color=type)) +
  facet_wrap(~model) +
    geom_point()


o <- list(tf_obs_exp = g_tfs_obs_exp,
          are_te_groups_more_similar_than_expected = g_are_te_groups_more_similar_than_expected,
     sig_more_related_binned_counts = g_counts_in_bin,
     lowdist_groups_are_small = g_lowdist_groups_are_small,
     examplary_dist_mats = g_exemplary_dist_mat, 
     exemplary_dist_jitter = g_exemplary_dist_jitter,
     similar_pairs_overview = g_related_pairs_scatter, 
     tf_te_groups_not_more_similar = g_tf_groups_not_more_similar)

write_rds(o,snakemake@output[["rds"]])

