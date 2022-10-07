# barchart for piRNA genes -----------------------------------------------------
g_lm_bar <- read_rds("results/plots/plot_pirna_genes_in_lms.easy_bar.rds") 

g_lm_bar <- g_lm_bar + scale_fill_manual(values=c(male="red",female="black"))

ggsave("~/Downloads/CSHL22_figs/pirna_simple_bar.svg",g_lm_bar, height = 7.6, width = 8.3)


#pirna_volc <- read_rds("results/plots/plot_pirna_genes_in_lms.volc.rds")
#read_rds("results/plots/plot_pirna_genes_in_lms.volc.rds")




#  stats about TE obs/exp per gene --------------------------------------------
g_intra_grp_dist_n <- intra_grp_dist_all$similar_pairs_overview$data %>% 
  arrange(in.dist - matched.dist.mean) %>%
  filter(similarity_class == "sig.") %>%
  inner_join(tx.related) %>%
  group_by(gene_symbol,regulator.type,family) %>%
  summarise(nTEs = paste(n_TEs,collapse = "/"),.groups = "drop" ) %>%
  mutate(family=ifelse(family == "Others","unchar.",family)) %>%
  mutate(class = paste(regulator.type,family,sep="/")) %>%
  count(class,regulator.type,family, sort = T) %>%
  slice_max(n,n=10) %>%
  mutate(class =fct_reorder(class,n)) %>%
  ggplot(aes(n,class)) +
  geom_col() +
  ylab("")

ggsave("~/Downloads/CSHL22_figs/sig_intra_grp_dist_n_tf_class.svg",
       g_intra_grp_dist_n, width = 7, height = 6)




# transcription regulators ----------------------------------
ggsave("~/Downloads/CSHL22_figs/gene_group_gsea_bar.svg",
       read_rds("results/plots/plot_gene_group_gsea.top10.rds")$plot + theme(axis.text.x = element_text(vjust=1)),width = 8, height = 8.4)



# tfrnai -----------------------------------------------------------------------
ggsave("~/Downloads/CSHL22_figs/tfrnai_coex_signature_gsea.bar.svg",read_rds("results/plots/plot_tfrnai_gsea.plot_list.rds")$bar, width = 4, height = 4.2)
