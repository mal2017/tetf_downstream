library(tidyverse)

#rew_fl <-"results/analysis/rewiring/candidate_rewiring_events.rds" 
rew_fl <- snakemake@input[["rew"]]

rew <- read_rds(rew_fl)

plot_mat <- function(x,y) {
    ggplot(x,aes_string(x=colnames(x)[2],y=colnames(x)[3])) +
      geom_point() +
      facet_wrap(~data_subset,ncol=1) +
      geom_smooth(method = "lm") +
    ggtitle(y)
  }

plot_box <-  function(x, y) {
  x %>% 
  ggplot(aes_string(x="data_subset",y=colnames(x)[3])) +
    geom_boxplot() +
    ggpubr::stat_compare_means()
}

rew <- rew %>%
  filter(padj < 0.1) %>%
  mutate(title = paste(gene_symbol,feature.y,nearby_gene,sep=">")) %>%
  mutate(gg=map2(data,title,plot_mat),
         gg2 = map2(data, title, plot_box))


oi <- rew %>%
  filter(p.value_focus_strains < 0.05 | p.value_other_strains < 0.05) %>%
  arrange(-(estimate1 - estimate2))
  #group_by(sex,gene_symbol,feature.y,nearby_gene) %>%
  #slice_head(n=1) %>%
  #filter(feature.y == "diver2") %>%
  #head(10) %>%

theme_set(theme_classic())
oi %>% 
  pull(gg) %>% .[[2]] +
  xlab("Hr39 expression") +
  ylab("Cyp6a2 expression") +
  facet_wrap(~data_subset)

  filter(oi, nearby_gene %in% c("FBgn0037140","FBgn0260397","FBgn0032615","FBgn0284237","FBgn0262509",
                                "FBgn0052226")) %>%
    pull(gg)

ggplot(filter(rew,padj < 0.1 & sign(rho_focus_strains)!=sign(rho_other_strains)),aes(rho_focus_strains,rho_other_strains,size=-log10(p.value),label=paste(gene_symbol,feature.y,nearby_gene))) + geom_label()

oi %>% filter(gene_symbol=="nau" & feature.y=="INE-1") %>% pull(nearby_gene) %>% walk(message)



oi %>% ggplot(aes(id,-log10(p.value),color=as.factor(sign(estimate1)))) +
  geom_jitter()

oi %>% group_by(id) %>% filter(n() >1) %>%
  filter(feature.y=="INE-1") %>%
  ggplot(aes(id,estimate1-estimate2,label=gene_symbol)) +
  geom_label()


oi %>% pull(nearby_gene) %>%
  walk(message)
