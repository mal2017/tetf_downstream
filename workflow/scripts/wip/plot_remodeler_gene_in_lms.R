library(tidyverse)
library(tidytext)

#mods_path <- "results/analysis/coexpression/filtered_models.tsv.gz"
mods_path <- snakemake@input[["mods"]]
mods <- read_tsv(mods_path)

mods %>% 
  filter(coef.quantile > 0.9) %>%
  filter(gene_symbol == "polybromo")

#pirna_path <- "results/resources/pirna_pathway.tsv"
oi_path <- snakemake@input[["oi"]]

oi_tbl <- read_tsv(oi_path)

# no filtering by coef extremeness
filtered_cts <- count(mods, model, gene_symbol,relationship) %>%
  pivot_wider(names_from = relationship, values_from = n, values_fill = 0) %>%
  mutate(nTotal = neg+pos) %>%
  mutate(is.piRNA = gene_symbol %in% pirna_tbl$gene_symbol)

mean_coefs <-  mods %>%
  group_by(model,gene_symbol) %>%
  summarise(mean_estimate.qnorm = mean(mean_estimate.qnorm), .groups = "drop") %>%
  mutate(is.piRNA = gene_symbol %in% pirna_tbl$gene_symbol)

##
# plots begin
#
# Note that these plots/analyses are performed without filtering for significantly
# extreme coefs - pairs are included if
# 1. TLS pkg says the coefs are significant
# 2. reproducible across salmon runs


# ------------------------------------------------------------------------------
# do we find anything special abt piRNA pathway genes in terms of number of assoc?
# the total number of reproducible TE associations is higher for piRNA genes
s_ntotal.ecdf <- filtered_cts %>%
  mutate(is.piRNA = as.factor(is.piRNA)) %>%
  split(.,.$model) %>%
  map_df(~broom::tidy(ks.test(nTotal~is.piRNA,data=.,alternative="two.sided")),.id="model")

g_ntotal.ecdf <- filtered_cts %>%
  ggplot(aes(nTotal,linetype=is.piRNA,color=model)) +
  stat_ecdf()

# ------------------------------------------------------------------------------
# do we find anything special abt piRNA pathway genes in terms of absolute coefs?
# The mean abs reproducible coef is higher for piRNA genes
s_mean.coef.ecdf <- mean_coefs %>% 
  #filtered_cts %>%
  mutate(is.piRNA = as.factor(is.piRNA)) %>%
  mutate(mean_estimate.qnorm.abs = abs(mean_estimate.qnorm)) %>%
  split(.,.$model) %>%
  map_df(~broom::tidy(ks.test(mean_estimate.qnorm.abs~is.piRNA,data=.,alternative="two.sided")),.id="model")

g_mean.coef.ecdf <- mean_coefs %>%
  ggplot(aes(abs(mean_estimate.qnorm),linetype=is.piRNA,color=model)) +
  stat_ecdf()


# ------------------------------------------------------------------------------
# just plot a nice volc as an overview
g_volcano <- full_join(filtered_cts,mean_coefs) %>%
  ggplot(aes((mean_estimate.qnorm),nTotal)) +
  geom_point(data = . %>% filter(!is.piRNA),color="grey") +
  geom_point(data = . %>% filter(is.piRNA)) +
  facet_wrap(~model)


# ------------------------------------------------------------------------------
# easy to follow piRNA figure
g_easy_bar <- filtered_cts %>%
  filter(is.piRNA) %>%
  group_by(model) %>%
  slice_max(nTotal, n=50) %>%
  #filter(gene_symbol %in% c("piwi","aub","armi","egg","Hen1","vret","shu","Nup54","wde","zuc","mael")) %>%
  mutate(gene_symbol=reorder_within(gene_symbol,nTotal,model)) %>%
  ggplot(aes(gene_symbol,nTotal,fill=model)) +
  geom_col(position="dodge") +
  facet_wrap(~model,scales = "free",ncol = 1) +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  scale_x_reordered()
  #geom_text(data = . %>% filter(is.piRNA) %>% count(model), aes(label=paste0(n,"/",nrow(oi)),x=10,y=10))


saveRDS(list(stat = s_mean.coef.ecdf, plot = g_mean.coef.ecdf),snakemake@output[["mean_coef_ecdf"]])
saveRDS(list(stat = s_ntotal.ecdf, plot = g_ntotal.ecdf),snakemake@output[["ntotal_ecdf"]])
saveRDS(g_volcano,snakemake@output[["volcano"]])
saveRDS(g_easy_bar,snakemake@output[["easy_bar"]])
