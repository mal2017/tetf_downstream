library(tidyverse)
library(tidytext)
library(patchwork)

#mods_path <- "upstream/final-models.collected-info.tsv.gz"
mods_path <- snakemake@input[["mods"]]
mods <- read_tsv(mods_path) #%>% filter(significant_x)

#pirna_path <- "results/resources/pirna_pathway.tsv"
pirna_path <- snakemake@input[["pirna"]]
pirna_tbl <- read_tsv(pirna_path)

cts <- count(mods, model, gene_symbol,relationship=if_else(significant_x,"sig","ns")) %>%
  pivot_wider(names_from = relationship, values_from = n, values_fill = 0) %>%
  mutate(Czech2013 = gene_symbol %in% filter(pirna_tbl,in.Czech13)$gene_symbol) %>%
  mutate(Handler2013 = gene_symbol %in% filter(pirna_tbl,in.Handler13)$gene_symbol) %>%
  mutate(other =map2_lgl(Czech2013, Handler2013, ~{!xor(.x,.y)}))


est <- mods %>%
  mutate(pct.var.expl = sumsq_anova_x/total_variance) %>%
  dplyr::select(model,gene_symbol,feature.y,estimate.qnorm, pct.var.expl) %>%
  mutate(Czech2013 = gene_symbol %in% filter(pirna_tbl,in.Czech13)$gene_symbol) %>%
  mutate(Handler2013 = gene_symbol %in% filter(pirna_tbl,in.Handler13)$gene_symbol) %>%
  mutate(other =map2_lgl(Czech2013, Handler2013, ~{!xor(.x,.y)}))

# ------------------------------------------------------------------------------
# stats computation
#-------------------------------------------------------------------------------
s_sig <- cts %>%
  mutate(across(contains("2013"), as.factor)) %>%
  split(.,.$model) %>%
  map_df(
    ~{
      bind_rows(Handler2013 = broom::tidy(wilcox.test(sig~Handler2013,data=.,alternative="two.sided")),
                Czech2013 = broom::tidy(wilcox.test(sig~Czech2013,data=.,alternative="two.sided")),.id="publication")
    },
    .id="model")

s_est <- est %>%
  mutate(estimate.qnorm = abs(estimate.qnorm)) %>%
  mutate(across(contains("2013"), as.factor)) %>%
  split(.,.$model) %>%
  map_df(
    ~{
      bind_rows(Handler2013 = broom::tidy(wilcox.test(estimate.qnorm~Handler2013,data=.,alternative="two.sided")),
                Czech2013 = broom::tidy(wilcox.test(estimate.qnorm~Czech2013,data=.,alternative="two.sided")),.id="publication")
    },
    .id="model")

s_var <- est %>%
  mutate(across(contains("2013"), as.factor)) %>%
  split(.,.$model) %>%
  map_df(
    ~{
      bind_rows(Handler2013 = broom::tidy(wilcox.test(pct.var.expl~Handler2013,data=.,alternative="two.sided")),
                Czech2013 = broom::tidy(wilcox.test(pct.var.expl~Czech2013,data=.,alternative="two.sided")),.id="publication")
    },
    .id="model")

stats_to_json <- bind_rows(pirna_var_exp_vs_others = s_var,
          pirna_coex.score_vs_others = s_est,
          pirna_n_sig_vs_others = s_sig, .id="stat_group") %>%
  nest(-model,-stat_group)

# ------------------------------------------------------------------------------
# begin plotting
# ------------------------------------------------------------------------------
toplot_continuous <- est %>%
  pivot_longer(-c(model,gene_symbol,feature.y,estimate.qnorm,pct.var.expl)) %>%
  filter(value) %>%
  mutate(name=fct_relevel(name,c("other","Czech2013","Handler2013"))) 
  
toplot_cts <- cts %>% 
  pivot_longer(-c(model,gene_symbol,sig,ns)) %>%
  filter(value) %>%
  mutate(name=fct_relevel(name,c("other","Czech2013","Handler2013")))

g_est <- toplot_continuous %>%
  ggplot(aes(model, abs(estimate.qnorm), fill=name)) +
  geom_boxplot(outlier.shape = NA) +
  coord_cartesian(ylim = c(0,0.25))+
  xlab("gene set") +
  ylab("abs(coexpression score)")

g_var <- toplot_continuous %>%
  ggplot(aes(model, pct.var.expl, fill=name)) +
  geom_boxplot(outlier.shape = NA) +
  coord_cartesian(ylim=c(0,0.02)) +
  xlab("gene set") +
  ylab("variance explained by gene expression") +
  scale_y_continuous(labels = scales::percent)

g_cts <- toplot_cts %>%
  ggplot(aes(model,sig,fill=name)) +
  geom_boxplot(outlier.shape = NA) +
  coord_cartesian(ylim = c(0,15))+
  xlab("gene set") +
  ylab("N. coexpressed TEs")

g <- g_var + g_est + g_cts + plot_layout(guides="collect")

# # ------------------------------------------------------------------------------
# # do we find anything special abt piRNA pathway genes in terms of number of assoc?
# # the total number of reproducible TE associations is higher for piRNA genes
# s_ntotal.ecdf <- cts %>%
#   mutate(across(contains("2013"), as.factor)) %>%
#   split(.,.$model) %>%
#   map_df(
#     ~{
#       bind_rows(Handler2013 = broom::tidy(ks.test(sig~Handler2013,data=.,alternative="two.sided")),
#                 Czech2013 = broom::tidy(ks.test(sig~Czech2013,data=.,alternative="two.sided")),.id="publication")
#       },
#     .id="model")
# 
# g_ntotal.ecdf <- est %>%
#   pivot_longer(c(Czech2013,Handler2013),names_to = "publication",values_to = "identified") %>%
#   #filter(publication == "Handler2013") %>%
#   mutate(gene.type = if_else(identified,"piRNA","other")) %>%
#   mutate(sex=ifelse(str_detect(model,"female"),"female","male")) %>%
#   ggplot(aes(abs(estimate.qnorm),linetype=gene.type,color=sex)) +
#   stat_ecdf() +
#   facet_wrap(~publication)

#saveRDS(list(coex.score = s_est, var.expl = s_var, nsig = s_sig), snakemake@output[["stats_rds"]])
ggsave(snakemake@output[["png"]],g)
saveRDS(g,snakemake@output[["rds"]])
jsonlite::write_json(stats_to_json, snakemake@output$json)