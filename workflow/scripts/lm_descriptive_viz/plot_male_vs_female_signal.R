library(tidyverse)
library(ggdensity)
library(patchwork)

mods_fl<- "upstream/final-models.collected-info.tsv.gz"
mods_fl <-snakemake@input[["mods"]]

dat <- read_tsv(mods_fl) %>% filter(significant_x)

dat2 <- dat %>% dplyr::select(model,
                      feature.x,feature.y,
                      valid,reproducible, significant_model,
                      estimate,
                      estimate.qnorm,
                      ftest_adjr2,
                      p.value_ftest_r2,
                      sumsq_anova_wolbachia,
                      sumsq_anova_overlap,
                      sumsq_anova_scaled.copies.y,
                      sumsq_anova_overlap.coex.gene,
                      sumsq_anova_x,
                      total_variance) %>%
  mutate(across(contains("p.value"),~{-log10(.x)})) %>%
  mutate(across(contains("sumsq"),~{.x/total_variance})) %>%
  pivot_longer(-c(model,feature.x,feature.y,valid,reproducible,significant_model)) %>%
  pivot_wider(names_from = model, values_from = c(value,valid,reproducible,significant_model))

probs <- c(0.3,0.6,0.9)

gs <- dat2 %>%
  #filter(name %in% c("sumsq_anova_overlap","sumsq_anova_scaled.copies.y")) %>%
  split(.,.$name) %>% 
  imap(.f = ~{
    ggplot(.x,aes(value_male,value_female)) +
      #geom_point() +
      geom_hdr_points() +
      ggtitle(.y)
  })

gs$sumsq_anova_overlap.coex.gene

write_rds(gs,snakemake@output[["rds"]])
ggsave(snakemake@output[["png"]],Reduce(`+`,gs))
