library(tidyverse)
library(ggdensity)
library(patchwork)
library(ggrastr)

#mods_fl<- "upstream/final-models.collected-info.tsv.gz"
mods_fl <-snakemake@input[["mods"]]

dat <- read_tsv(mods_fl) %>%
  group_by(feature.x, feature.y) %>%
  filter(any(significant_x)) %>%
  ungroup() 

dat2 <- dat %>% dplyr::select(model,
                      feature.x,feature.y,
                      valid,reproducible, significant_model,
                      #estimate,
                      estimate.qnorm,
                      #ftest_adjr2,
                      #p.value_ftest_r2,
                      #sumsq_anova_wolbachia,
                      #sumsq_anova_overlap,
                      #sumsq_anova_scaled.copies.y,
                      #sumsq_anova_overlap.coex.gene,
                      #sumsq_anova_x,
                      #total_variance
                      ) %>%
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
      #geom_point(size=0.1) +
      #ggrastr::rasterize(geom_point(size=0.001)) +
      #geom_hex(bins=60) +
      ggrastr::rasterise(ggdensity::geom_hdr_points(size=0.001),dpi=600) +
      #ggdensity::geom_hdr(method="freqpoly",) +
      #geom_smooth(method="lm", se=F) +
      #ggtitle(.y) +
      #coord_cartesian(xlim = c(-0.5,0.5), ylim=c(-0.5,0.5)) +
      ylim(c(-0.5,0.5)) + xlim(c(-0.5,0.5)) +
      ggpubr::stat_cor(method="spearman", size=1, label.x = -0.5, label.y=0.4) +
      theme(aspect.ratio = 1)
  })

#gs$estimate.qnorm

write_rds(gs,snakemake@output[["rds"]])
ggsave(snakemake@output[["png"]],Reduce(`+`,gs))
