library(tidyverse)

#mods_fl <- "upstream/final-models.collected-info.tsv.gz"
mods_fl <- snakemake@input[["mods"]]

mods <- read_tsv(mods_fl)

dat <- mods %>%
  filter(significant_model) %>%
  mutate(sumsq_anova_remaining = explained_variance*total_variance) %>%
  dplyr::select(model,feature.x,feature.y,total_variance,contains("sumsq"), sumsq_anova_remaining) %>%
  mutate(across(contains("sumsq"),~{100*.x/total_variance})) %>%
  pivot_longer(contains("sumsq"), names_to = "coef", values_to = "var.explained") %>%
  mutate(coef = str_remove(coef,"sumsq_anova_"))

dat <- dat %>% mutate(coef = ifelse(coef == "x", "gene expression", coef))

dat <- dat %>% 
  mutate(coef = fct_reorder(coef,var.explained,.fun = function(.x) {mean(.x,na.rm=T)})) %>%
  mutate(var.explained = var.explained + 1)

g <-  dat %>% ggplot(aes(coef,var.explained,fill=model)) +
  geom_boxplot(outlier.size = 0.01) +
  ylab("% explained variance")

saveRDS(g,snakemake@output[["rds"]])
ggsave(snakemake@output[["png"]],g)
