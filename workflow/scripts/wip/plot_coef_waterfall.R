# ------------------------------------------------------------------------------
# plot volcano
# this is intended to make the point about looking at extreme coefs
# rather than just reproducible ones
# IE many genes may be lowly correlated with TEs by chance, but not everything
# is extreme - and the extreme ones are the most interesting for
# other reasons - e.g. enrichment for TFs, etc
# ------------------------------------------------------------------------------

library(tidyverse)
library(ggdensity)

#mods_path <- "results/filtered_models.tsv.gz"
mods_path <- snakemake@input[["mods"]]

#theme_set(theme_classic(base_size = rel(5),base_line_size = rel(1.1)) +
#            theme(strip.text = element_text(size=rel(4)),strip.background.x = element_rect(size=rel(5))))


# not filtering by extreme coefs, will show them
mods <- read_tsv(mods_path)


g <- mods %>%
  group_by(model) %>%
  arrange(mean_estimate.qnorm) %>%
  mutate(rnk = row_number()) %>%
  mutate(extreme = ifelse(is.na(coef.quantile) | coef.quantile < 0.1,"sig.","sig.+extreme")) %>%
  ggplot(aes(rnk,mean_estimate.qnorm,label=feature.x,fill=extreme)) +
  geom_col(width = 1.1) +
  facet_wrap(~model,ncol=2, scales="free_x") +
  scale_fill_manual(values = c("sig."="grey","sig.+extreme"="red")) +
  theme(axis.text.x = element_blank()) +
  guides(fill=F) +
  coord_cartesian(ylim=c(-20,20)) +
  xlab("TE/gene pairs") +
  ylab("coef")

write_rds(g,snakemake@output[["rds"]])