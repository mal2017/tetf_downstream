library(tidyverse)

#mods_path <- "upstream/final-models.collected-info.tsv.gz"
mods_path <- snakemake@input[["mods"]]

mods <- vroom::vroom(mods_path)

dat <- list(total = mods %>% count(model),
     valid = mods %>% filter(valid) %>% count(model),
     reproducible = mods %>% filter(reproducible) %>% count(model),
     `signif. model` = mods %>% filter(significant_model) %>% count(model),
     `signif. gene coef` = mods %>% filter(reproducible & valid & significant_model & significant_x) %>% count(model)) %>%
  bind_rows(.id="filter")

g <- dat %>%
  mutate(filter = fct_reorder(filter,-n)) %>%
  ggplot(aes(filter,n,fill=model)) +
  geom_col(position = "dodge")

write_rds(g,snakemake@output[["rds"]])
ggsave(snakemake@output[["png"]],g)
