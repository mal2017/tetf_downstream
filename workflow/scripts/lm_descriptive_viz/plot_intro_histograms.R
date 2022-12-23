library(tidyverse)

mods <- ifelse(exists("snakmeake"), 
               snakemake@input$mods, 
               "upstream/final-models.collected-info.tsv.gz") %>%
  read_tsv()

g <- mods %>% 
  filter(significant_x) %>%
  ggplot(aes(estimate.qnorm)) +
  geom_histogram() #+facet_wrap(~model) # can facet in the fig script

write_rds(g,snakemake@output[["rds"]])
ggsave(snakemake@output[["png"]],g+facet_wrap(~model))