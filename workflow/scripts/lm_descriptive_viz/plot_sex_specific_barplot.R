library(tidyverse)

# Read in data
mods <- ifelse(exists("snakmeake"), 
    snakemake@input$mods, 
    "upstream/final-models.collected-info.tsv.gz") %>%
    read_tsv()

# dplyr pipeline that finds the  elements that are male-specific, female-specific, and shared by both
dat <- mods %>%
    filter(significant_x) %>%
    group_by(gene_symbol,feature.x,feature.y) %>%
    summarise(models = paste(model,collapse=","), .groups = "drop")

g <- dat %>%
  ggplot(aes(models)) +
  geom_bar() +
  scale_y_log10() +
  geom_text(data = . %>% count(models), aes(y=n, label=paste0("n=",n)), color="white",position = position_nudge(y=-0.2))

write_rds(g,snakemake@output[["rds"]])
ggsave(snakemake@output[["png"]],g)

#dat %>% count(models)

#dat %>% filter(models == "male,female") %>%
#  dplyr::select(gene_symbol) %>%
#  write_tsv("~/Downloads/test.tsv", col_names = F)

