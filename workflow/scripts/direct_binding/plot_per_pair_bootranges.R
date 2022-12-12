library(tidyverse)
library(rtracklayer)
library(plyranges)

df <- ifelse(exists("snakemake"),
                          snakemake@input[["rds"]],
                          "results/analysis/direct_binding/per_pair_bootranges.rds") %>%
  read_rds()

g <- df %>% 
  #filter(TF!="INE-1" & p < 0.05) %>%
  mutate(mean.expected = map_dbl(expected,mean)) %>%
  ggplot(aes(observed,mean.expected,color=p<0.05,label=paste(TF,TE,sep="/"))) +
  geom_point() +
  ggrepel::geom_text_repel(data = . %>% filter(p < 0.05),max.iter = 1000)

ggsave(snakemake@output[["png"]],g)
write_rds(g, snakemake@output[["rds"]])