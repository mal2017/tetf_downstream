library(tidyverse)
library(rtracklayer)
library(plyranges)

df <- ifelse(exists("snakemake"),
                          snakemake@input[["rds"]],
                          "results/analysis/direct_binding/per_pair_bootranges.rds") %>%
  read_rds()

g <- df %>% 
  #filter(TE!="INE-1" & p < 0.05) %>%
  mutate(mean.expected = map_dbl(expected,mean)) %>%
  ggplot(aes(observed,mean.expected,color=p<0.05,label=paste(TF,TE,sep="/"))) +
  geom_point() +
  ggrepel::geom_text_repel(data = . %>% filter(p < 0.05),max.iter = 1000) +
  ylab("mean(expected)")
  

ggsave(snakemake@output[["png"]],g)
write_rds(g, snakemake@output[["rds"]])


# we adjust p values here
df2 <- df %>% 
  mutate(padj = p.adjust(p,method="BH")) %>%
  mutate(mean.expected = map_dbl(expected,mean)) %>%
  dplyr::select(TF, TE, observed, mean.expected, p, padj)

df2 %>% 
  filter(p < 0.05) %>% filter(TF %in% c("pan","NfI")) %>%
  gt::gt()


df2 %>% 
  filter(padj < 0.1) %>%
  arrange(-(observed/mean.expected)) %>%
  gt::gt()

df2 %>% filter(padj < 0.1) %>%
  group_by(TE) %>%
  summarise(n=dplyr::n(), TFs=paste(TF,collapse=", ")) %>%
  arrange(-n) %>%
  gt::gt()
