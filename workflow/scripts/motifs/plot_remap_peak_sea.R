library(tidyverse)


df <- ifelse(exists("snakemake"), snakemake@input[["tsv"]],
             "results/analysis/motifs/remap_peak_sea.tsv.gz") %>%
  read_tsv()

g_bar <- df %>% 
  mutate(padj = p.adjust(PVALUE,method="bonferroni")) %>%
  group_by(DB) %>%
  summarise(`>0 de novo motif overrep. in REMAP peaks` = any(padj < 0.001)) %>%
  ggplot(aes(`>0 de novo motif overrep. in REMAP peaks`)) +
  geom_bar()

g <- df %>%
  #filter(QVALUE < 0.1) %>%
  arrange(-ENR_RATIO) %>%
  ggplot(aes(log2(ENR_RATIO),-log10(PVALUE),color=QVALUE < 0.1,label=DB)) +
  geom_point() +
  ggrepel::geom_text_repel(data= . %>% filter(QVALUE < 0.1)) +
  scale_color_grey(start = 0.6,end=0.4)

ggsave(snakemake@output[["nhits_png"]],g_bar)
ggsave(snakemake@output[["obs_exp_png"]],g)

write_rds(g_bar, snakemake@output[["nhits_rds"]])
write_rds(g, snakemake@output[["obs_exp_rds"]])