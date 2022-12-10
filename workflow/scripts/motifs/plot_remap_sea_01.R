library(tidyverse)


df <- read_tsv("results/analysis/motifs/remap_peak_sea.tsv.gz")


df %>% group_by(DB) %>%
  summarise(`>0 de novo motif overrep. in REMAP peaks` = any(QVALUE < 0.1)) %>%
  ggplot(aes(`>0 de novo motif overrep. in REMAP peaks`)) +
  geom_bar()


df %>%
  filter(QVALUE < 0.1) %>%
  arrange(-ENR_RATIO) %>%
  ggplot(aes(log2(ENR_RATIO),-log10(PVALUE),label=DB)) +
  geom_point() +
  ggrepel::geom_text_repel()


df %>%
  filter(QVALUE < 0.1) %>%
  ggplot(aes(`TP%`,`FP%`,size=-log10(PVALUE),label=DB)) +
  geom_point() +
  ggrepel::geom_text_repel()
