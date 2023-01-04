library(tidyverse)
library(rtracklayer)
library(plyranges)


lms <- ifelse(exists("snakemake"),snakemake@input[["lms"]],
              "upstream/final-models.collected-info.tsv.gz") %>%
  read_tsv()

x <- import("~/amarel-matt/ptera/subworkflows/modencode_dmel_chipseq/results/peaks/chip_macs/bdgcmp/TCF_WPP_ChIP.seq_ChIP_Rep2_SRR1164491_logFE.bdg",format = "bedgraph")

tes <- import("resources/Tidalbase_transposon_sequence.fasta") %>% names

mean.qp <- x %>% 
  filter(seqnames %in% tes) %>%
  as_tibble() %>%
  group_by(seqnames) %>%
  summarise(score = max(score))

lms %>% 
  filter(gene_symbol == "pan") %>%
  #filter(significant_x) %>%
  dplyr::select(feature.y, significant_x) %>%
  group_by(feature.y) %>%
  summarise(significant_x = any(significant_x)) %>%
  left_join(mean.qp,by=c(feature.y = "seqnames")) %>%
  ggplot(aes(significant_x,score)) +
  geom_boxplot(outlier.shape = NA) +
  #geom_jitter() +
  ggpubr::stat_compare_means(method = "wilcox.test") +
  xlab("")


lms %>% 
  filter(gene_symbol == "pan") %>%
  #filter(significant_x) %>%
  left_join(mean.qp,by=c(feature.y = "seqnames")) %>%
  ggplot(aes(significant_x,score)) +
  geom_boxplot(outlier.shape = NA) +
  #geom_jitter() +
  ggpubr::stat_compare_means(method = "wilcox.test") +
  xlab("TE coexpressed with pan?") +
  ylab("fold enrichment") +
  facet_wrap(~model)


lms %>% 
  filter(gene_symbol == "pan") %>%
  #filter(significant_x) %>%
  left_join(mean.qp,by=c(feature.y = "seqnames")) %>%
  ggplot(aes(score,estimate.qnorm, color=significant_x)) +
  geom_point() +
  facet_wrap(~model)



ggplot(aes(reorder(seqnames,score),score)) +
  geom_point() +
  theme(axis.text.x = element_text(angle=45, hjust=1))

x
