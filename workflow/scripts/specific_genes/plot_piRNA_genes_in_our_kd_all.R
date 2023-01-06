library(tidyverse)
library(GenomicRanges)

pirna <- ifelse(exists("snakemake"), snakemake@input[["pirna"]],
                "results/resources/pirna_pathway.tsv") %>%
  read_tsv()


res <- ifelse(exists("snakemake"), snakemake@input[["res"]],
              "results/analysis/deg/ourKD.de.grs.rds") %>%
  read_rds()

res2 <- res$adjusted %>% 
  map_df(as_tibble, .id="RNAi") %>%
  mutate(RNAi = fct_relevel(RNAi)) %>%
  mutate(RNAi= str_remove(RNAi,"knockdown2_")) %>%
  relocate(RNAi,feature)

select_pirna_oi <- . %>% filter(feature %in% pirna$gene_ID)

g <- res2 %>%
  left_join(pirna, by=c(feature="gene_ID")) %>%
  mutate(oi = feature %in% pirna$gene_ID) %>%
  filter(oi) %>% # possible plot all genes?
  arrange(log2FoldChange) %>%
  ggplot(aes(log2FoldChange,-log10(padj), color=padj<0.1,label=gene_symbol)) +
  geom_point(size=rel(0.2)) +
  ggrepel::geom_text_repel(data = . %>% filter(oi & padj < 0.1), max.iter = 10000, color="black") +
  facet_wrap(~RNAi, scales = "free") +
  theme(strip.text = element_text(size=rel(0.5)))


write_rds(g, snakemake@output[["rds"]])
ggsave(snakemake@output[["png"]], g)