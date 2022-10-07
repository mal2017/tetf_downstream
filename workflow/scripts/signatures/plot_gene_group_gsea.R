library(tidyverse)
library(clusterProfiler)
library(patchwork)
library(ggpmisc)

FONTSIZE = 7

#gsea_path <- "results/analysis/signatures/gene_group_gsea.tbl.rds"
gsea_path <- snakemake@input[["gene_group_gsea"]]

gsea_tbl <- read_rds(gsea_path)


s_gene_group_gsea <- gsea_tbl%>%
  dplyr::select(pathway,pval,padj,NES,ES) %>%
  mutate(label = paste0("NES=",round(NES,2),"\np=",format.pval(pval,3)))

g_top10_bar <- s_gene_group_gsea %>%
  mutate(pathway = ifelse(pathway == "Tx.related","AnimalTFDB 3.0",pathway)) %>%
  arrange(pval) %>%
  mutate(rnk = row_number()) %>%
  filter(str_detect(pathway,"piRNA") | rnk <= 5) %>%
  mutate(pathway2 = str_wrap(paste0(rnk,". ",pathway),width = 20)) %>%
  mutate(pathway2 = fct_reorder(pathway2,-rnk)) %>%
  ggplot(aes(-log10(pval),pathway2)) +
  geom_col() +
  geom_text(aes(label=pathway2),x=0.1,hjust=0,color="white",fontface="bold") +
  theme(axis.text.y = element_blank(),axis.ticks.y=element_blank()) +
  ylab("gene group")

oi <- c("Tx.related","piRNA",
        "C2H2 ZINC FINGER TRANSCRIPTION FACTORS",
        "OTHER DNA BINDING DOMAIN TRANSCRIPTION FACTORS")

g_tx_related_walk <- gsea_tbl %>%
  filter(pathway %in% oi) %>%
  dplyr::select(pathway,gsea.plt_tbl) %>% 
  unnest(gsea.plt_tbl) %>%
  arrange(pathway,x) %>%
  ggplot(aes(x,y)) +
  geom_path(color="tomato",size=rel(2)) +
  facet_wrap(~pathway, ncol=2,scales="free") +
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_text_npc(data=filter(s_gene_group_gsea,pathway %in% oi),aes(npcx=0.1,npcy=0.2,label=label),vjust="bottom",hjust="left",size=unit(FONTSIZE,"pt")) + 
  ylab("Enrichment Score") +
  xlab("rank") +
  scale_x_continuous(expand = expansion(0))

obar <- list(plot = g_top10_bar, stat = s_gene_group_gsea)
owalk <- list(plot = g_tx_related_walk, stat = s_gene_group_gsea)

write_rds(obar,snakemake@output[["top10"]])
write_rds(owalk,snakemake@output[["random_walk"]])