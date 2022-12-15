library(tidyverse)
library(clusterProfiler)

#limma_path <- "results/analysis/deg/s2rplus.res.tsv.gz"
limma_path <- snakemake@input[["limma"]]

#lkup_path <- "results/resources/gene_symbol_lookup.tsv.gz"
lkup_path <- snakemake@input[["lkup"]]

mods <- ifelse(exists("snakemake"),
               snakemake@input[["mods"]],
               "upstream/final-models.collected-info.tsv.gz") %>%
  read_tsv()

#gsea_pairs_path <- "results/analysis/signatures/s2rplus_te_gsea.pairs.tbl.rds"
gsea_pairs_path <- snakemake@input[["gsea_pairs"]]

# https://www.gsea-msigdb.org/gsea/doc/GSEAUserGuideFrame.html#:~:text=Therefore%2C%20by%20default%2C%20GSEA%20ignores,with%2010%2C000%20to%2020%2C000%20features.
# gsea ignores sizes < 15
gsea_pairs <- read_rds(gsea_pairs_path) #%>% filter(setSize >= 15)

lkup <- read_tsv(lkup_path)

limma <- read_tsv(limma_path)

# -----------------------
n_coex_hits <- mods %>% 
  filter(significant_x) %>%
  dplyr::select(gene_symbol,feature.y) %>% 
  distinct() %>%
  count(gene_symbol)


to_plot <- limma %>%
  dplyr::select(comparison) %>%
  distinct() %>%
  left_join(n_coex_hits, by=c(comparison = "gene_symbol")) %>% 
  mutate(n = replace_na(n,0)) %>%
  left_join(gsea_pairs,by=c(comparison="RNAi")) %>%
  mutate(any.predicted = n > 0) %>%
  mutate(at.least.15 = replace_na(setSize >=15 | !is.na(setSize),F)) %>%
  mutate(has.signature = padj < 0.1) %>%
  mutate(has.signature = replace_na(has.signature,F))

g_bar_hits <- to_plot %>% 
  count(has.signature,at.least.15) %>%
  mutate(out.of = sum(n)) %>%
  mutate(class = case_when(has.signature & at.least.15 ~ "coex. signature detected",
            !has.signature & at.least.15 ~ "signature not detected",
            !at.least.15 & !has.signature ~ "set size < 15 | signature not detected")) %>%
  mutate(class =str_wrap(class,20)) %>%
  mutate(class=fct_reorder(class,-n)) %>%
  ggplot(aes(class,n)) +
  geom_col() +
  xlab("")


#  gsea plt --------------------------------------------------------------------

to_anno <- to_plot %>% count(has.signature,at.least.15) %>%
  filter(at.least.15) %>%
  mutate(lab = if_else(has.signature,"KD enriched for coexpressed TEs:"," no enrichment:")) %>%
  mutate(lab=paste(lab,n))

# nudge_y = -1*sign(NES) * 0.001 * abs(NES)
g_nes_vs_p <- to_plot %>%
  drop_na() %>%
  #mutate(NES=replace_na(NES,0), pvalue= replace_na(pvalue,1)) %>%
  ggplot(aes(-log10(pvalue),NES,color=padj < 0.1)) +
  geom_point() +
  ggrepel::geom_text_repel(data = . %>% filter(padj < 0.1), aes(label=comparison),color="black",fontface="italic",max.overlaps = 20, max.iter = 10000)


#g_nes_vs_p <- g_nes_vs_p + annotate("text",label=paste(to_anno$lab,collapse = "\n"),x=4, y=0, hjust="center",vjust="center",size=1)


write_rds(list(bar = g_bar_hits, ne_vs_p = g_nes_vs_p), snakemake@output[["rds"]])

ggsave(snakemake@output[["png_bar"]], g_bar_hits)
ggsave(snakemake@output[["png_volc"]], g_nes_vs_p)

# tops <- limma %>% 
#   filter(str_detect(feature,"FBgn",negate = T)) %>%
#   filter(adj.P.Val < 0.1) %>%
#   count(comparison,direction = sign(logFC)) %>%
#   arrange(-n) %>%
#   head(10)
# 
# 
# limma %>% 
#   filter(comparison %in% tops$comparison) %>%
#   filter(!str_detect(feature,"FBgn")) %>%
#   ggplot(aes(logFC,-log10(P.Value))) +
#   geom_point(data=.%>% filter(adj.P.Val > 0.1),color="grey") +
#   geom_point(data=.%>% filter(adj.P.Val <= 0.1),color="red") +
#   facet_wrap(~comparison, ncol=5)
# 
# limma %>% 
#   filter(comparison == "pan") %>%
#   filter(!str_detect(feature,"FBgn")) %>%
#   ggplot(aes(logFC,-log10(P.Value))) +
#   geom_point(data=.%>% filter(adj.P.Val > 0.1),color="grey") +
#   geom_point(data=.%>% filter(adj.P.Val <= 0.1),color="red") +
#   facet_wrap(~comparison, ncol=5)
