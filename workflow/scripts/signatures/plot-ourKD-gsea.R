library(tidyverse)
library(clusterProfiler)
library(ggpmisc)
library(paletteer)

FONTSIZE=7

#gsea_tbl_path <- "results/analysis/signatures/ourKD_gsea.tbl.rds"
gsea_tbl_path <- snakemake@input[["gsea_tbl"]]

gsea_tbl <- read_rds(gsea_tbl_path) 

gsea_tbl2 <- gsea_tbl %>% 
  filter(kd == ID)

gsea_gg_tbl <- gsea_tbl %>%
  mutate(tissue = str_extract(comparison,"head|male_gonad|female_gonad")) %>%
  mutate(tissue = case_when(tissue == "male_gonad"~"testis",
                            tissue == "female_gonad"~"ovary",
                            T~tissue)) %>%
  mutate(driver = str_extract(comparison,"tj|aTub|Mef2.R")) %>%
  filter(p.adjust < 0.1) %>%
  dplyr::select(comparison,tissue,driver,kd,p.adjust,NES,gsea) %>%
  mutate(gsea.plt = pmap(list(gsea,kd,driver),
                         function(gsea,kd,driver) gseaplot(gsea,geneSetID = kd,title = paste0("UAS::",driver," ",kd,"-RNAi;\nSignature: TEs coex. w/ ",kd," in DGRP lines"), 
                                   by="runningScore")))

stat_res <- gsea_tbl2 %>%
  mutate(label = paste0("NES=",round(NES,2),"\npadj=",format.pval(p.adjust,3))) %>%
  dplyr::rename(RNAi="comparison")


to_plot <- gsea_gg_tbl %>% 
  dplyr::select(comparison,gsea.plt) %>% 
  deframe() %>%
  map(`$`,data) %>%
  map_df(as_tibble,.id="RNAi") %>%
  arrange(RNAi,x) %>%
  left_join(dplyr::select(stat_res,-gsea,-Description,-data)) %>%
  left_join(dplyr::select(gsea_gg_tbl,RNAi=comparison,tissue,driver)) %>%
  mutate(class = case_when(p.adjust < 0.1 & NES > 0~"pos",
                           p.adjust < 0.1 & NES < 0~"neg",
                           T~"n.s."))

to_plot <- to_plot %>%
  group_by(RNAi) %>%
  mutate(best.runningScore = runningScore[which.max(abs(runningScore))],
         nes.scale.factor = NES/best.runningScore,
         runningScore.nes = runningScore*nes.scale.factor)

# tissue <- c(vvl="head",
#      NfI="head",
#      CG16779.ovary="ovary",
#      CG16779.head="head",
#      Unr="head",
#      ct="ovary",
#      mamo="ovary") %>%
#   enframe(name = "RNAi",value="tissue")


gs <- to_plot %>%
  #left_join(tissue) %>%
  mutate(label=paste0("UAS::",driver," ",kd,"-RNAi (",tissue,")")) %>%
  mutate(significance = ifelse(p.adjust < 0.1,"sig.","n.s.")) %>%
  ggplot(aes(x,runningScore.nes,color=label,linetype=significance)) +
  geom_path(size=rel(1.1)) +
  #scale_color_manual(values=c(pos="red",neg="blue",n.s.="gray")) +
  #facet_wrap(~RNAi, ncol=1,scales="free") +
  geom_hline(yintercept = 0, linetype="dashed") +
  #ggrepel::geom_text_repel(max.overlaps = 10, max.iter = 100,color="black",fontface="italic",force_pull = 0.1,
  #                         data=. %>% 
  #                           filter(p.adjust < 0.1) %>% group_by(RNAi) %>% slice_max(abs(runningScore),n = 1, with_ties = F),
  #                         aes(label=str_wrap(RNAi,10)),vjust="bottom",hjust="left",size=unit(FONTSIZE-2,"pt")) + 
  ylab("NES") +
  xlab("rank") +
  #scale_x_continuous(expand = expansion(0)) +
  #scale_color_paletteer_d("ggsci::default_ucscgb") +
  scale_linetype_manual(values = c("n.s."="dotted","sig."="solid"))


o <- list(plot = gs, stats=stat_res)

saveRDS(o,snakemake@output[["rds"]])
ggsave(snakemake@output[["png"]],gs)




# 
# # -----------------------
# n_coex_hits <- filtered_mods %>% 
#   dplyr::select(gene_symbol,feature.y) %>% 
#   distinct() %>%
#   count(gene_symbol)
# 
# 
# to_plot <- limma %>%
#   dplyr::select(comparison) %>%
#   distinct() %>%
#   left_join(n_coex_hits, by=c(comparison = "gene_symbol")) %>% 
#   mutate(n = replace_na(n,0)) %>%
#   left_join(gsea_pairs,by=c(comparison="RNAi")) %>%
#   mutate(any.predicted = n > 0) %>%
#   mutate(at.least.15 = replace_na(setSize >=15 | !is.na(setSize),F)) %>%
#   mutate(has.signature = p.adjust < 0.1) %>%
#   mutate(has.signature = replace_na(has.signature,F))
# 
# g_bar_hits <- to_plot %>% 
#   count(has.signature,at.least.15) %>%
#   mutate(out.of = sum(n)) %>%
#   mutate(class = case_when(has.signature & at.least.15 ~ "coex. signature detected",
#             !has.signature & at.least.15 ~ "signature not detected",
#             !at.least.15 & !has.signature ~ "set size < 15 | signature not detected")) %>%
#   mutate(class =str_wrap(class,20)) %>%
#   mutate(class=fct_reorder(class,-n)) %>%
#   ggplot(aes(class,n)) +
#   geom_col() +
#   xlab("")
# 
# 
# #  gsea plt --------------------------------------------------------------------
# 
# # nudge_y = -1*sign(NES) * 0.001 * abs(NES)
# g_nes_vs_p <- to_plot %>%
#   #mutate(NES=replace_na(NES,0), pvalue= replace_na(pvalue,1)) %>%
#   ggplot(aes(-log10(pvalue),NES,color=p.adjust < 0.1)) +
#   geom_point() +
#   ggrepel::geom_text_repel(data = . %>% filter(p.adjust < 0.1), aes(label=comparison),color="black",fontface="italic",max.overlaps = 20, max.iter = 1000)
# 
# 
# write_rds(list(bar = g_bar_hits, ne_vs_p = g_nes_vs_p), snakemake@output[["rds"]])
# 
# # tops <- limma %>% 
# #   filter(str_detect(feature,"FBgn",negate = T)) %>%
# #   filter(adj.P.Val < 0.1) %>%
# #   count(comparison,direction = sign(logFC)) %>%
# #   arrange(-n) %>%
# #   head(10)
# # 
# # 
# # limma %>% 
# #   filter(comparison %in% tops$comparison) %>%
# #   filter(!str_detect(feature,"FBgn")) %>%
# #   ggplot(aes(logFC,-log10(P.Value))) +
# #   geom_point(data=.%>% filter(adj.P.Val > 0.1),color="grey") +
# #   geom_point(data=.%>% filter(adj.P.Val <= 0.1),color="red") +
# #   facet_wrap(~comparison, ncol=5)
# # 
# # limma %>% 
# #   filter(comparison == "pan") %>%
# #   filter(!str_detect(feature,"FBgn")) %>%
# #   ggplot(aes(logFC,-log10(P.Value))) +
# #   geom_point(data=.%>% filter(adj.P.Val > 0.1),color="grey") +
# #   geom_point(data=.%>% filter(adj.P.Val <= 0.1),color="red") +
# #   facet_wrap(~comparison, ncol=5)
