library(tidyverse)
library(clusterProfiler)
library(BiocGenerics)
library(patchwork)
library(ggpmisc)

FONTSIZE = 7

#gsea <- read_rds("results/analysis/signatures/s2rplus_te_gsea.all.tbl.rds")
gsea <- read_rds(snakemake@input[["all"]])

#gsea.res <- read_rds("results/analysis/signatures/s2rplus_te_gsea.pairs.tbl.rds")
gsea.res <- read_rds(snakemake@input[["pairs"]])


oi <- c("pan","ct","awd","mamo","Unr","NfI","vvl","CG16779")

gsea.plt <- gsea %>%
  filter(comparison %in% oi & comparison %in% filter(gsea.res,pvalue < 0.1)$RNAi) %>%
  mutate(gsea.plt = map2(gsea,comparison,
                         ~gseaplot(.x,geneSetID = .y,title = paste("Perturbation: ",.y,"RNAi in S2R+\n","Signature: TEs coex. w/",.y," in DGRP lines"), 
                                   by="runningScore")))

stat_res <- gsea.res %>%
  filter(RNAi %in% oi & pvalue < 0.1) %>%
  mutate(label = paste0("NES=",round(NES,2),"\np=",format.pval(pvalue,3)))

gs <- gsea.plt %>% dplyr::select(comparison,gsea.plt) %>% 
  deframe() %>%
  map(`$`,data) %>%
  map_df(as_tibble,.id="RNAi") %>%
  arrange(RNAi,x) %>%
  ggplot(aes(x,runningScore)) +
  geom_path(color="tomato",size=rel(2)) +
  facet_wrap(~RNAi, ncol=1,scales="free") +
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_text_npc(data=stat_res,aes(npcx=0.1,npcy=0.2,label=label),vjust="bottom",hjust="left",size=unit(FONTSIZE,"pt")) + 
  ylab("Enrichment Score") +
  xlab("rank") +
  scale_x_continuous(expand = expansion(0))


o <- list(plot = gs, stats=stat_res)

saveRDS(o,snakemake@output[["o"]])

