library(tidyverse)
library(rtracklayer)

#mods_path <- "upstream/final-models.collected-info.tsv.gz"
mods_path <-snakemake@input[["mods"]]
mods <- read_tsv(mods_path) %>% 
  filter(significant_x) %>%
  dplyr::select(feature.x=gene_symbol,feature.y) %>%
  mutate(is.coex=T) %>%
  distinct()

#deseq_gr_path <- "results/analysis/deg/ourKD.de.grs.rds"
deseq_gr_path <- snakemake@input[["deseq_gr"]]

res <- read_rds(deseq_gr_path)$adjusted %>%
  map_df(as_tibble,.id="RNAi")  %>%
  mutate(feature.x = str_extract(RNAi,"pan|NFI|vvl|ct|mamo|awd|CG16779|Unr")) %>%
  mutate(feature.x=if_else(feature.x == "NFI","NfI",feature.x))

res <- res %>%
  left_join(mods,by=c(feature.x="feature.x",feature="feature.y")) %>%
  mutate(class = case_when(padj >= 0.1 | is.na(padj) ~ "n.s. feature",
                           padj < 0.1 & is.na(is.coex) & str_detect(feature,"FBgn")~"gene",
                           padj < 0.1 & is.na(is.coex) & !str_detect(feature,"FBgn")~"other TE",
                          is.coex & padj < 0.1 ~ "coexpressed TE")) %>%
  mutate(class = fct_relevel(class,rev(c("coexpressed TE","other TE","gene","n.s. feature")))) %>%
  arrange(RNAi,class)


res <- res %>%
  mutate(tissue = str_extract(RNAi,"head|male_gonad|female_gonad")) %>%
  mutate(tissue = case_when(tissue == "male_gonad"~"testis",
                            tissue == "female_gonad"~"ovary",
                            T~tissue)) %>%
  mutate(driver = str_extract(RNAi,"tj|aTub|Mef2.R|C587")) %>%
  mutate(label=paste0("UAS::",driver,"; ",feature.x,"-RNAi (",tissue,")"))



plot_volc <- function(x,y) {
  x %>%
    #filter(class %in% c("coexpressed TE","n.s. feature")) %>%
  ggplot(aes(log2FoldChange,-log10(pvalue),color=class,size=class)) +
    geom_point() +
    #coord_cartesian(xlim=c(-10,10)) +
    scale_color_manual(values = c(`coexpressed TE`="red",`other TE`="blue",`gene`="darkgray",`n.s. feature`="lightgray")) +
    scale_size_manual(values = c(`coexpressed TE`=rel(3),`other TE`=rel(1),`gene`=rel(1),`n.s. feature`=rel(1)))+
    ggtitle(y)
}

gg <- (res %>% plot_volc("")) + facet_wrap(~label)

gg_tbl <- res %>% 
  nest(data=-c(feature.x,tissue,driver,label)) %>%
  mutate(gg=map2(data,str_wrap(label,width = 20),plot_volc))


gg_list <- gg_tbl %>% dplyr::select(label,gg) %>%
  deframe()

write_rds(gg,snakemake@output[["rds"]])

library(patchwork)
g <- (gg_list %>% Reduce(`+`,.)) + plot_layout(guides = "collect") & theme(legend.position = "bottom")

ggsave(snakemake@output[["png"]],g,width = 8.5, height = 11)

# ----------------------
# library(corrr)
# 
# res.wide <- res %>%
#   dplyr::select(label,score =log2FoldChange,feature) %>%
#   filter(!str_detect(feature,"FBgn")) %>%
#   pivot_wider(names_from = label,values_from = score)
# 
# cr <-  correlate(res.wide, method = "pearson")
# 
# autoplot(cr)
# rplot(cr)
# 
# res.wide %>%
#   ggplot(aes(`UAS::aTub; CG16779-RNAi (testis)`,`UAS::tj; CG16779-RNAi (ovary)`)) +
#   geom_point() +
#   theme(aspect.ratio = 1) +
#   ggpubr::stat_cor(method="pearson")
# 
# res.wide %>%
#   ggplot(aes(`UAS::tj; pan-RNAi (ovary)`,`UAS::aTub; pan-RNAi (testis)`)) +
#   geom_point() +
#   theme(aspect.ratio = 1) +
#   ggpubr::stat_cor(method="pearson")
# 
# res.wide %>%
#   ggplot(aes(`UAS::aTub; pan-RNAi (testis)`,`UAS::Mef2.R; Unr-RNAi (head)`)) +
#   geom_point() +
#   theme(aspect.ratio = 1) +
#   ggpubr::stat_cor(method="pearson")
# 
# res.wide %>%
#   ggplot(aes(`UAS::Mef2.R; CG16779-RNAi (head)`,`UAS::tj; CG16779-RNAi (ovary)`)) +
#   geom_point() +
#   theme(aspect.ratio = 1) +
#   ggpubr::stat_cor(method="pearson")
# 
#              