library(tidyverse)
library(rtracklayer)

#mods_path <- "results/analysis/coexpression/filtered_models.tsv.gz"
mods_path <-snakemake@input[["mods"]]
mods <- read_tsv(mods_path) %>% 
  dplyr::select(feature.x=gene_symbol,feature.y) %>%
  mutate(is.coex=T)

#deseq_gr_path <- "results/analysis/deg/ourKD.de.grs.rds"
deseq_gr_path <- snakemake@input[["deseq_gr"]]

res <- read_rds(deseq_gr_path)$corrected %>%
  map_df(as_tibble,.id="RNAi")  %>%
  mutate(feature.x = str_remove(RNAi,"\\..+"))

res <- res %>%
  left_join(mods,by=c(feature.x="feature.x",feature="feature.y")) %>%
  mutate(class = case_when(padj >= 0.1 ~ "n.s. feature",
                           padj < 0.1 & is.na(is.coex) & str_detect(feature,"FBgn")~"gene",
                           padj < 0.1 & is.na(is.coex) & !str_detect(feature,"FBgn")~"other TE",
                          is.coex & padj < 0.1 ~ "coexpressed TE")) %>%
  mutate(class = fct_relevel(class,rev(c("coexpressed TE","other TE","gene","n.s. feature")))) %>%
  arrange(RNAi,class)


plot_volc <- function(x,y) {
  x %>%
    #filter(class %in% c("coexpressed TE","n.s. feature")) %>%
  ggplot(aes(log2FoldChange,-log10(pvalue),color=class)) +
    geom_point() +
    coord_cartesian(xlim=c(-10,10)) +
    scale_color_manual(values = c(`coexpressed TE`="red",`other TE`="blue",`gene`="darkgray",`n.s. feature`="lightgray")) +
    ggtitle(paste0(y,"-RNAi"))
}


gg_tbl <- res %>% 
  nest(data=-RNAi) %>%
  mutate(gg=map2(data,RNAi,plot_volc))


gg_list <- gg_tbl %>% dplyr::select(RNAi,gg) %>%
  deframe()

write_rds(gg_list,snakemake@output[["rds"]])



             