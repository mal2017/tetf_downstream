library(tidyverse)
library(nullranges)

#each_gene_kmer_dist_path <- "results/analysis/direct_binding/each_gene_te_kmer_dist.tbl.rds"
each_gene_kmer_dist_path <- snakemake@input[["each_gene"]]

#within_tf_nr_path <- "results/analysis/direct_binding/within_tf_nullranges.rds"
within_tf_nr_path <- snakemake@input[["rds"]]

res <- read_rds(within_tf_nr_path)

# pval odds ratio scatter
g <- res %>%
  filter(estimate > 1) %>%
  ggplot(aes(estimate,-log10(p.value))) +
  geom_point(data =  . %>% filter(padj > 0.1),color="grey", size=rel(3)) +
  geom_point(data = . %>% filter(padj <= 0.1),color="red", size=rel(3)) +
  ggrepel::geom_text_repel(data =  . %>% filter(padj <=0.1 & estimate > 1),
                           aes(label=TF), max.overlaps = 20, fontface="italic") +
  xlab("odds ratio")
# 


####
# final tes in groups with sequence similarity
also_related_tes <- read_rds(each_gene_kmer_dist_path) %>%
  filter(padj < 0.1 & in.dist < matched.dist.mean) %>%
  group_by(feature.x) %>%
  slice_min(p.value,n=1) %>%
  left_join(res,.,by=c(TF="gene_symbol"),suffix=c(".overlap",".sim"))


g2 <- also_related_tes %>%
  filter(estimate > 1) %>% 
  ggplot(aes(estimate,-log10(p.value.overlap))) +
  geom_point(data =  . %>% filter(padj.overlap > 0.1),color="grey", size=rel(3)) +
  geom_point(data = . %>% filter(padj.overlap <= 0.1),color="red", size=rel(3)) +
  geom_point(data = . %>% filter(padj.sim <= 0.1 & in.dist < matched.dist.mean),size=rel(5),shape=21,color="blue") +
  ggrepel::geom_text_repel(data =  . %>% filter((padj.overlap <=0.1 & estimate > 1)),
                           aes(label=TF), max.overlaps = 20, fontface="italic") +
  
  xlab("odds ratio") + ylab("-log10(p)")


write_rds(g2,snakemake@output[["rds"]])




# 
# 
# gsea <- read_rds("results/analysis/signatures/s2rplus_te_gsea.pairs.tbl.rds")
# 
# gsea %>% filter(RNAi %in% also_related_tes$TF)
# 
# 
# library(gt)
# 
# for_table <- res %>%
#   filter(padj< 0.1 & estimate > 1) %>%
#   arrange(-estimate) %>%
#   mutate(OR.rank = row_number()) %>%
#   dplyr::select(OR.rank, TF, oddsRatio = estimate,conf.low,conf.high, p.value,padj,gg) %>%
#   #filter(OR.rank <= 10 | TF %in% c("Dref","pan")) %>%
#   mutate(`Prop. bound`=NA) 
# 
# gt(for_table) %>%
#   fmt_number(columns = c(oddsRatio,conf.low,conf.high),decimals = 2) %>%
#   fmt_scientific(columns = c(p.value,padj),decimals = 1) %>%
#   text_transform(locations = cells_body(columns = "Prop. bound"),
#                  fn=function(x) {
#                    for_table$gg %>% 
#                      map(.f=~{.x + 
#                          theme_prism()+
#                          theme(axis.text = element_text(size=rel(3.2))) +
#                          xlab("") + ylab("") +
#                          scale_x_discrete(labels = c("coex","ns"))
#                      }) %>%
#                      map(ggplot_image, height=px(75))
#                  }) %>%
#   cols_hide("gg")
# 
# 
# 

# 
# ggsave("~/Downloads/CSHL22_figs/chip_or_p_scatter.svg",g,width = 7, height = 6)
# 
# 
# res %>%
#   mutate(sig = ifelse(padj < 0.1,"padj<0.1","n.s.")) %>%
#   mutate(dir = ifelse(estimate > 1,"enr.","depl.")) %>%
#   ggplot(aes(dir,fill=sig)) +
#   geom_bar() +
#   theme_prism() +
#   xlab("") +
#   scale_fill_grey(start = 0.8, end=0.4) +
#   theme(aspect.ratio = 1.25)
# 
# 
# remap_simple[["CG10431"]] %>%
#   export("~/Downloads/CG10431_sites_on_tes.bed")
