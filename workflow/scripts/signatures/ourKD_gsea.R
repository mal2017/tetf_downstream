library(clusterProfiler)
library(tidyverse)
library(GenomicRanges)

gene_universe <- read_tsv("http://ftp.flybase.net/releases/FB2022_04/precomputed_files/genes/fbgn_fbtr_fbpp_expanded_fb_2022_04.tsv.gz", skip = 4)

lkup <- gene_universe %>% dplyr::select(gene_ID, gene_symbol) %>% distinct()

#lms_path <- "upstream/models.collected-info.tsv.gz"
lms_path <-snakemake@input[["mods"]]

lms <- read_tsv(lms_path)
lms <- lms %>% filter(significant_x)

# res_path <- "results/analysis/deg/ourKD.de.grs.rds"
res_path <- snakemake@input[["res"]]

res <- read_rds(res_path) %>%
  map_df(~map_df(.x,as_tibble,.id="comparison"),.id="adjustment")

res <- mutate(res, kd = str_extract(comparison,"(?<=knockdown2_).+?(?=_)"))

res <- res %>% mutate(kd=if_else(kd == "NFI","NfI",kd))

res <- res %>% left_join(lkup,by=c(kd="gene_symbol"))

t2g <- lms %>% 
  dplyr::select(gs_name=gene_symbol,ensembl_gene=feature.y) %>%
  distinct() %>%
  arrange(gs_name)

possibly_gsea <- possibly(function(.x,.y){GSEA(deframe(.x), TERM2GENE = filter(t2g,gs_name == .y),seed=2022,pvalueCutoff = 1,minGSSize = 5)},otherwise = NULL)

gsea.tbl <- res %>%
  filter(adjustment == 'adjusted') %>%
  dplyr::select(comparison,kd,feature,score=stat) %>%
  arrange(-score) %>%
  nest(c(feature,score)) %>%
  mutate(gsea = map2(data,kd,possibly_gsea)) %>%
  mutate(gsea.tidy = map(gsea,as_tibble)) %>%
  unnest(gsea.tidy) %>%
  mutate(p.adjust = p.adjust(p.adjust,method="BH"))


# "results/ourKD.gsea.tbl.rds"
saveRDS(gsea.tbl, snakemake@output[["rds"]])

# gsea.tbl %>%
#     dplyr::select(comparison,kd,gsea.tidy) %>%
#   unnest(gsea.tidy) %>%
#   filter(ID==kd) %>%
#   mutate(p.adjust = p.adjust(p.adjust,method="BH")) %>%
#   filter(p.adjust < 0.1)
#   group_by(comparison) %>%
#   mutate(pvRnk = dense_rank(-log10(pvalue)),
#          nesRnk = dense_rank(abs(NES))) %>%
#   ungroup() %>%
#   relocate(pvRnk,nesRnk) %>%
#  
#   filter(p.adjust < 0.1) %>%
#   dplyr::select(RNAi=comparison,TE.set=ID,setSize,NES,pvalue,p.adjust) %>%
#   arrange(NES) %>%
#   gt::gt()
#   
#  
# 
# gsea.tbl$gsea.plt
