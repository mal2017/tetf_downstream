library(clusterProfiler)
library(tidyverse)
library(furrr)
plan(multisession, workers = 2)

gene_universe <- read_tsv("http://ftp.flybase.net/releases/FB2022_04/precomputed_files/genes/fbgn_fbtr_fbpp_expanded_fb_2022_04.tsv.gz", skip = 4)

lkup <- gene_universe %>% dplyr::select(gene_ID, gene_symbol) %>% distinct()

#lms_path <- "results/analysis/coexpression/filtered_models.tsv.gz"
lms_path <-snakemake@input[["filtered_mods"]]

lms <- read_tsv(lms_path)

# res_path <- "results/analysis/deg/ourKD.de.grs.rds"
res_path <- snakemake@input[["res"]]
res <- read_rds(res_path) %>%
  map_df(~map_df(.x,as_tibble,.id="comparison"),.id="adjustment")

res <- mutate(res, kd = str_extract(comparison,".+(?=\\.)|.+$"))

res <- res %>% left_join(lkup,by=c(kd="gene_symbol"))

t2g <- lms %>% 
  #filter(pvalue < 0.05) %>%
  dplyr::select(gs_name=gene_symbol,ensembl_gene=feature.y) %>%
  distinct() %>%
  arrange(gs_name)

gsea.tbl <- res %>%
  filter(adjustment == 'corrected') %>%
  dplyr::select(comparison,kd,feature,score=stat) %>%
  arrange(-score) %>%
  nest(c(feature,score)) %>%
  #head(1) %>%
  mutate(gsea = map(data,~GSEA(deframe(.x), TERM2GENE = t2g,seed=2022,pvalueCutoff = 1))) %>%
  mutate(gsea.tidy = map(gsea,as_tibble)) #%>%
  #mutate(gsea.plt = map2(gsea,kd,~gseaplot(.x,geneSetID = .y,title = paste("TEs coex. w/",.y,"after",.y,"RNAi"), by="runningScore")))


# "results/ourKD.gsea.tbl.rds"
saveRDS(gsea.tbl, snakemake@output[["rds"]])


# gsea.tbl %>%
#   dplyr::select(comparison,kd,gsea.tidy) %>%
#   unnest(gsea.tidy) %>%
#   group_by(comparison) %>%
#   mutate(pvRnk = dense_rank(-log10(pvalue)),
#          nesRnk = dense_rank(abs(NES))) %>%
#   ungroup() %>%
#   relocate(pvRnk,nesRnk) %>%
#   filter(ID==kd)%>%
#   filter(p.adjust < 0.1) %>%
#   dplyr::select(RNAi=comparison,TE.set=ID,setSize,NES,pvalue,p.adjust) %>%
#   arrange(NES) %>%
#   gt::gt()
#   
#  
# 
# gsea.tbl$gsea.plt
