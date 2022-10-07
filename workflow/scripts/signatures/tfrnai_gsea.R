library(clusterProfiler)
library(tidyverse)
library(furrr)

threads <- 6
threads <-  snakemake@threads

plan(multisession, workers = threads)

#coex_path <- "results/analysis/coexpression/filtered_models.tsv.gz"
#res_path <- "results/analysis/deg/s2rplus.res.tsv.gz"

coex_path <- snakemake@input[["coex"]]
res_path <- snakemake@input[["deg"]]

lms <- read_tsv(coex_path)
res <- read_tsv(res_path)

# could prefilter to remove clearly uninteresting kds
# res <- res %>% group_by(comparison) %>% filter(sum(adj.P.Val < 0.1) > 25 )

t2g <- lms %>% 
  #filter(model == "female_model_01") %>%
  #filter(pvalue < 0.05) %>% # this pvalue is not for significance of the model - only for if the coefs are more extreme than expected
  dplyr::select(gs_name=gene_symbol,ensembl_gene=feature.y) %>%
  distinct() %>%
  arrange(gs_name)

gsea.tbl <- res %>%
  dplyr::select(comparison,feature,score=t) %>%
  arrange(-score) %>%
  nest(c(feature,score)) %>%
  mutate(gsea = future_map(data,~GSEA(deframe(.x), TERM2GENE = t2g,seed=2022,pvalueCutoff = 2),.options=furrr_options(seed=2022))) %>%
  mutate(gsea.tidy = future_map(gsea,as_tibble))

gsea.res <- gsea.tbl %>%
  dplyr::select(comparison,gsea.tidy) %>%
  mutate(kd = comparison) %>%
  unnest(gsea.tidy) %>%
  group_by(comparison) %>%
  mutate(pvRnk = dense_rank(-log10(pvalue)),
         nesRnk = dense_rank(abs(NES))) %>%
  ungroup() %>%
  relocate(pvRnk,nesRnk) %>%
  mutate(p.adjust = p.adjust(pvalue, method="BH")) %>%
  filter(ID==kd) %>%
  #filter(p.adjust < 0.1) %>%
  dplyr::select(RNAi=comparison,TE.set=ID,setSize,NES,p.adjust, pvalue) %>%
  arrange(NES)


write_rds(gsea.tbl,snakemake@output[["all"]])
write_rds(gsea.res,snakemake@output[["pairs"]])

