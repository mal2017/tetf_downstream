library(clusterProfiler)
library(tidyverse)
library(furrr)

threads <- 6
threads <-  snakemake@threads

plan(multisession, workers = threads)

#coex_path <- "upstream/final-models.collected-info.tsv.gz"
#res_path <- "results/analysis/deg/s2rplus.res.tsv.gz"

coex_path <- snakemake@input[["coex"]]
res_path <- snakemake@input[["deg"]]

lms <- read_tsv(coex_path) %>% filter(significant_x)
res <- read_tsv(res_path)

# could prefilter to remove clearly uninteresting kds
# res <- res %>% group_by(comparison) %>% filter(sum(adj.P.Val < 0.1) > 25 )

t2g <- lms %>% 
  dplyr::select(gs_name=gene_symbol,ensembl_gene=feature.y) %>%
  distinct() %>%
  group_by(gs_name) %>%
  ungroup() %>%
  arrange(gs_name)

possibly_gsea <-  possibly(function(.x,.y) {
  GSEA(deframe(.x), TERM2GENE = filter(t2g, gs_name == .y),seed=2022,pvalueCutoff = 2)
}, NULL)

gsea.tbl <- res %>%
  dplyr::select(comparison,feature,score=t) %>%
  arrange(-score) %>%
  nest(c(feature,score)) %>%
  #filter(comparison == "pan") %>%
  #mutate(gsea = future_map(data,possibly_gsea,.options=furrr_options(seed=2022))) %>%
  mutate(gsea = map2(data,comparison,possibly_gsea)) %>%
  mutate(gsea.tidy = map(gsea,as_tibble))

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
  #filter(ID==kd) %>%
  #filter(p.adjust < 0.1) %>%
  dplyr::select(RNAi=comparison,TE.set=ID,setSize,NES,pvalue) %>%
  mutate(padj = p.adjust(pvalue, method="BH")) %>%
  arrange(padj)


write_rds(gsea.tbl,snakemake@output[["all"]])
write_rds(gsea.res,snakemake@output[["pairs"]])

