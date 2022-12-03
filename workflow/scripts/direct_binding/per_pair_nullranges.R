library(tidyverse)
library(plyranges)
library(nullranges)
library(furrr)

threads <- ifelse(exists("snakemake"),snakemake@threads,4)
plan(multisession, workers = threads)

lms <- ifelse(exists("snakemake"),snakemake@input[["lms"]],
                   "upstream/final-models.collected-info.tsv.gz") %>%
  read_tsv() %>%
  filter(significant_x)

ins <- ifelse(exists("snakemake"),snakemake@input[["anno_ins"]],
                   "results/resources/annotated_fixed_insertions.gr.rds") %>%
  read_rds()

# 1. get fixed insertions to other script
# 2. For all TFs, and for all TEs
# 3. need bootstraps for sure. = maybe frame as bootstrap score rather than trying to do wilcox afterwards

# add1 is just the first TF, so this is the first index in our slice
tfs <- colnames(mcols(ins))[which(colnames(mcols(ins)) == "ADD1"):ncol(mcols(ins))]

tes <- unique(ins$repeat_element)

# gets contingency table for two sets of ranges
add_cont_tbl <-  function(mgr, TF) {
  if (is.null(mgr)) {
    return(NULL)
  }
  bound_and_target <- sum(mcols(focal(mgr))[[TF]])
  unbound_and_target <- sum(!mcols(focal(mgr))[[TF]])
  bound_and_other <- sum(mcols(mgr)[[TF]])
  unbound_and_other <- sum(!mcols(mgr)[[TF]])
  
  matrix(c(bound_and_target,unbound_and_target,bound_and_other,unbound_and_other),
         nrow = 2,
         dimnames = list(c("bound","unbound"),c("target","other")))
}

# gets matched nulls for a single TF, single TE
enrich_1te_1tf <- function(TF, TE, gr = ins) {
  #message(paste(TF,TE))
  set.seed(2022)
  matched <- matchRanges(filter(gr,repeat_element==TE),
                         filter(gr,repeat_element!=TE),
                         covar = ~ size + nearest.tss, method = "s")
  
  return(add_cont_tbl(matched,TF))
}

sl <- seqlengths(ins)

ins <- ins %>%
  mutate(.,chr=as.character(seqnames(.))) %>%
  mutate(.,dist.to.end = map2_dbl(chr,start,.f=~{
    len <- sl[.x]
    
    min(.y,len-.y)
  }))


res <- expand_grid(TF = tfs,TE = tes) %>%
  semi_join(lms, by=c(TF = "gene_symbol",TE="feature.y")) %>%
  #filter(TF == "Hr39" & TE == "17.6") %>%
  #filter(TF %in% c("pan")) %>%
  #head(2) %>% # for testing
  mutate(cont_mat = future_map2(TF,TE,possibly(enrich_1te_1tf,otherwise = NULL),.progress = T,.options = furrr_options(seed = TRUE)))

res <- lms %>%
  dplyr::select(TF = gene_symbol,TE=feature.y) %>%
  distinct() %>%
  mutate(is.coex=T) %>%
  left_join(res,.) %>%
  mutate(is.coex = replace_na(is.coex,F))

res <- res %>% 
  filter(!map_lgl(cont_mat,is.null)) %>%
  mutate(fish.test = map(cont_mat,~broom::tidy(fisher.test(.x,alternative = "greater")))) %>%
  unnest(fish.test,keep_empty = T) %>%
  bind_rows(filter(res,map_lgl(cont_mat,is.null))) %>%
  ungroup() %>%
  arrange(p.value)


write_rds(res,snakemake@output[["rds"]])


#   filter(padj < 0.1 & estimate > 1) %>%
#   arrange(padj)
# 
# #nullranges::plotCovariate(matched,covar = "GC")
# 
# res <- set_names(tfs,tfs) %>%
#   map(~add_cont_tbl(matched,.x)) %>%
#   enframe(name = "TF",value = "cont_mat")
# 
# res <- res %>%
#   mutate(fish.test = map(cont_mat,fisher.test)) %>%
#   mutate(fish.test.tidy = map(fish.test,broom::tidy)) %>%
#   unnest(fish.test.tidy) %>%
#   mutate(padj = p.adjust(p.value,method="BH")) %>%
#   mutate(is.coex = TF %in% (lms %>% filter(feature.y==te_oi) %>% pull(gene_symbol)))
# 
# res <- res %>% left_join(lms %>% filter(feature.y==te_oi), by=c(TF="gene_symbol"))
# 
# 
# # wtf -  Gal here is mislabeled by Remap22 -  it is actually H3K27me3... 
# # I traced it to this study: https://www.biorxiv.org/content/10.1101/127951v1.full.pdf
# res %>% 
#   filter(TF!="Gal" & TF %in% lms$gene_symbol) %>%
#   #filter(padj < 0.05 & estimate > 1) %>% 
#   arrange(-estimate) %>%
#   group_by(padj < 0.05 & estimate > 1,is.coex) %>%
#   tally()
# 
# 
# res %>% 
#   arrange(padj) %>%
#   filter(estimate > 1) %>%
#   pull(TF) %>%
#   walk(message)
# 
# 
# res %>%
#   #filter(is.coex) %>%
#   mutate(coef.quantile = replace_na(coef.quantile,0)) %>%
#   ggplot(aes(coef.quantile,-log10(p.value))) +
#   geom_point()
# 
# 
# # a bummer
# res %>%
#   filter(padj < 0.1) %>%
#   ggplot(aes(is.coex,log2(estimate))) +
#   geom_violin() + 
#   theme(aspect.ratio = 3) +
#   ggpubr::stat_compare_means()
