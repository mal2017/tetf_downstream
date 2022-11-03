library(tidyverse)
library(plyranges)
library(nullranges)

library(furrr)

plan(multisession, workers = 4)

#lms_path <- "results/resources/filtered_models.tsv.gz"
lms_path <- snakemake@input[["lms"]]

lms <- read_tsv(lms_path)

#anno_ins_path <- "results/resources/annotated_fixed_insertions.gr.rds"
anno_ins_path <- snakemake@input[["anno_ins"]]

# 1. splitto get fixed insertions to other script
# 2. For all TFs, and for all TEs
# 3. need bootstraps for sure. = maybe frame as bootstrap score rather than trying to do wilcox afterwards

tes2 <- read_rds(anno_ins_path)

tfs <- colnames(mcols(tes2))[which(colnames(mcols(tes2)) == "ADD1"):ncol(mcols(tes2))]


tes <- unique(tes2$repeat_element)


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


# for a single TF, single TE
enrich_1te_1tf <- function(TF, TE, gr = tes2) {
  #message(paste(TF,TE))
  set.seed(2022)
  matched <- matchRanges(filter(gr,repeat_element==TE),
                         filter(gr,repeat_element!=TE),
                         covar = ~ size + GC + s2_state, method = "s")
  
  return(add_cont_tbl(matched,TF))
}


res <- expand_grid(TF = tfs,TE = tes) %>%
  #head(2) %>% # for testing
  mutate(cont_mat = future_map2(TF,TE,possibly(enrich_1te_1tf,otherwise = NULL),.progress = T))

res <- lms %>%
  dplyr::select(TF = gene_symbol,TE=feature.y) %>%
  distinct() %>%
  mutate(is.coex=T) %>%
  left_join(res,.) %>%
  mutate(is.coex = replace_na(is.coex,F))

res <- res %>% 
  filter(!map_lgl(cont_mat,is.null)) %>%
  mutate(fish.test = map(cont_mat,~broom::tidy(fisher.test(.x)))) %>%
  unnest(fish.test,keep_empty = T) %>%
  bind_rows(filter(res,map_lgl(cont_mat,is.null))) %>%
  ungroup()


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
