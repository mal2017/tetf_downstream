library(tidyverse)
library(plyranges)
library(nullranges)
library(patchwork)
library(rtracklayer)
library(furrr)
library(ggprism)
library(tidyverse)

plan(multisession, workers = 4)

lms_path <- "results/resources/filtered_models.tsv.gz"
lms_path <- snakemake@input[["lms"]]

lms <- read_tsv(lms_path) %>%
  dplyr::rename(sex=model)

# -------- get TFs -------------------------------------------------------------

# remap_path <- "data/remap2022_nr_macs2_dm6_v1_0.bed.gz"
remap_path <- snakemake@input[["remap"]]

remap_simple <- import(remap_path) %>%
  mutate(ChIP = str_extract(name,".+(?=:)")) %>%
  split(.,.$ChIP) %>%
  GenomicRanges::reduce()

remap_simple <- remap_simple[names(remap_simple) %in% lms$gene_symbol]

#remap_simple <- remap_simple + 500

seqlevelsStyle(remap_simple) <- "NCBI"

# -------- chromatin -------------------------------------------------------------
#outside data downloaded from WUSTL GEP table browser aug 5, 2022 for dm6 (though most are lifted from r5).
#https://www.researchgate.net/figure/Chromatin-annotation-of-the-Drosophila-melanogaster-genome-a-A-9-state-model-of_fig6_279835262

#het_path <- "data/het_domains_r5todm6.bed.gz"
#s2_chromstate_path <- "data/chromstate_s2_9state_r5todm6.bed.gz"

het_path <- snakemake@input[["het"]]
s2_chromstate_path <- snakemake@input[["s2_chromstate"]]

het <- rtracklayer::import(het_path) %>%
  mutate(chromatin = str_extract(name,".+(?=_)"))

s2_chromstate <- rtracklayer::import(s2_chromstate_path) %>%
  mutate(s2_state = str_extract(name,".+(?=_)"))

seqlevelsStyle(het) <- "NCBI"
seqlevelsStyle(s2_chromstate) <- "NCBI"

# -------- get TEs -------------------------------------------------------------
# Which set of TE insertions to use is important to consider.. 
# One option is ISO-1 insertions via repeatmasker, as in the cell below (not run).

#insertions_path <- "data/insertions_by_strain.tsv.gz"
insertions_path <- snakemake@input[["insertions"]]

tes_p <- read_tsv(insertions_path) %>%
  filter(tidal_group == "DGRP_flies") %>%
  GRanges()

tes_p <- tes_p %>% filter(name %in% lms$feature.y)

tes_p1 <- tes_p %>%
  group_by(strain,name) %>%
  reduce_ranges() %>%
  ungroup()

strand(tes_p1) <- "*"

tes_p1.5 <- tes_p1 %>% 
  split(.,.$name) %>%
  lapply(FUN = function(x) {plyranges::mutate(x,nOL = count_overlaps(x,x))}) %>%
  GRangesList() %>%
  unlist()

n_strains <- length(unique(tes_p$strain))

# For INE-1, not all are fixed - most likely just small fragments aren't
fixed_in_dgrp <- tes_p1.5 %>% 
  filter(nOL == n_strains) %>% 
  filter(nOL == n_strains) %>% 
  split(.,.$name) %>% as.list()

# this just removes TEs not found in the reference at all
reference_hits <- tes_p %>% filter(source == "reference") %>% split(.,.$name) %>% lapply(reduce)
fixed_in_dgrp <- fixed_in_dgrp[names(fixed_in_dgrp) %in% names(reference_hits)]

# this only retains insertions that are 'fixed_in_dgrp' (ie nOL > x)
fixed_in_ref <- map(names(fixed_in_dgrp),.f=function(x){ 
  print(x)
  subsetByOverlaps(reference_hits[[x]],fixed_in_dgrp[[x]]) %>%
    mutate(name = x)
})

tes_p2 <-  plyranges::select(unlist(GRangesList(fixed_in_ref)), repeat_element = name) %>%
  split(.,.$repeat_element) %>%
  GenomicRanges::reduce() %>% 
  unlist() %>%
  mutate(.,repeat_element = names(.))

tes_p2 <- tes_p2[seqnames(tes_p2) %in% seqlevels(remap_simple)]

strand(tes_p2) <- "*" 

tes <- tes_p2

# -------- get GC -------------------------------------------------------------
library(BSgenome.Dmelanogaster.UCSC.dm6)

seqlevelsStyle(Dmelanogaster) <- "NCBI" 

tes$GC <- letterFrequency(getSeq(Dmelanogaster, tes), "GC", as.prob=T)[,1]


# ------------------------------------------------------------------------------
# annotate insertion
tes2 <- tes %>%
  plyranges::mutate(.,ix = 1: length(.)) %>%
  plyranges::join_overlap_left(het) %>%
  #plyranges::join_overlap_left(bg3_chromstate) %>%
  plyranges::join_overlap_left(s2_chromstate) %>%
  dplyr::select(repeat_element,GC, chromatin,s2_state,ix) %>%
  mutate_at(c("chromatin","s2_state"),replace_na,"other") %>%
  mutate(size = width) %>%
  group_by(ix) %>%
  arrange(chromatin) %>% # if any overlap with hetchrom, make sure we retain that info
  plyranges::slice(1) %>%
  ungroup() %>%
  sort()

# good spot to annotate with all TF overlaps
names(tes2) <- NULL

mcols(tes2) <- lapply(remap_simple, FUN=function(x)countOverlaps(tes2,x)) %>% 
  do.call(cbind,.,) %>%
  cbind(mcols(tes2),.)

te_oi <- "INE-1"

ins_oi <- tes2 %>% filter(repeat_element==te_oi)

matched <- matchRanges(ins_oi,
            tes2 %>% filter(repeat_element!=te_oi),
            covar = ~ size + GC + s2_state, method = "s")

nullranges::plotCovariate(matched,covar = "GC")

add_cont_tbl <-  function(mgr, TF) {
  if (is.null(mgr)) {
    return(NULL)
  }
  bound_and_coex <- sum(as.numeric(mcols(focal(mgr))[[TF]]) > 0)
  unbound_and_coex <- sum(as.numeric(mcols(focal(mgr))[[TF]]) == 0)
  bound_and_noncoex <- sum(as.numeric(mcols(mgr)[[TF]]) > 0)
  unbound_and_noncoex <- sum(as.numeric(mcols(mgr)[[TF]]) == 0)
  
  matrix(c(bound_and_coex,unbound_and_coex,bound_and_noncoex,unbound_and_noncoex),
         nrow = 2,
         dimnames = list(c("bound","unbound"),c("coex","noncoex")))
}


res <- names(remap_simple) %>%
  set_names(.,.) %>%
  map(~add_cont_tbl(matched,.x)) %>%
  enframe(name = "TF",value = "cont_mat")

res <- res %>%
  mutate(fish.test = map(cont_mat,fisher.test)) %>%
  mutate(fish.test.tidy = map(fish.test,broom::tidy)) %>%
  unnest(fish.test.tidy) %>%
  mutate(padj = p.adjust(p.value,method="BH")) %>%
  mutate(is.coex = TF %in% (lms %>% filter(feature.y==te_oi) %>% pull(gene_symbol)))

res <- res %>% left_join(lms %>% filter(feature.y==te_oi), by=c(TF="gene_symbol"))


res %>% 
  arrange(padj) %>%
  filter(estimate > 1) %>%
  pull(TF) %>%
  walk(message)

# wtf -  Gal here is mislabeled by Remap22 -  it is actually H3K27me3... 
# I traced it to this study: https://www.biorxiv.org/content/10.1101/127951v1.full.pdf
res %>% 
  filter(TF!="Gal" & TF %in% lms$gene_symbol) %>%
  #filter(padj < 0.05 & estimate > 1) %>% 
  arrange(-estimate) %>%
  group_by(padj < 0.05 & estimate > 1,is.coex) %>%
  tally()


res %>%
  #filter(is.coex) %>%
  mutate(coef.quantile = replace_na(coef.quantile,0)) %>%
  ggplot(aes(coef.quantile,-log10(p.value))) +
  geom_point()


# a bummer
res %>%
  filter(padj < 0.1) %>%
  ggplot(aes(is.coex,log2(estimate))) +
  geom_violin() + 
  theme(aspect.ratio = 3) +
  ggpubr::stat_compare_means()


# ----------get coefs ---------------------------------------------------------


sig_pairs_loose <- lms %>%
  filter(feature.y == "INE-1") %>%
  dplyr::select(sex,target=gene_symbol,userSet=feature.y,mean_estimate.qnorm) %>%
  dplyr::filter(target %in% names(remap_simple)) %>%
  dplyr::select(target, userSet) %>%
  distinct()


# -------------- do the test --------------------------------------------------
annotate_with_binding <- function(gr, lgr, maxgap=0L) {
  imap(lgr,~mutate(gr,n_olap = count_overlaps(gr,.x,maxgap = maxgap),TF=.y))
}

annot_with_coex <- function(gr,edges_to) {
  x <- mutate(gr, is.coex = repeat_element %in% edges_to)
  names(x) <- NULL
  x
}

get_matched_ranges <- function(gr, covar = ~ size + GC + s2_state, method="s", replace=F, seed=123, other_ins = other_ins) {
  set.seed(seed)
  foc <- gr
  pool <- other_ins
  if (length(foc) > 0 & length(pool) > 0) {
    return(matchRanges(focal = foc,
                       pool = pool,
                       method = method, 
                       replace = replace,
                       covar = covar))
  } else {
    return(NULL)
  }
  
}




get_pctages <- function(mat) {
  as_tibble(mat,rownames = "binding") %>% 
    pivot_longer(-binding, names_to = "coexpression",values_to = "n") %>% 
    group_by(coexpression) %>% 
    mutate(pct_bound  = n/sum(n)) %>%
    ungroup()
}

get_pct_plot <- function(tbl) {
  tbl %>% filter(binding == "bound") %>%
    ggplot(aes(coexpression,pct_bound)) +
    geom_col() +
    theme(aspect.ratio = 1)
}

tes_oi <- tes2 %>% filter(repeat_element == "INE-1")

other_ins <- tes2 %>% filter(repeat_element != "INE-1")

wip <- annotate_with_binding(tes_oi, as.list(remap_simple),maxgap = 1) %>%
  imap(~annot_with_coex(.x,unique(dplyr::filter(sig_pairs_loose,target == .y)$userSet)))


wip2 <- wip %>%
  enframe(name = "TF",value = "insertions") %>%
  mutate(mgr = map(insertions,possibly(get_matched_ranges,NULL))) %>%
  mutate(cont_tbl = map(mgr,add_cont_tbl)) %>%
  mutate(fisher = map(cont_tbl, possibly(fisher.test,NULL))) %>%
  mutate(fisher_tbl = map(fisher,broom::tidy)) %>%
  mutate(pct_tbl = map(cont_tbl, possibly(get_pctages,NULL))) %>%
  mutate(gg = map(pct_tbl, possibly(get_pct_plot,NULL)))


res <- wip2 %>%
  unnest(fisher_tbl) %>%
  filter(estimate > 1) %>%
  mutate(padj = p.adjust(p.value,method= "BH"))
