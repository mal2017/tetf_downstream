library(tidyverse)
library(ape)
library(furrr)

threads <- 2
threads <- snakemake@threads
plan(multisession, workers = threads)

# influenced by:
# Influenced by (this tutorial)[https://bioconductor.org/packages/devel/workflows/vignettes/fluentGenomics/inst/doc/fluentGenomics.html]

# get lms from snakemake
lms <- ifelse(exists("snakemake"),snakemake@input[["lms"]],"upstream/final-models.collected-info.tsv.gz") %>%
  read_tsv() %>%
  filter(significant_x)

tfs <- ifelse(exists("snakemake"),snakemake@input[["tfs"]],"resources/Drosophila_melanogaster_TF.txt") %>%
  read_tsv()

# this would be too easy if we compared to non-expressed TEs - as these are super weird
# so filter here to tonly consider TEs that make it into our coef set.
# I think this is the most stringent approach.
te.names <- unique(lms$feature.y)

# -------------------------- get dist ------------------------------------------
mash_dist_tbl <- ifelse(exists("snakemake"),snakemake@input[["mash_dist"]],"results/analysis/direct_binding/te_mash.txt") %>%
  read_tsv(col_names = c("seqA","seqB","mash_dist","pval","shared_hashes"))


kmer.dist <- mash_dist_tbl %>%
  dplyr::select(seqA,seqB,mash_dist) %>%
  pivot_wider(names_from = seqB,values_from = mash_dist,values_fill = 1) %>%
  column_to_rownames("seqA") %>%
  as.dist()

# --------------- bootstrapping ------------------------------------------------
helper_code <- "workflow/scripts/te_sequence_similarity/utils-kmer-dist.R"
stopifnot(file.exists(helper_code))
source(helper_code)

# random number consideration: https://www.r-bloggers.com/2020/09/future-1-19-1-making-sure-proper-random-numbers-are-produced-in-parallel-processing/
set.seed(2022)
boots0 <- lms %>% 
  dplyr::select(feature.x,gene_symbol,feature.y) %>%
  distinct() %>%
  group_by(gene_symbol,feature.x) %>%
  filter(n()>=2) %>% # doesn't make sense to get intra group distance for N=1
  nest(data=feature.y)

boots1 <- boots0 %>% 
  #head(n=2) %>% #for debuggin
  group_by(gene_symbol,feature.x) %>%
  summarise(TEs=map(data,pull,feature.y),n_TEs=map_dbl(TEs,length)) %>%
  ungroup() %>%
  mutate(in.dist = future_map_dbl(TEs, get_dists,.options = furrr_options(seed = TRUE))) 

nperm <- snakemake@params[["nperm"]]
bg <- 1:length(te.names) %>% set_names(.,as.character(.)) %>% map(get_boot_dists3,n=nperm)

get_pval <- function(d,n,np=nperm){
  sum(bg[[as.character(n)]] < d)/np
}

boots2 <- boots1 %>%
  mutate(pval = map2_dbl(in.dist,n_TEs,get_pval)) %>%
  mutate(padj = p.adjust(pval,method="BH"))

write_rds(kmer.dist,snakemake@output[["kmer_dist"]])
write_rds(boots2, snakemake@output[["each_gene"]])
write_rds(bg, snakemake@output[["null_model"]])

