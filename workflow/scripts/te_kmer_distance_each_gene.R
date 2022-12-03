library(kmer)
library(tidyverse)
library(ape)
library(furrr)

threads <- 2
threads <- snakemake@threads
plan(multisession, workers = threads)

# influenced by:
# Influenced by (this tutorial)[https://bioconductor.org/packages/devel/workflows/vignettes/fluentGenomics/inst/doc/fluentGenomics.html]

#te_fasta_path <- "data/Tidalbase_transposon_sequence.fasta"
te_fasta_path <- snakemake@input[["te_fasta"]]


#background_lms_path <- "results/analysis/coexpression/filtered_models.tsv.gz"
#lms_path <- "results/analysis/coexpression/filtered_models.tsv.gz"
lms_path <- snakemake@input[["lms"]]

# ---- get data ----------------------------------------------------------------
tes.bin <- ape::read.FASTA(te_fasta_path)

lms <- read_tsv(lms_path) %>% filter(significant_x)
#background_lms <- read_tsv(background_lms_path)

# this would be too easy if we compared to non-expressed TEs - as these are super weird
# so filter here to tonly consider TEs that make it into our coef set.
# I think this is the most stringent approach.
te.names <- names(tes.bin)[names(tes.bin) %in% lms$feature.y]
tes.bin <- tes.bin[te.names]

# -------------------------- get dist ------------------------------------------
# https://resources.qiagenbioinformatics.com/manuals/phylogenymodule/current/K_mer_based_distance_estimation.html#:~:text=K%2Dmer%20based%20distance%20estimation%20is%20an%20alternative%20to%20estimating,occuring%20in%20the%20two%20sequences.
#kmer.dist <- kdistance(tes.bin,k = 8, residues = "DNA",method = "euclidean")

#mash_dist_path <- "results/analysis/direct_binding/te_mash.txt"
mash_dist_path <- snakemake@input[["mash_dist"]]
mash_dist_tbl <- read_tsv(mash_dist_path,col_names = c("seqA","seqB","mash_dist","pval","shared_hashes"))


kmer.dist <- mash_dist_tbl %>%
  dplyr::select(seqA,seqB,mash_dist) %>%
  pivot_wider(names_from = seqB,values_from = mash_dist,values_fill = 1) %>%
  column_to_rownames("seqA") %>%
  as.dist()


# --------------- bootstrapping ------------------------------------------------

source("workflow/scripts/direct_binding/utils-kmer-dist.R")

# random number consideration: https://www.r-bloggers.com/2020/09/future-1-19-1-making-sure-proper-random-numbers-are-produced-in-parallel-processing/
set.seed(2022)
boots0 <- lms %>% 
  dplyr::select(model,feature.x,gene_symbol,feature.y) %>%
  distinct() %>%
  group_by(model,gene_symbol,feature.x) %>%
  filter(n()>=2) %>% # doesn't make sense to get intra group distance for N=1
  nest(data=feature.y)

boots1 <- boots0 %>% 
  #head(n=2) %>% #for debuggin
  group_by(model,gene_symbol,feature.x) %>%
  summarise(TEs=map(data,pull,feature.y),n_TEs=map_dbl(TEs,length)) %>%
  group_by(model) %>%
  mutate(in.dist = future_map_dbl(TEs, get_dists,.options = furrr_options(seed = TRUE))) 

boots2 <- boots1 %>% mutate(boot = future_map2(n_TEs,TEs,get_boot_dists2,.options = furrr_options(seed = TRUE)))

boots3 <- boots2 %>%
  unnest(boot) %>%  # unnest this to retrieve bootstrap reps used for plotting
  group_by(model, feature.x, gene_symbol, in.dist, n_TEs, TEs) %>%
  summarize(.,matched.dist.mean = mean(matched.dist), matched_boots = list(matched_boots), matched.dists = list(matched.dist)) %>% 
  ungroup()

boots <- boots3 %>%
  group_by(model) %>%
  mutate(stat.test = map2(in.dist,matched.dists,~broom::tidy(wilcox.test(x=unlist(.y),mu=.x)))) %>%
  unnest(stat.test) %>%
  mutate(padj = p.adjust(p.value,method = "BH")) %>%
  arrange(padj) %>%
  ungroup()

write_rds(kmer.dist,snakemake@output[["kmer_dist"]])
write_rds(boots, snakemake@output[["each_gene"]])

