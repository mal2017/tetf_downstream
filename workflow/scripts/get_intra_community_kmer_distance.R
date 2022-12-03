library(tidyverse)
library(furrr)


threads <- 2
threads <- snakemake@threads
plan(multisession, workers = threads)

#dist_path <-"results/analysis/direct_binding/te_kmer_dist_mat.rds"
dist_path <-snakemake@input[["dist"]]


#comms_path <-"results/analysis/direct_binding/coregulated_te_communities.communities_tbl.rds"
comms_path <-snakemake@input[["comms"]]

comms <- read_rds(comms_path)
kmer.dist <- read_rds(dist_path)


te.names <- comms$feature %>% unique()

source("workflow/scripts/direct_binding/utils-kmer-dist.R")

# random number consideration: https://www.r-bloggers.com/2020/09/future-1-19-1-making-sure-proper-random-numbers-are-produced-in-parallel-processing/
set.seed(2022)
boots0 <- comms %>%
  group_by(model,algorithm,community) %>%
  filter(n()>=2) %>% # doesn't make sense to get intra group distance for N=1
  nest(data=feature)

boots1 <- boots0 %>% 
  #head(n=2) %>% #for debuggin
  group_by(model,algorithm,community) %>%
  summarise(TEs=map(data,pull,feature),n_TEs=map_dbl(TEs,length)) %>%
  group_by(model) %>%
  mutate(in.dist = future_map_dbl(TEs, get_dists,.options = furrr_options(seed = TRUE))) 

boots2 <- boots1 %>% mutate(boot = map2(n_TEs,TEs,get_boot_dists2))

boots3 <- boots2 %>%
  unnest(boot) %>%  # unnest this to retrieve bootstrap reps used for plotting
  group_by(model, algorithm, community, in.dist, n_TEs, TEs) %>%
  summarize(.,matched.dist.mean = mean(matched.dist), matched_boots = list(matched_boots), matched.dists = list(matched.dist)) %>% 
  ungroup()

boots <- boots3 %>%
  mutate(stat.test = map2(in.dist,matched.dists,~broom::tidy(wilcox.test(x=unlist(.y),mu=.x)))) %>%
  unnest(stat.test) %>%
  group_by(model,algorithm) %>%
  mutate(padj = p.adjust(p.value,method = "BH")) %>%
  arrange(padj) %>%
  ungroup()


write_rds(boots, snakemake@output[["by_community"]])

# filter(boots, padj < 0.1 & in.dist < matched.dist.mean)
# 
# g_are_te_groups_more_similar_than_expected <- boots %>%
#   dplyr::select(model,algorithm,community,in.dist,matched.dist.mean) %>%
#   pivot_longer(-c(community,model,algorithm),names_to = "group",values_to = "D") %>%
#   ggplot(aes(group,D)) +
#   ggpubr::stat_compare_means() +
#   geom_boxplot() +
#   facet_grid(model~algorithm)
# 
# 
# 
# boots
# 
# 
# 
# te_fasta_path <- "data/Tidalbase_transposon_sequence.fasta"
# te_seqs <- Biostrings::readDNAStringSet(te_fasta_path)
# 
# coreg <- boots %>%
#   filter(algorithm == "leiden05") %>%
#   #filter(padj<0.1) %>%
#   arrange(-log10(p.value), in.dist/matched.dist.mean) %>%
#   unite(coreg,model,algorithm,community) %>%
#   dplyr::select(coreg,TEs) %>%
#   unnest(TEs) %>%
#   split(.,.$coreg)
# 
# 
# coreg$male_model_01_leiden05_2 %>%
#   pull(TEs) %>%
#   te_seqs[.] -> oi
# 
# Biostrings::writeXStringSet(oi,"~/Downloads/oi.fasta")
# 
# 
# library(rGADEM)
# library(BSgenome.Dmelanogaster.UCSC.dm6)
# gadem_res <- GADEM(oi,verbose = 1, genome = te_seqs ,seed = 2022,numGeneration = 10,)
# 
# consensus(gadem_res) %>% str_to_upper()
# 
# nOccurrences(gadem_res)
# 
# # could do something not unlike chromVar
# https://cran.r-project.org/web/packages/kmer/kmer.pdf
# # probably better to just use meme or homer if possible
# 
# # https://bioconductor.org/packages/3.14/bioc/vignettes/universalmotif/inst/doc/SequenceSearches.pdf