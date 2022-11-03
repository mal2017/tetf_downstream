library(tidyverse)
library(plyranges)
library(rtracklayer)

#lms_path <- "results/resources/filtered_models.tsv.gz"
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
  filter(nOL >= n_strains) %>% 
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
  group_by(ix) %>%
  plyranges::slice(1) %>% # retain only 1 locus in the case that multiple annotation regions overlap the same insertion
  ungroup() %>%
  sort()

# good spot to annotate with all TF overlaps
names(tes2) <- NULL

# could try jaccard dist -> clust -> check cluster concordance w/ similarity
chip_overlap_mat <- lapply(remap_simple, FUN=function(x){countOverlaps(tes2,x)>0}) %>% 
  do.call(cbind,.,)

# wtf -  Gal here is mislabeled by Remap22 -  it is actually H3K27me3... 
# I traced it to this study: https://www.biorxiv.org/content/10.1101/127951v1.full.pdf
chip_overlap_mat <- chip_overlap_mat[,colnames(chip_overlap_mat)!="Gal"]

mcols(tes2) <- cbind(mcols(tes2),chip_overlap_mat)

write_rds(tes2,snakemake@output[["rds"]])

