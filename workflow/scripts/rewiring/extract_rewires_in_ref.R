library(tidyverse)
library(GenomicRanges)
library(rtracklayer)
library(plyranges)
library(Biostrings)
library(BSgenome.Dmelanogaster.UCSC.dm6)

seqlevelsStyle(BSgenome.Dmelanogaster.UCSC.dm6) <- "NCBI"

ref_ins <- import("upstream/reference_insertions.bed")

# get candidate rewiring events involving INE-1, nautilus, and CG11902
rew <- read_rds("results/analysis/rewiring/candidate_rewiring_events.rds")

rew <- rew %>%
  filter(feature.y == "INE-1") %>%
  filter(gene_symbol %in% c("nau","CG11902")) %>%
  filter(padj < 0.1) %>%
  pull(id) %>%
  unique() %>%
  str_remove("INE-1\\.") %>%
  GRanges()

# get all INE-1s in the ref, annotate by whether it appears to be
# involved in rewiring with nautilus and CG11902
ine1 <- ref_ins %>% 
  filter(name == "INE-1") %>%
  mutate(.,rewiring_event = . %over% rew)

names(ine1) <- as.character(ine1) %>% paste0(ifelse(ine1$rewiring_event,"_RWE",""))

sq <- getSeq(BSgenome.Dmelanogaster.UCSC.dm6, ine1,)

export(sq,"~/work/ine1_phylo_v2/ine1_in_ref.fasta")
# now get sequence of all INE-1s in the host
# note which these are
# see if they cluster together by MASH, IQTREE, or if there are enriched motifs