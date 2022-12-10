library(regioneR)
library(rtracklayer)
library(tidyverse)

threads <- ifelse(exists("snakemake"),snakemake@threads,8)
options(mc.cores=threads)

remap <- readRDS(ifelse(exists("snakemake"),snakemake@input[["remap"]],"results/resources/remap.gr.rds"))

# each unique insertion -  retain only fixed (177 strains as far as TIDAL is concerned)
ins <- import(ifelse(exists("snakemake"),snakemake@input[["unique_ins"]],"results/resources/dgrp_tidal_insertions.unique.bb"))
ins <- ins[ins$score == max(ins$score)]

lms <- read_tsv(ifelse(exists("snakemake"),snakemake@input[["mods"]],"upstream/final-models.collected-info.tsv.gz")) %>%
  filter(significant_x)

remap <- remap[seqnames(remap) %in% seqnames(ins)]
ins <- ins[seqnames(ins) %in% seqlevels(remap)]


seqlevels(ins, pruning.mode="coarse") <- seqlevelsInUse(ins)
seqlevels(remap, pruning.mode="coarse") <- seqlevelsInUse(remap)

fa <- import(ifelse(exists("snakemake"),snakemake@input[["genome_fa"]],"resources/dmel-all-chromosome-r6.41.fasta.gz"))
names(fa) <- names(fa)%>% str_extract(".+(?=\\s+type)")

genome_df <- fa[seqlevels(ins),] %>% seqlengths() %>% 
  enframe(name = "seqnames",value = "end") %>%
  mutate(start=1) %>%
  dplyr::select(seqnames,start,end) %>%
  as.data.frame()

ins <- ins %>% split(.,.$name)

res2tibble <- function(x) {
  tibble(pval = x$pval, 
         z = x$zscore,
         alternative = x$alternative,
         observed = x$observed)
}

res <- expand_grid(TE=names(ins),TF=names(remap)) %>%
  semi_join(lms, by=c(TF = "gene_symbol",TE="feature.y")) %>%
  #head(5) %>%
  mutate(result = map2(TE,TF,.f=function(TE,TF) {
    set.seed(123)
    z <- overlapPermTest(A= remap[[TF]], B =ins[[TE]],verbose=T,
                         genome=genome_df, 
                         alternative = "greater",
                         non.overlapping =T,
                         per.chromosome =T,
                         mc.set.seed=FALSE,
                         nperm=10)
    
    res2tibble(z$numOverlaps)
  }))


res <- res %>% unnest(result)

write_tsv(res,snakemake@output[["tsv"]])






