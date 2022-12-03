library(tidyverse)
library(plyranges)
library(rtracklayer)

#fa_fl <- "resources/dmel-all-chromosome-r6.41.fasta.gz"
fa_fl <- snakemake@input[["genome_fa"]]
fa <- import(fa_fl)

names(fa) <- names(fa) %>% str_extract(".+?(?=\\s)")

#insertions_path <- "upstream/DGRP_flies.tsv.gz"
insertions_path <- snakemake@input[["ins"]]

ins <- read_tsv(insertions_path)

ins <- GRanges(ins)

ins <- ins %>% mutate(strain = str_extract(strain,"DGRP_\\d+"))

seqlengths(ins) <- seqlengths(fa)[seqlevels(ins)]

# Now get the reduced set in such a way that doesn't double count overlaps from the same strain
# strain 1: ------       --------
# strain 2:      ----------
# should yield...
#           -------------------- penetrance = 2
# and not penetrance=3

tes_all <- ins %>%
  group_by(name,strain) %>%
  reduce_ranges() %>%
  ungroup()

tes_penetrance <- tes_all %>%
  group_by(name) %>%
  reduce_ranges(score = plyranges::n_distinct(strain)) %>%
  ungroup()

ins <- mutate(ins,name=paste(name,strain,source,sep=".")) %>% plyranges::select(name,score)

export(ins,snakemake@output[["all_ins"]],format="bb")
export(tes_penetrance,snakemake@output[["penetrance"]],format="bb")