library(tidyverse)
library(rtracklayer)
library(Biostrings)
library(BSgenome)
library(GenomicRanges)

fac <-  ifelse(exists("snakemake"),snakemake@params[["tf"]],
               "pan")

# import bed file of repeatmasker features with rtracklayer
rpm <- ifelse(exists("snakemake"),snakemake@input[["rpm"]],
              "resources/plus-repeats.repeatmasked.fixednames.bed") %>%
  import()

strand(rpm) <- "*"

# import bed file of remap peaks with rtracklayer
grl <- ifelse(exists("snakemake"),snakemake@input[["bed"]],
              "results/resources/remap.gr.rds") %>%
    read_rds()

# keep only seqames 2R, 2L, 3R, 3L, X, Y, for each gr in the grl
# keep only peaks not on masked repeats
# and drop unused seqlevels
grl <- grl %>% 
    as.list() %>%
    map(~plyranges::filter(.x,seqnames %in% c("2R","2L","3R","3L","X","Y"))) %>%
    map(~plyranges::filter_by_non_overlaps(.x, rpm)) %>%
    GRangesList()

seqlevels(grl) <- seqlevelsInUse(grl)


# import genome fasta with rtrocklayer
fa <- ifelse(exists("snakemake"),snakemake@input[["fa"]],
             "resources/dmel-all-chromosome-r6.41.fasta.gz") %>%
    import()

# edit names of fasta to be only the first word before any whitespace
names(fa) <- str_split(names(fa),"\\s+") %>% map_chr(1)

gr <- grl %>% unlist() %>% plyranges::filter_by_non_overlaps(grl[[fac]]) %>% reduce()

# get sequence of grl from fa and save to DNAStringSet object at variable seqs
seqs <- BSgenome::getSeq(fa,gr)

names(seqs) <- as.character(gr)

# get output dir from snakemake
ofa <- ifelse(exists("snakemake"),snakemake@output[["fa"]],
               "results/analysis/motifs/inverted/remap_peaks/pan.fasta")

export(seqs, ofa, format = "fasta")

