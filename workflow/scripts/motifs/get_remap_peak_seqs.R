library(tidyverse)
library(rtracklayer)
library(Biostrings)
library(BSgenome)
library(GenomicRanges)

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

# get sequence of grl from fa and save to DNAStringSet object at variable seqs
seqs <- as.list(grl) %>% map(~{z <- BSgenome::getSeq(fa,.x);names(z) <- as.character(.x);z})


# get output dir from snakemake
odir <- ifelse(exists("snakemake"),snakemake@output[["odir"]],
               "results/analysis/motifs/remap_peaks/")

# create odir if it doesn't exist
dir.create(odir,recursive = T)

# save seqs to fasta files with one fasta file for each element in seqs
for (i in names(seqs)) {
    message(i)
    export(seqs[[i]],paste0(odir,"/",i,".fasta"))
}

