library(AnnotationDbi)
library(GenomicFeatures)
library(tidyverse)
library(BSgenome.Dmelanogaster.UCSC.dm6)

# flybase == NCBI
genome <- BSgenome.Dmelanogaster.UCSC.dm6
seqlevelsStyle(genome) <- "NCBI"

txdb <- ifelse(exists("snakemake"),snakemake@input[["gtf"]],
               "resources/dmel-all-r6.41.gtf.gz") %>%
  makeTxDbFromGFF(chrominfo = seqinfo(genome))

saveDb(txdb, file = snakemake@output[["txdb"]])
