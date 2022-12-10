library(AnnotationDbi)
library(GenomicFeatures)
library(tidyverse)

txdb <- ifelse(exists("snakemake"),snakemake@input[["gtf"]],
               "resources/dmel-all-r6.41.gtf.gz") %>%
  makeTxDbFromGFF()


saveDb(txdb, file = snakemake@output[["txdb"]])
