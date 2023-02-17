library(biomaRt)
library(tidyverse)
library("org.Dm.eg.db")

ensembl <- useEnsembl(biomart = "genes", dataset = "dmelanogaster_gene_ensembl")

#listFilters(ensembl)
#listAttributes(ensembl)

ids <- getBM(filters = "interpro",
    values = "IPR012934",
    attributes = c("ensembl_gene_id"),
    mart = ensembl)

ids <- ids %>% 
  as_tibble %>%
  distinct() %>%
  mutate(gene_symbol = mapIds(org.Dm.eg.db, 
                              keys=ensembl_gene_id, 
                              column="SYMBOL", 
                              keytype="ENSEMBL",
                              multiVals="first")) %>%
  distinct()

write_tsv(ids,snakemake@output$tsv)
