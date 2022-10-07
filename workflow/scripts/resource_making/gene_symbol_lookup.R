library(tidyverse)

gene_universe <- read_tsv("http://ftp.flybase.net/releases/FB2022_04/precomputed_files/genes/fbgn_fbtr_fbpp_expanded_fb_2022_04.tsv.gz", skip = 4)

lkup <- gene_universe %>% dplyr::select(gene_ID, gene_symbol, annotation_ID,gene_type) %>% distinct()

write_tsv(lkup,snakemake@output[["tsv"]])