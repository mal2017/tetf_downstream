library(magrittr)
library(tidyverse)

#tsvs <- Sys.glob("results/analysis/motifs/xstreme_per_tf/*/sea_out/sea.tsv")
tsvs <- snakemake@input[["tsvs"]] %>% paste0("/sea_out/sea.tsv")

tsvs <- tsvs %>% set_names(.,str_extract(.,"(?<=tf\\/).+(?=\\/sea_out)"))

res <- tsvs %>% 
  map_df(read_tsv, comment="#",.id = "te_group")

write_tsv(res, snakemake@output[["tsv"]])
