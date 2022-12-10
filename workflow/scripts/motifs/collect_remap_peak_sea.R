library(tidyverse)

#seas <- Sys.glob("results/analysis/motifs/sea_remap_peaks/*/sea.tsv")

# get list of snakemake inputs
seas <- snakemake@params[["seas"]] 

seas <- seas  %>%
  map_df(read_tsv, comment="#")


# mutate DB column to contain the name of the directory combined.meme is in
seas <- seas %>% 
  mutate(DB = str_extract(DB,regex("(?<=tf\\/).*?(?=\\/combined.meme)")))

# export to snakemake output as tsv
seas %>% 
  write_tsv(snakemake@output[["tsv"]])