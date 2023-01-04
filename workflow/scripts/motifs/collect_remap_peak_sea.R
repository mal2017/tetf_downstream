library(tidyverse)

#seas_fls <- Sys.glob("results/analysis/motifs/sea_remap_peaks/*/")

# get list of snakemake inputs
seas_fls <- snakemake@input[["seas"]] 

seas <- seas_fls %>% paste0("/sea.tsv")

seas <- seas %>% set_names(.,str_extract(.,"(?<=peaks\\/).+?(?=\\/)"))

seas <- seas  %>%
  #head() %>%
  map(read_tsv, comment="#")

# remove empty
seas <- seas[map_lgl(seas,~{nrow(.x) >= 1})]

seas <- seas %>% 
  map(~mutate(.x, ID=as.character(ID))) %>% 
  bind_rows(.id="peak_set")

# mutate DB column to contain the name of the directory combined.meme is in
seas <- seas %>% 
  mutate(te_group = str_extract(DB,regex("(?<=per_tf\\/).+(?=\\/combined.meme)"))) %>%
  #separate(ID, c("te_group","ID"),sep = "::") %>%
  relocate(te_group)

# export to snakemake output as tsv
seas %>% 
  write_tsv(snakemake@output[["tsv"]])



