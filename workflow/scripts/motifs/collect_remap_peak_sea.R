library(tidyverse)

#seas <- Sys.glob("results/analysis/motifs/sea_remap_peaks/*/")

# get list of snakemake inputs
seas <- snakemake@input[["seas"]] 

seas <- seas %>% paste0("/sea.tsv")

seas <- seas %>% set_names(.,str_extract(.,"(?<=peaks\\/).+?(?=\\/)"))

seas <- seas  %>%
  map_df(read_tsv, comment="#",.id="peak_set")


# mutate DB column to contain the name of the directory combined.meme is in
seas <- seas %>% 
  #mutate(te_group = str_extract(ID,regex(".*?(?=::)"))) %>%
  separate(ID, c("te_group","ID"),sep = "::") %>%
  relocate(te_group)

# export to snakemake output as tsv
seas %>% 
  write_tsv(snakemake@output[["tsv"]])



