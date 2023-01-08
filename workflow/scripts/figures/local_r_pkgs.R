library(tidyverse)
library(plyranges)
library(clusterProfiler)
library(fgsea)
library(DESeq2)
library(scater)
library(scran)
library(tesseract)
library(DescTools)
library(GenomicRanges)
library(ggpubr)
library(sessioninfo)

sessioninfo::external_info() %>% 
  as_tibble %>%
  pivot_longer(everything(),names_to = "lib", values_to = "version") %>% 
  write_tsv(snakemake@output$external_session_info)

sessioninfo::platform_info()  %>% 
  as_tibble %>% 
  pivot_longer(everything(),names_to = "item", values_to = "info") %>% 
  write_tsv(snakemake@output$platform_info)

sessioninfo::package_info("loaded")  %>%
  as_tibble %>% 
  write_tsv(snakemake@output$package_info)


