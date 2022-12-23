library(tidyverse)
library(SummarizedExperiment)
library(gridExtra)
library(ggpubr)

dds_path <- ifelse(exists("snakemake"),
       snakemake@input[["dds"]],"results/analysis/deg/ourKD.dds.list.rds")

dds <- read_rds(dds_path)

dds <- dds$adjusted

dat <- colData(dds) %>%
  as_tibble() %>%
  dplyr::select(target = knockdown, `BDSC stock #`=stock, 
                driver, `driver BDSC stock #`=driver_stock,
                tissue, batch) %>%
  mutate(tissue = str_replace_all(tissue,"_"," ")) %>%
  distinct()

to_italic <- dat$target %>% str_detect(.,"control") %>% {!.} %>% which %>% {.+1}

gtg<-ggtexttable(dat, rows=NULL)

gtg <- gtg %>% table_cell_font(row = to_italic,column = 1, face="italic")

saveRDS(gtg,snakemake@output[["rds"]])
