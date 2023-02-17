library(tidyverse)
library(clusterProfiler)

rnks <- ifelse(exists("snakemake"),snakemake@input$rds,"results/analysis/enrichment/rankings.rds") %>% read_rds()

rnks <- rnks %>%
  mutate(data = map(data, ~deframe(arrange(select(.x, gene_symbol, value),-value))))

zad_tbl <- ifelse(exists("snakemake"),snakemake@input$zad,"results/resources/zad_genes.tsv") %>% read_tsv()

zad_tbl <- zad_tbl %>% 
  mutate(gs_name="ZAD_ZNF") %>%
    dplyr::select(gs_name,ensembl_gene = gene_symbol)

pirna_tbl <- ifelse(exists("snakemake"),snakemake@input$pirna,"results/resources/pirna_pathway.tsv") %>% read_tsv()

pirna_tbl <- pirna_tbl %>% pivot_longer(cols = contains("in."),names_to = "gs_name", values_to = "is") %>%
  filter(is) %>%
  mutate(gs_name = str_remove(gs_name,"in.")) %>%
  dplyr::select(gs_name,ensembl_gene = gene_symbol)

grps <- read_tsv("http://ftp.flybase.net/releases/FB2022_04/precomputed_files/genes/gene_group_data_fb_2022_04.tsv.gz",skip=7)

t2g <- grps %>%
  dplyr::select(gs_name=FB_group_name,ensembl_gene=Group_member_FB_gene_symbol) %>%
  distinct() %>%
  drop_na() %>%
  group_by(gs_name) %>%
  filter(n()>10) %>%
  ungroup() %>%
  bind_rows(pirna_tbl) %>%
  bind_rows(zad_tbl)

#t2g %>% filter(str_detect(gs_name,"HIGH MOBIL")) %>% pull(gs_name) %>% unique %>% walk(message)

possibly_gsea <- possibly(function(.x,.y) {
  # note that per the function help page, the pv cutoff is for adjusted values
  # also - as written here it expects already unnested inputs
  st <- ifelse(str_detect(.y,"abs"),"pos","std")
  GSEA(.x, TERM2GENE = t2g,seed=2022,pvalueCutoff = 0.1,minGSSize = 10,pAdjustMethod = "BH", scoreType=st)
},otherwise = NULL)

res <- rnks %>%
  mutate(gsea = map2(data, metric, possibly_gsea)) %>%
  mutate(gsea.tidy = map(gsea, as_tibble))

write_rds(res, snakemake@output$rds)



