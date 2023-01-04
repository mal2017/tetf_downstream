library(tidyverse)
library(jsonlite)

# ------------------------------------------------------------------------------
# Read in data 
# ------------------------------------------------------------------------------
dat <- ifelse(exists("snakemake"), 
               snakemake@input$mods, 
               "upstream/final-models.collected-info.tsv.gz") %>%
  read_tsv()

limma_gsea <- ifelse(exists("snakemake"), 
              snakemake@input$gsea_pairs, 
              "results/analysis/signatures/s2rplus_te_gsea.pairs.tbl.rds") %>%
  read_rds()

limma <- ifelse(exists("snakemake"), 
                     snakemake@input$gsea_pairs, 
                     "results/analysis/deg/s2rplus.res.tsv.gz") %>%
  read_tsv()

pirna <- ifelse(exists("snakemake"), 
                snakemake@input$pirna, 
                "results/resources/pirna_pathway.tsv") %>%
  read_tsv()

# ------------------------------------------------------------------------------
# get data
# ------------------------------------------------------------------------------

n_tested <- limma_gsea$TE.set %>% unique() %>% length()

n_leading_edge_sig <- limma_gsea %>% filter(padj < 0.1) %>% nrow()


res <- list(n_tested=n_tested, n_leading_edge_sig=n_leading_edge_sig) %>% 
  enframe(name="statistic") %>%
  mutate(model = "all",
         stat_group = "s2rplus_gsea") %>%
  mutate(value=unlist(value))


# ------------------------------------------------------------------------------
# export 
# ------------------------------------------------------------------------------

to_write_json <- res %>%
  nest(-model,-stat_group) %>%
  arrange(model,stat_group)

# makes it wider so that each statistic is a named entry with a value pair
to_write_json <- to_write_json %>% mutate(data = map(data, ~pivot_wider(.x,names_from = "statistic", values_from = "value")))

write_json(to_write_json, snakemake@output$json, pretty=T)

