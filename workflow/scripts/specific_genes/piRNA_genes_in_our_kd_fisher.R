library(tidyverse)
library(GenomicRanges)

pirna <- ifelse(exists("snakemake"), snakemake@input[["pirna"]],
                "results/resources/pirna_pathway.tsv") %>%
  read_tsv()


res <- ifelse(exists("snakemake"), snakemake@input[["res"]],
              "results/analysis/deg/ourKD.de.grs.rds") %>%
  read_rds()


oi <-rev(c("knockdown2_CG16779_female_gonad_tj_control_female_gonad_tj",
        "knockdown2_NFI_female_head_Mef2.R_control_female_head_Mef2.R",
        "knockdown2_CG16779_male_gonad_aTub_control_male_gonad_aTub"))

res2 <- res$adjusted %>% 
  map_df(as_tibble, .id="RNAi") %>% 
  filter(RNAi %in% oi) %>%
  mutate(RNAi= str_remove(RNAi,"knockdown2_")) %>%
  relocate(RNAi,feature)

res2 %>% 
  mutate(is.DE = padj < 0.1,
         is.piRNA = feature %in% pirna$gene_ID) %>%
  dplyr::select(RNAi, is.DE, is.piRNA) %>%
  mutate(across(c("is.DE","is.piRNA"),replace_na,F)) %>%
  count(RNAi, is.piRNA, is.DE) %>%
  pivot_wider(names_from = is.piRNA, values_from = n,names_glue = "is.piRNA_{is.piRNA}") %>%
  mutate(is.DE=paste0("is.DE_",is.DE)) %>%
  arrange(RNAi,desc(is.DE)) %>%
  relocate(is.piRNA_TRUE,.after=is.DE) %>%
  nest(-RNAi) %>%
  mutate(data=map(data, column_to_rownames,"is.DE")) %>%
  mutate(fisher.res = map(data,~broom::tidy(fisher.test(.x)))) %>%
  unnest(fisher.res)


write_rds(g, snakemake@output[["rds"]])
ggsave(snakemake@output[["png"]], g)