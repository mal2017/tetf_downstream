library(tidyverse)
library(rtracklayer)
library(plyranges)

res_path <- "results/analysis/direct_binding/within_tf_regioner.tsv.gz"
res_path <- snakemake@input[["tsv"]]

res <- read_tsv(res_path)

res.sig <- res %>% filter(pval < 0.05 & z > 0)

g <-res.sig %>%
  mutate(TF=fct_infreq(TF)) %>%
  mutate(TE=fct_infreq(TE)) %>%
  ggplot(aes(TF,TE,size=-log10(pval))) +
  geom_point() +
  theme(axis.text.x  = element_text(angle=45,hjust=1))

write_rds(g, snakemake@output[["rds"]])

ggsave(snakemake@output[["png"]],g)



res %>%
  filter(TF == "nau")
  filter(TE=="INE-1")
