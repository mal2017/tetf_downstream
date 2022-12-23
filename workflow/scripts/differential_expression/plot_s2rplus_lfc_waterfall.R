library(tidyverse)

#res_path <- "results/analysis/deg/s2rplus.res.tsv.gz"
#lkup_path <- "results/resources/gene_symbol_lookup.tsv.gz"

res_path <- snakemake@input[["deg"]]

res <- read_tsv(res_path)
#lkup <- read_tsv(lkup_path)

#res <- res %>% left_join(lkup, by=c(feature = "gene_ID"))


res <- res %>% filter(!str_detect(feature,"FBgn"))

res <- res %>% group_by(comparison) %>%
  summarize(logFC = mean(logFC))%>% arrange(logFC) %>% 
  mutate(rnk = dense_rank(logFC))

g <- res %>%
  ggplot(aes(rnk,logFC)) +
  geom_col(width = 1.01) +
  geom_col(data =  . %>% filter(comparison == "pan"), width=1.01, color="red") +
  geom_text(data = . %>% filter(comparison == "pan"), 
            aes(label=comparison),nudge_x = 10, nudge_y = -0.1,
            color="red", fontface="italic")

saveRDS(g,snakemake@output[["rds"]])
ggsave(snakemake@output[["png"]],g)

