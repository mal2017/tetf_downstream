library(tidyverse)
library(tidytext)
#library(org.Dm.eg.db)

# https://support.bioconductor.org/p/129049/
#GOterms <- mapIds(org.Dm.eg.db, keys(org.Dm.eg.db, "GO"), "ENSEMBL", "GO", multiVals = "list")

#GOterms[c("GO:0034584","GO:0034587")]

#mods_path <- "upstream/final-models.collected-info.tsv.gz"
mods_path <- snakemake@input[["mods"]]
mods <- read_tsv(mods_path) #%>% filter(significant_x)

#pirna_path <- "results/resources/pirna_pathway.tsv"
pirna_path <- snakemake@input[["pirna"]]
pirna_tbl <- read_tsv(pirna_path)

cts <- count(mods, model, gene_symbol,relationship=if_else(significant_x,"sig","ns")) %>%
  pivot_wider(names_from = relationship, values_from = n, values_fill = 0) %>%
  mutate(Czech2013 = gene_symbol %in% filter(pirna_tbl,in.Czech13)$gene_symbol) %>%
  mutate(Handler2013 = gene_symbol %in% filter(pirna_tbl,in.Handler13)$gene_symbol)

##
# plots begin
#

# ------------------------------------------------------------------------------
# do we find anything special abt piRNA pathway genes in terms of number of assoc?
# the total number of reproducible TE associations is higher for piRNA genes
s_ntotal.ecdf <- cts %>%
  mutate(across(contains("2013"), as.factor)) %>%
  split(.,.$model) %>%
  map_df(
    ~{
      bind_rows(Handler2013 = broom::tidy(ks.test(sig~Handler2013,data=.,alternative="two.sided")),
                Czech2013 = broom::tidy(ks.test(sig~Czech2013,data=.,alternative="two.sided")),.id="publication")
      },
    .id="model")

g_ntotal.ecdf <- filtered_cts %>%
  pivot_longer(c(Czech2013,Handler2013),names_to = "publication",values_to = "identified") %>%
  filter(publication == "Handler2013") %>%
  mutate(gene.type = if_else(identified,"piRNA","other")) %>%
  mutate(sex=ifelse(str_detect(model,"female"),"female","male")) %>%
  ggplot(aes(sig,linetype=gene.type,color=sex)) +
  stat_ecdf() +
  facet_wrap(~publication)

saveRDS(list(stat = s_ntotal.ecdf, plot = g_ntotal.ecdf),snakemake@output[["ntotal_ecdf_rds"]])
ggsave(g_ntotal.ecdf,snakemake@output[["ntotal_ecdf_png"]])