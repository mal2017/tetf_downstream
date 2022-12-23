library(AUCell)
library(scater)
library(scran)
library(tidyverse)
library(paletteer)
library(patchwork)

pirna <- ifelse(exists("snakemake"),snakemake@input[["pirna"]],
                "results/resources/pirna_pathway.tsv") %>%
  read_tsv()

sce <- ifelse(exists("snakemake"),snakemake@input[["sce"]],
              "upstream/slaidina.sce.rds") %>%
  read_rds()

sce <- sce[,!is.na(sce$active.ident)]

sigs <- pirna %>%
  filter(gene_ID %in% rownames(sce)) %>%
  dplyr::select(gene_ID, in.Handler13, in.Czech13) %>%
  pivot_longer(c(in.Handler13, in.Czech13)) %>%
  filter(value) %>%
  mutate(name= str_remove(name,"in\\.")) %>%
  split(.,.$name) %>%
  map(~pull(.,gene_ID))


# --- signature
aggregated <- sumCountsAcrossFeatures(sce, sigs,
                                      exprs_values="reconstructed", average=T)

# pan info ------------------------
g_handler <- plotTSNE(sce, colour_by = I(aggregated["Handler13",]), point_size = 0.05, theme_size = 7)

g_czech <- plotTSNE(sce, colour_by = I(aggregated["Czech13",]), point_size = 0.05, theme_size = 7)



labs <- left_join(reducedDim(sce,"TSNE") %>% as_tibble(rownames="id"),
                  colData(sce) %>% as_tibble(rownames="id"), by="id") %>%
  group_by(active.ident) %>%
  summarise(across(c("V1","V2"),mean)) %>%
  rename(X=V1, Y=V2, colour_by = active.ident)

g_czech <- g_czech + geom_text(data = labs, aes(label=colour_by), size=rel(1.5))

g_handler <- g_handler + geom_text(data = labs, aes(label=colour_by), size=rel(1.5))


write_rds(list(czech = g_czech, handler =g_handler), snakemake@output[["rds"]])
ggsave(snakemake@output[["png"]], g_handler + g_czech)