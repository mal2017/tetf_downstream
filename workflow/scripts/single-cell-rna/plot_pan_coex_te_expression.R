library(AUCell)
library(scater)
library(scran)
library(tidyverse)
library(paletteer)
library(patchwork)

sce <- ifelse(exists("snakemake"),snakemake@input[["sce"]],
              "upstream/slaidina.sce.rds") %>%
  read_rds()

sce <- sce[,!is.na(sce$active.ident)]

sigs <- ifelse(exists("snakemake"),snakemake@input[["gsea"]],
               "~/work/tetf_downstream/results/analysis/signatures/ourKD_gsea.tbl.rds") %>%
  read_rds() %>%
  filter(str_detect(comparison,"female_gonad")) %>%
  dplyr::select(kd,core_enrichment) %>%
  mutate(core_enrichment = map(core_enrichment,~str_split(.x,"/")[[1]])) %>%
  unnest(core_enrichment) %>%
  split(.,.$kd) %>%
  map(~filter(.x,core_enrichment %in% rownames(sce))) %>%
  map(pull,core_enrichment)

# --- signature
aggregated <- sumCountsAcrossFeatures(sce, sigs,
                                      exprs_values="reconstructed", average=T)

# pan info ------------------------
g_tes <- plotTSNE(sce, colour_by = I(aggregated["pan",]), point_size = 0.05, theme_size = 7)

g_pan <- plotTSNE(sce, colour_by = "FBgn0085432",by_exprs_values = "reconstructed", point_size = 0.05, theme_size = 7)


labs <- left_join(reducedDim(sce,"TSNE") %>% as_tibble(rownames="id"),
                  colData(sce) %>% as_tibble(rownames="id"), by="id") %>%
  group_by(active.ident) %>%
  summarise(across(c("V1","V2"),mean)) %>%
  rename(X=V1, Y=V2, colour_by = active.ident)

g_tes <- g_tes + geom_text(data = labs, aes(label=colour_by), size=rel(1.5))

g_pan <- g_pan + geom_text(data = labs, aes(label=colour_by), size=rel(1.5))


write_rds(list(pan_core_enrich = g_tes, pan =g_pan), snakemake@output[["rds"]])
ggsave(snakemake@output[["png"]], g_pan + g_tes)