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

paletteer::palettes_d_names %>%
  filter(type=="qualitative") %>%
  arrange(-length) %>%
  filter(!novelty)

# pan info ------------------------

g <- plotTSNE(sce, colour_by = "active.ident",
         point_size = 0.05, theme_size = 7) +
  scale_color_paletteer_d("ggsci::default_igv") +
  guides(colour = guide_legend(override.aes = list(size=5)))
  theme(legend.key.size = unit(20,"pt"))

labs <- left_join(reducedDim(sce,"TSNE") %>% as_tibble(rownames="id"),
          colData(sce) %>% as_tibble(rownames="id"), by="id") %>%
  group_by(active.ident) %>%
  summarise(across(c("V1","V2"),mean)) %>%
  rename(X=V1, Y=V2, colour_by = active.ident)


g + geom_text(data = labs, aes(label=colour_by))

#write_rds(list(pan_core_enrich = g_tes, pan =g_pan), snakemake@output[["rds"]])
#ggsave(snakemake@output[["png"]], g_pan + g_tes)