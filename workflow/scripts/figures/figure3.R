library(plotgardener)
library(tidyverse)
library(patchwork)

theme_set(theme_classic() + theme(text = element_text(size=5),plot.title = element_text(hjust = 0.5)))

our_kds <- ifelse(exists("snakemake"),snakemake@input[["our_kds"]],"results/plots/ourKD_gsea.rds") %>%
  read_rds()

s2rnai <- ifelse(exists("snakemake"),snakemake@input[["s2rnai"]],"results/plots/plot_tfrnai_gsea.plot_list.rds") %>%
  read_rds() %>%
  .$ne_vs_p

s2rnai <- s2rnai + scale_color_grey(start = 0.6, end = 0.1) + theme(legend.position = c(1,1), legend.justification = c("right","top"))
our_kds


if (!interactive()) pdf(snakemake@output[["pdf"]],width = 7.5, height = 2)

pageCreate(width = 7.5, height = 2, default.units = "inches", showGuides = interactive())

pa <- plotGG(plot = our_kds$plot, x = 0.05, y=0.05, just = c("left","top"),width = 4, height=2)

pb <- plotGG(plot = s2rnai, x = 4.1, y=0.05, just = c("left","top"),width = 3, height=2)

plotText(label = "A", fontsize = 7,
         x = 0.125, y = 0.125, just = "center", default.units = "inches")

plotText(label = "B", fontsize = 7,
         x = 4.125, y = 0.125, just = "center", default.units = "inches")

if (!interactive()) dev.off()





