library(plotgardener)
library(tidyverse)
library(patchwork)
library(gridExtra)
library(ggpubr)

theme_set(theme_classic() + theme(text = element_text(size=5),plot.title = element_text(hjust = 0.5)))

prox <- ifelse(exists("snakemake"),snakemake@input[["proximity"]],
               "results/analysis/specific_genes/remap_peaks_near_pirna_genes_contingency.kd-chip-intersect.rds") %>%
  read_rds()

gtg <- ggtexttable(prox)


# edits
if (!interactive()) pdf(snakemake@output[["pdf"]],width = 7.5, height = 10)

pageCreate(width = 7.5, height = 10, default.units = "inches", showGuides = interactive())

pa <- plotGG(plot = gtg, x = -0.25, y=1.25, just = c("left","top"),width = 8, height=4.5)

plotText(label = "Figure 4 - Supplementary Table 1", fontsize = 12,
         x = 0.1, y = 9.9, just = c("left","bottom"), default.units = "inches")

if (!interactive()) dev.off()





