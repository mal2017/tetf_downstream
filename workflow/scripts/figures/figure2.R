library(plotgardener)
library(tidyverse)
library(patchwork)

theme_set(theme_classic() + theme(text = element_text(size=5),plot.title = element_text(hjust = 0.5)))

pirna_box <- ifelse(exists("snakemake"),snakemake@input[["pirna_box"]],"results/figs/pirna_genes_in_lms.rds") %>%
  read_rds()

tf_gsea <- ifelse(exists("snakemake"),snakemake@input[["gg_gsea_volc"]],"results/figs/gene_group_gsea.volc.rds") %>%
  read_rds()

pirna_box <- pirna_box & scale_fill_grey(start = 0, end = 0.9)
tf_gsea <- tf_gsea + scale_color_grey(start=0.6,end=0.1) + guides(color=F)

if (!interactive()) pdf(snakemake@output[["pdf"]],width = 7.5, height = 2)

pageCreate(width = 7.5, height = 2, default.units = "inches", showGuides = interactive())

pa <- plotGG(plot = pirna_box, x = 0.25, y=0.05, just = c("left","top"),width = 5, height=2)

pb <- plotGG(plot = tf_gsea, x = 5.25, y=0.05, just = c("left","top"),width = 2, height=2)

plotText(label = "A", fontsize = 7,
         x = 0.25, y = 0.25, just = "center", default.units = "inches")

plotText(label = "B", fontsize = 7,
         x = 5.125, y = 0.25, just = "center", default.units = "inches")

if (!interactive()) dev.off()





