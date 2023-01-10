library(plotgardener)
library(tidyverse)
library(patchwork)
library(gt)

fig_desc <- ifelse(exists("snakemake"),snakemake@params$plot_title,"Figure2 - Supplement 1")

theme_set(theme_classic() + theme(text = element_text(size=5),plot.title = element_text(hjust = 0.5)))

gg_gsea_table <- ifelse(exists("snakemake"),snakemake@input[["gg_gsea_table"]],"results/plots/gene_group_gsea.table.png")

gtg <- magick::image_read(gg_gsea_table) %>%
  magick::image_ggplot()

pirna_hist <- ifelse(exists("snakemake"),snakemake@input[["pirna_hist"]],"results/plots/pirna_score_hist.rds") %>%
  read_rds()

if (!interactive()) pdf(snakemake@output[["pdf"]],width = 7.5, height = 10)

pageCreate(width = 7.5, height = 10, default.units = "inches", showGuides = interactive())

#pa <- plotGG(plot = pirna_hist, x = 1.25, y=0.05, just = c("left","top"),width = 5, height=2)

pb <- plotGG(plot = gtg, x = 0.25, y=0.5, just = c("left","top"),width = 7.25, height=7)

#plotText(label = "A", fontsize = 7,
#         x = 1, y = 0.25, just = "center", default.units = "inches")

#plotText(label = "B", fontsize = 7,
#         x = 0.75, y = 2.5, just = "center", default.units = "inches")

plotText(label = fig_desc, fontsize = 12,
         x = 0.1, y = 9.9, just = c("left","bottom"), default.units = "inches")

if (!interactive()) dev.off()





