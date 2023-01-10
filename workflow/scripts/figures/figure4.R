library(plotgardener)
library(tidyverse)
library(patchwork)
library(gt)
library(gtExtras)
library(gridExtra)

fig_desc <- ifelse(exists("snakemake"),snakemake@params$plot_title,"Figure 4")

theme_set(theme_classic() + theme(text = element_text(size=5),plot.title = element_text(hjust = 0.5)))

our_kds <- ifelse(exists("snakemake"),snakemake@input[["our_kds"]],"results/plots/pirna_genes_in_our_kd.rds") %>%
  read_rds()

gtg_path <- ifelse(exists("snakemake"), snakemake@input[["motif_gt"]],
                   "results/plots/combined_motif_analysis_table.png")

#load(file = gtg_path)
gtg <- magick::image_read(gtg_path) %>%
  magick::image_ggplot()

slaidina <- ifelse(exists("snakemake"), snakemake@input[["slaidina"]],
                   "results/plots/pan_coex_te_expression_slaidina.rds") %>%
  read_rds()

slaidina <- slaidina$pan + slaidina$pan_core_enrich &
  theme(legend.key.size = unit(5,"pt"), 
        legend.text = element_text(size=unit(5,"pt")),
        legend.position = c(0.05,1.05), legend.justification = c(0,1), legend.title = element_blank())

slaidina <- slaidina + scale_size(range=0)

our_kds <- our_kds + guides(color="none") + scale_color_grey(start=0.6, end=0.3) + facet_wrap(~RNAi2, nrow=1)


if (!interactive()) pdf(snakemake@output[["pdf"]],width = 8.5, height = 11)

pageCreate(width = 8.5, height = 11, default.units = "inches", showGuides = interactive())

pa <- plotGG(plot = our_kds, x = 0.5, y=0.15, just = c("left","top"),width = 7.4, height=2)

pb <- plotGG(plot = gtg, x = 0.2125, y=2.2, just = c("left","top"),width = 2.75, height=3)

pc <- plotGG(plot=slaidina, x = 3.3, y=2.5, just = c("left","top"),width = 4.8, height=2.7)

plotText(label = "A", fontsize = 7,
         x = 0.5, y = 0.125, just = "center", default.units = "inches")

plotText(label = "B", fontsize = 7,
         x = 0.25, y = 2.25, just = "center", default.units = "inches")

plotText(label = "C", fontsize = 7,
         x = 3.25, y = 2.25, just = "center", default.units = "inches")

plotText(label = fig_desc, fontsize = 12,
         x = 0.1, y = 10.9, just = c("left","bottom"), default.units = "inches")

if (!interactive()) dev.off()





