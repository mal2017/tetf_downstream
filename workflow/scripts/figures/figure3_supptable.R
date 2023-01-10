library(plotgardener)
library(tidyverse)
library(patchwork)
library(gridExtra)

fig_desc <- ifelse(exists("snakemake"),snakemake@params$plot_title,"Figure 3 - Supplementary Table")

theme_set(theme_classic() + theme(text = element_text(size=5),plot.title = element_text(hjust = 0.5)))

kd_tab <- ifelse(exists("snakemake"),snakemake@input[["kd_info"]],"results/plots/kd_info.rds") %>%
  read_rds()


# edits ---------------

if (!interactive()) pdf(snakemake@output[["pdf"]],width = 7.5, height = 10)

pageCreate(width = 7.5, height = 10, default.units = "inches", showGuides = interactive())

pa <- plotGG(plot = kd_tab, x = -0.25, y=1.25, just = c("left","top"),width = 8, height=4.5)

plotText(label = fig_desc, fontsize = 12,
         x = 0.1, y = 9.9, just = c("left","bottom"), default.units = "inches")

if (!interactive()) dev.off()
