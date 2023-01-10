library(plotgardener)
library(tidyverse)
library(patchwork)

fig_desc <- ifelse(exists("snakemake"),snakemake@params$plot_title,"Figure 3")

theme_set(theme_classic() + theme(text = element_text(size=5),plot.title = element_text(hjust = 0.5)))

our_kds <- ifelse(exists("snakemake"),snakemake@input[["our_kds"]],"results/plots/ourKD_gsea.rds") %>%
  read_rds()

our_kds <- our_kds$plot + facet_wrap(~label, nrow=1)

if (!interactive()) pdf(snakemake@output[["pdf"]],width = 7.5, height = 10)

pageCreate(width = 7.5, height = 10, default.units = "inches", showGuides = interactive())

pa <- plotGG(plot = our_kds, x = 0.05, y=0.05, just = c("left","top"),width = 7, height=2)

plotText(label = fig_desc, fontsize = 12,
         x = 0.1, y = 9.9, just = c("left","bottom"), default.units = "inches")

if (!interactive()) dev.off()





