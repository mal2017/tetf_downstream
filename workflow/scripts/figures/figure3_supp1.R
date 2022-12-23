library(plotgardener)
library(tidyverse)
library(patchwork)

theme_set(theme_classic() + theme(text = element_text(size=5),plot.title = element_text(hjust = 0.5)))

exemplary <- ifelse(exists("snakemake"),snakemake@input[["exemplary_scatters"]],"results/plots/exemplary_scatters.rds") %>%
  read_rds()

# edits ---------------

if (!interactive()) pdf(snakemake@output[["pdf"]],width = 7.5, height = 10)

pageCreate(width = 7.5, height = 10, default.units = "inches", showGuides = interactive())

pa <- plotGG(plot = exemplary$pan + xlab("pan expression") + ylab("transposon expression"), x = 0.25, y=0.05, just = c("left","top"),width = 7, height=9.25)

plotText(label = "Figure 3 - Supplement 1", fontsize = 12,
         x = 0.1, y = 9.9, just = c("left","bottom"), default.units = "inches")

if (!interactive()) dev.off()
