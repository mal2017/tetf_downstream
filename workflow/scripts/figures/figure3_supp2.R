library(plotgardener)
library(tidyverse)
library(patchwork)

fig_desc <- ifelse(exists("snakemake"),snakemake@params$plot_title,"Figure 3 - Supplement 2")

theme_set(theme_classic() + theme(text = element_text(size=5),plot.title = element_text(hjust = 0.5)))

our_kd_box_coefs <- ifelse(exists("snakemake"),snakemake@input[["kd_gene_coefs_box"]],
                           "results/plots/kd_gene_coefs_boxplot.rds") %>%
  read_rds()

volcs <- ifelse(exists("snakemake"),snakemake@input[["volcs"]],
                    "results/plots/this_study_kd_deseq2.rds") %>%
  read_rds()

#waterfall <- ifelse(exists("snakemake"),snakemake@input[["waterfall"]],
#                              "results/plots/s2rplus_lfc_waterfall.rds") %>%
#  read_rds()


# edits ---------------

#waterfall <- waterfall + xlab("RNAi targets ranked by mean TE LFC")


if (!interactive()) pdf(snakemake@output[["pdf"]],width = 8.5, height = 11)

pageCreate(width = 8.5, height = 11, default.units = "inches", showGuides = interactive())

pa <- plotGG(plot = our_kd_box_coefs, x = 0.25, y=0.05, just = c("left","top"),width = 7.25, height=2)

pb <- plotGG(plot = volcs, x = 0.25, y=2.25, just = c("left","top"),width = 7.75, height=6)

#pc <- plotGG(plot = waterfall, x = 0.25, y=8.5, just = c("left","top"),width = 6.75, height=2)


plotText(label = "A", fontsize = 7,
         x = 0.25, y = 0.25, just = c("center"), default.units = "inches")

plotText(label = "B", fontsize = 7,
         x = 0.25, y = 2.25, just = c("center"), default.units = "inches")

#plotText(label = "C", fontsize = 7,
#         x = 0.25, y = 8.5, just = c("center"), default.units = "inches")


plotText(label = fig_desc, fontsize = 12,
         x = 0.1, y = 10.9, just = c("left","bottom"), default.units = "inches")

if (!interactive()) dev.off()
