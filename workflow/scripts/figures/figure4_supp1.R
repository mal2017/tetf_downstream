library(plotgardener)
library(tidyverse)
library(patchwork)

theme_set(theme_classic() + theme(text = element_text(size=5),plot.title = element_text(hjust = 0.5)))

our_kds <- ifelse(exists("snakemake"),snakemake@input[["our_kds_all"]],"results/plots/pirna_genes_in_our_kd_all.rds") %>%
  read_rds()

rda <- ifelse(exists("snakemake"), snakemake@input[["motifs"]],
                   "results/plots/combined_motif_analysis_table.rda")
load(rda)

slaidina <- ifelse(exists("snakemake"), snakemake@input[["pirna_slaidina"]],
                   "results/plots/pirna_gene_expression_slaidina.rds") %>%
  read_rds()

slaidina <- (slaidina$czech + ggtitle("Czech et al. 2013")) /
  (slaidina$handler + ggtitle("Handler et al. 2013")) &
  theme(legend.key.size = unit(5,"pt"), 
        legend.text = element_text(size=unit(5,"pt")),
        legend.position = c(0.05,1.05), legend.justification = c(0,1), legend.title = element_blank())

our_kds <- our_kds + guides(color="none") + scale_color_grey(start=0.6, end=0.3) + facet_wrap(~RNAi)

motif_alns <- tab %>%
  filter(padj < 0.1) %>%
  group_by(CONSENSUS) %>%
  slice_min(padj) %>%
  ungroup() %>%
  pull(gg) %>%
  Reduce(`+`,.)

if (!interactive()) pdf(snakemake@output[["pdf"]],width = 7.5, height = 10)

pageCreate(width = 7.5, height = 10, default.units = "inches", showGuides = interactive())

pa <- plotGG(plot = our_kds, x = 0.05, y=0.05, just = c("left","top"),width = 7, height=4)

pb <- plotGG(plot = slaidina, x = 0.15, y=4.5, just = c("left","top"),width = 2.65, height=4)

pc <- plotGG(plot=motif_alns, x = 3.15, y=4.5, just = c("left","top"),width = 4.15, height=4)

plotText(label = "A", fontsize = 7,
         x = 0.125, y = 0.125, just = "center", default.units = "inches")

plotText(label = "B", fontsize = 7,
         x = 0.125, y = 4.625, just = "center", default.units = "inches")

plotText(label = "C", fontsize = 7,
         x = 3.125, y = 4.625, just = "center", default.units = "inches")

plotText(label = "Figure 4 - Supplement 1", fontsize = 12,
         x = 0.1, y = 9.9, just = c("left","bottom"), default.units = "inches")

if (!interactive()) dev.off()





