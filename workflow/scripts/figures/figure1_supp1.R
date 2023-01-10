library(plotgardener)
library(tidyverse)
library(patchwork)

fig_desc <- ifelse(exists("snakemake"),snakemake@params$plot_title,"Figure 1 - Supplement 1")

theme_set(theme_classic() + theme(text = element_text(size=5),plot.title = element_text(hjust = 0.5)))

var_exp_box <- ifelse(exists("snakemake"),snakemake@input[["var_exp_box"]],"results/plots/variance_explained_overview_boxplot.rds") %>%
  read_rds()

alluvial <- ifelse(exists("snakemake"),snakemake@input[["alluvial"]],"results/plots/sig_coefs_alluvial.rds") %>%
  read_rds()

mf_corr <- ifelse(exists("snakemake"),snakemake@input[["mf_corr"]],"results/plots/male_vs_female_signal.rds") %>%
  read_rds()

sex_specific_bar <-  ifelse(exists("snakemake"),snakemake@input[["sex_specific_bar"]],"results/plots/sex_specific_barplot.rds") %>%
  read_rds()

filtering_barplot <-  ifelse(exists("snakemake"),snakemake@input[["filtering_barplot"]],"results/plots/n_models.rds") %>%
  read_rds()

coex_te_hist <-  ifelse(exists("snakemake"),snakemake@input[["coex_te_hist"]],"results/plots/plot_intro_histograms.rds") %>%
  read_rds()

consistency <-  ifelse(exists("snakemake"),snakemake@input[["consistency"]],"results/plots/consistency.rds") %>%
  read_rds()


# edits ---------------

filtering_barplot <- filtering_barplot + theme(legend.position = c(0.95,1.05), legend.justification = c(1,1), legend.title = element_blank())

alluvial <- alluvial + guides(fill="none")

coex_te_hist <- coex_te_hist + facet_wrap(~model) + xlab("coexpression score")

var_exp_box <- var_exp_box + theme(axis.text.x = element_text(angle=45, hjust=1))

mf_corr$estimate.qnorm <- mf_corr$estimate.qnorm + 
  xlab("male coexpression score") + ylab("female coexpression score")
  #ylim(c(-0.5, 0.5)) + xlim(c(-0.5,0.5)) 

if (!interactive()) pdf(snakemake@output[["pdf"]],width = 7.5, height = 10)

pageCreate(width = 7.5, height = 10, default.units = "inches", showGuides = interactive())

pa <- plotGG(plot = filtering_barplot, x = 0.25, y=0.05, just = c("left","top"),width = 2.75, height=2)

pb <- plotGG(plot = alluvial, x = 3.25, y=0.05, just = c("left","top"),width = 4.125, height=2)

pc <- plotGG(plot = var_exp_box, x = 0.25,  y=2.25, just = c("left","top"),width = 2.75, height=2)

pd <- plotGG(plot = coex_te_hist, x = 3.25,  y=2.25, just = c("left","top"),width = 4.125, height=2)

pe <- plotGG(plot = mf_corr$estimate.qnorm, x = 0.25,  y=4.3, just = c("left","top"),width = 2.7, height=2.7)

pf <- plotGG(plot = sex_specific_bar, x = 4.375-1.25,  y=4.5, just = c("left","top"),width = 4.25, height=2)

plotText(label = "A", fontsize = 7,
         x = 0.25, y = 0.25, just = "center", default.units = "inches")

plotText(label = "B", fontsize = 7,
         x = 3.25, y = 0.25, just = "center", default.units = "inches")

plotText(label = "C", fontsize = 7,
         x = 0.25, y = 2.25, just = "center", default.units = "inches")

plotText(label = "D", fontsize = 7,
         x = 3.25, y = 2.25, just = "center", default.units = "inches")

plotText(label = "E", fontsize = 7,
         x = 0.25, y = 4.5, just = "center", default.units = "inches")

plotText(label = "F", fontsize = 7,
         x = 4.5-1.25, y = 4.5, just = "center", default.units = "inches")

plotText(label = fig_desc, fontsize = 12,
         x = 0.1, y = 9.9, just = c("left","bottom"), default.units = "inches")

if (!interactive()) dev.off()