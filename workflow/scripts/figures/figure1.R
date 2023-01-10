library(plotgardener)
library(tidyverse)
library(png)
library(ComplexHeatmap)
library(patchwork)

fig_desc <- ifelse(exists("snakemake"),snakemake@params$plot_title,"Figure 1")

theme_set(theme_classic() + theme(text = element_text(size=5),plot.title = element_text(hjust = 0.5)))

cartoon_01 <- ifelse(exists("snakemake"),snakemake@input[["cartoon"]],"resources//simple_linmod_methods_v01.png") %>%
  readPNG(native = T)

intro_heats <- ifelse(exists("snakemake"),snakemake@input[["heatmaps"]],"results/plots/intro_heatmaps.gg-list.rds") %>%
  read_rds()

intro_scatter <- ifelse(exists("snakemake"),snakemake@input[["ncoex_scatter"]],"results/plots/intro_ncoex_scatter.rds") %>%
  read_rds()


# https://phanstiellab.github.io/plotgardener/articles/guides/incorporating_ggplots.html
# make these plottable as ggplot
intro_heats_male <- grid.grabExpr(draw(intro_heats$male_model_01))
intro_heats_female <- grid.grabExpr(draw(intro_heats$female_model_01))

# make this one plot
intro_scatter <- (Reduce(`+`,intro_scatter) + plot_layout(guides="collect")) & theme(legend.position = "bottom")

if (!interactive()) pdf(snakemake@output[["pdf"]],width = 7.5, height = 10)

pageCreate(width = 7.5, height = 10, default.units = "inches", showGuides = interactive())

f1a <- plotRaster(cartoon_01, x = 0.625, y=0.05, just = c("left","top"),width = 1.25, height=1.5)

f1b <- plotGG(plot = intro_heats_female, x = 2+(3/8), y=0.05, just = c("left","top"),width = 2.75-(3/8), height=2.9)
f1c <- plotGG(plot = intro_heats_male, x = 4.75, y=0.05, just = c("left","top"),width = 2.75, height=2.9)

f1d <- plotGG(intro_scatter, x = 0, y=1.56, just = c("left","top"),width = 2+(3/8), height=1.44)


plotText(label = "A", fontsize = 7,
         x = 0.25, y = 0.25, just = "center", default.units = "inches")

plotText(label = "B", fontsize = 7,
         x = 2.5, y = 0.25, just = "center", default.units = "inches")

plotText(label = "C", fontsize = 7,
         x = 5, y = 0.25, just = "center", default.units = "inches")


plotText(label = "D", fontsize = 7,
         x = 0.25, y = 1.625, just = "center", default.units = "inches")

plotText(label = fig_desc, fontsize = 12,
         x = 0.1, y = 9.9, just = c("left","bottom"), default.units = "inches")

if (!interactive()) dev.off()





