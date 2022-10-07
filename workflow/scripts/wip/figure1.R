library(plotgardener)
library(tidyverse)
library(png)

theme_set(theme_classic() + theme(text = element_text(size=7)))

# figure 1 - main --------------------------------------------------------------
# https://phanstiellab.github.io/plotgardener/articles/guides/incorporating_ggplots.html
make_gg <- . %>%  draw() %>% grid.grabExpr()

cartoon_01 <- readPNG("data/simple_linmod_methods_v01.png",native = T)
intro_heats <- read_rds("results/plots/intro_heatmaps.gg-list.rds")
intro_upset <- read_rds("results/plots/plot_intro_upsetplot.rds")
intro_scatter <- read_rds("results/plots/plot_intro_ncoex_scatter.rds")
intro_hists <- read_rds("results/plots/plot_intro_histograms.rds")

if (!interactive()) pdf("~/Downloads/figure1.pdf",width = 8.5, height = 11)

pageCreate(width = 8.5, height = 6, default.units = "inches", showGuides = interactive())

f1a <- plotRaster(cartoon_01, x = 0.6125, y=0.25, just = c("left","top"),width = 2.25, height=3)

f1b <- plotGG(intro_heats$combined, x = 3, y=0.1, just = c("left","top"),width = 5, height=3)

f1c <- plotGG(intro_scatter, x = 0.5, y=3.5, just = c("left","top"),width = 2.5, height=2.5)

f1d <- plotGG(intro_upset$plot, x = 3, y=3.5, just = c("left","top"),width = 5, height=2.5)


plotText(label = "A", fontsize = 7,
         x = 0.5, y = 0.25, just = "center", default.units = "inches")

plotText(label = "B", fontsize = 7,
         x = 3.25, y = 0.25, just = "center", default.units = "inches")

plotText(label = "C", fontsize = 7,
         x = 0.5, y = 3.5, just = "center", default.units = "inches")

plotText(label = "D", fontsize = 7,
         x = 3.25, y = 3.5, just = "center", default.units = "inches")


if (!interactive()) dev.off()





