library(tidyverse)
library(matrixStats)
library(ggplotify)
library(magick) # better for rasterization
library(ComplexHeatmap)
library(circlize)

#lm_path <- "upstream/final-models.collected-info.tsv.gz"
lm_path <- snakemake@input[["mods"]]

mods <- read_tsv(lm_path) %>% filter(significant_x)

mat <- mods %>% 
  dplyr::select(gene_symbol,feature.y) %>%
  mutate(coex = T) %>%
  distinct() %>% pivot_wider(names_from = feature.y, values_from = coex, values_fill = F) %>%
  column_to_rownames("gene_symbol") %>%
  as.matrix()

dist.genes <- dist(mat,method="binary")
dist.tes <- dist(t(mat),method="binary")

# https://scikit-learn.org/stable/auto_examples/cluster/plot_linkage_comparison.html
# https://stats.stackexchange.com/questions/195446/choosing-the-right-linkage-method-for-hierarchical-clustering
hc.genes <- hclust(dist.genes,method="ward.D")
hc.tes <- hclust(dist.tes,method="ward.D")


mode(mat) <- "integer"

col_fun = structure(c("white","black"), names = c("0","1"))

# classes_path <- "resources/Tidalbase_Dmel_TE_classifications_2015.txt"
classes_path <- snakemake@input[["classes"]]

te_anno <- read_tsv(classes_path) %>%
  dplyr::select(Flybase_name,repClass,repFamily) %>%
  distinct() %>%
  column_to_rownames("Flybase_name")

ha = HeatmapAnnotation(df=te_anno[colnames(mat),])

hms <- Heatmap(mat, top_annotation = ha,
        border = "black",
        cluster_rows = hc.genes, 
        cluster_columns = hc.tes,
        col = col_fun,
        column_names_gp = grid::gpar(fontsize = 5),
        show_row_names = F,
        name = "Coex.",
        use_raster = T,
        show_column_names = T)

write_rds(hms, snakemake@output[["rds"]])

pdf(snakemake@output[["pdf"]],onefile = T,width = 7.5)
hms
dev.off()
