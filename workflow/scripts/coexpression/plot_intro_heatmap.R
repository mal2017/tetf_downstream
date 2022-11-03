library(tidyverse)
library(matrixStats)
library(ggplotify)
library(magick) # better for rasterization
library(ComplexHeatmap)
library(circlize)

#lkup_path <- "results/resources/gene_symbol_lookup.tsv.gz"
lkup_path <- snakemake@input[["lkup"]]

#merged_mods_path <- "results/resources/merged_models.tsv.gz"
merged_mods_path <- snakemake@input[["merged_mods"]]

#filtered_mods_path <- "results/resources/filtered_models.tsv.gz"
filtered_mods_path <- snakemake@input[["filtered_mods"]]

#extreme_mods_path <- "results/resources/extreme_models.tsv.gz"
extreme_mods_path <- snakemake@input[["extreme_mods"]]

lkup <- read_tsv(lkup_path)

filtered_mods <- read_tsv(filtered_mods_path)
extreme_mods <- read_tsv(extreme_mods_path)
mods <- read_tsv(merged_mods_path)

mods <- filter(mods,feature.x %in% extreme_mods$feature.x)
filtered_mods <- filter(filtered_mods, feature.x %in% extreme_mods$feature.x)


mats <- mods %>%
  split(.,.$model) %>%
  map(~dplyr::select(.,feature.x,feature.y,mean_estimate.qnorm)) %>%
  map(distinct) %>%
  map(pivot_wider,names_from="feature.y", values_from = mean_estimate.qnorm, values_fill=0) %>%
  map(column_to_rownames,"feature.x") %>%
  map(as.matrix)

filtered_mats <- filtered_mods %>%
  split(.,.$model) %>%
  map(~dplyr::select(.,feature.x,feature.y,mean_estimate.qnorm)) %>%
  map(distinct) %>%
  map(pivot_wider,names_from="feature.y", values_from = mean_estimate.qnorm, values_fill=0) %>%
  map(column_to_rownames,"feature.x") %>%
  map(as.matrix)


uni_mat <- mods %>%
  filter(model %in% c("male_model_01","female_model_01")) %>%
  dplyr::select(.,feature.x,feature.y,mean_estimate.qnorm) %>%
  group_by(feature.x,feature.y) %>%
  summarise(mean_estimate.qnorm = mean(mean_estimate.qnorm), .groups = "drop") %>%
  pivot_wider(names_from="feature.y", values_from = mean_estimate.qnorm, values_fill=0) %>%
  column_to_rownames("feature.x") %>%
  as.matrix

filtered_uni_mat <- filtered_mods %>%
  filter(model %in% c("male_model_01","female_model_01")) %>%
  dplyr::select(.,feature.x,feature.y,mean_estimate.qnorm) %>%
  group_by(feature.x,feature.y) %>%
  summarise(mean_estimate.qnorm = mean(mean_estimate.qnorm), .groups = "drop") %>%
  pivot_wider(names_from="feature.y", values_from = mean_estimate.qnorm, values_fill=0) %>%
  column_to_rownames("feature.x") %>%
  as.matrix

mats[["combined"]] <- uni_mat
filtered_mats[["combined"]] <- filtered_uni_mat

mats <- mats %>% 
  map(~{.x[!rowAlls(.x,value=0),]})

filtered_mats <- filtered_mats %>% 
  map(~{.x[!rowAlls(.x,value=0),]})

# https://scikit-learn.org/stable/auto_examples/cluster/plot_linkage_comparison.html
# https://stats.stackexchange.com/questions/195446/choosing-the-right-linkage-method-for-hierarchical-clustering
hcs_from_mat <- . %>% dist() %>% hclust(method = "ward.D")

gene_hcs <- filtered_mats %>% #map(t) %>% 
  #map(abs) %>% 
  #map(scale,center=F) %>% 
  map(hcs_from_mat)

te_hcs <- filtered_mats %>% map(t) %>%
  #map(abs) %>% 
  #map(scale,center=F) %>% 
  map(hcs_from_mat)


col_fun = colorRamp2(colors = c("blue", "white", "red"),breaks = c(-10,0,10))


plot_hm <- function(m, t,g,row_cap = "", col_cap = "") {
  row_cap <- paste0(row_cap," (n=",nrow(m),")")
  col_cap <- paste0(col_cap," (n=",ncol(m),")")
  xh <- Heatmap(m, border = "black",
          cluster_rows = g, 
          na_col = "white",
          row_title_gp = gpar(fontsize=12), 
          column_title_gp = gpar(fontsize=12),
          cluster_columns = t,
          row_title = row_cap,
          column_title = col_cap,heatmap_legend_param = list(grid_width=unit(2,"mm"),
                                                             grid_height = unit(2,"mm"), 
                                                             labels_gp = gpar(fontsize = 12),
                                                             title_gp = gpar(fontsize = 12, fontface = "bold")),
          col = col_fun,
          show_row_names = F,name = "score",
          use_raster = T,
          show_column_names = F)
  return(xh)
}

hms <- list(combined = plot_hm(filtered_mats$combined, te_hcs$combined,gene_hcs$combined,row_cap = "genes", col_cap = "TEs"),
            male_model_01 = plot_hm(filtered_mats$male_model_01,  te_hcs$male_model_01,gene_hcs$male_model_01,row_cap = "genes", col_cap = "TEs"),
            female_model_01 = plot_hm(filtered_mats$female_model_01, te_hcs$female_model_01,gene_hcs$female_model_01, row_cap = "genes", col_cap = "TEs"))

n_hits_and_extreme <- filtered_mods %>%
  #filter(coef.quantile > 0.9) %>%
  group_by(feature.x,feature.y,gene_symbol) %>%
  slice_max(coef.quantile,n = 1, with_ties = F) %>%
  group_by(feature.x, gene_symbol) %>%
  dplyr::summarise(n=n(),nExtreme = sum(coef.quantile > 0.9),.groups = "drop") %>%
  column_to_rownames("feature.x") %>%
  .[rownames(filtered_mats$combined),]


hms$combined <- hms$combined +
  #Heatmap(n_hits_and_extreme$nExtreme, name="N extreme", width=unit(5,"mm")) +
  rowAnnotation(link = anno_mark(at = which(n_hits_and_extreme$nExtreme > 31), 
                                 labels = n_hits_and_extreme[n_hits_and_extreme$nExtreme > 31,"gene_symbol"], 
                                 link_width = unit(20,"mm"),
                                 labels_gp = gpar(fontsize = 12), padding = unit(1, "mm")))

hms$combined

write_rds(hms, snakemake@output[["heats"]])
write_rds(te_hcs, snakemake@output[["te_hcs"]])
write_rds(gene_hcs, snakemake@output[["gene_hcs"]])







# palette_len <- 255
# colors <- colorRampPalette( c("blue","white","red"))(palette_len)
# 
# myBreaks <- c(seq(-50,0, 
#                   length.out=ceiling(palette_len/2) + 1), 
#               seq(0.01, 50, 
#                   length.out=floor(palette_len/2)))
# 
# if(!interactive()) pdf(NULL)
# 
# xx <- pheatmap::pheatmap(mats$female_model_01,
#                          cluster_rows = gene_hcs$female_model_01,
#                          cluster_cols = te_hcs$female_model_01,
#                          color = colors, 
#                          breaks=myBreaks,
#                          treeheight_col = 25, treeheight_row = 25,
#                          show_rownames = F, show_colnames = F, scale = "none")
# if(!interactive()) pdf(NULL)
# xy <- pheatmap::pheatmap(mats$male_model_01,
#                          cluster_rows = gene_hcs$male_model_01,
#                          cluster_cols = te_hcs$male_model_01,
#                          color = colors, 
#                          breaks=myBreaks,
#                          treeheight_col = 25, treeheight_row = 25,
#                          show_rownames = F, show_colnames = F, scale = "none")
# if(!interactive()) pdf(NULL)
# combined <- pheatmap::pheatmap(mats$combined,
#                            cluster_rows = gene_hcs$combined,
#                          cluster_cols = te_hcs$combined,
#                          color = colors, 
#                          breaks=myBreaks,
#                          treeheight_col = 25, treeheight_row = 25,
#                          show_rownames = F, show_colnames = F, scale = "none")
# 
# gxx <- as.ggplot(xx)
# gxy <- as.ggplot(xy)
# gcombined <- as.ggplot(combined)

# write_rds(list(male_model_01 = gxy, female_model_01 = gxx, combined = gcombined), snakemake@output[["heats"]])


