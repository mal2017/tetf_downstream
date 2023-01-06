library(tidyverse)
library(matrixStats)
library(ggplotify)
library(magick) # better for rasterization
library(ComplexHeatmap)
library(circlize)


#merged_mods_path <- "upstream/final-models.collected-info.tsv.gz"
merged_mods_path <- snakemake@input[["mods"]]

merged_mods <- read_tsv(merged_mods_path)

filtered_mods <- merged_mods %>% filter(significant_x)

extreme_mods <- filtered_mods %>% filter(coef.quantile > .999)

mods <- filter(merged_mods,feature.x %in% extreme_mods$feature.x)

filtered_mods <- filter(filtered_mods, feature.x %in% extreme_mods$feature.x)

mats <- mods %>%
  split(.,.$model) %>%
  map(~dplyr::select(.,feature.x,feature.y,estimate.qnorm)) %>%
  map(distinct) %>%
  map(pivot_wider,names_from="feature.y", values_from = estimate.qnorm, values_fill=0) %>%
  map(column_to_rownames,"feature.x") %>%
  map(as.matrix)

filtered_mats <- filtered_mods %>%
  split(.,.$model) %>%
  map(~dplyr::select(.,feature.x, feature.y, estimate.qnorm)) %>%
  map(distinct) %>%
  map(pivot_wider,names_from="feature.y", values_from = estimate.qnorm, values_fill=0) %>%
  map(column_to_rownames,"feature.x") %>%
  map(as.matrix)

uni_mat <- mods %>%
  dplyr::select(.,feature.x, feature.y, estimate.qnorm) %>%
  group_by(feature.x, feature.y) %>%
  summarise(estimate.qnorm = mean(estimate.qnorm), .groups = "drop") %>%
  pivot_wider(names_from="feature.y", values_from = estimate.qnorm, values_fill=0) %>%
  column_to_rownames("feature.x") %>%
  as.matrix

filtered_uni_mat <- filtered_mods %>%
  dplyr::select(., feature.x, feature.y, estimate.qnorm) %>%
  group_by(feature.x, feature.y) %>%
  summarise(estimate.qnorm = mean(estimate.qnorm), .groups = "drop") %>%
  pivot_wider(names_from="feature.y", values_from = estimate.qnorm, values_fill=0) %>%
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

gene_hcs <- mats %>% #map(t) %>% 
  #map(abs) %>% 
  #map(scale,center=F) %>% 
  map(hcs_from_mat)

te_hcs <- mats %>% map(t) %>%
  #map(abs) %>% 
  #map(scale,center=F) %>% 
  map(hcs_from_mat)

col_fun = colorRamp2(colors = c("blue", "white", "red"),breaks = c(0.5*min(mats$combined),0,0.5*max(mats$combined)))

plot_hm <- function(m, t,g,row_cap = "", col_cap = "",plot_leg=T) {
  row_cap <- paste0(row_cap," (n=",nrow(m),")")
  col_cap <- paste0(col_cap," (n=",ncol(m),")")
  xh <- Heatmap(m, border = "black",
                show_heatmap_legend = plot_leg,
          cluster_rows = g, 
          column_dend_height = unit(5,"mm"),
          row_dend_width = unit(5,"mm"),
          na_col = "white",
          row_title_gp = gpar(fontsize=5), 
          column_title_gp = gpar(fontsize=5),
          cluster_columns = t,
          row_title = row_cap,
          column_title = col_cap,heatmap_legend_param = list(grid_width=unit(2,"mm"),
                                                             grid_height = unit(2,"mm"), 
                                                             labels_gp = gpar(fontsize = 5),
                                                             title_gp = gpar(fontsize = 5, fontface = "bold")),
          col = col_fun,
          show_row_names = F,name = "score",
          use_raster = T,
          show_column_names = F)
  return(xh)
}

hms <- list(combined = plot_hm(mats$combined, te_hcs$combined,gene_hcs$combined,row_cap = "genes", col_cap = "TEs"),
            male_model_01 = plot_hm(mats$male,  te_hcs$male, gene_hcs$male,row_cap = "genes", col_cap = "TEs"),
            female_model_01 = plot_hm(mats$female, te_hcs$female, gene_hcs$female, row_cap = "genes", col_cap = "TEs",plot_leg = F))

n_hits_and_extreme <- filtered_mods %>%
  #filter(coef.quantile > 0.9) %>%
  group_by(feature.x,feature.y,gene_symbol) %>%
  slice_max(coef.quantile,n = 1, with_ties = F) %>%
  group_by(feature.x, gene_symbol) %>%
  dplyr::summarise(n=n(),nExtreme = sum(coef.quantile > 0.99),.groups = "drop") %>%
  column_to_rownames("feature.x") %>%
  .[rownames(filtered_mats$combined),]

thresh <- 0.5

hms$combined <- hms$combined +
  #Heatmap(n_hits_and_extreme$nExtreme, name="N extreme", width=unit(5,"mm")) +
  rowAnnotation(link = anno_mark(at = which(n_hits_and_extreme$nExtreme > thresh*max(n_hits_and_extreme$nExtreme)), 
                                 labels = n_hits_and_extreme[n_hits_and_extreme$nExtreme > thresh*max(n_hits_and_extreme$nExtreme),"gene_symbol"], 
                                 link_width = unit(1,"mm"),
                                 labels_gp = gpar(fontsize = 5), padding = unit(2, "mm")))

hms$male_model_01<- hms$male_model_01 +
  #Heatmap(n_hits_and_extreme$nExtreme, name="N extreme", width=unit(5,"mm")) +
  rowAnnotation(link = anno_mark(at = which(n_hits_and_extreme$nExtreme > thresh*max(n_hits_and_extreme$nExtreme)), 
                                 labels = n_hits_and_extreme[n_hits_and_extreme$nExtreme > thresh*max(n_hits_and_extreme$nExtreme),"gene_symbol"], 
                                 link_width = unit(2,"mm"),
                                 labels_gp = gpar(fontsize = 5), padding = unit(2, "mm")))

hms$female_model_01 <- hms$female_model_01 +
  #Heatmap(n_hits_and_extreme$nExtreme, name="N extreme", width=unit(5,"mm")) +
  rowAnnotation(link = anno_mark(at = which(n_hits_and_extreme$nExtreme > thresh*max(n_hits_and_extreme$nExtreme)), 
                                 labels = n_hits_and_extreme[n_hits_and_extreme$nExtreme > thresh*max(n_hits_and_extreme$nExtreme),"gene_symbol"], 
                                 link_width = unit(2,"mm"),
                                 labels_gp = gpar(fontsize = 5), padding = unit(2, "mm")))



write_rds(hms, snakemake@output[["rds"]])

pdf(snakemake@output[["pdf"]],onefile = T)
hms
dev.off()




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


