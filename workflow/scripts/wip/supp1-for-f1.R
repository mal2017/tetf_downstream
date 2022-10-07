# Supp 1 (related to figure 1) -------------------------------------------------

if (!interactive()) pdf("~/Downloads/supp1.pdf",width = 6, height = 6)
pageCreate(width = 6, height = 6, default.units = "inches")



s1c <- plotGG(intro_hists$TEs + scale_fill_grey() + theme(legend.position = c(0.7,0.8)), x = 3.8, y=0.25, just = c("left","top"),width = 2.75*0.75, height=1.8*0.75)

f1d <- plotGG(intro_hists$genes + scale_fill_grey() + theme(legend.position = c(0.7,0.8)), x = 3.8, y=1.75, just = c("left","top"),width = 2.75*0.75, height=1.8*0.75)



# figure 2
read_rds("results/plots/plot_oi_s2rplus_te_gsea.rds")
read_rds("results/plots/plot_pirna_genes_in_lms.easy_bar.rds")

read_rds("results/plots/plot_pirna_genes_in_lms.volc.rds")

read_rds("results/plots/plot_pirna_genes_in_lms.mean_coef_ecdf.rds")

read_rds("results/plots/plot_pirna_genes_in_lms.ntotal_ecdf.rds")


read_rds("results/analysis/signatures/plot_gene_group_gsea.top10.rds")
read_rds("results/analysis/signatures/plot_gene_group_gsea.random_walk.rds")

read_rds("results/")
