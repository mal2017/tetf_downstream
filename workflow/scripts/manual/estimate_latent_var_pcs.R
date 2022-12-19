library(tidyverse)
library(SummarizedExperiment)
library(rlang)
library(sva)
library(DESeq2)

set.seed(1)

sefile <- "upstream/se.gene.0.rds"

gene_universe <- read_tsv("http://ftp.flybase.net/releases/FB2022_04/precomputed_files/genes/fbgn_fbtr_fbpp_expanded_fb_2022_04.tsv.gz", skip = 4)

allowed_genes <- gene_universe %>% filter(gene_type %in% c("protein_coding_gene"))

se <- read_rds(sefile)

x <- assay(se,"counts")

# set up filter
filt <- parse_expr( "(((!str_detect(rownames(x),'FBgn')) & rowSums(x > 1) > 10) | (rowSums(x > 10) > 100)) & rowSums(x == 0) < 0.3*ncol(x)" )
features_2_use <-  rownames(x[eval(filt),])

x1 <- x[features_2_use,]

# ------------------------- scaling ----------------------------------

rv <- rowVars(x1)

select <- order(rv, decreasing=TRUE)[seq_len(min(250, length(rv)))]

x2 <- t(scale(t(log2(x1+1))))

# comes out to 4
num.sv(x2[select,], mod = model.matrix(~sex + Strain,data = colData(se)),method="be",B=100,seed = 2022)

#read_tsv("http://ftp.flybase.net/releases/FB2022_05/precomputed_files/genes/gene_groups_HGNC_fb_2022_05.tsv.gz")