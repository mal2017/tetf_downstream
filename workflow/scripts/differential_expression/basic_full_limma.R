library(edgeR)
library(tidyverse)
library(SummarizedExperiment)

mods <- ifelse(exists("snakemake"),
               snakemake@input[["mods"]],
               "upstream/final-models.collected-info.tsv.gz") %>%
  read_tsv() %>%
  filter(significant_x)

#se_path <- "upstream/tfrnai.se.gene.0.rds"
se_path <- snakemake@input[["se"]]

#runselector_path <- "resources/full_tfrnai_srarunselector.txt"
runselector_path <- snakemake@input[["runselect"]]

#batch_path <- "resources/batch_data.tsv.gz"
batch_path <- snakemake@input[["batch"]]
tsv_results_path <- snakemake@output[["tsv"]]

gene_universe <- read_tsv("http://ftp.flybase.net/releases/FB2022_04/precomputed_files/genes/fbgn_fbtr_fbpp_expanded_fb_2022_04.tsv.gz", skip = 4)

allowed_genes <- gene_universe %>% filter(gene_type %in% c("protein_coding_gene"))

probable_contaminants <- gene_universe %>% filter(gene_type %in% c("tRNA_gene","rRNA_gene","mt_LSU_rRNA_gene","mt_LSU_rRNA_gene","RNase_MRP_RNA_gene","hpRNA_gene","H_ACA_box_snoRNA_gene","C_D_box_snoRNA_gene","sbRNA_gene", "antisense_lncRNA_gene",
                                                                   "pseudogene",
                                                                   "lncRNA_gene",
                                                                   "ncRNA_gene"))

lkup <- gene_universe %>% dplyr::select(gene_ID, gene_symbol) %>% distinct()

oi <- gene_universe %>% filter(gene_symbol %in% c("pan","ct","awd","mamo","Unr","NfI","vvl","CG16779","Dref")) %>%
  dplyr::select(gene_ID,gene_symbol) %>%
  distinct()

se <- read_rds(se_path)

#se <- map(c(oi$gene_symbol,"LacZ"),~str_detect(colnames(se),.)) %>%
#  map(as.matrix) %>%
#  do.call(cbind,.) %>%
#  rowAnys() %>%
#  se[,.]

cd <- colData(se) %>% as_tibble() %>% 
  mutate(Run = str_extract(sample_name,"SRR.+"))

runselector <- read_csv(runselector_path) %>% 
  dplyr::select(Run,GSM=`Sample Name`,RNAi_target_gene_name)

batch_info <- read_tsv(batch_path) %>% 
  dplyr::select(GSM,plate,well,gene_ID) %>%
  mutate(plate2 = str_extract(plate,".+(?=-)")) %>%
  mutate(plate3 = str_replace(plate,"-","\\."))

cd2 <- left_join(cd,runselector, by="Run") %>%
  left_join(batch_info,by="GSM") %>%
  drop_na() %>%
  column_to_rownames("names") %>%
  DataFrame()

stopifnot(rownames(cd2) == colnames(se))

colData(se) <- cd2

colData(se)$rnai <- colData(se)$rnai %>% as.factor() %>% relevel(ref = "LacZ")

se <- se[rownames(se) %in% allowed_genes$gene_ID | !str_detect(rownames(se),"FBgn") & !str_detect(rownames(se),"ERCC"),]

# REMOVE PER MDS PLOT
se <- se[,colnames(se)!="LacZ_SRR3488228"]

# we're interested in the knockdowns for TE correlated genes
se <- se[,se$rnai %in% mods$gene_symbol]

d0 <- DGEList(assay(se))

design <- model.matrix(~ plate + rnai , data = colData(se))

d0 <- calcNormFactors(d0,)
keep <- filterByExpr(d0, design,min.count=50, min.prop=0.75,)
d0 <- d0[keep, ]

v <- voom(d0,design,normalize.method = "quantile")

fit <- lmFit(v, design)

coefs <- unique(colData(se)[se$rnai!="LacZ",]$rnai) %>% 
  set_names(.,.) %>%
  map_chr(~sprintf("rnai%s",.))

eb <- eBayes(fit,trend = T,robust = T)

res <- map_df(coefs, 
              ~as_tibble(topTable(eb,number = Inf, coef = .,sort.by = "AveExpr"),rownames="feature"),
              .id="comparison")

write_tsv(res,tsv_results_path)