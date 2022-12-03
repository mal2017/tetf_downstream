library(tidyverse)
library(plyranges)
library(rtracklayer)
library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
library(GenomicFeatures)
library(DESeq2)
library(cocor)

# ------- Get params -----------------------------------------------------------
#tss_dist <- 5000 # 5000
tss_dist <- snakemake@params[["tss_dist"]]

#max_tes_per_gene <- 1 # 1
max_tes_per_gene <- snakemake@params[["max_tes_per_gene"]]

#coef_quantile_cutoff <- 0.95 #0.95
coef_quantile_cutoff <- snakemake@params[["coef_quantile_cutoff"]]

#min_pct_fixation <- 0.25 #0.25
min_pct_fixation <- snakemake@params[["min_pct_fixation"]]

#max_pct_fixation <- 0.75 #0.75
max_pct_fixation <- snakemake@params[["max_pct_fixation"]]

# ------- Get input files ------------------------------------------------------

#gtf_fl <- "resources/dmel-all-r6.41.gtf.gz"
gtf_fl <- snakemake@input[["gtf"]]
gtf <- GenomicFeatures::makeTxDbFromGFF(gtf_fl)


#mods_fl <- "upstream/final-models.collected-info.tsv.gz"
mods_fl <- snakemake@input[["mods"]]
lms <- read_tsv(mods_fl) %>% 
  filter(significant_x) #%>% 
  #filter(gene_symbol %in% c("CG16779","NfI","pan"))

#tfs_fl <- "resources/Drosophila_melanogaster_TF.txt"
tfs_fl <- snakemake@input[["tfs"]]
tfs <- read_tsv(tfs_fl)

#se_fl <- "upstream/se.gene.0.rds"
se_fl <- snakemake@input[["se"]]
se <- read_rds(se_fl)

# unique_insertions_fl <- "results/resources/dgrp_tidal_insertions.unique.bb" 
unique_insertions_fl <- snakemake@input[["unique_insertions"]]
ins0 <-import(unique_insertions_fl)

# all_ins_fl <- "results/resources/dgrp_tidal_insertions.bb"
all_ins_fl <-snakemake@input[["all_insertions"]]
all_ins<-import(all_ins_fl)

# ------- subset to get just the we'll be working with -------------------------
# mostly for memory purposes

lms2 <- lms %>% mutate(x_explained = sumsq_anova_x/total_variance) %>%
  filter(coef.quantile > coef_quantile_cutoff) %>%
  arrange(-x_explained) %>% 
  filter(feature.x %in% tfs$Ensembl) %>%
  group_by(model,gene_symbol) %>% slice_max(x_explained,n = max_tes_per_gene) %>% ungroup() %>%
  dplyr::select(feature.x,sex=model,gene_symbol,feature.y) %>%
  distinct()

ins <- ins0 %>% 
  filter(name %in% lms2$feature.y)

strains <- str_extract(all_ins$name,"DGRP_\\d+(?=\\.)") %>% unique()

tss <- genes(gtf) %>% resize(width=1,fix="start")

# get partially penetrant insertions
ins2 <- ins %>% filter(score < max_pct_fixation*max(score) & score > min_pct_fixation*max(score))

# get list of strains for each insertion
ins3 <- ins2 %>% 
  join_overlap_left(all_ins) %>% 
  filter(map2_lgl(name.y,name.x,str_detect)) %>%
  mutate(strain = str_extract(name.y,"DGRP_.+(?=\\.)")) %>%
  group_by(name.x) %>%
  reduce_ranges(strain = paste(strain,collapse=","))

# get genes whose 5' ends are withing `tss_dist` of a partially penetrant insertion
ins4 <- ins3 %>%
  plyranges::join_overlap_left(tss,maxgap = tss_dist) %>%
  filter(!is.na(gene_id)) %>%
  mutate(id = paste(name.x,as.character(.),sep=".")) %>%
  as_tibble() %>%
  dplyr::select(id,feature.y=name.x,strain,nearby_gene=gene_id) %>%
  left_join(lms2,.,by="feature.y") %>%
  mutate(strain = map(strain,~unlist(str_split(.x,","))))

# remove any instances of a "nearby_gene" that also happens to be correlated overall with
# the nearby TE insertion
ins4.5 <- ins4 %>% anti_join(lms2, by=c(nearby_gene = "feature.x",feature.y="feature.y"))

# transform and get subset of data we're actually going to use
dds <- DESeqDataSet(se,~1)
vds <- vst(dds)
vds <- vds[rownames(vds) %in% c(ins4$feature.x,ins4$nearby_gene,ins4$feature.y),se$Strain %in% strains]

# gets male or female expression data
get_mats <- function(strains,sex,gene1,gene2,ix) {
  message(ix)
  
  list(focus_strains = vds[c(gene1,gene2),vds$Strain %in% strains & vds$sex==sex],
       other_strains = vds[c(gene1,gene2),!vds$Strain %in% strains & vds$sex==sex]) %>%
    map(assay) %>%
    map(t) %>%
    enframe(name = "data_subset",value = "mat")
}

# get the actual expression data using above func
# We also 'slice_head' to not bother with multiple tests when the same TE is inserted
# near the same gene
ins5 <- ins4.5 %>%
  filter(nearby_gene %in% rownames(vds)) %>%
  group_by(feature.x,sex,gene_symbol,feature.y,nearby_gene) %>%
  slice_head(n=1) %>%
  ungroup() %>%
  mutate(ix = row_number()) %>%
  mutate(data = pmap(list(strain,sex,feature.x,nearby_gene,ix),get_mats)) %>%
  unnest(data) %>% 
  mutate(n = map_int(mat,nrow))

# calculate pearson corr + pvalue
ins6 <- ins5 %>%
  mutate(corr = map(mat,~broom::tidy(cor.test(.x[,1],.x[,2],method="pearson")))) %>%
  unnest(corr)

# convert the two mats of data to a single tibble to be more concise
ins7 <- ins6 %>%
  dplyr::select(id,sex,gene_symbol,feature.y,nearby_gene,data_subset,n,rho=estimate,p.value,mat) %>%
  mutate(across(contains("mat"),.fns = ~map(.x,as_tibble))) %>%
  pivot_wider(names_from = data_subset,values_from = c(rho,p.value,mat,n)) %>%
  mutate(data = map2(mat_focus_strains,mat_other_strains,~{bind_rows(focus_strains=.x,other_strains=.y,.id="data_subset")})) %>%
  dplyr::select(-contains("mat"))

# do a test that compares Fisher's z scores for each group of strains - 
# strains with an insertion and without
ins8 <- ins7 %>%
  drop_na() %>%
  mutate(cocor.test = pmap(list(rho_focus_strains,rho_other_strains,n_focus_strains,n_other_strains),
                           .f = function(w,x,y,z) as.htest(cocor.indep.groups(w,x,y,z)))) %>%
  mutate(fisherZ.tidy = map(cocor.test, ~{broom::tidy(.x$fisher1925)})) %>%
  unnest(fisherZ.tidy)

# adjust p values from Z test
ins9 <- ins8 %>% 
  mutate(padj = p.adjust(p.value,method="BH")) %>%
  arrange(-abs(rho_other_strains-rho_focus_strains))

write_rds(ins9,snakemake@output[["rds"]])
