library(tidyverse)
library(DESeq2)
library(Biostrings)
library(sva)
library(plyranges)

#se_path <- "upstream/kd.se.gene.0.rds"
se_path <- snakemake@input[["se"]]
gr_results_path <- snakemake@output[["grs"]]
dds_results_path <- snakemake@output[["dds"]]

gene_universe <- read_tsv("http://ftp.flybase.net/releases/FB2022_04/precomputed_files/genes/fbgn_fbtr_fbpp_expanded_fb_2022_04.tsv.gz", skip = 4)

allowed_genes <- gene_universe %>% filter(gene_type %in% c("protein_coding_gene"))

lkup <- gene_universe %>% dplyr::select(gene_ID, gene_symbol) %>% distinct()

se <- read_rds(se_path)

colData(se) <- colData(se) %>% 
  as_tibble(rownames = "name") %>%
  mutate(batch=case_when(str_detect(path_r1,"novo_")~"batch1",
                         str_detect(path_r1,"novo2")~"batch2",
                         str_detect(path_r1,"novo3")~"batch3")) %>%
  mutate(knockdown2 = paste(knockdown,tissue,driver,sep="_")) %>%
  mutate(driver_tissue = paste(driver,tissue,sep="_")) %>%
  column_to_rownames("name") %>%
  DataFrame()

se <- se [rownames(se) %in% allowed_genes$gene_ID | 
               !str_detect(rownames(se),"FBgn"),]

se <- se[,se$driver %in% c("aTub","tj","Mef2.R")]

# ----------- count filtering --------------------------------------------------
dds.pre <- DESeqDataSet(se,~ 1)

dds.pre <- estimateSizeFactors(dds.pre,type="poscounts")

keep <- rowSums(counts(dds.pre, normalized=TRUE) >= 10) >= 5 &
  rowSums(counts(dds.pre, normalized=TRUE) >= 1) >= 10

se <- se[keep,]

# annotation -------------------------------------------------------------------

elementMetadata(se) <- elementMetadata(se) %>%
  as_tibble() %>%
  left_join(lkup,by=c(gene_id="gene_ID")) %>%
  mutate(gene_symbol = if_else(is.na(gene_symbol),gene_id,gene_symbol)) %>%
  #left_join(gc_lkup %>% dplyr::select(gene_symbol,longest_tx_gc=gc,longest_tx_length=length), by="gene_symbol") %>%
  DataFrame()

# adjust for batch ---------- --------------------------------------------------
assay(se) <- assay(se) %>% round()
mat <-assay(se)

mode(mat) <- "integer"

adjusted <- ComBat_seq(counts = mat, 
                       batch = se$batch,
                       group = se$knockdown,
                       covar_mod = model.matrix( ~ driver_tissue , data=colData(se)),
                       shrink.disp = T,
                       shrink = T,gene.subset.n = 500)

mode(adjusted) <- "integer"

se2 <- se; assay(se2) <- adjusted

# differential expression on full set at once, raw and corrected----------------

ses <- list(raw = se, adjusted = se2)

run_deseq <- . %>%
  estimateSizeFactors(type="poscounts") %>%
  DESeq()

contrasts <- list(CG16779.head = c("knockdown2","CG16779_female_head_Mef2.R","control_female_head_Mef2.R"),
                  Unr = c("knockdown2","Unr_female_head_Mef2.R","control_female_head_Mef2.R"),
                  NfI = c("knockdown2","NFI_female_head_Mef2.R","control_female_head_Mef2.R"),
                  vvl = c("knockdown2","vvl_female_head_Mef2.R","control_female_head_Mef2.R"),
                  pan.testis = c("knockdown2","pan_male_gonad_aTub","control_male_gonad_aTub"),
                  pan.tj.ovary = c("knockdown2","pan_female_gonad_tj","control_female_gonad_tj"),
                  CG16779.testis = c("knockdown2","CG16779_male_gonad_aTub","control_male_gonad_aTub"),
                  CG16779.tj.ovary = c("knockdown2","CG16779_female_gonad_tj","control_female_gonad_tj"),
                  pan.head = c("knockdown2","pan_female_head_Mef2.R","control_female_head_Mef2.R"))

contrasts <- contrasts %>% set_names(.,map_chr(.,paste,collapse="_"))

get_res <- function(dds) {
  map(contrasts,~lfcShrink(dds,contrast = .x,format = "GRanges",type="normal")) %>%
    map(~{mutate(.x,feature=names(.))})
}

dds_list <- map(ses,~DESeqDataSet(.,~knockdown2)) %>%
  map(run_deseq)

res_list <- map(dds_list,get_res) 

# differential expression on split sets, no correction ----------------

dds.testis <- DESeqDataSet(se[,se$tissue == "male_gonad"],~ knockdown2)
dds.testis$knockdown2 <- relevel(dds.testis$knockdown2,ref="control_male_gonad_aTub")

dds.ovary_tj <- DESeqDataSet(se[,str_detect(se$tissue,"female_gonad") & se$driver == "tj"],~knockdown2)
dds.ovary_tj$knockdown2 <- relevel(dds.ovary_tj$knockdown2, ref = "control_female_gonad_tj")

dds.female_head <- DESeqDataSet(se[,str_detect(se$tissue,"female_head")],~knockdown2)
dds.female_head$knockdown2 <- relevel(dds.female_head$knockdown2, ref = "control_female_head_Mef2.R")

dds_list_indiv <- list(split.testis.raw=dds.testis,
                 split.ovary.tj.raw=dds.ovary_tj,
                 split.head.raw=dds.female_head)

dds_list_indiv <- map(dds_list_indiv, run_deseq)

res_list_indiv <- map(dds_list_indiv, ~{
  ds <- .x
  rn <- resultsNames(.x) %>% set_names(.,.) %>% .[2:length(.)] 
  map(rn,~lfcShrink(ds,coef = .x,type="normal",format="GRanges")) %>%
    map(~{mutate(.x,feature=names(.))})
})

#  combine all results ---------------------------------------------------------

res_list <- c(res_list,res_list_indiv)
dds_list <- c(dds_list,dds_list_indiv)

saveRDS(res_list,gr_results_path)
saveRDS(dds_list,dds_results_path)


# ---------------------------------------------------

# xx <- res_list %>% map(enframe) %>%
#   bind_rows(.id="approach") %>%
#   mutate(data = map(value,as_tibble)) %>%
#   dplyr::select(-value) %>%
#   unnest(data)
# 
# 
# xx %>%
#   filter(!str_detect(feature,"FBgn") & padj < 0.1) %>%
#   filter(str_detect(name,"CG16779|NFI")) %>%
#   #filter(str_detect(name,"pan")) %>%
#   ggplot(aes(log2FoldChange,-log10(pvalue))) +
#   geom_point() +
#   facet_wrap(~approach + name,scales = "free_y")
# 
# 
# xx %>%
#   filter(!str_detect(name,"Intercept")) %>%
#   unite(approach2,approach,name) %>%
#   dplyr::select(approach2,feature,log2FoldChange) %>%
#   dplyr::group_by(feature, approach2) %>%
#   dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
#   arrange(-n)
#   pivot_wider(names_from = "approach2", values_from = "log2FoldChange")
# 
# filter(xx,is.na(feature))


