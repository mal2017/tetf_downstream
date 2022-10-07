library(tidyverse)
library(DESeq2)
library(Biostrings)
library(sva)
#theme_set(theme_prism())

se_path <- snakemake@input[["se"]]
gr_results_path <- snakemake@output[["grs"]]
dds_results_path <- snakemake@output[["dds"]]

gene_universe <- read_tsv("http://ftp.flybase.net/releases/FB2022_04/precomputed_files/genes/fbgn_fbtr_fbpp_expanded_fb_2022_04.tsv.gz", skip = 4)

allowed_genes <- gene_universe %>% filter(gene_type %in% c("protein_coding_gene"))

probable_contaminants <- gene_universe %>% filter(gene_type %in% c("antisense_lncRNA_gene",
                                                                   "pseudogene",
                                                                   "lncRNA_gene",
                                                                   "ncRNA_gene",
                                                                   "tRNA_gene","rRNA_gene","mt_LSU_rRNA_gene",
                                                                   "mt_LSU_rRNA_gene","RNase_MRP_RNA_gene","hpRNA_gene",
                                                                   "H_ACA_box_snoRNA_gene","C_D_box_snoRNA_gene","sbRNA_gene"))

lkup <- gene_universe %>% dplyr::select(gene_ID, gene_symbol) %>% distinct()


se.all <- read_rds(se_path)

colData(se.all) <- colData(se.all) %>% as_tibble(rownames = "name") %>%
  mutate(batch=ifelse(str_detect(path_r1,"novo2"),"batch2","batch1")) %>%
  mutate(knockdown2 = paste(knockdown,tissue,sep="_")) %>%
  mutate(broad = if_else(str_detect(knockdown,'control'),"ctl","kd")) %>%
  mutate(batch2 = paste(batch,tissue,sep="_")) %>%
  column_to_rownames("name") %>%
  DataFrame()

se.cont <- se.all[rownames(se.all) %in% probable_contaminants$gene_ID,]

se <- se.all[rownames(se.all) %in% allowed_genes$gene_ID | !str_detect(rownames(se.all),"FBgn"),]

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

# including knockdown here doesn't do a good job of removing the
# batch effect.
# I think this approach is more conservative, as some differences
# between KDs and controls will not be retained, but
# I'd rather be conservative and not discover a TE change
# than make a false discovery
mod.mat <- model.matrix( ~  tissue, data=colData(se))

adjusted <- ComBat_seq(counts = assay(se), 
                       batch = sapply(se$batch, switch, "batch1" = 1, "batch2" = 2, USE.NAMES = F),
                       #group = sapply(se$tissue, switch, "female_head" = 1, "female_gonad" = 2, USE.NAMES = F),
                       covar_mod = mod.mat,
                       shrink = F, #full_mod = T,
                       shrink.disp = T)

se2 <- se

assay(se2) <- adjusted


# differential expression ------------------------------------------------------

ses <- list(raw = se,
            corrected = se2)

run_deseq <- . %>%
  DESeqDataSet(~ knockdown2) %>%
  estimateSizeFactors(type="poscounts") %>%
  DESeq()

contrasts <- list(CG16779.head = c("knockdown2","CG16779_female_head","control_female_head"),
                  Unr = c("knockdown2","Unr_female_head","control_female_head"),
                  NfI = c("knockdown2","NFI_female_head","control_female_head"),
                  vvl = c("knockdown2","vvl_female_head","control_female_head"),
                  ct  = c("knockdown2","ct_female_gonad","control_female_gonad"),
                  mamo = c("knockdown2","mamo_female_gonad","control_female_gonad"),
                  awd = c("knockdown2","awd_female_gonad","control_female_gonad"),
                  CG16779.ovary = c("knockdown2","CG16779_female_gonad","control_female_gonad"))

get_res <- function(dds) {
  #map(~lfcShrink(dds,contrast = dds,type = "normal"))%>%
  map(contrasts,~results(dds,contrast = .x,format = "GRanges")) %>%
    map(~plyranges::mutate(.,feature = names(.)))
  #map_df(as_tibble,rownames="feature",.id="comparison") %>%
  #mutate(kd = str_extract(comparison,".+(?=\\.)|.+$"))
}

dds_list <- map(ses,run_deseq) 

res_list <- map(dds_list,get_res)

saveRDS(res_list,gr_results_path)
saveRDS(dds_list,dds_results_path)
