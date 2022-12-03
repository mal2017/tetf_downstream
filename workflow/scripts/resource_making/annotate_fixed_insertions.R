library(tidyverse)
library(plyranges)
library(rtracklayer)

lms <- ifelse(exists("snakemake"),snakemake@input[["lms"]],
              "upstream/final-models.collected-info.tsv.gz") %>%
  read_tsv() %>%
  filter(significant_x) %>%
  dplyr::rename(sex=model)

# -------- get TFs -------------------------------------------------------------

remap_simple <- ifelse(exists("snakemake"),snakemake@input[["remap"]],
                     "resources/remap2022_nr_macs2_dm6_v1_0.bed.gz") %>%
  import() %>%
  mutate(ChIP = str_extract(name,".+(?=:)")) %>%
  split(.,.$ChIP) %>%
  GenomicRanges::reduce()

remap_simple <- remap_simple[names(remap_simple) %in% lms$gene_symbol]

seqlevelsStyle(remap_simple) <- "NCBI"

# -------- chromatin -------------------------------------------------------------
#outside data downloaded from WUSTL GEP table browser aug 5, 2022 for dm6 (though most are lifted from r5).
#https://www.researchgate.net/figure/Chromatin-annotation-of-the-Drosophila-melanogaster-genome-a-A-9-state-model-of_fig6_279835262

#het_path <- "resources/het_domains_r5todm6.bed.gz"
#s2_chromstate_path <- "resources/chromstate_s2_9state_r5todm6.bed.gz"

het_path <- snakemake@input[["het"]]
s2_chromstate_path <- snakemake@input[["s2_chromstate"]]

het <- rtracklayer::import(het_path) %>%
  mutate(chromatin = str_extract(name,".+(?=_)"))

s2_chromstate <- rtracklayer::import(s2_chromstate_path) %>%
  mutate(s2_state = str_extract(name,".+(?=_)"))

seqlevelsStyle(het) <- "NCBI"
seqlevelsStyle(s2_chromstate) <- "NCBI"

# -------- get TEs -------------------------------------------------------------
# Which set of TE insertions to use is important to consider.. 
# One option is ISO-1 insertions via repeatmasker, as in the cell below (not run).

#insertions_path <- "results/resources/dgrp_tidal_insertions.unique.bb"
insertions_path <- snakemake@input[["insertions"]]

#all_ins_path <- "results/resources/dgrp_tidal_insertions.bb"
all_ins_path <- snakemake@input[["all_ins"]]

ref_ins <- import(all_ins_path) %>%
  filter(str_detect(name,"reference"))

ref_ins <- ref_ins %>% split(.,str_extract(.$name,".+(?=\\.DGRP)"))

tes_p1 <- import(insertions_path) %>%
  mutate(nOL = score)

strand(tes_p1) <- "*"

n_strains <- max(tes_p1$nOL)

# For INE-1, not all are fixed - most likely just small fragments aren't
fixed_in_dgrp <- tes_p1 %>% 
  filter(nOL == n_strains) %>% 
  split(.,.$name) %>% as.list()

# this only retains insertions that are 'fixed_in_dgrp' (ie nOL > x)
fixed_in_ref <- map(names(fixed_in_dgrp),.f=function(x){ 
  print(x)
  subsetByOverlaps(ref_ins[[x]],fixed_in_dgrp[[x]]) %>%
    mutate(name = x)
})

tes_p2 <-  plyranges::select(unlist(GRangesList(fixed_in_ref)), repeat_element = name) %>%
  split(.,.$repeat_element) %>%
  GenomicRanges::reduce() %>% 
  unlist() %>%
  mutate(.,repeat_element = names(.))

tes_p2 <- tes_p2[seqnames(tes_p2) %in% seqlevels(remap_simple)]

strand(tes_p2) <- "*" 

tes <- subsetByOverlaps(tes_p1,tes_p2) %>% filter(nOL==n_strains)


# annotate by feature ----------------------------------------------------------
library(GenomicFeatures)
txdb <- GenomicFeatures::makeTxDbFromGFF(ifelse(exists("snakemake"),snakemake@input[["gtf"]],
                                                "resources/dmel-all-r6.41.gtf.gz"))

ex <-exons(txdb) %>% GenomicRanges::reduce() %>% unstrand()

intr <- intronsByTranscript(txdb) %>% unlist() %>% GenomicRanges::reduce() %>% unstrand()

fiveUTR <- fiveUTRsByTranscript(txdb) %>% unlist() %>% GenomicRanges::reduce() %>% unstrand()

threeUTR <- threeUTRsByTranscript(txdb) %>% unlist() %>% GenomicRanges::reduce() %>% unstrand()

prom2kbup <- GenomicFeatures::promoters(txdb,upstream = 2000,downstream = 0) %>% GenomicRanges::reduce() %>% unstrand()

genes <- genes(txdb)

down5kb <- flank_downstream(genes,width = 5000)

tes <- tes %>% 
  mutate(.,
         exonic = unstrand(.) %over% ex,
         intronic = unstrand(.) %over% intr,
         utr5 = unstrand(.) %over% fiveUTR,
         utr3 = unstrand(.) %over% threeUTR,
         prom = unstrand(.) %over% prom2kbup,
         dn5kb = unstrand(.) %over% down5kb)

tes <- tes %>% 
  plyranges::add_nearest_distance(name = "nearest.tss" ,promoters(txdb,upstream = 0,downstream = 0))

# -------- get GC -------------------------------------------------------------
library(BSgenome.Dmelanogaster.UCSC.dm6)

seqlevelsStyle(Dmelanogaster) <- "NCBI" 

tes$GC <- letterFrequency(getSeq(Dmelanogaster, tes), "GC", as.prob=T)[,1]


# ------------------------------------------------------------------------------
# annotate insertion
tes2 <- tes %>%
  plyranges::mutate(.,ix = 1: length(.)) %>%
  plyranges::join_overlap_left(het) %>%
  #plyranges::join_overlap_left(bg3_chromstate) %>%
  plyranges::join_overlap_left(s2_chromstate) %>%
  plyranges::select(repeat_element = name.x,nOL,GC,chromatin,s2_state,nearest.tss,exonic,intronic,utr5,utr3,prom,dn5kb,ix) %>%
  plyranges::mutate(chromatin = replace_na(chromatin,"other"), s2_state = replace_na(s2_state,"other")) %>%
  plyranges::mutate(size = width) %>%
  plyranges::group_by(ix) %>%
  arrange(chromatin) %>% # if any overlap with hetchrom, make sure we retain that info
  group_by(ix) %>%
  plyranges::slice(1) %>% # retain only 1 locus in the case that multiple annotation regions overlap the same insertion
  ungroup() %>%
  GenomicRanges::sort()

# good spot to annotate with all TF overlaps
names(tes2) <- NULL

chip_overlap_mat <- lapply(remap_simple, FUN=function(x){countOverlaps(tes2,x)>0}) %>% 
  do.call(cbind,.,)

# wtf -  Gal here is mislabeled by Remap22 -  it is actually H3K27me3... 
# I traced it to this study: https://www.biorxiv.org/content/10.1101/127951v1.full.pdf
chip_overlap_mat <- chip_overlap_mat[,colnames(chip_overlap_mat)!="Gal"]

mcols(tes2) <- cbind(mcols(tes2),chip_overlap_mat)

write_rds(tes2,snakemake@output[["rds"]])
write_rds(remap_simple,snakemake@output[["remap"]])

