library(tidyverse)
library(plyranges)
library(GenomicRanges)
library(GenomicFeatures)
library(rtracklayer)
library(AnnotationDbi)

# import saved txdb
txdb <- ifelse(exists("snakemake"),
    snakemake@input$txdb,
    "results/resources/txdb") %>%
    loadDb()

# import pirna gene ids
pirna_gene_ids <- ifelse(exists("snakemake"),
    snakemake@input$pirna,
    "results/resources/pirna_pathway.tsv") %>%
    read_tsv()

# import remap peaks as gr
remap <- ifelse(exists("snakemake"),
    snakemake@input$remap,
    "results/resources/remap.gr.rds") %>%
    readRDS()

#remap <- remap[c("pan","CG16779","NfI")]

remap <- remap %>%
  unlist() %>%
  mutate(.,ChIP = names(.)) 

all_genes <- genes(txdb) %>% 
  anchor_5p() %>%
  mutate(width=0)

pirna_genes <- all_genes %>% filter(gene_id %in% pirna_gene_ids$gene_ID)
other_genes <- all_genes %>% filter(!gene_id %in% pirna_gene_ids$gene_ID)

shared_seqs <- intersect(seqlevelsInUse(remap),seqlevelsInUse(pirna_genes))
seqlevels(pirna_genes, pruning.mode="coarse") <- shared_seqs
seqlevels(remap, pruning.mode="coarse") <- shared_seqs

gr <- remap %>% 
  split(.,.$ChIP) %>%
  as.list() %>%
  map(~join_nearest(pirna_genes, .x, suffix = c(".piRNA",".ChIP"), distance = T)) %>%
  GRangesList() %>%
  unlist()

  
dat <- gr %>%
  as_tibble() %>%
  mutate(ChIP2  = ifelse(ChIP %in% c("NfI","CG16779","pan"),ChIP,"other"))
  
bg <- dat$distance
# ------------------------------------------------------------------------------
# approach 1 - wilcox.test on distances
# ------------------------------------------------------------------------------
dat.test <- dat %>%
  nest(-ChIP,-ChIP2) %>%
  mutate(med.dist = map_dbl(data,~median(.x$distance)), med.sbg = median(bg)) %>%
  mutate(wilcox.tbl = map(data,~broom::tidy(wilcox.test(.x$distance,bg,alternative = "two.sided")))) %>%
  dplyr::select(-data) %>%
  unnest(wilcox.tbl) %>%
  mutate(padj = p.adjust(p.value,method="BH"))

dat  %>%
  ggplot(aes(fct_reorder(ChIP2, distance, .fun=mean),log10(distance+1))) +
  #geom_boxplot(outlier.shape = NA) +
  geom_violin(scale = "width") +
  geom_jitter(width = 0.25, size=0.3) +
  geom_text(data=filter(dat.test, ChIP2 %in% c("NfI","CG16779","pan")), aes(x=ChIP2,y=5, label=paste("padj=",format.pval(padj,digits = 1))))

# ------------------------------------------------------------------------------
# approach 2 - regioner, probably not appropriate.
# this basically says every TF significantly overlaps, probably because TFs
# bind promoters more than random spots in the genome
# ------------------------------------------------------------------------------
library(regioneR)
fa <- import(ifelse(exists("snakemake"),snakemake@input[["genome_fa"]],"resources/dmel-all-chromosome-r6.41.fasta.gz"))
names(fa) <- names(fa)%>% str_extract(".+(?=\\s+type)")

genome_df <- fa[seqlevels(remap),] %>% seqlengths() %>% 
  enframe(name = "seqnames",value = "end") %>%
  mutate(start=1) %>%
  dplyr::select(seqnames,start,end) %>%
  as.data.frame()

remap2 <- remap %>% split(.,.$ChIP)

res2tibble <- function(x) {
  tibble(pval = x$pval, 
         z = x$zscore,
         alternative = x$alternative,
         observed = x$observed)
}

res <- tibble(TF=unique(names(remap))) %>%
  filter(TF %in% c("pan","CG16779","NfI")) %>%
  mutate(result = map(TF,.f=function(TF) {
    set.seed(123)
    z <- overlapPermTest(A= remap2[[TF]], B = resize(pirna_genes, fix = "end",2000),verbose=T,
                         genome=genome_df, 
                         alternative = "greater",
                         non.overlapping =T,
                         per.chromosome =T,
                         mc.set.seed=FALSE,
                         nperm=10)
    
    res2tibble(z$numOverlaps)
  }))

# --------
c(sum(countOverlaps(resize(pirna_genes, fix="end",width=2000), remap2$pan) > 0),
  sum(countOverlaps( resize(other_genes, fix="end",width=2000), remap2$pan) > 0),
  sum(!countOverlaps(resize(pirna_genes, fix="end",width=2000), remap2$pan) > 0),
  sum(!countOverlaps(resize(other_genes, fix="end",width=2000), remap2$pan) > 0)) %>%
  matrix(nrow=2) %>%
  fisher.test()

pirna_genes

# -----