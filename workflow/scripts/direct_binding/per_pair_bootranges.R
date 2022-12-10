library(tidyverse)
library(plyranges)
library(nullranges)
library(GenomicFeatures)
library(AnnotationDbi)
library(progress)


L_s <- ifelse(exists("snakemake"),as.numeric(snakemake@params[["L_s"]]),1e5)
nseg <- ifelse(exists("snakemake"),as.numeric(snakemake@params[["nseg"]]),3)
R <- ifelse(exists("snakemake"),as.numeric(snakemake@params[["R"]]),3)
blockLength <- ifelse(exists("snakemake"),as.numeric(snakemake@params[["blockLength"]]),5e5)

print(L_s)
print(nseg)
print(R)
print(blockLength)


lms <- ifelse(exists("snakemake"),snakemake@input[["lms"]],
                   "upstream/final-models.collected-info.tsv.gz") %>%
                   read_tsv() #%>% filter(significant_x)

chip <- ifelse(exists("snakemake"),snakemake@input[["remap"]],
               "results/resources/remap.gr.rds") %>%
  read_rds()

anno_ins_path <- ifelse(exists("snakemake"),snakemake@input[["anno_ins"]],
                        "results/resources/annotated_fixed_insertions.gr.rds")


txdb <- ifelse(exists("snakemake"),snakemake@input[["txdb"]],
               "results/resources/txdb") %>%
  loadDb()

# needed to compute the sequence lengths of the chromosomes
fa <- rtracklayer::import(ifelse(exists("snakemake"),snakemake@input[["genome_fa"]],"results/resources/genome.fasta"))

# prep the insertions for bootstrapping
ins <- read_rds(anno_ins_path)
seqlengths(ins) <- seqlengths(fa)[seqlevels(ins)]
ins <- keepStandardChromosomes(ins, pruning.mode = "coarse")
ins <- sortSeqlevels(ins)
ins <- sort(ins)
names(ins) <- NULL

# in keeping with convention of the bootranges vignette,
# i add an 'iter' field to the granges of interest with value zero
# so i can split the granges by this field later if i ever combine this with boots
ins$iter <- 0

# generate segmentation by gene density
g <- genes(txdb)

seqlengths(g) <- seqlengths(fa)[seqlevels(g)]
g <- keepStandardChromosomes(g, pruning.mode = "coarse")
g <- sortSeqlevels(g)
g <- sort(g)
names(g) <- NULL

seg_cbs <- segmentDensity(g, n = nseg, L_s = L_s, type = "cbs")

# add1 is just the first TF, so this is the first index in our slice
tfs <- colnames(mcols(ins))[which(colnames(mcols(ins)) == "ADD1"):ncol(mcols(ins))]

tes <- unique(ins$repeat_element)

set.seed(5) # for reproducibility

# remove seqlevels that are not in use
seqlevels(ins) <- seqlevelsInUse(ins)

# this function is used to get the bootranges for each TF/TE pair
# and allow failure if it can't find ways to assign to matched segments
possibly_bootranges <- possibly(function(...) {
  pb$tick()$print()
  bootRanges(...)
},otherwise = NULL)

# get TE/TF pairs to check
res <- lms %>%
  filter(gene_symbol %in% names(chip)) %>%
  #filter(gene_symbol == "pan") %>%
  dplyr::select(TF=gene_symbol,TE=feature.y) %>%
  distinct() #%>% head(1)

pb <- progress_estimated(length(unique(res$TE)))

# get bootstraps for each TE
res <- res %>%
  dplyr::select(TE) %>%
  distinct() %>%
  #head(1) %>%
  mutate(ins.gr = map(TE,~{ins[ins$repeat_element == .x,c("repeat_element","iter")]})) %>%
  mutate(boot.gr = map(ins.gr, possibly_bootranges, blockLength, seg=seg_cbs, R=R,withinChrom=F)) %>%
  left_join(res,., by="TE") %>%
  drop_na()

# i previously added an 'iter' field to the granges of interest
# now i make a func to split ranges using this field to create a grl
possibly_split <- possibly( function(.x){split(.x,.x$iter)},NULL)

# first for the boots
pb <- progress_bar$new(total = nrow(res))
res <- res %>% mutate(boot.gr= map(boot.gr,~{pb$tick();possibly_split(.x)}))

# now for the ins
pb <- progress_bar$new(total = nrow(res))
res <- res %>% mutate(ins.gr= map(ins.gr,~{pb$tick();possibly_split(.x)}))

# overlaps with the chip data are the statistic for the bootstrap test
# this function is used to get the statistic for each TF/TE pair
possibly_get_stat <- possibly(function(gr,TF) {
  unname(countOverlaps(gr,chip[[TF]]))
},NULL)

# get the observed statistics for each TF/TE pair
pb <- progress_bar$new(total = nrow(res))
res <- res %>% mutate(expected = map2(boot.gr,TF,~{pb$tick();possibly_get_stat(.x,.y)}))
pb <- progress_bar$new(total = nrow(res))
res <- res %>% mutate(observed = map2_int(ins.gr,TF,~{pb$tick();possibly_get_stat(.x,.y)}))


# number of expected values that are greater than or equal to the observed
res <- res %>%
  mutate(p = map2_dbl(observed,expected,~{sum(.y>=.x)/length(.y)}))

write_rds(res,snakemake@output[["rds"]])
write_rds(seg_cbs,snakemake@output[["seg_rds"]])