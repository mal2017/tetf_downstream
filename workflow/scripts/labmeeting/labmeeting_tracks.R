library(plotgardener)
library(plotgardenerData)
library(rtracklayer)
library(GenomicRanges)
library(tidyverse)
library(plyranges)
library(GenomicFeatures)
library(BSgenome.Dmelanogaster.UCSC.dm6)

coex_w_dp <- read_tsv("results/resources/filtered_models.tsv.gz") %>%
  filter(gene_symbol == "Dp") %>%
  pull(feature.y) %>%
  unique()

dp_bw_path <- "~/work/221002_tetf_smartmap_02/results/macs2/Dp_embryo_ChIPseq_IP_Rep1_ENCSR120ICM/bdgcmp/Dp_embryo_ChIPseq_IP_Rep1_ENCSR120ICM_qpois.bw"

loci <- GRanges(c("4:1286935-1290988","2R:25259810-25264157","2R:25251704-25254715","3L:19859997-19862622"))

dp_bw_sig <- rtracklayer::import(dp_bw_path,which = loci)


tes <- import("data/my_repeatmask/genome.fasta.out.gff") %>%
  mutate(repeat_element = str_extract(Target,'(?<=Motif:).+(?=\\")')) %>%
  mutate(Name = make.unique(repeat_element)) %>%
  mutate(ID = paste0("mRNA:",Name),Parent=ID,type="mRNA") %>%
  filter(repeat_element %in% coex_w_dp)

seqlevelsStyle(tes) <- "UCSC"
seqlevelsStyle(dp_bw_sig) <- "UCSC"

params_track1 <- pgParams(
  assembly = "dm6",
  x = 0, just = c("left", "top"),
  fill = "red", linecolor = "red",scale=T,
  height=1,
  fontsize=12,
  width = 9.5, length =9.5, default.units = "inches"
)

params_feature1 <- pgParams(
  assembly = "dm6",
  fontsize=12,
  x = 0, just = c("left", "top"),
  fill = "black", linecolor = "black",scale=T,
  height=0.15,width=9.5,
  length = 9.5, default.units = "inches"
)


# 6.3, 11.1
png("~/Downloads/labmeeting_figs/tracks.png",width =9.5, height = 4.6,units="in", res=150)
pageCreate(showGuides = F,width = 9.5, height = 4.6)


# example 1 --------------------------------------------------------------------

ex1 <- plotSignal(
  chrom = "chr4",chromstart = 1286935, chromend = 1290988,
  data = dp_bw_sig,
  y = 0,
  params = params_track1
)

plotText(
  label = "qpois",just = c("center","top"),
  fontsize = 11, rot = 90, x=0, y=0.5)


plotText(fontcolor="red",fontface="bold",
  label = "Dp ChIP signal (ENCSR120ICM)",
  just = c("right","top"),
  fontsize = 12, x=9.25, y=0.125)


plotGenomeLabel(
  chrom = "chr4",chromstart = 1286935, chromend = 1290988,
  y = 1,scale = "Kb",
  params=params_feature1,
)

set.seed(22)
plotRanges(
  tes,limitLabel = T,
  chrom = "chr4",chromstart = 1286935, chromend = 1290988,
  params = params_feature1,
  y=1.2,
  fill="black"
)

# Tahre

plotText(
  label = "TAHRE",just = c("center","top"),
  fontsize = 12, x=6, y=1.1)

# example 2 --------------------------------------------------------------------


ex2 <- plotSignal(
  chrom = "chr2R",chromstart = 25251704, chromend = 25254715,
  data = dp_bw_sig,
  y = 1.5,
  params = params_track1
)

plotText(
  label = "qpois",just = c("center","top"),
  fontsize = 12, rot = 90, x=0, y=2)

plotText(fontcolor="red",fontface="bold",
  label = "Dp ChIP signal (ENCSR120ICM)",
  just = c("right","top"),
  fontsize = 12, x=9.25, y=1.625)

plotGenomeLabel(
  chrom = "chr2R",chromstart = 25251704, chromend = 25254715,
  y = 2.5,scale = "Kb",
  params=params_feature1,
)

plotRanges(
  tes,
  chrom = "chr2R",chromstart = 25251704, chromend = 25254715,
  params = params_feature1,
  y=2.7,
  fill="black"
)

plotGenes(
  params = params_feature1,
  y=2.9,
  chrom = "chr2R",chromstart = 25251704, chromend = 25254715,
)


plotText(
  label = "INE-1",just = c("center","top"),
  fontsize = 12, x=3.75, y=2.7)
# ine 1


# example 3 --------------------------------------------------------------------


ex3 <- plotSignal(
  chrom = "chr3L", chromstart = 19859997, chromend = 19862622,
  data = dp_bw_sig,
  y = 3,
  params = params_track1
)

plotText(
  label = "qpois",just = c("center","top"),
  fontsize = 12, rot = 90, x=0, y=3.5)

plotText(fontcolor="red",fontface="bold",
  label = "Dp ChIP signal (ENCSR120ICM)",
  just = c("right","top"),
  fontsize = 12, x=9.25, y=3.125)

plotGenomeLabel(
  chrom = "chr3L", chromstart = 19859997, chromend = 19862622,
  y = 4,scale = "Kb",
  params=params_feature1,
)

plotRanges(
  tes,
  chrom = "chr3L", chromstart = 19859997, chromend = 19862622,
  params = params_feature1,
  y=4.2,
  fill="black"
)

plotGenes(
  params = params_feature1,
  y=4.4,
  chrom = "chr3L", chromstart = 19859997, chromend = 19862622,
)
# Rt1b
plotText(
  label = "Rt1b",just = c("center","top"),
  fontsize = 12, x=5.2, y=4.25)



dev.off()

# see plot-within-grp-tf-nullranges.R
#chr4:1,286,935-1,290,988 Tahre
# chr2R:25,259,810-25,264,157 Tahre/het
# chr3L:27,857,995-27,862,810 INE-1
# chr4:40,948-44,551 INE-1
# chr4:661,822-663,324 not a repmap peak, but looks pretty good
# chr4:669,210-670,770 INE-1
# also_related_tes %>%
#   filter(estimate > 1 & padj.overlap <0.1) %>%
#   filter(TF == "Dp") %>%
#   pull(insertions) %>%
#   .[[1]] %>%
#   filter(is.coex & n_olap !=0) %>%
#   as_tibble() %>%
#   print(n=Inf)
#   #rtracklayer::export("~/Downloads/Dp_overlaps_TEs.gff")
#   .$repeat_element %>% table()