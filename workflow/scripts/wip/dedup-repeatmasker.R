library(plyranges)
library(tidyverse)
library(rtracklayer)

rpm <- import("data/my_repeatmask/genome.fasta.out.gff") %>%
  mutate(family = str_extract(Target,"(?<=Motif:).+(?=\"[[:space:]])"))

strand(rpm) <- "*"

rpm <- mutate(rpm, name = paste0(as.character(rpm),"_",family))

# first get the regions that have no overlaps with non-self ranged
only_1 <- rpm[countOverlaps(rpm) == 1]
multiple_overlaps <- rpm[countOverlaps(rpm) > 1]


# now get those that need to be resolved. For now we pick overlap of 0.5 or higher
hits <- findOverlaps(multiple_overlaps,drop.self=T, drop.redundant=T)

x <- multiple_overlaps[queryHits(hits)]
y <- multiple_overlaps[subjectHits(hits)]

needs_resolution_ix <- which(width(pintersect(x,y)) / pmin(width(x),width(y)) >= 0.5)

no_resolution_needed <- hits[-needs_resolution_ix,]

resolution_needed <- hits[needs_resolution_ix,]

resolution_net <- findConnectedComps(resolution_needed)

clustering <- as(resolution_net,"DataFrame") %>% 
  as_tibble() %>% 
  mutate(cluster=row_number()) %>% 
  unnest(X) %>% 
  dplyr::select(cluster,node=X)

multiple_overlaps2 <- mutate(multiple_overlaps, node=1:length(multiple_overlaps)) %>%
  as_tibble() %>%
  left_join(clustering, by="node") %>%
  GRanges()


resolved <- multiple_overlaps2 %>%filter(!is.na(cluster)) %>%
  mutate(width2=width) %>%
  group_by(cluster) %>%
  arrange(cluster,-width2) %>%
  slice(1) %>%
  ungroup()


no_resolution_needed_gr <- multiple_overlaps2 %>%filter(is.na(cluster))

unique(c(subjectHits(no_resolution_needed),queryHits(no_resolution_needed))) %in% no_resolution_needed_gr$node


# These funcs are adapted from a post by Michael Lawrence: https://support.bioconductor.org/p/85092/
setAs("Hits", "graph", function(from) {
  stopifnot(queryLength(from) == subjectLength(from))
  NEL <- split(subjectHits(from), queryHits(from))
  graph::graphNEL(c(names(NEL),unlist(NEL)) %>% unique(),edgeL =  NEL,edgemode='directed')
})


findConnectedComps <- function(x) {
  g <- as(x, "graph")
  ManyToOneGrouping(unname(RBGL::connectedComp(g)))
}












# I originally used a version of this function for my rotation in the Herranz lab
# https://github.com/mal2017/hlr/blob/master/scripts/utils.R
# Note that the score returned by Repeatmasker should not be used for comparison,
# see the following link for details why.
# https://www.repeatmasker.org/webrepeatmaskerhelp.html
# IMPORTANT
# this is good for getting a "representative" single peak,
# not for merging overlaps
resolve_overlaps4 <- function(loci0,min_score = 0, score_col="norm_score") {
  loci0 <- loci0[mcols(loci0)[[score_col]] >= min_score]
  
  if (isDisjoint(loci0)) {return(loci0)}
  
  only_1 <- countOverlaps(loci0,loci0) == 1
  scores <- mcols(loci0)[[score_col]]
  loci <- loci0[!only_1]
  
  findOverlaps(loci,loci) %>% 
    as_tibble() %>% 
    group_by(queryHits) %>% 
    summarize(hits = list(subjectHits)) %>%
    mutate(top = map_int(hits,.f = function(x) x[which.max(scores[x])])) %>%
    .[["top"]] %>%
    unique() %>%
    loci[.] %>%
    c(.,loci0[only_1]) %>%
    resolve_overlaps4(min_score = min_score, score_col = score_col)
  
}




reso <- mutate(rpm,width2=width) %>% resolve_overlaps4(score_col = "width2")

export(reso,"~/Downloads/deduped_tes.bed")
