library(plyranges)
library(tidyverse)
library(rtracklayer)

rpm <- import("data/my_repeatmask/genome.fasta.out.gff") %>%
  mutate(family = str_extract(Target,"(?<=Motif:).+(?=\"[[:space:]])"))

strand(rpm) <- "*"

rpm <- mutate(rpm, name = paste0(as.character(rpm),"_",family))

rpm <- mutate(rpm, width2 = width)

# first get the regions that have no overlaps with non-self ranged
only_1 <- rpm[countOverlaps(rpm) == 1]
multiple_overlaps <- rpm[countOverlaps(rpm) > 1] %>% 
  sort() %>%
  mutate(.,ix=1:length(.))


# now get those that need to be resolved. For now we pick overlap of 0.5 or higher
hits <- findOverlaps(multiple_overlaps,drop.self=T, drop.redundant=T)

x <- multiple_overlaps[queryHits(hits)]
y <- multiple_overlaps[subjectHits(hits)]

pct_olap <- width(pintersect(x,y)) / pmin(width(x),width(y))

needs_resolution_hits_ix <- which(pct_olap > 0.5)

resolution_needed_gr_ix <- hits[needs_resolution_hits_ix,] %>% {unique(c(subjectHits(.),queryHits(.)))}

resolution_needed_gr <- multiple_overlaps[resolution_needed_gr_ix,] %>% sort()
no_resolution_needed_gr <- multiple_overlaps[-resolution_needed_gr_ix,] %>% sort()

# sanity check: make sure we've not duplicated or removed anything
stopifnot(length(resolution_needed_gr) + length(no_resolution_needed_gr) == length(multiple_overlaps))

resolved <- resolution_needed_gr %>%
  sort() %>%
  reduce_ranges() %>%
  mutate(id=as.character(.)) %>%
  plyranges::join_overlap_left(multiple_overlaps[queryHits(resolution_needed)]) %>%
  group_by(id) %>%
  plyranges::slice(which.max(width2)[1]) %>%
  ungroup()

final <- c(only_1 %>% mutate(overlap_class = "no.overlap"),
  no_resolution_needed_gr %>% mutate(overlap_class = "below.overlap.cutoff"),
  resolved %>% mutate(overlap_class = "overlap.resolved"))


export(final,"~/Downloads/final.bed")

