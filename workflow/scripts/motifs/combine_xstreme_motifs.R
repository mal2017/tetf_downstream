library(universalmotif)
library(magrittr)
library(tidyverse)

# function renames individual slot in motif obj
rename_name <-  function(nm, TF, sep="::") {
  paste(TF,nm,sep=sep)
}

# function renames the name and altname slots of
# a motif obj
rename_motif <- function(mot,TF) {
  mot@altname %<>% rename_name(TF)
  mot@name <- (mot@altname) %>% str_split("-") %>% unlist %>% .[c(1,2)] %>% paste(collapse = "-")
  return(mot)
}

# renames a list of motif objects
rename_motifs <- function(mot_l, TF) {
  map(mot_l, rename_motif, TF)
}

#memes <- Sys.glob("results/analysis/motifs/xstreme_per_tf/*/combined.meme")

memes <- snakemake@input[["memes"]] %>% paste0("/combined.meme")

memes <- memes %>%
  set_names(.,str_extract(.,"(?<=tf\\/).+(?=\\/combined)"))

memes <- map(memes, read_meme)

# do the renamimg
memes <- memes %>% imap(rename_motifs)

# unlist
memes <- unlist(memes)

write_meme(memes,snakemake@output[["meme"]])



# 
# memes.filt <- memes[str_detect(names(memes),"Myc|Max|Hand|twi|E2f1|E2F2")]
# 
# 
# compr <- compare_motifs(memes.filt)
# compr <- 1-compr
# compr <- as.dist(compr)
# 
# labels <- attr(compr, "labels")
# 
# compr <- ape::as.phylo(hclust(compr))
# 
# family <- sapply(memes.filt, function(x) str_extract(x["name"],".+(?=::)"))
# family.unique <- unique(family)
# 
# family.annotations <- list()
# for (i in seq_along(family.unique)) {
#   family.annotations <- c(family.annotations,
#                           list(labels[family %in% family.unique[i]]))
# }
# names(family.annotations) <- family.unique
# 
# 
# compr <- ggtree::groupOTU(compr, family.annotations)
# 
# library(ggtree)
# tree <- ggtree(compr, layout = "circular") +
#   theme(legend.position = "bottom", legend.title = element_blank()) +
#   geom_tiplab()
# 
# tree





