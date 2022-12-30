library(tesseract)
library(tidyverse)
library(universalmotif)

eng <- tesseract("eng")
text <- tesseract::ocr("https://doi.org/10.1371/journal.pgen.1004591.s007", engine = eng)

hmg_motifs_df <- read_delim(text, skip = 8)

hmg_motifs_df <- hmg_motifs_df  %>%
  unite(dm3_location, c("Chr.","Location"),sep = ":")

hmg_motifs_df  <- hmg_motifs_df %>% dplyr::select(dm3_location, Sequence)

hmg_motifs_df <- hmg_motifs_df %>% mutate(nm = paste0("Archbold14::",Sequence))

motifs <- ifelse(exists("snakemake"),snakemake@input[["meme"]],
                 "results/analysis/motifs/combined_xstreme.meme") %>%
  read_meme()


names(motifs) <- motifs %>% map_chr( `@`, name)

degenerate_hmgs <- hmg_motifs_df %>%
  dplyr::select(nm,Sequence) %>%
  deframe() %>%
  map(universalmotif::create_motif)

degenerate_hmgs <- degenerate_hmgs %>% imap(~{.x@name <- .y; .x})

all_motifs <- c(degenerate_hmgs, motifs)

# default args, except for score.strat - the sampling distribution
# of pearson's rho is skewed, so averaging after FZT makes more sense
# I used the same approach to average the rho's for gene x gene correlation
# among salmon replicates
METHOD = "PCC"
SCORE.STRAT= "fzt"
compare_motifs2 <- function(m) {
  
  mat <- compare_motifs(motifs = m,
                                   method = METHOD,
                                   score.strat = SCORE.STRAT)
  
  p_df <- compare_motifs(motifs = m, 
                               compare.to = 1:length(all_motifs),
                               method = METHOD,
                               score.strat = SCORE.STRAT)
  
  p_df <- p_df[p_df$subject != p_df$target,]
  
  p_df$padj <- p.adjust(p_df$Pval, method="BH")
  
  return(list(p=p_df, sim = mat))
}


motif_comparison <- compare_motifs2(all_motifs)

motif_comparison$p <- motif_comparison$p %>%
  as_tibble() %>%
  #filter(str_detect(subject,"pan")|str_detect(target,"pan")) %>%
  mutate(motifs = map2(subject,target, ~{c(all_motifs[[.x]],all_motifs[[.y]])})) %>%
  mutate(gg = map(motifs, ~view_motifs(.x,method = METHOD, score.strat = SCORE.STRAT, text.size = 5)))

saveRDS(motif_comparison$p, snakemake@output[["motif_comparison"]])
saveRDS(motif_comparison$sim, snakemake@output[["motif_similarity"]])
write_tsv(hmg_motifs_df,snakemake@output[["archbold_motifs"]])


# # ------------------------------------------------------------------------------
motif_comparison$p %>% 
  as_tibble() %>%
  filter(padj < 0.1) %>%
   filter(!(str_detect(subject,"Arch") & str_detect(target,"Arch")))
# 
# # ------------------------------------------------------------------------------
# 
# 
# df %>%
#   filter(motifA!=motifB) %>%
#   filter(str_detect(motifA,"pan") & str_detect(motifB,"Arch|pan")) %>% 
#   arrange(-sim) %>%
#   mutate(motifBClass = str_extract(motifB,".+(?=::)")) %>%
#   mutate(motifA = fct_reorder(motifA,sim)) %>%
#   ggplot(aes(motifA,sim,color=motifBClass)) +
#   geom_point(position = ggplot2::position_jitterdodge()) +
#   theme(axis.text.x = element_text(angle=45,hjust=1))
# 
# # ------------------------------------------------------------------------------
# 
# compr <- 1-motif_comp_mat
# compr <- as.dist(compr)
# 
# labels <- attr(compr, "labels")
# 
# compr <- hclust(compr)
# 
# # see combine_xstreme_motifs.R for ideas about plotting as unrooted tree
# family <- sapply(all_motifs, function(x) str_extract(x["name"],".+(?=::)"))
# 
# pheatmap::pheatmap(motif_comp_mat,
#                    cluster_rows = compr, 
#                    cluster_cols = compr,
#                    scale = "none")
# 



