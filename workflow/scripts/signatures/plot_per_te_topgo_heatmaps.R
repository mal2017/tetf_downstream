library(tidyverse)
library(ComplexHeatmap)
library(GOSemSim)
library(furrr)
library(progressr)

threads <- 8
threads <- snakemake@threads
plan(multisession, workers = threads)
#topgo_fl <- "results/analysis/signatures/per_te_topgo.tsv.gz"
topgo_fl <- snakemake@input[["tsv"]]

df <- read_tsv(topgo_fl) %>% filter(Annotated < 250) %>%
  filter(pval < 0.05)

# --------- semantic similarity clustering -------------------------------------
# cluster the TEs by the semantic similarity of their enriched go terms
dmGO <- c("BP","MF","CC") %>%
  set_names(.,.) %>%
  map(~godata('org.Dm.eg.db', ont=.x))

sem_sim_df <- df %>%
  group_by(name,ont,dir) %>% 
  slice_min(pval,n=5) %>%
  slice_max((Significant+1)/(Expected + 1),n=5) %>%
  ungroup() %>%
  dplyr::select(name,ont,GO.ID,dir) %>% 
  group_by(name,ont,dir) %>% 
  summarise(terms = list(GO.ID),.groups = "drop") %>% 
  full_join(.,.,by=c("ont","dir"),suffix = c(".a",".b"))

tictoc::tic()
with_progress(
  {
    p <- progressor(steps=nrow(sem_sim_df))
    
    result <- sem_sim_df %>% 
      mutate(semsim = future_pmap_dbl(list(terms.a,terms.b,ont),
                                      .f= function(x,y,z){
                                        p()
                                        mgoSim(x,y, semData=dmGO[[z]],measure="Wang",combine="BMA")
                                        })) 
  }
)
tictoc::toc()


hcs.te <- result %>%
  tidyr::unite("grp",ont,dir) %>%
  dplyr::select(grp,name.a,name.b,semsim) %>%
  split(.,.$grp) %>%
  map(dplyr::select,-grp) %>%
  map(~pivot_wider(.x,names_from = name.b,values_from = semsim)) %>%
  map(~column_to_rownames(.x,"name.a")) %>%
  map(~{hclust(as.dist(1-.),method="ward.D")})


# ------------- term clustering -----------------------------------------------
# cluster the terms by their enrichment in each TEs gene sets
df2 <- df %>%
  unite("grp",ont,dir) %>%
  split(.,.$grp) %>%
  map(~{
      dplyr::select(.x, name, GO.ID, pval) %>%
      distinct() %>%
      filter(pval < 0.05) %>%
      pivot_wider(names_from = name, values_from = pval) %>%
      mutate(across(where(is.numeric),replace_na,1)) %>%
      column_to_rownames("GO.ID") %>%
      as.matrix()
  })

bin <- df2 %>% map(~{.x < 0.05})
# thanks gptChat
bin <- lapply(bin, function(x) matrix(as.integer(x), nrow = nrow(x), ncol = ncol(x),dimnames = dimnames(x)))

hcs.term <- bin %>% map(dist,method="binary") %>% map(hclust, method="ward.D")

col_fun = structure(c("white","black"), names = c("0","1"))

hms <- names(bin) %>% set_names(.,.) %>% 
  imap(~{
    Heatmap(bin[[.x]],
                   cluster_columns  = hcs.te[[.x]],
                   cluster_rows = hcs.term[[.x]],
                   border = "black",
                   col=col_fun,
                   column_names_gp = grid::gpar(fontsize = 5),
                   show_row_names = F,
                   name = .y,
                   use_raster = T,
                   show_column_names = T)
  })


write_rds(hms,snakemake@output[["rds"]])


pdf(snakemake@output[["pdf"]],onefile = T,width = 7.5)
hms
dev.off()
