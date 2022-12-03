library(tidyverse)
library(topGO)
library(org.Dm.eg.db)
library(furrr)
library(progressr)

#mods_fl <- "upstream/final-models.collected-info.tsv.gz"
mods_fl <- snakemake@input[["mods"]]
mods <- read_tsv(mods_fl)

# this will be easier to write a single selection expression for if i make it a 
# signed -log10 pval
# Then I select the model (between male and female) with best pval
# then I create all posssible combos of TE and gene and fill in 0-scores to capture the full universe of genes
ag <- mods %>%
  mutate(score = if_else(!significant_model,1,adj_p.value_anova_x)) %>%
  mutate(score =replace_na(score,1)) %>%
  mutate(score = -log10(score) * sign(estimate.qnorm)) %>%
  group_by(feature.x) %>%
  slice_max(abs(score),with_ties = F) %>%
  ungroup() %>%
  dplyr::select(feature.x,feature.y,score) %>%
  right_join(.,tidyr::expand(.,feature.x,feature.y)) %>%
  mutate(score =replace_na(score,0))%>%
  split(.,.$feature.y) %>%
  map(dplyr::select,feature.x,score) %>%
  map(deframe)


allGO2genes <- c("BP","MF","CC") %>% 
  set_names(.,.) %>%
  map(~annFUN.org(whichOnto=.x, mapping="org.Dm.eg.db", ID="ensembl"))

nodes <- 100

# selection functions for both directions
selection<- list(pos = function(x){ return(x > -log10(0.1))},
                 neg = function(x){ return(x < log10(0.1))},
                 both = function(x){ return(abs(x) > -log10(0.1))})

plan(multisession, workers = 8)

get_enrichment <-function(x,ont,dir) {
  #set.seed(101)
  GOdata <- new("topGOdata",
                    ontology=ont,
                    allGenes=x,
                    annot=annFUN.GO2genes,
                    GO2genes=allGO2genes[[ont]],
                    geneSel=selection[[dir]],
                    nodeSize=10)
  result <- runTest(GOdata, algorithm = "weight01", statistic = "fisher") # to use a continuous score with ks test would need to think about scoreOrder arg
      
  allRes <-GenTable(GOdata, pval=result, topNodes = nodes)

  allRes %>% 
    as_tibble() %>%
    mutate(pval= as.numeric(pval)) %>%
    arrange(pval)
}

# in case of failure of a single te (if it has no significant correlated genes, for example)
possibly_get_enrichment <- possibly(get_enrichment,otherwise = NA)

df <- ag %>% 
  enframe() %>%
  expand_grid(.,ont=c("BP","CC","MF")) %>%
  expand_grid(.,dir = c("pos","neg","both"))

result <- df %>%
    mutate(data = future_pmap(list(value,ont,dir),function(x,y,z) possibly_get_enrichment(x,y,z),.options = furrr_options(seed=T)))

result <- result %>% dplyr::select(-value) %>%
  unnest(data)

write_tsv(result,snakemake@output[["tsv"]])

#res %>% filter(!GO.ID %in% c("GO:0008150","GO:0009987")) %>%
#  group_by(name) %>% slice_min(pval,n=1,with_ties = F) %>% print(n=Inf)
