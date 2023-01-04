library(tidyverse)
library(plyranges)
library(GenomicRanges)
library(GenomicFeatures)
library(rtracklayer)
library(AnnotationDbi)
library(furrr)

# import saved txdb
threads <- ifelse(exists("snakemake"),
               snakemake@threads,4)

plan(multisession, workers = threads)

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

remap <- remap[c("pan","CG16779","NfI")]

remap <- remap %>%
  unlist() %>%
  mutate(.,ChIP = names(.)) 

# get first bp of all genes so I can set upstream distances on the fly
all_genes <- genes(txdb) %>% 
  anchor_5p() %>%
  mutate(width=0)

# remove genes way in the middle of hetchrom that don't have peaks near them - this makes the
# below test more realistic
all_genes <- all_genes %>% 
  add_nearest_distance(.,y=remap) %>%
  filter(distance < 1000)


pirna_genes <- all_genes %>% filter(gene_id %in% pirna_gene_ids$gene_ID)
other_genes <- all_genes %>% filter(!gene_id %in% pirna_gene_ids$gene_ID)

shared_seqs <- intersect(seqlevelsInUse(remap),seqlevelsInUse(pirna_genes))
seqlevels(pirna_genes, pruning.mode="coarse") <- shared_seqs
seqlevels(remap, pruning.mode="coarse") <- shared_seqs


get_cont_tbl <- function(d, fac, fix="end") {
  c(sum(countOverlaps(resize(pirna_genes, fix=fix,width=d), plyranges::filter(remap, ChIP==fac)) > 0),
    sum(countOverlaps( resize(other_genes, fix=fix,width=d), plyranges::filter(remap, ChIP==fac)) > 0),
    sum(countOverlaps(resize(pirna_genes, fix=fix,width=d), plyranges::filter(remap, ChIP==fac)) == 0),
    sum(countOverlaps(resize(other_genes, fix=fix,width=d), plyranges::filter(remap, ChIP==fac)) == 0)) %>%
    matrix(nrow=2,dimnames = list(c("piRNA","other"),c("bound","unbound")), byrow = F)
}

d <- c(500,1000, 2500, 5000, 10000, 25000)
# for all chips, see if they're enriched
res <- unique(remap$ChIP) %>%
  expand_grid(ChIP = ., 
            d = d) %>%
  mutate(contingency_table  = future_map2(d, ChIP, get_cont_tbl)) %>%
  mutate(fish.res = future_map(contingency_table, ~broom::tidy(fisher.test(.x)))) %>%
  unnest(fish.res)


res <- res %>% mutate(padj = p.adjust(p.value, method="BH"))

g <- all_genes %>% join_overlap_left_within(remap,maxgap=500) %>%
  as_tibble() %>%
  distinct() %>%
  group_by(gene_id) %>%
  summarise(n=sum(!is.na(ChIP))) %>%
  mutate(is.piRNA = gene_id %in% pirna_gene_ids$gene_ID) %>%
  ggplot(aes(is.piRNA,log10(n+1))) +
  geom_boxplot() +
  ggpubr::stat_compare_means()
  
ggsave(snakemake@output[["png"]],g)
saveRDS(g,snakemake@output[["rds"]])
saveRDS(res,snakemake@output[["res_rds"]])


# -------------------------------------------------------------------
# get list of piRNA pathway genes and the nearest peak to each
# -------------------------------------------------------------------

poi <- c("aub","piwi","rhi","del","AGO3","arx","tej","egg","vas")
pirna_genes %>%
  {x <- .; x@elementMetadata$distance <- NULL; x} %>%
  join_overlap_left(remap, maxgap=500) %>%
  #filter(distance == 0) %>%
  as_tibble() %>%
  mutate(bound= T) %>%
  distinct() %>%
  dplyr::select(gene_id, ChIP, bound) %>% pivot_wider(names_from = ChIP, values_from = bound) %>%
  mutate(across(where(is.logical), replace_na, F)) %>%
  dplyr::select(-`NA`) %>% 
  left_join(pirna_gene_ids, by=c(gene_id = "gene_ID")) %>%
  relocate(gene_symbol) %>% filter(!pan & (NfI | CG16779)) %>%
  print(n=Inf) %>%
  filter(gene_symbol %in% poi)

# ------------------------------------------------------------------------------
# approach Fishers 2
# the below is a similar analysis with slightly different coding convention
# using plyranges join_nearest. 
# ------------------------------------------------------------------------------

# gr <- remap %>%
#   split(.,.$ChIP) %>%
#   as.list() %>%
#   map(~join_nearest_left(all_genes, .x, suffix = c(".piRNA",".ChIP"), distance = T)) %>%
#   GRangesList() %>%
#   unlist()
# 
# 
# dat <- gr %>%
#   mutate(is.piRNA = gene_id %in% pirna_gene_ids$gene_ID) %>%
#   as_tibble()
# 
# 
# res <- dat %>%
#   count(ChIP,is.piRNA, distance < 500) %>%
#   pivot_wider(names_from = is.piRNA, values_from = n, names_glue = "is.piRNA_{is.piRNA}") %>%
#   mutate(across(contains("is.piRNA"), replace_na, 0)) %>%
#   arrange(ChIP, -`distance < 500`) %>%
#   relocate(is.piRNA_TRUE, .before = 'is.piRNA_FALSE') %>%
#   nest(-ChIP) %>%
#   mutate(data = map(data, column_to_rownames, "distance < 500")) %>%
#   mutate(fisher.res = map(data, ~broom::tidy(fisher.test(.x)))) %>%
#   unnest(fisher.res)
# 
# 
# res %>% mutate(padj = p.adjust(p.value, method="BH")) %>%
#   filter(ChIP== "pan") %>% pull(data)
#   filter(padj < 0.1 & estimate > 1) %>%
#   arrange(-estimate)
