library(tidyverse)

gene_universe <- read_tsv("http://ftp.flybase.net/releases/FB2022_04/precomputed_files/genes/fbgn_fbtr_fbpp_expanded_fb_2022_04.tsv.gz", skip = 4)

allowed_genes <- gene_universe %>% filter(gene_type %in% c("protein_coding_gene"))

probable_contaminants <- gene_universe %>% filter(gene_type %in% c("tRNA_gene","rRNA_gene","mt_LSU_rRNA_gene","mt_LSU_rRNA_gene","RNase_MRP_RNA_gene","hpRNA_gene","H_ACA_box_snoRNA_gene","C_D_box_snoRNA_gene","sbRNA_gene", "antisense_lncRNA_gene",
                                                                   "pseudogene",
                                                                   "lncRNA_gene",
                                                                   "ncRNA_gene"))

lkup <- gene_universe %>% dplyr::select(gene_ID, gene_symbol) %>% distinct()


#lm_paths <- Sys.glob("data/linear_models/*/*/lm.tidy.corrected.tsv.gz") 
lm_paths <- snakemake@input

# read in all models
all_models <-  lm_paths %>% 
  set_names(.,str_extract(.,"(?<=linear_models\\/).+/\\d")) %>%
  map_df(read_tsv,.id="model") %>%
  separate(model,into=c("model","rep"),sep = "/")

# won't work if using TLS
all_models <- all_models %>% group_by(model,rep) %>%
  mutate(padj = p.adjust(p.value,method="BH")) %>%
  ungroup() %>%
  mutate(significant = padj < 0.1)

# collapse replicates, but don't filter by reproducible yet
merged_models <- all_models %>%
  group_by(model,feature.x,feature.y) %>%
  #slice_max(abs(estimate.qnorm),n = 1,with_ties = F) %>%
  summarise(#min_lower.bound = min(lower.bound), 
            #max_upper.bound=max(upper.bound),
            mean_estimate = mean(estimate),
            max_pval = mean(p.value),
            mean_estimate.qnorm = mean(estimate.qnorm),
            relationship=unique(relationship),
            reproducible = length(unique(sign(estimate.qnorm)))==1 & all(significant),
            .groups = "drop")

# quantile options
# SLOW: https://stats.stackexchange.com/questions/50080/estimate-quantile-of-value-in-a-vector
# cume_dist: https://dplyr.tidyverse.org/reference/ranking.html
reproducible_models <- merged_models %>% 
  filter(reproducible) %>%
  group_by(model) %>%
  mutate(padj = p.adjust(max_pval, method="BH")) %>%
  ungroup() %>%
  left_join(lkup, by=c(feature.x = "gene_ID")) %>%
  mutate(gene_symbol = ifelse(is.na(gene_symbol),feature.x,gene_symbol)) %>%
  relocate(gene_symbol, .after = "feature.x") %>%
  group_by(model) %>%
  mutate(coef.quantile = cume_dist(abs(mean_estimate.qnorm))) %>%
  #mutate(z = scale(mean_estimate.qnorm)[,1]) %>% 
  #mutate(z.pvalue = 2*(1-pnorm(abs(z)))) %>%
  #mutate(z.padj = p.adjust(z.pvalue,method="BH")) %>%
  ungroup()

extreme_models <- reproducible_models %>%
  filter(coef.quantile > 0.9)


write_tsv(merged_models,snakemake@output[["merged_tsv"]])
write_tsv(reproducible_models,snakemake@output[["filtered_tsv"]])
write_tsv(extreme_models,snakemake@output[["extreme_tsv"]])
