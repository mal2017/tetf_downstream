library(tidyverse)

# Read in data
main <- ifelse(exists("snakemake"), 
              snakemake@input$mods, 
              "upstream/final-models.collected-info.tsv.gz") %>%
  read_tsv()

# Read in data
reps <- ifelse(exists("snakemake"), 
              snakemake@input$reps, 
              "upstream/final-models.collected-info.reps.tsv.gz") %>%
  read_tsv()

dat <- bind_rows(main=main, reps=reps, .id="dataset")

replicated <- ifelse(exists("snakemake"), 
              snakemake@input$replicated, 
              "results/analysis/independent_dataset/replicate_dataset_correlation.rds") %>%
  read_rds()

# bool values indicate significance in the main dataset and the replicate, independent dataset
replicated2 <- replicated %>% filter(result_set == "unfiltered") %>% dplyr::select(model, data) %>%
  unnest(data) %>%
  dplyr::select(model, feature.x, feature.y, starts_with("significant_x"))

# we want to be able to compute leading edge enrichments for rankings by the following:

# variance explained by each variable - abs and signed
var_exp <- dat %>% 
  mutate(across(contains("sumsq"),~{.x/total_variance},.names = "var_exp_{.col}")) %>%
  dplyr::rename_with(.cols = contains("var_exp_sumsq"),.fn = ~str_remove(.x,"sumsq_anova_")) %>%
  mutate(across(starts_with("var_exp_"),~{.x * sign(estimate.qnorm)}, .names="signed_{.col}")) %>%  
  dplyr::select(dataset,model,feature.x,gene_symbol,feature.y,contains("var_exp_x"))

# -log10 pvalue - abs and signed
logpv <- dat  %>%
  mutate(across(starts_with("p.value_anova_x"),~{-log10(.x)},.names = "log10_{.col}")) %>% 
  mutate(across(starts_with("log10"),~{.x * sign(estimate.qnorm)}, .names="signed_{.col}")) %>%
  dplyr::select(dataset,model,feature.x,gene_symbol,feature.y, contains("log10_"))

# coex scores - abs and signed
coex <- dat %>% 
  dplyr::select(dataset,model,feature.x,gene_symbol,feature.y,contains("estimate")) %>%
  mutate(across(starts_with("estimate"),abs, .names="abs_{.col}")) %>%
  dplyr::select(dataset,model,feature.x,gene_symbol,feature.y, contains("estimate"))

# combine these rankings into 1 tbl, containing both sexes and both datasets
dat <- full_join(coex,var_exp, by = c("dataset", "model", "feature.x", "gene_symbol", "feature.y")) %>%
  full_join(logpv, by = c("dataset", "model", "feature.x", "gene_symbol", "feature.y")) %>%
  left_join(replicated2,c("model", "feature.x", "feature.y")) %>%
  relocate(starts_with("significant_x"))

filtering_funs <- c(both = . %>% filter(significant_x.main & significant_x.rep),
                    main = . %>% filter(significant_x.main),
                    replicate = . %>% filter(significant_x.rep),
                    nofilt = . %>% return)

summary_funs <- list(mean = . %>% mean(na.rm=T))
                  #max = . %>% max(na.rm=T))

res <- filtering_funs %>%
  map_df(exec, dat,.id="filtering_approach") %>%
  nest(-dataset, -model, -filtering_approach)

res <- res %>% 
  filter(!(filtering_approach == "replicate" & dataset=="main")) %>%
  filter(!(filtering_approach == "main" & dataset=="reps"))

res2 <- res %>%
  mutate(stats = map(data, ~summarise(group_by(.x, gene_id = feature.x,gene_symbol),across(where(is.numeric), .fns=c(mean=mean), .names="{.fn}_{.col}", na.rm=T),.groups = "drop"))) %>%
  dplyr::select(-data)

res2 <- res2 %>%unnest(stats)
  
res2 <- res2 %>%
  pivot_longer(-c(dataset, model, filtering_approach,  gene_id, gene_symbol), 
               names_to = "metric", values_to = "value") %>%
  nest(gene_id, gene_symbol, value)

write_rds(res2, snakemake@output$rds)
