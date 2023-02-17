library(tidyverse)
library(clusterProfiler)
library(ggnewscale)

res <- ifelse(exists("snakemake"), snakemake@input$rds, "results/analysis/enrichment/fb_gg_gsea.rds") %>% read_rds

ranking_selector_funs <- list(main_male = . %>% filter(filtering_approach == "main" & dataset == "main" & model == "male" & metric == "mean_abs_estimate.qnorm"),
     main_female = . %>% filter(filtering_approach == "main" & dataset == "main" & model == "female" & metric == "mean_abs_estimate.qnorm"),
    indep_male = . %>% filter(filtering_approach == "replicate" & dataset == "reps" & model == "male" & metric == "mean_abs_estimate.qnorm"), 
     indep_female = . %>% filter(filtering_approach == "replicate" & dataset == "reps" & model == "female" & metric == "mean_abs_estimate.qnorm")) 

# sanity check
#pull(test,gsea)[[1]] %>% as_tibble() %>% filter(str_detect(ID,"C2H2")) %>% pull(core_enrichment) %>% str_count('/')
make_dotplot <- function(.x,.y) {
  dotplot(.x,showCategory=50, x="NES", size="str_count(core_enrichment,'/')") + 
  theme_linedraw() +
  labs(size="leading edge members") +
    ggtitle(.y)
}

ranking_selector_funs %>%
  map_df(exec,res, .id="label") %>%
  mutate(gg_dot = map2(gsea, label, make_dotplot)) %>%
  pull(gg_dot)




res %>%
  filter(metric == "mean_estimate.qnorm") %>%
  dplyr::select(model, filtering_approach, dataset, gsea.tidy) %>%
  unnest(gsea.tidy) %>%
  filter(ID == "ZAD_ZNF")
