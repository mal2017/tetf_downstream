library(tidyverse)
library(clusterProfiler)
library(ggnewscale)

# ideas

res <- ifelse(exists("snakemake"), snakemake@input$rds, "results/analysis/enrichment/fb_gg_gsea.rds") %>% read_rds

res %>% filter(filtering_approach == "main" & dataset == "main" & model == "female" & metric == "mean_abs_estimate.qnorm") %>%
  pull(gsea) %>%
  .[[1]] %>% 
  enrichplot::gseaplot2(geneSetID = "POLYBROMO-CONTAINING BRAHMA ASSOCIATED PROTEINS COMPLEX")

res %>% filter(filtering_approach == "main" & dataset == "main" & model == "female" & metric == "mean_abs_estimate.qnorm") %>%
  pull(gsea) %>%
  .[[1]] %>% 
  enrichplot::gseaplot2(geneSetID = "BRAHMA ASSOCIATED PROTEINS COMPLEX")


## C2H2 GSEA plot



res %>% filter(filtering_approach == "main" & dataset == "main" & model == "female" & metric == "mean_abs_estimate.qnorm") %>%
  pull(gsea) %>%
  .[[1]] %>% 
  enrichplot::gseaplot2(geneSetID = "C2H2 ZINC FINGER TRANSCRIPTION FACTORS")

res %>% filter(filtering_approach == "main" & dataset == "main" & model == "female" & metric == "mean_abs_estimate.qnorm") %>%
  pull(gsea) %>%
  .[[1]] %>% 
  enrichplot::gseaplot2(geneSetID = "ZAD_ZNF")


res %>% filter(filtering_approach == "main" & dataset == "main" & model == "male" & metric == "mean_abs_estimate.qnorm") %>%
  pull(gsea) %>%
  .[[1]] %>% 
  enrichplot::gseaplot2(geneSetID = "C2H2 ZINC FINGER TRANSCRIPTION FACTORS")

res %>% filter(filtering_approach == "main" & dataset == "main" & model == "male" & metric == "mean_abs_estimate.qnorm") %>%
  pull(gsea) %>%
  .[[1]] %>% 
  enrichplot::gseaplot2(geneSetID = "ZAD_ZNF")


res %>% filter(filtering_approach == "replicate" & dataset == "reps" & model == "male" & metric == "mean_abs_estimate.qnorm") %>%
  pull(gsea) %>%
  .[[1]] %>% 
  enrichplot::gseaplot2(geneSetID = "C2H2 ZINC FINGER TRANSCRIPTION FACTORS")

res %>% filter(filtering_approach == "replicate" & dataset == "reps" & model == "female" & metric == "mean_abs_estimate.qnorm") %>%
  pull(gsea) %>%
  .[[1]] %>% 
  enrichplot::gseaplot2(geneSetID = "C2H2 ZINC FINGER TRANSCRIPTION FACTORS")

res %>% filter(filtering_approach == "replicate" & dataset == "reps" & model == "female" & metric == "mean_abs_estimate.qnorm") %>%
  pull(gsea) %>%
  .[[1]] %>% 
  enrichplot::gseaplot2(geneSetID = "ZAD_ZNF")







# we want a plot that gives an overview of the main data set, a network plot, and an example of TFs being reproducible