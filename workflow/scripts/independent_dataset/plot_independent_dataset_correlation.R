library(tidyverse)
library(ggdensity)

# Read in data
dat <- ifelse(exists("snakemake"), 
               snakemake@input$rds, 
               "results/analysis/independent_dataset/replicate_dataset_correlation.rds") %>%
  read_rds()


plot_scatter <- function(data, label) {
  ggplot(data,aes(estimate.qnorm.main, estimate.qnorm.rep)) +
    ggdensity::geom_hdr_points() +
    geom_smooth(method="lm") +
    xlab("main results") +
    ylab("independent dataset") +
    annotate("text", -Inf, Inf, label = str_wrap(label,width = 20), hjust = 0, vjust = 1)
}


res <- dat %>% 
  filter(result_set!="unfiltered") %>%
  mutate(label = paste0(str_extract(method,"Pearson's|Spearman's",)," rho=",round(estimate,digits = 3),"; ",alternative," p=",format.pval(p.value,digits=3))) %>%
  mutate(label=str_replace(label,"=<","<")) %>%
  dplyr::relocate(label) %>%
  mutate(gg= map2(data, label, plot_scatter))


write_rds(res,snakemake@output$rds)
  
