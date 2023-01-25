library(tidyverse)
library(corrr)
library(jsonlite)

cols <- c("model","feature.x","gene_symbol","feature.y",
          "estimate.qnorm","significant_model","significant_x")

# Read in data
main <- ifelse(exists("snakemake"), 
               snakemake@input$mods, 
               "upstream/final-models.collected-info.tsv.gz") %>%
  read_tsv(col_select = all_of(cols))

reps <- ifelse(exists("snakemake"), 
               snakemake@input$indep, 
               "upstream/final-models.collected-info.reps.tsv.gz") %>%
  read_tsv(col_select = all_of(cols))

# combine main an indep results
dat <- full_join(main,reps,suffix=c(".main",".rep"),
                 by=c("model","feature.x","gene_symbol","feature.y"))

# correlation between each dataset's 'x' (gene expression coefs) for three filtering approaches
dat <- list(unfiltered = . %>% return,
     replicated = . %>% filter(significant_x.main & significant_x.rep),
     main_data = . %>% filter(significant_x.main)) %>%
  map_df(exec, dat,.id="result_set") %>%
  nest(-model,-result_set) %>%
  mutate(cor.test = map(data,~broom::tidy(cor.test(~ estimate.qnorm.main + estimate.qnorm.rep, data=.x)))) %>%
  unnest(cor.test)


dat %>%
  nest(-result_set,-model) %>%
  nest(-model) %>%
  jsonlite::write_json(snakemake@output$json, prettify=T)

write_rds(dat, snakemake@output[["rds"]])
