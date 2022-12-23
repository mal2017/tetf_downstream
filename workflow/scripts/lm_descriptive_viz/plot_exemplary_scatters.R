library(tidyverse)
library(ggdensity)

# Read in data
mods <- ifelse(exists("snakmeake"), 
               snakemake@input$mods, 
               "upstream/final-models.collected-info.tsv.gz") %>%
  read_tsv()

expression <- list(male = "upstream/male.0.expression.tsv.gz",
                   female =  "upstream/female.0.expression.tsv.gz") %>% 
  map_df(read_tsv,.id="sex") %>%
  pivot_longer(-c(feature,sex),names_to = "sample", values_to = "expression") %>%
  drop_na()

lkup <- mods %>% dplyr::select(feature=feature.x, gene_symbol) %>%
  distinct()

expression <- left_join(expression, lkup)

oi <- c("pan","NfI","CG16779")

dat <- full_join(expression %>% filter(gene_symbol %in% oi),
                 expression %>% filter(!str_detect(feature,"FBgn")),
                 by=c("sample","sex"), suffix = c(".x",".y")) %>%
  inner_join(filter(mods, significant_x), by = c("sex"="model", "feature.x", "feature.y"))


plot_scatter <- function(z) {
  ggplot(z,aes(expression.x, expression.y)) +
    #ggdensity::geom_hdr() + 
   #ggdensity::geom_hdr_points(size=0.01) +
    geom_point(size=0.1) +
    geom_smooth(method="lm", se=F, linewidth=0.5) +
    facet_wrap(sex~feature.y, scales="free") +
    theme(aspect.ratio = 1)
}

gs <- dat %>% split(.,.$gene_symbol) %>%
  map(plot_scatter)

write_rds(gs,snakemake@output[["rds"]])
