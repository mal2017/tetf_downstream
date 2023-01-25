library(tidyverse)

fl<-ifelse(exists("snakemake"),
           snakemake@input$rds, 
           "results/analysis/specific_genes/remap_peaks_near_pirna_genes_contingency.rds")

read_rds(fl) %>%
  filter(d == 1000) %>%
  mutate(contingency_table = map(contingency_table,~rownames_to_column(as.data.frame(.x)))) %>%
  unnest(contingency_table) %>%
  pivot_longer(c(bound,unbound), names_to="promoter", values_to = "n")%>%
  mutate(promoter = fct_rev(promoter)) %>%
  filter(rowname=="piRNA") %>%
  ggplot(aes(ChIP,n,fill=promoter)) +
  geom_col(position = "fill") +
  xlab("ChIP") + ylab("promoters bound") +
  scale_fill_grey(start = 0.8, end=0.4) +
  theme_classic() +
  scale_y_continuous(labels = scales::percent)
