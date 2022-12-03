library(tidyverse)

#rew_fl <-"results/analysis/rewiring/candidate_rewiring_events.rds" 
rew_fl <- snakemake@input[["rds"]]

rew <- read_rds(rew_fl)

g <- rew %>%
  ggplot(aes(rho_focus_strains, rho_other_strains)) +
  geom_point(data = . %>% filter(padj >= 0.1),color="gray") +
  geom_point(data = . %>% filter(padj < 0.1),color="red") +
  ggrepel::geom_text_repel(data = . %>% filter(rho_focus_strains > 0 & rho_other_strains < 0 & padj < 0.1),aes(label=paste(gene_symbol,feature.y,nearby_gene,sep=">")))

saveRDS(g,snakemake@output[["rds"]])
ggsave(snakemake@output[["png"]],g)


