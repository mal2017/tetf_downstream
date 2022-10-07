library(tidyverse)
library(tidytext)

#tfs_path <- "data/Drosophila_melanogaster_TF.txt"
#cofacs_path <- "data/Drosophila_melanogaster_TF_cofactors.txt"
#mods_path <- "results/analysis/coexpression/filtered_models.tsv.gz"
#lkup_path <- "results/resources/gene_symbol_lookup.tsv.gz"

lkup_path <- snakemake@input[["lkup"]]
mods_path <- snakemake@input[["mods"]]
tfs_path <- snakemake@input[["tfs"]]
cofacs_path <- snakemake@input[["cofacs"]]

lkup <- read_tsv(lkup_path)

mods <- read_tsv(mods_path)

tfs <- bind_rows(TF = read_tsv(tfs_path),
                 cofactor = read_tsv(cofacs_path),
                 .id="class") %>%
  dplyr::select(gene_id = Ensembl,gene_symbol=Symbol,family=Family,class)

# 
# out_of <- mods %>%
#   dplyr::select(feature.x) %>%
#   distinct() %>%
#   left_join(tfs,by=c(feature.x="gene_id")) %>%
#   mutate(across(where(is.character),replace_na,"other")) %>%
#     count(class,name = "out.of")


x1 <- mods %>%
  #filter(coef.quantile > 0.9) %>% 
  left_join(tfs,by=c(feature.x="gene_id",gene_symbol="gene_symbol")) %>%
  mutate(family = replace_na(family,"other"), class = replace_na(class,"other")) %>%
  group_by(gene_symbol,feature.y) %>%
  slice_max(coef.quantile,n=1,with_ties = F) %>% # when represented in both models, only use best
  ungroup() %>%
  select(model,gene_symbol,family,class,feature.y,coef.quantile) %>%
  mutate(bin = cut_interval(coef.quantile,length = 0.1))

# g <- x1 %>%
#   mutate(class = fct_relevel(class,rev(c("TF","cofactor","other")))) %>%
#   mutate(sex=str_extract(model,".*male")) %>%
#   ggplot(aes(bin,fill=class)) +
#   geom_bar(position="fill",width=1,color="black") +
#   #facet_wrap(~sex,nrow=2) +
#   scale_y_continuous(labels=scales::percent, expand = expansion()) +
#   scale_fill_manual(values=c(TF="red",cofactor="black",other="gray")) +
#   theme(axis.text.x = element_text(angle=45, hjust=1)) +
#   xlab("coex. score percentile") +
#   ylab("% coexpressed gene/TE pairs")


# -----------
all_genes <- lkup %>%
  filter(gene_type == "protein_coding_gene") %>%
  dplyr::select(gene_ID,gene_symbol) %>%
  distinct() %>%
  left_join(tfs,by=c(gene_symbol="gene_symbol",gene_ID="gene_id")) %>%
  mutate(across(where(is.character),replace_na,"other"))

to_plot <- x1 %>% 
  slice_max(bin) %>%
  count(gene_symbol,class) %>%
  left_join(all_genes,.) %>%
  mutate(n=replace_na(n,0)) %>%
  mutate(is.coex = n>0) %>%
  count(class,is.coex) %>%
  group_by(class) %>%
  mutate(prop = n/sum(n)) %>%
  filter(is.coex) %>%
  ungroup()
  
g <- to_plot %>%
  mutate(class = fct_reorder(class,prop)) %>%
  ggplot(aes(class,prop,fill=class)) +
  geom_col(color="black") +
  scale_y_continuous(labels=scales::percent, expand = expansion()) +
  scale_fill_manual(values=c(TF="red",cofactor="black",other="gray")) +
  ggtitle("% genes with >=1 coex. score > 90th percentile") +
  ylab("%") + xlab("") +
  guides(fill="none")

# to_plot <- x1  %>%
#   group_by(gene_symbol,family,class,bin) %>%
#   summarise(n=n(),.groups = "drop")  %>%
#   #filter(n >= 2) %>%
#   slice_max(bin)
# 
# to_label <- to_plot %>%
#   count(class,bin) %>%
#   slice_max(bin,n = 1) %>%
#   arrange(n) %>%
#   mutate(label = paste0(n," ",class,"s"))
# 
# 
# g <- to_plot%>%
#   count(bin,class) %>%
#   left_join(out_of) %>%
#   mutate(prop = n/out.of)
#   ggplot(aes(fct_infreq(class))) +
#   geom_bar() +
#   xlab("gene type") + 
#   annotate("text",x=2.5, y = 1000,label = paste(to_label$label,collapse = '\n')) +
#   ggtitle("genes with 90th percentile coex. scores")
  

write_rds(g,snakemake@output[["rds"]])


# how many