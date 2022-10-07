foi <- list(TF = "data/Drosophila_melanogaster_TF.txt",
            cofactor = "data/Drosophila_melanogaster_TF_cofactors.txt") %>%
  map_df(read_tsv, .id="type") %>%
  dplyr::select(type,gene_symbol=Symbol,gene_ID = Ensembl, family=Family) %>%
  arrange(desc(type)) %>%
  group_by(gene_ID) %>%
  slice_head(n=1) %>% ungroup()

oi <- c("pan","ct","awd","mamo","Unr","NfI","vvl","CG16779")
