library(tidyverse)
library(rtracklayer)

tfs <- ifelse(exists("snakemake"),snakemake@input[["tfs"]],
               "resources/Drosophila_melanogaster_TF.txt") %>%
  read_tsv()

lms <- ifelse(exists("snakemake"),snakemake@input[["mods"]],
    "upstream/final-models.collected-info.tsv.gz") %>% 
    read_tsv() %>%
    filter(significant_x) %>%
    filter(feature.x %in% tfs$Ensembl)

# finds number of TEs per TF
lms <- lms %>% 
    dplyr::select(gene_symbol,feature.y) %>%
    distinct() %>%
    group_by(gene_symbol) %>%
    summarise(n = n()) %>%
    arrange(desc(n)) %>%
    left_join(lms,.) %>%
    filter(n>5)



fa <- ifelse(exists("snakemake"),snakemake@input[["tes"]],
        "resources/Tidalbase_transposon_sequence.fasta") %>%
        import()

fa <- fa[names(fa) %in% lms$feature.y]

seqs <- lms %>% 
    dplyr::select(gene_symbol,feature.y) %>%
    distinct() %>%
    split(.,.$gene_symbol) %>%
    map(pull,feature.y) %>%
    map(.f = ~{list(coex = fa[.x],other = fa[!names(fa) %in% .x])})

odir0 <- ifelse(exists("snakemake"),snakemake@output[["odir"]],
               "results/analysis/motifs/consensus_tes_per_tf/")

for (s in names(seqs)) {
    message(s)
    sq <- seqs[[s]]
    
    for (k in names(sq)) {
        odir <- sprintf("%s/%s",odir0,s)
        dir.create(odir,recursive = T)
        export(sq[[k]],paste0(odir,"/",k,".fasta"))
    }
}
