library(tidyverse)
library(rtracklayer)
library(plyranges)
library(ggpointdensity)
library(ggdensity)
library(fgsea)

# ----------------------------------------------------------------------
# Stuff to note - just in case there's a weird bias like "TFs are higher expressed than average"
# or "TF quantification is lower certainty than average",
# I only consider a possible universe of genes that makes it into our final filtered LM set
# when I perform the gsea. The results are similar regardless.



#lms_path <- "results/analysis/coexpression/filtered_models.tsv.gz"
#tfs_path <- "data/Drosophila_melanogaster_TF.txt"
#cofacs_path <- "data/Drosophila_melanogaster_TF_cofactors.txt"
#pirna_path <- "results/resources/pirna_pathway.tsv"


lms_path <- snakemake@input[["lms"]]
tfs_path <- snakemake@input[["tfs"]]
cofacs_path <- snakemake@input[["cofacs"]]
pirna_path <- snakemake@input[["pirna"]]

lms <- read_tsv(lms_path)

tfs <- bind_rows(TF = read_tsv(tfs_path),
          cofactor = read_tsv(cofacs_path),
          .id="class") %>%
  dplyr::select(gene_id = Ensembl,gene_symbol=Symbol,family=Family,class) %>%
  filter(gene_id %in% lms$feature.x)

piRNAgenes <- read_tsv(pirna_path) %>%
  filter(gene_ID %in% lms$feature.x)

# grps <- read_tsv("http://ftp.flybase.net/releases/FB2022_04/precomputed_files/genes/gene_group_data_fb_2022_04.tsv.gz",skip=7) %>%
#   mutate(class=case_when(str_detect(FB_group_name,regex("kinase",ignore_case = T)) ~ "kinase",
#                          str_detect(FB_group_name,regex("GTPase",ignore_case=T)) ~ "GTPase",
#                          str_detect(FB_group_name,regex("transporter",ignore_case=T)) ~ "transporter",
#                          str_detect(FB_group_name,regex("peptidase",ignore_case=T)) ~ "peptidase",
#                          str_detect(FB_group_name,regex("nuclease",ignore_case=T)) ~ "nuclease",
#                          str_detect(FB_group_name,regex("ligase",ignore_case=T)) ~ "ligase")) %>%
#   dplyr::select(class,gene_id=Group_member_FB_gene_id) %>%
#   drop_na() %>%
#   distinct()


grps <- read_tsv("http://ftp.flybase.net/releases/FB2022_04/precomputed_files/genes/gene_group_data_fb_2022_04.tsv.gz",skip=7) %>%
  dplyr::select(class=FB_group_name,gene_id=Group_member_FB_gene_id) %>%
  drop_na() %>%
  distinct() %>%
  filter(gene_id %in% lms$feature.x)

# ------------------------------------------------------------------------------
# plot gsea for TFs
# ------------------------------------------------------------------------------

# not filtering by extreme coefs
to_plot <- lms %>% 
  group_by(feature.x) %>%
  summarise(coef = mean(abs(mean_estimate.qnorm)),n_tes = dplyr::n(),.groups = "drop") %>%
  mutate(Tx.related = ifelse(feature.x %in% tfs$gene_id,"TF/coact.","other"))

pathways <- list(Tx.related = unique(tfs %>% pull(gene_id)),piRNA=piRNAgenes$gene_ID) %>%
  c(.,
    grps %>% split(.,.$class) %>%
  map(pull,gene_id))

# exlcude small groups
pathways <- pathways[pathways %>% map_lgl(~{length(.)>=10})]

ranks <- to_plot %>%
  dplyr::select(feature.x,coef) %>%
  deframe

set.seed(2022)
gsea_res <- fgsea(pathways,ranks,scoreType="pos") %>% 
  as_tibble() %>%
  arrange(pval)

gsea_res_plt <- gsea_res %>%
  mutate(gsea.plt_tbl = map(pathway,~as_tibble(plotEnrichment(pathways[[.x]],ranks)$data)))



write_rds(gsea_res_plt,snakemake@output[["rds"]])


# 
# gsea_toplot %>%
#   mutate(gene.group)
#   ggplot(aes(x,y)) +
#   geom_path(color="tomato",size=rel(2)) +
#   geom_hline(yintercept = 0, linetype="dashed") +
#   #geom_text(data=labels,aes(label=label),vjust="top",hjust="right",size=rel(7)) + 
#   ylab("Enrichment Score") +
#   xlab("rank") +
#   scale_x_continuous(expand = expansion(0)) +
#   facet_wrap(~gene.group)

