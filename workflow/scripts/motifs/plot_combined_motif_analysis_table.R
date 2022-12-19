library(tidyverse)
library(gt)
library(gtExtras)
library(universalmotif)
library(webshot2)

memes <- ifelse(exists("snakemake"), snakemake@input[["memes"]], 
                "results/analysis/motifs/combined_xstreme.meme") %>% 
  read_meme() %>%
  set_names(., map_chr(., `@`, altname))



denovo <- ifelse(exists("snakemake"), snakemake@input[["denovo"]], 
                "results/analysis/motifs/combined_xstreme_results.tsv") %>% 
  read_tsv() %>%
  dplyr::select(-ALT_ID)

remap_enr <- ifelse(exists("snakemake"), snakemake@input[["remap_enr"]], 
                 "results/analysis/motifs/remap_peak_sea.tsv.gz") %>% 
  read_tsv()

archbold_compr <- ifelse(exists("snakemake"), snakemake@input[["archbold_compr"]], 
                    "results/analysis/motifs/archbold14_motif_comparison.rds") %>% 
  read_rds()

tab <- denovo %>%
  left_join(remap_enr, by=c("te_group","ID","CONSENSUS"), suffix = c(".TEs",".peaks")) %>%
  dplyr::select(te_group, peak_set, CONSENSUS, 
                starts_with("TP%"), starts_with("PVALUE"), ALT_ID) %>%
  mutate(across(starts_with("PVALUE"),p.adjust,method="BH",.names = "adj{.col}")) %>%
  filter(if_all(starts_with("adjPV"), ~{.x < 0.1})) %>%
  arrange(adjPVALUE.peaks) %>%
  mutate(motif = map(ALT_ID, ~memes[[.x]])) %>%
  mutate(logo  = map(motif, universalmotif::view_motifs))


for_table <- tab %>% 
  filter(te_group == "pan" & peak_set == "pan") %>%
  dplyr::select(logo, contains("TP"), contains("adjP")) %>%
  mutate(across(contains("%"), paste0, "%")) %>%
  mutate(across(contains("VALUE"), format.pval, 2)) %>%
  dplyr::rename_with(~str_squish(str_replace_all(.x,"TP|\\."," ")), .cols = contains("TP%")) %>%
  dplyr::rename_with(~str_replace_all(.x,"adjPVALUE\\.","padj "), .cols = contains("adj"))
  

gtg <- for_table %>%
  gt() %>%
  text_transform(
    locations = cells_body(columns = logo),
    fn = function(x) {
      map(for_table$logo, ggplot_image, height = px(75), aspect_ratio=2 )
    }
  ) %>%
  tab_header(title=md("SEA analysis"),
             subtitle = md("*de novo* motifs in *pan*-coexpressed TEs and *pan* REMAP peaks")) %>%
  cols_width(logo~px(200),
             everything()~px(100))

gtsave_extra(gtg, filename = snakemake@output[["png"]],zoom = 2, vwidth=600, vheight=1800)
save(gtg, for_table, file = snakemake@output[["rda"]])