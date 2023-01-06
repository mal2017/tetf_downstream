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
  mutate(ALT_ID = paste(ID,ALT_ID,sep='-'))

archbold_compr <- ifelse(exists("snakemake"), snakemake@input[["archbold_compr"]], 
                    "results/analysis/motifs/archbold14_motif_comparison.rds") %>% 
  read_rds()

# make a joinable column
archbold_compr <- archbold_compr %>%
  mutate(CONSENSUS = str_remove(target,"pan::\\d+-"))

tab <- denovo %>%
  left_join(archbold_compr) %>%
  filter(padj < 0.1) %>%
  dplyr::select(te_group, CONSENSUS, EVALUE, estimated_fdr, padj, ALT_ID, gg) %>%
  mutate(motif = map2(te_group, ALT_ID, ~memes[[paste0(.x,"::",.y)]])) %>%
  mutate(logo  = map(motif, universalmotif::view_motifs))

for_table <- tab %>% 
  group_by(CONSENSUS) %>%
  slice_min(padj) %>%
  ungroup() %>%
  filter(padj <0.1) %>%
  filter(te_group == "pan") %>%
  dplyr::select(logo, `E-value`=EVALUE, `est. FDR` = estimated_fdr, padj) %>%
  #mutate(across(contains("%"), paste0, "%")) %>%
  mutate(across(contains("VALUE")|contains("padj")|contains("FDR"), format.pval, 2))

gtg <- for_table %>%
  arrange(`est. FDR`) %>%
  gt() %>%
  text_transform(
    locations = cells_body(columns = logo),
    fn = function(x) {
      map(for_table$logo, ggplot_image, height = px(75), aspect_ratio=2 )
    }
  ) %>%
  tab_header(#title=md("SEA analysis"),
             title = md("*de novo* motifs"),
             subtitle = md("found in TEs coexpressed with *pan*")) %>%
  cols_width(logo~px(200),
             everything()~px(100)) %>%
  tab_spanner(label = "STREME", columns = c(`E-value`,`est. FDR`)) %>%
  tab_spanner(label='similarity to HMG motif', columns = c(`padj`)) %>%
  cols_align("center")

gtsave_extra(gtg, filename = snakemake@output[["png"]],zoom = 2, vwidth=600, vheight=1800)
save(gtg, for_table, tab, file = snakemake@output[["rda"]])