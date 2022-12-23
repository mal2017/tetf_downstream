library(tidyverse)
library(ggpmisc)
library(gt)
library(gtExtras)

FONTSIZE = 7L

theme_set(theme_classic() + theme(axis.title = element_text(size=rel(5)), axis.text = element_text(size=rel(2))))

gsea_path <- ifelse(exists("snakemake"),snakemake@input[["gene_group_gsea"]],
                    "results/analysis/signatures/gene_group_gsea.tbl.rds")

gsea_tbl <- read_rds(gsea_path) %>%
  mutate(pathway = ifelse(pathway == "Tx.related","AnimalTFDB 3.0",pathway))

plot_gsea <- function(d) {
  ggplot(d,aes(x,y)) +
  geom_path(linewidth=rel(1.1)) +
  geom_hline(yintercept = 0, linetype="dashed") +
  #ggrepel::geom_text_repel(max.overlaps = 10, max.iter = 100,color="black",fontface="italic",force_pull = 0.1,
  #                         data=. %>% 
  #                           filter(p.adjust < 0.1) %>% group_by(RNAi) %>% slice_max(abs(runningScore),n = 1, with_ties = F),
  #                         aes(label=str_wrap(RNAi,10)),vjust="bottom",hjust="left",size=unit(FONTSIZE-2,"pt")) + 
  ylab("NES") +
  xlab("rank")
}

plts_df <- gsea_tbl %>%
  filter(padj < 0.1) %>%
  mutate(gg = map(gsea.plt_tbl,.f=plot_gsea)) %>%
  arrange(pval) %>%
  dplyr::select(pathway,size,pval,padj,NES,ES,gg) %>%
  mutate(enrichment = NA)

gtg <- plts_df %>%
  dplyr::select(-gg, -pval, -ES) %>%
  mutate(padj = map(padj,format.pval,3),
         NES = map(NES, round, 3)) %>%
  gt() %>%
  text_transform(
    locations = cells_body(columns = enrichment),
    fn = function(x) {
      map(plts_df$gg, ggplot_image, height = px(75), aspect_ratio=2 )
    }
  ) %>%
  tab_header(title=md("Flybase gene group enrichment"), subtitle = md("among genes coexpressed with TEs")) %>%
  cols_width(enrichment~px(400),
             pathway~px(400),
             everything()~px(200))

gtsave_extra(gtg, filename = snakemake@output[["png"]],zoom = 2, vwidth=1400, vheight=3600)
save(gtg, plts_df, file = snakemake@output[["rda"]])