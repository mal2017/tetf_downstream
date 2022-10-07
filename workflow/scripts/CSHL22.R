library(tidyverse)
library(ggpmisc)
library(ggprism)
library(ComplexHeatmap)
library(grid)
library(magick)
library(patchwork)
library(magrittr)
library(paletteer)
library(tidygraph)
library(ggraph)

FONTSIZE = 12

theme_set(theme_prism() + theme(axis.text.x=element_text(hjust=0.5,vjust=1,size = FONTSIZE),
                                axis.text.y = element_text(size=FONTSIZE),
                                axis.title = element_text(size=FONTSIZE),plot.title = element_text(size=FONTSIZE)))

lkup <- read_tsv("results/resources/gene_symbol_lookup.tsv.gz")

tfs <- read_tsv("data/Drosophila_melanogaster_TF.txt")

cofacs <- read_tsv("data/Drosophila_melanogaster_TF_cofactors.txt")

tx.related <- bind_rows(TF=tfs,cofactor=cofacs,
          .id="regulator.type") %>%
  dplyr::select(regulator.type,gene_symbol=Symbol,gene_ID=Ensembl,family=Family) %>%
  group_by(gene_ID) %>%
  slice_head(n=1) %>%
  ungroup()

dir.create("~/Downloads/CSHL22_figs")

# -----------------------------
# coexpression
# -----------------------------

# heatmap ----------------------------------------------------------------------
intro_heats <- read_rds("results/plots/intro_heatmaps.gg-list.rds")
svg("~/Downloads/CSHL22_figs/heatmap.svg",width = 7.5, height = 6.1)
draw(intro_heats$combined)
dev.off()


# scatter showing genes are either pos or neg ----------------------------------
intro_scatter <- read_rds("results/plots/plot_intro_ncoex_scatter.rds")
intro_scatter <- (intro_scatter$tes + ylab("N neg. coex. genes") + xlab("N pos. coex. genes"))  / 
                  (intro_scatter$genes + ylab("N neg. coex. TEs") + xlab("N pos. coex. TEs")) + 
  plot_layout(guides = 'collect') & 
  theme(legend.position = "bottom")
ggsave("~/Downloads/CSHL22_figs/ncoex_scatter.svg",intro_scatter, height = 6.5, width = 4.5)


# ------------------------------------------------------------------------------
# te communities
# ------------------------------------------------------------------------------

# community obs/expected ---------------------------------------------------------------
comm_obs_exp <- read_rds("results/analysis/direct_binding/by_community_te_kmer_dist.tbl.rds") %>%
  filter(algorithm == "leiden20" & model == "female_model_01") %>%
  mutate(sig = padj < 0.1 & in.dist < matched.dist.mean) %>%
  dplyr::select(community,observed = in.dist,expected = matched.dist.mean,sig)

cols <- comm_obs_exp %>% filter(sig) %>% mutate(color = as.character(paletteer::paletteer_d("ggsci::default_igv",n = n()))) %>%
  dplyr::select(community,color) %>%
  deframe()

g_comm_obs_exp <- comm_obs_exp %>%
  pivot_longer(c(-community,-sig), names_to = "metric", values_to = "mash.dist") %>%
  mutate(metric = fct_relevel(metric,"observed")) %>%
  ggplot(aes(metric,mash.dist)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color=community),size=3,width = 0.05) +
  geom_line(aes(group=community,color=sig)) +
  #scale_color_manual(values = c(`TRUE` = "red",`FALSE`="black")) +
  scale_color_manual(values = cols) +
  guides(color="none") +
  xlab("")

ggsave("~/Downloads/CSHL22_figs/te_comm_obs_exp.svg",
       g_comm_obs_exp, width = 5.45, height = 3.1)

# communities ---------------------------------------------------------------
te_comms <- read_rds("results/analysis/direct_binding/coregulated_te_communities.communities_tbl.rds") %>% filter(model=="female_model_01" & algorithm=="leiden20")
set.seed(2021)
g_te_comms <- read_rds("results/analysis/direct_binding/coregulated_te_communities.igraph_list.rds")$female_model_01 %>%
  as_tbl_graph() %>%
  left_join(te_comms,by=c(name="feature")) %>%
  left_join(comm_obs_exp,by="community") %>%
  ggraph(layout = "fr") +
  geom_edge_link2(check_overlap = T,alpha=0.2) +
  geom_node_point(size=rel(4),color="white") +
  geom_node_point(size=rel(4),aes(fill=community),shape=21,alpha=0.3) +
  geom_node_label(data=. %>% filter(sig),aes(label=name,fill=community),check_overlap = T,family="bold") +
  #scale_fill_paletteer_d("ggsci::default_igv") +
  scale_fill_manual(values=cols) +
  guides(fill="none") +
  theme_graph()

ggsave("~/Downloads/CSHL22_figs/te_comms.svg",
       g_te_comms, width = 11, height = 6)


# ------------------------------------------------------------------------------
#  TE obs/exp per gene 
# ------------------------------------------------------------------------------
intra_grp_dist_all <- read_rds("results/plots/plot_each_gene_te_kmer_dist.rds")


ggsave("~/Downloads/CSHL22_figs/tf_obs_exp.svg",
       intra_grp_dist_all$tf_obs_exp +
         ylab("observed mash dist") +
         xlab("expected mash dist") +
         ggrepel::geom_text_repel(data = . %>% filter(in.dist< 0.25 & gene_symbol %in% tx.related$gene_symbol),aes(label=paste0("italic('",gene_symbol,"')")),
                                  max.overlaps = 15,max.iter = 1000,
                                  size=FONTSIZE-8,parse=T),
       width = 5.45, height = 3.1)

# ------------------------------------------------------------------------------
# leading edge -> tx regulators
# ------------------------------------------------------------------------------

# gsea plots -------------------------------------------------------------------
gsea_gene_grp <- read_rds("results/analysis/signatures/gene_group_gsea.tbl.rds") %>%
  mutate(pathway = ifelse(pathway == "Tx.related","AnimalTFDB 3.0",str_wrap(pathway,20)))

plot_gene_grp_gsea <- function(grp, tbl = gsea_gene_grp) {
  s_grp<- tbl%>%
    dplyr::select(pathway,pval,padj,NES,ES) %>%
    mutate(label = paste0("NES=",round(NES,2),"\np=",format.pval(pval,3))) %>%
    filter(pathway %in% grp)
  
  tbl %>%
    filter(pathway %in% grp) %>%
    dplyr::select(pathway,gsea.plt_tbl) %>% 
    unnest(gsea.plt_tbl) %>%
    arrange(pathway,x) %>%
    ggplot(aes(x,y)) +
    geom_path(color="tomato",size=rel(2)) +
    facet_wrap(~pathway, ncol=1,scales="free") +
    geom_hline(yintercept = 0, linetype="dashed") +
    geom_text_npc(data=s_grp,aes(npcx=0.4,npcy=0.4,label=label),vjust="center",hjust="center",size=6) + 
    ylab("Enrichment Score") +
    xlab("rank") +
    scale_x_continuous(expand = expansion(0))
}


(plot_gene_grp_gsea(c("AnimalTFDB 3.0","C2H2 ZINC FINGER\nTRANSCRIPTION\nFACTORS")) + facet_wrap(~pathway,ncol=2)) %>%
  ggsave("~/Downloads/CSHL22_figs/animaltdb_gsea.svg", ., width = 9.5, height = 3.7)

# -- top decile n bar
read_rds("results/plots/plot_tx_related_in_lms.rds") %>%
  ggsave("~/Downloads/CSHL22_figs/tx_related_in_lms_bar.svg",.,width = 5.6, height = 4.2)

# ------------------------------------------------------------------------------
# binding
# ------------------------------------------------------------------------------
ggsave("~/Downloads/CSHL22_figs/with_grp_nullranges.svg",read_rds("results/plots/plot_within_tf_grp_nullranges.rds"),height = 4.2, width = 5.2)


# ------------------------------------------------------------------------------
#piRNA pathway
# ------------------------------------------------------------------------------

# ecdf plots showing shift for piRNA genes -------------------------------------
pirna_genes_in_lms.mean_coef_ecdf <- read_rds("results/plots/plot_pirna_genes_in_lms.mean_coef_ecdf.rds")
pirna_genes_in_lms.ntotal_ecdf <- read_rds("results/plots/plot_pirna_genes_in_lms.ntotal_ecdf.rds")

pirna_genes_in_lms.ntotal_ecdf$plot %<>% {. + ylab("prop")}
pirna_genes_in_lms.mean_coef_ecdf$plot %<>% {. + ylab("prop") + xlab("abs(coex. score)")}

get_pv_lab_from_ks_tbl <- .  %>%
  mutate(sex=ifelse(str_detect(model,"female"),"female","male")) %>%
  mutate(lab = paste0(sex," p=",format.pval(p.value,digits = 1))) 

pirna_genes_in_lms.mean_coef_ecdf$stat %<>% get_pv_lab_from_ks_tbl()
pirna_genes_in_lms.ntotal_ecdf$stat %<>% get_pv_lab_from_ks_tbl()

pirna_genes_in_lms.ntotal_ecdf$plot %<>% {.+ 
    annotate("text_npc",npcx = 0.5, npcy = 0.5,label = paste(c("2-sided KS test",pirna_genes_in_lms.ntotal_ecdf$stat$lab),collapse = "\n"),size=6) +
    scale_color_manual(values=c(male="red",female="black"))}
#scale_color_paletteer_d("nord::aurora")

pirna_genes_in_lms.mean_coef_ecdf$plot %<>% {.+ 
    annotate("text_npc",npcx = 0.5, npcy = 0.5,label = paste(c("2-sided KS test",pirna_genes_in_lms.mean_coef_ecdf$stat$lab),collapse = "\n"),size=6) +
    scale_color_manual(values=c(male="red",female="black"))}


(pirna_genes_in_lms.mean_coef_ecdf$plot +
    pirna_genes_in_lms.ntotal_ecdf$plot + xlab("N coex. TEs") +
    plot_layout(guides = 'collect') &
    theme(plot.margin = margin())
    #guides(color=guide_legend(nrow=2, byrow=TRUE), linetype=guide_legend(nrow=2, byrow=TRUE)) #
    #theme(legend.position = "bottom",axis.title = element_text(size=12),plot.margin = margin(r=10),
    #      legend.text = element_text(size=12))
  ) %>%
ggsave("~/Downloads/CSHL22_figs/pirna_ecdfs.svg",.,
       height = 3.2, width = 9)


# ------------------------------------------------------------------------------
# knockdowns
# ------------------------------------------------------------------------------


(read_rds("results/plots/plot_tfrnai_gsea.plot_list.rds")$ne_vs_p + 
    guides(color="none") + 
    scale_color_manual(values=c(`TRUE`="red",`FALSE`='gray'))) %>%
ggsave("~/Downloads/CSHL22_figs/tfrnai_coex_signature_gsea.nes.svg", .,
       width = 9, height =5)


# our KDs ----------------------------------------------------------------------

(read_rds("results/plots/plot_ourKD_gsea.plot_list.rds")$plot + 
   #guides(linetype="none") + 
    guides(color=guide_legend(nrow=4, byrow=TRUE),
           linetype=guide_legend(nrow=2,byrow=T))+
    theme(legend.position = "bottom",
          legend.text = element_text(size=8),
          legend.key.size = unit(5,"pt"),
          plot.margin = margin(r=10))) %>%
  ggsave("~/Downloads/CSHL22_figs/our_kd_gsea.svg",plot = .,width = 5, height = 4.2)



(read_rds("results/plots/plot_this_study_kd_deseq2.gg_list.rds")$NfI + theme(legend.position = c(0,0.9),legend.justification = "left") + xlab("log2(NfI-RNAi / GFP)")) %>%
  ggsave("~/Downloads/CSHL22_figs/our_kd_NfI.svg",.,width = 5, height = 4.2)
