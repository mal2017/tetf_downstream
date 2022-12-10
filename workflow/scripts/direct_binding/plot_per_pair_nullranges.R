library(tidyverse)
library(rtracklayer)
library(plyranges)

#nullranges_path <- "results/analysis/direct_binding/within_tf_nullranges.rds"
nullranges_path <- snakemake@input[["rds"]]

res <- read_rds(nullranges_path) %>%
  arrange(p.value) %>%
  ungroup()

res.sig <- res %>% filter(p.value < 0.05 & estimate > 1)

# res.sig %>%
#   mutate(TE = fct_infreq(TE)) %>%
#   ggplot(aes(TE)) +
#   geom_bar()

#res.sig %>%
#  mutate(TF = fct_infreq(TF)) %>%
#  ggplot(aes(TF)) +
#  geom_bar() + theme(axis.text.x = element_text(angle=45,hjust=1))

g <-res.sig %>%
  mutate(TF=fct_infreq(TF)) %>%
  mutate(TE=fct_infreq(TE)) %>%
  ggplot(aes(TF,TE,size=-log10(p.value))) +
  geom_point() +
  theme(axis.text.x  = element_text(angle=45,hjust=1))

write_rds(g, snakemake@output[["rds"]])

ggsave(snakemake@output[["png"]],g)


#res.sig %>%
#  filter(TF == "pan") %>%
#  dplyr::select(-cont_mat,-is.coex) %>%
#  gt::gt()


#res.sig %>%
#  filter(TE == "INE-1") %>%
#  dplyr::select(-cont_mat,-is.coex) %>%
#  gt::gt()



# res %>%
#   filter(estimate > 1 & TE=="INE-1") %>%
#   ggplot(aes(estimate,-log10(p.value))) +
#   geom_point() +
#   geom_point(data= . %>% filter(padj < 0.1),shape=21,color="red",size=4) +
#   ggrepel::geom_text_repel(data= . %>% filter((padj < 0.1)),aes(label=TF),max.overlaps = 50)
# 









# 
# 
# 
# 
# remap <- import("data/remap2022_nr_macs2_dm6_v1_0.bed.gz")
# 
# dre4 <- remap %>% filter(str_detect(name,"dre4")) %>%
#   reduce() %>%
#   as_tibble()
# 
# 
# fixed.dre4 <- read_rds("results/resources/annotated_fixed_insertions.gr.rds") %>%
#   filter(repeat_element == "INE-1" & dre4)
# 
# fixed.notdre4 <- read_rds("results/resources/annotated_fixed_insertions.gr.rds") %>%
#   filter(repeat_element == "INE-1" & !dre4)
# 
# 
# fixed.dre4 %>%
#   mutate(name=ix) %>%
#   export("~/Downloads/fixed_ine1_dre4.bed")
# 
# fixed.dre4.tbl <- fixed %>%
#   as_tibble()
# 
# fixed.dre4.tbl%>%
#   ggplot(aes(seqnames)) +
#   geom_bar()
