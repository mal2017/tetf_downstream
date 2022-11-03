library(tidyverse)
library(rtracklayer)
library(plyranges)

#nullranges_path <- "results/analysis/direct_binding/per_pair_nullranges.rds"
nullranges_path <- snakemake@input[["nullranges_path"]]

res <- read_rds(nullranges_path) %>%
  arrange(p.value) %>%
  mutate(padj = p.adjust(p.value,method="BH")) %>%
  ungroup()

res.sig <- res %>% filter(padj < 0.1 & estimate > 1)

res.sig %>%
  mutate(TE = fct_infreq(TE)) %>%
  ggplot(aes(TE)) +
  geom_bar()


res %>%
  filter(estimate > 1 & TE=="INE-1") %>%
  ggplot(aes(estimate,-log10(p.value))) +
  geom_point() +
  geom_point(data= . %>% filter(padj < 0.1),shape=21,color="red",size=4) +
  ggrepel::geom_text_repel(data= . %>% filter((padj < 0.1)),aes(label=TF),max.overlaps = 50)


remap <- import("data/remap2022_nr_macs2_dm6_v1_0.bed.gz")

dre4 <- remap %>% filter(str_detect(name,"dre4")) %>%
  reduce() %>%
  as_tibble()


fixed.dre4 <- read_rds("results/resources/annotated_fixed_insertions.gr.rds") %>%
  filter(repeat_element == "INE-1" & dre4)

fixed.notdre4 <- read_rds("results/resources/annotated_fixed_insertions.gr.rds") %>%
  filter(repeat_element == "INE-1" & !dre4)


fixed.dre4 %>%
  mutate(name=ix) %>%
  export("~/Downloads/fixed_ine1_dre4.bed")

fixed.dre4.tbl <- fixed %>%
  as_tibble()

fixed.dre4.tbl%>%
  ggplot(aes(seqnames)) +
  geom_bar()
