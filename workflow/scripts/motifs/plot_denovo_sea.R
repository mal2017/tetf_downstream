library(tidyverse)


path <- Sys.glob("results/analysis/motifs/xstreme_per_tf/*/sea_out/sea.tsv")
path <- snakemake@input[["seas"]]

df <- map_df(path,read_tsv,comment = "#")


df <- mutate(df, DB = str_extract(DB,regex("(?<=tf\\/).*?(?=\\/streme_out)")))


df %>%
  filter(QVALUE < 0.1) %>%
  count(DB)%>%
  arrange(n)

g <- df %>%
  #filter(QVALUE < 0.1) %>%
  arrange(-ENR_RATIO) %>%
  ggplot(aes(log2(ENR_RATIO),-log10(PVALUE),color=QVALUE < 0.1,label=DB)) +
  geom_point() +
  ggrepel::geom_text_repel(data= . %>% filter(QVALUE < 0.1)) +
  scale_color_grey(start = 0.6,end=0.4)

g + theme_classic()

#read_tsv("results/analysis/motifs/xstreme_per_tf/Aef1/sea_out/sea.tsv")
