library(tidyverse)
library(ggalluvial)

#mods_path <- "upstream/final-models.collected-info.tsv.gz"
mods_path <- snakemake@input[["mods"]]

mods <- vroom::vroom(mods_path)

dat <- mods %>%
  filter(significant_model) %>%
  count(model,across(starts_with("adj_p.value_anova"), ~{.<0.1}))

# true in a col denotes that that term was significant
dat <- dat %>% 
  rename_with(~str_remove(.,".*p.value_anova_"), starts_with("adj_p.value_anova"))

dat <- dat %>%
  mutate(across(where(is.logical),~if_else(.x,"sig.","n.s.")))

dat2 <- dat %>%
  mutate(combo = row_number()) %>%
  pivot_longer(-c(n,combo,model))

dat2 <- dat2 %>%
  mutate(name=fct_relevel(name,c("overlap","wolbachia","scaled.copies.y","x")))

g <- ggplot(dat2,
       aes(x=name,stratum=value,alluvium=combo,y = n,label=value,fill=value)) +
  scale_x_discrete(expand = expansion()) +
  geom_flow() +
  geom_stratum() +
  geom_text(stat = "stratum",min.y=10000) +
  facet_wrap(~model, scales = "free")


# g <- ggplot(dat,
#       aes(axis1=wolbachia,axis2 = overlap, axis3 = scaled.copies.y,axis4=x,
#           y = n)) +
#  scale_x_discrete(limits = c("wolbachia","overlap","scaled.copies.y","x"), expand = c(.2, .05)) +
#  geom_alluvium(width = 1/4) +
#  geom_stratum() +
# geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
#  facet_wrap(~model, scales = "free")

write_rds(g,snakemake@output[["rds"]])
ggsave(snakemake@output[["png"]],g)