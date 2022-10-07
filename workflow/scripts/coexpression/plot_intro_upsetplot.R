library(tidyverse)
library(ComplexUpset)

#mods_path <- "results/analysis/coexpression/filtered_models.tsv.gz"
mods_path <- snakemake@input[["mods"]]

mods <- read_tsv(mods_path)

to_plot <- mods %>% 
  filter(model %in% c("male_model_01","female_model_01")) %>%
  mutate(sex = ifelse(model == "male_model_01","M","F")) %>%
  unite(pair,feature.x,feature.y,sep = "~") %>%
  mutate(centile = ifelse(coef.quantile >0.9,">90 centile",NA)) %>%
  dplyr::select(sex, pair,relationship, centile) %>%
  mutate(class = str_squish(paste(sex,relationship)), extreme = ifelse(!is.na(centile),str_squish(paste(sex,centile)),NA)) %>%
  pivot_wider(names_from = sex, values_from = c(class,extreme), names_sep = ".") %>%
  #mutate(class = case_when(relationship.F == relationship.M ~paste0("both ",relationship.M),
  #                         relationship.F != relationship.M ~paste0("M ",relationship.M,"; F ", relationship.F),
  #                         is.na(relationship.M)~paste0("F ",relationship.F," only"),
  #                         is.na(relationship.F)~paste0("M ",relationship.M," only"),
  #                         T~"other")) %>%
  group_by(pair) %>%
  summarize(class = pmap(list(class.F,class.M,extreme.F,extreme.M),.f=function(w,x,y,z){c(w,x,y,z) %>%.[!is.na(.)]}), .groups = "drop")


to_plot2 <- to_plot %>%
  unnest(class) %>%
  mutate(ph = T) %>%
  pivot_wider(names_from = "class", values_from = ph,values_fill = F)


g <- upset(to_plot2,colnames(to_plot2)[-1],
           set_sizes = F,
           encode_sets = T,height_ratio = 0.25,
           mode="distinct",
           base_annotations = list('Size'=intersection_size(
             mode="distinct",
             size=0,
             counts=F)),
           wrap = F,matrix = intersection_matrix(geom=geom_point(size=1)), min_size=50)


#g <- ggplot(to_plot,aes(class)) +
#  geom_bar() +
#  scale_x_upset() +
#  scale_y_log10()

# export -----------------------------------------------------------------------

write_rds(list(plot=g,data=to_plot2),snakemake@output[["rds"]])

# plots  -----------------------------------------------------------------------
# g <- mods %>% 
#   count(model,type) %>%
#   mutate(type = fct_reorder(type,-n)) %>%
#   ggplot(aes(n,type,fill=model)) +
#   geom_col(position="dodge")
# 
# 
# mods2 %>%
#   group_by(model,type,feature.x) %>%
#   dplyr::select(model,type,family,gene_symbol.x) %>%
#   distinct() %>%
#   group_by(model,type) %>%
#   tally() %>%
#   mutate(type = fct_reorder(type,-n)) %>%
#   #filter(type!="other") %>%
#   ggplot(aes(type,n,fill=model)) +
#   geom_col(position="dodge")
# 
# 
# 
# 
# mods2 %>%
#    filter(type!="other") %>%
#    unite(col = "family2",family,type,sep=".",remove = F) %>%
#    dplyr::select(model,type,family,family2, gene_symbol.x) %>%
#   distinct() %>%
#    count(type,model,family,family2,sort = T) %>%
#   filter(n >1) %>%
#   ggplot(aes(reorder(family2,n),n,fill=model)) +
#   geom_col(position="dodge") +
#   coord_flip()
  



