library(tidyverse)

# Read in data
mods <- ifelse(exists("snakmeake"), 
               snakemake@input$mods, 
               "upstream/final-models.collected-info.tsv.gz") %>%
  read_tsv()

mods2 <- mods %>% 
  filter(significant_model) %>%
  dplyr::select(model,feature.x, gene_symbol, feature.y, contains("adj") & (contains("wolbachia") | contains("copies")))

# convert padj to bool sig/nsig
mods2 <- mods2 %>% 
  mutate(across(contains("adj"),~{.<0.1}))

# gets mode of lgl vec
get_lgl_mode <- function(x) {
  names(which.max(table(x)))
}

# get the mode of each TE's results for the selected coefs
mods2 <- mods2 %>%
  group_by(model, feature.y) %>%
  mutate(across(contains("adj"),get_lgl_mode,.names = "{.col}_mode"))

# still grouped - new var if the significance of each entry matches the
# mode of its group
mods2 <- mods2 %>%
  mutate(wolbachia_consistent = adj_p.value_anova_wolbachia==adj_p.value_anova_wolbachia,
         copies_consistent = adj_p.value_anova_scaled.copies.y == adj_p.value_anova_scaled.copies.y_mode)

mods2 <- mods2 %>%
  summarise(n=n(), across(contains("consistent"),sum)) %>%
  mutate(across(contains("consistent"), ~{.x/n})) %>%
  arrange(wolbachia_consistent) %>%
  pivot_longer(contains("consistent"),names_to = "term", values_to = "consistent results per TE") %>%
  mutate(term = str_remove(term,"_consistent"))

g <- mods2 %>%
  ggplot(aes(term,`consistent results per TE`)) +
  geom_jitter() +
  facet_wrap(~model) +
  scale_y_continuous(labels = scales::percent)

write_rds(g,snakemake@output[["rds"]])
ggsave(snakemake@output[["png"]],g)
