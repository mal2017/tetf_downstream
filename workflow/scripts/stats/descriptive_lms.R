library(tidyverse)
library(jsonlite)

# Read in data
dat <- ifelse(exists("snakemake"), 
               snakemake@input$mods, 
               "upstream/final-models.collected-info.tsv.gz") %>%
  read_tsv()


# add new columns that are the result of dividing each column with "sumsq_" in the name by the value of `total_variance`
dat %>% 
    mutate(across(c(contains("sumsq_")), ~ . / total_variance,.names = "var_exp_{.col}")) -> dat

# ---------------------------------------------------------------------------------------------------------
# get basic stats for significant models
# ---------------------------------------------------------------------------------------------------------

# for models where significant_model == TRUE, calculate sd, mean, median, min, max, for all columns startining with 'var_exp_'
# and all columns with 'r2' but not 'p.value' or 'padj' in the name
sig_mod_stats <- dat %>% 
    filter(significant_model == TRUE) %>%
    bind_rows(.,mutate(., model="all")) %>%
    dplyr::select(-contains("p.value")) %>%
    group_by(model) %>%
    summarise(across(contains("var_exp_") | contains("r2"), .fns = funs(sd, mean, median, min, max, .args = list(na.rm=T))))

sig_mod_stats <- sig_mod_stats %>% 
  pivot_longer(-model, names_to = "statistic", values_to = "value")

# ---------------------------------------------------------------------------------------------------------
# get basic stats for significant x feature models
# ---------------------------------------------------------------------------------------------------------

where_x_sig_mod_stats <- dat %>% 
  filter(significant_model & significant_x) %>%
  bind_rows(.,mutate(., model="all")) %>%
  dplyr::select(-contains("p.value")) %>%
  group_by(model) %>%
  summarise(across(contains("var_exp_") | contains("r2"), .fns = funs(sd, mean, median, min, max, .args = list(na.rm=T)))) %>%
  pivot_longer(-model, names_to = "statistic", values_to = "value")

# ------------------------------------------------------------------------------
# percentage significant terms among valid models
# ------------------------------------------------------------------------------
get_pct_sig<- . %>%
  group_by(model) %>%
  filter(significant_model) %>%
  summarize(across(starts_with("adj"),~{sum(.<0.1, na.rm = T)/n()}))
  

significant_term_pct <- bind_rows(dat %>% get_pct_sig(),
          dat %>% mutate(model = "all") %>% get_pct_sig()) %>%
  rename_with(~str_remove(.x, "adj_p.value_anova_|adj_p.value_ftest_")) %>%
  pivot_longer(-model, names_to = "statistic", values_to = "value") %>%
  mutate(statistic = paste0("pct_sig_",statistic))

# ------------------------------------------------------------------------------
# N
# ------------------------------------------------------------------------------
significant_models <- dat %>% 
  filter(significant_model) %>% 
  dplyr::select(feature.x,feature.y) %>%
  distinct()

# calculate the total number of unique feature.y entries
dat %>% dplyr::select(feature.y) %>%
  distinct() %>%
  nrow() %>%
  list(n_features.y_tested = . ) -> n_features.y

# calculate the total number of unique feature.x entries
dat %>% dplyr::select(feature.x) %>%
  distinct() %>%
  nrow() %>%
  list(n_features.x_tested =.)-> n_features.x

#
n_valid_nonrepro <- dat %>% 
  filter(valid & !reproducible) %>%
  dplyr::anti_join(significant_models) %>%
  dplyr::select(feature.x, feature.y) %>%
  distinct() %>%
  nrow() %>%
  list(n_valid_nonrepro= .)
  

# Calculate total number of unique pairs assessed
dat %>% dplyr::select(feature.x, feature.y) %>%
  distinct() %>%
  nrow()  %>%
  list(n_total_pairs_tested = .)-> n_pairs

n_retained <- dat %>% filter(significant_model) %>%
  dplyr::select(feature.x, feature.y) %>%
  distinct() %>%
  nrow() %>%
  list(n_retained = .)

n_removed <- (n_pairs$n_total_pairs_tested - n_retained$n_retained) %>%
  list(n_removed = .)

basic_n <- c(n_features.y, n_features.x, n_pairs, n_removed, n_retained, n_valid_nonrepro) %>% 
  enframe(name="statistic") %>%
  mutate(model = "all",
         stat_group = "basic_n") %>%
  mutate(value=unlist(value))

# -----------------------------
# combine stats for json writing
# ---------------

to_write_json <- bind_rows(basic_n = basic_n,
                           sig_mod_stats = sig_mod_stats,
                           where_x_sig_mod_stats = where_x_sig_mod_stats,
                           significant_term_pct = significant_term_pct,
                           .id="stat_group") %>%
  nest(-model,-stat_group) %>%
  arrange(model,stat_group)

# makes it wider so that each statistic is a named entry with a value pair
to_write_json <- to_write_json %>% mutate(data = map(data, ~pivot_wider(.x,names_from = "statistic", values_from = "value")))

write_json(to_write_json, snakemake@output$json, pretty=T)

