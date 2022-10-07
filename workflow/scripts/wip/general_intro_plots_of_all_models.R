# establishing normality is reasonable
ggplot(mods,aes(sample=mean_estimate.qnorm)) +
  stat_qq() +
  stat_qq_line() +
  facet_wrap(~model)



# filtered_models %>%
#   group_by(model,feature.x,relationship) %>%
#   mutate(n=n()) %>%
#   ggplot(aes(n,mean_estimate.qnorm)) +
#   geom_point() +
#   facet_grid(model~relationship)
# 
# 
# 
# 
# filtered_models %>%
#   dplyr::select(model,mean_estimate.qnorm,feature.x,feature.y) %>%
#   pivot_wider(names_from = model, values_from = mean_estimate.qnorm,values_fill = 0) %>%
#   arrange(-female_model_01)
#   ggplot(aes(female_model_01, male_model_01)) +
#   geom_point()
# 
# filtered_models %>%
#   count(feature.x,model, relationship) %>%
#   pivot_wider(names_from = relationship, values_from = n,values_fill = 0) %>%
#   ggplot(aes(neg,pos)) +
#   facet_wrap(~model) +
#   geom_point()
#   
# 
# filtered_models %>%
#   filter(feature.x %in% c("FBgn0267033","FBgn0085432","FBgn0037698","FBgn0004198","FBgn0000150","FBgn0042696","FBgn0263352","FBgn0086680")) %>%
#   count(feature.x,model, relationship) %>%
#   pivot_wider(names_from = model, values_from = n,values_fill = 0) %>%
#   ggplot(aes(male_model_01,female_model_01)) +
#   facet_wrap(~relationship) +
#   geom_point() +
#   geom_label(aes(label=feature.x))