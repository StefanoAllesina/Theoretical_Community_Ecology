plot_pred_obs <- function(output){
  # assign communities and take mean observed
  cm <- apply((output$observed>0) * 1, 1, paste, collapse = "")
  tmp <- output$predicted %>% as_tibble() %>% add_column(community = cm) %>% 
    mutate(experiment = row_number()) %>% 
    pivot_longer(cols = -c(experiment, community), names_to = "strain", values_to = "predicted") %>% 
    inner_join(output$observed %>% as_tibble() %>% add_column(community = cm) %>% 
                 mutate(experiment = row_number()) %>% 
                 pivot_longer(cols = -c(experiment, community), names_to = "strain", values_to = "observed"),
               by = join_by(experiment, strain, community))
  ggplot(tmp %>% filter(observed > 0), 
         aes(x = predicted, y = observed, colour = strain)) + 
      geom_point(alpha = 0.5, shape = 4) + 
      geom_abline(intercept = 0, slope = 1, linetype = 2) + 
      theme_bw() + 
    geom_point(data = tmp %>% filter(observed > 0) %>% group_by(strain, community) %>% mutate(
      predicted = mean(predicted), observed = mean(observed)), 
      shape = 16, size = 2)
  
}