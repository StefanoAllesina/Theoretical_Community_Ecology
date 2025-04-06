library(tidyverse)

plot_results_boxplot <- function(E, Epred, inout = NULL){
  if (is.null(inout)) inout <- rep(FALSE, nrow(Epred))
  # input E Epred (optional) inout
  spnames <- colnames(E)
  #spnames <- sapply(spnames, substr, 1, 6)
  colnames(E) <- spnames
  colnames(Epred) <- spnames
  presence <- E > 0
  community <- apply(presence, 1, function(x) paste(spnames[x], collapse = "-"))
  tmp1 <- E %>% as_tibble() %>% 
    mutate(ID = row_number()) %>% add_column(community = community, inout = inout) %>% 
    pivot_longer(names_to = "species", values_to = "observed", cols = -c("ID", "community", "inout"))
  tmp2 <- Epred %>% as_tibble() %>% 
    mutate(ID = row_number()) %>% add_column(community = community, inout = inout) %>% 
    pivot_longer(names_to = "species", values_to = "predicted", cols = -c("ID", "community", "inout"))
  tmp <- inner_join(tmp1, tmp2, by = c("ID", "community", "inout", "species"))
  tmp <- tmp %>% filter(observed > 0)
  pl <- tmp %>% 
    ggplot(aes(x = species, y = observed, colour = species)) +
    facet_wrap(~community, scales = "free_y") + 
    #facet_wrap(~community) + 
    geom_boxplot(aes(alpha = 0.6 - 0.2 * inout, fill = species)) + 
    stat_summary(fun = mean, geom="point", shape=2, size=5) +
    geom_point() + 
    scale_alpha(range = c(0.05,0.5)) + 
    geom_point(data = tmp %>% dplyr::select(species, predicted, community) %>% distinct(), 
             aes(x = species, y = predicted), colour = "black", size = 3, shape = 1) + 
    theme_bw() + theme(legend.position = "none") + scale_y_log10() + 
    scale_fill_brewer(palette = "Set1") + 
    scale_colour_brewer(palette = "Set1")
  return(pl)
}
