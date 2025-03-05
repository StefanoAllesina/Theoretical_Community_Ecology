library(tidyverse)
library(deSolve)

# out1 (and optionally out2) are the ouput of ode from deSolve
plot_dynamics <- function(out1, out2){
  toplot <- out1 %>% 
    as.data.frame() %>% 
    pivot_longer(cols = -time) %>% 
    add_column(ts = "a")
  if (!is.null(out2)){
    toplot <- toplot %>% bind_rows(
      out2 %>% 
        as.data.frame() %>% 
        pivot_longer(cols = -time) %>% 
        add_column(ts = "b")
    )}
  pl <- toplot %>% 
    ggplot(aes(x = time, y = value, colour = name, linetype = ts)) + 
    geom_line(linewidth = 1) + 
    theme_bw() + 
    theme(legend.position = "none") + 
    scale_y_sqrt()
  return(pl)
}
