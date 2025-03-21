mu <- 10

get_roots <- function(rho, mu) {
  # roots of the polynomial
  # c0 + c1 x + c2 x^2 + c3 x^3
  c3 <- -rho / mu
  c2 <- rho
  c1 <- -1 - rho / mu
  c0 <- rho
  polyroot(c(c0, c1, c2, c3))
}

rho <- seq(0.3, 0.6, by = 0.001)
r <- t(sapply(rho, get_roots, mu = mu))
r <- as.data.frame(r)
colnames(r) <- c("r1", "r2", "r3")
r <- cbind(r, rho)
r <-
  r %>% pivot_longer(names_to = "root",
                     values_to = "equilibrium",
                     cols = -rho)
r <- r %>% filter(abs(Im(equilibrium)) < 10 ^ -6)
r <-
  r %>% mutate(equilibrium = Re(equilibrium)) %>% filter(equilibrium > 0)
pl <-
  ggplot(r, aes(x = rho, y = equilibrium, colour = root)) + geom_line() + theme_bw() + theme(legend.position = "none")
