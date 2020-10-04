library(deSolve)
THRESH <- 10^-9

cr_sep <- function(t, x, pars){
  x[x < THRESH] <- 0
  v <- pars$v
  s <- pars$s
  B <- pars$B
  dxdt <- (v * x) * (s - B %*% x)
  return(list(dxdt))
}

n <- 6
# find feasible, stable parameterization
success <- FALSE
while(success == FALSE){
  B <- matrix(runif(n * n), n, n)
  B <- B + t(B)
  # make stable
  eB <- eigen(B, symmetric = TRUE, only.values = TRUE)$values
  B <- B - diag(rep(min(Re(eB) * 1.1), n))
  # efficiencies
  v <- runif(n)
  # choose feasible equilibrium
  xstar <- runif(n)
  # choose s such that xstar is equilibrium
  s <- as.numeric(B %*% xstar)
  if (all(s > 0)) success <- TRUE
}
pars <- list(s = s, v = v, B = B)
# integrate dynamics
x0 <- runif(n) + 0.05
out <- ode(y = x0, times = seq(0, 1000, by = 0.1), 
           func = cr_sep, 
           parms = pars, method = "ode45")
# Lyapunov function
lyap <- function(out, pars){
  V <- numeric(0)
  s <- pars$s
  B <- pars$B
  for (i in 1:nrow(out)){
    x <- out[i, -1]
    V <- c(V,
           2* sum(x * s) - as.numeric(t(x) %*% B %*% x)
           )
  }
  return(V)
}

V <- lyap(out, pars)
