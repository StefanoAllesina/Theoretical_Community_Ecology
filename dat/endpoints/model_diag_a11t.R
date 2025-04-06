# Approximate B
# B = D(1/d) + a1 1^t 
# where a is a scalar
# then
# B^(-1) = D + D1 1^t D / (1 + a 1 D 1^t)
#        = D + d d^t / (1 + a 1^t d)

# build the matrix from parameters
get_B <- function(pars){
  n <- length(pars) -1
  d <- pars[1:n]
  a <- pars[n + 1]
  ones <- rep(1, n)
  return(diag(1/d) + a * ones %o% ones)
}

# get inverse from parameters
get_B_inv <- function(pars){
  n <- length(pars) -1
  d <- pars[1:n]
  a <- pars[n + 1]
  ones <- rep(1, n)
  return(diag(d) - a * (d) %o% (d) / (1 + a * sum(d)))
}

# predictions
get_pred <- function(x, d, a){
  p <- (x > 0) * 1
  d_tilde <- p * d
  theta <- 1 + a * sum(d_tilde)
  return(d_tilde/ theta)
}

# matrix of predictions given parameters
get_X <- function(E, pars){
  n <- length(pars) -1
  d <- pars[1:n]
  a <- pars[n + 1]
  ones <- rep(1, n)
  return(t(apply(E, 1, get_pred, d = d, a = a)))
}
