# Approximate
# B = D(1/d) + D(1/d) (v w^t) D(1/d)
# then
# B^(-1) = D + (v w^t) / (1 + 1^t(v w / d))
# where w is a unit vector (not enforced right now)

# build the matrix from parameters
get_B <- function(pars){
  n <- length(pars) /3
  d <- pars[1:n]
  b <- pars[n + 1:n]
  v <- pars[2 * n + 1:n]
  #v <- v / sqrt(sum(v^2))
  return(diag(1/d) + (b/d) %o% (v/d))
}

# get inverse from parameters
get_B_inv <- function(pars){
  n <- length(pars) /3
  d <- pars[1:n]
  b <- pars[n + 1:n]
  v <- pars[2 * n + 1:n]
  #v <- v / sqrt(sum(v^2))
  return(diag(d) - (b) %o% (v) / (1 + sum(b * v / d)))
}

# predictions
get_pred <- function(x, d, b, v){
  p <- (x > 0) * 1
  d_tilde <- p * d
  b_tilde <- p * b
  theta <- sum(p * v) / (1 + sum(b_tilde * v / d))
  return(d_tilde - b_tilde * theta)
}

# matrix of predictions given parameters
get_X <- function(E, pars){
  n <- ncol(E)
  d <- pars[1:n]
  b <- pars[n + 1:n]
  v <- pars[2 * n + 1:n]
  #v <- v / sqrt(sum(v^2))
  return(t(apply(E, 1, get_pred, d = d, b = b, v = v)))
}
