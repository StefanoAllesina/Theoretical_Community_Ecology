# Approximate
# B = D(1/d) + D(1/d) (v v^t) D(1/d)
# then
# B^(-1) = D + (v v^t) / (1 + 1^t(v^2 / d))

# build the matrix from parameters
get_B <- function(pars){
  n <- length(pars) * 0.5
  d <- pars[1:n]
  b <- pars[n + 1:n]
  return(diag(1/d) + (b/d) %o% (b/d))
}

# get inverse from parameters
get_B_inv <- function(pars){
  n <- length(pars) * 0.5
  d <- pars[1:n]
  b <- pars[n + 1:n]
  return(diag(d) - (b) %o% (b) / (1 + sum(b^2 / d)))
}

# predictions
get_pred <- function(x, d, b){
  p <- (x > 0) * 1
  d_tilde <- p * d
  b_tilde <- p * b
  theta <- sum(b_tilde) / (1 + sum(b_tilde^2 / d))
  return(d_tilde - b_tilde * theta)
}

# matrix of predictions given parameters
get_X <- function(E, pars){
  n <- ncol(E)
  d <- pars[1:n]
  b <- pars[n + 1:n]
  return(t(apply(E, 1, get_pred, d = d, b = b)))
}
