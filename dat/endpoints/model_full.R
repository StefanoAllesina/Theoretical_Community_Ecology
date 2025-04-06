# build the matrix from parameters
get_B <- function(pars){
  n <- sqrt(length(pars))
  Binv <- matrix(pars, n, n)
  # avoid problems with singular matrices
  try({return(solve(Binv))}, silent = TRUE)
  return(diag(rep(1,n)))
}

# get inverse from parameters
get_B_inv <- function(pars){
  n <- sqrt(length(pars))
  Binv <- matrix(pars, n, n)
  return(Binv)
}

# predictions
get_pred <- function(x, B){
  Bk <- B[x>0, x>0, drop = FALSE]
  # avoid problems with singular matrices
  x[x > 0] <- -1
  try({
    x[x < 0] <- rowSums(solve(Bk))
  },silent = TRUE)
  return(x)
}

# matrix of predictions given parameters
get_X <- function(E, pars){
  B <- get_B(pars)
  return(t(apply(E, 1, get_pred, B = B)))
}
