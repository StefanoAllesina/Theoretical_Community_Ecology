# Maximize likelihood when sampling from lognormal with mean x_i^(k) and variance sigma_i
goal_type <- "LikLN"

goal_function <- function(X, E, Evar, return_pred = FALSE, use_penalization = TRUE){
  if (return_pred) return(X)
  tominimize <- function(sdlog, obs, pred){
    -sum(dlnorm(obs, meanlog = log(pred), sdlog = sdlog, log = TRUE))
  }
  maximize_likelihood_column <- function(x, e){
    x <- x[e > 0]
    e <- e[e > 0]
    penalization <- sum(x < 0)
    x[x < 0] <- 1.0
    # find the best standard deviation
    tmp <- optimise(tominimize, interval = c(0,5), obs = e, pred = x)
    return(tmp$objective + penalization * 1000)
  }
  negative_log_likelihoods <- sapply(1:ncol(E), function(k) maximize_likelihood_column(X[,k], E[,k]))
  return(sum(negative_log_likelihoods))
}
