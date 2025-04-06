goal_type <- "SSQ"

goal_function <- function(X, E, Evar, return_pred = FALSE, use_penalization = TRUE){
  if (return_pred) return(X)
  penalization <- 0
  if (use_penalization){
    penalization <- sum(X < 0) * max(E) * 100000 - sum(X[X < 0])
    X <- abs(X)
  }
  ssq <- sum((X - E)^2)
  return(ssq + penalization)
}
