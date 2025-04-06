goal_type <- "WLS"

goal_function <- function(X, E, Evar, return_pred = FALSE, use_penalization = TRUE){
  if (return_pred) return(X)
  vE <- as.vector(E)
  vX <- as.vector(X)
  vVar <- as.vector(Evar)
  # remove zeros
  vX <- vX[vE > 0]
  vVar <- vVar[vE > 0]
  vE <- vE[vE > 0]
  penalization <- 0
  if (use_penalization){
    # penalize for negative predictions
    penalization <- sum(vX < 0) * max(E) * 100000 - sum(vX[vX < 0])
    vX <- abs(vX)
  }
  # compute weighted sum of squares
  res_sq <- (vX - vE)^2
  WLS <- sum(res_sq / vVar)
  return(WLS + penalization)
}
