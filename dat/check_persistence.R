# Find a persistent GLV that has a feasibile but unstable equilibrium
find_persistent <- function(n, includebound = FALSE) {
  success <- FALSE
  while (!success) {
    A <- -matrix(abs(rnorm(n ^ 2)), n, n)
    r <- abs(rnorm(n)) + 0.1
    # A <- -matrix(sample(1:15, n^2), n, n)
    # r <- -rowSums(A)
    n <- nrow(A)
    # check feasible not stable
    if (abs(det(A)) < 10^-10) {
      xs <- rep(0, n)
    } else {
      xs <- solve(A, -r)
      l1 <- -1
      if (all(xs > 0)) {
        l1 <- max(Re(eigen(diag(xs) %*% A)$values))
      }
      if (l1 > 0) {
        boundary <- matrix(0, 2 ^ n - 2, n)
        for (i in 1:(2 ^ n - 2)) {
          status <- as.integer(intToBits(i)[1:n])
          Ak <- A[status > 0, status > 0, drop = FALSE]
          rk <- r[status > 0]
          if(abs(det(Ak)) > 10^-10){
          xk <- solve(Ak, -rk)
          if (all(xk > 0))
            boundary[i, status > 0] <- xk
          }
        }
        # all boundary equil are positive
        allp <- all(boundary >= 0)
        if (allp) {
          all_components <- round(t(A %*% t(boundary) + r), 14)
          tooptim <- function(p) {
            p <- exp(p) + 0.01
            mysum <- rowSums(all_components %*% diag(p))
            return(-sum(mysum[mysum < 0]))
          }
          tmp <- optim(rep(1, n), tooptim)
          if (tmp$value == 0)
            success <- TRUE
        }
      }
    }
  }
  if (includebound) return(list(A = A, r = r, boundary = boundary))
  return(list(A = A, r = r))
}