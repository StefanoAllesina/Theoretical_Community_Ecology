n <- 10
k <- 12

P <- matrix(runif(n * k), n, k)
B <- -diag(runif(n))
v <- runif(k)
Z <- matrix(0, k, k)

A <- cbind(B, -P)
A <- rbind(A, cbind(diag(v) %*% t(P), Z))