library(deSolve)

f <- function(t, x, parms){
  return(list(c(
    1 - x[2],
    -1 + x[1]
  )))
}

tt <- seq(0, 10, by = 0.01)
out <- ode(y = c(1/2,6), func = f, times = tt, parms = NULL, method = "ode45")

x <- out[,2]
y <- out[,3]
