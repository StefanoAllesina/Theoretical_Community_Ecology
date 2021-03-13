# goal: 
# given a contact network (with few nodes), build the graph described in
# Fig 1 of Virus Spread in Networks, by Van Mieghem et al. 2009
# https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=4549746
#

library(igraph)

build_SIS_network <- function(A){
  n <- nrow(A) 
  states <- 0:(2^n-1)
  vectors <- sapply(states, function(x) as.numeric(intToBits(x)[1:n]))
  labels <- apply(vectors, 2, function(x) paste(x, collapse = ""))
  nstates <- length(states)
  edges <- matrix("", 0, 2)
  # edgelist
  for (i in 1:(nstates-1)){
    for (j in (i+1):nstates){
     vi <- vectors[,i] 
     vj <- vectors[,j]
     # if they are one infection apart
     diff <- sum(abs(vi - vj))
     if (diff == 1){
       if (sum(vj - vi) == 1){
         # case 1: j contains i
         # add an edge from j to i (recovery)
         edges <- rbind(edges, c(labels[j], labels[i]))
         # now find who is being infected
         infected <- which.max(vj - vi)
         # check whether this can be infected by a node infected in i
         for (l in 1:n){
           if (vi[l] == 1 & A[l, infected] == 1){
             edges <- rbind(edges, c(labels[i], labels[j]))
             break
           }
         }
       } else{
         # case 2: i contains j
         # add an edge from i to j (recovery)
         edges <- rbind(edges, c(labels[i], labels[j]))
         # now find who is being infected
         infected <- which.max(vi - vj)
         # check whether this can be infected by a node infected in i
         for (l in 1:n){
           if (vj[l] == 1 & A[l, infected] == 1){
             edges <- rbind(edges, c(labels[j], labels[i]))
             break
           }
         }
       }
     }
    }
  }
  g_transitions <- graph_from_edgelist(edges, directed = TRUE)
  return(g_transitions)
}

set.seed(3)
n <- 3
# draw a small contact network
A <- matrix(0, n, n)
# star network
A[1,] <- 1
# make symmetric
A <- A + t(A)
# remove self-loops
diag(A) <- 0
# make binary
A <- (A > 0) * 1

# graph of social network
g_social <- graph_from_adjacency_matrix(A)
g_transitions <- build_SIS_network(A)