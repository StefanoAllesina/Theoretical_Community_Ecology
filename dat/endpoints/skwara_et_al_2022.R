library(tidyverse)
source("dat/endpoints/plot_results.R") # for plotting results
source("dat/endpoints/plot_pred_obs.R") # for plotting results

get_P <- function(E, pars){
  X <- get_X(E, pars)
  Bt <- t(get_B(pars))
  return(X %*% Bt)
}

find_best_approx_numerically <- function(pars, E, Evar, P, emmaxit, momentum){
  initialpars <- pars
  tmp <- list(par = pars)
  minimize_discrep <- function(pars, E, Evar, P){
    X <- P %*% t(get_B_inv(pars))
    return(goal_function(X = X, E = E, Evar = Evar, return_pred = FALSE, use_penalization = FALSE))
  }
  tmp <- optim(par = tmp$par, fn = minimize_discrep, E = E, Evar = Evar, P = P, 
               method = "Nelder-Mead",
               control = list(maxit = emmaxit, trace = FALSE))
  return(tmp$par * (1 - momentum) + momentum * initialpars)
}

EM_step <- function(pars, E, Evar, emmaxit = 500, momentum = 0){
  # compute P for current pars
  P <- get_P(E, pars)
  # find best pars for the current P
  pars <- find_best_approx_numerically(pars = pars, E = E, Evar = Evar, P = P, emmaxit = emmaxit, momentum = momentum)
  return(pars)
}

minimize_given_pars <- function(pars, E, Evar, D, replicates, return_pred = FALSE, use_penalization = TRUE){
  X <- get_X(pars = pars, E = D)
  X <- X[replicates, ]
  return(goal_function(X = X, E = E, Evar = Evar, return_pred = return_pred, use_penalization = use_penalization))
}

brute_force <- function(pars, E, Evar, D, replicates){
  tmp <- list(par = pars)
  for (i in 1:10){
    tmp <- optim(par = tmp$par, fn = minimize_given_pars, E = E, Evar = Evar, D = D, replicates = replicates,
                 method = "BFGS",
                 control = list(maxit = 2500, trace = FALSE))
    tmp <- optim(par = tmp$par, fn = minimize_given_pars, E = E, Evar = Evar, D = D, replicates = replicates,
                 method = "Nelder-Mead",
                 control = list(maxit = 5000, trace = FALSE))
    print(tmp$value)
  }
  return(tmp$par)
}

compute_variances <- function(E){
  comm <- apply(E, 1, function(x) paste0((x>0)*1, collapse = ""))
  Evar <- E %>% as_tibble() %>% add_column(comm = comm)
  Evar <- Evar %>% group_by(comm) %>% mutate_if(is.numeric, var) %>% 
    ungroup() %>% dplyr::select(-comm) %>% as.matrix()
  return(Evar)
}

single_run <- function(pars, E, Evar, D, replicates, skipEM, outfile, optimstep = 500){
  X <- get_X(D, pars)
  X <- X[replicates, ]
  # store goal function here (to be minimized)
  res <- goal_function(X = X, E = E, Evar = Evar, return_pred = FALSE, use_penalization = FALSE)
  # store best parameters found so far
  bestpars <- pars
  bestres <- res[1]
  if (!skipEM){
    # this is the EM
    for (i in 1:25){
      pars <- EM_step(pars = pars, E = E, Evar = Evar, emmaxit = optimstep + i * optimstep / 10, momentum = 0.)
      res <- c(res, minimize_given_pars(pars = pars, E = E, Evar = Evar, D = D, replicates = replicates, return_pred = FALSE, use_penalization = FALSE))
      P <- get_P(E, pars)
      if (res[i+1] < bestres){
        bestpars <- pars
        bestres <- res[i+1]
      }
      print(paste("iterative step", i, "->", bestres))
      P <- get_P(E, pars)
    }
  }
  print("numerical search")
  # now optimize numerically
  pars <- brute_force(pars = bestpars, E = E, Evar = Evar, D = D, replicates = replicates)
  res <- c(res, minimize_given_pars(pars = pars, E = E, Evar = Evar,  D = D, replicates = replicates, return_pred = FALSE))
  Epred <- minimize_given_pars(pars = pars, E = E, Evar = Evar,  D = D, replicates = replicates, return_pred = TRUE)
  # save results
  output <- list(
    data_name = outfile,
    observed = E,
    predicted = Epred,
    variances = Evar,
    B = get_B(pars),
    pars = pars,
    goal_function = minimize_given_pars(pars = pars, E = E, Evar = Evar, D = D, replicates = replicates, return_pred = FALSE),
    goal_type = goal_type
  )
  return(output)
}

run_model <- function(datafile, # location of the data
                     model = "diag_a11t", # one of "diag_a11t", "diag_vvt", "diag_vwt", or "full"
                     goalf = "SSQ", # one of SSQ, WLS
                     pars = NULL, # pre-computed parameters; otherwise use the identity matrix
                     skipEM = FALSE, # skip the iterative algorithm
                     plot_results = FALSE # plot the results
                     ){
  # load goal function
  if (goalf == "SSQ") source("dat/endpoints/ssq.R")
  if (goalf == "WLS") source("dat/endpoints/wls.R")
  if (goalf == "LikLN") source("dat/endpoints/LikLN.R")
  # load model structure
  print(paste0("dat/endpoints/model_", model, ".R"))
  source(paste0("dat/endpoints/model_", model, ".R"))
  # read the data
  E <- read.csv(datafile) %>% as.matrix()
  n <- ncol(E)
  # base name for output
  outfile <- tools::file_path_sans_ext(basename(datafile))
  print(outfile)
  # reorder the matrix by community
  community <- apply((E > 0) * 1, 1, paste, collapse = "")
  E <- E[order(community),]
  community <- community[order(community)]
  # compute variances
  Evar <- compute_variances(E)
  # compute design matrix D, and identify replicates
  D <- unique((E > 0) * 1)
  replicates <- as.numeric(as.factor(community))
  # if parameters are not provided, use the identity matrix
  if(is.null(pars)){
    if (model == "diag_a11t") pars <- c(rep(1,n), 0)
    if (model == "diag_vvt") pars <- c(rep(1,n), rep(0, n))
    if (model == "diag_vwt") pars <- c(rep(1,n), rep(0, 2 * n))
    if (model == "full") pars <- as.vector(diag(rep(1,n)))
  }
  output <- single_run(pars = pars, E = E, Evar = Evar, 
                       D = D, replicates = replicates, 
                       skipEM = skipEM, outfile = outfile)
  #if (plot_results) (plot_results_boxplot(output$observed, output$predicted))
  #save(output, file = paste0("dat/endpoints/results/", model, "_", outfile,"_", goal_type, ".Rdata"))
  return(output)
}

run_model_LOO <- function(datafile, # location of the data
                      model = "diag_a11t", # one of "diag_a11t", "diag_vvt", "diag_vwt", or "full"
                      goalf = "SSQ", # one of SSQ, WLS
                      LOO_row_num, # row number of the community to leave out: all communities of the same type will be left out as well
                      pars = NULL, # pre-computed parameters; otherwise use the identity matrix
                      skipEM = FALSE, # skip the iterative algorithm
                      plot_results = TRUE # plot the results
){
  # load goal function
  if (goalf == "SSQ") source("dat/endpoints/ssq.R")
  if (goalf == "WLS") source("dat/endpoints/wls.R")
  if (goalf == "LikLN") source("dat/endpoints/LikLN.R")
  # load model structure
  source(paste0("dat/endpoints/model_", model, ".R"))
  # read the data
  E <- read.csv(datafile) %>% as.matrix()
  n <- ncol(E)
  #E <- E / mean(E[E>0])
  # base name for output
  outfile <- tools::file_path_sans_ext(basename(datafile))
  # reorder the matrix by community
  community <- apply((E > 0) * 1, 1, paste, collapse = "")
  E <- E[order(community),]
  community <- community[order(community)]
  D <- unique((E > 0) * 1)
  replicates <- as.integer(as.factor(community))
  # compute variances
  Evar <- compute_variances(E)
  # now identify community to leave out
  labelcomm <- apply((E>0) * 1, 1, paste, collapse = "")
  infit <- labelcomm != labelcomm[LOO_row_num]
  Einfit <- E[infit == TRUE, , drop = FALSE]
  Evarinfit <- Evar[infit == TRUE, , drop = FALSE]
  comminfit <- apply((Einfit > 0) * 1, 1, paste, collapse = "")
  Dinfit <- unique((Einfit > 0) * 1)
  replicatesinfit <- as.numeric(as.factor(comminfit))
  
  # first run: identity matrix
  if (model == "diag_a11t") pars <- c(rep(1,n), 0)
  if (model == "diag_vvt") pars <- c(rep(1,n), rep(0, n))
  if (model == "diag_vwt") pars <- c(rep(1,n), rep(0, 2 * n))
  if (model == "full") pars <- as.vector(diag(rep(1,n)))
  
  output <- single_run(pars = pars, 
                       D = Dinfit, replicates = replicatesinfit,
                       E = Einfit, 
                       Evar = Evarinfit, skipEM = skipEM, outfile = outfile)
  
  pars <- output$pars
  
  # now that we have the best parameters, get predictions for the whole data
  source("dat/endpoints/model_full.R")
  Epred <- minimize_given_pars(pars = as.vector(solve(output$B)), E = E, D = D, replicates = replicates,
                               Evar = Evar, return_pred = TRUE)
  
  # save results
  output <- list(
    data_name = outfile,
    observed = E,
    predicted = Epred,
    variances = Evar,
    B = output$B,
    pars = pars,
    infit = infit,
    goal_function = minimize_given_pars(pars = as.vector(solve(output$B)), E = E, D = D, replicates = replicates,
                                        Evar = Evar),
    goal_type = goal_type
  )
  
  
  if (plot_results) show(plot_results_boxplot(output$observed, output$predicted, output$infit))
  save(output, file = paste0("dat/endpoints/results/", model, "_", outfile,"_", goal_type, "_LOO_", LOO_row_num, ".Rdata"))
  return(output)
}

