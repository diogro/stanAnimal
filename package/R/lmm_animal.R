#' Fit linear mixed animal model
#' 
#' Fits a univariate or multivariate linear mixed animal model using a know relatedness matrix
#' @param Y response variable, univariate or multivariate
#' @param X fixed model matrix.
#' @param A relatedness matrix.
#' @param iter number of MCMC iterations.
#' @param warmup number of warmup iterations, must be smaller than iter
#' @param chains number of chains to run.
#' @param control list of stan control arguments.
#' @param cmd_stan Logical, use cmd_stan. If false, uses rstan.
#' @param ... additional parameters for stan sampling
#' @export
#' @author Diogo Melo
lmm_animal = function(Y, X, A, iter = 3000, warmup = floor(iter/2), 
                      chains = 4, control = list(adapt_delta = 0.99), 
                      cmd_stan = TRUE, ...){
  if(cmd_stan)
    if(!"cmdstanr" %in% installed.packages()[,"Package"]){
      warning("cmdstanr not installed, falling back to rstan.\nInstall cmdstanr using remotes::install_github(\"stan-dev/cmdstanr)\".\nNext, run cmdstanr::install_cmdstan().")
      cmd_stan = FALSE
    }
      
  p = ncol(Y)
  if(is.null(p)) p = 1
  if(p == 1){
    stan_data = list(J = ncol(X),
                     N = length(Y),
                     X = X,
                     Y = as.numeric(Y),
                     A = as.matrix(A))
    fit = rstan::sampling(stanmodels$animalModelUni, data = stan_data, 
                          iter = iter, warmup = warmup, chains = chains, control = control, ...)
    if(!is.null(colnames(X))){
      names(fit)[grep("beta", names(fit))] = colnames(X)
    }
    fit = new("animalUni", fit)
  } else{
    stan_data = list(K = ncol(Y),
                     J = ncol(X),
                     N = nrow(Y),
                     X = X,
                     Y = Y,
                     A = as.matrix(A))
    if(cmd_stan)
      fit = cstan(model_code = stanmodels$animalModel@model_code, 
                  data = stan_data, 
                  iter = iter, warmup = warmup, chains = chains , threads=1 , 
                  control=control , ...)
    else
      fit = rstan::sampling(stanmodels$animalModel, 
                            data = stan_data, 
                            iter = iter, warmup = warmup, chains = chains, 
                            control = control, ...)
    if(!is.null(colnames(X))){
      names(fit)[grep("beta", names(fit))] = colnames(X)
    }
  }
  return(fit)
}

setClass("animalUni",    contains = "stanfit")
setClass("animalMulti", contains = "stanfit")