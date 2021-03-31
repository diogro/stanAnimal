#' Fit linear mixed animal model
#' 
#' Fits a univariate or multivariate linear mixed animal model using a know relatedness matrix
#' @param Y response variable, univariate or multivariate
#' @param X fixed model matrix.
#' @param A relatedness matrix.
#' @param lkj_prior Parameter of the prior on the correlations. Larger means more shrinkage.
#' @param beta_prior Parameter of the beta prior on the h2. 
#' @param iter number of MCMC iterations.
#' @param warmup number of warmup iterations, must be smaller than iter
#' @param chains number of chains to run.
#' @param control list of stan control arguments.
#' @param cmd_stan Logical, use cmd_stan. If false, uses rstan.
#' @param ... additional parameters for stan sampling
#' @export
#' @author Diogo Melo
lmm_animal = function(Y, X, A, lkj_prior = 2, beta_prior = c(2, 4), 
                      iter = 3000, warmup = floor(iter/2), 
                      chains = 4, cores = chains, model_file = NULL, control = list(adapt_delta = 0.99), 
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
                     A = as.matrix(A), 
                     lkj_prior =  lkj_prior, 
                     beta_prior= beta_prior)
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
                     A = as.matrix(A), 
                     lkj_prior =  lkj_prior,
                     beta_prior= beta_prior)
    if(cmd_stan){
      if(is.null(model_file))
        model_code = stanmodels$animalModel@model_code
      else
        model_code = readChar(model_file, 1e5)
      fit = cstan(model_code = model_code, 
                  data = stan_data, 
                  iter = iter, warmup = warmup, chains = chains , cores = cores, 
                  control=control , ...)
    }
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