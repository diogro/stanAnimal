#' Fit linear mixed animal model
#' 
#' Fits a univariate or multivariate linear mixed animal model using a know relatedness matrix
#' @param Y response variable, univariate or multivariate
#' @param X fixed model matrix
#' @param A relatedness matrix
#' @param iter number of MCMC iterations
#' @param warmup number of warmup iterations, must be smaller than iter
#' @param chains number of chains to run
#' @param control list of stan control arguments
#' @param ... additional parameters for stan sampling
#' @author Diogo Melo
#' @export
lmm_animal = function(Y, X, A, iter = 3000, warmup = 2000, 
                      chains = 4, control = list(adapt_delta = 0.99), ...){
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
    fit = rstan::sampling(stanmodels$animalModel, data = stan_data, 
                          iter = iter, warmup = warmup, chains = chains, control = control, ...)
    if(!is.null(colnames(X))){
      names(fit)[grep("beta", names(fit))] = colnames(X)
    }
  }
  return(fit)
}

setClass("animalUni",    contains = "stanfit")
setClass("animalMulti", contains = "stanfit")