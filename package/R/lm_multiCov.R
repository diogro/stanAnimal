#' Fit a bayesian covariance matrix
#' 
#' Fits a multivariate covariance matrix for the residuals a linear model. 
#' @param Y response variable, univariate or multivariate
#' @param X Fixed-effects model matrix. Continuous columns are scaled internally; returned beta draws use original X units.
#' @param LKJ parameter of the LKJ prior on the correlations. Higher values result in lower correlations. A LKJ=1 implies uniform prior on the correlations. Default is 4, which tends to give reasonable amount of shrinkage.
#' @param iter number of MCMC iterations
#' @param warmup number of warmup iterations, must be smaller than iter
#' @param chains number of chains to run
#' @param control list of stan control arguments
#' @param ... additional parameters for Stan sampling
#' @details Continuous predictors are divided by their sample standard deviation inside Stan. Constant and 0/1 indicator columns are not scaled. Returned \code{beta} coefficients are backtransformed to original predictor and response units; \code{beta_scaled} draws are also available.
#' @author Diogo Melo
#' @export
lmm_multiCov = function(Y, X, LKJ = 4, iter = 3000, warmup = 2000, 
                      chains = 4, control = list(adapt_delta = 0.99), ...){
  p = ncol(Y)
  if(is.null(p)) p = 1
  if(p == 1){
    stop("This function is for multivariate responses")
  } else{
    stan_data = list(K = ncol(Y),
                     J = ncol(X),
                     N = nrow(Y),
                     X = X,
                     Y = Y,
                     c = LKJ)
    model <- .get_rstan_model(stanmodels$multiNormalCov@model_code,
                  "multiNormalCov_runtime", "multiNormalCov")
    fit = rstan::sampling(model, data = stan_data,
                          iter = iter, warmup = warmup, chains = chains, control = control, ...)
    fit <- .rename_fixed_effect_draws(fit, X, ncol(Y))
  }
  return(fit)
}
