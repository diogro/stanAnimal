#' Fit a Bayesian phenotypic covariance matrix
#'
#' Estimates the phenotypic covariance matrix \code{P} of one or more traits
#' after accounting for fixed effects, for data with no relatedness or genomic
#' information. It is the counterpart of \code{\link{lmm_animal}} without the
#' genetic component, and uses the same priors, scaling and outputs.
#' @param Y Numeric response vector or matrix, with individuals in rows.
#' @param X Finite numeric fixed-effects matrix with the same rows as Y and full column rank. Continuous columns are scaled internally; returned beta draws use original X units.
#' @param lkj_prior Parameter of the LKJ prior on the correlation matrix. Higher values result in lower correlations, and 1 is uniform. Default is 2.
#' @param iter number of MCMC iterations.
#' @param warmup number of warmup iterations, must be smaller than iter.
#' @param chains number of chains to run.
#' @param cores number of chains to run in parallel (cmdstanr only).
#' @param control list of stan control arguments.
#' @param cmd_stan Logical, use cmdstanr. If FALSE, uses rstan.
#' @param model_file Optional path to a custom Stan file replacing the default model.
#' @param ... additional parameters for stan sampling.
#' @return With \code{cmd_stan = TRUE}, a \code{CmdStanMCMC} object from cmdstanr, usable with posterior and bayesplot (see \code{\link{lmm_animal}}); use \code{\link{animal_draws}} for \code{beta} draws named after the columns of \code{X}. With \code{cmd_stan = FALSE}, a \code{stanfit} object.
#' @details The model estimates \code{P} as scales times an LKJ-distributed correlation matrix, with half-normal priors on the scales of the standardized responses. Returned quantities (original response scale) are \code{P}, \code{corrP}, \code{sigma_P} and the fixed effects \code{beta[trait,j]}.
#'
#' Continuous predictors are divided by their sample standard deviation inside Stan. Constant and 0/1 indicator columns are not scaled. The returned \code{beta} coefficients are backtransformed to original predictor and response units; \code{beta_scaled} draws are also available. The response is centered and standardized.
#' @author Diogo Melo
#' @export
lmm_multiCov = function(Y, X, lkj_prior = 2, iter = 3000, warmup = floor(iter/2),
                        chains = 4, cores = chains, control = list(adapt_delta = 0.99),
                        cmd_stan = TRUE, model_file = NULL, ...){
  YX <- .check_response_design(Y, X)
  Y <- YX$Y
  X <- YX$X

  if(cmd_stan && !requireNamespace("cmdstanr", quietly = TRUE)){
    warning("cmdstanr not installed, falling back to rstan.\nInstall cmdstanr using remotes::install_github(\"stan-dev/cmdstanr)\".\nNext, run cmdstanr::install_cmdstan().")
    cmd_stan = FALSE
  }

  stan_data = list(K = ncol(Y),
                   J = ncol(X),
                   N = nrow(Y),
                   X = X,
                   Y = Y,
                   lkj_prior = lkj_prior)
  model_code <- .read_animal_model_code(model_file, "phenotypicModel.stan")

  if(cmd_stan){
    fit = cstan(model_code = model_code,
                data = stan_data,
                iter = iter, warmup = warmup, chains = chains, cores = cores,
                control = control, ...)
    attr(fit, "x_names") <- colnames(X)
  } else {
    model <- .get_rstan_model(
      model_code,
      "phenotypicModel_runtime",
      if (is.null(model_file)) "phenotypicModel" else NULL
    )
    fit = rstan::sampling(model, data = stan_data,
                          iter = iter, warmup = warmup, chains = chains,
                          control = control, ...)
    fit <- .rename_fixed_effect_draws(fit, X, ncol(Y))
  }
  return(fit)
}
