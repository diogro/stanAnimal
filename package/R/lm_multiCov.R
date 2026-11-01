#' Fit a Bayesian phenotypic covariance matrix
#'
#' Estimates the phenotypic covariance matrix \code{P} of one or more traits
#' after accounting for fixed effects, for data with no relatedness or genomic
#' information. It is the counterpart of \code{\link{lmm_animal}} without the
#' genetic component, and uses the same priors, scaling and outputs.
#' @param Y Numeric response vector or matrix, with individuals in rows. Every trait must vary: constant traits cannot be standardized and are an error.
#' @param X Finite numeric fixed-effects matrix with the same rows as Y and full column rank. The columns must be able to represent a constant (an intercept column, or the 0/1 indicators of every level of a factor), so that the mean of the response can be reported. Continuous columns are centered and scaled internally; returned beta draws use original X units.
#' @param lkj_prior Shape parameter of the LKJ prior on the correlation matrix. \code{1} is uniform; values above \code{1} shrink correlations toward zero and values below \code{1} favor larger correlations. The default \code{2} is recommended for many traits, because at high dimensionality the uniform prior concentrates mass near \code{+/-1} and overestimates the magnitude of correlations.
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
#' The likelihood is computed by whitening the residuals with the Cholesky factor of \code{P} in a single triangular solve, at \code{O(N K^2 + K^3)} per evaluation rather than a per-individual loop, so the model scales to dozens of traits.
#'
#' Continuous predictors are centered and divided by their sample standard deviation inside Stan. Constant and 0/1 indicator columns are not centered or scaled. Standardized slopes have \code{normal(0, 1)} priors; the intercept contribution has a wider \code{normal(0, 2.5)} prior. Returned \code{beta} coefficients are backtransformed to original predictor and response units, including restoring the response mean to the intercept; \code{beta_scaled} draws are also available.
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
                   lkj_prior = lkj_prior,
                   x_mean_weight = YX$x_mean_weight)
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
    fit = .rstan_sample(model, stan_data, iter, warmup, chains, control, ...)
    fit <- .rename_fixed_effect_draws(fit, X, ncol(Y))
  }
  return(fit)
}
