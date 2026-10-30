#' Fit linear mixed animal model
#' 
#' Fits a univariate or multivariate linear mixed animal model using a know relatedness matrix
#' @param Y response variable, univariate or multivariate
#' @param X Finite numeric fixed-effects matrix with the same rows as Y and full column rank. Continuous columns are scaled internally; returned beta draws use original X units.
#' @param A Finite symmetric positive-semidefinite relationship matrix matching the rows of Y.
#' @param lkj_prior Parameter of the prior on the correlations. Larger means more shrinkage.
#' @param beta_prior Parameters of the beta prior on h2. Used only by the univariate model; the multivariate model parametrizes G and E directly and has no prior on h2.
#' @param iter number of MCMC iterations.
#' @param warmup number of warmup iterations, must be smaller than iter
#' @param chains number of chains to run.
#' @param control list of stan control arguments.
#' @param cmd_stan Logical, use cmdstanr. If FALSE, uses rstan.
#' @param blend_A If TRUE, replace A with 0.99 * A + 0.01 * I before fitting, making singular A positive definite.
#' @param residual_structure Residual covariance for multivariate fits: "diagonal" or "unstructured" (full covariance between traits). With a diagonal residual, environmental covariances between traits are absorbed into the genetic covariance, so "unstructured" is recommended when estimating G.
#' @param ... additional parameters for stan sampling
#' @return With \code{cmd_stan = TRUE}, a \code{CmdStanMCMC} object from cmdstanr: use \code{fit$draws()} (a \code{posterior} draws array) and \code{fit$summary()}, and pass these, or the fit itself, to bayesplot (e.g. \code{bayesplot::mcmc_trace(fit$draws("h2"))}, \code{bayesplot::nuts_params(fit)}). Fixed effects are \code{beta[j]} (univariate) or \code{beta[trait,j]} for column \code{j} of \code{X}; use \code{\link{animal_draws}} to get draws named after the columns of \code{X}. With \code{cmd_stan = FALSE}, a \code{stanfit} object whose \code{beta} draws are named after the columns of \code{X}.
#' @details The multivariate model estimates the additive genetic covariance matrix \code{G} and the residual covariance \code{E} directly, each as scales times an LKJ-distributed correlation matrix (\code{lkj_prior}), with half-normal priors on the scales of the standardized responses. Breeding values are marginalized using the eigendecomposition of \code{A}, so the likelihood needs no per-individual random effects and only requires \code{A} to be positive semidefinite. Returned quantities (original response scale) are \code{G}, \code{E}, \code{P = G + E}, \code{corrG}, \code{sigma_G}, \code{sigma_E} and \code{h2 = diag(G) / diag(P)}.
#'
#' Continuous predictors are divided by their sample standard deviation inside Stan. Constant and 0/1 indicator columns are not scaled. The returned \code{beta} coefficients are backtransformed to original predictor and response units; \code{beta_scaled} draws are also available. The response is centered and standardized in all models (including univariate), with a standard normal prior on \code{beta_scaled} and a free total variance on the standardized scale.
#' @export
#' @author Diogo Melo
lmm_animal = function(Y, X, A, lkj_prior = 2, beta_prior = c(2, 4), 
                      iter = 3000, warmup = floor(iter/2), 
                      chains = 4, cores = chains, model_file = NULL, control = list(adapt_delta = 0.99), 
                      cmd_stan = TRUE, blend_A = FALSE,
                      residual_structure = c("diagonal", "unstructured"), ...){
  residual_structure <- match.arg(residual_structure)
  if (is.null(dim(Y))) Y <- matrix(Y, ncol = 1L)
  Y <- as.matrix(Y)
  X <- as.matrix(X)
  A <- as.matrix(A)

  if (!is.numeric(Y) || any(!is.finite(Y)))
    stop("Y must contain only finite numeric values.")
  if (!is.numeric(X) || any(!is.finite(X)))
    stop("X must contain only finite numeric values.")
  if (!is.logical(blend_A) || length(blend_A) != 1L || is.na(blend_A))
    stop("blend_A must be TRUE or FALSE.")

  n <- nrow(Y)
  p <- ncol(X)
  if (n < 2L || p < 1L || nrow(X) != n)
    stop("Y and X must have the same number of rows, and X must have columns.")
  if (qr(X)$rank < p)
    stop("X must have full column rank.")
  if (!is.numeric(A) || any(!is.finite(A)) || !identical(dim(A), c(n, n)))
    stop("A must be a finite numeric square matrix matching the rows of Y.")

  a_scale <- max(1, max(abs(A)))
  if (max(abs(A - t(A))) > 1e-8 * a_scale)
    stop("A must be symmetric.")
  A <- (A + t(A)) / 2
  eigenvalues <- eigen(A, symmetric = TRUE, only.values = TRUE)$values
  eigen_tol <- 1e-8 * max(1, max(abs(eigenvalues)))
  if (min(eigenvalues) < -eigen_tol)
    stop("A must be positive semidefinite.")

  if (blend_A) {
    A <- 0.99 * A + 0.01 * diag(n)
    eigenvalues <- eigen(A, symmetric = TRUE, only.values = TRUE)$values
    eigen_tol <- 1e-8 * max(1, max(abs(eigenvalues)))
  }
  if (min(eigenvalues) <= eigen_tol)
    stop("A must be positive definite for Stan; set blend_A = TRUE to regularize it.")

  if(cmd_stan)
    if(!requireNamespace("cmdstanr", quietly = TRUE)){
      warning("cmdstanr not installed, falling back to rstan.\nInstall cmdstanr using remotes::install_github(\"stan-dev/cmdstanr)\".\nNext, run cmdstanr::install_cmdstan().")
      cmd_stan = FALSE
    }
      
  if(ncol(Y) == 1L){
    stan_data = list(J = ncol(X),
                     N = nrow(Y),
                     X = X,
                     Y = as.numeric(Y),
                     A = A,
                     lkj_prior =  lkj_prior, 
                     beta_prior= beta_prior)
    if(cmd_stan){
      model_code = if(is.null(model_file))
        stanmodels$animalModelUni@model_code else
          .read_animal_model_code(model_file, "animalModelUni.stan")
      fit = cstan(model_code = model_code,
                  data = stan_data,
                  iter = iter, warmup = warmup, chains = chains, cores = cores,
                  control = control, ...)
      attr(fit, "x_names") <- colnames(X)
    }
    else {
      model_code <- if (is.null(model_file))
        stanmodels$animalModelUni@model_code else
          .read_animal_model_code(model_file, "animalModelUni.stan")
      model <- .get_rstan_model(
        model_code,
        "animalModelUni_runtime",
        if (is.null(model_file)) "animalModelUni" else NULL
      )
      fit = rstan::sampling(model, data = stan_data,
                            iter = iter, warmup = warmup, chains = chains, control = control, ...)
    }
    if(!cmd_stan)
      fit <- .rename_fixed_effect_draws(fit, X, ncol(Y))
  } else{
    stan_data = list(K = ncol(Y),
                     J = ncol(X),
                     N = nrow(Y),
                     X = X,
                     Y = Y,
                     A = A,
                     lkj_prior = lkj_prior,
                     unstructured_E = as.integer(residual_structure == "unstructured"))
    model_code <- .read_animal_model_code(model_file)
    if(cmd_stan){
      fit = cstan(model_code = model_code,
                  data = stan_data,
                  iter = iter, warmup = warmup, chains = chains , cores = cores,
                  control=control , ...)
      attr(fit, "x_names") <- colnames(X)
    } else {
      model <- .get_rstan_model(
        model_code,
        "animalModelG_runtime",
        if (is.null(model_file)) "animalModelG" else NULL
      )
      fit = rstan::sampling(model,
                            data = stan_data,
                            iter = iter, warmup = warmup, chains = chains,
                            control = control, ...)
      fit <- .rename_fixed_effect_draws(fit, X, ncol(Y))
    }
  }
  return(fit)
}

.beta_names <- function(x_names, n_traits) {
  if (n_traits == 1L)
    return(list(old = paste0("beta[", seq_along(x_names), "]"), new = paste0("beta[", x_names, "]")))
  grid <- expand.grid(trait = seq_len(n_traits), j = seq_along(x_names))
  list(old = paste0("beta[", grid$trait, ",", grid$j, "]"),
       new = paste0("beta[", grid$trait, ",", x_names[grid$j], "]"))
}

.rename_fixed_effect_draws <- function(fit, X, n_traits) {
  if (is.null(colnames(X)))
    return(fit)

  beta_indices <- grep("^beta\\[", names(fit))
  beta_names <- .beta_names(colnames(X), n_traits)$new
  if (length(beta_indices) == length(beta_names))
    names(fit)[beta_indices] <- beta_names
  fit
}

#' Posterior draws with named fixed effects
#'
#' Returns draws from a \code{lmm_animal} fit made with \code{cmd_stan = TRUE},
#' with the \code{beta} variables named after the columns of \code{X}
#' (\code{beta[name]} for one trait, \code{beta[trait,name]} for several).
#' If \code{X} had no column names the draws are returned unchanged.
#' @param fit A \code{CmdStanMCMC} object returned by \code{lmm_animal}.
#' @param variables Optional character vector of variables to extract, using
#'   Stan's names (e.g. \code{"beta"}, \code{"h2"}). Default is all variables.
#' @param ... Further arguments passed to \code{fit$draws()}, e.g. \code{format}.
#' @return A \code{posterior} draws object, usable with bayesplot.
#' @export
animal_draws <- function(fit, variables = NULL, ...) {
  draws <- fit$draws(variables = variables, ...)
  x_names <- attr(fit, "x_names")
  if (is.null(x_names))
    return(draws)
  vars <- posterior::variables(draws)
  n_traits <- if (any(grepl("^beta\\[[0-9]+,[0-9]+\\]$", vars))) {
    max(as.integer(sub("^beta\\[([0-9]+),.*", "\\1", grep("^beta\\[[0-9]+,", vars, value = TRUE))))
  } else 1L
  map <- .beta_names(x_names, n_traits)
  hit <- match(vars, map$old)
  vars[!is.na(hit)] <- map$new[hit[!is.na(hit)]]
  posterior::variables(draws) <- vars
  draws
}

.read_animal_model_code <- function(model_file, default = "animalModelG.stan") {
  if (!is.null(model_file))
    return(readChar(model_file, file.info(model_file)$size))

  model_filename <- default
  model_path <- system.file("stan", model_filename, package = "stanAnimal")
  if (!nzchar(model_path)) {
    candidates <- file.path(c("inst/stan", "package/inst/stan"), model_filename)
    model_path <- candidates[file.exists(candidates)][1L]
  }
  if (is.na(model_path) || !file.exists(model_path))
    stop("Could not locate the requested Stan model source.")
  paste(readLines(model_path, warn = FALSE), collapse = "\n")
}

.rstan_model_cache <- new.env(parent = emptyenv())

.get_rstan_model <- function(model_code, model_name, cache_key = NULL) {
  if (!is.null(cache_key) && exists(cache_key, envir = .rstan_model_cache,
                                    inherits = FALSE)) {
    cached <- get(cache_key, envir = .rstan_model_cache, inherits = FALSE)
    if (identical(cached$code, model_code))
      return(cached$model)
    model_name <- basename(tempfile(paste0(model_name, "_")))
  } else if (is.null(cache_key)) {
    model_name <- basename(tempfile(paste0(model_name, "_")))
  }

  model <- rstan::stan_model(model_code = model_code, model_name = model_name)
  if (!is.null(cache_key))
    assign(cache_key, list(code = model_code, model = model),
           envir = .rstan_model_cache)
  model
}
