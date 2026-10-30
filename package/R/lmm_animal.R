#' Fit linear mixed animal model
#' 
#' Fits a univariate or multivariate linear mixed animal model using a know relatedness matrix
#' @param Y response variable, univariate or multivariate
#' @param X Finite numeric fixed-effects matrix with the same rows as Y and full column rank. Continuous columns are scaled internally; returned beta draws use original X units.
#' @param A Finite symmetric positive-semidefinite relationship matrix matching the rows of Y.
#' @param lkj_prior Parameter of the prior on the correlations. Larger means more shrinkage.
#' @param beta_prior Parameter of the beta prior on the h2. 
#' @param iter number of MCMC iterations.
#' @param warmup number of warmup iterations, must be smaller than iter
#' @param chains number of chains to run.
#' @param control list of stan control arguments.
#' @param cmd_stan Logical, use cmd_stan. If false, uses rstan.
#' @param blend_A If TRUE, replace A with 0.99 * A + 0.01 * I before fitting, making singular A positive definite.
#' @param residual_structure Residual covariance for multivariate fits: "diagonal" or "unstructured".
#' @param ... additional parameters for stan sampling
#' @details Continuous predictors are divided by their sample standard deviation inside Stan. Constant and 0/1 indicator columns are not scaled. The returned \code{beta} coefficients are backtransformed to original predictor and response units; \code{beta_scaled} draws are also available. The univariate model uses the response SD to set the scaled coefficient prior.
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
    if(!"cmdstanr" %in% installed.packages()[,"Package"]){
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
      if(is.null(model_file))
        model_code = stanmodels$animalModelUni@model_code
      else
        model_code = readChar(model_file, 1e5)
      fit = cstan(model_code = model_code,
                  data = stan_data,
                  iter = iter, warmup = warmup, chains = chains, cores = cores,
                  control = control, ...)
    }
    else {
      model_code <- if (is.null(model_file))
        stanmodels$animalModelUni@model_code else
          .read_animal_model_code(model_file, "diagonal")
      model <- .get_rstan_model(
        model_code,
        "animalModelUni_runtime",
        if (is.null(model_file)) "animalModelUni" else NULL
      )
      fit = rstan::sampling(model, data = stan_data,
                            iter = iter, warmup = warmup, chains = chains, control = control, ...)
    }
    fit <- .rename_fixed_effect_draws(fit, X, ncol(Y))
    fit = new("animalUni", fit)
  } else{
    stan_data = list(K = ncol(Y),
                     J = ncol(X),
                     N = nrow(Y),
                     X = X,
                     Y = Y,
                     A = A,
                     lkj_prior =  lkj_prior,
                     beta_prior= beta_prior)
    if(cmd_stan){
      if(is.null(model_file) && residual_structure == "diagonal")
        model_code = stanmodels$animalModel@model_code
      else
        model_code = .read_animal_model_code(model_file, residual_structure)
      fit = cstan(model_code = model_code, 
                  data = stan_data, 
                  iter = iter, warmup = warmup, chains = chains , cores = cores, 
                  control=control , ...)
    } else {
      model_code <- if (is.null(model_file) && residual_structure == "diagonal")
        stanmodels$animalModel@model_code else
          .read_animal_model_code(model_file, residual_structure)
      cache_key <- if (is.null(model_file))
        paste0("animalModel_", residual_structure) else NULL
      model <- .get_rstan_model(
        model_code,
        paste0("animalModel_", residual_structure, "_runtime"),
        cache_key
      )
      fit = rstan::sampling(model,
                            data = stan_data, 
                            iter = iter, warmup = warmup, chains = chains, 
                            control = control, ...)
    }
    fit <- .rename_fixed_effect_draws(fit, X, ncol(Y))
  }
  return(fit)
}

.rename_fixed_effect_draws <- function(fit, X, n_traits) {
  if (is.null(colnames(X)))
    return(fit)

  beta_indices <- grep("^beta\\[", names(fit))
  if (n_traits == 1L) {
    beta_names <- colnames(X)
  } else {
    beta_names <- unlist(lapply(seq_len(ncol(X)), function(j) {
      paste0("beta[", seq_len(n_traits), ",", colnames(X)[j], "]")
    }), use.names = FALSE)
  }
  if (length(beta_indices) == length(beta_names))
    names(fit)[beta_indices] <- beta_names
  fit
}

.read_animal_model_code <- function(model_file, residual_structure) {
  if (!is.null(model_file))
    return(readChar(model_file, file.info(model_file)$size))

  model_filename <- if (residual_structure == "unstructured")
    "animalModelCorR.stan" else "animalModel.stan"
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

setClass("animalUni",    contains = "stanfit")
setClass("animalMulti", contains = "stanfit")