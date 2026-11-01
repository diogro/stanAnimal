# Wrapped so that the tests can pretend a package is not installed.
.is_installed <- function(package) requireNamespace(package, quietly = TRUE)

# Which Stan engine to use. Both cmdstanr and rstan are suggested packages, so
# a fit falls back to rstan when cmdstanr is not installed, and errors with the
# installation instructions when neither is available.
.resolve_stan_backend <- function(cmd_stan) {
  if (cmd_stan) {
    if (.is_installed("cmdstanr")) return(TRUE)
    if (!.is_installed("rstan"))
      stop("fitting a model needs cmdstanr or rstan: install one of them with ",
           "remotes::install_github(\"stan-dev/cmdstanr\") (then ",
           "cmdstanr::install_cmdstan()) or install.packages(\"rstan\").",
           call. = FALSE)
    warning("cmdstanr is not installed, falling back to rstan. Install it with ",
            "remotes::install_github(\"stan-dev/cmdstanr\"), then run ",
            "cmdstanr::install_cmdstan().", call. = FALSE)
    return(FALSE)
  }
  if (!.is_installed("rstan"))
    stop("cmd_stan = FALSE needs rstan: install.packages(\"rstan\").",
         call. = FALSE)
  FALSE
}

#' Fit linear mixed animal model
#' 
#' Fits a univariate or multivariate linear mixed animal model using a known relatedness matrix
#' @param Y response variable, univariate or multivariate. Every trait must vary: constant traits cannot be standardized and are an error.
#' @param X Finite numeric fixed-effects matrix with the same rows as Y and full column rank. The columns must be able to represent a constant (an intercept column, or the 0/1 indicators of every level of a factor), so that the mean of the response can be reported. Continuous columns are centered and scaled internally; returned beta draws use original X units.
#' @param A Finite symmetric positive-definite relationship matrix matching the rows of Y. A matrix that is only positive semidefinite (duplicated individuals, for example) is rejected; set \code{blend_A = TRUE} to regularize it. When \code{A} and \code{Y} both have row names, they are used to put \code{A} in the order of \code{Y}; names that disagree are an error. Without names the rows are matched by position.
#' @param lkj_prior Shape parameter of the LKJ prior on the correlation matrices of the multivariate model. \code{1} is uniform over correlation matrices, values above \code{1} shrink correlations toward zero, and values below \code{1} favor larger correlations. The default \code{2} is recommended for many traits: at high dimensionality the uniform prior concentrates mass near \code{+/-1} and overestimates the magnitude of genetic correlations, a well-known problem in quantitative genetics.
#' @param beta_prior Parameters of the beta prior on h2, used only by the univariate model with \code{parametrization = "h2"}. The default \code{c(1, 1)} is flat on h2.
#' @param parametrization Univariate models only. \code{"variance"} (default) estimates the genetic and residual standard deviations with half-normal priors and derives h2. \code{"h2"} estimates a total standard deviation and h2 with a beta prior (\code{beta_prior}). If the variance version shows low E-BFMI, a warning suggests switching to \code{"h2"}. Ignored for multivariate fits.
#' @param iter number of MCMC iterations.
#' @param warmup number of warmup iterations, must be smaller than iter
#' @param chains number of chains to run.
#' @param cores number of chains to run in parallel (cmdstanr only).
#' @param model_file Optional path to a custom Stan file replacing the default model.
#' @param control Named list of sampler controls, passed to whichever engine runs: \code{adapt_delta}, \code{adapt_engaged}, \code{max_treedepth}, \code{step_size} (rstan's \code{stepsize} is accepted), \code{init_buffer}, \code{term_buffer}, \code{window} and \code{metric}. Any other name is an error, so a control is never silently ignored.
#' @param cmd_stan Logical, use cmdstanr. If FALSE, uses rstan. Both engines are suggested packages: when cmdstanr is not installed the fit falls back to rstan with a warning, and when neither is installed the fit errors.
#' @param blend_A If TRUE, replace A with 0.99 * A + 0.01 * I before fitting, making singular A positive definite.
#' @param residual_structure Residual covariance for multivariate fits: "unstructured" (default, full covariance between traits) or "diagonal". With a diagonal residual, environmental covariances between traits are absorbed into the genetic covariance and inflate heritability, so "unstructured" is recommended when estimating G.
#' @param ... additional parameters for stan sampling
#' @return With \code{cmd_stan = TRUE}, a \code{CmdStanMCMC} object from cmdstanr: use \code{fit$draws()} (a \code{posterior} draws array) and \code{fit$summary()}, and pass these, or the fit itself, to bayesplot (e.g. \code{bayesplot::mcmc_trace(fit$draws("h2"))}, \code{bayesplot::nuts_params(fit)}). Fixed effects are \code{beta[j]} (univariate) or \code{beta[trait,j]} for column \code{j} of \code{X}; use \code{\link{animal_draws}} to get draws named after the columns of \code{X}. With \code{cmd_stan = FALSE}, a \code{stanfit} object whose \code{beta} draws are named after the columns of \code{X}.
#' @details The models estimate the additive genetic and residual standard deviations (\code{sigma_G}, \code{sigma_E}) directly, with half-normal priors on the standardized response scale, and report \code{h2} as a derived quantity; there is no prior on h2 itself. The multivariate model estimates the additive genetic covariance matrix \code{G} and the residual covariance \code{E} directly, each as scales times an LKJ-distributed correlation matrix (\code{lkj_prior}), with half-normal priors on the scales of the standardized responses. Breeding values are marginalized using the eigendecomposition of \code{A}, so the likelihood needs no per-individual random effects; only the eigenvalues of \code{A} are used, which is why a singular \code{A} can be fitted once it is regularized with \code{blend_A = TRUE}. Returned quantities (original response scale) are \code{G}, \code{E}, \code{P = G + E}, \code{corrG}, \code{sigma_G}, \code{sigma_E} and \code{h2 = diag(G) / diag(P)}.
#'
#' The per-individual likelihoods are whitened by the residual Cholesky factor and rotated into a single eigendecomposition of \code{E^{-1/2} G E^{-1/2}}, so each evaluation costs \code{O(N K^2 + K^3)} rather than \code{O(N K^3)}. This lets the multivariate model scale to dozens of traits without a Cholesky factorization per individual. When the eigenvalues of that matrix are nearly coincident (for example at initialization, or when \code{G} is roughly proportional to \code{E}) the model automatically falls back to a per-individual Cholesky evaluation, which has a well-behaved gradient.
#'
#' The cost that grows with the number of individuals is the dense relationship matrix: \code{A} is passed in full and each chain eigendecomposes it before sampling (\code{O(N^3)} time and \code{O(N^2)} memory per chain). That cost is paid once, so with many traits sampling dominates and the decomposition is a small share of the fit, while with few traits and short runs it can be most of it. The models are intended for many traits on hundreds to a few thousand individuals. For cohorts of tens of thousands, where that matrix is the bottleneck, an engine that uses the sparse inverse of \code{A} (such as ASReml) is the better tool.
#'
#' The LKJ shape \code{lkj_prior} (default 2) deliberately shrinks correlations toward zero. With many traits, an unstructured correlation matrix fit under a flat prior tends to overestimate the magnitude of correlations, especially genetic correlations between traits with limited information, so a shape above 1 acts as a mild regularizer; increase it for stronger shrinkage. Note that fitting many traits remains expensive: the likelihood is cheap, but an unstructured \code{G} for \code{K = 40} traits has on the order of 1700 parameters and NUTS needs many leapfrog steps, so a single fit can take an hour or more. Fewer iterations, a larger \code{lkj_prior}, or \code{residual_structure = "diagonal"} all make fits faster.
#'
#' Continuous predictors are centered and divided by their sample standard deviation inside Stan. Constant and 0/1 indicator columns are not centered or scaled. Standardized slopes have \code{normal(0, 1)} priors; the intercept contribution has a wider \code{normal(0, 2.5)} prior. Returned \code{beta} coefficients are backtransformed to original predictor and response units, including restoring the response mean to the intercept; \code{beta_scaled} draws are also available.
#' @export
#' @author Diogo Melo
lmm_animal = function(Y, X, A, lkj_prior = 2, beta_prior = c(1, 1),
                      parametrization = c("variance", "h2"),
                      iter = 3000, warmup = floor(iter/2), 
                      chains = 4, cores = chains, model_file = NULL, control = list(adapt_delta = 0.99), 
                      cmd_stan = TRUE, blend_A = FALSE,
                      residual_structure = c("unstructured", "diagonal"), ...){
  residual_structure <- match.arg(residual_structure)
  parametrization <- match.arg(parametrization)
  YX <- .check_response_design(Y, X)
  Y <- YX$Y
  X <- YX$X
  A <- as.matrix(A)

  if (!is.logical(blend_A) || length(blend_A) != 1L || is.na(blend_A))
    stop("blend_A must be TRUE or FALSE.")

  n <- nrow(Y)
  A <- .check_relationship_matrix(A, Y)

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

  cmd_stan <- .resolve_stan_backend(cmd_stan)
      
  if(ncol(Y) == 1L){
    stan_data = list(J = ncol(X),
                     N = nrow(Y),
                     X = X,
                     Y = as.numeric(Y),
                     A = A,
                     x_mean_weight = YX$x_mean_weight)
    uni_file <- if(parametrization == "h2") "animalModelUniH2.stan" else "animalModelUni.stan"
    if(parametrization == "h2")
      stan_data$beta_prior <- beta_prior
    model_code <- .read_animal_model_code(model_file, uni_file)
    if(cmd_stan){
      fit = cstan(model_code = model_code,
                  data = stan_data,
                  iter = iter, warmup = warmup, chains = chains, cores = cores,
                  control = control, ...)
      attr(fit, "x_names") <- colnames(X)
    }
    else {
      model <- .get_rstan_model(
        model_code,
        paste0(sub("\\.stan$", "", uni_file), "_runtime"),
        if (is.null(model_file)) sub("\\.stan$", "", uni_file) else NULL
      )
      fit = .rstan_sample(model, stan_data, iter, warmup, chains, control, ...)
      fit <- .rename_fixed_effect_draws(fit, X, ncol(Y))
    }
    if(parametrization == "variance")
      .warn_low_bfmi(fit)
  } else{
    stan_data = list(K = ncol(Y),
                     J = ncol(X),
                     N = nrow(Y),
                     X = X,
                     Y = Y,
                     A = A,
                     lkj_prior = lkj_prior,
                     x_mean_weight = YX$x_mean_weight,
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
      fit = .rstan_sample(model, stan_data, iter, warmup, chains, control, ...)
      fit <- .rename_fixed_effect_draws(fit, X, ncol(Y))
    }
  }
  return(fit)
}

.check_response_design <- function(Y, X) {
  if (is.null(dim(Y))) Y <- matrix(Y, ncol = 1L)
  Y <- as.matrix(Y)
  X <- as.matrix(X)
  if (!is.numeric(Y) || any(!is.finite(Y)))
    stop("Y must contain only finite numeric values.")
  if (!is.numeric(X) || any(!is.finite(X)))
    stop("X must contain only finite numeric values.")
  p <- ncol(X)
  if (nrow(Y) < 2L || p < 1L || nrow(X) != nrow(Y))
    stop("Y and X must have the same number of rows, and X must have columns.")
  if (qr(X)$rank < p)
    stop("X must have full column rank.")
  # The models standardize the responses by their standard deviation, so a
  # constant trait would divide by zero inside Stan
  spread <- apply(Y, 2L, stats::sd)
  constant <- spread <= 0
  if (any(constant)) {
    labels <- colnames(Y) %||% rep("", ncol(Y))
    blank <- is.na(labels) | !nzchar(labels)
    labels[blank] <- paste("column", which(blank))
    stop("every trait in Y must vary, so that it can be standardized: ",
         paste(labels[constant], collapse = ", "), call. = FALSE)
  }
  weights <- .intercept_weights(X)
  if (is.null(weights))
    stop("X cannot represent a constant, so the mean of the response cannot be ",
         "reported: include an intercept column, or the 0/1 indicator columns ",
         "of every level of a factor.", call. = FALSE)
  list(Y = Y, X = X, x_mean_weight = weights)
}

# Weights w with X w equal to a column of ones, when the constant is in the
# column space of X, and NULL when it is not. The models centre the responses,
# so the response mean is reported back through these weights: an intercept
# column gets all of it, the indicators of a factor share it.
.intercept_weights <- function(X) {
  weights <- qr.solve(X, rep(1, nrow(X)))
  if (anyNA(weights)) return(NULL)
  fitted <- drop(X %*% weights)
  if (max(abs(fitted - 1)) > 1e-8 * max(1, max(abs(fitted)))) return(NULL)
  weights
}

# Checks that A matches the rows of Y and returns it with its individuals in the
# order of Y. When both A and Y are named the names are used to align A; names
# that disagree are an error. Without names the rows are matched by position, so
# only the dimensions are checked.
.check_relationship_matrix <- function(A, Y) {
  Y <- as.matrix(Y)
  n <- nrow(Y)
  if (!is.numeric(A) || any(!is.finite(A)) || is.null(dim(A)))
    stop("A must be a finite numeric square matrix matching the rows of Y.",
         call. = FALSE)
  if (!identical(dim(A), c(n, n)))
    stop("A must be a square matrix matching the rows of Y: A is ", nrow(A),
         " x ", ncol(A), " but Y has ", n, " rows.", call. = FALSE)

  y_ids <- rownames(Y)
  a_ids <- rownames(A)
  if (is.null(y_ids) || is.null(a_ids))
    return(A)
  if (is.null(colnames(A)) || !identical(a_ids, colnames(A)))
    stop("A must have identical row and column names.", call. = FALSE)
  if (anyDuplicated(a_ids))
    stop("A has duplicated individual names.", call. = FALSE)

  missing_in_A <- setdiff(y_ids, a_ids)
  missing_in_Y <- setdiff(a_ids, y_ids)
  if (length(missing_in_A) || length(missing_in_Y))
    stop("the individuals in A do not match the rows of Y: ",
         length(missing_in_A), " rows of Y are not in A and ",
         length(missing_in_Y), " individuals in A are not in Y. Rename or ",
         "reorder them so the names agree, or drop the names from both to ",
         "match by position.", call. = FALSE)

  if (identical(y_ids, a_ids))
    return(A)
  A[y_ids, y_ids, drop = FALSE]
}

.warn_low_bfmi <- function(fit, threshold = 0.3) {
  bfmi <- tryCatch(
    if (inherits(fit, "stanfit")) rstan::get_bfmi(fit) else
      fit$diagnostic_summary(quiet = TRUE)$ebfmi,
    error = function(e) NULL)
  if (length(bfmi) && any(bfmi < threshold, na.rm = TRUE))
    warning(sum(bfmi < threshold, na.rm = TRUE), " of ", length(bfmi),
            " chains had E-BFMI below ", threshold, ". The variance ",
            "parametrization can sample poorly when the genetic variance is ",
            "weakly identified; consider parametrization = \"h2\".",
            call. = FALSE)
  invisible(fit)
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
  if (n_traits == 1L && any(grepl("^beta\\[1,[0-9]+\\]$", vars)))
    map$old <- paste0("beta[1,", seq_along(x_names), "]")
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

.rstan_sample <- function(model, stan_data, iter, warmup, chains, control, ...) {
  # init = 0 starts at the prior mode (identity correlations, unit scales),
  # avoiding degenerate cholesky_factor_corr draws that can abort a chain at
  # high trait counts. Overridden by an explicit `init` in `...`.
  control <- .rstan_control(.check_control(control))
  args <- list(object = model, data = stan_data, iter = iter, warmup = warmup,
               chains = chains, control = control, init = 0)
  do.call(rstan::sampling, utils::modifyList(args, list(...)))
}
