#' Fit an animal model with a Gibbs sampler
#'
#' Fits a univariate or multivariate linear animal model with fixed effects and relatedness
#' matrix using a single-chain Gibbs sampler. This is intended for quick
#' simulations and checks; its priors differ from those in \code{lmm_animal}.
#' @param Y Numeric response vector or matrix, with individuals in rows. Every
#'   trait must vary, and is checked as in \code{\link{lmm_animal}}.
#' @param X Optional numeric full-column-rank fixed-effects design matrix, with
#'   the same rows as \code{Y}. The default \code{NULL} fits an intercept-only
#'   model: one column of ones, named \code{"(Intercept)"}. Unlike the Stan
#'   models, it needs no column representing a constant: the sampler does not
#'   shift the response mean through the design.
#' @param A Symmetric positive-semidefinite relationship matrix, with rows and
#'   columns in the same individual order as the rows of \code{Y}. When both
#'   \code{A} and \code{Y} have row names they are used to put \code{A} in the
#'   order of \code{Y}; names that disagree are an error.
#' @param iter Total number of Gibbs iterations, a whole number.
#' @param warmup Number of initial iterations to discard, a whole number smaller
#'   than \code{iter}.
#' @param thin Keep every \code{thin}th iteration after burn-in, a whole number;
#'   \code{warmup + thin} must not exceed \code{iter}.
#' @return A list of posterior draws: \code{G} and \code{R} are covariance
#'   arrays, \code{beta} contains fixed-effect matrix draws, \code{h2}
#'   contains trait-wise variance ratios, and \code{iterations} gives the
#'   retained iteration numbers. Quantities are indexed by trait first and draw
#'   last, from named dimensions, the layout that
#'   \code{\link{back_transform}} returns as well.
#' @details Fixed effects have an improper flat prior. Genetic and residual
#'   covariance matrices have inverse-Wishart priors with degrees of freedom
#'   \code{K + 1} and scale matrix \code{(K + 1) I / 2}, where \code{K} is
#'   the number of traits. These differ from the priors used by
#'   \code{lmm_animal}; the resulting posterior draws are not directly
#'   interchangeable.
#' @export
lmm_animal_gibbs <- function(Y, X = NULL, A, iter = 1000, warmup = floor(iter / 2), thin = 1) {
  .check_iterations(iter, warmup, thin)
  # No mean is shifted through the design here, so X needs no constant column.
  YX <- .check_response_design(Y, X, require_constant = FALSE)
  Y <- YX$Y
  X <- YX$X

  n <- nrow(Y)
  k <- ncol(Y)
  p <- ncol(X)
  A <- .check_relationship_matrix(as.matrix(A), Y)
  A <- .check_symmetric(A)
  .check_positive_semidefinite(A)

  eig_A <- eigen(A, symmetric = TRUE)
  lambda <- pmax(eig_A$values, max(1e-12, max(abs(eig_A$values)) * 1e-8))

  y_centered <- sweep(Y, 2L, colMeans(Y), "-")
  y_cov <- crossprod(y_centered) / (n - 1)

  keep <- seq.int(warmup + thin, iter, by = thin)
  n_keep <- length(keep)
  Q <- eig_A$vectors
  Ys <- crossprod(Q, Y)
  Xs <- crossprod(Q, X)
  beta_cov <- solve(crossprod(Xs))
  beta <- beta_cov %*% crossprod(Xs, Ys)

  V <- diag(k) / 2
  nu <- k + 1
  prior_scale <- nu * V
  initial_var <- pmax(diag(y_cov) / 2, 1e-6)
  G <- R <- diag(initial_var, nrow = k)
  u <- matrix(0, nrow = n, ncol = k)

  out_G <- array(NA_real_, c(k, k, n_keep),
                 dimnames = list(colnames(Y), colnames(Y), NULL))
  out_R <- array(NA_real_, c(k, k, n_keep),
                 dimnames = list(colnames(Y), colnames(Y), NULL))
  out_beta <- array(NA_real_, c(k, p, n_keep),
                    dimnames = list(colnames(Y), colnames(X), NULL))
  out_h2 <- matrix(NA_real_, nrow = k, ncol = n_keep,
                   dimnames = list(colnames(Y), NULL))
  next_keep <- 1L

  for (iteration in seq_len(iter)) {
    Lr <- t(chol(R))
    residual <- Ys - Xs %*% beta
    W <- t(forwardsolve(Lr, t(residual)))
    Gw <- forwardsolve(Lr, t(forwardsolve(Lr, G)))
    eig_G <- eigen((Gw + t(Gw)) / 2, symmetric = TRUE)
    s <- outer(lambda, pmax(eig_G$values, 0))
    shrink <- s / (s + 1)
    Uw <- (W %*% eig_G$vectors) * shrink +
      matrix(stats::rnorm(n * k), nrow = n, ncol = k) * sqrt(shrink)
    u <- Uw %*% t(eig_G$vectors) %*% t(Lr)

    beta_mean <- beta_cov %*% crossprod(Xs, Ys - u)
    beta <- beta_mean + t(chol(beta_cov)) %*%
      matrix(stats::rnorm(p * k), nrow = p, ncol = k) %*% chol(R)

    G <- .rinvwishart(nu + n, prior_scale +
                        crossprod(sweep(u, 1L, sqrt(lambda), "/")))
    residual <- Ys - Xs %*% beta - u
    R <- .rinvwishart(nu + n, prior_scale + crossprod(residual))

    if (iteration == keep[next_keep]) {
      out_G[, , next_keep] <- G
      out_R[, , next_keep] <- R
      out_beta[, , next_keep] <- t(beta)
      out_h2[, next_keep] <- diag(G) / (diag(G) + diag(R))
      next_keep <- next_keep + 1L
      if (next_keep > n_keep) break
    }
  }

  list(G = out_G, R = out_R, beta = out_beta, h2 = out_h2,
       iterations = keep)
}

.rinvwishart <- function(df, scale) {
  k <- nrow(scale)
  wishart <- matrix(stats::rWishart(1L, df, solve(scale))[, , 1L],
                    nrow = k, ncol = k)
  draw <- solve(wishart)
  (draw + t(draw)) / 2
}