# Random covariance matrices for simulations.
#
# The construction is the one in evolqg::RandomMatrix(): draw a correlation
# matrix, scale it by a vector of variances. The partial-correlation method is
# a port of that function, with the one- and two-trait cases fixed (the
# original indexes 2:1 and 3:2 there), and the LKJ method reuses the same
# generator as the prior draws.

# Random correlation matrix from uniformly drawn partial correlations.
# `tolerance` is the threshold below which a partial correlation is treated as
# zero, as `ke` in evolqg::RandomMatrix().
.random_correlation_partial <- function(k, tolerance = 1e-3) {
  if (k == 1L) return(matrix(1, 1, 1))

  correlation <- matrix(0, k, k)
  b <- matrix(0, k, k)
  b[lower.tri(b, diag = TRUE)] <- 1

  # First column: correlations with trait 1, uniform on (-1, 1).
  correlation[2:k, 1] <- stats::runif(k - 1L, -1, 1)
  b[2:k, 1] <- correlation[2:k, 1]
  for (i in 2:k) b[i, 2:i] <- sqrt(1 - correlation[i, 1]^2)

  # Later columns: draw a partial correlation, then turn it into the
  # corresponding correlation and update the partial-correlation rows below it.
  if (k >= 3L) {
    for (i in seq.int(3L, k)) {
      for (j in seq.int(2L, i - 1L)) {
        b1 <- as.vector(b[j, seq_len(j - 1L)] %*% b[i, seq_len(j - 1L)])
        b2 <- b[j, j] * b[i, j]
        lower <- b1 - b2
        upper <- b1 + b2
        correlation[i, j] <- if (b2 < tolerance) b1 else
          lower + (upper - lower) * stats::runif(1)
        cos_angle <- (correlation[i, j] - b1) / b2
        if (is.finite(cos_angle)) {
          if (cos_angle > 1) {
            b[i, seq.int(j + 1L, k)] <- 0
          } else if (cos_angle < -1) {
            b[i, j] <- -b[i, j]
            b[i, seq.int(j + 1L, k)] <- 0
          } else {
            b[i, j] <- b[i, j] * cos_angle
            b[i, seq.int(j + 1L, k)] <-
              b[i, seq.int(j + 1L, k)] * sqrt(1 - cos_angle^2)
          }
        }
      }
    }
  }

  correlation <- correlation + t(correlation) + diag(1, k)
  permutation <- sample(k)
  correlation[permutation, permutation]
}

# Relative eigenvalue variance, the integration index returned by
# evolqg::CalcEigenVar(rel = TRUE): the observed eigenvalue variance divided by
# the maximum for a matrix of the same size and total variance. Zero
# eigenvalues are kept: a rank-one matrix has a relative variance of 1, and
# dropping the zeros would call it 0.
.relative_eigen_variance <- function(S) {
  eigenvalues <- pmax(eigen(S, symmetric = TRUE, only.values = TRUE)$values, 0)
  k <- length(eigenvalues)
  if (k < 2L) return(0)
  average <- mean(eigenvalues)
  observed <- sum((eigenvalues - average)^2) / k
  maximum <- (k - 1) * sum(eigenvalues)^2 / k^2
  observed / maximum
}

# A symmetric matrix with the same correlations as `R` but with a positive
# leading eigenvector: flip the sign of the traits that load negatively on PC1.
.orient_positive <- function(R) {
  axis <- eigen(R, symmetric = TRUE)$vectors[, 1L]
  signs <- ifelse(axis < 0, -1, 1)
  R * tcrossprod(signs)
}

# A random correlation matrix with a positive size axis and a chosen
# integration level. The random part is a sign-oriented LKJ correlation shrunk
# toward independence, so it keeps a positive PC1 while the other components
# are free to have negative correlations. Two regimes share that random part:
#
#   weak integration:    R = (1 - alpha) * I + alpha * random
#   strong integration:  R = rho * ones + (1 - rho) * random
#
# The first shrinks the whole matrix toward independence, which reaches targets
# below the integration of the random part without flattening it into
# equicorrelation. The second adds the equal-loading size factor; at rho = 1
# the matrix is all ones, whose relative eigenvalue variance is 1, so every
# target above the random part is reachable. Both regimes have diagonal 1 and
# a positive PC1, and the two meet where alpha = 1 and rho = 0.
.factor_correlation <- function(k, integration, lkj_prior = 2,
                                residual = 0.3) {
  if (k == 1L || integration <= 0) return(diag(1, k))

  identity <- diag(1, k)
  random <- (1 - residual) * identity +
    residual * .orient_positive(.rlkj(k, lkj_prior))

  baseline <- .relative_eigen_variance(random)
  if (integration <= baseline) {
    # Shrinking toward identity scales every eigenvalue deviation from 1 by
    # alpha, so the relative variance scales by alpha^2 and the root is exact.
    alpha <- sqrt(integration / baseline)
    return((1 - alpha) * identity + alpha * random)
  }

  ones <- matrix(1, k, k)
  strength <- stats::uniroot(
    function(rho)
      .relative_eigen_variance(rho * ones + (1 - rho) * random) - integration,
    c(0, 1 - 1e-9), tol = 1e-10
  )$root
  strength * ones + (1 - strength) * random
}

# `integration` is NULL (the default range), one target, or a range to draw
# from.
.check_integration <- function(integration) {
  if (is.null(integration)) return(c(0.1, 0.8))
  if (!is.numeric(integration) || !length(integration) %in% c(1L, 2L) ||
      any(!is.finite(integration)) || any(integration < 0) ||
      any(integration >= 1))
    stop("`integration` must be one or two numbers in [0, 1).", call. = FALSE)
  if (length(integration) == 2L && integration[1L] > integration[2L])
    stop("`integration` must be an increasing range.", call. = FALSE)
  as.numeric(integration)
}

#' Random covariance matrices
#'
#' Draws random covariance matrices for simulations: a random correlation
#' matrix scaled by a vector of variances. Correlations come from the LKJ
#' distribution used as a prior by \code{\link{lmm_animal}}, from the
#' partial-correlation construction of \code{evolqg::RandomMatrix()}, or from a
#' two-regime model that combines a positive size axis with a sign-oriented
#' random LKJ residual.
#'
#' @param n_traits Number of traits, a whole number of at least 1.
#' @param n_matrices Number of matrices to draw, a whole number of at least 1.
#' @param variances Optional vector of variances, one per trait, with finite
#'   positive values. The default draws them uniformly between \code{min_var}
#'   and \code{max_var}.
#' @param min_var,max_var Bounds of the uniform distribution of variances when
#'   \code{variances} is \code{NULL}: finite positive numbers with
#'   \code{min_var} not above \code{max_var}.
#' @param method Correlation model. \code{"lkj"} (default) draws a correlation
#'   matrix from the LKJ distribution with shape \code{lkj_prior};
#'   \code{"partial"} draws partial correlations uniformly, reproducing
#'   \code{evolqg::RandomMatrix()}; \code{"factor"} combines a positive size
#'   factor with a sign-oriented random LKJ residual, so PC1 is a size axis and
#'   draws differ at a fixed integration level.
#' @param lkj_prior LKJ shape for \code{method = "lkj"} and for the random
#'   residual of \code{method = "factor"}, a finite number above 0. The default
#'   \code{2} matches the prior of \code{\link{lmm_animal}}.
#' @param integration For \code{method = "factor"}, the target integration
#'   level, measured as the relative eigenvalue variance of
#'   \code{evolqg::CalcEigenVar(rel = TRUE)} and so between 0 and 1. A single
#'   number targets that level for every matrix; two numbers draw the target
#'   uniformly between them. The default \code{NULL} draws from 0.1 to 0.8,
#'   spanning weak to strong integration. Supplying it with another method is
#'   an error.
#' @param trait_names Optional trait names, one per trait; the default is
#'   \code{trait1}, \code{trait2}, ...
#' @return With \code{n_matrices = 1}, a named symmetric positive-definite
#'   matrix. With more than one matrix, a list of such matrices.
#' @details The variances are independent of the correlations, and the
#'   correlation matrix \eqn{R} is scaled as \eqn{R \odot \sqrt{v v^T}}, so the
#'   diagonal is \eqn{v}. With the default \code{min_var = max_var = 1} the
#'   result is a correlation matrix, as in \code{evolqg::RandomMatrix()}; using
#'   those defaults for both \code{G} and \code{E} gives a heritability of 0.5
#'   for every trait. Use \code{variances}, or a wider
#'   \code{min_var}/\code{max_var}, to vary the heritabilities.
#'
#'   The LKJ method is the default because it matches the prior family of the
#'   models and stays well conditioned: correlations are shrunk toward zero by
#'   \code{lkj_prior}, and a five-trait matrix typically has a condition number
#'   in the tens. The partial-correlation method reproduces
#'   \code{evolqg::RandomMatrix()} and is the one to use when the simulation
#'   should span a realistic range of integration: at five traits its first
#'   principal component can explain up to about 90\% of the variance and the
#'   mean squared correlation can approach 0.8, while LKJ(2) stays well below
#'   that. The cost is conditioning: the strongly integrated draws often have
#'   condition numbers in the hundreds or more, whereas real size-dominated
#'   matrices have a single positive axis and can reach the same integration
#'   with much smaller condition numbers. Use the partial method to span the
#'   range, or deliberately as a hard case for an estimator, not as a neutral
#'   default.
#'
#'   Unlike \code{evolqg::RandomMatrix()}, this function also handles one and
#'   two traits with \code{method = "partial"}.
#'
#'   With \code{method = "factor"} the random part is
#'   \eqn{R_{res} = (1 - \delta) I + \delta R_{LKJ}}, where \eqn{R_{LKJ}} is a
#'   random LKJ correlation whose leading eigenvector has been oriented to be
#'   positive and \eqn{\delta = 0.3}. Orienting keeps PC1 positive while leaving
#'   the remaining components free to carry negative correlations. For targets
#'   below the integration of \eqn{R_{res}}, the matrix is
#'   \eqn{(1 - \alpha) I + \alpha R_{res}}, which shrinks the random structure
#'   toward independence without flattening it into equicorrelation. For higher
#'   targets the matrix is \eqn{\rho 11^T + (1 - \rho) R_{res}}, which adds the
#'   equal-loading size factor; at \eqn{\rho = 1} it is all ones, whose relative
#'   eigenvalue variance is 1. Both regimes have diagonal 1 and are solved to
#'   match the requested integration. Because \eqn{R_{LKJ}} is redrawn for every
#'   matrix, draws at a fixed integration differ in their lower components. This
#'   is the method to use for size-dominated morphometric simulations: it spans
#'   the full integration range while staying far better conditioned than
#'   bending an LKJ spectrum directly (at five traits the median condition
#'   number is about 4 at 0.1, 16 at 0.5 and 53 at 0.8). At two traits the
#'   target fixes the single correlation, so draws at a fixed integration are
#'   identical up to the sign orientation.
#' @seealso \code{\link{animal_prior_draws}}, \code{\link{size_shape}}
#' @examples
#' set.seed(1)
#' G <- random_covariance(3, min_var = 0.5, max_var = 1.5)
#' round(G, 2)
#' round(cov2cor(G), 2)
#' eigen(G, symmetric = TRUE, only.values = TRUE)$values
#'
#' # Reproduce the evolqg construction
#' random_covariance(2, n_matrices = 3, method = "partial")
#'
#' # A strongly integrated size axis: PC1 dominates and is positive
#' S <- random_covariance(5, method = "factor", integration = 0.8)
#' R <- cov2cor(S)
#' round(R, 2)
#' eigen(R, symmetric = TRUE)$vectors[, 1]
#' @export
random_covariance <- function(n_traits, n_matrices = 1, variances = NULL,
                              min_var = 1, max_var = 1,
                              method = c("lkj", "partial", "factor"),
                              lkj_prior = 2, integration = NULL,
                              trait_names = NULL) {
  n_traits <- .check_count(n_traits, "n_traits")
  n_matrices <- .check_count(n_matrices, "n_matrices")
  method <- match.arg(method)
  min_var <- .check_positive(min_var, "min_var")
  max_var <- .check_positive(max_var, "max_var")
  if (min_var > max_var)
    stop("`min_var` must not exceed `max_var`.", call. = FALSE)
  lkj_prior <- .check_positive(lkj_prior, "lkj_prior")

  if (method != "factor" && !is.null(integration))
    stop("`integration` is only used by method = \"factor\".", call. = FALSE)
  if (method == "factor") {
    integration <- .check_integration(integration)
    targets <- if (length(integration) == 1L) {
      rep(integration, n_matrices)
    } else {
      stats::runif(n_matrices, integration[1L], integration[2L])
    }
  } else {
    targets <- rep(NA_real_, n_matrices)
  }

  if (!is.null(variances)) {
    if (!is.numeric(variances) || length(variances) != n_traits ||
        any(!is.finite(variances)) || any(variances <= 0))
      stop("`variances` must be one finite positive number per trait (",
           n_traits, ").", call. = FALSE)
    variances <- as.numeric(variances)
  }

  if (is.null(trait_names)) {
    trait_names <- paste0("trait", seq_len(n_traits))
  } else {
    if (length(trait_names) != n_traits)
      stop("`trait_names` must have one name per trait (", n_traits, ").",
           call. = FALSE)
    trait_names <- as.character(trait_names)
  }

  draw_one <- function(target) {
    v <- variances %||% stats::runif(n_traits, min_var, max_var)
    correlation <- switch(method,
                          lkj = .rlkj(n_traits, lkj_prior),
                          partial = .random_correlation_partial(n_traits),
                          factor = .factor_correlation(n_traits, target,
                                                       lkj_prior))
    out <- correlation * sqrt(outer(v, v))
    dimnames(out) <- list(trait_names, trait_names)
    out
  }

  if (n_matrices == 1L) return(draw_one(targets[1L]))
  lapply(seq_len(n_matrices), function(i) draw_one(targets[i]))
}
