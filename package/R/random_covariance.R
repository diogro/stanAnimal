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

#' Random covariance matrices
#'
#' Draws random covariance matrices for simulations: a random correlation
#' matrix scaled by a vector of variances. Correlations come either from the
#' LKJ distribution used as a prior by \code{\link{lmm_animal}}, or from the
#' partial-correlation construction of \code{evolqg::RandomMatrix()}.
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
#'   \code{evolqg::RandomMatrix()}.
#' @param lkj_prior LKJ shape for \code{method = "lkj"}, a finite number above
#'   0. The default \code{2} matches the prior of \code{\link{lmm_animal}}.
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
#'   \code{evolqg::RandomMatrix()}, but it samples correlations near \eqn{\pm
#'   1} and frequently produces nearly singular matrices: at five traits the
#'   median condition number is in the hundreds. Use it to reproduce that
#'   construction, not as a neutral simulation default.
#'
#'   Unlike \code{evolqg::RandomMatrix()}, this function also handles one and
#'   two traits with \code{method = "partial"}.
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
#' @export
random_covariance <- function(n_traits, n_matrices = 1, variances = NULL,
                              min_var = 1, max_var = 1,
                              method = c("lkj", "partial"), lkj_prior = 2,
                              trait_names = NULL) {
  n_traits <- .check_count(n_traits, "n_traits")
  n_matrices <- .check_count(n_matrices, "n_matrices")
  method <- match.arg(method)
  min_var <- .check_positive(min_var, "min_var")
  max_var <- .check_positive(max_var, "max_var")
  if (min_var > max_var)
    stop("`min_var` must not exceed `max_var`.", call. = FALSE)
  lkj_prior <- .check_positive(lkj_prior, "lkj_prior")

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

  draw_one <- function() {
    v <- variances %||% stats::runif(n_traits, min_var, max_var)
    correlation <- if (method == "lkj") .rlkj(n_traits, lkj_prior)
                   else .random_correlation_partial(n_traits)
    out <- correlation * sqrt(outer(v, v))
    dimnames(out) <- list(trait_names, trait_names)
    out
  }

  if (n_matrices == 1L) return(draw_one())
  lapply(seq_len(n_matrices), function(i) draw_one())
}
