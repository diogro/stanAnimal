#' Share of variance in the leading principal component
#'
#' Diagnostic for the size and shape reparametrization: when a single component
#' dominates, the traits are mostly one size gradient plus noise, and fitting
#' size and shape scores is a much better behaved model than fitting all of the
#' raw traits.
#'
#' @param Y Numeric matrix of traits, individuals in rows.
#' @return The share of the total variance taken by the first principal
#'   component of the correlation matrix of \code{Y}, between \code{1 / ncol(Y)}
#'   and 1 (1 for a single trait).
#' @details The correlation matrix is used so that traits measured on different
#'   scales do not dominate through their units alone. A share near one means
#'   the traits are one gradient, which is what \code{\link{size_shape}} is for;
#'   a share close to \code{1 / ncol(Y)} means the traits carry independent
#'   information and no reparametrization is needed.
#' @seealso \code{\link{size_shape}}, \code{\link{animal_design}}
#' @examples
#' gradient <- rnorm(50)
#' correlated <- cbind(a = gradient + rnorm(50, sd = 0.2),
#'                     b = 2 * gradient + rnorm(50, sd = 0.2),
#'                     c = -gradient + rnorm(50, sd = 0.2))
#' pc_dominance(correlated)                    # close to 1
#' pc_dominance(matrix(rnorm(150), 50, 3))     # close to 1 / 3
#' @export
pc_dominance <- function(Y) {
  Y <- as.matrix(Y)
  if (!is.numeric(Y) || any(!is.finite(Y)))
    stop("Y must be a finite numeric matrix of traits.", call. = FALSE)
  k <- ncol(Y)
  if (is.null(k) || k < 1L)
    stop("Y must have at least one column of traits.", call. = FALSE)
  if (k == 1L)
    return(1)
  if (nrow(Y) < 3L)
    stop("pc_dominance needs at least three individuals.", call. = FALSE)
  if (any(apply(Y, 2L, stats::sd) == 0))
    stop("every trait in Y must vary.", call. = FALSE)
  eigenvalues <- eigen(stats::cor(Y), symmetric = TRUE, only.values = TRUE)$values
  eigenvalues[1L] / sum(eigenvalues)
}

#' Size and shape scores for many correlated traits
#'
#' Rewrites a set of correlated traits as log size plus orthogonal shape axes,
#' the reparametrization used to fit many allometric traits: raw traits are
#' usually dominated by a single size gradient, so their genetic covariance is a
#' long, thin ridge that is slow to sample. Size and shape span the same space,
#' so nothing is lost, but the sampler works on well-separated axes.
#'
#' @param Y Numeric matrix of traits, individuals in rows and traits in columns.
#'   Log-transformed measurements are the usual input.
#' @param n_shape Number of shape axes to keep, between 1 and
#'   \code{ncol(Y) - 1}. The default keeps all of them, so the transformation is
#'   a bijection; fewer axes give a lower-rank (singular) covariance in the
#'   original trait space.
#' @return A list of class \code{size_shape} with
#'   \item{scores}{matrix with the size and shape scores, one column each, named
#'     \code{size}, \code{shape1}, \code{shape2}, ...;}
#'   \item{loadings}{\code{ncol(Y)} by \code{1 + n_shape} matrix \code{B} such
#'     that, up to a constant absorbed by the model intercepts,
#'     \code{Y = scores \%*\% t(B)};}
#'   \item{traits}{names of the model traits (\code{size}, \code{shape1}, ...);}
#'   \item{raw_traits}{names of the original traits;}
#'   \item{centre}{column means removed from the shape axes.}
#' @details Size is the mean of the traits, \code{rowMeans(Y)}, i.e. the log
#'   geometric mean for log measurements. Shape is \code{Y - size}, which sums
#'   to zero within each individual and therefore lives in an
#'   \code{ncol(Y) - 1} dimensional subspace; its principal components are the
#'   shape axes. The leading axis of variation is size, and the shape axes are
#'   ordered by the variance they explain.
#'
#'   The transformation is linear, so genetic and residual covariance matrices
#'   transform as \eqn{B M B^T}; use \code{\link{back_transform}} to do this for
#'   the draws of a fit. Signs of the shape axes are arbitrary.
#' @seealso \code{\link{back_transform}}, \code{\link{animal_design}}
#' @examples
#' size_gradient <- rnorm(20)
#' Y <- cbind(a = size_gradient + rnorm(20),
#'            b = 2 * size_gradient + rnorm(20),
#'            c = -size_gradient + rnorm(20))
#' ss <- size_shape(Y)
#' head(ss$scores)
#' ss$loadings
#' @export
size_shape <- function(Y, n_shape = NULL) {
  Y <- as.matrix(Y)
  if (!is.numeric(Y) || any(!is.finite(Y)))
    stop("Y must be a finite numeric matrix of traits.", call. = FALSE)
  k <- ncol(Y)
  if (is.null(k) || k < 2L)
    stop("size_shape needs at least two traits.", call. = FALSE)
  if (nrow(Y) < 3L)
    stop("size_shape needs at least three individuals.", call. = FALSE)
  if (is.null(n_shape))
    n_shape <- k - 1L
  if (!is.numeric(n_shape) || length(n_shape) != 1L || is.na(n_shape) ||
      n_shape != floor(n_shape) || n_shape < 1L || n_shape > k - 1L)
    stop("n_shape must be an integer between 1 and ncol(Y) - 1 (", k - 1L, ").",
         call. = FALSE)
  n_shape <- as.integer(n_shape)

  size <- rowMeans(Y)
  shape <- Y - size
  pca <- stats::prcomp(shape)
  axes <- seq_len(n_shape)
  B_shape <- pca$rotation[, axes, drop = FALSE]
  B <- cbind(size = 1, B_shape)

  raw_traits <- colnames(Y) %||% paste0("trait", seq_len(k))
  scores <- cbind(size, pca$x[, axes, drop = FALSE])
  colnames(scores) <- c("size", paste0("shape", axes))
  rownames(B) <- raw_traits
  rownames(scores) <- rownames(Y)

  structure(list(scores = scores, loadings = B, traits = colnames(scores),
                 raw_traits = raw_traits, centre = colMeans(shape)),
            class = "size_shape")
}

#' @export
print.size_shape <- function(x, ...) {
  cat("<size_shape>", nrow(x$scores), "individuals,",
      length(x$traits), "model traits from", length(x$raw_traits),
      "raw traits\n")
  cat("model traits:", paste(x$traits, collapse = ", "), "\n")
  cat("raw traits:  ", paste(x$raw_traits, collapse = ", "), "\n")
  invisible(x)
}

#' Back-transform posterior draws to the original trait scale
#'
#' Takes a fit of the size and shape scores made with
#' \code{\link{lmm_animal}} and returns draws of the genetic and residual
#' covariance matrices, the phenotypic covariance, the genetic correlations and
#' the heritabilities on the scale of the original traits.
#'
#' @param fit A fit returned by \code{\link{lmm_animal}} (a
#'   \code{CmdStanMCMC} object), or a \code{posterior} draws object, or a matrix
#'   of draws, containing the variables \code{G} and \code{E}, or \code{P} for a
#'   phenotypic fit from \code{\link{lmm_multiCov}}.
#' @param transform A \code{size_shape} object from \code{\link{size_shape}}.
#' @return A \code{draws_matrix} on the scale of the original traits. For an
#'   animal fit it holds \code{G[i,j]}, \code{E[i,j]}, \code{P[i,j]},
#'   \code{corrG[i,j]} and \code{corrP[i,j]} (below the diagonal) and
#'   \code{h2[i]}; for a phenotypic fit it holds \code{P[i,j]},
#'   \code{corrP[i,j]} (below the diagonal) and \code{sigma_P[i]}.
#' @details For each draw the covariance matrices are mapped as
#'   \eqn{M \to B M B^T}, which is exact because size and shape are a linear
#'   reparametrization of the traits. Heritabilities and correlations are then
#'   recomputed from the transformed matrices: unlike the covariance, they are
#'   not invariant under a change of traits, so the h2 of a shape axis is not
#'   the h2 of any original trait.
#'
#'   The fixed effects \code{beta} stay on the model scale, since a coefficient
#'   for size or a shape axis does not correspond to a single original trait.
#' @seealso \code{\link{size_shape}}, \code{\link{animal_draws}}
#' @examples
#' ss <- size_shape(cbind(a = rnorm(20), b = rnorm(20), c = rnorm(20)))
#' k <- length(ss$raw_traits)
#' names_mat <- function(var) sprintf(paste0(var, "[%d,%d]"),
#'                                    rep(seq_len(k), k), rep(seq_len(k), each = k))
#' # a real fit supplies G and E; here they are fixed at the identity
#' m <- cbind(matrix(rep(as.vector(diag(k)), 4), 4, byrow = TRUE),
#'            matrix(rep(as.vector(diag(k)), 4), 4, byrow = TRUE))
#' colnames(m) <- c(names_mat("G"), names_mat("E"))
#' post <- back_transform(m, ss)
#' posterior::summarise_draws(
#'   posterior::subset_draws(post, variable = c("h2[1]", "G[1,1]", "corrG[1,2]")))
#' @export
back_transform <- function(fit, transform) {
  if (!inherits(transform, "size_shape"))
    stop("transform must be a size_shape object from size_shape().",
         call. = FALSE)
  B <- transform$loadings
  k_model <- ncol(B)
  k <- nrow(B)

  draws <- if (inherits(fit, "CmdStanFit") || (is.list(fit) && is.function(fit$draws))) {
    variables <- if (is.function(fit$metadata)) fit$metadata()$stan_variables else NULL
    wanted <- intersect(c("G", "E", "P"), variables %||% c("G", "E", "P"))
    if (length(wanted)) fit$draws(wanted) else fit$draws()
  } else fit
  m <- if (inherits(draws, "draws")) unclass(posterior::as_draws_matrix(draws))
       else as.matrix(draws)
  if (!is.matrix(m) || is.null(colnames(m)))
    stop("draws must be a fit, a posterior draws object or a named matrix.",
         call. = FALSE)

  in_g <- .matrix_variable_names("G", k_model)
  in_e <- .matrix_variable_names("E", k_model)
  in_p <- .matrix_variable_names("P", k_model)
  has_components <- all(in_g %in% colnames(m)) && all(in_e %in% colnames(m))
  if (!has_components && !all(in_p %in% colnames(m)))
    stop("draws must contain G and E (lmm_animal) or P (lmm_multiCov) for ",
         k_model, " model traits, as named by the model functions for the ",
         "size_shape scores.", call. = FALSE)

  out_g <- .matrix_variable_names("G", k)
  out_e <- .matrix_variable_names("E", k)
  out_p <- .matrix_variable_names("P", k)
  lower <- upper.tri(diag(k))
  triangle <- which(lower, arr.ind = TRUE)

  if (has_components) {
    out <- matrix(NA_real_, nrow(m), k * k * 3L + 2L * nrow(triangle) + k,
                  dimnames = list(NULL,
                                  c(out_g, out_e, out_p,
                                    sprintf("corrG[%d,%d]", triangle[, 1], triangle[, 2]),
                                    sprintf("corrP[%d,%d]", triangle[, 1], triangle[, 2]),
                                    sprintf("h2[%d]", seq_len(k)))))
    for (i in seq_len(nrow(m))) {
      G <- B %*% matrix(m[i, in_g], k_model) %*% t(B)
      E <- B %*% matrix(m[i, in_e], k_model) %*% t(B)
      P <- G + E
      out[i, ] <- c(G, E, P, cov2cor(G)[lower], cov2cor(P)[lower],
                    diag(G) / diag(P))
    }
  } else {
    out <- matrix(NA_real_, nrow(m), k * k + nrow(triangle) + k,
                  dimnames = list(NULL,
                                  c(out_p,
                                    sprintf("corrP[%d,%d]", triangle[, 1], triangle[, 2]),
                                    sprintf("sigma_P[%d]", seq_len(k)))))
    for (i in seq_len(nrow(m))) {
      P <- B %*% matrix(m[i, in_p], k_model) %*% t(B)
      out[i, ] <- c(P, cov2cor(P)[lower], sqrt(diag(P)))
    }
  }
  posterior::as_draws_matrix(out)
}

# Names of a matrix-valued Stan variable in column-major order, as cmdstanr
# reports them.
.matrix_variable_names <- function(name, k) {
  sprintf("%s[%d,%d]", name, rep(seq_len(k), k), rep(seq_len(k), each = k))
}
