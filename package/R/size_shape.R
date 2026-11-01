#' Share of variance explained by the first principal component
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
#' pc1_variance_share(correlated)                 # close to 1
#' pc1_variance_share(matrix(rnorm(150), 50, 3))  # close to 1 / 3
#' @export
pc1_variance_share <- function(Y) {
  Y <- as.matrix(Y)
  if (!is.numeric(Y) || any(!is.finite(Y)))
    stop("Y must be a finite numeric matrix of traits.", call. = FALSE)
  k <- ncol(Y)
  if (is.null(k) || k < 1L)
    stop("Y must have at least one column of traits.", call. = FALSE)
  if (k == 1L)
    return(1)
  if (nrow(Y) < 3L)
    stop("pc1_variance_share needs at least three individuals.", call. = FALSE)
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
#'   shape axes. Size is often close to the leading axis of variation of a
#'   morphometric data set, but it need not be: it is the unweighted mean of the
#'   traits, while the principal components are shaped by their covariance
#'   structure. The shape axes are ordered by the residual variance they
#'   explain.
#'
#'   The transformation is linear, so genetic and residual covariance matrices
#'   transform as \eqn{B M B^T}; use \code{\link{back_transform}} to do this for
#'   the draws of a fit. Signs of the shape axes are arbitrary.
#'
#'   The prior is not invariant under this change of traits, and that is usually
#'   an asset. The LKJ prior that \code{\link{lmm_animal}} puts on the
#'   correlation matrix acts on correlations between size and shape instead of
#'   correlations between raw traits. The shape axes are principal components of
#'   \code{Y - size}, so they are uncorrelated in the sample and their genetic
#'   correlations are expected to be small; a larger \code{lkj_prior} than for
#'   the same traits fitted raw (stronger shrinkage towards zero correlation) is
#'   therefore reasonable, and it can also make the model cheaper to sample.
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
#'   phenotypic fit from \code{\link{lm_phen_cov}}.
#' @param transform A \code{size_shape} object from \code{\link{size_shape}}.
#' @return A list of class \code{size_shape_draws} holding one element per
#'   quantity, in the same layout as \code{\link{lmm_animal_gibbs}}: traits
#'   first and draws last, with every dimension named after the traits. For an
#'   animal fit it holds \code{G}, \code{E}, \code{P}, \code{corrG} and
#'   \code{corrP} as \code{traits} by \code{traits} by \code{draws} arrays,
#'   and \code{h2} as a \code{traits} by \code{draws} matrix. For a
#'   phenotypic fit it holds \code{P} and \code{corrP} as arrays and
#'   \code{sigma_P} as a matrix. \code{traits} holds the trait names.
#' @details For each draw the covariance matrices are mapped as
#'   \eqn{M \to B M B^T}, which is exact because size and shape are a linear
#'   reparametrization of the traits. Heritabilities and correlations are then
#'   recomputed from the transformed matrices: unlike the covariance, they are
#'   not invariant under a change of traits, so the h2 of a shape axis is not
#'   the h2 of any original trait.
#'
#'   The fixed effects \code{beta} stay on the model scale, since a coefficient
#'   for size or a shape axis does not correspond to a single original trait.
#'
#'   Draws are the last dimension, as in \code{\link{lmm_animal_gibbs}}, so
#'   posterior summaries are one \code{apply()} away:
#'   \code{apply(post$G, c(1, 2), mean)} is the posterior mean of \code{G} on
#'   the trait scale, and \code{apply(post$corrG, c(1, 2), quantile, 0.05)} its
#'   5\% quantile.
#' @seealso \code{\link{size_shape}}, \code{\link{rename_beta}},
#'   \code{\link{to_draws_list}}
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
#' dim(post$G)                    # traits x traits x draws
#' apply(post$G, c(1, 2), mean)   # posterior mean G on the trait scale
#' apply(post$h2, 1, mean)
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
    stop("draws must contain G and E (lmm_animal) or P (lm_phen_cov) for ",
         k_model, " model traits, as named by the model functions for the ",
         "size_shape scores.", call. = FALSE)

  traits <- transform$raw_traits
  if (has_components) {
    components <- list(G = .transform_draws(m, "G", B, traits),
                       E = .transform_draws(m, "E", B, traits))
    components$P <- components$G + components$E
    components$corrG <- .corr_draws(components$G)
    components$corrP <- .corr_draws(components$P)
    components$h2 <- .diag_draws(components$G) / .diag_draws(components$P)
  } else {
    components <- list(P = .transform_draws(m, "P", B, traits))
    components$corrP <- .corr_draws(components$P)
    components$sigma_P <- sqrt(.diag_draws(components$P))
  }
  structure(c(components, list(traits = traits)), class = "size_shape_draws")
}

#' Convert the back-transformed draws to posterior draws
#'
#' Reshapes the arrays returned by \code{\link{back_transform}} into
#' \code{\link[posterior]{draws_matrix}} objects, one per quantity, so that
#' the \pkg{posterior} and \pkg{bayesplot} tools can be used on them.
#'
#' @param x A \code{size_shape_draws} object from
#'   \code{\link{back_transform}}.
#' @return A named list of \code{draws_matrix} objects, one per quantity of
#'   \code{x}, with the variables named as the models report them:
#'   \code{G[i,j]}, \code{E[i,j]} and \code{P[i,j]} for every element of the
#'   covariance matrices, \code{corrG[i,j]} and \code{corrP[i,j]} for the
#'   correlations above the diagonal only, and \code{h2[i]} or
#'   \code{sigma_P[i]} for the vectors. The trait names are kept in
#'   \code{attr(, "traits")}.
#' @details Correlations are reported once per pair, for \code{i < j}, because
#'   the two halves of the matrix are identical and a draws object has no
#'   notion of symmetry.
#'
#'   To get everything back in a single draws object, bind the elements with
#'   \code{do.call(posterior::bind_draws, unname(draws))}. The names of the
#'   list have to be dropped, because \code{bind_draws()} takes the draws
#'   objects as unnamed arguments.
#' @seealso \code{\link{back_transform}}
#' @examples
#' ss <- size_shape(cbind(a = rnorm(20), b = rnorm(20), c = rnorm(20)))
#' k <- length(ss$raw_traits)
#' m <- cbind(matrix(rep(as.vector(diag(k)), 3), 3, byrow = TRUE),
#'            matrix(rep(as.vector(diag(k)), 3), 3, byrow = TRUE))
#' colnames(m) <- c(sprintf("G[%d,%d]", rep(seq_len(k), k), rep(seq_len(k), each = k)),
#'                  sprintf("E[%d,%d]", rep(seq_len(k), k), rep(seq_len(k), each = k)))
#' draws <- to_draws_list(back_transform(m, ss))
#' names(draws)
#' posterior::summarise_draws(draws$corrG)
#' combined <- do.call(posterior::bind_draws, unname(draws))
#' posterior::variables(combined)[1:4]
#' @export
to_draws_list <- function(x) {
  if (!inherits(x, "size_shape_draws"))
    stop("x must be a size_shape_draws object from back_transform().",
         call. = FALSE)

  components <- setdiff(names(x), "traits")
  out <- lapply(components, function(component) {
    values <- x[[component]]
    if (length(dim(values)) == 3L) {
      k <- dim(values)[1]
      n <- dim(values)[3]
      if (startsWith(component, "corr")) {
        pairs <- which(upper.tri(diag(k)), arr.ind = TRUE)
        pairs <- pairs[order(pairs[, 1], pairs[, 2]), , drop = FALSE]
        index <- cbind(rep(pairs[, 1], each = n), rep(pairs[, 2], each = n),
                       rep(seq_len(n), times = nrow(pairs)))
        draws <- matrix(values[index], n, nrow(pairs))
        variables <- sprintf("%s[%d,%d]", component, pairs[, 1], pairs[, 2])
      } else {
        draws <- t(matrix(values, nrow = k * k))
        variables <- sprintf("%s[%d,%d]", component, rep(seq_len(k), k),
                             rep(seq_len(k), each = k))
      }
    } else {
      draws <- t(values)
      variables <- sprintf("%s[%d]", component, seq_along(x$traits))
    }
    colnames(draws) <- variables
    posterior::as_draws_matrix(draws)
  })
  names(out) <- components
  attr(out, "traits") <- x$traits
  out
}

#' @export
print.size_shape_draws <- function(x, ...) {
  components <- setdiff(names(x), "traits")
  dims <- dim(x[[components[1L]]])
  cat("<size_shape_draws>", dims[length(dims)], "draws,",
      length(x$traits), "traits\n")
  cat("components:", paste(components, collapse = ", "), "\n")
  if (length(x$traits) > 6L) {
    cat("traits:    ", paste(utils::head(x$traits, 6L), collapse = ", "), "...\n")
  } else {
    cat("traits:    ", paste(x$traits, collapse = ", "), "\n")
  }
  invisible(x)
}

# B M B' for every draw of a matrix-valued variable, as trait x trait x draw.
.transform_draws <- function(draws, name, B, traits) {
  k_model <- ncol(B)
  parsed <- array(t(draws[, .matrix_variable_names(name, k_model), drop = FALSE]),
                  c(k_model, k_model, nrow(draws)))
  k <- nrow(B)
  out <- array(apply(parsed, 3L, function(M) B %*% M %*% t(B)), c(k, k, nrow(draws)))
  dimnames(out) <- list(traits, traits, NULL)
  out
}

# cov2cor() of every draw of a trait x trait x draw array.
.corr_draws <- function(A) {
  out <- array(apply(A, 3L, stats::cov2cor), dim(A), dimnames = dimnames(A))
  out
}

# Diagonal of every draw, as a trait x draw matrix.
.diag_draws <- function(A) {
  k <- dim(A)[1]
  n <- dim(A)[3]
  index <- cbind(rep(seq_len(k), times = n), rep(seq_len(k), times = n),
                 rep(seq_len(n), each = k))
  matrix(A[index], k, n, dimnames = list(dimnames(A)[[1]], NULL))
}

# Names of a matrix-valued Stan variable in column-major order, as cmdstanr
# reports them.
.matrix_variable_names <- function(name, k) {
  sprintf("%s[%d,%d]", name, rep(seq_len(k), k), rep(seq_len(k), each = k))
}
