#' Genomic relationship matrix from genotype data
#'
#' Computes the VanRaden (2008) genomic relationship matrix (GRM) from biallelic
#' marker genotypes, for use as the \code{A} argument of \code{\link{lmm_animal}}.
#'
#' @param M Numeric genotype matrix with individuals in rows and markers in
#'   columns, coded as the number of copies (0, 1, 2) of the counted allele.
#'   Missing genotypes are allowed (\code{NA}); non-integer values such as
#'   imputed dosages are also accepted.
#' @param p Optional vector of allele frequencies of the counted allele, one per
#'   column of \code{M}. If \code{NULL} (default), they are estimated from
#'   \code{M} as \code{colMeans(M, na.rm = TRUE) / 2}. Supplying base-population
#'   frequencies is recommended when they are known.
#' @param maf Markers with minor allele frequency below this value are dropped.
#'   Monomorphic markers (and, when \code{p} is estimated, markers with a single
#'   observed genotype class) are always dropped.
#' @param max_missing Markers with a proportion of missing genotypes above this
#'   value are dropped.
#' @param blend Weight \code{w} in \code{(1 - w) * G + w * I}. Default 0 returns
#'   the raw GRM; a small value (e.g. 0.01) makes a singular GRM, such as one
#'   with fewer markers than individuals or duplicated genotypes, positive
#'   definite.
#' @return A symmetric positive-semidefinite matrix with the individuals' row
#'   names of \code{M} on both dimensions. The number of markers used and their
#'   scaling constant are stored in the attributes \code{n_markers} and
#'   \code{scale}.
#' @details The GRM is \code{G = Z \%*\% t(Z) / (2 * sum(p * (1 - p)))}, where
#'   \code{Z = M - 2p}. Missing genotypes are replaced by their expected value
#'   \code{2p}, so they contribute zero to the relationship. With frequencies
#'   estimated from the data, the diagonal averages about one plus the
#'   inbreeding coefficient, and the off-diagonals are relative to the
#'   sample average.
#' @references VanRaden, P.M. (2008). Efficient methods to compute genomic
#'   predictions. Journal of Dairy Science 91: 4414-4423.
#' @export
#' @examples
#' set.seed(1)
#' M <- matrix(rbinom(50 * 500, 2, 0.3), nrow = 50)
#' G <- genomic_relationship(M, blend = 0.01)
#' dim(G)
genomic_relationship <- function(M, p = NULL, maf = 0.01, max_missing = 0.1,
                                 blend = 0) {
  M <- .check_numeric_matrix(M, "M", allow_na = TRUE)
  if (nrow(M) < 2L || ncol(M) < 1L)
    stop("`M` must have at least two individuals and one marker.",
         call. = FALSE)
  if (any(M < 0 | M > 2, na.rm = TRUE))
    stop("`M` must have genotype values between 0 and 2.", call. = FALSE)
  maf <- .check_in_range(maf, "maf", 0, 0.5, open_high = TRUE)
  max_missing <- .check_in_range(max_missing, "max_missing", 0, 1)
  blend <- .check_in_range(blend, "blend", 0, 1, open_high = TRUE)

  constant <- rep(FALSE, ncol(M))
  if (is.null(p)) {
    p <- colMeans(M, na.rm = TRUE) / 2
    constant <- apply(M, 2L, function(x) length(unique(x[!is.na(x)])) < 2L)
  } else {
    if (!is.numeric(p) || length(p) != ncol(M) || any(!is.finite(p)) ||
        any(p < 0 | p > 1))
      stop("`p` must be a vector of allele frequencies in [0, 1], one per ",
           "column of `M`.", call. = FALSE)
  }

  missing_rate <- colMeans(is.na(M))
  keep <- is.finite(p) & pmin(p, 1 - p) > maf & pmin(p, 1 - p) > 0 &
    missing_rate <= max_missing & !constant
  if (!any(keep))
    stop("no markers left after filtering: check `maf`, `max_missing` and the ",
         "genotype coding.", call. = FALSE)

  M <- M[, keep, drop = FALSE]
  p <- p[keep]
  Z <- sweep(M, 2L, 2 * p, "-")
  Z[is.na(Z)] <- 0

  scale <- 2 * sum(p * (1 - p))
  G <- tcrossprod(Z) / scale
  G <- (G + t(G)) / 2
  if (blend > 0)
    G <- (1 - blend) * G + blend * diag(nrow(G))

  dimnames(G) <- list(rownames(M), rownames(M))
  attr(G, "n_markers") <- sum(keep)
  attr(G, "scale") <- scale
  G
}
