#' Response and design matrices from a model formula
#'
#' Builds the \code{Y} and \code{X} arguments of \code{\link{lmm_animal}},
#' \code{\link{lm_phen_cov}} and \code{\link{lmm_animal_gibbs}} from a formula
#' and a data frame, so a model written the usual way (\code{y ~ sex + age}, or
#' \code{cbind(x, y, z) ~ sex + age}) can be fitted without building the
#' matrices by hand.
#'
#' @param formula Model formula. The left-hand side is a single response or
#'   \code{cbind()} of several traits; the right-hand side lists the fixed
#'   effects.
#' @param data Data frame holding the variables in \code{formula}. When
#'   \code{pedigree} or \code{genotypes} are used, its row names must be the
#'   individual IDs.
#' @param na.action How to handle missing values, passed to
#'   \code{\link[stats]{model.frame}}; the default drops incomplete rows.
#' @param pedigree Optional pedigree, used to compute the additive relationship
#'   matrix \code{A} with \code{\link{makeA}}. A data frame
#'   whose first three columns are the individual, its sire and its dam, with
#'   \code{0}, \code{NA} or \code{""} for unknown parents. The pedigree may
#'   contain ancestors that are not in \code{data}.
#' @param genotypes Optional genotype matrix, used to compute a genomic
#'   relationship matrix (see \code{\link{genomic_relationship}}): individuals
#'   in rows, with their IDs as row names, and biallelic markers in columns,
#'   coded 0, 1, 2.
#' @param ... Further arguments passed to \code{\link{genomic_relationship}},
#'   such as \code{p}, \code{maf} or \code{blend}.
#' @return A list with
#'   \item{Y}{response matrix, one column per trait, named after the
#'     response(s) in the formula;}
#'   \item{X}{fixed-effects design matrix, with columns named after the terms in
#'     the formula;}
#'   \item{A}{relationship matrix, present only when \code{pedigree} or
#'     \code{genotypes} is given, restricted to the individuals in \code{Y} and
#'     in the same row and column order.}
#'   Rows are in the order of \code{data}. \code{attr(, "individuals")} holds
#'   the row names of \code{data} for the rows kept, so naming the rows with the
#'   individual IDs lets a relationship matrix \code{A} be matched to
#'   \code{Y}: \code{\link{lmm_animal}} puts a named \code{A} in the row order of
#'   \code{Y} and errors if the names disagree.
#' @details Random-effect terms such as \code{(1 | animal)} are ignored, so
#'   formulas written for \code{lme4} can be reused unchanged: relatedness is
#'   supplied through \code{A}, which is passed to the model directly.
#'
#'   The result satisfies what the model functions expect: \code{Y} is numeric
#'   and \code{X} has full column rank. Factor predictors are dummy coded and
#'   the intercept is kept unless the formula removes it, so the columns of
#'   \code{X} match the \code{beta} draws reported by a fit. Removing the
#'   intercept is only allowed when the remaining columns still sum to a
#'   constant, as the 0/1 indicators of a factor do, because the models centre
#'   the responses and report their mean back through the design; a formula such
#'   as \code{y ~ 0 + x} is an error.
#'
#'   When \code{pedigree} or \code{genotypes} is given, \code{A} is computed for
#'   all individuals in the pedigree or genotype matrix and then subset to the
#'   rows of \code{Y}, so ancestors without records and individuals with missing
#'   phenotypes are handled: every row of \code{Y} must appear in the pedigree
#'   or in the row names of the genotype matrix.
#'
#'   If the response has three or more traits and its first principal component
#'   explains at least 80\% of the variance, a warning suggests fitting size and
#'   shape scores instead (see \code{\link{size_shape}});
#'   \code{\link{pc1_variance_share}} reports that share.
#' @seealso \code{\link{lmm_animal}}, \code{\link{makeA}},
#'   \code{\link{genomic_relationship}}, \code{\link{rename_beta}}
#' @examples
#' dat <- data.frame(trait = rnorm(6), x = rnorm(6),
#'                   sex = rep(c(0, 1), each = 3), id = 1:6)
#' rownames(dat) <- dat$id
#' design <- animal_design(trait ~ sex + x, dat)
#' design$Y
#' design$X
#'
#' # formulas written for lme4 work as they are; the random term is ignored
#' animal_design(trait ~ sex + (1 | id), dat)
#'
#' # relatedness from a pedigree (first three columns: id, sire, dam)
#' ped <- data.frame(id = 1:6, sire = c(0, 0, 1, 1, 4, 4),
#'                   dam = c(0, 0, 2, 2, 3, 5))
#' design <- animal_design(trait ~ sex, dat, pedigree = ped)
#' round(design$A, 2)
#'
#' # relatedness from genotypes (individuals in rows, markers in columns)
#' M <- matrix(rbinom(6 * 40, 2, 0.3), 6,
#'             dimnames = list(as.character(1:6), NULL))
#' design <- animal_design(trait ~ sex, dat, genotypes = M, blend = 0.01)
#' dim(design$A)
#' @export
animal_design <- function(formula, data, na.action = stats::na.omit,
                          pedigree = NULL, genotypes = NULL, ...) {
  if (!inherits(formula, "formula"))
    stop("formula must be a model formula, such as y ~ sex + age.", call. = FALSE)
  if (length(formula) < 3L)
    stop("formula must have a response on the left-hand side, such as y ~ sex.",
         call. = FALSE)
  if (missing(data) || is.null(data))
    stop("data must be a data frame holding the variables in formula.",
         call. = FALSE)
  if (!is.null(pedigree) && !is.null(genotypes))
    stop("pass either pedigree or genotypes, not both.", call. = FALSE)

  formula <- .drop_random_terms(formula)
  frame <- stats::model.frame(formula, data = data, na.action = na.action,
                              drop.unused.levels = TRUE)
  response <- stats::model.response(frame)
  if (is.null(response))
    stop("formula must have a response on the left-hand side.", call. = FALSE)

  individuals <- rownames(frame)
  if (is.null(dim(response))) {
    Y <- matrix(response, ncol = 1L,
                dimnames = list(individuals,
                                deparse(formula[[2L]], width.cutoff = 500L)))
  } else {
    Y <- as.matrix(response)
    rownames(Y) <- individuals
    if (is.null(colnames(Y)))
      colnames(Y) <- paste0("trait", seq_len(ncol(Y)))
  }
  if (!is.numeric(Y))
    stop("the response must be numeric: a single trait or cbind() of traits.",
         call. = FALSE)
  if (any(!is.finite(Y)))
    stop("Y must contain only finite numeric values.", call. = FALSE)

  if (nrow(frame) < 2L)
    stop("at least two complete rows are needed to build Y and X.", call. = FALSE)

  if (ncol(Y) >= 3L && nrow(Y) >= 3L) {
    pc1_share <- pc1_variance_share(Y)
    if (pc1_share >= 0.8)
      warning("the first principal component of Y explains ",
              round(100 * pc1_share), "% of the variance. If these are raw ",
              "measurements, consider fitting size and shape scores instead ",
              "(see ?size_shape).", call. = FALSE)
  }

  X <- stats::model.matrix(stats::terms(frame), data = frame)
  if (qr(X)$rank < ncol(X))
    stop("the fixed-effects design matrix is not full rank; remove redundant ",
         "or collinear terms from the formula.", call. = FALSE)
  if (is.null(.intercept_weights(X)))
    stop("the fixed-effects design matrix cannot represent a constant, so the ",
         "mean of the response cannot be reported: keep the intercept, or the ",
         "0/1 indicator columns of every level of a factor.", call. = FALSE)

  out <- list(Y = Y, X = X)
  if (!is.null(pedigree) || !is.null(genotypes)) {
    if (is.null(individuals))
      stop("with pedigree or genotypes, the rows of data must be named with ",
           "the individual IDs.", call. = FALSE)
    if (!is.null(pedigree)) {
      A <- makeA(pedigree)
      source <- "pedigree"
    } else {
      A <- genomic_relationship(genotypes, ...)
      source <- "genotypes"
    }
    out$A <- .check_relationship_matrix(
      .subset_relationship_matrix(A, individuals, source), Y)
  } else if (length(list(...))) {
    stop("... is only used when genotypes is given.", call. = FALSE)
  }

  structure(out,
            individuals = individuals,
            na.action = attr(frame, "na.action"))
}

# Keeps the individuals of Y, in the order of Y, from a named relationship
# matrix that may cover a larger set (ancestors, individuals without records).
.subset_relationship_matrix <- function(A, ids, source) {
  A_ids <- rownames(A)
  if (is.null(A_ids))
    stop("the ", source, " must carry individual identifiers as row names to ",
         "match the rows of Y.", call. = FALSE)
  missing <- setdiff(ids, A_ids)
  if (length(missing))
    stop(length(missing), " rows of Y are not in the ", source, " (e.g. ",
         paste(utils::head(missing, 3L), collapse = ", "), "). Set the ",
         "individual IDs as row names of data with rownames(data) <- data$id.",
         call. = FALSE)
  out <- A[ids, ids, drop = FALSE]
  # keep informational attributes such as n_markers and scale
  extra <- setdiff(names(attributes(A)), c("dim", "dimnames"))
  for (attribute in extra)
    attr(out, attribute) <- attr(A, attribute)
  out
}

# Drops random-effect terms, e.g. (1 | animal), so lme4-style formulas can be
# reused; relatedness is passed to the model as A instead.
.drop_random_terms <- function(formula) {
  term_labels <- attr(stats::terms(formula), "term.labels")
  keep <- !grepl("|", term_labels, fixed = TRUE)
  if (all(keep))
    return(formula)
  stats::reformulate(term_labels[keep],
                     response = deparse(formula[[2L]], width.cutoff = 500L),
                     intercept = attr(stats::terms(formula), "intercept") == 1L,
                     env = environment(formula) %||% parent.frame())
}
