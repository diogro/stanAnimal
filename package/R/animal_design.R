#' Response and design matrices from a model formula
#'
#' Builds the \code{Y} and \code{X} arguments of \code{\link{lmm_animal}},
#' \code{\link{lmm_multiCov}} and \code{\link{lmm_animal_gibbs}} from a formula
#' and a data frame, so a model written the usual way (\code{y ~ sex + age}, or
#' \code{cbind(x, y, z) ~ sex + age}) can be fitted without building the
#' matrices by hand.
#'
#' @param formula Model formula. The left-hand side is a single response or
#'   \code{cbind()} of several traits; the right-hand side lists the fixed
#'   effects.
#' @param data Data frame holding the variables in \code{formula}.
#' @param na.action How to handle missing values, passed to
#'   \code{\link[stats]{model.frame}}; the default drops incomplete rows.
#' @return A list with
#'   \item{Y}{response matrix, one column per trait, named after the
#'     response(s) in the formula;}
#'   \item{X}{fixed-effects design matrix, with columns named after the terms in
#'     the formula.}
#'   Rows are in the order of \code{data}. \code{attr(, "individuals")} holds
#'   the row names of \code{data} for the rows kept, so naming the rows with the
#'   individual IDs lets a relationship matrix \code{A} be subset and reordered
#'   to match.
#' @details Random-effect terms such as \code{(1 | animal)} are ignored, so
#'   formulas written for \code{lme4} can be reused unchanged: relatedness is
#'   supplied through \code{A}, which is passed to the model directly.
#'
#'   The result satisfies what the model functions expect: \code{Y} is numeric
#'   and \code{X} has full column rank. Factor predictors are dummy coded and
#'   the intercept is kept unless the formula removes it, so the columns of
#'   \code{X} match the \code{beta} draws reported by a fit.
#' @seealso \code{\link{lmm_animal}}, \code{\link{animal_draws}}
#' @examples
#' dat <- data.frame(trait = rnorm(6), x = rnorm(6),
#'                   sex = rep(c(0, 1), each = 3), id = 1:6)
#' design <- animal_design(trait ~ sex + x, dat)
#' design$Y
#' design$X
#'
#' # formulas written for lme4 work as they are; the random term is ignored
#' animal_design(trait ~ sex + (1 | id), dat)
#' @export
animal_design <- function(formula, data, na.action = stats::na.omit) {
  if (!inherits(formula, "formula"))
    stop("formula must be a model formula, such as y ~ sex + age.", call. = FALSE)
  if (length(formula) < 3L)
    stop("formula must have a response on the left-hand side, such as y ~ sex.",
         call. = FALSE)
  if (missing(data) || is.null(data))
    stop("data must be a data frame holding the variables in formula.",
         call. = FALSE)

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
  X <- stats::model.matrix(stats::terms(frame), data = frame)
  if (qr(X)$rank < ncol(X))
    stop("the fixed-effects design matrix is not full rank; remove redundant ",
         "or collinear terms from the formula.", call. = FALSE)

  structure(list(Y = Y, X = X),
            individuals = individuals,
            na.action = attr(frame, "na.action"))
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
