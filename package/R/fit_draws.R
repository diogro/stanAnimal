# Posterior draws shaped by trait.
#
# The models report covariance matrices one element at a time, as `G[1,1]`,
# `G[2,1]`, ... Turning those into a matrix means matching and reshaping names
# by hand, and every analysis then repeats the same helper. `fit_draws()`
# does it once: it returns trait by trait by draw arrays with named dimensions,
# and `draws_mean()` / `draws_quantile()` summarise one of them in the same
# shape, without the draw dimension.

.model_components <- c("G", "E", "P", "corrG", "corrE", "corrP",
                       "h2", "sigma_G", "sigma_E", "sigma_P", "beta")

.matrix_components <- c("G", "E", "P", "corrG", "corrE", "corrP")

.vector_components <- c("h2", "sigma_G", "sigma_E", "sigma_P")

#' Posterior draws of a fit, shaped by trait
#'
#' Extracts the draws of an animal model and reshapes every parameter so that
#' its first dimensions are the traits. Covariance matrices come out as trait
#' by trait by draw arrays with named dimensions, ready to average, take
#' correlations of, or compare with a simulated truth, instead of one Stan
#' variable per element such as \code{G[1,1]}.
#'
#' @param fit A fit from \code{\link{lmm_animal}}, \code{\link{lm_phen_cov}} or
#'   \code{\link{lmm_animal_gibbs}}, a \code{posterior} draws object, or a named
#'   matrix of draws such as \code{fit$draws(format = "draws_matrix")}. Objects
#'   from \code{\link{back_transform}} also work.
#' @param traits Optional trait names, one per trait. The default uses the
#'   column names of the response when the fit kept them, and
#'   \code{trait1}, \code{trait2}, ... otherwise.
#' @param terms Optional names for the columns of \code{X} in the \code{beta}
#'   draws, in their order in the fit. The default uses the column names of
#'   \code{X} when the fit kept them, and the names the fit reports otherwise.
#' @return A list of class \code{fit_draws}, with one element per quantity
#'   the fit holds: \code{G}, \code{E}, \code{P}, \code{corrG}, \code{corrE},
#'   \code{corrP} as trait by trait by draw arrays, \code{h2}, \code{sigma_G},
#'   \code{sigma_E} and \code{sigma_P} as trait by draw matrices, and
#'   \code{beta} as a trait by term by draw array. Every dimension
#'   is named, and draws are last, the same layout as
#'   \code{\link{lmm_animal_gibbs}} and \code{\link{back_transform}}.
#'   Quantities the fit does not report are derived when they follow from the
#'   others: \code{P = G + E}, \code{corrG}, \code{corrE} and \code{corrP} from
#'   the covariances, \code{sigma_G}, \code{sigma_E} and \code{sigma_P} from
#'   their diagonals, and \code{h2} from \code{G} and \code{P}.
#' @details The Gibbs sampler calls its residual covariance \code{R}; it is
#'   returned as \code{E} so that the two engines give the same names. The
#'   priors differ between the engines, so the draws themselves are not
#'   interchangeable.
#' @seealso \code{\link{draws_mean}}, \code{\link{draws_quantile}},
#'   \code{\link{animal_parameter_draws}}
#' @examples
#' k <- 2
#' vars <- sprintf("G[%d,%d]", rep(1:k, k), rep(1:k, each = k))
#' m <- cbind(matrix(stats::rnorm(40), 10, 4, dimnames = list(NULL, vars)),
#'            matrix(stats::rnorm(40), 10, 4,
#'                   dimnames = list(NULL, sub("G", "E", vars))))
#' d <- fit_draws(m, traits = c("t1", "t2"))
#' names(d)
#' dim(d$G)                       # traits x traits x draws
#' round(draws_mean(d, "G"), 2)   # posterior mean G, as a matrix
#' round(draws_mean(d, "corrG"), 2)
#' @export
fit_draws <- function(fit, traits = NULL, terms = NULL) {
  if (inherits(fit, "size_shape_draws"))
    return(.complete_components(unclass(fit)[setdiff(names(fit), "traits")],
                                traits = traits %||% fit$traits,
                                terms = terms))
  if (.is_gibbs_draws(fit))
    return(.complete_components(
      list(G = fit$G, E = fit$R %||% fit$E, h2 = fit$h2, beta = fit$beta),
      traits = traits %||% dimnames(fit$G)[[1L]], terms = terms))
  m <- .draws_matrix(fit)
  vars <- colnames(m)
  k <- .n_traits_from_draws(vars)
  traits <- .draw_traits(fit, traits, k)
  components <- list(
    G = .matrix_draws(m, "G", k, traits),
    E = .matrix_draws(m, "E", k, traits),
    P = .matrix_draws(m, "P", k, traits),
    corrG = .matrix_draws(m, "corrG", k, traits),
    corrE = .matrix_draws(m, "corrE", k, traits),
    corrP = .matrix_draws(m, "corrP", k, traits),
    h2 = .vector_draws(m, "h2", k, traits),
    sigma_G = .vector_draws(m, "sigma_G", k, traits),
    sigma_E = .vector_draws(m, "sigma_E", k, traits),
    sigma_P = .vector_draws(m, "sigma_P", k, traits),
    beta = .beta_draws(m, k, traits, terms, attr(fit, "x_names"))
  )
  .complete_components(components, traits = traits, terms = terms)
}

#' @export
print.fit_draws <- function(x, ...) {
  cat("<fit_draws>", attr(x, "n_draws"), "draws,",
      length(attr(x, "traits")), "trait(s)\n")
  cat("traits:", paste(attr(x, "traits"), collapse = ", "), "\n")
  for (name in names(x)) {
    dims <- dim(x[[name]])
    cat(sprintf("  %-8s %s\n", name, paste(dims, collapse = " x ")))
  }
  invisible(x)
}

#' Posterior means or quantiles of a fit, as matrices
#'
#' Summarises the draws of an animal model without the draw dimension, so each
#' covariance matrix is returned as a trait by trait matrix and each trait
#' vector as a named vector. These are the shapes used to compare estimates with
#' a simulated truth, compute heritabilities, or look at a correlation matrix.
#'
#' @inheritParams fit_draws
#' @param x A fit, or an object from \code{\link{fit_draws}}.
#' @param parameter Name of the quantity to summarise: \code{"G"}, \code{"E"},
#'   \code{"P"}, \code{"corrG"}, \code{"corrE"}, \code{"corrP"}, \code{"h2"},
#'   \code{"sigma_G"}, \code{"sigma_E"}, \code{"sigma_P"} or \code{"beta"}.
#'   The default summarises every quantity the fit holds.
#' @param probs For \code{draws_quantile}, one or more probabilities between 0
#'   and 1; the default is the 95\% equal-tailed interval.
#' @return For one \code{parameter}, a named matrix or vector: \code{k} by
#'   \code{k} for a covariance or correlation matrix, \code{k} by the number of
#'   terms for \code{beta}, and a vector named after the traits for the others.
#'   For \code{draws_quantile} the last dimension is the quantile, named after
#'   the probability. When \code{parameter} is \code{NULL}, a named list with one
#'   summary per quantity.
#' @details The summaries are computed from the draws, so the posterior mean of
#'   \code{G} is the average of the draws and not the covariance of any average.
#'   \code{draws_quantile} returns the same shapes as \code{draws_mean}, with a
#'   final dimension for the probabilities.
#' @seealso \code{\link{fit_draws}}
#' @examples
#' k <- 2
#' vars <- sprintf("G[%d,%d]", rep(1:k, k), rep(1:k, each = k))
#' m <- cbind(matrix(stats::rnorm(40), 10, 4, dimnames = list(NULL, vars)),
#'            matrix(stats::rnorm(40), 10, 4,
#'                   dimnames = list(NULL, sub("G", "E", vars))))
#' draws_mean(m, "G", traits = c("t1", "t2"))
#' round(draws_mean(m, "G", traits = c("t1", "t2")), 2)
#' draws_quantile(m, "h2", traits = c("t1", "t2"))
#' @export
draws_mean <- function(x, parameter = NULL, traits = NULL) {
  d <- .as_fit_draws(x, traits)
  .summarise_components(d, parameter, .component_mean)
}

#' @rdname draws_mean
#' @export
draws_quantile <- function(x, parameter = NULL, probs = c(0.025, 0.975),
                           traits = NULL) {
  probs <- .check_probs(probs)
  d <- .as_fit_draws(x, traits)
  .summarise_components(d, parameter, function(values)
    .component_quantile(values, probs))
}

.as_fit_draws <- function(x, traits = NULL) {
  if (inherits(x, "fit_draws")) {
    if (is.null(traits)) return(x)
    return(.relabel_fit_draws(x, traits))
  }
  fit_draws(x, traits = traits)
}

.relabel_fit_draws <- function(d, traits) {
  k <- length(attr(d, "traits"))
  if (length(traits) != k)
    stop("`traits` must have one name per trait (", k, ").", call. = FALSE)
  traits <- as.character(traits)
  for (name in intersect(names(d), c(.matrix_components, .vector_components))) {
    dimnames(d[[name]]) <- if (name %in% .matrix_components)
      list(traits, traits, NULL) else list(traits, NULL)
  }
  if (!is.null(d$beta)) dimnames(d$beta)[[1L]] <- traits
  attr(d, "traits") <- traits
  d
}

.summarise_components <- function(d, parameter, fun) {
  wanted <- .component_names(d, parameter)
  out <- lapply(stats::setNames(wanted, wanted), function(name) fun(d[[name]]))
  if (length(out) == 1L && !is.null(parameter)) out[[1L]] else out
}

.component_names <- function(d, parameter) {
  if (is.null(parameter)) return(names(d))
  parameter <- as.character(parameter)
  unknown <- setdiff(parameter, names(d))
  if (length(unknown))
    stop("`parameter` is not in the fit: ", paste(unknown, collapse = ", "),
         ". Available: ", paste(names(d), collapse = ", "), ".", call. = FALSE)
  parameter
}

.component_mean <- function(values) {
  dims <- dim(values)
  if (length(dims) == 3L) {
    out <- apply(values, c(1L, 2L), mean)
    dimnames(out) <- dimnames(values)[1:2]
    out
  } else {
    out <- rowMeans(values)
    names(out) <- rownames(values)
    out
  }
}

.component_quantile <- function(values, probs) {
  labels <- paste0(100 * probs, "%")
  dims <- dim(values)
  if (length(dims) == 3L) {
    out <- array(NA_real_, c(dims[1L], dims[2L], length(probs)),
                 dimnames = list(dimnames(values)[[1L]], dimnames(values)[[2L]],
                                 labels))
    for (q in seq_along(probs))
      out[, , q] <- apply(values, c(1L, 2L), stats::quantile,
                          probs = probs[q], names = FALSE)
    out
  } else {
    q <- apply(values, 1L, stats::quantile, probs = probs, names = FALSE)
    matrix(q, nrow = dims[1L], ncol = length(probs),
           dimnames = list(rownames(values), labels))
  }
}

.check_probs <- function(probs) {
  if (!is.numeric(probs) || !length(probs) || any(!is.finite(probs)) ||
      any(probs < 0 | probs > 1))
    stop("`probs` must be one or more numbers between 0 and 1.", call. = FALSE)
  sort(as.numeric(probs))
}

# The draws as a matrix, with one column per Stan variable, whatever the engine.
.draws_matrix <- function(fit) {
  if (inherits(fit, "draws")) {
    m <- as.matrix(posterior::as_draws_matrix(fit))
  } else if (inherits(fit, "stanfit")) {
    m <- as.matrix(fit)
  } else if (inherits(fit, "CmdStanFit") ||
             (is.list(fit) && is.function(fit$draws))) {
    draws <- fit$draws(format = "draws_matrix")
    m <- if (inherits(draws, "draws"))
      as.matrix(posterior::as_draws_matrix(draws)) else as.matrix(draws)
  } else if (is.matrix(fit)) {
    m <- fit
  } else {
    stop("`fit` must be a fit from `lmm_animal()`, `lm_phen_cov()` or ",
         "`lmm_animal_gibbs()`, or a posterior draws object.", call. = FALSE)
  }
  if (!is.numeric(m) || is.null(colnames(m)))
    stop("`fit` must hold named numeric draws.", call. = FALSE)
  m
}

.is_gibbs_draws <- function(x) {
  is.list(x) && !inherits(x, "draws") && is.null(x$draws) &&
    is.array(x$G) && length(dim(x$G)) == 3L
}

.n_traits_from_draws <- function(vars) {
  pairs <- grep("^[A-Za-z_][A-Za-z0-9_]*\\[[0-9]+,[0-9]+\\]$", vars,
                value = TRUE)
  if (length(pairs)) {
    trait <- as.integer(sub("^[^[]*\\[([0-9]+),.*$", "\\1", pairs))
    return(max(trait))
  }
  vectors <- grep("^(h2|sigma_G|sigma_E|sigma_P)\\[[0-9]+\\]$", vars,
                  value = TRUE)
  if (length(vectors)) {
    trait <- as.integer(sub("^[^[]*\\[([0-9]+)\\]$", "\\1", vectors))
    return(max(trait))
  }
  1L
}

.draw_traits <- function(fit, traits, k) {
  if (!is.null(traits)) {
    if (length(traits) != k)
      stop("`traits` must have one name per trait (", k, ").", call. = FALSE)
    return(as.character(traits))
  }
  y_names <- attr(fit, "y_names")
  if (!is.null(y_names) && length(y_names) == k) return(as.character(y_names))
  paste0("trait", seq_len(k))
}

.matrix_draws <- function(m, name, k, traits) {
  vars <- .matrix_variable_names(name, k)
  if (!all(vars %in% colnames(m))) return(NULL)
  array(t(m[, vars, drop = FALSE]), c(k, k, nrow(m)),
        dimnames = list(traits, traits, NULL))
}

.vector_draws <- function(m, name, k, traits) {
  vars <- if (k == 1L) name else sprintf("%s[%d]", name, seq_len(k))
  if (!all(vars %in% colnames(m))) return(NULL)
  matrix(t(m[, vars, drop = FALSE]), k, nrow(m),
         dimnames = list(traits, NULL))
}

# Beta draws as trait by term by draw. The engines report the fixed effects
# either by index (`beta[2,1]`, `beta[1]`) or already named (`beta[2,age]`,
# `beta[age]`), so both spellings are parsed, and the columns of X name the
# terms when the fit kept them.
.beta_draws <- function(m, k, traits, terms, x_names = NULL) {
  vars <- grep("^beta\\[", colnames(m), value = TRUE)
  if (!length(vars)) return(NULL)
  pieces <- strsplit(sub("^beta\\[", "", sub("\\]$", "", vars)), ",", fixed = TRUE)
  if (length(pieces[[1L]]) == 2L) {
    trait <- vapply(pieces, `[`, character(1), 1L)
    term <- vapply(pieces, `[`, character(1), 2L)
    if (any(!grepl("^[0-9]+$", trait)))
      stop("`beta` draw names must start with the trait index.", call. = FALSE)
    trait <- as.integer(trait)
    levels <- unique(term)
    labels <- .term_labels(levels, terms, x_names)
    out <- array(NA_real_, c(k, length(labels), nrow(m)),
                 dimnames = list(traits, labels, NULL))
    for (v in seq_along(vars))
      out[trait[v], match(term[v], levels), ] <- m[, vars[v]]
    out
  } else {
    labels <- .term_labels(vapply(pieces, `[`, character(1), 1L), terms,
                           x_names)
    out <- array(NA_real_, c(1L, length(labels), nrow(m)),
                 dimnames = list(traits[1L], labels, NULL))
    out[1L, , ] <- t(m[, vars, drop = FALSE])
    out
  }
}

# Term labels in the order the fit reports them, with the columns of X taking
# over when the fit reports indices.
.term_labels <- function(reported, terms, x_names) {
  levels <- unique(reported)
  if (!is.null(terms)) {
    if (length(terms) != length(levels))
      stop("`terms` must have one name per column of `X` (", length(levels),
           ").", call. = FALSE)
    return(as.character(terms))
  }
  if (!is.null(x_names) && all(grepl("^[0-9]+$", levels))) {
    index <- as.integer(levels)
    if (length(x_names) < max(index))
      stop("`terms` was not given and the column names of `X` do not match ",
           "the fixed effects in the fit.", call. = FALSE)
    return(as.character(x_names[index]))
  }
  levels
}

.components_n_traits <- function(components) {
  for (name in intersect(c("G", "E", "P", "corrG", "beta"),
                         names(components)))
    return(dim(components[[name]])[1L])
  for (name in intersect(.vector_components, names(components)))
    return(dim(components[[name]])[1L])
  1L
}

# Bring any set of components to the documented layout and derive what follows
# from what is there.
.complete_components <- function(components, traits = NULL, terms = NULL) {
  components <- components[!vapply(components, is.null, logical(1L))]
  if (!length(components))
    stop("`fit` does not hold the draws of an animal model: expected G and E, ",
         "P, or h2.", call. = FALSE)
  k <- .components_n_traits(components)
  traits <- traits %||% paste0("trait", seq_len(k))
  if (length(traits) != k)
    stop("`traits` must have one name per trait (", k, ").", call. = FALSE)
  traits <- as.character(traits)

  for (name in intersect(names(components), .matrix_components))
    dimnames(components[[name]]) <- list(traits, traits, NULL)
  for (name in intersect(names(components), .vector_components))
    dimnames(components[[name]]) <- list(traits, NULL)
  if (!is.null(components$beta)) {
    if (!is.null(terms)) {
      if (dim(components$beta)[2L] != length(terms))
        stop("`terms` must have one name per column of `X` (",
             dim(components$beta)[2L], ").", call. = FALSE)
      dimnames(components$beta)[[2L]] <- as.character(terms)
    }
    dimnames(components$beta)[[1L]] <- traits
  }

  if (is.null(components$P) && !is.null(components$G) && !is.null(components$E))
    components$P <- components$G + components$E
  if (is.null(components$corrG) && !is.null(components$G))
    components$corrG <- .corr_draws(components$G)
  if (is.null(components$corrE) && !is.null(components$E))
    components$corrE <- .corr_draws(components$E)
  if (is.null(components$corrP) && !is.null(components$P))
    components$corrP <- .corr_draws(components$P)
  if (is.null(components$sigma_G) && !is.null(components$G))
    components$sigma_G <- sqrt(.diag_draws(components$G))
  if (is.null(components$sigma_E) && !is.null(components$E))
    components$sigma_E <- sqrt(.diag_draws(components$E))
  if (is.null(components$sigma_P) && !is.null(components$P))
    components$sigma_P <- sqrt(.diag_draws(components$P))
  if (is.null(components$h2) && !is.null(components$G) && !is.null(components$P))
    components$h2 <- .diag_draws(components$G) / .diag_draws(components$P)

  components <- components[intersect(.model_components, names(components))]
  n_draws <- dim(components[[1L]])[length(dim(components[[1L]]))]
  structure(components, class = "fit_draws", traits = traits,
            terms = if (!is.null(components$beta)) dimnames(components$beta)[[2L]],
            n_draws = n_draws)
}
