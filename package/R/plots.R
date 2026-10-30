#' Prior draws for heritabilities and genetic correlations
#'
#' Samples the priors that \code{\link{lmm_animal}} places on the genetic and
#' residual covariance matrices and summarises them as heritabilities and
#' genetic correlations, which do not depend on the scale of the data.
#'
#' @param n_traits Number of traits.
#' @param n Number of prior draws.
#' @param lkj_prior LKJ shape used in the fit.
#' @param unstructured_E If \code{FALSE}, the residual covariance is diagonal.
#' @param trait_names Optional trait names.
#' @return A long data frame with columns \code{type} (\code{"Heritability"} or
#'   \code{"Genetic correlation"}), \code{parameter} and \code{value}.
#' @details Standard deviations are half-normal(0, 1) and correlation matrices
#'   are LKJ(\code{lkj_prior}), as in the Stan models. Univariate fits use the
#'   default \code{parametrization = "variance"} prior.
#' @export
animal_prior_draws <- function(n_traits, n = 4000, lkj_prior = 2,
                               unstructured_E = TRUE, trait_names = NULL) {
  trait_names <- .trait_names(trait_names, n_traits)
  h2 <- matrix(NA_real_, n, n_traits, dimnames = list(NULL, trait_names))
  pairs <- .trait_pairs(trait_names)
  cors <- matrix(NA_real_, n, nrow(pairs), dimnames = list(NULL, pairs$parameter))
  for (i in seq_len(n)) {
    sG <- abs(stats::rnorm(n_traits))
    sE <- abs(stats::rnorm(n_traits))
    Og <- .rlkj(n_traits, lkj_prior)
    Gv <- sG^2
    h2[i, ] <- Gv / (Gv + sE^2)
    if (n_traits > 1L)
      cors[i, ] <- Og[cbind(pairs$i, pairs$j)]
  }
  .long_draws(h2, if (n_traits > 1L) cors, "prior")
}

#' Heritability and genetic correlation draws from a fit
#'
#' @param fit A fit returned by \code{lmm_animal} with \code{cmd_stan = TRUE}, or
#'   a \code{posterior} draws object containing \code{h2} and, for several
#'   traits, \code{corrG}.
#' @param trait_names Optional trait names.
#' @return A long data frame with columns \code{type}, \code{parameter} and
#'   \code{value}.
#' @export
animal_parameter_draws <- function(fit, trait_names = NULL) {
  draws <- fit
  if (inherits(fit, "CmdStanFit") || (is.list(fit) && is.function(fit$draws))) {
    vars <- if (is.function(fit$metadata)) fit$metadata()$stan_variables else c("h2", "corrG")
    draws <- fit$draws(intersect(c("h2", "corrG"), vars))
  }
  df <- as.data.frame(posterior::as_draws_df(draws))
  h2_cols <- grep("^h2(\\[[0-9]+\\])?$", names(df), value = TRUE)
  if (!length(h2_cols)) stop("No h2 draws found.")
  trait_names <- .trait_names(trait_names, length(h2_cols))
  h2 <- as.matrix(df[h2_cols]); colnames(h2) <- trait_names
  cors <- NULL
  if (length(h2_cols) > 1L) {
    pairs <- .trait_pairs(trait_names)
    src <- sprintf("corrG[%d,%d]", pairs$i, pairs$j)
    if (!all(src %in% names(df))) stop("No corrG draws found.")
    cors <- as.matrix(df[src]); colnames(cors) <- pairs$parameter
  }
  .long_draws(h2, cors, "posterior")
}

#' Posterior versus prior for heritabilities and genetic correlations
#'
#' @inheritParams animal_parameter_draws
#' @param lkj_prior,unstructured_E,n_prior Passed to \code{\link{animal_prior_draws}}
#'   (\code{n_prior} is \code{n}).
#' @param prior Optional prior draws in the format returned by
#'   \code{animal_prior_draws}; computed if \code{NULL}.
#' @return A ggplot object. \code{plot_prior_posterior} shows overlaid densities
#'   in one panel per parameter, and \code{plot_forest} shows means and
#'   \code{prob} intervals of the posterior and the prior.
#' @details Where prior and posterior overlap, the data add little information.
#' @export
plot_prior_posterior <- function(fit, trait_names = NULL, lkj_prior = 2,
                                 unstructured_E = TRUE, n_prior = 4000,
                                 prior = NULL, ncol = 4) {
  d <- .prior_posterior(fit, trait_names, lkj_prior, unstructured_E, n_prior, prior)
  ggplot2::ggplot(d, ggplot2::aes(.data$value, fill = .data$source, colour = .data$source)) +
    ggplot2::geom_density(alpha = 0.35, adjust = 1.2) +
    ggplot2::facet_wrap(ggplot2::vars(.data$parameter), scales = "free_y", ncol = ncol) +
    ggplot2::scale_fill_manual(values = .source_colours$fill) +
    ggplot2::scale_colour_manual(values = .source_colours$line) +
    ggplot2::labs(x = NULL, y = "Density", fill = NULL, colour = NULL,
                  title = "Posterior versus prior") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(legend.position = "top", axis.text.y = ggplot2::element_blank())
}

#' @rdname plot_prior_posterior
#' @param prob Probability mass of the highest density interval.
#' @param ncol Number of facet columns for \code{plot_prior_posterior}.
#' @export
plot_forest <- function(fit, trait_names = NULL, lkj_prior = 2,
                        unstructured_E = TRUE, n_prior = 4000,
                        prior = NULL, prob = 0.95) {
  d <- .prior_posterior(fit, trait_names, lkj_prior, unstructured_E, n_prior, prior)
  groups <- split(d, list(d$type, d$parameter, d$source), drop = TRUE)
  s <- do.call(rbind, lapply(groups, function(g) {
    hdi <- .hdi(g$value, prob)
    data.frame(type = g$type[1], parameter = g$parameter[1], source = g$source[1],
               mean = mean(g$value), lower = hdi[1], upper = hdi[2])
  }))
  s$parameter <- factor(s$parameter, levels = rev(levels(d$parameter)))
  ggplot2::ggplot(s, ggplot2::aes(.data$mean, .data$parameter, colour = .data$source)) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
    ggplot2::geom_pointrange(ggplot2::aes(xmin = .data$lower, xmax = .data$upper),
                             position = ggplot2::position_dodge(width = 0.6)) +
    ggplot2::facet_wrap(ggplot2::vars(.data$type), scales = "free", ncol = 2) +
    ggplot2::scale_colour_manual(values = .source_colours$line) +
    ggplot2::labs(x = sprintf("Estimate and %g%% HPD interval", 100 * prob),
                  y = NULL, colour = NULL,
                  title = "Posterior and prior intervals") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(legend.position = "top")
}

.source_colours <- list(fill = c(Posterior = "firebrick", Prior = "grey60"),
                        line = c(Posterior = "firebrick", Prior = "grey40"))

.prior_posterior <- function(fit, trait_names, lkj_prior, unstructured_E,
                             n_prior, prior) {
  post <- animal_parameter_draws(fit, trait_names)
  if (is.null(prior)) {
    k <- length(unique(post$parameter[post$type == "Heritability"]))
    prior <- animal_prior_draws(k, n_prior, lkj_prior, unstructured_E,
                                unique(post$parameter[post$type == "Heritability"]))
  }
  d <- rbind(post, prior)
  d$source <- factor(ifelse(d$source == "posterior", "Posterior", "Prior"),
                     levels = c("Posterior", "Prior"))
  d$type <- factor(d$type, levels = c("Heritability", "Genetic correlation"))
  d$parameter <- factor(d$parameter, levels = unique(post$parameter))
  d[order(d$type), ]
}

.trait_names <- function(trait_names, k) {
  if (is.null(trait_names)) return(if (k == 1L) "h2" else paste0("trait", seq_len(k)))
  if (length(trait_names) != k) stop("trait_names must have one name per trait.")
  as.character(trait_names)
}

.trait_pairs <- function(trait_names) {
  k <- length(trait_names)
  idx <- which(upper.tri(diag(k)), arr.ind = TRUE)
  idx <- idx[order(idx[, 1], idx[, 2]), , drop = FALSE]
  data.frame(i = idx[, 1], j = idx[, 2],
             parameter = paste(trait_names[idx[, 1]], trait_names[idx[, 2]], sep = " ~ "))
}

.long_draws <- function(h2, cors, source) {
  out <- data.frame(type = "Heritability",
                    parameter = rep(colnames(h2), each = nrow(h2)),
                    value = as.vector(h2))
  if (!is.null(cors))
    out <- rbind(out, data.frame(type = "Genetic correlation",
                                 parameter = rep(colnames(cors), each = nrow(cors)),
                                 value = as.vector(cors)))
  out$source <- source
  out
}

# Shortest interval containing `prob` of the draws.
.hdi <- function(x, prob) {
  x <- sort(x)
  n <- length(x)
  m <- max(1L, floor(prob * n))
  w <- x[seq(m + 1L, n)] - x[seq_len(n - m)]
  if (!length(w)) return(range(x))
  i <- which.min(w)
  c(x[i], x[i + m])
}

# LKJ correlation matrix via the onion-free vine method (Lewandowski et al. 2009).
.rlkj <- function(k, eta) {
  if (k == 1L) return(matrix(1, 1, 1))
  P <- matrix(0, k, k)
  S <- diag(k)
  for (m in seq_len(k - 1L)) {
    beta <- eta + (k - 1L - m) / 2
    for (i in seq(m + 1L, k)) {
      P[m, i] <- 2 * stats::rbeta(1, beta, beta) - 1
      p <- P[m, i]
      if (m > 1L)
        for (l in rev(seq_len(m - 1L)))
          p <- p * sqrt((1 - P[l, i]^2) * (1 - P[l, m]^2)) + P[l, i] * P[l, m]
      S[m, i] <- S[i, m] <- p
    }
  }
  S
}

#' Pairs plot of the data with posterior draws of G ellipses
#'
#' For each pair of traits, the lower triangle shows the data with the
#' \code{level} ellipses of a sample of posterior G draws (thin lines), the
#' posterior mean G (thick line) and the observed phenotypic covariance
#' (dashed), all centred on the data mean. The upper triangle reports the
#' posterior mean genetic correlation and its HPD interval.
#'
#' @param fit A multivariate fit from \code{lmm_animal} (\code{cmd_stan = TRUE})
#'   or a \code{posterior} draws object containing \code{G} and \code{corrG}.
#' @param Y Matrix of the observed traits (the \code{Y} given to the fit).
#' @param trait_names Optional trait names; default are the column names of
#'   \code{Y}.
#' @param n_show Number of posterior draws of G to draw.
#' @param level Probability level of the ellipses.
#' @param prob Probability mass of the HPD interval of the correlations.
#' @return A ggplot object.
#' @export
plot_G_ellipses <- function(fit, Y, trait_names = NULL, n_show = 100,
                            level = 0.95, prob = 0.95) {
  Y <- as.matrix(Y)
  k <- ncol(Y)
  if (k < 2L) stop("At least two traits are needed.")
  trait_names <- trait_names %||% colnames(Y) %||% paste0("trait", seq_len(k))
  if (length(trait_names) != k) stop("trait_names must have one name per trait.")

  draws <- fit
  if (inherits(fit, "CmdStanFit") || (is.list(fit) && is.function(fit$draws)))
    draws <- fit$draws("G")
  df <- as.data.frame(posterior::as_draws_df(draws))
  g_names <- outer(seq_len(k), seq_len(k), function(i, j) sprintf("G[%d,%d]", i, j))
  if (!all(g_names %in% names(df))) stop("G draws for all traits not found.")
  G <- array(unlist(df[as.vector(g_names)]), c(nrow(df), k, k))

  mean_G <- apply(G, c(2, 3), mean)
  corr <- t(apply(G, 1, function(g) stats::cov2cor(g)[upper.tri(diag(k))]))
  if (!is.matrix(corr)) corr <- matrix(corr, ncol = k * (k - 1) / 2)
  pairs <- .trait_pairs(trait_names)
  shown <- sample(dim(G)[1], min(n_show, dim(G)[1]))
  P_obs <- stats::cov(Y)
  centre <- colMeans(Y)
  lab <- factor(trait_names, levels = trait_names)

  lower <- list(points = list(), lines = list())
  upper <- list(); diag_l <- list(); ranges <- list()
  for (c_ in seq_len(k)) {
    r <- range(Y[, c_])
    ranges[[c_]] <- r + c(-1, 1) * 0.35 * diff(r)
  }
  for (row in seq_len(k)) for (col in seq_len(k)) {
    row_f <- lab[row]; col_f <- lab[col]
    mid <- c(mean(ranges[[col]]), mean(ranges[[row]]))
    if (row == col) {
      diag_l[[length(diag_l) + 1L]] <- data.frame(
        row = row_f, col = col_f, x = mid[1], y = mid[2], label = trait_names[row])
    } else if (row < col) {
      m <- which(pairs$i == row & pairs$j == col)
      hdi <- .hdi(corr[, m], prob)
      upper[[length(upper) + 1L]] <- data.frame(
        row = row_f, col = col_f, x = mid[1], y = mid[2],
        label = sprintf("genetic correlation\n%.2f\n%g%% HPD [%.2f, %.2f]",
                        mean(corr[, m]), 100 * prob, hdi[1], hdi[2]))
    } else {
      idx <- c(col, row)
      mk <- function(S, group, draw = 1L)
        data.frame(row = row_f, col = col_f, group = group, draw = draw,
                   .ellipse(S[idx, idx], centre[idx], level))
      lower$points[[length(lower$points) + 1L]] <- data.frame(
        row = row_f, col = col_f, x = Y[, col], y = Y[, row])
      lower$lines <- c(lower$lines,
        lapply(shown, function(s) mk(G[s, , ], "G draws", s)),
        list(mk(mean_G, "Posterior mean G"), mk(P_obs, "Observed P")))
    }
  }
  bind <- function(l) do.call(rbind, l)
  points <- bind(lower$points); lines <- bind(lower$lines)
  lines$id <- interaction(lines$row, lines$col, lines$group, lines$draw)
  txt <- rbind(cbind(bind(upper), size = 3.2, face = "plain"),
               cbind(bind(diag_l), size = 5, face = "bold"))
  # Fix the axis ranges of each trait so every panel in a row/column shares them.
  blank <- do.call(rbind, lapply(seq_len(k), function(i) do.call(rbind, lapply(
    seq_len(k), function(j) data.frame(
      row = lab[i], col = lab[j],
      x = ranges[[j]], y = ranges[[i]])))))

  ggplot2::ggplot() +
    ggplot2::geom_blank(data = blank, ggplot2::aes(.data$x, .data$y)) +
    ggplot2::geom_point(data = points, ggplot2::aes(.data$x, .data$y),
                        colour = "grey30", size = 0.8) +
    ggplot2::geom_path(data = lines[lines$group == "G draws", ],
                       ggplot2::aes(.data$x, .data$y, group = .data$id),
                       colour = "firebrick", alpha = 0.12) +
    ggplot2::geom_path(data = lines[lines$group != "G draws", ],
                       ggplot2::aes(.data$x, .data$y, group = .data$id,
                                    colour = .data$group, linetype = .data$group),
                       linewidth = 1) +
    ggplot2::geom_text(data = txt, ggplot2::aes(.data$x, .data$y, label = .data$label,
                                                size = I(.data$size), fontface = .data$face)) +
    ggplot2::facet_grid(rows = ggplot2::vars(.data$row), cols = ggplot2::vars(.data$col),
                        scales = "free") +
    ggplot2::scale_colour_manual(values = c("Posterior mean G" = "firebrick",
                                            "Observed P" = "royalblue")) +
    ggplot2::scale_linetype_manual(values = c("Posterior mean G" = "solid",
                                              "Observed P" = "dashed")) +
    ggplot2::labs(x = NULL, y = NULL, colour = NULL, linetype = NULL) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(legend.position = "top", panel.border = ggplot2::element_rect(
      colour = "grey70", fill = NA))
}

.ellipse <- function(S, centre, level = 0.95) {
  e <- eigen(S, symmetric = TRUE)
  a <- seq(0, 2 * pi, length.out = 100)
  circle <- cbind(cos(a), sin(a)) * sqrt(stats::qchisq(level, 2))
  out <- sweep(circle %*% diag(sqrt(pmax(e$values, 0))) %*% t(e$vectors), 2, centre, "+")
  data.frame(x = out[, 1], y = out[, 2])
}

`%||%` <- function(a, b) if (is.null(a)) b else a
