# A fit is a list with a draws() method plus the names of X and of the traits,
# the shape cmdstanr returns. The draws are built by hand, so the tests check
# the reshaping against known values instead of against a sampler.
named_columns <- function(x, names) {
  colnames(x) <- names
  x
}

matrix_names <- function(name, k) {
  sprintf(paste0(name, "[%d,%d]"), rep(seq_len(k), k), rep(seq_len(k), each = k))
}

mock_fit <- function(k = 2, n = 20, terms = c("Int", "x"), seed = 1,
                     traits = paste0("trait", seq_len(k))) {
  set.seed(seed)
  G <- E <- array(0, c(k, k, n))
  for (i in seq_len(n)) {
    G[, , i] <- crossprod(matrix(stats::rnorm(k * k), k)) + diag(k)
    E[, , i] <- (crossprod(matrix(stats::rnorm(k * k), k)) + diag(k)) * 2
  }
  P <- G + E
  corrG <- array(apply(G, 3L, stats::cov2cor), dim(G))
  beta <- array(stats::rnorm(k * length(terms) * n),
                c(k, length(terms), n))
  by_draw <- function(A) t(apply(A, 3L, as.vector))
  diagonal <- function(A, fun) t(apply(A, 3L, function(M) fun(diag(M))))
  h2 <- t(sapply(seq_len(n), function(i) diag(G[, , i]) / diag(P[, , i])))

  m <- cbind(
    named_columns(by_draw(G), matrix_names("G", k)),
    named_columns(by_draw(E), matrix_names("E", k)),
    named_columns(by_draw(P), matrix_names("P", k)),
    named_columns(by_draw(corrG), matrix_names("corrG", k)),
    named_columns(diagonal(G, sqrt), sprintf("sigma_G[%d]", seq_len(k))),
    named_columns(diagonal(E, sqrt), sprintf("sigma_E[%d]", seq_len(k))),
    named_columns(t(apply(beta, 3L, as.vector)),
                  sprintf("beta[%d,%d]", rep(seq_len(k), length(terms)),
                          rep(seq_along(terms), each = k))),
    named_columns(h2, sprintf("h2[%d]", seq_len(k)))
  )
  structure(list(draws = function(variables = NULL, format = NULL, ...)
    posterior::as_draws_matrix(m)),
    class = "CmdStanMCMC", x_names = terms, y_names = traits,
    truth = list(G = G, E = E, P = P, corrG = corrG, beta = beta, m = m))
}

test_that("fit_draws reshapes covariance and vector draws by trait", {
  fit <- mock_fit()
  d <- fit_draws(fit)
  expect_s3_class(d, "fit_draws")
  expect_equal(names(d), c("G", "E", "P", "corrG", "corrE", "corrP", "h2",
                           "sigma_G", "sigma_E", "sigma_P", "beta"))
  expect_equal(dim(d$G), c(2L, 2L, 20L))
  expect_equal(dimnames(d$G), list(c("trait1", "trait2"),
                                   c("trait1", "trait2"), NULL))
  expect_equal(dim(d$h2), c(2L, 20L))
  expect_equal(rownames(d$h2), c("trait1", "trait2"))
  expect_equal(dim(d$beta), c(2L, 2L, 20L))
  expect_equal(dimnames(d$beta), list(c("trait1", "trait2"), c("Int", "x"), NULL))

  truth <- attr(fit, "truth")
  expect_equal(unname(d$G), truth$G)
  expect_equal(unname(d$P), truth$P)
  expect_equal(unname(d$corrG), truth$corrG)
  expect_equal(d$corrE[, , 3], stats::cov2cor(matrix(d$E[, , 3], 2)),
               ignore_attr = TRUE)
  expect_equal(unname(d$h2[, 1]), diag(truth$G[, , 1]) / diag(truth$P[, , 1]))
  expect_equal(unname(d$sigma_G[, 1]), sqrt(diag(truth$G[, , 1])))
  expect_equal(unname(d$beta[, , 1]), unname(truth$beta[, , 1]))
})

test_that("draws_mean and draws_quantile drop the draw dimension", {
  d <- fit_draws(mock_fit(), traits = c("a", "b"))
  expect_equal(dim(draws_mean(d, "G")), c(2L, 2L))
  expect_equal(dimnames(draws_mean(d, "G")), list(c("a", "b"), c("a", "b")))
  expect_equal(draws_mean(d, "G")[1, 2], mean(d$G[1, 2, ]))
  expect_equal(unname(draws_mean(d, "h2")), unname(apply(d$h2, 1L, mean)))
  expect_equal(dim(draws_mean(d, "beta")), c(2L, 2L))
  all_means <- draws_mean(d)
  expect_equal(names(all_means), names(d))

  q <- draws_quantile(d, "G", probs = c(0.1, 0.9))
  expect_equal(dim(q), c(2L, 2L, 2L))
  expect_equal(dimnames(q)[[3]], c("10%", "90%"))
  expect_equal(q[1, 1, 1], unname(stats::quantile(d$G[1, 1, ], 0.1)))
  q1 <- draws_quantile(d, "h2", probs = 0.5)
  expect_equal(dim(q1), c(2L, 1L))
  expect_equal(q1[, 1], apply(d$h2, 1L, stats::median))
})

test_that("fit_draws validates its inputs", {
  d <- fit_draws(mock_fit())
  expect_error(draws_mean(d, "R"), "`parameter`")
  expect_error(draws_quantile(d, "G", probs = 2), "`probs`")
  expect_error(fit_draws(d, traits = "one name"), "one name per trait")
  expect_error(fit_draws(list(a = 1)), "`fit` must")
  expect_error(fit_draws(matrix(1:4, 2)), "named numeric draws")
  expect_output(print(d), "<fit_draws>")
})

test_that("fit_draws handles a univariate fit", {
  m <- cbind(h2 = stats::rnorm(10), sigma_G = abs(stats::rnorm(10)),
             sigma_E = abs(stats::rnorm(10)), beta = stats::rnorm(10))
  colnames(m) <- c("h2", "sigma_G", "sigma_E", "beta[1]")
  fit <- structure(list(draws = function(...) posterior::as_draws_matrix(m)),
                   class = "CmdStanMCMC", x_names = "(Intercept)")
  d <- fit_draws(fit, traits = "y")
  expect_equal(dim(d$h2), c(1L, 10L))
  expect_equal(dim(d$beta), c(1L, 1L, 10L))
  expect_equal(dimnames(d$beta), list("y", "(Intercept)", NULL))
  expect_equal(unname(draws_mean(d, "h2")), mean(m[, "h2"]))
})

test_that("fit_draws reads the Gibbs sampler and keeps its layout", {
  k <- 2; n <- 5
  G <- array(0, c(k, k, n), dimnames = list(c("a", "b"), c("a", "b"), NULL))
  R <- G
  for (i in seq_len(n)) {
    G[, , i] <- matrix(c(1, 0.2, 0.2, 1), k)
    R[, , i] <- diag(2)
  }
  beta <- array(stats::rnorm(k * 2 * n), c(k, 2, n),
                dimnames = list(c("a", "b"), c("Int", "x"), NULL))
  gibbs <- list(G = G, R = R, beta = beta, h2 = matrix(0.5, k, n),
                iterations = seq_len(n))
  d <- fit_draws(gibbs)
  expect_equal(names(d)[1:2], c("G", "E"))
  expect_equal(d$E, R)
  expect_equal(dimnames(d$beta)[[2]], c("Int", "x"))
  expect_equal(draws_mean(d, "P")[1, 2], mean(G[1, 2, ] + R[1, 2, ]))
})

test_that("fit_draws accepts back_transform output", {
  Y <- cbind(a = stats::rnorm(20), b = stats::rnorm(20), c = stats::rnorm(20))
  ss <- size_shape(Y)
  k <- 3
  m <- cbind(matrix(rep(as.vector(diag(k)), 4), 4, byrow = TRUE),
             matrix(rep(as.vector(diag(k)), 4), 4, byrow = TRUE))
  colnames(m) <- c(sprintf("G[%d,%d]", rep(1:k, k), rep(1:k, each = k)),
                   sprintf("E[%d,%d]", rep(1:k, k), rep(1:k, each = k)))
  post <- back_transform(m, ss)
  d <- fit_draws(post)
  expect_equal(dim(d$G), c(3L, 3L, 4L))
  expect_equal(attr(d, "traits"), ss$raw_traits)
  expect_equal(dimnames(d$G)[[1]], ss$raw_traits)
})
