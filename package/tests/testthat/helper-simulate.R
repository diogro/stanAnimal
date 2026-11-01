skip_if_no_cmdstan <- function() {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("cmdstanr")
  ok <- tryCatch({cmdstanr::cmdstan_version(); TRUE}, error = function(e) FALSE)
  if (!ok) testthat::skip("CmdStan is not installed")
}

skip_if_no_rstan <- function() {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("rstan")
  testthat::skip_if(Sys.getenv("STANIMAL_SKIP_RSTAN") == "true",
                    "rstan tests disabled")
}

# Relationship matrix with `n_fam` families of `size` full sibs (r = 0.5).
make_A <- function(n_fam, size) {
  block <- matrix(0.5, size, size)
  diag(block) <- 1
  kronecker(diag(n_fam), block)
}

rmvn <- function(n, Sigma) {
  matrix(stats::rnorm(n * nrow(Sigma)), n) %*% chol(Sigma)
}

# Simulates K traits with an intercept, a continuous predictor (mean 50, sd 10)
# and a 0/1 predictor. Returns the truth on the original scale.
simulate_animal <- function(n_fam = 25, size = 8, G = NULL, E = NULL,
                            seed = 1, genetic = TRUE) {
  set.seed(seed)
  n <- n_fam * size
  K <- if (!is.null(G)) nrow(G) else nrow(E)
  A <- make_A(n_fam, size)
  X <- cbind(Int = 1, x = stats::rnorm(n, 50, 10), z = stats::rbinom(n, 1, 0.5))
  B <- rbind(c(10, -5, 2)[seq_len(K)], c(0.05, 0.1, -0.02)[seq_len(K)],
             c(1, -0.5, 0.3)[seq_len(K)])
  u <- if (genetic) t(chol(A)) %*% rmvn(n, G) else 0
  Y <- X %*% B + u + rmvn(n, E)
  list(Y = Y, X = X, A = A, B = B, G = G, E = E, P = if (genetic) G + E else E)
}

G3 <- matrix(c(1.0, 0.5, 0.2,
               0.5, 0.8, 0.1,
               0.2, 0.1, 0.6), 3)
E3 <- matrix(c(1.0, 0.3, 0.0,
               0.3, 0.9, 0.2,
               0.0, 0.2, 0.7), 3)

# Posterior mean of a matrix-valued variable such as "G" from a CmdStanMCMC.
post_matrix <- function(fit, name, K) {
  s <- fit$summary(name, "mean")
  matrix(s$mean[match(sprintf("%s[%d,%d]", name, rep(seq_len(K), K),
                              rep(seq_len(K), each = K)), s$variable)], K)
}

expect_clean_fit <- function(fit, rhat_max = 1.05) {
  ds <- fit$diagnostic_summary(quiet = TRUE)
  expect_equal(sum(ds$num_divergent), 0)
  s <- fit$summary()
  expect_lt(max(s$rhat, na.rm = TRUE), rhat_max)
}
