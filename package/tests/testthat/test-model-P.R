skip_if_no_cmdstan()

fit_P <- function(sim, ...) {
  lmm_multiCov(sim$Y, sim$X, iter = 1000, chains = 2, cores = 2, seed = 1,
               refresh = 0, show_messages = FALSE, ...)
}

test_that("recovers P, corrP and original-scale beta", {
  sim <- simulate_animal(n_fam = 1, size = 300, E = G3, seed = 41,
                         genetic = FALSE)
  fit <- suppressWarnings(fit_P(sim))
  expect_s3_class(fit, "CmdStanMCMC")
  expect_clean_fit(fit)

  Pm <- post_matrix(fit, "P", 3)
  expect_lt(max(abs(Pm - G3)), 0.25)
  # close to the OLS residual covariance
  ols <- crossprod(stats::resid(stats::lm(sim$Y ~ sim$X - 1))) / (300 - 3)
  expect_lt(max(abs(Pm - ols)), 0.1)
  expect_lt(max(abs(post_matrix(fit, "corrP", 3) - stats::cov2cor(G3))), 0.12)

  b <- fit$summary("beta", "mean")
  est <- matrix(b$mean[match(sprintf("beta[%d,%d]", rep(1:3, 3), rep(1:3, each = 3)),
                             b$variable)], 3)
  expect_lt(max(abs(est[, 2:3] - t(sim$B)[, 2:3])), 0.15)
  expect_lt(max(abs(est[, 1] - sim$B[1, ])), 2)
})

test_that("P is consistent with sigma_P and corrP", {
  sim <- simulate_animal(n_fam = 1, size = 150, E = G3, seed = 42,
                         genetic = FALSE)
  fit <- suppressWarnings(fit_P(sim))
  d <- fit$draws(c("P", "corrP", "sigma_P"), format = "draws_matrix")
  expect_equal(as.numeric(d[, "P[1,2]"]),
               as.numeric(d[, "corrP[1,2]"] * d[, "sigma_P[1]"] * d[, "sigma_P[2]"]),
               tolerance = 1e-6)
  expect_equal(as.numeric(d[, "P[3,3]"]), as.numeric(d[, "sigma_P[3]"]^2),
               tolerance = 1e-6)
})

test_that("a single trait works", {
  sim <- simulate_animal(n_fam = 1, size = 200, E = matrix(1.5), seed = 43,
                         genetic = FALSE)
  fit <- suppressWarnings(fit_P(sim))
  expect_clean_fit(fit)
  expect_lt(abs(fit$summary("P", "mean")$mean - 1.5), 0.4)
  expect_true(all(c("beta[Int]", "beta[x]", "beta[z]") %in%
                    posterior::variables(animal_draws(fit, "beta"))))
})

test_that("animal_draws names multi-trait beta", {
  sim <- simulate_animal(n_fam = 1, size = 100, E = G3, seed = 44,
                         genetic = FALSE)
  fit <- suppressWarnings(fit_P(sim))
  expect_true("beta[2,x]" %in% posterior::variables(animal_draws(fit, "beta")))
})

test_that("a custom model file is used", {
  sim <- simulate_animal(n_fam = 1, size = 60, E = G3[1:2, 1:2], seed = 45,
                         genetic = FALSE)
  f <- tempfile(fileext = ".stan")
  writeLines("parameters { real x; } model { x ~ normal(0, 1); }", f)
  fit <- lmm_multiCov(sim$Y, sim$X, model_file = f, iter = 400, chains = 1,
                      refresh = 0, show_messages = FALSE,
                      control = list(adapt_delta = 0.8))
  expect_equal(setdiff(fit$metadata()$stan_variables, "lp__"), "x")
})
