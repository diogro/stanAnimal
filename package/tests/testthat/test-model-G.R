skip_if_no_cmdstan()

fit_G <- function(sim, ...) {
  suppressWarnings(lmm_animal(sim$Y, sim$X, sim$A, iter = 1000, chains = 2,
                              cores = 2, seed = 1, refresh = 0,
                              show_messages = FALSE, ...))
}

test_that("unstructured G/E model recovers G, E and original-scale beta", {
  sim <- simulate_animal(n_fam = 30, size = 8, G = G3, E = E3, seed = 21)
  fit <- fit_G(sim)
  expect_s3_class(fit, "CmdStanMCMC")
  expect_clean_fit(fit)

  expect_lt(max(abs(post_matrix(fit, "G", 3) - G3)), 0.6)
  expect_lt(max(abs(post_matrix(fit, "E", 3) - E3)), 0.6)
  expect_lt(max(abs(post_matrix(fit, "P", 3) - (G3 + E3))), 0.6)

  # beta[k, j] is trait k, column j of X, on the original scale
  b <- fit$summary("beta", "mean")
  est <- matrix(b$mean[match(sprintf("beta[%d,%d]", rep(1:3, 3), rep(1:3, each = 3)),
                             b$variable)], 3)
  expect_lt(max(abs(est[, 2:3] - t(sim$B)[, 2:3])), 0.25)
  expect_lt(max(abs(est[, 1] - sim$B[1, ])), 3)
})

test_that("outputs are internally consistent", {
  sim <- simulate_animal(n_fam = 20, size = 6, G = G3, E = E3, seed = 22)
  fit <- fit_G(sim)
  d <- fit$draws(c("G", "E", "P", "corrG", "sigma_G", "sigma_E", "h2"),
                 format = "draws_matrix")
  vars <- colnames(d)
  sel <- function(x) d[, grep(paste0("^", x, "\\["), vars), drop = FALSE]
  expect_equal(unname(as.matrix(sel("P"))),
               unname(as.matrix(sel("G")) + as.matrix(sel("E"))), tolerance = 1e-6)
  expect_true(all(sel("h2") > 0 & sel("h2") < 1))
  expect_true(all(abs(sel("corrG")) <= 1 + 1e-8))
  expect_true(all(sel("sigma_G") > 0))
  expect_true(all(sel("sigma_E") > 0))
  # sigma_G^2 are the diagonal elements of G
  g11 <- d[, "G[1,1]"]
  expect_equal(as.numeric(g11), as.numeric(d[, "sigma_G[1]"]^2), tolerance = 1e-6)
})

test_that("animal_draws names beta after the columns of X", {
  sim <- simulate_animal(n_fam = 20, size = 6, G = G3, E = E3, seed = 23)
  fit <- fit_G(sim)
  v <- posterior::variables(animal_draws(fit, "beta"))
  expect_true(all(c("beta[1,Int]", "beta[2,x]", "beta[3,z]") %in% v))
  expect_equal(length(v), 9)
})

test_that("diagonal residual structure runs and returns a diagonal E", {
  sim <- simulate_animal(n_fam = 20, size = 6, G = G3, E = diag(diag(E3)),
                         seed = 24)
  fit <- fit_G(sim, residual_structure = "diagonal")
  E <- post_matrix(fit, "E", 3)
  expect_equal(E[upper.tri(E)], rep(0, 3), tolerance = 1e-8)
  expect_lt(max(abs(E - diag(diag(E3)))), 0.6)
})

test_that("two traits and a constant non-unit column are backtransformed", {
  G2 <- G3[1:2, 1:2]; E2 <- E3[1:2, 1:2]
  sim <- simulate_animal(n_fam = 25, size = 8, G = G2, E = E2, seed = 25)
  sim$X[, "Int"] <- 5
  fit <- fit_G(sim)
  expect_clean_fit(fit)
  b <- fit$summary("beta", "mean")
  est <- matrix(b$mean[match(sprintf("beta[%d,%d]", rep(1:2, 3), rep(1:3, each = 2)),
                             b$variable)], 2)
  # the response intercept 10 is now carried by the constant column of 5
  expect_equal(est[1, 1], sim$B[1, 1] / 5, tolerance = 0.5)
  expect_lt(max(abs(est[, 2:3] - t(sim$B)[, 2:3])), 0.25)
})
