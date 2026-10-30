skip_if_no_rstan()

sim <- simulate_animal(n_fam = 20, size = 8, G = G3, E = E3, seed = 51)

test_that("rstan path fits the multivariate G model with renamed beta", {
  fit <- suppressWarnings(lmm_animal(sim$Y, sim$X, sim$A, cmd_stan = FALSE,
                                     iter = 600, chains = 2, cores = 1,
                                     seed = 1, refresh = 0))
  expect_s4_class(fit, "stanfit")
  expect_true(all(c("beta[1,Int]", "beta[2,x]", "beta[3,z]") %in% names(fit)))
  s <- rstan::summary(fit, "G")$summary
  expect_equal(nrow(s), 9)
  expect_lt(max(s[, "Rhat"]), 1.1)
  Gm <- matrix(s[, "mean"], 3)
  expect_lt(max(abs(Gm - G3)), 0.8)
})
