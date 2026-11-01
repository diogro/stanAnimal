skip_if_no_rstan()

test_that("rstan path fits the P model with renamed beta", {
  sim <- simulate_animal(n_fam = 1, size = 150, E = G3, seed = 53,
                         genetic = FALSE)
  fit <- suppressWarnings(lm_phen_cov(sim$Y, sim$X, cmd_stan = FALSE,
                                       iter = 600, chains = 2, cores = 1,
                                       seed = 1, refresh = 0))
  expect_s4_class(fit, "stanfit")
  expect_true(all(c("beta[1,Int]", "beta[2,x]", "beta[3,z]") %in% names(fit)))
  s <- rstan::summary(fit, "P")$summary
  expect_lt(max(abs(matrix(s[, "mean"], 3) - G3)), 0.5)
})
