skip_if_no_rstan()

sim <- simulate_animal(n_fam = 20, size = 8, G = matrix(0.8), E = matrix(0.8),
                       seed = 52)

test_that("rstan path fits both univariate parametrizations", {
  for (p in c("variance", "h2")) {
    fit <- suppressWarnings(lmm_animal(sim$Y, sim$X, sim$A, cmd_stan = FALSE,
                                       parametrization = p, iter = 600,
                                       chains = 2, cores = 1, seed = 1,
                                       refresh = 0))
    expect_s4_class(fit, "stanfit")
    expect_true(all(c("beta[Int]", "beta[x]", "beta[z]") %in% names(fit)))
    h2 <- rstan::extract(fit, "h2")[[1]]
    expect_true(all(h2 > 0 & h2 < 1))
  }
})
