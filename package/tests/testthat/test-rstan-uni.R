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

test_that("the rstan path forwards step_size as stepsize", {
  set.seed(64)
  n <- 40
  A <- diag(n)
  X <- cbind(Intercept = 1, x = stats::rnorm(n))
  y <- stats::rnorm(n)
  fit <- suppressWarnings(lmm_animal(
    y, X, A, cmd_stan = FALSE, iter = 200, warmup = 100, chains = 1, refresh = 0,
    control = list(adapt_engaged = FALSE, step_size = 0.3)))
  stepsize <- rstan::get_sampler_params(fit, inc_warmup = FALSE)[[1]][, "stepsize__"]
  expect_equal(unique(stepsize), 0.3)
})
