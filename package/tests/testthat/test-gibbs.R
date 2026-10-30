test_that("univariate sampler recovers fixed effects and h2", {
  sim <- simulate_animal(n_fam = 40, size = 8, G = matrix(1), E = matrix(1),
                         seed = 11)
  fit <- lmm_animal_gibbs(sim$Y, sim$X, sim$A, iter = 2000, warmup = 500)
  expect_named(fit, c("G", "R", "beta", "h2", "iterations"))
  expect_equal(dim(fit$beta), c(1, 3, 1500))
  est <- unname(apply(fit$beta, 2, mean))
  expect_lt(max(abs(est - sim$B[, 1])[2:3]), 0.2)
  expect_equal(mean(fit$h2), 0.5, tolerance = 0.25)
})

test_that("multivariate sampler returns G, R, beta and h2 of the right size", {
  sim <- simulate_animal(n_fam = 30, size = 8, G = G3, E = E3, seed = 12)
  fit <- lmm_animal_gibbs(sim$Y, sim$X, sim$A, iter = 1500, warmup = 500)
  expect_equal(dim(fit$G), c(3, 3, 1000))
  expect_equal(dim(fit$R), c(3, 3, 1000))
  expect_equal(dim(fit$beta), c(3, 3, 1000))
  expect_equal(dim(fit$h2), c(3, 1000))
  Gm <- apply(fit$G, c(1, 2), mean)
  expect_lt(max(abs(Gm - G3)), 1)
  Bm <- apply(fit$beta, c(1, 2), mean)
  expect_lt(max(abs(t(Bm) - sim$B)[2:3, ]), 0.2)
})

test_that("thinning reduces the number of stored draws", {
  sim <- simulate_animal(n_fam = 10, size = 5, G = matrix(1), E = matrix(1))
  fit <- lmm_animal_gibbs(sim$Y, sim$X, sim$A, iter = 400, warmup = 100,
                          thin = 3)
  expect_equal(dim(fit$beta)[3], 100)
})
