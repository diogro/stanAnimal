skip_if_no_cmdstan()

sim <- simulate_animal(n_fam = 30, size = 8, G = matrix(0.8), E = matrix(0.8),
                       seed = 31)

fit_uni <- function(...) {
  suppressWarnings(lmm_animal(sim$Y, sim$X, sim$A, iter = 1000, chains = 2,
                              cores = 2, seed = 1, refresh = 0,
                              show_messages = FALSE, ...))
}

check_uni <- function(fit) {
  expect_s3_class(fit, "CmdStanMCMC")
  expect_clean_fit(fit)
  b <- fit$summary("beta", "mean")$mean
  expect_length(b, 3)
  # slope and indicator effects on the original scale, intercept loosely
  expect_lt(max(abs(b[2:3] - sim$B[2:3, 1])), 0.5)
  expect_lt(abs(b[1] - sim$B[1, 1]), 3)
  h2 <- fit$draws("h2", format = "draws_matrix")
  expect_true(all(h2 > 0 & h2 < 1))
  expect_lt(abs(mean(h2) - 0.5), 0.25)
  expect_true(all(c("beta[Int]", "beta[x]", "beta[z]") %in%
                    posterior::variables(animal_draws(fit, "beta"))))
}

test_that("variance parametrization recovers beta and h2", {
  check_uni(fit_uni())
})

test_that("h2 parametrization recovers beta and h2", {
  check_uni(fit_uni(parametrization = "h2"))
})

test_that("h2 parametrization accepts an informative beta prior", {
  flat <- fit_uni(parametrization = "h2")
  strong <- fit_uni(parametrization = "h2", beta_prior = c(5, 50))
  h_flat <- mean(flat$draws("h2", format = "draws_matrix"))
  h_strong <- mean(strong$draws("h2", format = "draws_matrix"))
  expect_lt(h_strong, h_flat)
  expect_lt(h_strong, 0.2)
})

test_that("the E-BFMI warning is raised only for the variance parametrization", {
  skip_on_cran()
  mock <- structure(list(diagnostic_summary = function(...) list(ebfmi = 0.1)),
                    class = "CmdStanMCMC")
  expect_warning(stanAnimal:::.warn_low_bfmi(mock), "parametrization")
})
