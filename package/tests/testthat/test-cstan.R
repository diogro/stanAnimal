test_that("cstan fits a Stan program and returns a CmdStanMCMC", {
  skip_if_no_cmdstan()
  code <- "
  data { int<lower=0> N; vector[N] y; }
  parameters { real mu; real<lower=0> sigma; }
  model { y ~ normal(mu, sigma); mu ~ normal(0, 10); sigma ~ normal(0, 5); }"
  set.seed(1)
  y <- stats::rnorm(50, 3, 2)
  fit <- stanAnimal:::cstan(code, list(N = 50, y = y), chains = 2, cores = 2,
                            iter = 800, refresh = 0, show_messages = FALSE,
                            seed = 1)
  expect_s3_class(fit, "CmdStanMCMC")
  expect_equal(posterior::ndraws(fit$draws()), 800)
  expect_lt(abs(fit$summary("mu", "mean")$mean - 3), 0.8)
})

test_that("cstan passes adapt_delta and max_treedepth from control", {
  skip_if_no_cmdstan()
  code <- "parameters { real x; } model { x ~ normal(0, 1); }"
  fit <- stanAnimal:::cstan(code, list(), chains = 1, iter = 200,
                            control = list(adapt_delta = 0.8, max_treedepth = 5),
                            refresh = 0, show_messages = FALSE)
  expect_equal(fit$metadata()$adapt_delta, 0.8)
  expect_equal(fit$metadata()$max_treedepth, 5)
})
