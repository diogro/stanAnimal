test_that(".beta_names uses trait indices only for several traits", {
  expect_equal(stanimal:::.beta_names(c("Int", "x"), 1)$new,
               c("beta[Int]", "beta[x]"))
  expect_equal(stanimal:::.beta_names(c("Int", "x"), 1)$old,
               c("beta[1]", "beta[2]"))
  nm <- stanimal:::.beta_names(c("Int", "x"), 2)
  expect_equal(nm$old, c("beta[1,1]", "beta[2,1]", "beta[1,2]", "beta[2,2]"))
  expect_equal(nm$new, c("beta[1,Int]", "beta[2,Int]", "beta[1,x]", "beta[2,x]"))
})

test_that("rename_beta renames beta and leaves other variables alone", {
  skip_if_not_installed("posterior")
  draws <- posterior::as_draws_array(array(
    stats::rnorm(10 * 2 * 3), c(10, 2, 3),
    dimnames = list(NULL, NULL, c("beta[1]", "beta[2]", "h2"))))
  fit <- structure(list(draws = function(variables = NULL, ...) draws),
                   x_names = c("Int", "x"))
  out <- rename_beta(fit)
  expect_equal(posterior::variables(out), c("beta[Int]", "beta[x]", "h2"))
  attr(fit, "x_names") <- NULL
  expect_equal(posterior::variables(rename_beta(fit)),
               c("beta[1]", "beta[2]", "h2"))
})

test_that("rename_beta renames multi-trait beta", {
  skip_if_not_installed("posterior")
  vars <- c("beta[1,1]", "beta[2,1]", "beta[1,2]", "beta[2,2]")
  draws <- posterior::as_draws_array(array(
    stats::rnorm(5 * 1 * 4), c(5, 1, 4), dimnames = list(NULL, NULL, vars)))
  fit <- structure(list(draws = function(variables = NULL, ...) draws),
                   x_names = c("Int", "x"))
  expect_equal(posterior::variables(rename_beta(fit)),
               c("beta[1,Int]", "beta[2,Int]", "beta[1,x]", "beta[2,x]"))
})

test_that(".warn_low_bfmi warns only below the threshold", {
  mk <- function(e) structure(list(diagnostic_summary = function(...)
    list(ebfmi = e)), class = "CmdStanMCMC")
  expect_warning(stanimal:::.warn_low_bfmi(mk(c(0.1, 0.2, 0.9, 0.9))), "h2")
  expect_no_warning(stanimal:::.warn_low_bfmi(mk(c(0.8, 0.9))))
})
