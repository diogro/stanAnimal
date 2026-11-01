test_that("the Stan backend falls back to rstan with a usable message", {
  local_mocked_bindings(
    .is_installed = function(package) identical(package, "rstan"),
    .package = "stanimal")
  expect_warning(backend <- stanimal:::.resolve_stan_backend(TRUE),
                 'install.packages\\("cmdstanr", repos = c\\("https://stan-dev[.]r-universe[.]dev"')
  expect_false(backend)
})

test_that("the Stan backend explains what to install when neither engine is here", {
  local_mocked_bindings(
    .is_installed = function(...) FALSE,
    .package = "stanimal")
  expect_error(stanimal:::.resolve_stan_backend(TRUE), "cmdstanr or rstan")
  expect_error(stanimal:::.resolve_stan_backend(FALSE),
               "cmd_stan = FALSE needs rstan")
})

test_that("the Stan backend picks an installed engine", {
  skip_if_not_installed("cmdstanr")
  skip_if_not_installed("rstan")
  expect_true(stanimal:::.resolve_stan_backend(TRUE))
  expect_false(stanimal:::.resolve_stan_backend(FALSE))
})
