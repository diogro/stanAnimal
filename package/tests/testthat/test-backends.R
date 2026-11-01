test_that("the Stan backend falls back to rstan with a usable message", {
  local_mocked_bindings(
    .is_installed = function(package) identical(package, "rstan"),
    .package = "stanAnimal")
  expect_warning(backend <- stanAnimal:::.resolve_stan_backend(TRUE),
                 'remotes::install_github\\("stan-dev/cmdstanr"\\)')
  expect_false(backend)
})

test_that("the Stan backend explains what to install when neither engine is here", {
  local_mocked_bindings(
    .is_installed = function(...) FALSE,
    .package = "stanAnimal")
  expect_error(stanAnimal:::.resolve_stan_backend(TRUE), "cmdstanr or rstan")
  expect_error(stanAnimal:::.resolve_stan_backend(FALSE),
               "cmd_stan = FALSE needs rstan")
})

test_that("the Stan backend picks an installed engine", {
  skip_if_not_installed("cmdstanr")
  skip_if_not_installed("rstan")
  expect_true(stanAnimal:::.resolve_stan_backend(TRUE))
  expect_false(stanAnimal:::.resolve_stan_backend(FALSE))
})
