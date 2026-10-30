test_that("cache directory follows the option and is created", {
  d <- file.path(tempfile("cache"), "nested")
  withr_old <- options(stanAnimal.cache_dir = d)
  on.exit(options(withr_old))
  expect_equal(stan_cache_dir(), normalizePath(d))
  expect_true(dir.exists(d))
  writeLines("x", file.path(d, "a.txt"))
  expect_equal(stanAnimal_clear_cache(), 1)
  expect_length(list.files(d), 0)
})

test_that("compiled models are reused across calls", {
  skip_if_no_cmdstan()
  d <- tempfile("cache")
  old <- options(stanAnimal.cache_dir = d)
  on.exit(options(old))
  code <- "parameters { real x; } model { x ~ normal(0, 1); }"
  run <- function() system.time(suppressMessages(
    stanAnimal:::cstan(code, list(), chains = 1, iter = 100, refresh = 0,
                       show_messages = FALSE)))[["elapsed"]]
  first <- run()
  exes <- list.files(d, full.names = TRUE)
  exes <- exes[!grepl("\\.stan$|\\.hpp$|\\.o$", exes) & !dir.exists(exes)]
  expect_length(exes, 1)
  second <- run()
  expect_lt(second, first)
  expect_length(list.files(d, pattern = "\\.stan$"), 1)
})
