make_draws <- function(k, n = 200) {
  vars <- c(sprintf("h2[%d]", seq_len(k)),
            if (k > 1) c("corrG[1,2]", "corrG[2,1]"))
  posterior::as_draws_array(array(stats::runif(n * 2 * length(vars)),
                                  c(n, 2, length(vars)),
                                  dimnames = list(NULL, NULL, vars)))
}

test_that("prior draws have valid ranges and names", {
  d <- animal_prior_draws(3, n = 300, trait_names = c("a", "b", "c"))
  expect_setequal(unique(d$parameter),
                  c("a", "b", "c", "a ~ b", "a ~ c", "b ~ c"))
  expect_true(all(d$value > -1 & d$value < 1))
  expect_true(all(d$value[d$type == "Heritability"] > 0))
  lkj <- stanAnimal:::.rlkj(4, 2)
  expect_true(all(eigen(lkj)$values > 0))
  expect_equal(diag(lkj), rep(1, 4))
})

test_that("LKJ(1) pairwise correlations are uniform-ish and eta shrinks", {
  set.seed(1)
  r2 <- vapply(1:2000, function(i) stanAnimal:::.rlkj(2, 2)[1, 2], 0)
  expect_lt(abs(var(r2) - 1 / 5), 0.03) # Var = 1 / (2 * eta + 1)
})

test_that("parameter draws are extracted from draws objects", {
  expect_error(animal_parameter_draws(make_draws(3)[, , 1:3]), "corrG")
  d <- animal_parameter_draws(make_draws(1))
  expect_equal(unique(d$parameter), "h2")
})

test_that("plots build", {
  skip_if_not_installed("ggplot2")
  vars <- c("h2[1]", "h2[2]", "corrG[1,2]", "corrG[2,1]")
  fit <- posterior::as_draws_array(array(
    stats::runif(200 * 2 * 4), c(200, 2, 4), dimnames = list(NULL, NULL, vars)))
  expect_s3_class(plot_prior_posterior(fit, c("a", "b"), n_prior = 200), "ggplot")
  p <- plot_forest(fit, c("a", "b"), n_prior = 200)
  expect_s3_class(p, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p))
})

test_that("plot_G_ellipses builds", {
  set.seed(1)
  Y <- matrix(stats::rnorm(60 * 3), 60, 3, dimnames = list(NULL, c("a", "b", "c")))
  vars <- c(outer(1:3, 1:3, function(i, j) sprintf("G[%d,%d]", i, j)))
  G <- replicate(50, c(crossprod(matrix(stats::rnorm(30, sd = 0.3), 10, 3)) + diag(0.3, 3)))
  d <- posterior::as_draws_array(array(t(G), c(50, 1, 9), dimnames = list(NULL, NULL, vars)))
  p <- plot_G_ellipses(d, Y, n_show = 10)
  expect_s3_class(p, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p))
  expect_error(plot_G_ellipses(d, Y[, 1, drop = FALSE]), "two traits")
})

test_that("plot_corr builds and keeps matrix labels", {
  M <- stats::cor(mtcars[, 1:4])
  p <- plot_corr(M, "t")
  expect_s3_class(p, "ggplot")
  expect_equal(levels(p$data$Var1), colnames(M))
  expect_no_error(ggplot2::ggplot_build(p))
  expect_s3_class(plot_corr(unname(M)), "ggplot")
})

make_G_beta_draws <- function(k, J, B, n = 40) {
  vars <- c(sprintf("G[%d,%d]", rep(seq_len(k), k), rep(seq_len(k), each = k)),
            sprintf("beta[%d,%d]", rep(seq_len(k), J), rep(seq_len(J), each = k)))
  arr <- array(0, c(n, 2, length(vars)), dimnames = list(NULL, NULL, vars))
  values <- rep(c(as.vector(diag(k) + 0.1), as.vector(B)), each = n)
  for (chain in 1:2) arr[, chain, ] <- values
  posterior::as_draws_array(arr)
}

test_that("plot_G_ellipses plots residuals when X is given", {
  set.seed(4)
  n <- 40; k <- 2; J <- 2
  Y <- cbind(a = stats::rnorm(n, 10, 2), b = stats::rnorm(n, 5, 1))
  X <- cbind(1, stats::rnorm(n))
  B <- matrix(c(10, 5, 0.5, -0.5), k, J)
  d <- make_G_beta_draws(k, J, B)

  p <- plot_G_ellipses(d, Y, X = X)
  built <- ggplot2::ggplot_build(p)$data
  points <- built[[which(vapply(built, function(layer)
    all(layer$geom == "point"), logical(1)))]]
  expected <- Y - X %*% t(B)
  expect_equal(sort(points$x), sort(expected[, 1]))
  expect_equal(sort(points$y), sort(expected[, 2]))
  expect_true("Residual P" %in% unique(p$layers[[4]]$data$group))

  # without X the raw traits are plotted and the covariance is the observed one
  p_raw <- plot_G_ellipses(d, Y)
  expect_true("Observed P" %in% unique(p_raw$layers[[4]]$data$group))
})

test_that("plot_G_ellipses validates X", {
  Y <- cbind(a = stats::rnorm(30), b = stats::rnorm(30))
  d <- make_G_beta_draws(2, 2, matrix(c(1, 1, 0.2, -0.2), 2, 2), n = 10)
  expect_error(plot_G_ellipses(d, Y, X = cbind(1, 2, 3)), "one row per row")
  expect_error(plot_G_ellipses(d, Y, X = matrix("a", 30, 2)), "finite numeric")
  expect_error(plot_G_ellipses(d, Y, X = cbind(1, stats::rnorm(30),
                                                stats::rnorm(30))), "fixed effects")
})
