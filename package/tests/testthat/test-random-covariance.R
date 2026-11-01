test_that("random_covariance returns a named positive-definite covariance", {
  set.seed(1)
  S <- random_covariance(4, variances = c(1, 2, 3, 4))
  expect_equal(dim(S), c(4L, 4L))
  expect_equal(diag(S), c(1, 2, 3, 4), ignore_attr = TRUE)
  expect_equal(S, t(S), ignore_attr = TRUE)
  expect_equal(dimnames(S), list(paste0("trait", 1:4), paste0("trait", 1:4)))
  expect_gt(min(eigen(S, symmetric = TRUE, only.values = TRUE)$values), 0)
  expect_equal(unname(diag(cov2cor(S))), rep(1, 4))
})

test_that("random_covariance supports every correlation model", {
  for (method in c("lkj", "partial", "factor")) {
    set.seed(2)
    S <- random_covariance(5, min_var = 0.5, max_var = 1.5, method = method)
    expect_true(all(diag(S) >= 0.5 & diag(S) <= 1.5))
    expect_gt(min(eigen(S, symmetric = TRUE, only.values = TRUE)$values), 0)
  }
})

test_that("random_covariance draws several matrices", {
  set.seed(3)
  S <- random_covariance(3, n_matrices = 4, min_var = 1, max_var = 2)
  expect_length(S, 4L)
  expect_true(all(vapply(S, function(x) all(dim(x) == 3L), logical(1))))
  expect_true(all(vapply(S, function(x) all(diag(x) >= 1 & diag(x) <= 2),
                         logical(1))))
  expect_false(identical(S[[1L]], S[[2L]]))

  fixed <- random_covariance(3, n_matrices = 2, variances = c(1, 2, 3))
  expect_equal(diag(fixed[[1L]]), c(1, 2, 3), ignore_attr = TRUE)
  expect_equal(diag(fixed[[2L]]), c(1, 2, 3), ignore_attr = TRUE)
})

test_that("random_covariance handles one and two traits with partial correlations", {
  set.seed(4)
  S1 <- random_covariance(1, method = "partial")
  expect_equal(dim(S1), c(1L, 1L))
  expect_equal(unname(S1[1, 1]), 1)

  S2 <- random_covariance(2, method = "partial")
  expect_equal(dim(S2), c(2L, 2L))
  expect_equal(unname(diag(S2)), c(1, 1))
  expect_equal(S2[1, 2], S2[2, 1])
  expect_true(abs(S2[1, 2]) <= 1)
})

test_that("random_covariance uses the LKJ shape and trait names", {
  set.seed(5)
  S <- random_covariance(3, lkj_prior = 1, trait_names = c("a", "b", "c"))
  expect_equal(dimnames(S), list(c("a", "b", "c"), c("a", "b", "c")))
  expect_gt(min(eigen(S, symmetric = TRUE, only.values = TRUE)$values), 0)
})

test_that("factor method targets integration with a positive size axis", {
  relvar <- stanimal:::.relative_eigen_variance
  for (target in c(0, 0.3, 0.8)) {
    set.seed(6)
    S <- random_covariance(6, method = "factor", integration = target)
    R <- cov2cor(S)
    expect_equal(unname(diag(R)), rep(1, 6))
    expect_equal(R, t(R), ignore_attr = TRUE)
    expect_gt(min(eigen(R, symmetric = TRUE, only.values = TRUE)$values), 0)
    expect_equal(relvar(R), target, tolerance = 1e-6)
    if (target > 0) {
      axis <- eigen(R, symmetric = TRUE)$vectors[, 1]
      expect_true(all(axis > 0) || all(axis < 0))
    }
  }
})

test_that("factor method spans integration but keeps the size axis", {
  set.seed(7)
  S <- random_covariance(5, n_matrices = 20, method = "factor")
  r2 <- vapply(S, function(x) stanimal:::.relative_eigen_variance(cov2cor(x)),
               numeric(1))
  expect_true(all(r2 >= 0.1 & r2 <= 0.8))
  expect_gt(stats::var(r2), 0)
  axes <- vapply(S, function(x) {
    axis <- eigen(cov2cor(x), symmetric = TRUE)$vectors[, 1]
    if (sum(axis) < 0) axis <- -axis
    axis
  }, numeric(5))
  expect_equal(axes, matrix(axes[, 1], 5, 20), ignore_attr = TRUE,
               tolerance = 1e-8)
})

test_that("factor loadings fix the size pattern", {
  loadings <- c(0.5, 1, 1.5, 2)
  S <- random_covariance(4, method = "factor", loadings = loadings,
                         integration = 0.5, variances = rep(1, 4))
  R <- cov2cor(S)
  ratios <- R[upper.tri(R)] / outer(loadings, loadings)[upper.tri(R)]
  expect_equal(ratios, rep(ratios[1], length(ratios)), tolerance = 1e-8)
  axis <- eigen(R, symmetric = TRUE)$vectors[, 1]
  expect_true(all(axis > 0) || all(axis < 0))
})

test_that("random_covariance validates its inputs", {
  expect_error(random_covariance(0), "`n_traits`")
  expect_error(random_covariance(3, n_matrices = 0), "`n_matrices`")
  expect_error(random_covariance(3, variances = c(1, 1)), "`variances`")
  expect_error(random_covariance(3, variances = c(1, -1, 1)), "`variances`")
  expect_error(random_covariance(3, variances = c(1, Inf, 1)), "`variances`")
  expect_error(random_covariance(3, min_var = 2, max_var = 1), "`min_var`")
  expect_error(random_covariance(3, min_var = 0), "`min_var`")
  expect_error(random_covariance(3, lkj_prior = 0), "`lkj_prior`")
  expect_error(random_covariance(3, trait_names = "a"), "`trait_names`")
  expect_error(random_covariance(3, method = "none"), "should be one of")
  expect_error(random_covariance(3, method = "factor", integration = 1),
               "`integration`")
  expect_error(random_covariance(3, method = "factor", integration = "x"),
               "`integration`")
  expect_error(random_covariance(3, method = "factor",
                                 integration = c(0.5, 0.2)), "increasing")
  expect_error(random_covariance(3, method = "factor", loadings = c(1, 0, 1)),
               "`loadings`")
  expect_error(random_covariance(3, method = "factor", loadings = 1),
               "`loadings`")
  expect_error(random_covariance(3, loadings = rep(1, 3)),
               "only used by method")
  expect_error(random_covariance(3, integration = 0.5),
               "only used by method")
})
