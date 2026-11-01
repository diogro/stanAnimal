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

test_that("random_covariance supports both correlation models", {
  for (method in c("lkj", "partial")) {
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
})
