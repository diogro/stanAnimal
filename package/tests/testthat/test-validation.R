sim <- simulate_animal(n_fam = 5, size = 4, G = G3, E = E3)
Y <- sim$Y
X <- sim$X
A <- sim$A

# The checks run before any Stan code, so these are fast.
test_that("lmm_animal validates Y and X", {
  Yna <- Y; Yna[1, 1] <- NA
  expect_error(lmm_animal(Yna, X, A), "Y must contain only finite")
  Xinf <- X; Xinf[2, 2] <- Inf
  expect_error(lmm_animal(Y, Xinf, A), "X must contain only finite")
  expect_error(lmm_animal(Y, X[-1, ], A), "same number of rows")
  expect_error(lmm_animal(Y, cbind(X, X[, 2]), A), "full column rank")
})

test_that("lmm_animal validates A and flags", {
  expect_error(lmm_animal(Y, X, A[-1, -1]), "square matrix")
  As <- A; As[1, 2] <- 5
  expect_error(lmm_animal(Y, X, As), "symmetric")
  Ai <- diag(nrow(A)); Ai[1, 1] <- -1
  expect_error(lmm_animal(Y, X, Ai), "positive semidefinite")
  Asing <- A; Asing[2, ] <- Asing[1, ]; Asing[, 2] <- Asing[, 1]
  expect_error(lmm_animal(Y, X, Asing), "positive definite")
  expect_error(lmm_animal(Y, X, A, blend_A = NA), "blend_A")
  expect_error(lmm_animal(Y, X, A, residual_structure = "bad"), "should be one of")
  expect_error(lmm_animal(Y, X, A, parametrization = "bad"), "should be one of")
})

test_that("constant traits are rejected before the model is fitted", {
  # univariate: a constant response has no standard deviation to divide by
  expect_error(lmm_animal(rep(5, nrow(Y)), X, A),
               "every trait in Y must vary")
  # multivariate: one constant column is enough
  expect_error(lmm_animal(cbind(Y, flat = 1), X, A),
               "every trait in Y must vary.*: flat")
  # partially named and unnamed responses still say which column it was
  expect_error(lmm_animal(cbind(Y, 1), X, A), "column 4")
  expect_error(lmm_animal(cbind(unname(Y[, 1]), 1), X, A), "column 2")
  # the phenotypic model shares the check, and so does the rstan path
  expect_error(lmm_multiCov(cbind(Y, flat = 1), X), "flat")
  expect_error(lmm_animal(cbind(Y, flat = 1), X, A, cmd_stan = FALSE), "flat")
})

test_that("X must be able to represent a constant", {
  no_intercept <- cbind(x = X[, "x"])
  expect_error(lmm_animal(Y, no_intercept, A), "cannot represent a constant")
  expect_error(lmm_multiCov(Y, no_intercept), "cannot represent a constant")
  expect_error(lmm_animal(Y, no_intercept, A, cmd_stan = FALSE),
               "cannot represent a constant")
  # factor indicators span the constant, so they are fine
  expect_silent(stanAnimal:::.check_response_design(
    Y, stats::model.matrix(~ 0 + factor(rep(c(0, 1), length.out = nrow(Y))))))
})

test_that("lmm_multiCov validates Y and X", {
  Yna <- Y; Yna[1, 1] <- NA
  expect_error(lmm_multiCov(Yna, X), "Y must contain only finite")
  expect_error(lmm_multiCov(Y, X[-1, ]), "same number of rows")
  expect_error(lmm_multiCov(Y, cbind(X, X[, 2])), "full column rank")
})

test_that("lmm_animal_gibbs validates its inputs", {
  expect_error(lmm_animal_gibbs(Y, X, A[-1, -1]), "square matrix")
  expect_error(lmm_animal_gibbs(Y, cbind(X, X[, 2]), A), "full column rank")
  expect_error(lmm_animal_gibbs(Y, X, A, iter = 10, warmup = 20), "iter")
  Ai <- diag(nrow(A)); Ai[1, 1] <- -1
  expect_error(lmm_animal_gibbs(Y, X, Ai), "positive semidefinite")
  expect_error(lmm_animal_gibbs(cbind(Y[, 1], 1), X, A), "nonzero variance")
})

test_that("custom model files are read and missing defaults error", {
  f <- tempfile(fileext = ".stan")
  writeLines("parameters { real x; } model { x ~ normal(0, 1); }", f)
  expect_match(stanAnimal:::.read_animal_model_code(f), "normal\\(0, 1\\)")
  expect_match(stanAnimal:::.read_animal_model_code(NULL, "animalModelG.stan"),
               "lkj_corr_cholesky")
  expect_match(stanAnimal:::.read_animal_model_code(NULL, "phenotypicModel.stan"),
               "sigma_P")
  expect_error(stanAnimal:::.read_animal_model_code(NULL, "nope.stan"),
               "Could not locate")
})

test_that("all shipped Stan files are discoverable", {
  for (f in c("animalModelG.stan", "animalModelUni.stan",
              "animalModelUniH2.stan", "phenotypicModel.stan"))
    expect_true(nzchar(stanAnimal:::.read_animal_model_code(NULL, f)))
})

test_that("an invalid LKJ shape is rejected before Stan runs", {
  expect_error(lmm_animal(Y, X, A, lkj_prior = 0), "`lkj_prior`")
  expect_error(lmm_animal(Y, X, A, lkj_prior = -1), "`lkj_prior`")
  expect_error(lmm_animal(Y, X, A, lkj_prior = c(1, 2)), "`lkj_prior`")
  expect_error(lmm_multiCov(Y, X, lkj_prior = 0), "`lkj_prior`")
  expect_error(lmm_multiCov(Y, X, lkj_prior = NA_real_), "`lkj_prior`")
})
