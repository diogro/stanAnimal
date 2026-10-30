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
