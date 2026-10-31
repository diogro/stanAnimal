sim <- simulate_animal(n_fam = 5, size = 4, G = G3, E = E3)
Y <- sim$Y
X <- sim$X
A <- sim$A
ids <- paste0("id", seq_len(nrow(Y)))

test_that("A must match the dimensions of Y", {
  expect_error(lmm_animal(Y, X, A[-1, -1]), "square matrix")
  expect_error(lmm_animal_gibbs(Y, X, A[-1, -1]), "square matrix")
  expect_error(stanAnimal:::.check_relationship_matrix(as.vector(A), Y),
               "square matrix")
})

test_that("A is put in the row order of Y by name", {
  An <- A
  rownames(An) <- colnames(An) <- ids
  Yn <- Y
  rownames(Yn) <- ids
  expect_equal(stanAnimal:::.check_relationship_matrix(An, Yn), An)
  permuted <- rev(seq_len(nrow(A)))
  expect_equal(stanAnimal:::.check_relationship_matrix(An[permuted, permuted], Yn),
               An)
})

test_that("names that disagree between A and Y are an error", {
  An <- A
  rownames(An) <- colnames(An) <- ids
  Ybad <- Y
  rownames(Ybad) <- c("someone_else", ids[-1])
  expect_error(stanAnimal:::.check_relationship_matrix(An, Ybad), "do not match")
  expect_error(lmm_animal(Ybad, X, An), "do not match")
  expect_error(lmm_animal_gibbs(Ybad, X, An), "do not match")
})

test_that("A must have matching row and column names", {
  An <- A
  rownames(An) <- ids
  colnames(An) <- rev(ids)
  Yn <- Y
  rownames(Yn) <- ids
  expect_error(stanAnimal:::.check_relationship_matrix(An, Yn),
               "row and column names")
  Ad <- An
  rownames(Ad) <- c(ids[1], ids[1], ids[-c(1, 2)])
  colnames(Ad) <- rownames(Ad)
  expect_error(stanAnimal:::.check_relationship_matrix(Ad, Yn), "duplicated")
})

test_that("unnamed A and Y fall back to matching by position", {
  Yn <- Y
  rownames(Yn) <- ids
  expect_silent(stanAnimal:::.check_relationship_matrix(unname(A), Yn))
  expect_silent(stanAnimal:::.check_relationship_matrix(A, Y))
})
