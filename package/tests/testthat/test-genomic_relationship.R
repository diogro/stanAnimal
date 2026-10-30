test_that("matches an independent VanRaden calculation", {
  set.seed(1)
  M <- matrix(rbinom(40 * 300, 2, 0.35), 40)
  G <- genomic_relationship(M, maf = 0)
  p <- colMeans(M) / 2
  Z <- sweep(M, 2, 2 * p)
  expect_equal(unname(G), unname(tcrossprod(Z) / (2 * sum(p * (1 - p)))),
               tolerance = 1e-10, ignore_attr = TRUE)
  expect_equal(dim(G), c(40, 40))
  expect_true(isSymmetric(unname(G)))
  expect_equal(attr(G, "n_markers"), 300)
})

test_that("keeps row names", {
  set.seed(1)
  M <- matrix(rbinom(10 * 50, 2, 0.4), 10, dimnames = list(letters[1:10], NULL))
  G <- genomic_relationship(M)
  expect_equal(rownames(G), letters[1:10])
  expect_equal(colnames(G), letters[1:10])
})

test_that("missing genotypes are handled and constant markers dropped", {
  set.seed(2)
  M <- matrix(rbinom(30 * 100, 2, 0.4), 30)
  M[sample(length(M), 50)] <- NA
  M[, 1] <- 1
  G <- genomic_relationship(M)
  expect_false(anyNA(G))
  expect_lt(attr(G, "n_markers"), 100)
})

test_that("rare and heavily missing markers are filtered", {
  set.seed(3)
  M <- matrix(rbinom(30 * 100, 2, 0.4), 30)
  M[, 2] <- 0
  M[1:10, 3] <- NA
  G <- genomic_relationship(M, maf = 0.05, max_missing = 0.1)
  expect_equal(attr(G, "n_markers"), 98)
})

test_that("user-supplied allele frequencies are used", {
  set.seed(4)
  M <- matrix(rbinom(30 * 100, 2, 0.4), 30)
  G1 <- genomic_relationship(M, p = rep(0.4, 100), maf = 0)
  G2 <- genomic_relationship(M, maf = 0)
  expect_false(isTRUE(all.equal(unname(G1), unname(G2))))
})

test_that("blending toward the identity makes the matrix positive definite", {
  set.seed(5)
  M <- matrix(rbinom(50 * 20, 2, 0.4), 50)
  G <- genomic_relationship(M, blend = 0.05)
  expect_gt(min(eigen(unname(G), symmetric = TRUE)$values), 0)
})

test_that("bad inputs are rejected", {
  M <- matrix(rbinom(20 * 30, 2, 0.4), 20)
  expect_error(genomic_relationship("a"), "numeric")
  expect_error(genomic_relationship(M[1, , drop = FALSE]), "at least two")
  expect_error(genomic_relationship(M + 1), "between 0 and 2")
  expect_error(genomic_relationship(M, maf = 0.6), "maf")
  expect_error(genomic_relationship(M, max_missing = 2), "max_missing")
  expect_error(genomic_relationship(M, blend = 1), "blend")
  expect_error(genomic_relationship(M, p = rep(0.5, 3)), "p must")
  expect_error(genomic_relationship(matrix(1, 5, 5)), "No markers left")
})
