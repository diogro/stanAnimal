set.seed(11)
n <- 40
size <- rnorm(n)
Y4 <- cbind(a = size + rnorm(n, sd = 0.3),
            b = 2 * size + rnorm(n, sd = 0.3),
            c = -size + rnorm(n, sd = 0.3),
            d = size + rnorm(n, sd = 0.5))
colnames(Y4) <- c("a", "b", "c", "d")

names_mat <- function(var, k) {
  sprintf(paste0(var, "[%d,%d]"), rep(seq_len(k), k), rep(seq_len(k), each = k))
}

test_that("size_shape returns size plus orthogonal shape axes", {
  ss <- size_shape(Y4)
  expect_s3_class(ss, "size_shape")
  expect_equal(colnames(ss$scores), c("size", "shape1", "shape2", "shape3"))
  expect_equal(ss$raw_traits, colnames(Y4))
  expect_equal(dim(ss$loadings), c(4L, 4L))
  expect_equal(unname(ss$scores[, "size"]), unname(rowMeans(Y4)))
  expect_equal(cor(ss$scores[, -1]), diag(3), ignore_attr = TRUE, tolerance = 1e-8)
  expect_output(print(ss), "size_shape")
})

test_that("size and shape reconstruct the traits up to their means", {
  ss <- size_shape(Y4)
  reconstructed <- ss$scores %*% t(ss$loadings)
  expect_equal(scale(reconstructed, scale = FALSE), scale(Y4, scale = FALSE),
               ignore_attr = TRUE)
})

test_that("n_shape keeps fewer shape axes", {
  ss <- size_shape(Y4, n_shape = 2)
  expect_equal(colnames(ss$scores), c("size", "shape1", "shape2"))
  expect_equal(dim(ss$loadings), c(4L, 3L))
  expect_equal(unname(ss$loadings[, "size"]), rep(1, 4))
})

test_that("size_shape validates its input", {
  expect_error(size_shape(Y4[, 1]), "at least two traits")
  expect_error(size_shape(matrix("a", 5, 2)), "numeric")
  expect_error(size_shape(Y4, n_shape = 0), "n_shape")
  expect_error(size_shape(Y4, n_shape = 4), "n_shape")
  expect_error(size_shape(Y4, n_shape = 1.5), "n_shape")
})

test_that("back_transform maps G and E and recomputes corrG and h2", {
  ss <- size_shape(Y4)
  k <- 4
  G <- matrix(0.3, k, k); diag(G) <- 1
  E <- matrix(0.1, k, k); diag(E) <- 2
  m <- cbind(matrix(rep(as.vector(G), 3), 3, byrow = TRUE),
             matrix(rep(as.vector(E), 3), 3, byrow = TRUE))
  colnames(m) <- c(names_mat("G", k), names_mat("E", k))

  post <- back_transform(m, ss)
  B <- ss$loadings
  exp_G <- unname(B %*% G %*% t(B))
  exp_E <- unname(B %*% E %*% t(B))
  exp_P <- exp_G + exp_E

  expect_s3_class(post, "size_shape_draws")
  expect_named(post, c("G", "E", "P", "corrG", "corrP", "h2", "traits"))
  expect_equal(dim(post$G), c(k, k, 3L))
  expect_identical(dimnames(post$G), list(ss$raw_traits, ss$raw_traits, NULL))
  expect_equal(unname(post$G[, , 1]), exp_G)
  expect_equal(unname(post$E[, , 1]), exp_E)
  expect_equal(unname(post$P[, , 1]), exp_P)
  expect_equal(unname(post$corrG[, , 1]), unname(cov2cor(exp_G)))
  expect_equal(unname(post$corrP[, , 1]), unname(cov2cor(exp_P)))
  expect_equal(unname(post$h2[, 1]), unname(diag(exp_G) / diag(exp_P)))
  expect_output(print(post), "size_shape_draws")
})

test_that("pc1_variance_share detects a single gradient", {
  set.seed(7)
  gradient <- rnorm(50)
  correlated <- cbind(a = gradient + rnorm(50, sd = 0.2),
                      b = 2 * gradient + rnorm(50, sd = 0.2),
                      c = -gradient + rnorm(50, sd = 0.2))
  expect_gt(pc1_variance_share(correlated), 0.95)
  expect_lt(pc1_variance_share(matrix(rnorm(150), 50, 3)), 0.7)
  expect_equal(pc1_variance_share(rnorm(50)), 1)
  expect_error(pc1_variance_share(matrix("a", 5, 2)), "numeric")
  expect_error(pc1_variance_share(cbind(1:5, rep(1, 5))), "must vary")
  expect_error(pc1_variance_share(matrix(0, 5, 3)), "must vary")
})

test_that("renamed functions replace the old exports", {
  ns <- asNamespace("stanAnimal")
  expect_true(exists("pc1_variance_share", envir = ns, inherits = FALSE))
  expect_true(exists("lm_phen_cov", envir = ns, inherits = FALSE))
  expect_false(exists("pc_dominance", envir = ns, inherits = FALSE))
  expect_false(exists("lmm_multiCov", envir = ns, inherits = FALSE))
})

test_that("animal_design warns about a dominant principal component", {
  set.seed(8)
  gradient <- rnorm(30)
  dominant <- data.frame(a = gradient + rnorm(30, sd = 0.1),
                         b = 2 * gradient + rnorm(30, sd = 0.1),
                         c = -gradient + rnorm(30, sd = 0.1))
  expect_warning(animal_design(cbind(a, b, c) ~ 1, dominant), "size_shape")
  independent <- data.frame(a = rnorm(30), b = rnorm(30), c = rnorm(30))
  expect_silent(animal_design(cbind(a, b, c) ~ 1, independent))
})

test_that("back_transform handles a phenotypic fit (P only)", {
  ss <- size_shape(Y4)
  k <- 4
  P <- matrix(0.2, k, k); diag(P) <- 1.5
  m <- matrix(rep(as.vector(P), 2), 2, byrow = TRUE)
  colnames(m) <- names_mat("P", k)

  post <- back_transform(m, ss)
  B <- ss$loadings
  exp_P <- unname(B %*% P %*% t(B))

  expect_named(post, c("P", "corrP", "sigma_P", "traits"))
  expect_equal(dim(post$P), c(k, k, 2L))
  expect_equal(dim(post$sigma_P), c(k, 2L))
  expect_equal(unname(post$P[, , 1]), exp_P)
  expect_equal(unname(post$corrP[, , 1]), unname(cov2cor(exp_P)))
  expect_equal(unname(post$sigma_P[, 1]), sqrt(diag(exp_P)))
})

test_that("back_transform works with fewer shape axes than traits", {
  ss <- size_shape(Y4, n_shape = 2)
  k_model <- ncol(ss$loadings)
  k <- nrow(ss$loadings)
  P <- matrix(0.2, k_model, k_model); diag(P) <- 1.5
  m <- matrix(rep(as.vector(P), 2), 2, byrow = TRUE)
  colnames(m) <- names_mat("P", k_model)

  post <- back_transform(m, ss)
  expect_equal(dim(post$P), c(k, k, 2L))
  expect_identical(dimnames(post$P), list(ss$raw_traits, ss$raw_traits, NULL))
  expect_equal(unname(post$P[, , 1]), unname(ss$loadings %*% P %*% t(ss$loadings)))
})

test_that("back_transform validates its input", {
  ss <- size_shape(Y4)
  expect_error(back_transform(matrix(0, 1, 1), list()), "size_shape object")
  m <- matrix(0, 2, 2, dimnames = list(NULL, c("G[1,1]", "E[1,1]")))
  expect_error(back_transform(m, ss), "must contain G and E")
  expect_error(back_transform(matrix(0, 2, 2), ss), "named matrix")
})

test_that("to_draws_list reshapes the arrays into posterior draws", {
  ss <- size_shape(Y4)
  k <- 4
  G <- matrix(0.3, k, k); diag(G) <- 1
  E <- matrix(0.1, k, k); diag(E) <- 2
  m <- cbind(matrix(rep(as.vector(G), 3), 3, byrow = TRUE),
             matrix(rep(as.vector(E), 3), 3, byrow = TRUE))
  colnames(m) <- c(names_mat("G", k), names_mat("E", k))
  post <- back_transform(m, ss)

  draws <- to_draws_list(post)
  expect_named(draws, c("G", "E", "P", "corrG", "corrP", "h2"))
  expect_identical(attr(draws, "traits"), ss$raw_traits)
  expect_s3_class(draws$G, "draws_matrix")
  expect_equal(dim(draws$G), c(3L, k * k))
  expect_identical(posterior::variables(draws$G), names_mat("G", k))
  expect_equal(as.numeric(as.matrix(draws$G)[1, ]), as.numeric(post$G[, , 1]))

  # correlations appear once per pair, above the diagonal
  pairs <- which(upper.tri(diag(k)), arr.ind = TRUE)
  pairs <- pairs[order(pairs[, 1], pairs[, 2]), , drop = FALSE]
  expect_identical(posterior::variables(draws$corrG),
                   sprintf("corrG[%d,%d]", pairs[, 1], pairs[, 2]))
  expect_equal(as.numeric(as.matrix(draws$corrG)[1, ]),
               as.numeric(post$corrG[cbind(pairs[, 1], pairs[, 2], rep(1L, nrow(pairs)))]))

  expect_identical(posterior::variables(draws$h2), sprintf("h2[%d]", seq_len(k)))
  expect_equal(as.numeric(as.matrix(draws$h2)[1, ]), as.numeric(post$h2[, 1]))

  # the draws objects can be summarised on their own or bound together
  expect_s3_class(posterior::summarise_draws(draws$corrG), "data.frame")
  # bind_draws() takes the draws objects as unnamed arguments
  expect_s3_class(do.call(posterior::bind_draws, unname(draws)), "draws_matrix")
})

test_that("to_draws_list handles a phenotypic fit and validates its input", {
  ss <- size_shape(Y4)
  k <- 4
  P <- matrix(0.2, k, k); diag(P) <- 1.5
  m <- matrix(rep(as.vector(P), 2), 2, byrow = TRUE)
  colnames(m) <- names_mat("P", k)
  post <- back_transform(m, ss)

  draws <- to_draws_list(post)
  expect_named(draws, c("P", "corrP", "sigma_P"))
  expect_identical(posterior::variables(draws$sigma_P),
                   sprintf("sigma_P[%d]", seq_len(k)))
  expect_equal(as.numeric(as.matrix(draws$sigma_P)[1, ]), as.numeric(post$sigma_P[, 1]))
  expect_error(to_draws_list(list()), "size_shape_draws")
})
