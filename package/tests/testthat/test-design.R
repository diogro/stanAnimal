dat <- data.frame(trait = rnorm(8),
                  x = rnorm(8),
                  x2 = rnorm(8),
                  sex = rep(c(0, 1), 4),
                  id = 1:8)

test_that("animal_design builds Y and X from a formula", {
  d <- animal_design(trait ~ sex + x, dat)
  expect_true(is.matrix(d$Y))
  expect_equal(dim(d$Y), c(8L, 1L))
  expect_equal(colnames(d$Y), "trait")
  expect_equal(unname(d$Y[, 1]), dat$trait)
  expect_equal(colnames(d$X), c("(Intercept)", "sex", "x"))
  expect_equal(unname(d$X[, 3]), dat$x)
  expect_equal(attr(d, "individuals"), as.character(1:8))
})

test_that("animal_design handles several traits with cbind", {
  d <- animal_design(cbind(trait, x2) ~ sex, dat)
  expect_equal(colnames(d$Y), c("trait", "x2"))
  expect_equal(dim(d$Y), c(8L, 2L))
  expect_equal(unname(d$Y[, 2]), dat$x2)
  expect_equal(colnames(d$X), c("(Intercept)", "sex"))
})

test_that("animal_design ignores random-effect terms", {
  d <- animal_design(trait ~ sex + (1 | id), dat)
  expect_equal(colnames(d$X), c("(Intercept)", "sex"))
  expect_equal(dim(d$Y), c(8L, 1L))
  d2 <- animal_design(trait ~ (1 | id), dat)
  expect_equal(colnames(d2$X), "(Intercept)")
})

test_that("animal_design drops incomplete rows and reports the ones kept", {
  dat_na <- dat
  dat_na$x[3] <- NA
  d <- animal_design(trait ~ sex + x, dat_na)
  expect_equal(nrow(d$Y), 7L)
  expect_equal(attr(d, "individuals"), as.character(c(1, 2, 4, 5, 6, 7, 8)))
})

test_that("animal_design dummy codes factors", {
  d <- animal_design(trait ~ factor(sex), dat)
  expect_equal(colnames(d$X), c("(Intercept)", "factor(sex)1"))
})

test_that("animal_design validates its input", {
  expect_error(animal_design(dat$trait, dat), "model formula")
  expect_error(animal_design(~sex, dat), "left-hand side")
  expect_error(animal_design(trait ~ sex, data = NULL), "data must be")
  expect_error(animal_design(trait ~ sex + I(1 - sex), dat), "full rank")
  chars <- dat
  chars$trait <- letters[1:8]
  expect_error(animal_design(trait ~ sex, chars), "must be numeric")
  infin <- dat
  infin$trait[2] <- Inf
  expect_error(animal_design(trait ~ sex, infin), "finite")
})

ped6 <- data.frame(id = 1:6, sire = c(0, 0, 1, 1, 4, 4),
                   dam = c(0, 0, 2, 2, 3, 5))

named_dat <- function(pheno = rnorm(6)) {
  d <- data.frame(trait = pheno, x = rnorm(6), id = 1:6)
  rownames(d) <- d$id
  d
}

test_that("animal_design builds A from a pedigree", {
  d <- animal_design(trait ~ x, named_dat(), pedigree = ped6)
  expect_equal(unname(d$A), unname(makeA(ped6)))
  expect_equal(rownames(d$A), as.character(1:6))
  expect_equal(colnames(d$A), rownames(d$A))
  expect_true(isSymmetric(d$A))
  # full sibs are related by 0.5, and 5 and 6 are inbred
  expect_equal(unname(d$A["3", "4"]), 0.5)
  expect_equal(unname(diag(d$A)), c(1, 1, 1, 1, 1.25, 1.375))
})

test_that("the pedigree may contain ancestors without records", {
  ped <- data.frame(id = 1:8,
                    sire = c(0, 0, 0, 0, 1, 3, 5, 5),
                    dam = c(0, 0, 0, 0, 2, 4, 6, 7))
  dd <- data.frame(trait = rnorm(4), id = 5:8)
  rownames(dd) <- dd$id
  d <- animal_design(trait ~ 1, dd, pedigree = ped)
  expect_equal(dim(d$A), c(4L, 4L))
  expect_equal(rownames(d$A), as.character(5:8))
})

test_that("A follows the rows kept when phenotypes are missing", {
  pheno <- rnorm(6)
  pheno[3] <- NA
  d <- animal_design(trait ~ 1, named_dat(pheno), pedigree = ped6)
  expect_equal(nrow(d$Y), 5L)
  expect_equal(rownames(d$A), c("1", "2", "4", "5", "6"))
})

test_that("animal_design builds A from genotypes", {
  set.seed(2)
  M <- matrix(rbinom(6 * 60, 2, 0.3), 6,
              dimnames = list(as.character(1:6), NULL))
  d <- animal_design(trait ~ x, named_dat(), genotypes = M)
  expect_equal(d$A, genomic_relationship(M))
  d2 <- animal_design(trait ~ x, named_dat(), genotypes = M, blend = 0.01)
  expect_equal(d2$A, genomic_relationship(M, blend = 0.01))
})

test_that("genotypes may cover individuals without records", {
  set.seed(3)
  M <- matrix(rbinom(9 * 60, 2, 0.3), 9,
              dimnames = list(as.character(1:9), NULL))
  d <- animal_design(trait ~ 1, named_dat(), genotypes = M)
  ref <- genomic_relationship(M)[as.character(1:6), as.character(1:6)]
  expect_equal(unname(d$A), unname(ref), ignore_attr = TRUE)
  expect_equal(rownames(d$A), as.character(1:6))
})

test_that("animal_design validates pedigrees and genotypes", {
  dat <- named_dat()
  set.seed(4)
  M <- matrix(rbinom(6 * 60, 2, 0.3), 6,
              dimnames = list(as.character(1:6), NULL))
  expect_error(animal_design(trait ~ 1, dat, pedigree = ped6, genotypes = M),
               "either pedigree or genotypes")
  expect_error(animal_design(trait ~ 1, dat, blend = 0.1),
               "only used when genotypes")

  expect_error(animal_design(trait ~ 1, dat, pedigree = ped6[, 1:2]),
               "three columns")
  expect_error(animal_design(trait ~ 1, dat, pedigree = rbind(ped6, ped6[1, ])),
               "duplicated")
  bad <- dat
  rownames(bad) <- c("x1", as.character(2:6))
  expect_error(animal_design(trait ~ 1, bad, pedigree = ped6),
               "not in the pedigree")
  expect_error(animal_design(trait ~ 1, bad, genotypes = M),
               "not in the genotypes")
})
