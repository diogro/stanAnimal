test_that("ratones_morpho has the documented shape", {
  expect_s3_class(ratones_morpho, "data.frame")
  expect_identical(dim(ratones_morpho), c(329L, 46L))
  expect_identical(rownames(ratones_morpho), as.character(ratones_morpho$ID))
  expect_identical(colnames(ratones_morpho)[1:11],
                   c(".id", "ID", "SEX", "GER", "LIN", "MADRE", "PADRE",
                     "AGE", "P49", "line", "selection"))
  expect_length(names(ratones_morpho)[12:46], 35L)
  expect_false(anyNA(ratones_morpho))
  expect_setequal(unique(ratones_morpho$selection),
                  c("control", "downwards", "upwards"))
  expect_true(all(ratones_morpho$SEX %in% c("F", "M")))
})

test_that("ratones_pedigree has the documented shape and covers the mice", {
  expect_s3_class(ratones_pedigree, "data.frame")
  expect_identical(names(ratones_pedigree), c("id", "sire", "dam"))
  expect_identical(dim(ratones_pedigree), c(525L, 3L))
  expect_equal(anyDuplicated(ratones_pedigree$id), 0L)
  expect_true(all(ratones_morpho$ID %in% ratones_pedigree$id))
  # founders and unknown parents are NA
  expect_true(all(is.na(ratones_pedigree$sire) == is.na(ratones_pedigree$dam)))
  expect_gt(sum(is.na(ratones_pedigree$sire)), 0)
})

test_that("the bundled pedigree builds a usable relationship matrix", {
  skip_if_not_installed("nadiv")
  dat <- ratones_morpho[c(1:20, 300:329), ]
  dat$trait <- dat[["IS_PM"]]
  d <- animal_design(trait ~ SEX + AGE, dat, pedigree = ratones_pedigree)
  expect_identical(dim(d$A), c(nrow(dat), nrow(dat)))
  expect_identical(rownames(d$A), as.character(dat$ID))
  expect_identical(colnames(d$A), rownames(d$A))
  expect_true(isSymmetric(d$A))
  # every measured mouse is related to itself, and the pedigree has ancestors
  expect_equal(unname(diag(d$A)), rep(1, nrow(dat)))
  expect_true(any(d$A[upper.tri(d$A)] > 0))
  expect_gt(min(eigen(d$A, symmetric = TRUE, only.values = TRUE)$values), 0)
})

test_that("the morphometrics are moderately integrated", {
  traits <- names(ratones_morpho)[12:46]
  Y <- log(as.matrix(ratones_morpho[, traits]))
  dominance <- pc_dominance(Y)
  expect_gt(dominance, 0.3)
  expect_lt(dominance, 0.8)
  expect_gt(mean(cor(Y)[upper.tri(diag(ncol(Y)))]), 0.2)
})
