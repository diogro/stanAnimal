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
