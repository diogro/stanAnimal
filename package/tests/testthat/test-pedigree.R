ped6 <- data.frame(id = 1:6, sire = c(0, 0, 1, 1, 4, 4),
                   dam = c(0, 0, 2, 2, 3, 5))

# The tabular method, worked out by hand for ped6: individuals 5 and 6 are
# inbred through the full-sib mating of 4 with 3 and of 5 with 4.
expected_ped6 <- structure(
  c(1, 0, 0.5, 0.5, 0.5, 0.5,
    0, 1, 0.5, 0.5, 0.5, 0.5,
    0.5, 0.5, 1, 0.5, 0.75, 0.625,
    0.5, 0.5, 0.5, 1, 0.75, 0.875,
    0.5, 0.5, 0.75, 0.75, 1.25, 1,
    0.5, 0.5, 0.625, 0.875, 1, 1.375),
  dim = c(6L, 6L),
  dimnames = list(as.character(1:6), as.character(1:6)))

test_that("makeA computes the tabular relationship matrix", {
  A <- makeA(ped6)
  expect_equal(A, expected_ped6)
  expect_true(isSymmetric(A))
  expect_equal(unname(diag(A)), c(1, 1, 1, 1, 1.25, 1.375))
  expect_gt(min(eigen(A, symmetric = TRUE, only.values = TRUE)$values), 0)
})

test_that("the order of the pedigree does not matter", {
  shuffled <- ped6[c(5, 1, 6, 3, 2, 4), ]
  A <- makeA(shuffled)
  expect_equal(unname(A[rownames(expected_ped6), colnames(expected_ped6)]),
               unname(expected_ped6))
})

test_that("unknown parents can be coded 0, NA or empty", {
  from_na <- data.frame(id = 1:6, sire = c(NA, NA, 1, 1, 4, 4),
                        dam = c(NA, NA, 2, 2, 3, 5))
  from_empty <- ped6
  from_empty$sire <- as.character(from_empty$sire)
  from_empty$dam <- as.character(from_empty$dam)
  from_empty$sire[1:2] <- ""
  from_empty$dam[1:2] <- ""
  expect_equal(makeA(from_na), expected_ped6)
  expect_equal(makeA(from_empty), expected_ped6)
})

test_that("character IDs are kept", {
  ped <- data.frame(id = letters[1:4], sire = c(NA, NA, "a", "c"),
                    dam = c(NA, NA, "b", "b"))
  A <- makeA(ped)
  expect_identical(rownames(A), letters[1:4])
  expect_identical(colnames(A), letters[1:4])
  expect_equal(unname(A["c", "d"]), 0.75)
  expect_equal(unname(A["a", "d"]), 0.25)
})

test_that("a parent that is not in the pedigree is added as a founder", {
  ped <- data.frame(id = 1:3, sire = c(NA, 1, 9), dam = c(NA, NA, NA))
  expect_warning(A <- makeA(ped), "added as founders: 9")
  expect_identical(rownames(A), c("1", "2", "3", "9"))
  expect_equal(unname(A["3", "9"]), 0.5)
  expect_equal(unname(A["3", "1"]), 0)
})

test_that("makeA validates the pedigree", {
  loop <- data.frame(id = 1:2, sire = c(1, 1), dam = c(NA, 2))
  expect_error(makeA(loop), "loop")
  repeated <- data.frame(id = c(1, 1), sire = c(0, 0), dam = c(0, 0))
  expect_error(makeA(repeated), "duplicated")
  expect_error(makeA(data.frame(id = 1:2, sire = c(0, 0))),
               "three columns")
  no_id <- data.frame(id = c(1, NA), sire = c(0, 0), dam = c(0, 0))
  expect_error(makeA(no_id), "missing ID")
})
