# Renames in 0.4.0 replaced the old exports rather than deprecating them, so a
# script written against an older version fails at the call. This keeps the old
# names from creeping back as aliases.

test_that("renamed functions replace the old exports", {
  ns <- asNamespace("stanAnimal")
  expect_true(exists("rename_beta", envir = ns, inherits = FALSE))
  expect_true(exists("pc1_variance_share", envir = ns, inherits = FALSE))
  expect_true(exists("lm_phen_cov", envir = ns, inherits = FALSE))
  expect_false(exists("animal_draws", envir = ns, inherits = FALSE))
  expect_false(exists("pc_dominance", envir = ns, inherits = FALSE))
  expect_false(exists("lmm_multiCov", envir = ns, inherits = FALSE))
})

test_that("the old names are not exported as aliases either", {
  exports <- getNamespaceExports("stanAnimal")
  expect_false(any(c("animal_draws", "pc_dominance", "lmm_multiCov") %in% exports))
})
