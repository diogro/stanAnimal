# The README is not part of the built package, so these tests only run from a
# checkout of the repository (as in devtools::test()) and skip otherwise.

readme_source <- function() {
  candidates <- c(testthat::test_path("..", "..", "..", "README.md"),
                  testthat::test_path("..", "..", "README.md"))
  candidates <- candidates[file.exists(candidates)]
  if (!length(candidates)) return(NULL)
  normalizePath(candidates[[1L]])
}

# R code blocks of the README, with the line where each one starts. Blocks that
# install software are instructions rather than examples, so they are skipped.
readme_code_blocks <- function(path) {
  lines <- readLines(path, warn = FALSE)
  fences <- grep("^```", lines)
  out <- list()
  for (k in seq_len(length(fences) %/% 2L)) {
    open <- fences[2L * k - 1L]
    close <- fences[2L * k]
    if (!grepl("^```r\\s*$", lines[open])) next
    code <- lines[seq.int(open + 1L, close - 1L)]
    if (any(grepl("pkg_install\\(|install\\.packages\\(|install_cmdstan\\(", code)))
      next
    out[[length(out) + 1L]] <- list(line = open + 1L, code = paste(code, collapse = "\n"))
  }
  out
}

# The data and objects the snippets expect, plus short fits so the blocks are
# cheap: the snippets never pass iter or chains themselves.
readme_context <- function() {
  env <- new.env(parent = globalenv())
  short_animal <- function(...)
    stanAnimal::lmm_animal(..., iter = 200, warmup = 100, chains = 2, cores = 2,
                           refresh = 0, show_messages = FALSE)
  short_multiCov <- function(...)
    stanAnimal::lmm_multiCov(..., iter = 200, warmup = 100, chains = 2, cores = 2,
                             refresh = 0, show_messages = FALSE)
  assign("lmm_animal", short_animal, envir = env)
  assign("lmm_multiCov", short_multiCov, envir = env)

  local({
    set.seed(42)
    n <- 60
    ids <- as.character(seq_len(n))
    # founders 1-12, then 48 offspring of a sire and a dam from different pairs
    offspring <- seq.int(13L, n)
    ped <- data.frame(id = seq_len(n),
                      sire = c(rep(NA_integer_, 12L), 1L + (offspring - 13L) %% 6L),
                      dam = c(rep(NA_integer_, 12L), 7L + (offspring - 13L) %% 6L))
    pedigree <- ped
    family <- rep(seq_len(12), each = 5)
    A0 <- outer(family, family, "==") * 0.5
    diag(A0) <- 1
    dimnames(A0) <- list(ids, ids)
    A <- A0
    M <- matrix(stats::rbinom(n * 200, 2, 0.3), n, dimnames = list(ids, NULL))
    dat <- data.frame(t1 = stats::rnorm(n), t2 = stats::rnorm(n), t3 = stats::rnorm(n),
                      age = stats::rnorm(n, 50, 10), sex = stats::rbinom(n, 1, 0.5),
                      site = sample(c("north", "south"), n, TRUE), id = seq_len(n))
    rownames(dat) <- dat$id
    morphology <- dat[, c("age", "sex", "site", "id")]
    log_traits <- matrix(log(stats::rlnorm(n * 5, 3, 0.4)), n, 5,
                         dimnames = list(ids, paste0("m", seq_len(5))))
    X <- cbind(Intercept = 1, age = dat$age, sex = dat$sex)
    B <- cbind(c(1, 0.1, 0.5), c(-1, -0.1, 0.3))
    Y <- X %*% B +
      t(chol(A0)) %*% matrix(stats::rnorm(n * 2), n) %*% diag(c(0.6, 0.6)) +
      matrix(stats::rnorm(n * 2), n) %*% diag(c(0.5, 0.5))
    colnames(Y) <- c("t1", "t2")
    y <- Y[, 1]
    assign("ped", ped, envir = env)
    assign("pedigree", pedigree, envir = env)
    assign("A", A, envir = env)
    assign("M", M, envir = env)
    assign("dat", dat, envir = env)
    assign("morphology", morphology, envir = env)
    assign("log_traits", log_traits, envir = env)
    assign("X", X, envir = env)
    assign("Y", Y, envir = env)
    assign("y", y, envir = env)
  })

  env
}

test_that("every example in the README runs", {
  path <- readme_source()
  if (is.null(path))
    skip("README.md is not available outside a checkout of the repository")
  skip_if_no_cmdstan()
  skip_if_not_installed("bayesplot")
  skip_if_not_installed("posterior")

  blocks <- readme_code_blocks(path)
  expect_gt(length(blocks), 10L)

  env <- readme_context()
  for (block in blocks) {
    failure <- tryCatch({
      utils::capture.output(eval(parse(text = block$code), envir = env))
      NULL
    }, error = function(e) conditionMessage(e))
    expect_true(is.null(failure),
                info = sprintf("README code block at line %d: %s", block$line,
                               if (is.null(failure)) "" else failure))
  }
})
