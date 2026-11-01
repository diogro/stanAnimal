# Input checks shared across the package.
#
# Every check follows the same shape: one helper per kind of value, the argument
# name in backticks, and `call. = FALSE` so the message names the argument the
# caller passed rather than the internal call that raised it. Public functions
# validate at the top, before any work, and reuse these helpers instead of
# repeating conditions inline.

`%||%` <- function(a, b) if (is.null(a)) b else a

# A matrix of finite numbers, with NA allowed only where the data can have
# missing values (marker genotypes). Vectors become one-column matrices.
.check_numeric_matrix <- function(x, arg, allow_na = FALSE) {
  x <- as.matrix(x)
  if (!is.numeric(x) || any(is.infinite(x)) || (!allow_na && anyNA(x)))
    stop("`", arg, "` must contain only finite numeric values",
         if (allow_na) ", or NA for missing", ".", call. = FALSE)
  x
}

.check_flag <- function(x, arg) {
  if (!is.logical(x) || length(x) != 1L || is.na(x))
    stop("`", arg, "` must be TRUE or FALSE.", call. = FALSE)
  x
}

# A single whole number, at least `min`: counts, numbers of draws, iterations.
.check_whole <- function(x, arg, min = 1) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < min || x != trunc(x))
    stop("`", arg, "` must be one whole number, at least ", min, ".",
         call. = FALSE)
  # beyond this the coercion to integer is NA, which callers would use as a size
  if (x > .Machine$integer.max)
    stop("`", arg, "` must be at most ", .Machine$integer.max, ".", call. = FALSE)
  as.integer(x)
}

.check_count <- function(x, arg) .check_whole(x, arg, 1)

# A single probability, strictly inside (0, 1).
.check_prob <- function(x, arg) .check_in_range(x, arg, 0, 1, TRUE, TRUE)

# A single finite number above zero: scale parameters such as the LKJ shape.
.check_positive <- function(x, arg) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0)
    stop("`", arg, "` must be one finite number above 0.", call. = FALSE)
  as.numeric(x)
}

# A single number inside a range, with either end optionally excluded.
.check_in_range <- function(x, arg, low, high, open_low = FALSE,
                            open_high = FALSE) {
  ok <- is.numeric(x) && length(x) == 1L && is.finite(x)
  if (ok) ok <- if (open_low) x > low else x >= low
  if (ok) ok <- if (open_high) x < high else x <= high
  if (!ok)
    stop("`", arg, "` must be one number in ", if (open_low) "(" else "[", low,
         ", ", high, if (open_high) ")" else "]", ".", call. = FALSE)
  as.numeric(x)
}

# Iterations of a sampler: whole numbers, at least one draw retained.
.check_iterations <- function(iter, warmup, thin = NULL) {
  iter <- .check_whole(iter, "iter", 1)
  warmup <- .check_whole(warmup, "warmup", 0)
  if (warmup >= iter)
    stop("`warmup` must be smaller than `iter` (", iter, ").", call. = FALSE)
  if (!is.null(thin)) {
    thin <- .check_whole(thin, "thin", 1)
    if (warmup + thin > iter)
      stop("`iter` must leave at least one draw after `warmup`: ",
           warmup + thin, " is above ", iter, ".", call. = FALSE)
  }
  invisible(list(iter = iter, warmup = warmup, thin = thin))
}

.check_sampler_size <- function(chains, cores) {
  chains <- .check_count(chains, "chains")
  cores <- .check_count(cores, "cores")
  if (cores > chains)
    stop("`cores` must not exceed `chains`: there is one core per chain at ",
         "most, and `cores` is ", cores, " for ", chains, " chains.",
         call. = FALSE)
  invisible(list(chains = chains, cores = cores))
}
