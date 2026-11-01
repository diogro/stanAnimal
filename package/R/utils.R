# Small internal helpers shared across files.

`%||%` <- function(a, b) if (is.null(a)) b else a

# A single probability, strictly inside (0, 1).
.check_prob <- function(x, arg) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0 || x >= 1)
    stop("`", arg, "` must be one number strictly between 0 and 1.", call. = FALSE)
  as.numeric(x)
}

# A single whole number, at least 1: counts, numbers of draws, numbers of traits.
.check_count <- function(x, arg) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < 1 || x != trunc(x))
    stop("`", arg, "` must be one whole number, at least 1.", call. = FALSE)
  # beyond this the coercion to integer is NA, which callers would use as a size
  if (x > .Machine$integer.max)
    stop("`", arg, "` must be at most ", .Machine$integer.max, ".", call. = FALSE)
  as.integer(x)
}

# A single finite number above zero: scale parameters such as the LKJ shape.
.check_positive <- function(x, arg) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0)
    stop("`", arg, "` must be one finite number above 0.", call. = FALSE)
  as.numeric(x)
}
