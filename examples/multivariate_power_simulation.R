# Parallel simulation checker for genetic-correlation recovery in the
# marginalized multivariate animal model.
# Run from the project root with: Rscript check_multivariate_power.R
# Set STANANIMAL_POWER_METHOD=gibbs to fit the same datasets with lmm_animal_gibbs.

method <- Sys.getenv("STANANIMAL_POWER_METHOD", "stan")
if (!method %in% c("stan", "gibbs"))
  stop("STANANIMAL_POWER_METHOD must be 'stan' or 'gibbs'.")

if (method == "stan" && !requireNamespace("cmdstanr", quietly = TRUE))
  stop("Install cmdstanr before running this checker.")
if (!requireNamespace("posterior", quietly = TRUE))
  stop("Install posterior before running this checker.")
if (.Platform$OS.type == "windows")
  stop("This checker uses parallel::mclapply and currently requires Linux or macOS.")

suppressPackageStartupMessages(library(rstan))
source("package/R/cmdstan_support.r")
source("package/R/lmm_animal.R")
source("package/R/lmm_animal_gibbs.R")

n_traits <- 4L
sample_sizes <- c(48L, 96L, 192L)
n_replicates <- as.integer(Sys.getenv("STANANIMAL_POWER_REPS", "40"))
iterations <- as.integer(Sys.getenv("STANANIMAL_POWER_ITER", "1200"))
warmup <- as.integer(Sys.getenv("STANANIMAL_POWER_WARMUP", as.character(iterations / 2)))
default_workers <- parallel::detectCores()
if (is.na(default_workers)) default_workers <- 1L
workers <- as.integer(Sys.getenv(
  "STANANIMAL_POWER_WORKERS", as.character(min(32L, default_workers))))
seed <- as.integer(Sys.getenv("STANANIMAL_POWER_SEED", "20261030"))
adapt_delta <- as.numeric(Sys.getenv("STANANIMAL_POWER_ADAPT_DELTA", "0.99"))
output_dir <- Sys.getenv("STANANIMAL_POWER_OUTPUT",
                         if (method == "gibbs") "results/multivariate_power_gibbs" else "results/multivariate_power")
chains <- if (method == "stan") 4L else 1L
gibbs_thin <- as.integer(Sys.getenv("STANANIMAL_POWER_GIBBS_THIN", "1"))

if (anyNA(c(n_replicates, iterations, warmup, workers, seed)) ||
    n_replicates < 1L || iterations < 2L || warmup < 1L ||
    warmup >= iterations || workers < 1L || length(adapt_delta) != 1L ||
    !is.finite(adapt_delta) || adapt_delta <= 0.8 || adapt_delta >= 1)
  stop("Check STANANIMAL_POWER_REPS, ITER, WARMUP, WORKERS, and SEED.")

# A lower-triangular factor defines an SPD correlation matrix with high,
# moderate, weak/negative, and near-zero pairwise correlations.
G_factor <- rbind(
  c(0.75, 0, 0, 0),
  c(0.50, 0.40, 0, 0),
  c(0.25, 0.35, 0.30, 0),
  c(-0.12, 0.15, 0.40, 0.30)
)
G_correlation <- cov2cor(G_factor %*% t(G_factor))
G_sd <- c(0.80, 0.90, 0.75, 0.85)
E_sd <- c(0.90, 1.00, 0.85, 1.05)
G_true <- diag(G_sd) %*% G_correlation %*% diag(G_sd)
E_correlation <- outer(seq_len(n_traits), seq_len(n_traits),
                       function(i, j) 0.30 ^ abs(i - j))
E_true <- diag(E_sd) %*% E_correlation %*% diag(E_sd)
correlation_truth <- cov2cor(G_true)
h2_truth <- diag(G_true) / diag(G_true + E_true)

if (min(eigen(G_correlation, symmetric = TRUE, only.values = TRUE)$values) <= 0)
  stop("The genetic correlation truth must be positive definite.")

make_relationship <- function(n, family_size = 4L) {
  if (n %% family_size != 0L)
    stop("Sample size must be divisible by family_size.")
  A <- diag(n)
  for (start in seq.int(1L, n, by = family_size)) {
    members <- start:(start + family_size - 1L)
    A[members, members] <- 0.5
    diag(A)[members] <- 1
  }
  A
}

simulate_one <- function(n, replicate) {
  A <- make_relationship(n)
  X <- cbind(
    intercept = 1,
    treatment = stats::rbinom(n, 1, 0.5),
    age = stats::rnorm(n, mean = 50, sd = 10)
  )
  beta_true <- rbind(
    c(1.0, 0.30, 0.04),
    c(-0.5, -0.25, -0.03),
    c(0.7, 0.20, 0.02),
    c(0.2, 0.35, 0.01)
  )

  breeding_values <- t(chol(A)) %*%
    matrix(stats::rnorm(n * n_traits), n, n_traits) %*% chol(G_true)
  residuals <- matrix(stats::rnorm(n * n_traits), n, n_traits) %*% chol(E_true)
  Y <- X %*% t(beta_true) + breeding_values + residuals
  colnames(Y) <- paste0("trait", seq_len(n_traits))

  trait_pairs <- combn(n_traits, 2L, simplify = FALSE)
  if (method == "stan") {
    fit <- lmm_animal(
      Y = Y,
      X = X,
      A = A,
      lkj_prior = 2,
      iter = iterations,
      warmup = warmup,
      chains = 4,
      cores = 1,
      control = list(adapt_delta = adapt_delta),
      cmd_stan = TRUE,
      residual_structure = "unstructured"
    )
    draws <- posterior::as_draws_matrix(fit$draws(c("corrG", "h2")))
    corr_draws <- function(i, j) draws[, sprintf("corrG[%d,%d]", i, j)]
    h2_draws <- function(trait) draws[, sprintf("h2[%d]", trait)]
    fit_summary <- fit$summary(variables = "corrG")
    sampler <- data.frame(
      max_rhat_corrG = max(fit_summary$rhat, na.rm = TRUE),
      min_ess_bulk_corrG = min(fit_summary$ess_bulk, na.rm = TRUE),
      divergences = sum(fit$diagnostic_summary(quiet = TRUE)$num_divergent)
    )
  } else {
    fit <- lmm_animal_gibbs(Y, X, A, iter = iterations, warmup = warmup, thin = gibbs_thin)
    corr_array <- apply(fit$G, 3L, stats::cov2cor)
    corr_draws <- function(i, j) corr_array[(j - 1L) * n_traits + i, ]
    h2_draws <- function(trait) fit$h2[trait, ]
    ess <- vapply(trait_pairs, function(pair)
      posterior::ess_bulk(corr_draws(pair[1L], pair[2L])), numeric(1))
    sampler <- data.frame(
      max_rhat_corrG = NA_real_,
      min_ess_bulk_corrG = min(ess),
      divergences = 0L
    )
  }

  pair_rows <- lapply(trait_pairs, function(pair) {
    trait_i <- pair[1L]
    trait_j <- pair[2L]
    values <- corr_draws(trait_i, trait_j)
    lower <- unname(stats::quantile(values, 0.025))
    upper <- unname(stats::quantile(values, 0.975))
    truth <- correlation_truth[trait_i, trait_j]
    excludes_zero <- lower > 0 || upper < 0
    data.frame(
      n = n,
      replicate = replicate,
      trait_i = trait_i,
      trait_j = trait_j,
      truth = truth,
      estimate = mean(values),
      lower_95 = lower,
      upper_95 = upper,
      covered_95 = lower <= truth && truth <= upper,
      excludes_zero = excludes_zero,
      detects_direction = excludes_zero && sign(mean(values)) == sign(truth),
      stringsAsFactors = FALSE
    )
  })

  h2_rows <- lapply(seq_len(n_traits), function(trait) {
    values <- h2_draws(trait)
    lower <- unname(stats::quantile(values, 0.025))
    upper <- unname(stats::quantile(values, 0.975))
    data.frame(
      n = n,
      replicate = replicate,
      trait = trait,
      truth = h2_truth[trait],
      estimate = mean(values),
      covered_95 = lower <= h2_truth[trait] && h2_truth[trait] <= upper,
      stringsAsFactors = FALSE
    )
  })

  list(
    correlations = do.call(rbind, pair_rows),
    heritability = do.call(rbind, h2_rows),
    diagnostics = cbind(data.frame(n = n, replicate = replicate), sampler)
  )
}

# Precompile the exact source once before forking. lmm_animal/cstan writes this
# same code hash into the configured cache and reuses the executable per fit.
if (method == "stan") {
  model_code <- .read_animal_model_code(NULL)
  stan_file <- cmdstanr::write_stan_file(model_code, dir = stan_cache_dir())
  message("Precompiling marginalized multivariate animal model...")
  invisible(cmdstanr::cmdstan_model(stan_file))
}

tasks <- expand.grid(
  n = sample_sizes,
  replicate = seq_len(n_replicates),
  KEEP.OUT.ATTRS = FALSE
)
RNGkind("L'Ecuyer-CMRG")
set.seed(seed)
workers <- min(workers, nrow(tasks))
message(sprintf("Running %d %s fits across %d workers (%d chains per fit).",
                nrow(tasks), method, workers, chains))

fits <- parallel::mclapply(seq_len(nrow(tasks)), function(index) {
  task <- tasks[index, ]
  tryCatch(
    simulate_one(task$n, task$replicate),
    error = function(error) list(error = data.frame(
      n = task$n,
      replicate = task$replicate,
      message = conditionMessage(error),
      stringsAsFactors = FALSE
    ))
  )
}, mc.cores = workers, mc.set.seed = TRUE, mc.preschedule = FALSE)

failed <- vapply(fits, function(result) !is.null(result$error), logical(1))
if (any(failed)) {
  failures <- do.call(rbind, lapply(fits[failed], `[[`, "error"))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(failures, file.path(output_dir, "failures.csv"), row.names = FALSE)
  stop(sum(failed), " simulation fits failed; see failures.csv.")
}

correlations <- do.call(rbind, lapply(fits, `[[`, "correlations"))
heritability <- do.call(rbind, lapply(fits, `[[`, "heritability"))
diagnostics <- do.call(rbind, lapply(fits, `[[`, "diagnostics"))

group_key <- interaction(correlations$n, correlations$trait_i,
                         correlations$trait_j, drop = TRUE)
correlation_summary <- do.call(rbind, lapply(split(correlations, group_key), function(group) {
  true_value <- group$truth[1L]
  error <- group$estimate - true_value
  data.frame(
    n = group$n[1L],
    trait_i = group$trait_i[1L],
    trait_j = group$trait_j[1L],
    truth = true_value,
    mean_estimate = mean(group$estimate),
    bias = mean(error),
    rmse = sqrt(mean(error^2)),
    coverage_95 = mean(group$covered_95),
    power_direction = if (abs(true_value) >= 0.1) mean(group$detects_direction) else NA_real_,
    false_positive_rate = if (abs(true_value) < 0.05) mean(group$excludes_zero) else NA_real_,
    mean_interval_width = mean(group$upper_95 - group$lower_95),
    replicates = nrow(group)
  )
}))

h2_key <- interaction(heritability$n, heritability$trait, drop = TRUE)
h2_summary <- do.call(rbind, lapply(split(heritability, h2_key), function(group) {
  error <- group$estimate - group$truth[1L]
  data.frame(
    n = group$n[1L],
    trait = group$trait[1L],
    truth = group$truth[1L],
    mean_estimate = mean(group$estimate),
    bias = mean(error),
    rmse = sqrt(mean(error^2)),
    coverage_95 = mean(group$covered_95),
    replicates = nrow(group)
  )
}))

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(correlations, file.path(output_dir, "correlation_draw_summaries.csv"), row.names = FALSE)
utils::write.csv(correlation_summary, file.path(output_dir, "correlation_power_summary.csv"), row.names = FALSE)
utils::write.csv(heritability, file.path(output_dir, "heritability_draw_summaries.csv"), row.names = FALSE)
utils::write.csv(h2_summary, file.path(output_dir, "heritability_summary.csv"), row.names = FALSE)
utils::write.csv(diagnostics, file.path(output_dir, "sampler_diagnostics.csv"), row.names = FALSE)
utils::write.csv(data.frame(
  sample_sizes = paste(sample_sizes, collapse = ","),
  replicates = n_replicates,
  iterations = iterations,
  warmup = warmup,
  method = method,
  chains = chains,
  parallel_workers = workers,
  adapt_delta = adapt_delta,
  seed = seed,
  stringsAsFactors = FALSE
), file.path(output_dir, "run_config.csv"), row.names = FALSE)

print(correlation_summary, row.names = FALSE)
print(h2_summary, row.names = FALSE)
message("Results written to: ", normalizePath(output_dir))