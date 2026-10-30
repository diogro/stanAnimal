# Helpers for fitting Stan models with cmdstanr.
# Adapted from rethinking (https://github.com/rmcelreath/rethinking).

# Fit a model from Stan code using cmdstanr and return the CmdStanMCMC object.
# The result works directly with posterior (fit$draws()) and bayesplot
# (e.g. bayesplot::mcmc_trace(fit$draws("beta")), bayesplot::nuts_params(fit)).
cstan <- function(model_code, data = list(), chains = 1, cores = 1, iter = 1000,
                  warmup = floor(iter / 2), threads = 1,
                  control = list(adapt_delta = 0.95), cpp_options = list(), ...) {
  if (threads > 1) cpp_options[["stan_threads"]] <- TRUE

  # File names are hashes of the code, so cmdstanr reuses the compiled
  # executable across R sessions and recompiles only when the code changes.
  stan_file <- cmdstanr::write_stan_file(model_code, dir = stan_cache_dir())
  mod <- cmdstanr::cmdstan_model(stan_file, cpp_options = cpp_options)

  adapt_delta <- if (is.null(control[["adapt_delta"]])) 0.95 else
    as.numeric(control[["adapt_delta"]])
  max_treedepth <- if (is.null(control[["max_treedepth"]])) 11 else
    as.numeric(control[["max_treedepth"]])

  args <- list(data = data,
               chains = chains,
               parallel_chains = cores,
               iter_warmup = warmup,
               iter_sampling = iter - warmup,
               adapt_delta = adapt_delta,
               max_treedepth = max_treedepth)
  if (threads > 1) args$threads_per_chain <- threads
  args <- utils::modifyList(args, list(...))
  do.call(mod$sample, args)
}

#' Location of the compiled Stan model cache
#'
#' Models fit with \code{cmd_stan = TRUE} are compiled once and the executables
#' are reused in later R sessions. The cache lives in the directory given by
#' \code{options(stanAnimal.cache_dir = )} or the environment variable
#' \code{STANANIMAL_CACHE_DIR}, and defaults to the user cache directory
#' returned by \code{tools::R_user_dir("stanAnimal", "cache")}. Executables are
#' specific to the CmdStan version and platform; after upgrading CmdStan, call
#' \code{stanAnimal_clear_cache()}.
#' @return \code{stan_cache_dir()} returns the (created) cache directory.
#'   \code{stanAnimal_clear_cache()} deletes its contents and returns the
#'   number of files removed, invisibly.
#' @rdname stan_cache
#' @export
stan_cache_dir <- function() {
  dir <- getOption("stanAnimal.cache_dir")
  if (is.null(dir) || !nzchar(dir))
    dir <- Sys.getenv("STANANIMAL_CACHE_DIR", unset = "")
  if (!nzchar(dir))
    dir <- tools::R_user_dir("stanAnimal", "cache")
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  normalizePath(dir)
}

#' @rdname stan_cache
#' @export
stanAnimal_clear_cache <- function() {
  files <- list.files(stan_cache_dir(), full.names = TRUE, all.files = TRUE,
                      no.. = TRUE)
  unlink(files, recursive = TRUE)
  invisible(length(files))
}
