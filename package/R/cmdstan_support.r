# Helpers for fitting Stan models with cmdstanr.
# Adapted from rethinking (https://github.com/rmcelreath/rethinking).

# Sampler controls the package forwards to either engine. Names follow cmdstanr;
# rstan spells step_size "stepsize". Anything else is an error, so that a
# control accepted by one engine is never silently dropped by the other.
.stan_controls <- c("adapt_delta", "adapt_engaged", "max_treedepth", "step_size",
                    "init_buffer", "term_buffer", "window", "metric")

.check_control <- function(control) {
  if (is.null(control)) return(list())
  if (!is.list(control) || (length(control) && is.null(names(control))))
    stop("control must be a named list of sampler controls.", call. = FALSE)
  # rstan calls step_size "stepsize"; accept it and settle on a single name.
  names(control)[names(control) == "stepsize"] <- "step_size"
  dup <- unique(names(control)[duplicated(names(control))])
  if (length(dup))
    stop("control sets ", paste(dup, collapse = ", "), " more than once.",
         call. = FALSE)
  unknown <- setdiff(names(control), .stan_controls)
  if (length(unknown))
    stop("unsupported sampler controls: ", paste(unknown, collapse = ", "),
         ". Supported: ", paste(.stan_controls, collapse = ", "), ".",
         call. = FALSE)
  control
}

# The same controls in the spelling rstan uses.
.rstan_control <- function(control) {
  names(control)[names(control) == "step_size"] <- "stepsize"
  control
}

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

  control <- .check_control(control)

  args <- list(data = data,
               chains = chains,
               parallel_chains = cores,
               iter_warmup = warmup,
               iter_sampling = iter - warmup,
               adapt_delta = 0.95,
               max_treedepth = 11,
               init = 0)
  if (threads > 1) args$threads_per_chain <- threads
  args <- utils::modifyList(args, control)
  args <- utils::modifyList(args, list(...))
  do.call(mod$sample, args)
}

#' Location of the compiled Stan model cache
#'
#' Models fit with \code{cmd_stan = TRUE} are compiled once and the executables
#' are reused in later R sessions. The cache lives in the directory given by
#' \code{options(stanimal.cache_dir = )} or the environment variable
#' \code{STANIMAL_CACHE_DIR}, and defaults to the user cache directory
#' returned by \code{tools::R_user_dir("stanimal", "cache")}. Executables are
#' specific to the CmdStan version and platform; after upgrading CmdStan, call
#' \code{stanimal_clear_cache()}.
#' @return \code{stan_cache_dir()} returns the (created) cache directory.
#'   \code{stanimal_clear_cache()} deletes its contents and returns the
#'   number of files removed, invisibly.
#' @rdname stan_cache
#' @export
stan_cache_dir <- function() {
  dir <- getOption("stanimal.cache_dir")
  if (is.null(dir) || !nzchar(dir))
    dir <- Sys.getenv("STANIMAL_CACHE_DIR", unset = "")
  if (!nzchar(dir))
    dir <- tools::R_user_dir("stanimal", "cache")
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  normalizePath(dir)
}

#' @rdname stan_cache
#' @export
stanimal_clear_cache <- function() {
  files <- list.files(stan_cache_dir(), full.names = TRUE, all.files = TRUE,
                      no.. = TRUE)
  unlink(files, recursive = TRUE)
  invisible(length(files))
}
