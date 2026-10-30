# Helpers for fitting Stan models with cmdstanr.
# Adapted from rethinking (https://github.com/rmcelreath/rethinking).

# Fit a model from Stan code using cmdstanr and return the CmdStanMCMC object.
# The result works directly with posterior (fit$draws()) and bayesplot
# (e.g. bayesplot::mcmc_trace(fit$draws("beta")), bayesplot::nuts_params(fit)).
cstan <- function(model_code, data = list(), chains = 1, cores = 1, iter = 1000,
                  warmup = floor(iter / 2), threads = 1,
                  control = list(adapt_delta = 0.95), cpp_options = list(), ...) {
  if (threads > 1) cpp_options[["stan_threads"]] <- TRUE

  stan_file <- cmdstanr::write_stan_file(model_code, dir = tempdir())
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
