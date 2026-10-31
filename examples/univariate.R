# Univariate animal model example
#
# lmm_animal fits univariate models under two parametrisations:
#   "variance": half-normal priors on the genetic and residual standard deviations,
#               h2 derived as s_G^2 / (s_G^2 + s_E^2)
#   "h2":       half-normal prior on the total standard deviation and a Beta prior
#               on h2 (beta_prior, c(1, 1) is flat)
# Both are fitted here and compared. Run from the project root:
#   Rscript examples/univariate.R

# Install (if needed) and load all dependencies in one go
pkgs <- c("rstan", "nadiv", "bayesplot", "ggplot2",
          "cran/pedtools")  # pedtools is no longer on CRAN -> install from GitHub
if(!requireNamespace("pak", quietly = TRUE)) install.packages("pak")
pak::pkg_install(pkgs, upgrade = FALSE)
invisible(lapply(sub("^.*/", "", pkgs), library, character.only = TRUE))

# cmdstanr is not on CRAN; install it only if missing
if(!requireNamespace("cmdstanr", quietly = TRUE)) pak::pkg_install("stan-dev/cmdstanr")
library(cmdstanr)

rstan_options(auto_write = TRUE)
options(mc.cores = 4)

library(stanAnimal)

## Simulation
ped = randomPed(1500, 50, seed = 2)
ped = data.frame(ID = ped$ID, sire = ped$FIDX, dam = ped$MIDX)
ped = prepPed(ped)
A <- as.matrix(nadiv::makeA(ped))
L_A = chol(A)

n = nrow(A)

set.seed(41)
h2_true = 0.6
P = 1
G = h2_true * P        # additive genetic variance
E = (1 - h2_true) * P  # residual variance

a = sqrt(G) * (t(L_A) %*% rnorm(n))  # Var(a) = G A

beta = matrix(c(1, 0.3), 2, 1, dimnames = list(c("Intercept", "sex"), "trait"))
X = cbind(Intercept = 1, sex = sample(c(0, 1), n, replace = TRUE))
e = rnorm(n, sd = sqrt(E))

Y = X %*% beta + a + e
colnames(Y) = "trait"

# Simulated truth, used as the reference line in the interval plots
truth = c(h2 = h2_true, sigma_G = sqrt(G), sigma_E = sqrt(E))
round(truth, 3)

## Model fit: variance parametrisation
# Half-normal(0, 1) priors on the standardized sigma_G and sigma_E; no prior on h2
fit_var = lmm_animal(Y, X, A, parametrization = "variance",
                     iter = 3000, warmup = 2000, cores = 4, chains = 4)

## Model fit: h2 parametrisation
# Half-normal(0, 1) prior on the standardized total sd and h2 ~ Beta(beta_prior)
fit_h2 = lmm_animal(Y, X, A, parametrization = "h2", beta_prior = c(1, 1),
                    iter = 3000, warmup = 2000, cores = 4, chains = 4)

# Both parametrisations report the same quantities on the original scale:
# h2, sigma_G, sigma_E (derived in generated quantities), plus the fixed effects
fit_var$summary(c("h2", "sigma_G", "sigma_E", "beta"))
fit_h2$summary(c("h2", "sigma_G", "sigma_E", "beta"))

# Fixed effects are named after the columns of X
beta_draws_var = animal_draws(fit_var, "beta")
beta_draws_h2 = animal_draws(fit_h2, "beta")
posterior::summarise_draws(beta_draws_var)
posterior::summarise_draws(beta_draws_h2)

# Posterior means compared with the simulated truth
post_means = function(fit) {
  s = fit$summary(c("h2", "sigma_G", "sigma_E"), "mean")
  setNames(s$mean, s$variable)
}
rbind(variance = post_means(fit_var), h2 = post_means(fit_h2), truth = truth)

# Fixed effects, original scale
round(setNames(fit_var$summary("beta", "mean")$mean, colnames(X)), 3)
t(beta)

## Diagnostics
# E-BFMI is where the two parametrisations can differ: when the genetic variance
# is weakly identified the variance version mixes poorly in the neck between the
# two variance components, and lmm_animal warns and suggests parametrization = "h2"
rbind(variance = fit_var$diagnostic_summary(quiet = TRUE)$ebfmi,
      h2 = fit_h2$diagnostic_summary(quiet = TRUE)$ebfmi)
fit_var$diagnostic_summary()
fit_h2$diagnostic_summary()

## Sampler diagnostics (divergences, treedepth, E-BFMI)
fit_var$summary(c("h2", "sigma_G", "sigma_E"),
                posterior::default_convergence_measures())
fit_h2$summary(c("h2", "sigma_G", "sigma_E"),
               posterior::default_convergence_measures())

## Visualization
# All figures are written to a single PDF
dir.create("examples/figures", recursive = TRUE, showWarnings = FALSE)
pdf("examples/figures/univariate_plots.pdf", width = 10, height = 7)
bayesplot::color_scheme_set("blue")

# Traces for both parametrisations
print(bayesplot::mcmc_trace(fit_var$draws(c("sigma_G", "sigma_E", "h2"))))
print(bayesplot::mcmc_trace(fit_h2$draws(c("sigma_G", "sigma_E", "h2"))))

# Convergence and sampler checks
sampler_checks = function(fit) {
  all_draws = fit$draws()
  print(bayesplot::mcmc_rhat(posterior::summarise_draws(all_draws, "rhat")$rhat))
  print(bayesplot::mcmc_neff(posterior::summarise_draws(all_draws, "ess_bulk")$ess_bulk /
                               posterior::ndraws(all_draws)))
  print(bayesplot::mcmc_nuts_energy(bayesplot::nuts_params(fit)))
}
sampler_checks(fit_var)
sampler_checks(fit_h2)

# Intervals with the simulated truth. mcmc_recover_intervals() matches the named
# true values to the draws and marks them with hollow triangles.
interval_plot = function(fit, pars, truth) {
  bayesplot::mcmc_recover_intervals(
    fit$draws(pars), true = stats::setNames(truth, pars),
    prob = 0.8, prob_outer = 0.99, point_est = "mean")
}
print(interval_plot(fit_var, c("h2", "sigma_G", "sigma_E"), truth))
print(interval_plot(fit_h2, c("h2", "sigma_G", "sigma_E"), truth))

# Fixed effects
print(bayesplot::mcmc_intervals(beta_draws_var))

# Prior versus posterior for the variance parametrisation
# (the prior animal_prior_draws uses: half-normal(0, 1) scales on both SDs)
print(plot_prior_posterior(fit_var, trait_names = "trait", ncol = 1))
print(plot_forest(fit_var, trait_names = "trait"))

# Prior versus posterior for h2 under both parametrisations
# variance: h2 is derived from two half-normal SDs, so its prior is not a Beta
# h2:       h2 ~ Beta(1, 1), i.e. flat
h2_draws = rbind(
  data.frame(h2 = as.vector(posterior::as_draws_matrix(fit_var$draws("h2"))),
             parametrization = "variance"),
  data.frame(h2 = as.vector(posterior::as_draws_matrix(fit_h2$draws("h2"))),
             parametrization = "h2"))
print(ggplot2::ggplot(h2_draws,
                      ggplot2::aes(h2, colour = parametrization, fill = parametrization)) +
        ggplot2::geom_density(alpha = 0.2) +
        ggplot2::geom_vline(xintercept = h2_true, linetype = "dashed") +
        ggplot2::stat_function(fun = stats::dbeta, args = list(shape1 = 1, shape2 = 1),
                               linetype = "dotted", colour = "grey30",
                               inherit.aes = FALSE) +
        ggplot2::labs(x = "h2", y = "Density",
                      title = "Posterior h2 under both parametrisations (dashed = truth, dotted = Beta(1,1) prior)"))

# The same comparison for all three reported quantities
compare = rbind(
  data.frame(parametrization = "variance",
             as.data.frame(as.matrix(
               posterior::as_draws_matrix(fit_var$draws(c("h2", "sigma_G", "sigma_E")))))),
  data.frame(parametrization = "h2",
             as.data.frame(as.matrix(
               posterior::as_draws_matrix(fit_h2$draws(c("h2", "sigma_G", "sigma_E")))))))
compare_long = rbind(
  data.frame(parameter = "h2", value = compare$h2, parametrization = compare$parametrization),
  data.frame(parameter = "sigma_G", value = compare$sigma_G, parametrization = compare$parametrization),
  data.frame(parameter = "sigma_E", value = compare$sigma_E, parametrization = compare$parametrization))
compare_long$parameter = factor(compare_long$parameter,
                                levels = c("h2", "sigma_G", "sigma_E"))
truth_df = data.frame(parameter = factor(names(truth), levels = levels(compare_long$parameter)),
                      truth = unname(truth))
print(ggplot2::ggplot(compare_long,
                      ggplot2::aes(value, colour = parametrization, fill = parametrization)) +
        ggplot2::geom_density(alpha = 0.2) +
        ggplot2::geom_vline(data = truth_df, ggplot2::aes(xintercept = truth),
                            linetype = "dashed", inherit.aes = FALSE) +
        ggplot2::facet_wrap(~ parameter, scales = "free") +
        ggplot2::labs(x = NULL, y = "Density",
                      title = "Posterior under both parametrisations (dashed = simulated truth)"))

dev.off()
