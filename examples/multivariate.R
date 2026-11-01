# Multivariate animal model example
#
# Simulates correlated traits on a random pedigree with a known genetic
# covariance matrix G, fits the marginalized multivariate model with
# lmm_animal(), and compares the posterior draws of G, E, P, corrG and h2 (and
# the fixed effects) with the simulated truth. Run from the project root:
#   Rscript examples/multivariate.R

# Install (if needed) and load all dependencies in one go
pkgs <- c("mvtnorm", "rstan", "evolqg", "corrplot", "patchwork", "bayesplot",
          "cran/pedtools")  # pedtools is no longer on CRAN -> install from GitHub
if(!requireNamespace("pak", quietly = TRUE)) install.packages("pak")
pak::pkg_install(pkgs, upgrade = FALSE)
invisible(lapply(sub("^.*/", "", pkgs), library, character.only = TRUE))

# cmdstanr is not on CRAN; install it only if missing
if(!requireNamespace("cmdstanr", quietly = TRUE)) pak::pkg_install("stan-dev/cmdstanr")
library(cmdstanr)

rstan_options(auto_write = TRUE)
options(mc.cores = 4)

library(stanimal)

## Simulation
ped = randomPed(1500, 50, seed = 2)
ped = data.frame(id = ped$ID, sire = ped$FIDX, dam = ped$MIDX)
A <- makeA(ped)
L_A = chol(A)

n_traits = 5
traits = paste0("trait", seq_len(n_traits))

set.seed(41)
# random_covariance(method = "partial") is the package port of
# evolqg::RandomMatrix(LKJ = FALSE), with the one- and two-trait cases fixed.
(G <- random_covariance(n_traits, method = "partial", trait_names = traits))
E = diag(n_traits)

n = nrow(A)
p = nrow(G)
a = t(L_A) %*% matrix(rnorm(n*p), n, p) %*% chol(G)
cor(a)

beta = matrix(c(rep(1, n_traits),
                rnorm(n_traits, 0, 1),
                rnorm(n_traits, 1, 1)), 3, n_traits, byrow = TRUE)
rownames(beta) = c("Intercept", "sex", "Z")
colnames(beta) = traits

sex = sample(c(0, 1), nrow(a), replace = TRUE)
Z = rnorm(nrow(a))
Intercept = rep(1, nrow(a))
X = cbind(Intercept, sex, Z)
e = rmvnorm(nrow(a), sigma = E)

Y = X %*% beta + a + e
colnames(Y) = traits

## Model fit
stan_model = lmm_animal(Y, X, A, lkj_prior = 2,
                        iter = 3000, warmup = 2000,
                        cores = 4, chains = 4)

# The fit is a CmdStanMCMC object: summaries and diagnostics use its methods,
# and draws work with posterior and bayesplot.

## Summaries
stan_model$summary(c("h2", "G", "E", "corrG"))
stan_model$summary(c("h2", "G", "E", "corrG"),
                   posterior::default_convergence_measures())

# All draws as trait x trait x draw arrays, with named dimensions
draws = fit_draws(stan_model)
draws

## Sampler diagnostics (divergences, treedepth, E-BFMI)
stan_model$diagnostic_summary()

# Fixed effects are named after the columns of X: beta[trait, column]
beta_draws = rename_beta(stan_model, "beta")
posterior::summarise_draws(beta_draws)

# Residual correlations are not a model output; fit_draws() derives them
# from the E draws, and the summary functions give matrices back
round(draws_mean(stan_model, "corrE"), 2)
round(draws_quantile(stan_model, "corrE", probs = c(0.05, 0.95)), 2)

# Posterior means as matrices, compared with the simulated truth
round(draws_mean(stan_model, "G"), 2); round(G, 2)
round(draws_mean(stan_model, "E"), 2); round(E, 2)
round(draws_mean(stan_model, "corrP"), 2); round(cov2cor(G + E), 2)

# Fixed effects, original scale: rows are traits, as in t(beta)
round(draws_mean(stan_model, "beta"), 2)
t(beta)

# Pairs of traits, for the interval plots below
pairs = which(lower.tri(diag(n_traits)), arr.ind = TRUE)

## Visualization
# All figures are written to a single PDF
dir.create("examples/figures", recursive = TRUE, showWarnings = FALSE)
pdf("examples/figures/multivariate_plots.pdf", width = 10, height = 7)
bayesplot::color_scheme_set("blue")

# Traces and densities
print(bayesplot::mcmc_trace(stan_model$draws("G")))
print(bayesplot::mcmc_trace(stan_model$draws("E")))
print(bayesplot::mcmc_dens_overlay(stan_model$draws("h2")))

# Convergence and sampler checks
all_draws = stan_model$draws()
print(bayesplot::mcmc_rhat(posterior::summarise_draws(all_draws, "rhat")$rhat))
print(bayesplot::mcmc_neff(posterior::summarise_draws(all_draws, "ess_bulk")$ess_bulk /
                             posterior::ndraws(all_draws)))
print(bayesplot::mcmc_nuts_energy(bayesplot::nuts_params(stan_model)))

# Fixed effects
print(bayesplot::mcmc_intervals(beta_draws))

# Intervals with the simulated truth. mcmc_recover_intervals() matches the named
# true values to the draws and marks them with hollow triangles.
interval_plot = function(pars, truth) {
  bayesplot::mcmc_recover_intervals(
    stan_model$draws(pars), true = stats::setNames(truth, pars),
    prob = 0.8, prob_outer = 0.99, point_est = "mean")
}
print(interval_plot(sprintf("h2[%d]", seq_len(n_traits)),
                    diag(G) / (diag(G) + diag(E))))
print(interval_plot(sprintf("G[%d,%d]", seq_len(n_traits), seq_len(n_traits)), diag(G)))
print(interval_plot(sprintf("E[%d,%d]", seq_len(n_traits), seq_len(n_traits)), diag(E)))
print(interval_plot(sprintf("G[%d,%d]", pairs[, 1], pairs[, 2]), G[pairs]))
print(interval_plot(sprintf("corrG[%d,%d]", pairs[, 1], pairs[, 2]), cov2cor(G)[pairs]))

# Prior versus posterior for heritabilities and genetic correlations
# (where they overlap, the data add little information)
print(plot_prior_posterior(stan_model, trait_names = traits, lkj_prior = 2,
                           ncol = min(3, n_traits)))
print(plot_forest(stan_model, trait_names = traits, lkj_prior = 2))

# Residuals after the fixed effects (X), with posterior draws of the G ellipses
# and genetic correlations
print(plot_G_ellipses(stan_model, Y, n_show = 100, X = X))

# Realized breeding values, estimated and simulated genetic correlations
par(mfrow = c(1, 3))
corrplot.mixed(cor(a), upper = "ellipse", mar = c(0, 0, 2, 0),
               title = "Realized Breeding Values Correlations")
corrplot.mixed(draws_mean(stan_model, "corrG"), upper = "ellipse", mar = c(0, 0, 2, 0),
               title = "Estimated Genetic Correlations")
corrplot.mixed(cov2cor(G), upper = "ellipse", mar = c(0, 0, 2, 0),
               title = "True Genetic Correlations")
par(mfrow = c(1, 1))

# Heatmaps of the estimated and simulated genetic correlations
print(plot_corr(draws_mean(stan_model, "corrG"), "Estimated") +
        plot_corr(cov2cor(G), "True"))

# Posterior distribution of the eigenvalue variance of G
G_draws <- draws$G                        # traits x traits x draws
eig_var <- apply(G_draws, 3, CalcEigenVar)
eig_var_true <- CalcEigenVar(G)
print(ggplot2::ggplot(data.frame(value = as.vector(eig_var)), ggplot2::aes(value)) +
        ggplot2::geom_density(fill = "grey70", colour = NA) +
        ggplot2::geom_vline(xintercept = eig_var_true, colour = "red", linetype = "dashed") +
        ggplot2::labs(x = "Eigenvalue variance of G", y = "Density",
                      title = "Posterior eigenvalue variance of G (dashed line = true value)"))

dev.off()

# Same plot on screen, after the PDF device is closed
print(plot_G_ellipses(stan_model, Y, n_show = 100, X = X))
