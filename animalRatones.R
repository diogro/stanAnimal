list_pkgs <- c("dplyr", "bayesplot", "posterior", "cmdstanr")
new_pkgs <- list_pkgs[!(list_pkgs %in% installed.packages()[,"Package"])]
if(length(new_pkgs) > 0){ install.packages(new_pkgs) }

library(stanAnimal)
ratonesdf <- read.table("https://raw.githubusercontent.com/diogro/evofencom/refs/heads/main/Tutoriais/ratones.tsv", header = TRUE)

# lmm_multiCov estimates the phenotypic covariance matrix P after fixed effects;
# there is no relatedness, so no G or breeding values.
traits = names(dplyr::select(ratonesdf, IS_PM:NA_BR))
Y = as.matrix(ratonesdf[, traits])
X = model.matrix(~ SEX + AGE + LIN, data = ratonesdf)

cov_fit = lmm_multiCov(Y, X, lkj_prior = 2, chains = 4, iter = 800,
                       control = list(adapt_delta = 0.99))

# Summaries and sampler diagnostics
cov_fit$summary(c("P[1,1]", "P[2,2]", "P[3,3]", "corrP[1,2]", "corrP[1,3]", "corrP[2,3]"))
cov_fit$diagnostic_summary()

# Fixed effects, named after the columns of X
beta_draws = animal_draws(cov_fit, "beta")
posterior::summarise_draws(beta_draws)

# Posterior mean P, compared with the covariance of the residuals
k = ncol(Y)
P_hat = matrix(cov_fit$summary("P", "mean")$mean, k)
dimnames(P_hat) = list(traits, traits)
round(P_hat, 3)
round(cov(residuals(lm(Y ~ X - 1))), 3)

library(evolqg)
MatrixCompare(P_hat, cov(residuals(lm(Y ~ X - 1))))


# Plots, all written to a single PDF
pdf("animalRatones_plots.pdf", width = 10, height = 7)
bayesplot::color_scheme_set("blue")
print(bayesplot::mcmc_trace(cov_fit$draws(c("P[1,1]", "P[2,2]", "P[3,3]"))))
print(bayesplot::mcmc_intervals(beta_draws))
print(bayesplot::mcmc_intervals(
  cov_fit$draws(), pars = c("P[1,1]", "P[2,2]", "P[3,3]"),
  prob = 0.8, prob_outer = 0.99, point_est = "mean"))
print(bayesplot::mcmc_intervals(
  cov_fit$draws(), pars = c("corrP[1,2]", "corrP[1,3]", "corrP[2,3]"),
  prob = 0.8, prob_outer = 0.99, point_est = "mean"))
print(plot_corr(cov2cor(P_hat), "Posterior mean phenotypic correlations"))
dev.off()
