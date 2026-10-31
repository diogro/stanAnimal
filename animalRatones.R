list_pkgs <- c("dplyr", "bayesplot", "posterior", "cmdstanr")
new_pkgs <- list_pkgs[!(list_pkgs %in% installed.packages()[,"Package"])]
if(length(new_pkgs) > 0){ install.packages(new_pkgs) }

library(stanAnimal)
ratonesdf <- read.table("https://raw.githubusercontent.com/diogro/evofencom/refs/heads/main/Tutoriais/ratones.tsv", header = TRUE)

# lmm_multiCov estimates the phenotypic covariance matrix P after fixed effects;
# there is no relatedness, so no G or breeding values.
all_traits = names(ratonesdf)[12:46]
traits = names(dplyr::select(ratonesdf, IS_PM:PM_ZI))

# animal_design() builds Y (one column per trait) and X from the data frame,
# so the model can be written as a formula. With many traits the cbind()
# response is assembled from the column names.
ratones_formula = stats::reformulate(c("SEX", "AGE", "LIN"),
                                     response = sprintf("cbind(%s)",
                                                        paste(traits, collapse = ", ")))
ratones = animal_design(ratones_formula, data = ratonesdf)
Y = ratones$Y
X = ratones$X
cat(nrow(Y), "individuals,", ncol(Y), "traits,", ncol(X), "fixed effects\n")

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
