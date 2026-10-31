list_pkgs <- c("dplyr", "bayesplot", "posterior", "cmdstanr", "patchwork")
new_pkgs <- list_pkgs[!(list_pkgs %in% installed.packages()[,"Package"])]
if(length(new_pkgs) > 0){ install.packages(new_pkgs) }
invisible(lapply(sub("^.*/", "", list_pkgs), library, character.only = TRUE))


library(stanAnimal)
ratonesdf <- read.table("https://raw.githubusercontent.com/diogro/evofencom/refs/heads/main/Tutoriais/ratones.tsv", header = TRUE)

# lmm_multiCov estimates the phenotypic covariance matrix P after fixed effects;
# there is no relatedness, so no G or breeding values.
all_traits = names(ratonesdf)[12:46]
traits = all_traits[1:35]  # subset to keep the example small (50s on full matrix)

# Size and shape reparametrization: the measurements share a single size
# gradient, so the model is fit to log size and orthogonal shape axes, and every
# draw is transformed back to the log-trait scale below. Keeping all k - 1 shape
# axes makes the reparametrization exact; a smaller n_shape gives a lower-rank,
# and much cheaper, model.
log_traits = log(as.matrix(ratonesdf[, traits]))
ss = size_shape(log_traits)
ratonesdf = cbind(ratonesdf, ss$scores)
cat(length(traits), "measured traits ->", ncol(ss$scores), "model traits\n")

# animal_design() builds Y (one column per model trait) and X from the data
# frame, so the model can be written as a formula. With many traits the cbind()
# response is assembled from the column names.
ratones_formula = stats::reformulate(c("SEX", "AGE", "LIN"),
                                     response = sprintf("cbind(%s)",
                                                        paste(ss$traits, collapse = ", ")))
ratones = animal_design(ratones_formula, data = ratonesdf)
Y = ratones$Y
X = ratones$X
cat(nrow(Y), "individuals,", ncol(Y), "model traits,", ncol(X), "fixed effects\n")

cov_fit = lmm_multiCov(Y, X, lkj_prior = 6, iter = 800, 
                       warmup = 400, cores = 4,  chains = 4)

# Summaries and sampler diagnostics (model scale: size and shape axes)
cov_fit$summary(c("sigma_P[1]", "P[1,1]", "P[2,2]", "corrP[1,2]"))
cov_fit$diagnostic_summary()

# Fixed effects, named after the columns of X (model scale: one coefficient per
# size or shape trait)
beta_draws = animal_draws(cov_fit, "beta")
posterior::summarise_draws(beta_draws)

# Back-transform every draw to the log-trait scale
post = back_transform(cov_fit, ss)
k = length(traits)
P_names = sprintf("P[%d,%d]", rep(1:k, k), rep(1:k, each = k))
post_draws = as.matrix(post)
posterior::summarise_draws(posterior::subset_draws(
  post, variable = c("sigma_P[1]", "sigma_P[2]", "corrP[1,2]", "corrP[1,3]", "corrP[2,3]")))

# Posterior mean P on the log-trait scale, compared with the residual covariance
P_hat = matrix(colMeans(post_draws[, P_names]), k, dimnames = list(traits, traits))
round(P_hat, 3)
round(cov(residuals(lm(log_traits ~ X - 1))), 3)

library(evolqg)
MatrixCompare(P_hat, cov(residuals(lm(log_traits ~ X - 1))))


# Plots, all written to a single PDF
pdf("animalRatones_plots.pdf", width = 10, height = 7)
bayesplot::color_scheme_set("blue")
print(bayesplot::mcmc_trace(cov_fit$draws(c("P[1,1]", "P[2,2]", "P[3,3]"))))
print(bayesplot::mcmc_intervals(beta_draws))
print(bayesplot::mcmc_intervals(
  posterior::subset_draws(post, variable = c("sigma_P[1]", "sigma_P[2]", "sigma_P[3]")),
  prob = 0.8, prob_outer = 0.99, point_est = "mean"))
print(bayesplot::mcmc_intervals(
  posterior::subset_draws(post, variable = c("corrP[1,2]", "corrP[1,3]", "corrP[2,3]")),
  prob = 0.8, prob_outer = 0.99, point_est = "mean"))
print(plot_corr(cov2cor(P_hat), "Posterior mean phenotypic correlations") +
plot_corr(cov2cor(cov(residuals(lm(log_traits ~ X - 1)))), "Residual correlations"))
dev.off()
